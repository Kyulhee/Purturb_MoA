"""
Run 12: Cross-CT Epistasis Transfer — Full Ablation + Coverage + R2/Pearson
============================================================================
Based on run_11, with key additions:
  1. Coverage (90% CI) measurement for UQ calibration
  2. R2/Pearson reported alongside Spearman rho
  3. Full ablation: A1-A7 + B1(CPA)
  4. MC-only UQ (ICM violation removed per run_11 finding)
  5. Residual rank transfer as primary (A7 confirmed best)
"""
import sys, io, os, json, time, copy
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np, torch, torch.nn as nn, torch.nn.functional as F
from sklearn.metrics import r2_score
from scipy.stats import spearmanr, pearsonr
import scanpy as sc, warnings
warnings.filterwarnings("ignore")
torch.manual_seed(42); np.random.seed(42)
RESULTS_DIR = os.path.dirname(os.path.abspath(__file__))
def P(*a,**k): print(*a,**k,flush=True)
def fmt(v): return f"{v:.3f}" if isinstance(v,(int,float,np.floating)) else str(v)

# ===== Data (from cached h5ad — no gears dependency) =====
def load_replogle(dp='outputs/analysis/run_04/data/gears_data',mx=200):
    P("  Loading K562 h5ad..."); ak=sc.read_h5ad(os.path.join(dp,'replogle_k562_essential','perturb_processed.h5ad'))
    P("  Loading RPE1 h5ad..."); ar=sc.read_h5ad(os.path.join(dp,'replogle_rpe1_essential','perturb_processed.h5ad'))
    sh=set(ak.obs['condition'].unique())&set(ar.obs['condition'].unique()); sh.discard('ctrl')
    gn=sorted(set(ak.var_names)&set(ar.var_names)); kp=sh|{'ctrl'}
    ak=ak[ak.obs['condition'].isin(kp),gn].copy(); ar=ar[ar.obs['condition'].isin(kp),gn].copy()
    ak=_sub(ak,mx); ar=_sub(ar,mx); ak.obs['cell_type']='K562'; ar.obs['cell_type']='RPE1'
    ac=ak.concatenate(ar,batch_key='batch'); P(f"  Combined: {ac.shape}")
    return ac,sh

def load_norman(dp='outputs/analysis/run_04/data/gears_data'):
    return sc.read_h5ad(os.path.join(dp,'norman','perturb_processed.h5ad'))

def _sub(ad,mx):
    idx=[]
    for c in ad.obs['condition'].unique():
        m=np.where(ad.obs['condition']==c)[0]
        idx.extend(np.random.choice(m,min(len(m),mx),replace=False)if len(m)>mx else m)
    return ad[sorted(idx)].copy()

def prep(ad,nh=500):
    sc.pp.filter_genes(ad,min_cells=50); sc.pp.normalize_total(ad,1e4); sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad,n_top_genes=nh,flavor='seurat'); ad=ad[:,ad.var.highly_variable].copy()
    X=ad.X.toarray().astype(np.float32)if hasattr(ad.X,'toarray')else ad.X.astype(np.float32)
    P(f"  HVG: {ad.shape}"); return ad,X

def g2c(cu):
    r={}
    for c in cu:
        if c=='ctrl': continue
        if '+ctrl' in c: r[c.replace('+ctrl','')]=c
        elif 'ctrl+' in c: r[c.replace('ctrl+','')]=c
        elif '+' not in c: r[c]=c
    return r

# ===== Model =====
class Enc(nn.Module):
    def __init__(s,ng,np_,zd,nct):
        super().__init__(); s.zd=zd
        s.xe=nn.Sequential(nn.Linear(ng,256),nn.ReLU(),nn.Dropout(.1),nn.Linear(256,128),nn.ReLU())
        s.pe=nn.Embedding(np_,zd); s.zx=nn.Linear(128+nct,zd*2); s.zt=nn.Linear(128+zd,zd*2); s.ztx=nn.Linear(128+zd,zd*2)
    def forward(s,x,pid,ctoh):
        h=s.xe(x); zi=s.pe(pid)
        a=s.zx(torch.cat([h,ctoh],-1)); b=s.zt(torch.cat([h,zi],-1)); c=s.ztx(torch.cat([h,zi],-1))
        return(a[:,:s.zd],a[:,s.zd:]),(b[:,:s.zd],b[:,s.zd:]),(c[:,:s.zd],c[:,s.zd:])

class Dec(nn.Module):
    def __init__(s,zd,ng,dr=.1):
        super().__init__()
        s.d=nn.Sequential(nn.Linear(3*zd,256),nn.ReLU(),nn.Dropout(dr),nn.Linear(256,128),nn.ReLU(),nn.Dropout(dr),nn.Linear(128,ng))
    def forward(s,zx,zt,ztx): return s.d(torch.cat([zx,zt,ztx],-1))

def rp(m,lv): return m+torch.exp(.5*lv)*torch.randn_like(m)

def icm_reg(zm,cts):
    ut=torch.unique(cts)
    if len(ut)<2: return((zm.var(0)-1)**2).mean()
    loss=torch.tensor(0.,device=zm.device)
    for i in range(len(ut)):
        for j in range(i+1,len(ut)):
            zi,zj=zm[cts==ut[i]],zm[cts==ut[j]]
            loss+=(zi.mean(0)-zj.mean(0)).pow(2).sum()
            ns=min(50,zi.shape[0],zj.shape[0])
            if ns>5:
                si,sj=zi[:ns],zj[:ns]; sig=1.
                xx=torch.exp(-torch.cdist(si,si).pow(2)/(2*sig)).mean()
                yy=torch.exp(-torch.cdist(sj,sj).pow(2)/(2*sig)).mean()
                xy=torch.exp(-torch.cdist(si,sj).pow(2)/(2*sig)).mean()
                loss+=xx+yy-2*xy
    return loss

def train_fcr(X,pi,ci,np_,nct,zd=8,icm=False,ep=150,bs=512,iw=10.,mc_drop=0.1):
    ng=X.shape[1]; enc=Enc(ng,np_,zd,nct); dec=Dec(zd,ng,dr=mc_drop)
    opt=torch.optim.Adam(list(enc.parameters())+list(dec.parameters()),lr=1e-3)
    xt=torch.FloatTensor(X); pt=torch.LongTensor(pi); ct=torch.LongTensor(ci)
    ch=F.one_hot(ct,nct).float()
    dl=torch.utils.data.DataLoader(torch.utils.data.TensorDataset(xt,pt,ct,ch),batch_size=bs,shuffle=True,drop_last=True)
    for e in range(ep):
        el=0
        for bx,bp,bc,bch in dl:
            opt.zero_grad()
            (zxm,zxlv),(ztm,ztlv),(ztxm,ztxlv)=enc(bx,bp,bch)
            xr=dec(rp(zxm,zxlv),rp(ztm,ztlv),rp(ztxm,ztxlv))
            rl=F.mse_loss(xr,bx,reduction='sum')
            kl=sum(-.5*torch.sum(1+lv-m.pow(2)-lv.exp())for m,lv in[(zxm,zxlv),(ztm,ztlv),(ztxm,ztxlv)])
            loss=rl+.5*kl
            if icm and nct>1: loss=loss+iw*icm_reg(ztxm,bc)
            loss.backward(); opt.step(); el+=loss.item()
        if(e+1)%30==0: P(f"    Ep {e+1}/{ep}: loss={el/len(dl):.1f}")
    return enc,dec

# ===== CPA Baseline (B1) =====
class CPAEnc(nn.Module):
    def __init__(s,ng,np_,nct,zd=8):
        super().__init__()
        s.xe=nn.Sequential(nn.Linear(ng,256),nn.ReLU(),nn.Linear(256,128),nn.ReLU())
        s.pe=nn.Embedding(np_,zd); s.ce=nn.Embedding(nct,zd)
        s.mu=nn.Linear(128+zd+zd,zd); s.lv=nn.Linear(128+zd+zd,zd)
    def forward(s,x,pid,ctid):
        h=s.xe(x); pi=s.pe(pid); ci=s.ce(ctid)
        inp=torch.cat([h,pi,ci],-1); return s.mu(inp),s.lv(inp)

class CPADec(nn.Module):
    def __init__(s,zd,ng):
        super().__init__()
        s.d=nn.Sequential(nn.Linear(zd,256),nn.ReLU(),nn.Linear(256,128),nn.ReLU(),nn.Linear(128,ng))
    def forward(s,z): return s.d(z)

def train_cpa(X,pi,ci,np_,nct,zd=8,ep=150,bs=512):
    ng=X.shape[1]; enc=CPAEnc(ng,np_,nct,zd); dec=CPADec(zd,ng)
    opt=torch.optim.Adam(list(enc.parameters())+list(dec.parameters()),lr=1e-3)
    xt=torch.FloatTensor(X); pt=torch.LongTensor(pi); ctt=torch.LongTensor(ci)
    dl=torch.utils.data.DataLoader(torch.utils.data.TensorDataset(xt,pt,ctt),batch_size=bs,shuffle=True,drop_last=True)
    for e in range(ep):
        el=0
        for bx,bp,bc in dl:
            opt.zero_grad(); mu,lv=enc(bx,bp,bc)
            z=mu+torch.exp(.5*lv)*torch.randn_like(mu); xr=dec(z)
            rl=F.mse_loss(xr,bx,reduction='sum')
            kl=-0.5*torch.sum(1+lv-mu.pow(2)-lv.exp())
            loss=rl+.5*kl; loss.backward(); opt.step(); el+=loss.item()
        if(e+1)%30==0: P(f"    CPA Ep {e+1}/{ep}: loss={el/len(dl):.1f}")
    return enc,dec

# ===== Epistasis Scores (3-Formula + MC Dropout UQ + Coverage) =====
def compute_epi_scores(enc,dec,ad,X,pim,ctm,nct,zd,n_mc=30):
    enc.eval(); dec.eval()
    conds=ad.obs['condition'].values; cts=ad.obs['cell_type'].values; ctn=sorted(ad.obs['cell_type'].unique())
    ctrl_m={}
    for ct in ctn:
        m=(conds=='ctrl')&(cts==ct)
        if m.sum()>=5: ctrl_m[ct]=X[m].mean(0)
    pd_={}
    sps=sorted([c for c in np.unique(conds) if c!='ctrl'and('+'not in c or c.endswith('+ctrl')or c.startswith('ctrl+'))])
    for ct in ctn:
        ci=ctm[ct]
        for pn in sps:
            if pn not in pim: continue
            m=(conds==pn)&(cts==ct)
            if m.sum()<10: continue
            nu=min(50,m.sum()); xp=torch.FloatTensor(X[m][:nu])
            pt=torch.full((nu,),pim[pn],dtype=torch.long)
            ch=F.one_hot(torch.full((nu,),ci,dtype=torch.long),nct).float()
            with torch.no_grad(): (zxm,_),(ztm,_),(ztxm,_)=enc(xp,pt,ch)
            if pn not in pd_: pd_[pn]={}
            pd_[pn][ct]={'z_tx':ztxm.mean(0),'z_x':zxm.mean(0),'z_t':ztm.mean(0),'y_mean':X[m].mean(0)}

    scores={}
    for ct in ctn:
        if ct not in ctrl_m: continue
        cm=ctrl_m[ct]; scores[ct]={}
        for pn in sps:
            if pn not in pd_ or ct not in pd_[pn]: continue
            d=pd_[pn][ct]; y_obs=d['y_mean']
            zxr=d['z_x'].unsqueeze(0); ztr=d['z_t'].unsqueeze(0); ztxr=d['z_tx'].unsqueeze(0)
            with torch.no_grad(): y_pred=dec(zxr,ztr,ztxr)[0].numpy()
            # MC Dropout (n_mc samples) — MC-only UQ
            dec.train(); mcs=[]
            with torch.no_grad():
                for _ in range(n_mc): mcs.append(dec(zxr,ztr,ztxr)[0].numpy())
            dec.eval(); mcs=np.array(mcs)  # (n_mc, ng)
            mc_var=float(np.var(mcs,axis=0).mean())
            mc_mean=mcs.mean(axis=0)
            pred_err=float(np.abs(y_obs-y_pred).mean())
            # Coverage: 90% CI — what fraction of genes have y_obs within MC 5-95 percentile
            mc_lo=np.percentile(mcs,5,axis=0)  # (ng,)
            mc_hi=np.percentile(mcs,95,axis=0)  # (ng,)
            in_ci=((y_obs>=mc_lo)&(y_obs<=mc_hi)).astype(float)
            coverage_90=float(in_ci.mean())
            # 3 formulas
            e_add=pred_err
            e_mult=float(np.abs(y_obs-y_pred).mean()/(np.abs(y_pred).mean()+1e-8))
            pn_=float(np.abs(y_obs*cm-y_pred*cm).mean())/(float(np.abs(y_obs*cm).mean())+1e-8)
            # Self-combo residual (A7)
            r_add=y_obs-(2*y_pred-cm); r_add_mag=float(np.abs(r_add).mean())
            # R2 and Pearson per perturbation
            r2_pert=float(r2_score(y_obs,y_pred))
            pr_pert,_=pearsonr(y_obs,y_pred)
            scores[ct][pn]={'add':e_add,'mult':e_mult,'prod':pn_,'pred_err':pred_err,
                            'mc_var':mc_var,'r_add_mag':r_add_mag,'coverage_90':coverage_90,
                            'r2':r2_pert,'pearson':pr_pert}
    return scores,pd_

# ===== Cross-CT Transfer (RQ3 Core) =====
def cross_ct_transfer(scores,src='K562',tgt='RPE1'):
    if src not in scores or tgt not in scores: return{'error':'missing CT'}
    sh=set(scores[src].keys())&set(scores[tgt].keys())
    if len(sh)<10: return{'error':f'too few: {len(sh)}'}
    r={'n_shared':len(sh)}
    # Per-formula rho, pearson, R2
    for f in['add','mult','prod','pred_err','mc_var','r_add_mag','coverage_90']:
        sv=[scores[src][p].get(f)for p in sh]; tv=[scores[tgt][p].get(f)for p in sh]
        # Skip if key missing in any entry (e.g. CPA has no coverage_90)
        if None in sv or None in tv: continue
        rho,pv=spearmanr(sv,tv); r[f'rho_{f}']=float(rho); r[f'p_{f}']=float(pv)
    # Composite (MC-only: best weight search)
    br=-2; bw=None
    for wm in np.arange(0,1.01,.2):
        for we in np.arange(0,1.01-wm,.2):
            wi=1-wm-we
            sv=[wm*scores[src][p]['mc_var']+wi*scores[src][p]['pred_err']+we*scores[src][p]['r_add_mag']for p in sh]
            tv=[wm*scores[tgt][p]['mc_var']+wi*scores[tgt][p]['pred_err']+we*scores[tgt][p]['r_add_mag']for p in sh]
            rho,_=spearmanr(sv,tv)
            if rho>br: br=rho; bw={'mc':float(wm),'err':float(wi),'r_add':float(we)}
    r['rho_composite']=float(br); r['composite_w']=bw
    # Top-k overlap
    sl=sorted(sh)
    for k in[10,20,50]:
        if k>len(sl): continue
        st=set(sorted(sl,key=lambda p:-scores[src][p]['add'])[:k])
        tt=set(sorted(sl,key=lambda p:-scores[tgt][p]['add'])[:k])
        ol=len(st&tt)/k; rnd=k/len(sl)
        r[f'top{k}_overlap']=float(ol); r[f'top{k}_random']=float(rnd); r[f'top{k}_imp']=float(ol/max(rnd,1e-8))
    # Aggregate R2 and Pearson (skip if not available, e.g. CPA)
    src_r2=[scores[src][p].get('r2')for p in sh]; tgt_r2=[scores[tgt][p].get('r2')for p in sh]
    src_pr=[scores[src][p].get('pearson')for p in sh]; tgt_pr=[scores[tgt][p].get('pearson')for p in sh]
    if None not in src_r2 and None not in tgt_r2:
        r['mean_r2_src']=float(np.mean(src_r2)); r['mean_r2_tgt']=float(np.mean(tgt_r2))
        r2_rho,_=spearmanr(src_r2,tgt_r2); r['rho_r2_cross']=float(r2_rho)
    if None not in src_pr and None not in tgt_pr:
        r['mean_pearson_src']=float(np.mean(src_pr)); r['mean_pearson_tgt']=float(np.mean(tgt_pr))
        pr_rho,_=spearmanr(src_pr,tgt_pr); r['rho_pearson_cross']=float(pr_rho)
    return r

# ===== A1 vs A2 Ablation =====
def ablation_icm(si,sn,src='K562',tgt='RPE1'):
    ti=cross_ct_transfer(si,src,tgt); tn=cross_ct_transfer(sn,src,tgt)
    r={'A1_icm':ti,'A2_noicm':tn}
    for f in['add','mult','prod','pred_err','mc_var','r_add_mag','composite']:
        ri=ti.get(f'rho_{f}',0); rn=tn.get(f'rho_{f}',0)
        r[f'imp_{f}']=float(ri/rn)if abs(rn)>1e-8 else(float('inf')if abs(ri)>1e-8 else 0.)
    for k in[10,20,50]:
        r[f'top{k}_icm']=ti.get(f'top{k}_overlap',0); r[f'top{k}_noicm']=tn.get(f'top{k}_overlap',0)
    return r

# ===== A7: Residual Rank Transfer =====
def transfer_res_rank(scores,src='K562',tgt='RPE1'):
    if src not in scores or tgt not in scores: return{'error':'missing CT'}
    sh=set(scores[src].keys())&set(scores[tgt].keys())
    if len(sh)<10: return{'error':'too few'}
    sv=[scores[src][p]['r_add_mag']for p in sh]; tv=[scores[tgt][p]['r_add_mag']for p in sh]
    rho,pv=spearmanr(sv,tv); pr,ppv=pearsonr(sv,tv)
    sl=sorted(sh); tk={}
    for k in[10,20,50]:
        if k>len(sl): continue
        st=set(sorted(sl,key=lambda p:-scores[src][p]['r_add_mag'])[:k])
        tt=set(sorted(sl,key=lambda p:-scores[tgt][p]['r_add_mag'])[:k])
        tk[f'top{k}_overlap']=float(len(st&tt)/k)
    return{'rho':float(rho),'pearson':float(pr),'p':float(pv),'n':len(sh),**tk,'method':'A7_residual_rank'}

# ===== A3: Single Formula =====
def transfer_single_formula(scores,src='K562',tgt='RPE1',formula='add'):
    """A3: Use only additive formula (no 3-formula sensitivity)"""
    if src not in scores or tgt not in scores: return{'error':'missing CT'}
    sh=set(scores[src].keys())&set(scores[tgt].keys())
    if len(sh)<10: return{'error':'too few'}
    sv=[scores[src][p][formula]for p in sh]; tv=[scores[tgt][p][formula]for p in sh]
    rho,pv=spearmanr(sv,tv)
    return{'rho':float(rho),'p':float(pv),'n':len(sh),'method':f'A3_single_{formula}'}

# ===== A4: Trivial Decomposition =====
def transfer_trivial(scores,src='K562',tgt='RPE1'):
    """A4: Use pred_err directly (no residual decomposition, trivial baseline)"""
    if src not in scores or tgt not in scores: return{'error':'missing CT'}
    sh=set(scores[src].keys())&set(scores[tgt].keys())
    if len(sh)<10: return{'error':'too few'}
    sv=[scores[src][p]['pred_err']for p in sh]; tv=[scores[tgt][p]['pred_err']for p in sh]
    rho,pv=spearmanr(sv,tv)
    return{'rho':float(rho),'p':float(pv),'n':len(sh),'method':'A4_trivial_decomp'}

# ===== RQ2: Holdout UQ with Coverage =====
def holdout_eval(enc,dec,ad,X,pim,ctm,nct,zd,n_ho=50,n_mc=30):
    enc.eval(); dec.eval()
    conds=ad.obs['condition'].values; cts=ad.obs['cell_type'].values; ctn=sorted(ad.obs['cell_type'].unique())
    sps=[c for c in np.unique(conds) if c!='ctrl'and('+'not in c or c.endswith('+ctrl')or c.startswith('ctrl+'))]
    np.random.shuffle(sps); ho=sps[:n_ho]
    czx={}
    for ct in ctn:
        ci=ctm[ct]; m=(conds=='ctrl')&(cts==ct)
        if m.sum()<5: continue
        nu=min(50,m.sum()); xc=torch.FloatTensor(X[m][:nu])
        pt=torch.full((nu,),pim['ctrl'],dtype=torch.long)
        ch=F.one_hot(torch.full((nu,),ci,dtype=torch.long),nct).float()
        with torch.no_grad():(a,_,_)=enc(xc,pt,ch); czx[ct]=a[0].mean(0)

    res=[]
    for pn in ho:
        if pn not in pim: continue
        for ct in ctn:
            ci=ctm[ct]; m=(conds==pn)&(cts==ct)
            if m.sum()<10 or ct not in czx: continue
            yt=X[m].mean(0); nu=min(30,m.sum())
            xp=torch.FloatTensor(X[m][:nu]); pt=torch.full((nu,),pim[pn],dtype=torch.long)
            ch=F.one_hot(torch.full((nu,),ci,dtype=torch.long),nct).float()
            with torch.no_grad():(zxm,_),(ztm,_),(ztxm,_)=enc(xp,pt,ch)
            zxc=czx[ct].unsqueeze(0).expand(nu,-1)
            with torch.no_grad(): yp=dec(zxc,ztm,ztxm)
            ypm=yp.mean(0).numpy()
            # MC Dropout for UQ + Coverage
            dec.train(); mcp=[]
            with torch.no_grad():
                zs=czx[ct].unsqueeze(0); zsm=ztm.mean(0).unsqueeze(0); ztm2=ztxm.mean(0).unsqueeze(0)
                for _ in range(n_mc): mcp.append(dec(zs,zsm,ztm2).numpy()[0])
            dec.eval(); mcp=np.array(mcp)
            mc_var=float(np.var(mcp,axis=0).mean())
            # Coverage: 90% CI
            mc_lo=np.percentile(mcp,5,axis=0); mc_hi=np.percentile(mcp,95,axis=0)
            coverage=float(((yt>=mc_lo)&(yt<=mc_hi)).mean())
            res.append({'pert':pn,'ct':ct,'error':float(np.abs(yt-ypm).mean()),
                        'r2':float(r2_score(yt,ypm)),'mc_var':mc_var,
                        'coverage_90':coverage,
                        'pearson':float(pearsonr(yt,ypm)[0])})
    return res

def uq_holdout(hr):
    if len(hr)<5: return{'error':'too few'}
    e=np.array([r['error']for r in hr]); mc=np.array([r['mc_var']for r in hr])
    cov=np.array([r['coverage_90']for r in hr]); r2=np.array([r['r2']for r in hr])
    pr=np.array([r['pearson']for r in hr])
    rmc,pmc=spearmanr(mc,e)
    # MC-only UQ (no ICM violation combination per run_11 finding)
    return{'n':len(hr),'mean_r2':float(np.mean(r2)),'mean_pearson':float(np.mean(pr)),
           'rho_mc':float(rmc),'p_mc':float(pmc),
           'mean_coverage_90':float(np.mean(cov)),'std_coverage_90':float(np.std(cov)),
           'coverage_in_range':bool(0.85<=np.mean(cov)<=0.95)}

# ===== AL Simulation (RQ4) =====
def al_sim(scores,src='K562',tgt='RPE1'):
    if src not in scores or tgt not in scores: return{'error':'missing CT'}
    sh=sorted(set(scores[src].keys())&set(scores[tgt].keys()))
    if len(sh)<20: return{'error':'too few'}
    tgt_epi={p:scores[tgt][p]['add']for p in sh}
    se=sorted(tgt_epi.values(),reverse=True); thr=se[int(len(se)*.3)]
    is_s={p:tgt_epi[p]>=thr for p in sh}; ns=sum(is_s.values())
    np.random.seed(42); ro=list(np.random.permutation(sh))
    uqo=sorted(sh,key=lambda p:-(scores[src][p]['mc_var']+scores[src][p]['pred_err']))
    io=sorted(sh,key=lambda p:-scores[src][p]['icm_viol']if'icm_viol'in scores[src][p]else 0)
    eo=sorted(sh,key=lambda p:-scores[src][p]['add'])
    oo=sorted(sh,key=lambda p:-tgt_epi[p])
    def cum(o):
        f=0;c=[]
        for x in o: f+=is_s[x]; c.append(f/max(ns,1))
        return c
    rc,uc,ic,ec,oc=cum(ro),cum(uqo),cum(io),cum(eo),cum(oo)
    tk={}
    for k in[5,10,20]:
        if k<=len(sh):
            tk[k]={'rand':rc[k-1],'uq':uc[k-1],'epi':ec[k-1],'oracle':oc[k-1],
                    'imp_uq':uc[k-1]/max(rc[k-1],1e-8),'imp_epi':ec[k-1]/max(rc[k-1],1e-8)}
    s20=set(eo[:min(20,len(eo))]); te=sorted(tgt_epi.values(),reverse=True); tt=te[int(len(te)*.3)]
    ts={k for k,v in tgt_epi.items()if v>=tt}
    return{'n':len(sh),'n_strong':ns,'top_k':tk,
           'transfer_overlap':float(len(s20&ts)/max(len(s20),1))}

# ===== Norman Epistasis Precision =====
def norman_precision(ad_n,X_n):
    conds=ad_n.obs['condition'].values; gc=g2c(np.unique(conds))
    sps,dps=[],[]
    for c in np.unique(conds):
        if c=='ctrl': continue
        ps=c.split('+')
        if len(ps)==1: sps.append(c)
        elif len(ps)==2:
            if ps[0]=='ctrl'or ps[1]=='ctrl': sps.append(c)
            else: dps.append(c)
    cm_=ad_n[conds=='ctrl'].X; cm=cm_.toarray().mean(0)if hasattr(cm_,'toarray')else cm_.mean(0)
    ac=sorted(ad_n.obs['condition'].unique()); pm={c:i for i,c in enumerate(ac)}; np_=len(ac)
    ci=np.zeros(len(conds),dtype=np.int64); pi=np.array([pm[c]for c in conds],dtype=np.int64)
    P("  Training Norman FCR..."); enc,dec=train_fcr(X_n,pi,ci,np_,1,zd=8,icm=False,ep=100)
    enc.eval(); dec.eval()
    se={}
    for g,cn in gc.items():
        if cn not in pm: continue
        m=conds==cn
        if m.sum()<5: continue
        nu=min(50,m.sum()); xp=torch.FloatTensor(X_n[m][:nu])
        pt=torch.full((nu,),pm[cn],dtype=torch.long); ch=F.one_hot(torch.zeros(nu,dtype=torch.long),1).float()
        with torch.no_grad():(zxm,_),(ztm,_),(ztxm,_)=enc(xp,pt,ch)
        se[g]={'z_tx':ztxm.mean(0),'z_x':zxm.mean(0),'z_t':ztm.mean(0),'y_mean':X_n[m].mean(0)}

    res=[]
    for dp in dps:
        p1,p2=dp.split('+')
        if p1 not in se or p2 not in se: continue
        dm=conds==dp
        if dm.sum()<5: continue
        yr=X_n[dm].mean(0); yA=se[p1]['y_mean']; yB=se[p2]['y_mean']
        yad=yA+yB-cm; rad=yr-yad
        with torch.no_grad():
            ztc=se[p1]['z_tx']+se[p2]['z_tx']
            yp=dec(se[p1]['z_x'].unsqueeze(0),se[p1]['z_t'].unsqueeze(0),ztc.unsqueeze(0))[0].numpy()
        # MC Dropout for coverage on Norman
        dec.train(); mcs=[]
        with torch.no_grad():
            for _ in range(30): mcs.append(dec(se[p1]['z_x'].unsqueeze(0),se[p1]['z_t'].unsqueeze(0),ztc.unsqueeze(0))[0].numpy()[0])
        dec.eval(); mcs=np.array(mcs)
        mc_lo=np.percentile(mcs,5,axis=0); mc_hi=np.percentile(mcs,95,axis=0)
        cov_90=float(((yr>=mc_lo)&(yr<=mc_hi)).mean())
        res.append({'dp':dp,'r_add_mag':float(np.abs(rad).mean()),'pred_res_mag':float(np.abs(yr-yp).mean()),
                     'coverage_90':cov_90})

    if len(res)<5: return{'error':f'too few doubles: {len(res)}'}
    em=np.array([r['r_add_mag']for r in res]); pm_=np.array([r['pred_res_mag']for r in res])
    cv=np.array([r['coverage_90']for r in res])
    rho,pv=spearmanr(pm_,em); thr=np.percentile(em,70); is_e=em>=thr; ne=int(is_e.sum())
    pc={}
    for kf in[.1,.2,.3]:
        k=max(1,int(len(res)*kf)); tk=np.argsort(-pm_)[:k]; p_=float(is_e[tk].mean())
        r_=float(ne/len(res)); pc[f'prec_top{int(kf*100)}']={'precision':p_,'random':r_,'imp':float(p_/max(r_,1e-8))}
    return{'n_doubles':len(res),'n_epi':ne,'trivial_prec':float(ne/len(res)),
           'rho_pred_actual':float(rho),'p_pred_actual':float(pv),
           'mean_coverage_90':float(np.mean(cv)),'precision':pc}

# ===== CPA Epistasis (B1) =====
def compute_cpa_epi(cpa_enc,cpa_dec,ad,X,pim,ctm,nct,zd):
    """B1: CPA epistasis scores — single z, no factorized decomposition"""
    cpa_enc.eval(); cpa_dec.eval()
    conds=ad.obs['condition'].values; cts=ad.obs['cell_type'].values; ctn=sorted(ad.obs['cell_type'].unique())
    pd_={}
    sps=sorted([c for c in np.unique(conds) if c!='ctrl'and('+'not in c or c.endswith('+ctrl')or c.startswith('ctrl+'))])
    for ct in ctn:
        ci=ctm[ct]
        for pn in sps:
            if pn not in pim: continue
            m=(conds==pn)&(cts==ct)
            if m.sum()<10: continue
            nu=min(50,m.sum()); xp=torch.FloatTensor(X[m][:nu])
            pt=torch.full((nu,),pim[pn],dtype=torch.long); ctt=torch.full((nu,),ci,dtype=torch.long)
            with torch.no_grad(): mu,lv=cpa_enc(xp,pt,ctt); yp=cpa_dec(mu)
            if pn not in pd_: pd_[pn]={}
            pd_[pn][ct]={'y_mean':X[m].mean(0),'y_pred':yp.mean(0).numpy()}
    scores={}
    for ct in ctn:
        scores[ct]={}
        for pn in sps:
            if pn not in pd_ or ct not in pd_[pn]: continue
            d=pd_[pn][ct]; y_obs=d['y_mean']; y_pred=d['y_pred']
            err=float(np.abs(y_obs-y_pred).mean())
            scores[ct][pn]={'add':err,'mult':float(err/(np.abs(y_pred).mean()+1e-8)),
                            'prod':err,'pred_err':err,'mc_var':0.0,'r_add_mag':err}
    return scores

# ===== Main =====
def main():
    t0=time.time()
    P("="*70+"\nRun 12: Cross-CT Epistasis Transfer — Full Ablation + Coverage\n"+"="*70)

    # Phase 1: Load Replogle
    P("\n[Phase 1] Loading Replogle K562+RPE1...")
    ad,sh=load_replogle(); ad,X=prep(ad); ng=X.shape[1]
    ctn=sorted(ad.obs['cell_type'].unique()); ctm={n:i for i,n in enumerate(ctn)}; nct=len(ctn)
    ap=sorted(ad.obs['condition'].unique()); pim={n:i for i,n in enumerate(ap)}; np_=len(ap)
    P(f"  Genes: {ng}, Perts: {np_}, CTs: {nct}, Shared: {len(sh)}")
    ci=np.array([ctm[c]for c in ad.obs['cell_type'].values],dtype=np.int64)
    pi=np.array([pim[c]for c in ad.obs['condition'].values],dtype=np.int64)

    # ===== A1: Full model (ICM + residual rank transfer) =====
    P("\n[Config A1] FCR + ICM...")
    e1,d1=train_fcr(X,pi,ci,np_,nct,icm=True,iw=10.)

    # ===== A2: No ICM =====
    P("\n[Config A2] FCR baseline (no ICM)...")
    e2,d2=train_fcr(X,pi,ci,np_,nct,icm=False)

    # ===== B1: CPA Baseline =====
    P("\n[Config B1] CPA baseline...")
    cpa_enc,cpa_dec=train_cpa(X,pi,ci,np_,nct)

    # ===== RQ2: Holdout UQ with Coverage =====
    P("\n[RQ2] Holdout UQ + Coverage (A1 ICM model)...")
    hr=holdout_eval(e1,d1,ad,X,pim,ctm,nct,8,n_ho=50,n_mc=30)
    uq=uq_holdout(hr)
    P(f"  n={uq.get('n','N/A')}, rho_mc={fmt(uq.get('rho_mc','N/A'))}, coverage_90={fmt(uq.get('mean_coverage_90','N/A'))}")
    P(f"  mean_r2={fmt(uq.get('mean_r2','N/A'))}, mean_pearson={fmt(uq.get('mean_pearson','N/A'))}")
    P(f"  coverage_in_range={uq.get('coverage_in_range','N/A')}")

    P("\n[RQ2] Holdout UQ + Coverage (A2 no-ICM model)...")
    hr2=holdout_eval(e2,d2,ad,X,pim,ctm,nct,8,n_ho=50,n_mc=30)
    uq2=uq_holdout(hr2)
    P(f"  n={uq2.get('n','N/A')}, rho_mc={fmt(uq2.get('rho_mc','N/A'))}, coverage_90={fmt(uq2.get('mean_coverage_90','N/A'))}")

    # ===== Epistasis Scores =====
    P("\n[Phase 2] Computing epistasis scores (A1 ICM)...")
    si,_=compute_epi_scores(e1,d1,ad,X,pim,ctm,nct,8,n_mc=30)
    P("\n[Phase 2b] Computing epistasis scores (A2 no-ICM)...")
    sn,_=compute_epi_scores(e2,d2,ad,X,pim,ctm,nct,8,n_mc=30)

    # ===== RQ3: Cross-CT Transfer =====
    P("\n[RQ3] Cross-CT Epistasis Transfer...")
    ti=cross_ct_transfer(si,'K562','RPE1')
    tn=cross_ct_transfer(sn,'K562','RPE1')
    P(f"  A1(ICM):  rho_add={fmt(ti.get('rho_add','N/A'))}, rho_prod={fmt(ti.get('rho_prod','N/A'))}")
    P(f"  A2(noICM): rho_add={fmt(tn.get('rho_add','N/A'))}, rho_prod={fmt(tn.get('rho_prod','N/A'))}, rho_composite={fmt(tn.get('rho_composite','N/A'))}")
    P(f"  A1 coverage_90 rho={fmt(ti.get('rho_coverage_90','N/A'))}, A2 coverage_90 rho={fmt(tn.get('rho_coverage_90','N/A'))}")
    P(f"  A1 mean_r2_src={fmt(ti.get('mean_r2_src','N/A'))}, mean_r2_tgt={fmt(ti.get('mean_r2_tgt','N/A'))}")
    P(f"  A1 mean_pearson_src={fmt(ti.get('mean_pearson_src','N/A'))}, mean_pearson_tgt={fmt(ti.get('mean_pearson_tgt','N/A'))}")

    # ===== Ablation: A1 vs A2 =====
    P("\n[Ablation] A1(ICM) vs A2(no-ICM)...")
    abl=ablation_icm(si,sn,'K562','RPE1')
    for f in['add','mult','prod','pred_err','mc_var','r_add_mag','composite']:
        P(f"  imp_{f}: {fmt(abl.get(f'imp_{f}','N/A'))}")
    P(f"  ICM improvement on epistasis transfer: {fmt(abl.get('imp_add','N/A'))}× (expect ~1.0, negative result)")

    # ===== A7: Residual Rank Transfer =====
    P("\n[A7] Residual Rank Transfer...")
    a7=transfer_res_rank(si,'K562','RPE1')
    P(f"  rho={fmt(a7.get('rho','N/A'))}, pearson={fmt(a7.get('pearson','N/A'))}, p={fmt(a7.get('p','N/A'))}")
    for k in[10,20,50]:
        if f'top{k}_overlap' in a7: P(f"  top{k}_overlap={fmt(a7[f'top{k}_overlap'])}")

    # ===== A3: Single Formula (additive only) =====
    P("\n[A3] Single Formula Ablation...")
    a3_add=transfer_single_formula(si,'K562','RPE1','add')
    a3_mult=transfer_single_formula(si,'K562','RPE1','mult')
    a3_prod=transfer_single_formula(si,'K562','RPE1','prod')
    P(f"  add_only:  rho={fmt(a3_add.get('rho','N/A'))}")
    P(f"  mult_only: rho={fmt(a3_mult.get('rho','N/A'))}")
    P(f"  prod_only: rho={fmt(a3_prod.get('rho','N/A'))}")

    # ===== A4: Trivial Decomposition =====
    P("\n[A4] Trivial Decomposition (pred_err only, no residual)...")
    a4=transfer_trivial(si,'K562','RPE1')
    P(f"  rho={fmt(a4.get('rho','N/A'))}, p={fmt(a4.get('p','N/A'))}")

    # ===== B1: CPA Baseline =====
    P("\n[B1] CPA Epistasis Evaluation...")
    cpa_scores=compute_cpa_epi(cpa_enc,cpa_dec,ad,X,pim,ctm,nct,8)
    cpa_transfer=cross_ct_transfer(cpa_scores,'K562','RPE1')
    P(f"  CPA rho_add={fmt(cpa_transfer.get('rho_add','N/A'))}, rho_prod={fmt(cpa_transfer.get('rho_prod','N/A'))}")
    P(f"  CPA rho_composite={fmt(cpa_transfer.get('rho_composite','N/A'))}")

    # ===== 3-Formula Sensitivity =====
    P("\n[3-Formula Sensitivity] Cross-CT rho across formulas...")
    for f in['add','mult','prod','pred_err','mc_var','r_add_mag','coverage_90']:
        ri=ti.get(f'rho_{f}',0); rn=tn.get(f'rho_{f}',0)
        P(f"  {f}: A1={fmt(ri)}, A2={fmt(rn)}")

    # ===== RQ4: Active Learning =====
    P("\n[RQ4] Active Learning Simulation...")
    al=al_sim(si,'K562','RPE1')
    P(f"  n={al.get('n','N/A')}, n_strong={al.get('n_strong','N/A')}")
    P(f"  transfer_overlap={fmt(al.get('transfer_overlap','N/A'))}")
    if'top_k'in al:
        for k,v in al['top_k'].items():
            P(f"  top{k}: rand={fmt(v.get('rand','N/A'))}, uq={fmt(v.get('uq','N/A'))}, epi={fmt(v.get('epi','N/A'))}, oracle={fmt(v.get('oracle','N/A'))}")
            P(f"         imp_uq={fmt(v.get('imp_uq','N/A'))}×, imp_epi={fmt(v.get('imp_epi','N/A'))}×")
    # A5/A6 derived from AL results (no separate training needed)
    a5_a6={'A5_note':'w/o UQ+AL = random baseline from AL sim','A6_note':'w/o AL = random baseline from AL sim'}
    if'top_k'in al and 10 in al['top_k']:
        a5_a6['A5_rand_top10']=al['top_k'][10]['rand']
        a5_a6['A6_rand_top10']=al['top_k'][10]['rand']
    P(f"  A5(w/o UQ+AL): random baseline only — top10={fmt(a5_a6.get('A5_rand_top10','N/A'))}")
    P(f"  A6(w/o AL):    random baseline only — top10={fmt(a5_a6.get('A6_rand_top10','N/A'))}")

    # ===== Norman Epistasis Precision =====
    P("\n[RQ1] Norman Epistasis Precision + Coverage...")
    try:
        ad_n=load_norman(); ad_n,X_n=prep(ad_n,nh=500)
        nprec=norman_precision(ad_n,X_n)
        P(f"  n_doubles={nprec.get('n_doubles','N/A')}, n_epi={nprec.get('n_epi','N/A')}")
        P(f"  trivial_prec={fmt(nprec.get('trivial_prec','N/A'))}, rho_pred_actual={fmt(nprec.get('rho_pred_actual','N/A'))}")
        P(f"  mean_coverage_90={fmt(nprec.get('mean_coverage_90','N/A'))}")
        if'precision'in nprec:
            for kf,v in nprec['precision'].items():
                P(f"  {kf}: prec={fmt(v.get('precision','N/A'))}, rand={fmt(v.get('random','N/A'))}, imp={fmt(v.get('imp','N/A'))}×")
    except Exception as ex:
        P(f"  Norman FAILED: {ex}"); nprec={'error':str(ex)}

    # ===== Summary =====
    dt=time.time()-t0
    P("\n"+"="*70+"\nSUMMARY\n"+"="*70)
    P(f"RQ2 (UQ):          rho_mc={fmt(uq.get('rho_mc','N/A'))}, coverage_90={fmt(uq.get('mean_coverage_90','N/A'))}, in_range={uq.get('coverage_in_range','N/A')}")
    P(f"RQ3 (Cross-CT):    A1 rho_add={fmt(ti.get('rho_add','N/A'))}, A7 rho={fmt(a7.get('rho','N/A'))}")
    P(f"RQ3 (ICM ablation):imp_add={fmt(abl.get('imp_add','N/A'))}× (negative if <1)")
    P(f"RQ3 (CPA B1):      rho_add={fmt(cpa_transfer.get('rho_add','N/A'))}")
    P(f"RQ3 (A3 single):   add={fmt(a3_add.get('rho','N/A'))}, mult={fmt(a3_mult.get('rho','N/A'))}, prod={fmt(a3_prod.get('rho','N/A'))}")
    P(f"RQ3 (A4 trivial):  rho={fmt(a4.get('rho','N/A'))}")
    P(f"RQ4 (AL):          top10 imp_uq={fmt(al.get('top_k',{}).get(10,{}).get('imp_uq','N/A'))}×")
    if not isinstance(nprec,dict)or'error'not in nprec:
        P(f"RQ1 (Norman):      prec_top20={fmt(nprec.get('precision',{}).get('prec_top20',{}).get('precision','N/A'))}, cov={fmt(nprec.get('mean_coverage_90','N/A'))}")
    P(f"Runtime: {dt:.1f}s")

    # ===== PASS/FAIL vs Success Criteria =====
    P("\n--- Success Criteria Check ---")
    rho_a7=a7.get('rho',0); rho_uq=uq.get('rho_mc',0)
    cov_90=uq.get('mean_coverage_90',0)
    al_imp=al.get('top_k',{}).get(10,{}).get('imp_uq',0)if'top_k'in al else 0
    P(f"  Cross-CT rho (A7):  {fmt(rho_a7)} {'PASS'if rho_a7>0.4 else'PARTIAL'if rho_a7>0.2 else'FAIL'} (target>0.4)")
    P(f"  U-Error rho (MC):   {fmt(rho_uq)} {'PASS'if rho_uq>0.6 else'FAIL'} (target>0.6)")
    P(f"  Coverage (90% CI):  {fmt(cov_90)} {'PASS'if 0.85<=cov_90<=0.95 else'FAIL'} (target 0.85-0.95)")
    P(f"  AL improvement:     {fmt(al_imp)}× {'PASS'if al_imp>2 else'PARTIAL'if al_imp>1.5 else'FAIL'} (target>2×)")
    if not isinstance(nprec,dict)or'error'not in nprec:
        np_val=nprec.get('precision',{}).get('prec_top20',{}).get('precision',0)
        P(f"  Norman Precision:   {fmt(np_val)} {'PASS'if np_val>0.6 else'FAIL'} (target>0.6)")

    # ===== Save JSON =====
    results={'run':'run_12','timestamp':time.strftime('%Y-%m-%d %H:%M:%S'),'runtime_s':dt,
             'RQ2_UQ':uq,'RQ2_UQ_noICM':uq2,
             'RQ3_A1_ICM':ti,'RQ3_A2_noICM':tn,'RQ3_ablation_A1vsA2':abl,
             'RQ3_A7_res_rank':a7,'RQ3_A3_single':{'add':a3_add,'mult':a3_mult,'prod':a3_prod},
             'RQ3_A4_trivial':a4,'RQ3_B1_CPA':cpa_transfer,
             'RQ4_AL':al,'A5_A6_derived':a5_a6,
             'RQ1_Norman':nprec if not isinstance(nprec,dict)or'error'not in nprec else str(nprec)}
    jp=os.path.join(RESULTS_DIR,'run_12_results.json')
    with open(jp,'w',encoding='utf-8')as f: json.dump(results,f,indent=2,default=str)
    P(f"\nResults saved to {jp}")

if __name__=='__main__': main()