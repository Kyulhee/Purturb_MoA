"""
Run 11: Cross-CT Epistasis Transfer Prediction + Full Ablation + UQ
====================================================================
  1. A1 vs A2: ICM vs no-ICM epistasis transfer (critical ablation)
  2. A7: Residual rank transfer, A8: UQ rank transfer
  3. Cross-CT top-k overlap, 3-formula sensitivity
  4. Norman epistasis Precision
  5. ICM transfer improvement ratio (primary RQ3 metric)
"""
import sys, io, os, json, time
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import numpy as np, torch, torch.nn as nn, torch.nn.functional as F
from sklearn.metrics import r2_score
from scipy.stats import spearmanr
import scanpy as sc, warnings
warnings.filterwarnings("ignore")
torch.manual_seed(42); np.random.seed(42)
RESULTS_DIR = os.path.dirname(os.path.abspath(__file__))
def P(*a,**k): print(*a,**k,flush=True)
def fmt(v): return f"{v:.3f}" if isinstance(v,(int,float,np.floating)) else str(v)

# ===== Data =====
def load_replogle(dp='outputs/analysis/run_04/data/gears_data',mx=200):
    from gears import PertData
    P("  Loading K562..."); pk=PertData(dp); pk.load(data_name='replogle_k562_essential')
    P("  Loading RPE1..."); pr=PertData(dp); pr.load(data_name='replogle_rpe1_essential')
    ak,ar=pk.adata,pr.adata; sh=set(ak.obs['condition'].unique())&set(ar.obs['condition'].unique()); sh.discard('ctrl')
    gn=sorted(set(ak.var_names)&set(ar.var_names)); kp=sh|{'ctrl'}
    ak=ak[ak.obs['condition'].isin(kp),gn].copy(); ar=ar[ar.obs['condition'].isin(kp),gn].copy()
    ak=_sub(ak,mx); ar=_sub(ar,mx); ak.obs['cell_type']='K562'; ar.obs['cell_type']='RPE1'
    ac=ak.concatenate(ar,batch_key='batch'); P(f"  Combined: {ac.shape}")
    return ac,sh

def load_norman(dp='outputs/analysis/run_04/data/gears_data'):
    from gears import PertData; pn=PertData(dp); pn.load(data_name='norman'); return pn.adata

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

def train_fcr(X,pi,ci,np_,nct,zd=8,icm=False,ep=150,bs=512,iw=10.):
    ng=X.shape[1]; enc=Enc(ng,np_,zd,nct); dec=Dec(zd,ng)
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

# ===== Epistasis Scores (3-Formula) =====
def compute_epi_scores(enc,dec,ad,X,pim,ctm,nct,zd):
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
            # MC Dropout (30 samples)
            dec.train(); mcp=[]
            with torch.no_grad():
                for _ in range(30): mcp.append(dec(zxr,ztr,ztxr)[0].numpy())
            dec.eval(); mc_var=float(np.var(mcp,axis=0).mean())
            pred_err=float(np.abs(y_obs-y_pred).mean())
            icm_viol=float(np.abs(d['z_tx'].numpy()).mean())
            # 3 formulas
            e_add=pred_err
            e_mult=float(np.abs(y_obs-y_pred).mean()/(np.abs(y_pred).mean()+1e-8))
            pn_=float(np.abs(y_obs*cm-y_pred*cm).mean())/(float(np.abs(y_obs*cm).mean())+1e-8)
            # Self-combo residual (A7)
            r_add=y_obs-(2*y_pred-cm); r_add_mag=float(np.abs(r_add).mean())
            scores[ct][pn]={'add':e_add,'mult':e_mult,'prod':pn_,'icm_viol':icm_viol,
                            'pred_err':pred_err,'mc_var':mc_var,'r_add_mag':r_add_mag}
    return scores,pd_

# ===== Cross-CT Transfer (RQ3 Core) =====
def cross_ct_transfer(scores,src='K562',tgt='RPE1'):
    if src not in scores or tgt not in scores: return{'error':'missing CT'}
    sh=set(scores[src].keys())&set(scores[tgt].keys())
    if len(sh)<10: return{'error':f'too few: {len(sh)}'}
    r={}
    # Per-formula rho
    for f in['add','mult','prod']:
        sv=[scores[src][p][f]for p in sh]; tv=[scores[tgt][p][f]for p in sh]
        rho,pv=spearmanr(sv,tv); r[f'rho_{f}']=float(rho); r[f'p_{f}']=float(pv)
    # Signal-specific transfer
    for sig in['icm_viol','pred_err','mc_var','r_add_mag']:
        sv=[scores[src][p][sig]for p in sh]; tv=[scores[tgt][p][sig]for p in sh]
        rho,pv=spearmanr(sv,tv); r[f'rho_{sig}']=float(rho); r[f'p_{sig}']=float(pv)
    # Composite: best weighted combination
    br=-2; bw=None
    for wm in np.arange(0,1.01,.2):
        for wi in np.arange(0,1.01-wm,.2):
            we=1-wm-wi
            sv=[wm*scores[src][p]['mc_var']+wi*scores[src][p]['icm_viol']+we*scores[src][p]['pred_err']for p in sh]
            tv=[wm*scores[tgt][p]['mc_var']+wi*scores[tgt][p]['icm_viol']+we*scores[tgt][p]['pred_err']for p in sh]
            rho,_=spearmanr(sv,tv)
            if rho>br: br=rho; bw={'mc':float(wm),'icm':float(wi),'err':float(we)}
    r['rho_composite']=float(br); r['composite_w']=bw
    # Top-k overlap
    sl=sorted(sh)
    for k in[10,20,50]:
        if k>len(sl): continue
        st=set(sorted(sl,key=lambda p:-scores[src][p]['add'])[:k])
        tt=set(sorted(sl,key=lambda p:-scores[tgt][p]['add'])[:k])
        ol=len(st&tt)/k; rnd=k/len(sl)
        r[f'top{k}_overlap']=float(ol); r[f'top{k}_random']=float(rnd); r[f'top{k}_imp']=float(ol/max(rnd,1e-8))
    r['n_shared']=len(sh)
    return r

# ===== Ablation: A1 vs A2 =====
def ablation_icm(si,sn,src='K562',tgt='RPE1'):
    ti=cross_ct_transfer(si,src,tgt); tn=cross_ct_transfer(sn,src,tgt)
    r={'A1_icm':ti,'A2_noicm':tn}
    for f in['add','mult','prod','icm_viol','pred_err','mc_var','r_add_mag','composite']:
        ri=ti.get(f'rho_{f}',0); rn=tn.get(f'rho_{f}',0)
        r[f'imp_{f}']=float(ri/rn)if abs(rn)>1e-8 else(float('inf')if abs(ri)>1e-8 else 0.)
    for k in[10,20,50]:
        r[f'top{k}_icm']=ti.get(f'top{k}_overlap',0); r[f'top{k}_noicm']=tn.get(f'top{k}_overlap',0)
    return r

# ===== A7/A8 =====
def transfer_res_rank(scores,src='K562',tgt='RPE1'):
    if src not in scores or tgt not in scores: return{'error':'missing CT'}
    sh=set(scores[src].keys())&set(scores[tgt].keys())
    if len(sh)<10: return{'error':'too few'}
    sv=[scores[src][p]['r_add_mag']for p in sh]; tv=[scores[tgt][p]['r_add_mag']for p in sh]
    rho,pv=spearmanr(sv,tv); sl=sorted(sh); tk={}
    for k in[10,20,50]:
        if k>len(sl): continue
        st=set(sorted(sl,key=lambda p:-scores[src][p]['r_add_mag'])[:k])
        tt=set(sorted(sl,key=lambda p:-scores[tgt][p]['r_add_mag'])[:k])
        tk[f'top{k}_overlap']=float(len(st&tt)/k)
    return{'rho':float(rho),'p':float(pv),'n':len(sh),**tk,'method':'A7_residual_rank'}

def transfer_uq_rank(scores,src='K562',tgt='RPE1'):
    if src not in scores or tgt not in scores: return{'error':'missing CT'}
    sh=set(scores[src].keys())&set(scores[tgt].keys())
    if len(sh)<10: return{'error':'too few'}
    sv=[scores[src][p]['mc_var']+scores[src][p]['pred_err']for p in sh]
    tv=[scores[tgt][p]['mc_var']+scores[tgt][p]['pred_err']for p in sh]
    rho,pv=spearmanr(sv,tv); sl=sorted(sh); tk={}
    for k in[10,20,50]:
        if k>len(sl): continue
        st=set(sorted(sl,key=lambda p:-(scores[src][p]['mc_var']+scores[src][p]['pred_err']))[:k])
        tt=set(sorted(sl,key=lambda p:-(scores[tgt][p]['mc_var']+scores[tgt][p]['pred_err']))[:k])
        tk[f'top{k}_overlap']=float(len(st&tt)/k)
    return{'rho':float(rho),'p':float(pv),'n':len(sh),**tk,'method':'A8_uq_rank'}

# ===== RQ2 Holdout =====
def holdout_eval(enc,dec,ad,X,pim,ctm,nct,zd,n_ho=50):
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
            dec.train(); mcp=[]
            with torch.no_grad():
                zs=czx[ct].unsqueeze(0); zsm=ztm.mean(0).unsqueeze(0); ztm2=ztxm.mean(0).unsqueeze(0)
                for _ in range(20): mcp.append(dec(zs,zsm,ztm2).numpy()[0])
            dec.eval()
            res.append({'pert':pn,'ct':ct,'error':float(np.abs(yt-ypm).mean()),
                        'r2':float(r2_score(yt,ypm)),'mc_var':float(np.var(mcp,axis=0).mean()),
                        'ztx_var':float(ztxm.var(0).mean().item())})
    return res

def uq_holdout(hr):
    if len(hr)<5: return{'error':'too few'}
    e=np.array([r['error']for r in hr]); mc=np.array([r['mc_var']for r in hr]); zv=np.array([r['ztx_var']for r in hr])
    rmc,pmc=spearmanr(mc,e); rzv,pzv=spearmanr(zv,e)
    br,bw=-1,{'mc':.5,'ztx':.5}
    for wm in np.arange(0,1.01,.1):
        wz=1-wm; u=wm*mc+wz*zv; rho,_=spearmanr(u,e)
        if rho>br: br=rho; bw={'mc':float(wm),'ztx':float(wz)}
    return{'n':len(hr),'mean_r2':float(np.mean([r['r2']for r in hr])),
           'rho_mc':float(rmc),'rho_ztx':float(rzv),'rho_combined':float(br),'best_w':bw}

# ===== AL Simulation (RQ4) =====
def al_sim(scores,src='K562',tgt='RPE1'):
    if src not in scores or tgt not in scores: return{'error':'missing CT'}
    sh=sorted(set(scores[src].keys())&set(scores[tgt].keys()))
    if len(sh)<20: return{'error':'too few'}
    tgt_epi={p:scores[tgt][p]['add']for p in sh}
    se=sorted(tgt_epi.values(),reverse=True); thr=se[int(len(se)*.3)]
    is_s={p:tgt_epi[p]>=thr for p in sh}; ns=sum(is_s.values())
    # Orderings
    np.random.seed(42); ro=list(np.random.permutation(sh))
    uqo=sorted(sh,key=lambda p:-(scores[src][p]['mc_var']+scores[src][p]['pred_err']))
    io=sorted(sh,key=lambda p:-scores[src][p]['icm_viol'])
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
            tk[k]={'rand':rc[k-1],'uq':uc[k-1],'icm':ic[k-1],'epi':ec[k-1],'oracle':oc[k-1],
                    'imp_uq':uc[k-1]/max(rc[k-1],1e-8),'imp_epi':ec[k-1]/max(rc[k-1],1e-8)}
    # Transfer overlap: top-20 source vs target
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
        res.append({'dp':dp,'r_add_mag':float(np.abs(rad).mean()),'pred_res_mag':float(np.abs(yr-yp).mean())})

    if len(res)<5: return{'error':f'too few doubles: {len(res)}'}
    em=np.array([r['r_add_mag']for r in res]); pm_=np.array([r['pred_res_mag']for r in res])
    rho,pv=spearmanr(pm_,em); thr=np.percentile(em,70); is_e=em>=thr; ne=int(is_e.sum())
    pc={}
    for kf in[.1,.2,.3]:
        k=max(1,int(len(res)*kf)); tk=np.argsort(-pm_)[:k]; p_=float(is_e[tk].mean())
        r_=float(ne/len(res)); pc[f'prec_top{int(kf*100)}']={'precision':p_,'random':r_,'imp':float(p_/max(r_,1e-8))}
    return{'n_doubles':len(res),'n_epi':ne,'trivial_prec':float(ne/len(res)),
           'rho_pred_actual':float(rho),'p_pred_actual':float(pv),'precision':pc}

# ===== Main =====
def main():
    t0=time.time()
    P("="*70+"\nRun 11: Cross-CT Epistasis Transfer + Full Ablation\n"+"="*70)

    # Phase 1: Load Replogle
    P("\n[Phase 1] Loading Replogle K562+RPE1...")
    ad,sh=load_replogle(); ad,X=prep(ad); ng=X.shape[1]
    ctn=sorted(ad.obs['cell_type'].unique()); ctm={n:i for i,n in enumerate(ctn)}; nct=len(ctn)
    ap=sorted(ad.obs['condition'].unique()); pim={n:i for i,n in enumerate(ap)}; np_=len(ap)
    P(f"  Genes: {ng}, Perts: {np_}, CTs: {nct}, Shared: {len(sh)}")
    ci=np.array([ctm[c]for c in ad.obs['cell_type'].values],dtype=np.int64)
    pi=np.array([pim[c]for c in ad.obs['condition'].values],dtype=np.int64)

    # Train A1 (ICM) and A2 (no ICM)
    P("\n[Config A2] FCR baseline (no ICM)...")
    e2,d2=train_fcr(X,pi,ci,np_,nct,icm=False)
    P("\n[Config A1] FCR + ICM...")
    e1,d1=train_fcr(X,pi,ci,np_,nct,icm=True,iw=10.)

    summary={'dataset':'Replogle K562+RPE1','n_shared_perts':len(sh),'n_genes':ng,'nct':nct}

    # ===== RQ2: Holdout UQ =====
    P("\n[RQ2] Holdout UQ (ICM model)...")
    hr=holdout_eval(e1,d1,ad,X,pim,ctm,nct,8,n_ho=50)
    uq=uq_holdout(hr)
    P(f"  n={uq.get('n','N/A')}, rho_mc={fmt(uq.get('rho_mc','N/A'))}, rho_combined={fmt(uq.get('rho_combined','N/A'))}")
    summary['rq2_holdout']=uq

    # ===== Epistasis Scores =====
    P("\n[Phase 2] Computing epistasis scores (ICM model)...")
    si,_=compute_epi_scores(e1,d1,ad,X,pim,ctm,nct,8)
    P("\n[Phase 2b] Computing epistasis scores (no-ICM model)...")
    sn,_=compute_epi_scores(e2,d2,ad,X,pim,ctm,nct,8)

    # ===== RQ3: Cross-CT Transfer =====
    P("\n[RQ3] Cross-CT Epistasis Transfer...")
    ti=cross_ct_transfer(si,'K562','RPE1')
    tn=cross_ct_transfer(sn,'K562','RPE1')
    P(f"  A1(ICM):  rho_add={fmt(ti.get('rho_add','N/A'))}, rho_prod={fmt(ti.get('rho_prod','N/A'))}, rho_composite={fmt(ti.get('rho_composite','N/A'))}")
    P(f"  A2(noICM): rho_add={fmt(tn.get('rho_add','N/A'))}, rho_prod={fmt(tn.get('rho_prod','N/A'))}, rho_composite={fmt(tn.get('rho_composite','N/A'))}")
    summary['rq3_A1_icm']=ti; summary['rq3_A2_noicm']=tn

    # ===== A1 vs A2 Ablation =====
    P("\n[Ablation] A1 vs A2 (ICM Transfer Improvement)...")
    abl=ablation_icm(si,sn,'K562','RPE1')
    P(f"  ICM improvement (add): {fmt(abl.get('imp_add','N/A'))}x")
    P(f"  ICM improvement (composite): {fmt(abl.get('imp_composite','N/A'))}x")
    for k in[10,20,50]:
        P(f"  Top-{k} overlap: ICM={fmt(abl.get(f'top{k}_icm','N/A'))}, noICM={fmt(abl.get(f'top{k}_noicm','N/A'))}")
    summary['ablation_A1_vs_A2']=abl

    # ===== A7: Residual Rank Transfer =====
    P("\n[A7] Residual Rank Transfer...")
    a7=transfer_res_rank(si,'K562','RPE1')
    P(f"  rho={fmt(a7.get('rho','N/A'))}, p={fmt(a7.get('p','N/A'))}")
    summary['A7_residual_rank']=a7

    # ===== A8: UQ Rank Transfer =====
    P("\n[A8] UQ Rank Transfer...")
    a8=transfer_uq_rank(si,'K562','RPE1')
    P(f"  rho={fmt(a8.get('rho','N/A'))}, p={fmt(a8.get('p','N/A'))}")
    summary['A8_uq_rank']=a8

    # ===== 3-Formula Sensitivity =====
    P("\n[3-Formula] Cross-CT sensitivity...")
    for f in['add','mult','prod']:
        P(f"  {f}: ICM rho={fmt(ti.get(f'rho_{f}','N/A'))}, noICM rho={fmt(tn.get(f'rho_{f}','N/A'))}")
    summary['formula_sensitivity']={
        'ICM':{f:ti.get(f'rho_{f}',0)for f in['add','mult','prod']},
        'noICM':{f:tn.get(f'rho_{f}',0)for f in['add','mult','prod']}
    }

    # ===== RQ4: Active Learning =====
    P("\n[RQ4] Active Learning...")
    alr=al_sim(si,'K562','RPE1')
    if'error'not in alr:
        for k,v in alr.get('top_k',{}).items():
            P(f"  Top-{k}: rand={v['rand']:.3f}, uq={v['uq']:.3f}, epi={v['epi']:.3f}, imp_uq={v['imp_uq']:.2f}x, imp_epi={v['imp_epi']:.2f}x")
        P(f"  Transfer overlap: {fmt(alr.get('transfer_overlap','N/A'))}")
    summary['rq4_al']=alr

    # ===== Norman Precision =====
    P("\n[Norman] Epistasis Precision...")
    try:
        ad_n=load_norman(); ad_n,X_n=prep(ad_n)
        np_=norman_precision(ad_n,X_n)
        P(f"  n_doubles={np_.get('n_doubles','N/A')}, trivial_prec={fmt(np_.get('trivial_prec','N/A'))}")
        P(f"  rho_pred_actual={fmt(np_.get('rho_pred_actual','N/A'))}")
        for k,v in np_.get('precision',{}).items():
            P(f"  {k}: precision={fmt(v.get('precision','N/A'))}, random={fmt(v.get('random','N/A'))}, imp={fmt(v.get('imp','N/A'))}x")
        summary['norman_precision']=np_
    except Exception as ex:
        P(f"  Norman failed: {ex}")
        summary['norman_precision']={'error':str(ex)}

    # ===== Summary =====
    el=time.time()-t0
    P("\n"+"="*70+"\nRUN 11 SUMMARY\n"+"="*70)
    P(f"  Dataset: Replogle K562+RPE1, {len(sh)} shared perts, {ng} genes")
    P(f"  RQ2: rho_combined={fmt(uq.get('rho_combined','N/A'))} (target>0.6)")
    P(f"  RQ3 A1(ICM):  rho_add={fmt(ti.get('rho_add','N/A'))}, rho_composite={fmt(ti.get('rho_composite','N/A'))}")
    P(f"  RQ3 A2(noICM): rho_add={fmt(tn.get('rho_add','N/A'))}, rho_composite={fmt(tn.get('rho_composite','N/A'))}")
    P(f"  ICM improvement (add): {fmt(abl.get('imp_add','N/A'))}x (target>1.5x)")
    P(f"  A7 residual rank rho: {fmt(a7.get('rho','N/A'))}")
    P(f"  A8 uq rank rho: {fmt(a8.get('rho','N/A'))}")
    if 'top_k' in alr:
        for k,v in alr['top_k'].items():
            P(f"  AL Top-{k}: imp_epi={v['imp_epi']:.2f}x (target>2x)")
    if 'norman_precision' in summary and 'error' not in summary['norman_precision']:
        np_=summary['norman_precision']
        P(f"  Norman: trivial_prec={fmt(np_.get('trivial_prec','N/A'))}, rho={fmt(np_.get('rho_pred_actual','N/A'))}")
    P(f"\n  Elapsed: {el:.0f}s")

    # Save results
    with open(os.path.join(RESULTS_DIR,'run_11_results.json'),'w') as f:
        json.dump(summary,f,indent=2,default=lambda o:float(o)if isinstance(o,np.floating)else int(o)if isinstance(o,np.integer)else str(o))
    P(f"\n  Saved to {RESULTS_DIR}/run_11_results.json")

if __name__=='__main__': main()