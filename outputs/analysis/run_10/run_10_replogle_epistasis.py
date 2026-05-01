"""
Run 10: Replogle Multi-Cell-Type Epistasis + UQ Pipeline
=========================================================
Key differences from run_09c (Norman):
  - Replogle has single-gene perturbations only (no double-KO)
  - RQ2: Holdout evaluation (non-circular!)
  - Cross-cell-type epistasis consistency
  - RQ1/RQ3: Multi-cell-type ICM validation
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

def P(*a, **k): print(*a, **k, flush=True)
def fmt(v): return f"{v:.3f}" if isinstance(v,(int,float,np.floating)) else str(v)

# ---- Data Loading (from run_07) ----
def load_and_combine(data_path='outputs/analysis/run_04/data/gears_data', max_cells=200):
    from gears import PertData
    P("  Loading K562..."); pd_k = PertData(data_path); pd_k.load(data_name='replogle_k562_essential')
    P("  Loading RPE1..."); pd_r = PertData(data_path); pd_r.load(data_name='replogle_rpe1_essential')
    ak, ar = pd_k.adata, pd_r.adata
    shared = set(ak.obs['condition'].unique()) & set(ar.obs['condition'].unique()); shared.discard('ctrl')
    genes = sorted(set(ak.var_names) & set(ar.var_names))
    keep = shared | {'ctrl'}
    ak = ak[ak.obs['condition'].isin(keep), genes].copy()
    ar = ar[ar.obs['condition'].isin(keep), genes].copy()
    ak = _subsample(ak, max_cells); ar = _subsample(ar, max_cells)
    ak.obs['cell_type']='K562'; ar.obs['cell_type']='RPE1'
    ac = ak.concatenate(ar, batch_key='batch')
    P(f"  Combined: {ac.shape}, CTs: {ac.obs['cell_type'].unique()}")
    return ac, shared

def _subsample(ad, mx):
    idx=[]
    for c in ad.obs['condition'].unique():
        m=np.where(ad.obs['condition']==c)[0]
        idx.extend(np.random.choice(m,min(len(m),mx),replace=False) if len(m)>mx else m)
    return ad[sorted(idx)].copy()

def preprocess(ad, n_hvg=500):
    sc.pp.filter_genes(ad, min_cells=50); sc.pp.normalize_total(ad,1e4); sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, n_top_genes=n_hvg, flavor='seurat')
    ad = ad[:,ad.var.highly_variable].copy()
    X = ad.X.toarray().astype(np.float32) if hasattr(ad.X,'toarray') else ad.X.astype(np.float32)
    P(f"  After HVG: {ad.shape}"); return ad, X

# ---- Model ----
class Enc(nn.Module):
    def __init__(s, ng, np_, zd, nct):
        super().__init__(); s.zd=zd
        s.xe=nn.Sequential(nn.Linear(ng,256),nn.ReLU(),nn.Dropout(.1),nn.Linear(256,128),nn.ReLU())
        s.pe=nn.Embedding(np_,zd); s.cte=nn.Embedding(nct,zd)
        s.zx=nn.Linear(128+nct,zd*2); s.zt=nn.Linear(128+zd,zd*2); s.ztx=nn.Linear(128+zd,zd*2)
    def forward(s,x,pid,ctoh):
        h=s.xe(x); zt_in=s.pe(pid)
        a=s.zx(torch.cat([h,ctoh],-1)); b=s.zt(torch.cat([h,zt_in],-1)); c=s.ztx(torch.cat([h,zt_in],-1))
        return (a[:,:s.zd],a[:,s.zd:]),(b[:,:s.zd],b[:,s.zd:]),(c[:,:s.zd],c[:,s.zd:])

class Dec(nn.Module):
    def __init__(s,zd,ng,dr=.1):
        super().__init__()
        s.d=nn.Sequential(nn.Linear(3*zd,256),nn.ReLU(),nn.Dropout(dr),nn.Linear(256,128),nn.ReLU(),nn.Dropout(dr),nn.Linear(128,ng))
    def forward(s,zx,zt,ztx): return s.d(torch.cat([zx,zt,ztx],-1))

def reparam(m,lv): return m+torch.exp(.5*lv)*torch.randn_like(m)

def icm_reg(zm,cts):
    ut=torch.unique(cts)
    if len(ut)<2: return ((zm.var(0)-1)**2).mean()
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

def train_fcr(X,pi,ci,np_,nct,zd=8,icm=False,epochs=150,bs=512,icm_w=10.):
    ng=X.shape[1]; enc=Enc(ng,np_,zd,nct); dec=Dec(zd,ng)
    opt=torch.optim.Adam(list(enc.parameters())+list(dec.parameters()),lr=1e-3)
    xt=torch.FloatTensor(X); pt=torch.LongTensor(pi); ct=torch.LongTensor(ci)
    coh=F.one_hot(ct,nct).float()
    dl=torch.utils.data.DataLoader(torch.utils.data.TensorDataset(xt,pt,ct,coh),batch_size=bs,shuffle=True,drop_last=True)
    for ep in range(epochs):
        el=0
        for bx,bp,bc,bch in dl:
            opt.zero_grad()
            (zxm,zxlv),(ztm,ztlv),(ztxm,ztxlv)=enc(bx,bp,bch)
            xr=dec(reparam(zxm,zxlv),reparam(ztm,ztlv),reparam(ztxm,ztxlv))
            rl=F.mse_loss(xr,bx,reduction='sum')
            kl=sum(-.5*torch.sum(1+lv-m.pow(2)-lv.exp()) for m,lv in[(zxm,zxlv),(ztm,ztlv),(ztxm,ztxlv)])
            loss=rl+.5*kl
            if icm and nct>1: loss=loss+icm_w*icm_reg(ztxm,bc)
            loss.backward(); opt.step(); el+=loss.item()
        if(ep+1)%30==0: P(f"    Ep {ep+1}/{epochs}: loss={el/len(dl):.1f}")
    return enc,dec

# ---- Phase 2: Holdout RQ2 ----
def holdout_eval(enc,dec,ad,X,pim,ctm,nct,zd,n_ho=50):
    enc.eval(); dec.eval()
    conds=ad.obs['condition'].values; cts=ad.obs['cell_type'].values
    ctn=sorted(ad.obs['cell_type'].unique())
    sps=[c for c in np.unique(conds) if c!='ctrl' and ('+' not in c or c.endswith('+ctrl') or c.startswith('ctrl+'))]
    np.random.shuffle(sps); ho=sps[:n_ho]
    P(f"  Holdout: {len(ho)} perts")
    ctrl_zx={}
    for ct in ctn:
        ci=ctm[ct]; m=(conds=='ctrl')&(cts==ct)
        if m.sum()<5: continue
        nu=min(50,m.sum()); xc=torch.FloatTensor(X[m][:nu])
        pt=torch.full((nu,),pim['ctrl'],dtype=torch.long)
        ch=F.one_hot(torch.full((nu,),ci,dtype=torch.long),nct).float()
        with torch.no_grad(): (zxm,_),_,_=enc(xc,pt,ch)
        ctrl_zx[ct]=zxm.mean(0)
    res=[]
    for pn in ho:
        if pn not in pim: continue
        for ct in ctn:
            ci=ctm[ct]; m=(conds==pn)&(cts==ct)
            if m.sum()<10 or ct not in ctrl_zx: continue
            yt=X[m].mean(0); nu=min(30,m.sum())
            xp=torch.FloatTensor(X[m][:nu]); pt=torch.full((nu,),pim[pn],dtype=torch.long)
            ch=F.one_hot(torch.full((nu,),ci,dtype=torch.long),nct).float()
            with torch.no_grad(): (zxm,_),(ztm,_),(ztxm,_)=enc(xp,pt,ch)
            zxc=ctrl_zx[ct].unsqueeze(0).expand(nu,-1)
            with torch.no_grad(): yp=dec(zxc,ztm,ztxm)
            ypm=yp.mean(0).numpy()
            dec.train(); mcp=[]
            with torch.no_grad():
                zs=ctrl_zx[ct].unsqueeze(0); zsm=ztm.mean(0).unsqueeze(0); ztm2=ztxm.mean(0).unsqueeze(0)
                for _ in range(20): mcp.append(dec(zs,zsm,ztm2).numpy()[0])
            dec.eval(); mcv=float(np.var(mcp,axis=0).mean())
            ztv=float(ztxm.var(0).mean().item())
            err=float(np.abs(yt-ypm).mean()); r2=float(r2_score(yt,ypm))
            res.append({'pert':pn,'ct':ct,'error':err,'r2':r2,'mc_var':mcv,'ztx_var':ztv})
    return res

# ---- Phase 3: Compose Pairs ----
def compose_pairs(enc,dec,ad,X,pim,ctm,nct,zd,npairs=200):
    enc.eval(); dec.eval()
    conds=ad.obs['condition'].values; cts=ad.obs['cell_type'].values
    ctn=sorted(ad.obs['cell_type'].unique())
    sps=sorted([c for c in np.unique(conds) if c!='ctrl' and ('+' not in c or c.endswith('+ctrl') or c.startswith('ctrl+'))])
    P(f"  Single perts: {len(sps)}")
    # Encode
    pzt,pxp,ptp,pme={},{},{},{}
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
            if pn not in pzt: pzt[pn]={}; pxp[pn]={}; ptp[pn]={}; pme[pn]={}
            pzt[pn][ct]=ztxm.mean(0); pxp[pn][ct]=zxm.mean(0); ptp[pn][ct]=ztm.mean(0)
            pme[pn][ct]=X[m].mean(0)
    vp=[p for p in sps if p in pzt and len(pzt[p])==len(ctn)]
    P(f"  Valid for pairing: {len(vp)}")
    # Generate pairs
    np.random.shuffle(vp); pairs=[]
    for i in range(min(len(vp),25)):
        for j in range(i+1,min(len(vp),25)): pairs.append((vp[i],vp[j]))
    np.random.shuffle(pairs); pairs=pairs[:npairs]
    P(f"  Pairs: {len(pairs)}")
    # z_tx stats for OOD
    zts={}
    for pn in vp:
        az=np.array([pzt[pn][ct].numpy() for ct in ctn if ct in pzt[pn]])
        if len(az)>0: zts[pn]={'mean':az.mean(0),'std':az.std(0)+1e-6}
    res=[]
    for p1,p2 in pairs:
        for ct in ctn:
            if ct not in pzt.get(p1,{}) or ct not in pzt.get(p2,{}): continue
            cm=(conds=='ctrl')&(cts==ct); cmn=X[cm].mean(0)
            ztc=pzt[p1][ct]+pzt[p2][ct]
            zxr=pxp[p1][ct].unsqueeze(0); ztr=ptp[p1][ct].unsqueeze(0)
            with torch.no_grad(): yp=dec(zxr,ztr,ztc.unsqueeze(0))[0].numpy()
            ya=pme[p1][ct]; yb=pme[p2][ct]; yadd=ya+yb-cmn
            r=yp-yadd
            # OOD
            oods=[float((np.abs(ztc.numpy()-s['mean'])/s['std']).mean()) for s in zts.values()]
            ood=float(np.mean(oods))
            # MC
            dec.train(); mcp=[]
            with torch.no_grad():
                for _ in range(15): mcp.append(dec(zxr,ztr,ztc.unsqueeze(0)).numpy()[0])
            dec.eval(); mcv=float(np.var(mcp,axis=0).mean())
            r2v=r2_score(yp,yadd) if not np.isnan(r2_score(yp,yadd)) else -999.
            res.append({'p1':p1,'p2':p2,'ct':ct,'y_pred':yp,'y_add':yadd,'res':r,
                        'res_mag':float(np.abs(r).mean()),'ood':ood,'mc_var':mcv,'r2_add':float(r2v)})
    return res

# ---- Phase 4: UQ ----
def uq_holdout(hr):
    if len(hr)<5: return {'error':'too few'}
    e=np.array([r['error'] for r in hr])
    mc=np.array([r['mc_var'] for r in hr]); zv=np.array([r['ztx_var'] for r in hr])
    rmc,pmc=spearmanr(mc,e); rzv,pzv=spearmanr(zv,e)
    br,bw=-1,{'mc':.5,'ztx':.5}
    for wm in np.arange(0,1.01,.1):
        wz=1-wm; u=wm*mc+wz*zv; rho,_=spearmanr(u,e)
        if rho>br: br=rho; bw={'mc':float(wm),'ztx':float(wz)}
    return {'n':len(hr),'mean_r2':float(np.mean([r['r2'] for r in hr])),
            'rho_mc':float(rmc),'p_mc':float(pmc),'rho_ztx':float(rzv),
            'rho_combined':float(br),'best_w':bw}

def uq_composed(cr):
    if len(cr)<5: return {'error':'too few'}
    e=np.array([r['res_mag'] for r in cr])
    od=np.array([r['ood'] for r in cr]); mc=np.array([r['mc_var'] for r in cr])
    rod,pod=spearmanr(od,e); rmc,pmc=spearmanr(mc,e)
    br,bw=-1,{'ood':.5,'mc':.5}
    for wo in np.arange(0,1.01,.1):
        wm=1-wo; u=wo*od+wm*mc; rho,_=spearmanr(u,e)
        if rho>br: br=rho; bw={'ood':float(wo),'mc':float(wm)}
    return {'n':len(cr),'mean_res':float(e.mean()),'rho_ood':float(rod),'p_ood':float(pod),
            'rho_mc':float(rmc),'rho_combined':float(br),'best_w':bw}

# ---- Phase 5: Cross-CT Epistasis Consistency ----
def cross_ct_epi(cr):
    pd={}
    for r in cr:
        k=(r['p1'],r['p2'])
        if k not in pd: pd[k]={}
        pd[k][r['ct']]=r
    ctn=sorted(set(r['ct'] for r in cr))
    if len(ctn)<2: return {'error':'need 2 CTs'}
    sp=[k for k,v in pd.items() if len(v)==len(ctn)]
    P(f"  Shared pairs: {len(sp)}")
    if len(sp)<5: return {'error':'too few'}
    rv={ct:[] for ct in ctn}
    for pk in sp:
        for ct in ctn: rv[ct].append(pd[pk][ct]['res_mag'])
    rho,p=spearmanr(rv[ctn[0]],rv[ctn[1]])
    da=[]
    for pk in sp:
        r0,r1=pd[pk][ctn[0]]['res'],pd[pk][ctn[1]]['res']
        da.append(float((np.sign(r0)==np.sign(r1)).mean()))
    return {'n_shared':len(sp),'cross_ct_rho':float(rho),'cross_ct_p':float(p),
            'mean_dir_agree':float(np.mean(da))}

# ---- Phase 6: AL ----
def al_sim(cr):
    ctn=sorted(set(r['ct'] for r in cr)); sct=ctn[0]
    sr=[r for r in cr if r['ct']==sct]
    if len(sr)<10: return {'error':'too few'}
    n=len(sr); gt={r['p1']+'+'+r['p2']:r['res_mag'] for r in sr}
    ap=list(gt.keys())
    effs=sorted(gt.values(),reverse=True); thr=effs[int(n*.3)]
    istr={c:gt[c]>=thr for c in ap}; ns=sum(istr.values())
    np.random.seed(42); ro=list(np.random.permutation(ap))
    ol={r['p1']+'+'+r['p2']:r['ood'] for r in sr}
    oo=sorted(ap,key=lambda c:-ol.get(c,0))
    ml={r['p1']+'+'+r['p2']:r['mc_var'] for r in sr}
    mo=sorted(ap,key=lambda c:-ml.get(c,0))
    eo=sorted(ap,key=lambda c:-gt[c])
    def cumrec(o): f=0;c=[];[c.append((f:=f+istr.get(x,False))/max(ns,1)) for x in o];return c
    rc,oc,mc2,ec=cumrec(ro),cumrec(oo),cumrec(mo),cumrec(eo)
    tk={}
    for k in [5,10,20]:
        if k<=n: tk[k]={'rand':rc[k-1],'ood':oc[k-1],'mc':mc2[k-1],'epi':ec[k-1],
                         'imp_ood':oc[k-1]/max(rc[k-1],1e-8),'imp_epi':ec[k-1]/max(rc[k-1],1e-8)}
    tr=None
    if len(ctn)>1:
        tct=ctn[1]; tr2={r['p1']+'+'+r['p2']:r['res_mag'] for r in cr if r['ct']==tct}
        if len(tr2)>5:
            s20=set(eo[:min(20,len(eo))])
            te=sorted(tr2.values(),reverse=True); tt=te[int(len(te)*.3)]
            ts={k for k,v in tr2.items() if v>=tt}
            tr={'src':sct,'tgt':tct,'top20_overlap':float(len(s20&ts)/max(len(s20),1))}
    return {'n':n,'n_strong':ns,'sct':sct,'top_k':tk,'transfer':tr}

# ---- Phase 7: RQ1 + RQ3 ----
def eval_rq1(enc,ad,X,pim,ctm,nct,zd):
    enc.eval(); ctn=sorted(ad.obs['cell_type'].unique())
    conds=ad.obs['condition'].values; cts=ad.obs['cell_type'].values; r=[]
    for c in np.unique(conds):
        if c=='ctrl' or c not in pim: continue
        zp=[]
        for ct in ctn:
            ci=ctm[ct]; m=(conds==c)&(cts==ct)
            if m.sum()<10: zp.append(None); continue
            nu=min(100,m.sum()); xt=torch.FloatTensor(X[m][:nu])
            ch=F.one_hot(torch.full((nu,),ci,dtype=torch.long),nct).float()
            pt=torch.full((nu,),pim[c],dtype=torch.long)
            with torch.no_grad(): _,_,(ztxm,_)=enc(xt,pt,ch)
            zp.append(ztxm.mean(0).numpy())
        v=[z for z in zp if z is not None]
        if len(v)==2:
            co=np.corrcoef(v[0],v[1])[0,1]
            cs=np.dot(v[0],v[1])/(np.linalg.norm(v[0])*np.linalg.norm(v[1])+1e-8)
            r.append({'pert':c,'corr':float(co),'cos':float(cs)})
    if r: return {'mean_corr':float(np.mean([x['corr'] for x in r])),
                  'mean_cos':float(np.mean([x['cos'] for x in r])),'n':len(r)}
    return None

def eval_rq3(enc,dec,ad,X,pim,ctm,nct,zd,sct='K562',tct='RPE1'):
    enc.eval(); dec.eval()
    conds=ad.obs['condition'].values; cts=ad.obs['cell_type'].values
    ctn=sorted(ad.obs['cell_type'].unique()); si=ctm[sct]; ti=ctm[tct]
    r=[]
    for c in np.unique(conds):
        if c=='ctrl' or c not in pim: continue
        ms=(conds==c)&(cts==sct); mt=(conds==c)&(cts==tct)
        if ms.sum()<10 or mt.sum()<10: continue
        ns=min(100,ms.sum()); xs=torch.FloatTensor(X[ms][:ns])
        cs=F.one_hot(torch.full((ns,),si,dtype=torch.long),nct).float()
        ps=torch.full((ns,),pim[c],dtype=torch.long)
        with torch.no_grad(): (zxs,_),(zts,_),(ztxs,_)=enc(xs,ps,cs)
        nt=min(100,mt.sum()); xt=torch.FloatTensor(X[mt][:nt])
        ct2=F.one_hot(torch.full((nt,),ti,dtype=torch.long),nct).float()
        pt2=torch.full((nt,),pim[c],dtype=torch.long)
        with torch.no_grad(): (zxt,_),(ztt,_),(ztxt,_)=enc(xt,pt2,ct2)
        with torch.no_grad():
            xp=dec(zxt.mean(0).unsqueeze(0),ztt.mean(0).unsqueeze(0),ztxs.mean(0).unsqueeze(0))
            xo=dec(zxt.mean(0).unsqueeze(0),ztt.mean(0).unsqueeze(0),ztxt.mean(0).unsqueeze(0))
        act=X[mt].mean(0); pp=xp[0].numpy(); oo=xo[0].numpy()
        r.append({'r2_t':float(r2_score(act,pp)),'corr_t':float(np.corrcoef(act,pp)[0,1]),
                  'r2_o':float(r2_score(act,oo)),'corr_o':float(np.corrcoef(act,oo)[0,1]),
                  'ztx_cc':float(np.corrcoef(ztxs.mean(0).numpy(),ztxt.mean(0).numpy())[0,1])})
    if r: return {'mean_r2_t':float(np.mean([x['r2_t'] for x in r])),
                  'mean_corr_t':float(np.mean([x['corr_t'] for x in r])),
                  'mean_r2_o':float(np.mean([x['r2_o'] for x in r])),
                  'mean_corr_o':float(np.mean([x['corr_o'] for x in r])),
                  'mean_ztx_cc':float(np.mean([x['ztx_cc'] for x in r])),'n':len(r)}
    return None

# ---- Main ----
def main():
    t0=time.time()
    P("="*70+"\nRun 10: Replogle Multi-CT Epistasis+UQ\n"+"="*70)
    P("\n[Phase 1] Loading Replogle K562+RPE1...")
    ad,sp=load_and_combine(); ad,X=preprocess(ad)
    ng=X.shape[1]
    ctn=sorted(ad.obs['cell_type'].unique()); ctm={n:i for i,n in enumerate(ctn)}; nct=len(ctn)
    ap=sorted(ad.obs['condition'].unique()); pim={n:i for i,n in enumerate(ap)}; np_=len(ap)
    P(f"  Genes: {ng}, Perts: {np_}, CTs: {nct}")
    ci=np.array([ctm[c] for c in ad.obs['cell_type'].values],dtype=np.int64)
    pi=np.array([pim[c] for c in ad.obs['condition'].values],dtype=np.int64)

    # Train both configs
    P("\n[Config 1] FCR baseline (no ICM)...")
    e1,d1=train_fcr(X,pi,ci,np_,nct,icm=False)
    P("\n[Config 2] FCR + ICM...")
    e2,d2=train_fcr(X,pi,ci,np_,nct,icm=True,icm_w=10.)

    summary={'dataset':'Replogle K562+RPE1','n_shared_perts':len(sp),'n_genes':ng,'nct':nct}

    # RQ1
    P("\n[RQ1] z_tx Invariance...")
    rq1_1=eval_rq1(e1,ad,X,pim,ctm,nct,8)
    rq1_2=eval_rq1(e2,ad,X,pim,ctm,nct,8)
    if rq1_1: P(f"  Baseline: corr={rq1_1['mean_corr']:.4f}, cos={rq1_1['mean_cos']:.4f}")
    if rq1_2: P(f"  ICM:      corr={rq1_2['mean_corr']:.4f}, cos={rq1_2['mean_cos']:.4f}")
    summary['rq1_base_corr']=rq1_1['mean_corr'] if rq1_1 else None
    summary['rq1_icm_corr']=rq1_2['mean_corr'] if rq1_2 else None

    # RQ3
    P("\n[RQ3] Zero-Shot Transfer K562->RPE1...")
    rq3_1=eval_rq3(e1,d1,ad,X,pim,ctm,nct,8)
    rq3_2=eval_rq3(e2,d2,ad,X,pim,ctm,nct,8)
    if rq3_1: P(f"  Baseline: R2={rq3_1['mean_r2_t']:.4f}, corr={rq3_1['mean_corr_t']:.4f}, ztx_cc={rq3_1['mean_ztx_cc']:.4f}")
    if rq3_2: P(f"  ICM:      R2={rq3_2['mean_r2_t']:.4f}, corr={rq3_2['mean_corr_t']:.4f}, ztx_cc={rq3_2['mean_ztx_cc']:.4f}")
    summary['rq3_base_r2']=rq3_1['mean_r2_t'] if rq3_1 else None
    summary['rq3_icm_r2']=rq3_2['mean_r2_t'] if rq3_2 else None

    # Phase 2: Holdout RQ2 (ICM model)
    P("\n[RQ2-Holdout] Non-circular UQ (ICM model)...")
    hr=holdout_eval(e2,d2,ad,X,pim,ctm,nct,8,n_ho=50)
    P(f"  Holdout results: {len(hr)}")
    uq_h=uq_holdout(hr)
    P(f"  Mean R2: {uq_h.get('mean_r2','N/A')}")
    P(f"  rho_mc: {fmt(uq_h.get('rho_mc','N/A'))}, rho_ztx: {fmt(uq_h.get('rho_ztx','N/A'))}")
    P(f"  rho_combined: {fmt(uq_h.get('rho_combined','N/A'))}")
    summary['rq2_holdout']=uq_h

    # Phase 3+4: Compose pairs + UQ
    P("\n[Phase 3+4] Compose pairs + UQ (ICM model)...")
    cr=compose_pairs(e2,d2,ad,X,pim,ctm,nct,8,npairs=200)
    P(f"  Composed: {len(cr)}")
    if cr:
        rm=np.mean([r['res_mag'] for r in cr])
        P(f"  Mean residual mag: {rm:.4f}")
    uq_c=uq_composed(cr)
    P(f"  rho_ood: {fmt(uq_c.get('rho_ood','N/A'))}, rho_mc: {fmt(uq_c.get('rho_mc','N/A'))}")
    P(f"  rho_combined: {fmt(uq_c.get('rho_combined','N/A'))}")
    summary['rq2_composed']=uq_c

    # Phase 5: Cross-CT epistasis
    P("\n[Phase 5] Cross-CT Epistasis Consistency...")
    cce=cross_ct_epi(cr)
    P(f"  Shared pairs: {cce.get('n_shared','N/A')}")
    P(f"  Cross-CT rho: {fmt(cce.get('cross_ct_rho','N/A'))}")
    P(f"  Direction agreement: {fmt(cce.get('mean_dir_agree','N/A'))}")
    summary['cross_ct_epi']=cce

    # Phase 6: AL
    P("\n[Phase 6] Active Learning...")
    alr=al_sim(cr)
    if 'error' not in alr:
        for k,v in alr.get('top_k',{}).items():
            P(f"  Top-{k}: rand={v['rand']:.3f}, ood={v['ood']:.3f}, epi={v['epi']:.3f}, imp_ood={v['imp_ood']:.2f}x, imp_epi={v['imp_epi']:.2f}x")
        if alr.get('transfer'):
            P(f"  Transfer: top20_overlap={alr['transfer']['top20_overlap']:.3f}")
    summary['al']=alr

    # Summary
    el=time.time()-t0
    P("\n"+"="*70+"\nRUN 10 SUMMARY\n"+"="*70)
    P(f"  Dataset: Replogle K562+RPE1, {len(sp)} shared perts, {ng} genes")
    if rq1_1 and rq1_2: P(f"  RQ1: baseline corr={rq1_1['mean_corr']:.4f} -> ICM={rq1_2['mean_corr']:.4f}")
    if rq3_1 and rq3_2: P(f"  RQ3: baseline R2={rq3_1['mean_r2_t']:.4f} -> ICM={rq3_2['mean_r2_t']:.4f}")
    P(f"  RQ2 (holdout): rho_combined={fmt(uq_h.get('rho_combined','N/A'))} (target>0.6)")
    P(f"  RQ2 (composed): rho_ood={fmt(uq_c.get('rho_ood','N/A'))}")
    P(f"  Cross-CT epi rho: {fmt(cce.get('cross_ct_rho','N/A'))}")
    P(f"  Cross-CT dir agree: {fmt(cce.get('mean_dir_agree','N/A'))}")
    P(f"\n  Elapsed: {el:.0f}s")

    with open(os.path.join(RESULTS_DIR,'run_10_results.json'),'w') as f:
        json.dump(summary,f,indent=2,default=lambda o:float(o) if isinstance(o,np.floating) else int(o) if isinstance(o,np.integer) else str(o))
    P(f"\n  Saved to {RESULTS_DIR}/run_10_results.json")

if __name__=='__main__': main()
