import numpy as np
import pandas as pd
from scipy.stats import rankdata
from .config import CANDIDATE_GENES, COMPARISON_GENES

SEED=20260907

def association(x,y):
    x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float)
    if len(x)!=len(y) or len(x)<3 or not np.isfinite(x).all() or not np.isfinite(y).all():raise ValueError('Invalid association input')
    xc=x-x.mean();yc=y-y.mean();ss=xc@xc
    if ss<=0:return np.nan,np.nan
    beta=(xc@yc)/ss;residual=yc-beta*xc;se=np.sqrt((residual@residual)/(len(x)-2)/ss)
    return float(beta),float(beta/se) if se>0 else np.nan

def permute_labels(y,rng):return np.asarray(y)[rng.permutation(len(y))]

def permutation(x,y,n=5000,seed=SEED+1):
    rng=np.random.default_rng(seed);beta,t=association(x,y)
    ts=np.array([association(x,permute_labels(y,rng))[1] for _ in range(n)])
    if not np.isfinite(ts).all():raise ValueError('Invalid permutation statistics')
    return dict(status='estimated',seed=seed,iterations=n,observed_beta=beta,observed_abs_t=abs(t),empirical_p=float((1+(abs(ts)>=abs(t)).sum())/(n+1))),ts

def bootstrap(x,y,n=2000,seed=SEED+2):
    rng=np.random.default_rng(seed);x=np.asarray(x);y=np.asarray(y);betas=[]
    for _ in range(n):
        idx=rng.integers(0,len(x),len(x));betas.append(association(x[idx],y[idx])[0])
    valid=np.asarray(betas);valid=valid[np.isfinite(valid)]
    if len(valid)<.95*n:raise ValueError('Too many degenerate bootstrap samples')
    return dict(status='estimated',seed=seed,iterations=n,valid_iterations=len(valid),ci_low=float(np.quantile(valid,.025)),ci_high=float(np.quantile(valid,.975)),median_beta=float(np.median(valid))),betas

def leave_one_out(z):
    if not set(CANDIDATE_GENES).issubset(z.index):raise ValueError('Missing candidate')
    return {g:z.loc[[h for h in CANDIDATE_GENES if h!=g]].mean(axis=0) for g in CANDIDATE_GENES}

def matched_random(expression,z,patient_ids,y,observed_x,n=5000,seed=SEED):
    excluded=set(CANDIDATE_GENES)|set(COMPARISON_GENES)
    mean=expression.mean(axis=1);var=expression.var(axis=1,ddof=0)
    eligible=(np.expm1(expression*np.log(2))>=.1).mean(axis=1)>=.2
    eligible &= var>1e-8
    pool=[g for g in expression.index if eligible.loc[g] and g not in excluded and np.isfinite(z.loc[g]).all()]
    if len(pool)<100:raise ValueError('Insufficient eligible random pool')
    names=pool+list(CANDIDATE_GENES)
    coords=np.column_stack([rankdata(mean.loc[names])/len(names),rankdata(var.loc[names])/len(names)])
    neighbors=[];distances=[]
    for i in range(9):
        d=np.sqrt(((coords[:len(pool)]-coords[len(pool)+i])**2).sum(axis=1))
        neighbors.append(np.argsort(d,kind='stable')[:min(200,len(pool))]);distances.append(d)
    rng=np.random.default_rng(seed);matrix=z.loc[pool,patient_ids].to_numpy();rows=[]
    for iteration in range(n):
        chosen=[];ds=[]
        for i,options in enumerate(neighbors):
            available=[int(k) for k in options if k not in chosen]
            selected=int(rng.choice(available));chosen.append(selected);ds.append(distances[i][selected])
        genes=[pool[k] for k in chosen]
        beta,t=association(matrix[chosen].mean(axis=0),y)
        rows.append(dict(iteration=iteration+1,genes=';'.join(genes),beta=beta,abs_t=abs(t),mean_matching_distance=float(np.mean(ds)),max_matching_distance=float(max(ds))))
    frame=pd.DataFrame(rows);beta,t=association(observed_x,y)
    if not np.isfinite(frame.abs_t).all():raise ValueError('Degenerate random signature')
    summary=dict(status='estimated',seed=seed,iterations=n,pool_n=len(pool),neighbors_per_slot=min(200,len(pool)),observed_beta=beta,observed_abs_t=abs(t),
                 empirical_p=float((1+(frame.abs_t>=abs(t)).sum())/(n+1)),percentile=float(100*(frame.abs_t<abs(t)).mean()),
                 matching_mean_distance=float(frame.mean_matching_distance.mean()),matching_max_distance=float(frame.max_matching_distance.max()))
    return frame,summary
