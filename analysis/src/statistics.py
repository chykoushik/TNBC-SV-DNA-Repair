import warnings
import numpy as np
import pandas as pd
import scipy.stats as st
import statsmodels.api as sm
from statsmodels.stats.diagnostic import het_breuschpagan, linear_reset
from statsmodels.stats.outliers_influence import variance_inflation_factor

def bh(pvalues):
    p=np.asarray(pvalues,dtype=float);out=np.full(len(p),np.nan);valid=np.isfinite(p)
    if np.any((p[valid]<0)|(p[valid]>1)):raise ValueError('Invalid P value')
    order=np.argsort(p[valid]);ranked=p[valid][order];adjusted=np.minimum.accumulate((ranked*len(ranked)/np.arange(1,len(ranked)+1))[::-1])[::-1]
    subset=np.empty(len(ranked));subset[order]=np.minimum(1,adjusted);out[valid]=subset
    return out

def unavailable(outcome,predictor,label,n,reason):
    return dict(model=label,outcome=outcome,predictor=predictor,N=n,status='unavailable',reason=reason,
                beta=np.nan,se=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,r_squared=np.nan,
                hc3_se=np.nan,hc3_ci_low=np.nan,hc3_ci_high=np.nan,hc3_p=np.nan,
                spearman_rho=np.nan,spearman_p=np.nan,pearson_r=np.nan,pearson_p=np.nan)

def fit_model(data,outcome='HRD_total',predictor='candidate_mean_z',covariates=(),label='Model 0',binary=False):
    columns=list(dict.fromkeys([outcome,predictor,*covariates]))
    if not set(columns).issubset(data):return unavailable(outcome,predictor,label,0,'missing columns'),None,[]
    d=data[columns].replace([np.inf,-np.inf],np.nan).dropna()
    if len(d)<10:return unavailable(outcome,predictor,label,len(d),'fewer than 10 complete patients'),None,[]
    x=pd.get_dummies(d[[predictor,*covariates]],columns=[c for c in covariates if c in ['stage','PAM50']],drop_first=True,dtype=float)
    x=sm.add_constant(x.astype(float),has_constant='add');y=d[outcome].astype(float)
    if np.linalg.matrix_rank(x)<x.shape[1] or y.nunique()<2:return unavailable(outcome,predictor,label,len(d),'rank deficient/constant outcome'),None,[]
    if covariates and len(d)<10*x.shape[1]:return unavailable(outcome,predictor,label,len(d),'fewer than 10 complete patients per estimated parameter'),None,[]
    vif={c:float(variance_inflation_factor(x.values,i)) for i,c in enumerate(x.columns) if c!='const'}
    if binary:
        if set(y.unique())!={0.,1.} or min(y.value_counts())<10:return unavailable(outcome,predictor,label,len(d),'binary class has fewer than 10 patients'),None,[]
        scale=float(x[predictor].std(ddof=0));x[predictor]/=scale
        try:
            with warnings.catch_warnings(record=True) as caught:
                fitted=sm.Logit(y,x).fit(disp=False)
            if caught or not fitted.mle_retvals['converged']:raise ValueError('Logistic fit warning/nonconvergence')
        except Exception as exc:return unavailable(outcome,predictor,label,len(d),str(exc)),None,[]
        ci=fitted.conf_int().loc[predictor]
        result=dict(model=label,outcome=outcome,predictor=predictor,N=len(d),status='estimated',reason='',beta=float(fitted.params[predictor]),
                    se=float(fitted.bse[predictor]),ci_low=float(ci.iloc[0]),ci_high=float(ci.iloc[1]),p=float(fitted.pvalues[predictor]),
                    odds_ratio=float(np.exp(fitted.params[predictor])),or_ci_low=float(np.exp(ci.iloc[0])),or_ci_high=float(np.exp(ci.iloc[1])),
                    predictor_sd_for_OR=scale,r_squared=float(fitted.prsquared),vif=vif)
        return result,fitted,[]
    fitted=sm.OLS(y,x).fit();robust=fitted.get_robustcov_results(cov_type='HC3',use_t=True);idx=list(x.columns).index(predictor)
    ci=fitted.conf_int().loc[predictor];hc=robust.conf_int()[idx]
    pear=st.pearsonr(d[predictor],y);spear=st.spearmanr(d[predictor],y)
    result=dict(model=label,outcome=outcome,predictor=predictor,N=len(d),status='estimated',reason='',beta=float(fitted.params[predictor]),se=float(fitted.bse[predictor]),
                ci_low=float(ci.iloc[0]),ci_high=float(ci.iloc[1]),p=float(fitted.pvalues[predictor]),r_squared=float(fitted.rsquared),
                hc3_se=float(robust.bse[idx]),hc3_ci_low=float(hc[0]),hc3_ci_high=float(hc[1]),hc3_p=float(robust.pvalues[idx]),
                spearman_rho=float(spear.statistic),spearman_p=float(spear.pvalue),pearson_r=float(pear.statistic),pearson_p=float(pear.pvalue),
                vif=vif,coefficient_stability_warning=any(v>5 for v in vif.values()))
    influence=fitted.get_influence();cooks=influence.cooks_distance[0]
    result.update(breusch_pagan_p=float(het_breuschpagan(fitted.resid,x)[1]),shapiro_residual_p=float(st.shapiro(fitted.resid).pvalue),
                  reset_quadratic_p=float(linear_reset(fitted,power=2,use_f=True).pvalue),cooks_threshold=4/len(d),influential_n=int((cooks>4/len(d)).sum()))
    diagnostics=[dict(row_index=int(i),fitted=float(f),residual=float(r),cooks_distance=float(c),leverage=float(h),influential=bool(c>4/len(d))) for i,f,r,c,h in zip(d.index,fitted.fittedvalues,fitted.resid,cooks,influence.hat_matrix_diag)]
    return result,fitted,diagnostics
