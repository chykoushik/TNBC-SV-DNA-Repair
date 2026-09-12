import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import statsmodels.api as sm

def render(root,data,tables,models,null_summary):
    dest=root/'outputs/v2/figures/stage2';dest.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    def save(fig,name):
        fig.tight_layout();fig.savefig(dest/(name+'.pdf'));fig.savefig(dest/(name+'.png'),dpi=300);plt.close(fig)
    def unavailable(ax,title):
        ax.set_axis_off();ax.set_title(title);ax.text(.5,.5,'Not estimable\nIndependent genomic outcome data unavailable locally',ha='center',va='center',transform=ax.transAxes)
    def scatter(ax,x,y,xlabel,ylabel,model):
        ax.scatter(x,y,s=20,alpha=.7,color='#28668a',edgecolors='none')
        if model is not None:
            grid=np.linspace(min(x),max(x),150);pred=model.get_prediction(sm.add_constant(grid)).summary_frame()
            ax.plot(grid,pred['mean'],color='#333333');ax.fill_between(grid,pred.mean_ci_lower,pred.mean_ci_upper,color='#28668a',alpha=.18)
        ax.set_xlabel(xlabel);ax.set_ylabel(ylabel)
        if ylabel=='HRD total':ax.set_ylim(bottom=0)
    fig,ax=plt.subplots(figsize=(6,4.5));d=data[['candidate_mean_z','HRD_total']].dropna()
    if models.get('primary') is not None:scatter(ax,d.candidate_mean_z,d.HRD_total,'Candidate mean z-score','HRD total',models['primary'])
    else:unavailable(ax,'Figure 1: HRD total versus candidate program')
    save(fig,'figure1_hrd_scatter')
    fig,ax=plt.subplots(figsize=(7,5));rows=[]
    for row in tables['primary']+tables['secondary']:
        if row.get('status')!='estimated' or row['outcome']=='WGD':continue
        subset=data[[row['outcome'],'candidate_mean_z']].dropna();scale=subset.candidate_mean_z.std(ddof=0)/subset[row['outcome']].std(ddof=0)
        rows.append((row['outcome'],row['beta']*scale,row['ci_low']*scale,row['ci_high']*scale))
    if rows:
        for i,(name,b,l,u) in enumerate(rows):ax.errorbar(b,i,xerr=[[b-l],[u-b]],fmt='o',color='#28668a')
        ax.set_yticks(range(len(rows)),[r[0] for r in rows]);ax.axvline(0,color='gray',lw=.8);ax.set_xlabel('Standardized effect (outcome SD per predictor SD), 95% CI')
        ax.set_title('Figure 2: genomic outcomes; WGD OR reported separately')
    else:unavailable(ax,'Figure 2: genomic outcome effects')
    save(fig,'figure2_outcome_forest')
    for key,title,name in [('single','Figure 3: individual candidate genes','figure3_single_genes'),('logo','Figure 5: leave-one-gene-out','figure5_leave_one_gene_out')]:
        fig,ax=plt.subplots(figsize=(7,5));rows=[r for r in tables[key] if r['status']=='estimated']
        if rows:
            for i,r in enumerate(rows):ax.errorbar(r['beta'],i,xerr=[[r['beta']-r['ci_low']],[r['ci_high']-r['beta']]],fmt='o',color='#28668a')
            labels=[r.get('gene',r.get('omitted_gene',r['predictor']))+(f" (q={r['fdr']:.3g})" if 'fdr' in r else '') for r in rows]
            ax.set_yticks(range(len(rows)),labels);ax.axvline(0,color='gray',lw=.8);ax.set_xlabel('HRD units per signature unit (95% CI)');ax.set_title(title)
        else:unavailable(ax,title)
        save(fig,name)
    fig,ax=plt.subplots(figsize=(6,4.5));null=tables['null']
    if len(null):
        ax.hist(null.abs_t,bins=45,color='#8fb2c7');ax.axvline(null_summary['observed_abs_t'],color='#a34a27',label='Observed candidate');ax.set_xlim(left=0)
        ax.set_xlabel('Absolute unadjusted t statistic');ax.set_ylabel('Random gene-set count');ax.legend()
    else:unavailable(ax,'Figure 4: matched random gene-set null')
    save(fig,'figure4_random_null')
    fig,ax=plt.subplots(figsize=(6,4.5));d=data[['candidate_mean_z','proliferation']].dropna()
    if models.get('proliferation') is not None:
        scatter(ax,d.candidate_mean_z,d.proliferation,'Candidate mean z-score','Proliferation mean z-score',models['proliferation'])
        adjusted=next(r for r in tables['proliferation'] if r['model']=='HRD adjusted for proliferation')
        annotation=f"HRD-adjusted candidate beta={adjusted['beta']:.3g}; P={adjusted['p']:.3g}" if adjusted['status']=='estimated' else 'HRD-adjusted association unavailable: no local HRD outcomes'
        ax.set_title('Figure 6: candidate program and proliferation');fig.text(.5,.01,annotation,ha='center',fontsize=8);fig.subplots_adjust(bottom=.2)
    else:unavailable(ax,'Figure 6: proliferation control')
    save(fig,'figure6_proliferation')
    model=models.get('primary')
    if model is not None:
        fig,axes=plt.subplots(1,3,figsize=(12,4));axes[0].scatter(model.fittedvalues,model.resid,s=15);axes[0].axhline(0,color='gray');axes[0].set(xlabel='Fitted HRD',ylabel='Residual')
        sm.qqplot(model.resid,line='45',fit=True,ax=axes[1]);axes[2].stem(np.arange(len(model.resid)),model.get_influence().cooks_distance[0]);axes[2].axhline(4/len(model.resid),color='gray');axes[2].set(xlabel='Model patient index',ylabel="Cook's distance")
        save(fig,'supplement_primary_diagnostics')
