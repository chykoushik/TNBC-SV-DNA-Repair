import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import statsmodels.api as sm

def make_figures(root,metabric,meta_results,geo,geo_summary,single,sc_summary,evidence):
    dest=root/'outputs/v2/figures/stage3';dest.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    def save(fig,name):
        fig.tight_layout();fig.savefig(dest/(name+'.pdf'),bbox_inches='tight');fig.savefig(dest/(name+'.png'),dpi=300,bbox_inches='tight');plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    d=metabric[['candidate_mean_z','altered_gene_fraction']].dropna();axes[0].scatter(d.candidate_mean_z,d.altered_gene_fraction,s=12,alpha=.55)
    fit=sm.OLS(d.altered_gene_fraction,sm.add_constant(d.candidate_mean_z)).fit();grid=np.linspace(d.candidate_mean_z.min(),d.candidate_mean_z.max(),150);pred=fit.get_prediction(sm.add_constant(grid)).summary_frame()
    axes[0].plot(grid,pred['mean'],color='black');axes[0].fill_between(grid,pred.mean_ci_lower,pred.mean_ci_upper,alpha=.2)
    adjusted=next(r for r in meta_results if r['endpoint']=='altered_gene_fraction_adjusted')
    axes[0].set(xlabel='Frozen candidate mean z-score',ylabel='Fraction of assayed genes altered',ylim=(-.025,1.025),title=f"METABRIC TNBC (N={len(d)})\nAge/proliferation-adjusted P={adjusted['p']:.3g}")
    for ax,endpoint in zip(axes[1:],['OS','RFS']):
        r=next(r for r in meta_results if r['endpoint']==endpoint and r['gene']=='composite')
        if r['status']=='estimated':
            hr=r['hazard_ratio'];ax.errorbar(hr,0,xerr=[[hr-r['hr_ci_low']],[r['hr_ci_high']-hr]],fmt='o',color='#28668a');ax.axvline(1,color='gray');ax.set_yticks([])
            ax.set(xlabel='Hazard ratio per signature SD (95% CI)',title=f"{endpoint}: N={r['N']}, events={r['events']}\nP={r['p']:.3g}; q={r['fdr']:.3g}")
        else:ax.set_axis_off();ax.text(.5,.5,endpoint+' unavailable',ha='center',transform=ax.transAxes)
    fig.suptitle('Partial genomic replication only: gene-level CNA is not HRD scar',y=1.04)
    save(fig,'figure1_metabric_validation')
    fig,ax=plt.subplots(figsize=(6,4))
    if 'candidate_mean_z' not in geo:
        counts=geo.pCR.value_counts();ax.bar(['Residual disease','pCR'],[counts.get(0,0),counts.get(1,0)],color=['#879aa4','#467e99'])
        ax.set(ylabel='Patients with observed response',title=f"GSE25066: {geo_summary['clinical_tnbc_n']} clinical TNBC")
        ax.text(.5,1.02,'Signature association unavailable: GPL96 annotation missing',ha='center',transform=ax.transAxes,fontsize=9)
    else:
        d=geo[['candidate_mean_z','pCR']].dropna();ax.scatter(d.candidate_mean_z,d.pCR,s=15,alpha=.4)
        ax.set(xlabel='Prespecified reduced-platform mean z-score',ylabel='Observed pCR (1) / residual disease (0)')
    save(fig,'figure2_gse25066_response')
    fig,ax=plt.subplots(figsize=(6,4.5));d=single[single.subtype=='TNBC'].pivot(index='donor_sample_id',columns='cell_class',values='candidate_mean_z')
    paired=d[['malignant_epithelial','nonmalignant_combined']].dropna()
    for donor,row in paired.iterrows():ax.plot([0,1],row.values,marker='o',alpha=.65,lw=1,label=donor)
    ax.set_xticks([0,1],['Malignant epithelial','Non-malignant combined']);ax.set_ylabel('Frozen nine-gene pseudobulk mean z-score')
    r=sc_summary['primary_result'];ax.set_title(f"TNBC paired donor/sample units N={r['N']}\nMean difference={r['effect']:.3f}; Wilcoxon P={r['p']:.3g}")
    save(fig,'figure3_single_cell_program')
    for name,title in [('figure4_depmap_volcano','DepMap genome-wide dependency volcano'),('figure5_dependency_forest','Top dependency effect sizes')]:
        fig,ax=plt.subplots(figsize=(6,4));ax.set_axis_off();ax.set_title(title)
        ax.text(.5,.5,'Unavailable: required DepMap files are not local\nNo genome-wide dependency tests performed',ha='center',va='center',transform=ax.transAxes)
        save(fig,name)
    columns=[('TCGA HRD','TCGA_HRD','HRD_total'),('TCGA FGA','TCGA_CIN','FGA'),('METABRIC CNA','METABRIC','altered_gene_fraction'),
             ('METABRIC CNA adjusted','METABRIC','altered_gene_fraction_adjusted'),('GSE25066 pCR','GSE25066','pCR'),
             ('TNBC malignant specificity','single_cell','malignant_vs_nonmalignant_combined'),('DepMap','DepMap','genome_wide_CRISPR'),('Drug response','drug','compound_wide_response')]
    genes=['composite',*__import__('src.config',fromlist=['CANDIDATE_GENES']).CANDIDATE_GENES]
    codes={'unavailable':0,'does_not_support':1,'supports':2};matrix=np.zeros((len(genes),len(columns)));labels={}
    for i,g in enumerate(genes):
        for j,(_,domain,endpoint) in enumerate(columns):
            selected=evidence[(evidence.gene==g)&(evidence.domain==domain)&(evidence.endpoint==endpoint)]
            if len(selected):
                r=selected.iloc[0];matrix[i,j]=codes[r.support]
                if r.support!='unavailable':labels[i,j]=f"q={r.fdr:.2g}" if pd.notna(r.fdr) else f"P={r.p:.2g}"
    fig,ax=plt.subplots(figsize=(13,7));ax.imshow(matrix,cmap=ListedColormap(['#ededed','#ddb788','#6caaa4']),vmin=0,vmax=2,aspect='auto')
    ax.set_xticks(range(len(columns)),[c[0] for c in columns],rotation=35,ha='right');ax.set_yticks(range(len(genes)),genes)
    for (i,j),label in labels.items():ax.text(j,i,label,ha='center',va='center',fontsize=8)
    ax.set_title('Integrated evidence: categorical results, not an evidence score')
    ax.legend(handles=[Patch(color=c,label=l) for c,l in [('#6caaa4','Supports at threshold'),('#ddb788','Does not support at threshold'),('#ededed','Unavailable / not tested')]],loc='upper center',bbox_to_anchor=(.5,-.35),ncol=3)
    save(fig,'figure6_integrated_evidence')
