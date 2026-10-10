from pathlib import Path
import json, shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle
from PIL import Image
import fitz
R=Path(__file__).resolve().parent;E=R.parent/'experiments';OLD=R.parents[1]/'figures'
C=dict(blue='#146B9D',teal='#007F83',purple='#7851A9',gold='#B87916',red='#BC4057',ink='#183047',muted='#556B7D',line='#D4E0E8')
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7.5,'ytick.labelsize':7.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'legend.frameon':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','lines.linewidth':1.4})

def export(fig,name):
    fig.savefig(R/f'{name}.pdf',facecolor='white');fig.savefig(R/f'{name}.svg',facecolor='white')
    fig.savefig(R/f'{name}.png',dpi=600,facecolor='white')
    im=Image.open(R/f'{name}.png').convert('RGB');im.save(R/f'{name}.tiff',compression='tiff_lzw',dpi=(600,600));im.save(R/f'{name}.jpeg',quality=100,subsampling=0,dpi=(600,600));plt.close(fig)

def schematic():
    fig=plt.figure(figsize=(170/25.4,155/25.4));ax=fig.add_axes([.02,.02,.96,.96]);ax.set(xlim=(0,100),ylim=(0,100));ax.axis('off')
    def txt(x,y,s,size=8,color=None,weight='normal',**kw):ax.text(x,y,s,fontsize=size,color=color or C['ink'],weight=weight,va='top',**kw)
    def rect(x,y,w,h,color,fill=None):ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.35,rounding_size=1.1',edgecolor=color,facecolor=fill or 'white',linewidth=.7))
    stages=[('TARGET','Who, when, care?','blue'),('OBSERVATION','What is recorded?','teal'),('IDENTIFICATION','What is compatible?','purple'),('INFERENCE','What is covered?','gold'),('DECISION','Which loss matters?','red')]
    for i,(a,b,c) in enumerate(stages):
        x=.4+20*i;rect(x,88,18.7,10,C[c],C[c]);txt(x+1.3,95.5,a,7.5,'white','bold');txt(x+1.3,92,b,6.5,'white')
        if i<4:ax.annotate('',xy=(x+20,93),xytext=(x+19,93),arrowprops=dict(arrowstyle='->',lw=.6,color=C['muted']))
    rect(.4,45,48.2,39.5,C['line'],'#F3F8FB');rect(51.3,45,48.2,39.5,C['line'],'#F8F5FB')
    txt(2.5,82,'a  Synthetic study population',10,C['blue'],'bold')
    for x,color,label in [(4,C['blue'],'Source'),(29,C['teal'],'Target')]:
        ax.add_patch(Rectangle((x,65),11,9,facecolor='white',edgecolor=color,lw=1))
        txt(x+3.8,73,'+',14,color,'bold')
        for k in range(3):ax.add_patch(Rectangle((x+1.4+3*k,66.2),1.6,1.4,facecolor=color))
        txt(x,63.4,label,8,color,'bold')
    ax.annotate('',xy=(28,69.5),xytext=(17,69.5),arrowprops=dict(arrowstyle='->',lw=1.3,color=C['muted']))
    txt(3,59.5,'Paired reference\nand routine proxy',7.5);txt(28,59.5,'Routine proxy\n+ optional audit',7.5)
    txt(3,52,'Two locked rules; fixed care; t = 0.20',7.6,weight='bold')
    txt(3,48.5,'G+  rule 1 only: 2/3      G−  rule 2 only: 1/3',7.1)
    txt(53.5,82,'b  One observed law, three worlds',10,C['purple'],'bold')
    txt(54,76.5,'Same proxy rates: 0.31 and 0.40',8)
    xmin,xmax=58,95
    xpos=lambda v:xmin+(v+.02)/.08*(xmax-xmin)
    ax.plot([xpos(0)]*2,[57,72],color=C['muted'],ls='--',lw=.7)
    for y,w,v,m,c in [(70,'B',-1/72,'v','red'),(65,'C',0,'s','muted'),(60,'A',1/18,'o','teal')]:
        txt(54,y+1.2,w,8,weight='bold');ax.plot(xpos(v),y,m,color=C[c],ms=5);txt(98,y+3.5,f'{v:+.4f}',7.4,ha='right')
    ax.plot([xpos(-1/72),xpos(1/18)],[54.5]*2,color=C['purple'],lw=3)
    txt(54,51,'Sharp population set [−1/72, 1/18]',7.7,weight='bold');txt(54,47.5,'Identification, without a confidence level',7)
    rect(.4,13,48.2,28.3,C['line'],'#F3F9F8');rect(51.3,13,48.2,28.3,C['line'],'#FCF8F1')
    txt(2.5,39,'c  Test the observation assumptions',9.5,C['teal'],'bold')
    txt(3,34,'Shared vs group-specific label channels',7.8)
    txt(3,29.5,'Random vs proxy-selected reference audit',7.8)
    txt(3,25,'Accurate vs fallible reference recording',7.8)
    txt(3,19.5,'Separate inversion → joint constraints',8,C['teal'],'bold')
    txt(3,15.7,'12,000 reused evidence records',7.5)
    txt(53.5,39,'d  Compare what new data can resolve',9.3,C['gold'],'bold')
    txt(54,34,'Random, accurate audits in both groups',7.7)
    txt(54,29.5,'Equal new budgets: 200 to 12,800',7.8)
    txt(54,25,'G+ allocation: 1/3, 1/2 or 2/3',7.8)
    txt(54,19.5,'Coverage • sign certification • regret',7.8,C['gold'],'bold')
    txt(54,15.7,'108,000 simulations; 3,000 per cell',7.5)
    txt(1,9,'Match the remedy to the source of uncertainty',10,weight='bold')
    txt(1,4.9,'Unobserved outcomes, discarded constraints and limited samples require different responses.',7.6)
    export(fig,'figure_1')

def jointplot():
    j=json.loads((E/'joint_summary.json').read_text());data=j['summary']
    fig,axs=plt.subplots(1,2,figsize=(170/25.4,84/25.4),gridspec_kw={'width_ratios':[1.25,1]});fig.subplots_adjust(left=.11,right=.98,bottom=.20,top=.78,wspace=.62)
    for i,(world,rec) in enumerate([(w,r) for w in ('B','C') for r in ('accurate','fallible')]):
        for k,(m,c,marker) in enumerate([('oracle','blue','o'),('bounded','teal','s')]):
            x=next(d for d in data if d['world']==world and d['recording']==rec and d['method']==m);yy=i+(k-.5)*.22
            val=x['width_reduction'];axs[0].errorbar(val['mean'],yy,xerr=1.96*val['mcse'],fmt=marker,color=C[c],ms=4,capsize=2,label=m.capitalize() if i==0 else None)
    axs[0].set_yticks(range(4),['B / accurate','B / fallible','C / accurate','C / fallible']);axs[0].invert_yaxis();axs[0].set_xlim(0,.045);axs[0].set_xlabel('Paired reduction in envelope width');axs[0].set_title('a  Finite samples',loc='left',weight='bold');axs[0].legend(loc='lower left',bbox_to_anchor=(-.07,1.08),ncol=2,fontsize=7)
    old_summary=json.loads((OLD.parent/'extensions/bounded_correction/summary_main.json').read_text())
    old_pop=next(x for x in old_summary['population_outer_envelopes_not_sharp'] if x['world']=='B' and x['recording']=='accurate' and x['method']=='oracle')
    new_pop=next(x for x in j['population'] if x['world']=='B' and x['recording']=='accurate' and x['method']=='oracle' and x['grid']==256)
    pop=[('Separate outer',old_pop['lower'],old_pop['upper'],C['muted'],'o'),('Joint cell outer',new_pop['lower'],new_pop['upper'],C['blue'],'s'),('Exact joint set',-1/72,-1/72,C['red'],'D')]
    for i,(label,lo,hi,c,m) in enumerate(pop):axs[1].plot([lo,hi],[i,i],color=c,lw=2);axs[1].plot((lo+hi)/2,i,m,ms=4,color=c)
    axs[1].axvline(0,color=C['muted'],ls='--',lw=.7);axs[1].set_yticks(range(3),[x[0] for x in pop]);axs[1].invert_yaxis();axs[1].set_xlim(-.018,.007);axs[1].set_xlabel('Net-benefit difference');axs[1].set_title('b  Population B / oracle',loc='left',weight='bold');axs[1].ticklabel_format(axis='x',style='plain');axs[1].set_xticks([-.015,0])
    # Read the plotted cell interval from source, not a display constant.
    p=next(x for x in j['population'] if x['world']=='B' and x['recording']=='accurate' and x['method']=='oracle' and x['grid']==256)
    assert abs(p['lower']-pop[1][1])<1e-9 and abs(p['upper']-pop[1][2])<1e-9
    fig.text(.01,.035,'a: 3,000 paired records per cell; bars ±1.96 Monte Carlo SE.  b: population objects, no confidence level.',fontsize=6.6,color=C['muted'])
    export(fig,'figure_6')

def budgetplot():
    data=json.loads((E/'budget_summary.json').read_text())['summary'];fig,axs=plt.subplots(2,2,figsize=(170/25.4,132/25.4));fig.subplots_adjust(left=.12,right=.98,bottom=.17,top=.86,wspace=.30,hspace=.57)
    configs=[('A','correct_sign',100,'a  Positive world A','Correct sign certified (%)'),('B','correct_sign',100,'b  Negative world B','Correct sign certified (%)'),('B','regret',1,'c  Decision loss in world B','Mean net-benefit regret'),('C','wrong_sign',100,'d  False sign in tie world C','Wrong sign certified (%)')]
    for ax,(w,metric,mul,title,ylab) in zip(axs.flat,configs):
        for alloc,c,mark,ls in [(1/3,'red','^',':'),(.5,'purple','s','--'),(2/3,'blue','o','-')]:
            rows=sorted([x for x in data if x['world']==w and abs(x['allocation']-alloc)<1e-8],key=lambda x:x['budget']);x=np.arange(4)
            if metric=='regret':y=np.array([r[metric]['mean'] for r in rows]);err=np.array([1.96*r[metric]['mcse'] for r in rows])
            else:
                y=np.array([r[metric]['rate'] for r in rows]);err=np.array([[r[metric]['rate']-r[metric]['low'] for r in rows],[r[metric]['high']-r[metric]['rate'] for r in rows]])
            ax.errorbar(x,y*mul,yerr=err*mul,color=C[c],marker=mark,ls=ls,ms=3.8,capsize=2,label={1/3:'1/3 in G+',.5:'1/2 in G+',2/3:'2/3 in G+'}[alloc])
        ax.set_title(title,loc='left',weight='bold');ax.set_xticks(range(4),['200','800','3,200','12,800']);ax.set_xlabel('Total new reference outcomes');ax.set_ylabel(ylab);ax.set_ylim(bottom=0)
        if metric=='correct_sign':ax.set_ylim(0,105)
    handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.52,.99),ncol=3,fontsize=8)
    fig.text(.01,.035,'Equal per-case costs; random, accurate reference ascertainment. Bars show Monte Carlo uncertainty.\nNo source/proxy data enter any allocation. Regret uses an explicit two-rule minimax choice.',fontsize=7,color=C['muted'])
    export(fig,'figure_7')

def inherited():
    for name in ['figure_2','figure_3','figure_4','figure_5','figure_S1','figure_S2']:
        src=fitz.open(OLD/f'{name}.pdf');p=src[0];width=170/25.4*72;scale=width/p.rect.width;dst=fitz.open();new=dst.new_page(width=width,height=p.rect.height*scale);new.show_pdf_page(new.rect,src,0);dst.save(R/f'{name}.pdf');pix=new.get_pixmap(dpi=600,alpha=False);pix.save(R/f'{name}.png');dst.close();src.close()
        im=Image.open(R/f'{name}.png').convert('RGB');im.save(R/f'{name}.tiff',compression='tiff_lzw',dpi=(600,600));im.save(R/f'{name}.jpeg',quality=100,subsampling=0,dpi=(600,600))
        shutil.copy2(OLD/f'{name}.svg',R/f'{name}.svg')
    (R/'INHERITED_FIGURES.md').write_text('Figures 2–5 and S1–S2 preserve previous verified figure contents and source SVG, resized as vector PDFs to 170 mm width. TIFF/JPEG are rendered directly from these vectors at 600 dpi. Original editable plotting sources and data remain in the computational companion. New figures 1, 6 and 7 are generated by build_figures.py.\n',encoding='utf-8')

if __name__=='__main__':
    schematic();jointplot();budgetplot();inherited()
    manifest=[]
    for p in sorted(R.glob('*.tiff')):
        im=Image.open(p);manifest.append(dict(file=p.name,pixels=im.size,dpi=[float(v) for v in im.info.get('dpi')],bytes=p.stat().st_size));assert p.stat().st_size<10_000_000
    (R/'figure_inventory.json').write_text(json.dumps(manifest,indent=2));print('Nine figure sets exported')
