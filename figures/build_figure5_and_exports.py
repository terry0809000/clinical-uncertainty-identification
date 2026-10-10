from pathlib import Path
from fractions import Fraction as F
import json,shutil,subprocess,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'figures'; OLD=ROOT/'legacy'
PROJ=ROOT.parents[1]
new=json.loads((ROOT/'extensions/group_channels/group_channels_results.json').read_text())['rows']
sim=json.loads((ROOT/'extensions/bounded_correction/summary_main.json').read_text())
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
fig,(ax,bx)=plt.subplots(1,2,figsize=(7.1,4.0),gridspec_kw={'width_ratios':[1,1.16]})
fig.subplots_adjust(left=.14,right=.98,bottom=.20,top=.82,wspace=.38)
colors=['#7854A1','#087E80']
labels=['No audit','A  exact G−','B  exact G−','C  exact G+']
look=['none','A_minus','B_minus','C_plus']
for j,ch in enumerate(['shared','group_specific']):
    for i,a in enumerate(look):
        r=next(x for x in new if x['risk_cap']=='1/10' and x['channel_cap']=='1/20' and x['channel']==ch and x['audit']==a)
        y=3-i+(.12 if j==0 else -.12); lo,hi=float(F(r['lower'])),float(F(r['upper']))
        ax.plot([lo,hi],[y,y],lw=2.5,color=colors[j],solid_capstyle='butt')
        ax.plot([lo,hi],[y,y],ls='none',marker='|' if j==0 else 'o',ms=4,color=colors[j],label=['Shared channel','Group-specific'][j] if i==0 else None)
ax.axvline(0,color='#798894',ls='--',lw=.7)
ax.set(yticks=[3,2,1,0],yticklabels=labels,xlim=(-.095,.14),xticks=[-.08,0,.08],xlabel='Net-benefit difference')
ax.set_title('a  Structural information',loc='left',fontweight='bold',fontsize=9)
ax.legend(loc='lower left',bbox_to_anchor=(-.10,1.08),fontsize=7,frameon=False)
ax.tick_params(axis='y',length=0); ax.spines['left'].set_visible(False)
methods=['naive','oracle','bounded','wrong_selection_bound']; cols=['#526678','#1478B5','#087E80','#BB4057']; marks=['o','s','D','^']
for k,(world,ref) in enumerate([('B','accurate'),('B','fallible'),('C','accurate'),('C','fallible')]):
    for j,m in enumerate(methods):
        row=next(x for x in sim['cells'] if x['seed']=='pooled' and x['world']==world and x['recording']==ref and x['method']==m)
        z=row['metrics']['wrong_sign']; val=100*z['proportion']; lo,hi=np.array(z['mc_95_interval'])*100
        bx.errorbar(val,3-k+(1.5-j)*.13,xerr=[[val-lo],[hi-val]],fmt=marks[j],color=cols[j],ms=3.6,lw=.85,capsize=1.8,label=['Naive','Oracle-known','Truth-containing bounds','False selection bounds'][j] if k==0 else None)
bx.set(yticks=[3,2,1,0],yticklabels=['B / accurate','B / fallible','C / accurate','C / fallible'],xlim=(-3,100),xticks=[0,25,50,75,100],xlabel='Wrong nonzero sign (%)')
bx.set_title('b  Audit mechanism information',loc='left',fontweight='bold',fontsize=9)
bx.grid(axis='x',color='#DEE6E9',lw=.55); bx.set_axisbelow(True); bx.tick_params(axis='y',length=0); bx.spines['left'].set_visible(False)
bx.legend(loc='lower left',bbox_to_anchor=(-.12,1.07),ncol=2,fontsize=6.8,frameon=False,columnspacing=.8,handletextpad=.4)
fig.text(.025,.086,'a: exact population sets under original sensitivity caps.  b: 3,000 repetitions per cell.',fontsize=7,color='#526678')
fig.text(.025,.045,'Oracle and bounded corrections made no sign calls. Error bars: exact 95% Monte Carlo intervals.',fontsize=7,color='#526678')
for ext in ('pdf','svg','eps','png'): fig.savefig(OUT/f'figure_5.{ext}',dpi=600,facecolor='white')
plt.close(fig)

# Preserve originals, copy styles/results, repair Figure2 bottom note clearance.
mapping={'figure2_identification':'figure_2','figure3_simulation':'figure_3','figure4_stress_test':'figure_4','figureS1_sensitivity_surface':'figure_S1','figureS2_threshold_criterion':'figure_S2'}
for stem,dest in mapping.items():
    for ext in ('pdf','svg','png'):
        shutil.copy2(OLD/'figures'/f'{stem}.{ext}',OUT/f'{dest}.{ext}')
    source=(OLD/'figures'/f'{stem}.py').read_text(encoding='utf-8')
    source=source.replace('Path(__file__).resolve().parents[1] / "simulation"','Path(__file__).resolve().parents[1] / "legacy" / "simulation"')
    source=source.replace('Path(__file__).resolve().parents[1] / "analytic"','Path(__file__).resolve().parents[1] / "legacy" / "analytic"')
    source=source.replace('OUT.parent / "analytic"', 'Path(__file__).resolve().parents[1] / "legacy" / "analytic"')
    if stem=='figure2_identification':
        source=source.replace('(6.69, 3.32)','(6.69, 3.6)').replace('0.23, 0.29, 0.60','0.30, 0.29, 0.52').replace('0.23, 0.43, 0.60','0.30, 0.43, 0.52')
        source=source.replace('("pdf", "svg", "png")','("pdf", "svg", "eps", "png")').replace('figure2_identification.{ext}','figure_2.{ext}')
        f=OUT/f'{stem}.py'; f.write_text(source,encoding='utf-8'); subprocess.run([sys.executable,str(f)],check=True)
    else: (OUT/f'{stem}.py').write_text(source,encoding='utf-8')

inventory=[]
for name in ['figure_1','figure_2','figure_3','figure_4','figure_5','figure_S1','figure_S2']:
    with Image.open(OUT/f'{name}.png') as im:
        rgb=im.convert('RGB'); rgb.save(OUT/f'{name}.tiff',compression='tiff_lzw',dpi=(600,600)); rgb.save(OUT/f'{name}.jpeg',quality=96,subsampling=0,dpi=(600,600))
        inventory.append({'figure':name,'pixels':list(im.size),'dpi':600,'tiff':'RGB LZW','editable':'SVG + plotting source','vector':'PDF'})
(OUT/'figure_inventory.json').write_text(json.dumps(inventory,indent=2),encoding='utf-8')
print('7 figures exported with TIFF/JPEG and vector sources')
