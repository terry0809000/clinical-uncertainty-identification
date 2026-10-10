"""Editable graphical overview. Geometric icons; no generated patient imagery."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Circle,Rectangle,FancyArrowPatch
from PIL import Image
OUT=Path(__file__).resolve().parent
C={'ink':'#14283F','blue':'#1478B5','teal':'#087E80','purple':'#7854A1','orange':'#B97512','red':'#BB4057','line':'#CEDAE4','gray':'#526678'}
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none'})
fig,ax=plt.subplots(figsize=(7.1,5.55)); fig.subplots_adjust(0,0,1,1); ax.set(xlim=(0,1),ylim=(0,1)); ax.axis('off')
def text(x,y,s,size=8,color='ink',bold=False,ha='left',**kw):
    ax.text(x,y,s,fontsize=size,color=C.get(color,color),weight='bold' if bold else 'normal',ha=ha,va='center',**kw)
def box(x,y,w,h,face='white',edge='line',radius=.015):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f'round,pad=0.003,rounding_size={radius}',facecolor=C.get(face,face),edgecolor=C.get(edge,edge),lw=.75))
def arrow(x,y,xx,yy,color='gray'):
    ax.add_patch(FancyArrowPatch((x,y),(xx,yy),arrowstyle='-|>',mutation_scale=8,lw=1,color=C.get(color,color)))
def person(x,y,color):
    ax.add_patch(Circle((x,y),.009,fc=C[color],ec='none'))
    ax.add_patch(FancyBboxPatch((x-.007,y-.031),.014,.020,boxstyle='round,pad=.002',fc=C[color],ec='none'))
def hospital(x,y,color):
    ax.add_patch(Rectangle((x,y),.062,.056,fc='white',ec=C[color],lw=1))
    ax.add_patch(Rectangle((x+.025,y+.034),.013,.004,fc=C[color]))
    ax.add_patch(Rectangle((x+.0295,y+.029),.004,.014,fc=C[color]))
    for xx in (x+.011,x+.026,x+.041):
        ax.add_patch(Rectangle((xx,y+.012),.008,.009,fc=C[color]))

text(.035,.965,'What evidence makes a clinical prediction actionable?',11.1,bold=True)
text(.035,.929,'A critical narrative synthesis with exact and finite-sample synthetic experiments',7.6,color='gray')
stages=[('TARGET','Who, when, under what care?','blue'),('OBSERVATION','What was measured?','teal'),('IDENTIFICATION','What is compatible?','purple'),('INFERENCE','What does an interval cover?','orange'),('DECISION','Which consequence matters?','red')]
for i,(title,sub,col) in enumerate(stages):
    x=.035+i*.191
    box(x,.820,.168,.080,face=col,edge=col,radius=.009)
    text(x+.010,.872,title,8,color='white',bold=True)
    # Two-line subtitles fit final-size typography.
    subtitles=['Who, when,\nunder what care?','What was\nmeasured?','What is\ncompatible?','What does an\ninterval cover?','Which consequence\nmatters?']
    text(x+.010,.841,subtitles[i],6.9,color='white',linespacing=1.25)
    if i<4: arrow(x+.171,.861,x+.188,.861)

box(.035,.330,.448,.458,face='#F4F8FC')
box(.515,.330,.450,.458,face='#F8F5FC')
text(.050,.761,'a  DEFINE THE OBSERVATION PROCESS',8.4,bold=True,color='blue')
hospital(.061,.649,'blue'); text(.135,.689,'Source setting',8.1,bold=True)
text(.135,.650,'Paired reference\n+ routine proxy',7.1,color='gray',linespacing=1.2)
hospital(.300,.649,'teal'); text(.380,.689,'Target',8.1,bold=True)
text(.380,.660,'Proxy only',7.1,color='gray')
arrow(.279,.678,.294,.678)
for i in range(4): person(.080+i*.031,.603,'blue')
text(.224,.591,'Fixed care policy',7.2,color='gray')
text(.057,.539,'Two locked rules disagree',8.2,bold=True)
box(.058,.446,.183,.062,face='#E4F1FA',edge='#B6D6E8',radius=.007)
box(.265,.446,.195,.062,face='#EEE7F8',edge='#D1C0E3',radius=.007)
text(.149,.481,r'$G_+$: rule 1 only',7.6,color='blue',bold=True,ha='center')
text(.362,.481,r'$G_-$: rule 2 only',7.6,color='purple',bold=True,ha='center')
text(.149,.461,r'Weight $2/3$',6.7,ha='center',color='gray')
text(.362,.461,r'Weight $1/3$',6.7,ha='center',color='gray')
text(.057,.408,'Routine proxy ≠ clinical reference',8.0,color='red',bold=True)
text(.057,.366,'The endpoint and its ascertainment are part of the target.',7.0,color='gray')

text(.533,.761,'b  HOLD OBSERVED DATA FIXED',8.4,bold=True,color='purple')
text(.533,.719,'Compatible worlds can reverse the ranking',8.1,bold=True)
text(.533,.687,'Reference-outcome net-benefit difference at t = 0.20',7.1,color='gray')
x0=.700; scale=2.5
ax.plot([.569,.914],[.529,.529],color=C['line'],lw=.8)
ax.plot([x0,x0],[.523,.657],color=C['gray'],lw=.65,ls='--')
for y,val,label,col,mark in ((.638,-1/72,'B  negative','red','v'),(.591,0,'C  tie','gray','s'),(.544,1/18,'A  positive','teal','o')):
    text(.540,y,label,6.9,color=col)
    ax.plot(x0+val*scale,y,marker=mark,ms=5.1,color=C[col])
    text(.937,y,f'{val:+.3f}',7.0,color=col,ha='right')
text(x0,.503,'0',6.8,color='gray',ha='center')
text(.533,.459,'Identification set',7.6,bold=True)
ax.plot([x0-scale/72,x0+scale/18],[.422,.422],color=C['purple'],lw=4,solid_capstyle='butt')
text(.533,.384,'Same proxy data · different compatible reference risks',7.0,color='gray')
text(.533,.353,'Population calculation, not a confidence interval',7.0,color='gray')

text(.035,.293,'c  TEST WHICH ADDITIONAL EVIDENCE CHANGES THE CONCLUSION',8.5,bold=True)
items=[(.035,'01  Exact identification','Shared versus group-specific\nlabel channels; audit information','purple'),(.350,'02  Finite samples','Coverage, width and sign calls\nwith Monte Carlo uncertainty','blue'),(.665,'03  Audit validity','Selection and reference error;\nknown or bounded correction','teal')]
for x,title,sub,col in items:
    box(x,.152,.299,.106,face='white',edge=col,radius=.008)
    text(x+.012,.231,title,8.1,color=col,bold=True)
    text(x+.012,.188,sub,7.2,color='gray',linespacing=1.35)
box(.035,.038,.929,.079,face='#EAF3F5',edge='#CDDFE4',radius=.009)
text(.051,.087,'Interpretation',8.1,bold=True,color='teal')
text(.218,.087,'Match the next study to the missing information.',8.5,bold=True)
text(.218,.058,'More proxy observations alone need not resolve a reference-outcome decision.',7.3,color='gray')
for ext in ('pdf','svg','eps','png'):
    fig.savefig(OUT/f'figure_1.{ext}',dpi=600,facecolor='white')
plt.close(fig)
with Image.open(OUT/'figure_1.png') as im:
    im.convert('RGB').save(OUT/'figure_1.tiff',compression='tiff_lzw',dpi=(600,600))
    im.convert('RGB').save(OUT/'figure_1.jpeg',quality=96,subsampling=0,dpi=(600,600))
print('Figure 1 exported: PDF SVG EPS PNG TIFF JPEG')
