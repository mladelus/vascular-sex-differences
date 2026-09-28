import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams['font.family']='DejaVu Sans'
W,H=6.69,5.1
fig=plt.figure(figsize=(W,H)); ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,100); ax.set_ylim(0,100); ax.axis('off')
rows=[('A','Healthy capillary endothelium','#1E88E5','#E3F2FD'),
      ('B','Arteries (GTEx bulk RNA-seq)','#8E24AA','#F3E5F5'),
      ('C','Kidney capillaries (KPMP)','#E53935','#FFEBEE'),
      ('D','Sex chromosomes and hormone receptors','#2E7D32','#E8F5E9')]
cols=[(7,26),(30.5,49.5),(54,73),(77,99)]
rowy=[(75.6,99.2),(50.8,74.4),(26,49.6),(0.8,24.8)]   # (bottom, top) of each band incl. title
def band(i):
    L,t,c,l=rows[i]; y0,y1=rowy[i]
    ax.add_patch(FancyBboxPatch((0.8,y0),98.9,y1-y0,boxstyle='round,pad=0,rounding_size=1.2',fc=l,ec='none',alpha=0.55,zorder=0))
    ax.add_patch(FancyBboxPatch((0.8,y0),4.2,y1-y0,boxstyle='round,pad=0,rounding_size=1.2',fc=c,ec='none',zorder=1))
    ax.text(2.9,(y0+y1)/2,L,color='white',fontsize=12,fontweight='bold',ha='center',va='center',zorder=2)
    ax.text(7,y1-2.7,t,color=c,fontsize=7.8,fontweight='bold',ha='left',va='center')
def card(x0,x1,y0,y1,title,sub,res,c,final=False,l='#fff',k=1.0):
    fc=c if final else 'white'
    ax.add_patch(FancyBboxPatch((x0,y0),x1-x0,y1-y0,boxstyle='round,pad=0,rounding_size=1.3',
                 fc=fc,ec=c,lw=1.1 if not final else 0,zorder=3))
    cx=(x0+x1)/2
    if final:
        ax.text(cx,(y0+y1)/2,title,color='white',fontsize=6.7,fontweight='bold',ha='center',va='center',zorder=4,linespacing=1.3)
        return
    lines=[(title,dict(fontsize=6.5,fontweight='bold',color='#212121'))]
    if sub: lines.append((sub,dict(fontsize=5.8,color='#616161')))
    if res: lines.append((res,dict(fontsize=5.9,color=c,fontstyle='italic',fontweight='bold')))
    n=sum(s.count('\n')+1 for s,_ in lines); step=2.25*6.3/H*k
    y=(y0+y1)/2+step*(n-1)/2
    for s,kw in lines:
        k=s.count('\n')+1
        ax.text(cx,y-step*(k-1)/2,s,ha='center',va='center',zorder=4,linespacing=1.15,**kw); y-=step*k
def arrow(x0,y0,x1,y1,c):
    ax.add_patch(FancyArrowPatch((x0,y0),(x1,y1),arrowstyle='-|>',mutation_scale=8,color=c,lw=1.1,zorder=2,shrinkA=0,shrinkB=0))
for i in range(4): band(i)
# Row A
c='#1E88E5'; b,t=rowy[0]; t-=5.4; b+=1.6; m=(b+t)/2
card(*cols[0],b,t,'Tabula Sapiens 2.0','6 W · 8 M\nmany organs','no sex difference',c)
card(*cols[1],b,t,'Heart Cell Atlas','7 W · 7 M\nmyocardium','hypoxia, TGF-β\nhigher in women',c)
card(*cols[2],m+0.6,t,'Independent LV','15 W · 14 M','not replicated',c,k=0.8)
card(*cols[2],b,m-0.6,'GTEx left ventricle','137 W · 294 M','not replicated',c,k=0.8)
card(*cols[3],b,t,'No reproducible\nsex difference in\nhealthy capillaries',None,None,c,final=True)
arrow(cols[0][1],m,cols[1][0],m,c); arrow(cols[1][1],m,cols[2][0],(m+t)/2+0.4,c); arrow(cols[1][1],m,cols[2][0],(b+m)/2-0.4,c)
arrow(cols[2][1],(m+t)/2+0.4,cols[3][0],m+1.5,c); arrow(cols[2][1],(b+m)/2-0.4,cols[3][0],m-1.5,c)
# Row B
c='#8E24AA'; b,t=rowy[1]; t-=5.4; b+=1.6; m=(b+t)/2
card(*cols[0],b,t,'GTEx arteries','coronary 94 W · 146 M\naorta 152 W · 279 M\ntibial 206 W · 445 M',None,c)
card(*cols[1],b,t,'Adjusted models','age, ischemic time,\nmanner of death, RIN',None,c)
card(*cols[2],b,t,'Sensitivity','self-contained test (fry),\nimmune cells, ventilator\ndeaths, under vs over 50',None,c)
card(*cols[3],b,t,'Coronary, under 50:\nlower interferon-γ and\ninflammatory signaling\nin women',None,None,c,final=True)
for k in range(3): arrow(cols[k][1],m,cols[k+1][0],m,c)
# Row C
c='#E53935'; b,t=rowy[2]; t-=5.4; b+=1.6; m=(b+t)/2
card(*cols[0],b,t,'KPMP kidney atlas','healthy 24 · CKD 39\ndonors',None,c)
card(*cols[1],b,t,'Sex × CKD interaction','peritubular and\nglomerular capillaries',None,c)
card(*cols[2],b,t,'Robustness','self-contained test (fry),\nimmune-RNA adjustment',None,c)
card(*cols[3],b,t,'Modest sex × CKD\ninteraction in\nglomerular capillary\ninterferon-γ signaling',None,None,c,final=True)
for k in range(3): arrow(cols[k][1],m,cols[k+1][0],m,c)
# Row D
c='#2E7D32'; b,t=rowy[3]; t-=5.4; b+=1.6; m=(b+t)/2
card(*cols[0],b,t,'3 single-nucleus\ncohorts','+ KPMP single cells',None,c)
card(*cols[1],b,t,'X–Y paralogue pairs','KDM6A/UTY,\nKDM5C/KDM5D','X copy higher in women;\nY copy 52–74% in men',c)
card(*cols[2],m+0.6,t,'Loss of Y rare','1.6% (kidney)',None,c,k=0.8)
card(*cols[2],b,m-0.6,'Hormone receptors','AR, ESR1, PGR; no aromatase',None,c,k=0.8)
card(*cols[3],b,t,'In all 3 cohorts:\nlower combined\nX + Y dosage of\nchromatin regulators\nin women',None,None,c,final=True)
arrow(cols[0][1],m,cols[1][0],m,c); arrow(cols[1][1],m,cols[2][0],(m+t)/2+0.4,c); arrow(cols[1][1],m,cols[2][0],(b+m)/2-0.4,c)
arrow(cols[2][1],(m+t)/2+0.4,cols[3][0],m+1.5,c); arrow(cols[2][1],(b+m)/2-0.4,cols[3][0],m-1.5,c)
fig.savefig('fig1_study_design.png',dpi=300,facecolor='white')
fig.savefig('fig1_study_design.pdf',facecolor='white')
fig.savefig('fig1_study_design.tiff',dpi=600,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
