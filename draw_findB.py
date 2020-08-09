import pandas as pd
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib
import matplotlib.pyplot as plt
import sys

matplotlib.rcParams['xtick.direction'] = 'in'
matplotlib.rcParams['xtick.top'] = True
matplotlib.rcParams['xtick.minor.visible'] = True
matplotlib.rcParams['xtick.major.size'] = 10
matplotlib.rcParams['xtick.minor.size'] = 5
matplotlib.rcParams['ytick.direction'] = 'in'
matplotlib.rcParams['ytick.right'] = True
matplotlib.rcParams['ytick.minor.visible'] = True
matplotlib.rcParams['ytick.major.size'] = 10
matplotlib.rcParams['ytick.minor.size'] = 5
                                    
df = pd.read_hdf(sys.argv[1], key='cats')

new_cats = [(df['nb']==2) & (df['nsv']==0),
            (df['nb']==2) & (df['nsv']==1),
            (df['nb']==2) & (df['nsv']==2),
            (df['nb']==3) & (df['nsv']==0),
            (df['nb']==3) & (df['nsv']==1),
            (df['nb']==4) & (df['nsv']==0)]

old_cats = [(df['nb']==2),
            (df['nb']==2) & (df['nsv']==100),
            (df['nb']==2) & (df['nsv']==100),
            (df['nb']==3),
            (df['nb']==3) & (df['nsv']==100),
            (df['nb']==4)]

name_cats = [r'$n_{b} = 2$, $n_{SV}=0$',
             r'$n_{b} = 2$, $n_{SV}=1$',
             r'$n_{b} = 2$, $n_{SV}=2$',
             r'$n_{b} = 3$, $n_{SV}=0$',
             r'$n_{b} = 3$, $n_{SV}=1$',
             r'$n_{b} = 4$, $n_{SV}=0$']

new_count = np.array([float(len(np.where(cat)[0]))/float(len(df)) for cat in new_cats])
old_count = np.array([float(len(np.where(cat)[0]))/float(len(df)) for cat in old_cats])
x = np.arange(len(name_cats))  # the label locations
print new_count, old_count, x
width = 0.35  # the width of the bars

fig, ax = plt.subplots(figsize=(12, 6), dpi=200, facecolor='w', edgecolor='k')
rects1 = ax.bar(x - width/2, old_count, width, label='Only jets')
rects2 = ax.bar(x + width/2, new_count, width, label='Jets and SV')

ax.set_ylabel('Number of Entries')
ax.set_title('Acceptance and 4-body mass per topological region')
ax.set_xticks(x)
ax.set_xticklabels(name_cats)
ax.set_ylim(ymin=0., ymax=0.55)
ax.legend(loc=2)


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('%.2f' % height,
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
        

autolabel(rects1)
autolabel(rects2)

#fig.tight_layout()

ax2 = inset_axes(ax,
                 width='100%',
                 height='100%',
                 bbox_to_anchor=(0.39, 0.7, 0.11, 0.3),
                 bbox_transform=ax.transAxes)
ax3 = inset_axes(ax,
                 width='100%',
                 height='100%',
                 bbox_to_anchor=(0.55, 0.7, 0.11, 0.3),
                 bbox_transform=ax.transAxes)
ax4 = inset_axes(ax,
                 width='100%',
                 height='100%',
                 bbox_to_anchor=(0.72, 0.7, 0.11, 0.3),
                 bbox_transform=ax.transAxes)
ax5 = inset_axes(ax,
                 width='100%',
                 height='100%',
                 bbox_to_anchor=(0.88, 0.7, 0.11, 0.3),
                 bbox_transform=ax.transAxes)

ax2.hist(df[new_cats[2]]['hm'],bins=45,range=[50,140],histtype='step',color='tab:red',linewidth=2)
ax2.set_xlabel('Mass [GeV]')
ax2.set_ylabel('Entries/2 GeV')
ax2.set_ylim(ymax=250)
ax2.text(0.15,0.8, '$\sigma = {:03.1f}\,$ GeV'.format(df[new_cats[2]]['hm'].std()), transform=ax2.transAxes, fontsize=8)
ax3.hist(df[new_cats[3]]['hm'],bins=45,range=[50,140],histtype='step',color='tab:red',linewidth=2)
ax3.set_xlabel('Mass [GeV]')
ax3.set_ylim(ymax=570)
ax3.text(0.15,0.8, '$\sigma = {:03.1f}\,$ GeV'.format(df[new_cats[3]]['hm'].std()), transform=ax3.transAxes, fontsize=8)
ax4.hist(df[new_cats[4]]['hm'],bins=45,range=[50,140],histtype='step',color='tab:red',linewidth=2)
ax4.set_xlabel('Mass [GeV]')
ax4.set_ylim(ymax=500)
ax4.text(0.15,0.8, '$\sigma = {:03.1f}\,$ GeV'.format(df[new_cats[4]]['hm'].std()), transform=ax4.transAxes, fontsize=8)
ax5.hist(df[new_cats[5]]['hm'],bins=45,range=[50,140],histtype='step',color='tab:red',linewidth=2)
ax5.set_xlabel('Mass [GeV]')
ax5.set_ylim(ymax=310)
ax5.text(0.15,0.8, '$\sigma = {:03.1f}\,$ GeV'.format(df[new_cats[5]]['hm'].std()), transform=ax5.transAxes, fontsize=8)
ax.text(0.04, 0.83, "ATLAS", fontweight='bold', fontstyle='italic', verticalalignment='bottom', transform=ax.transAxes)
ax.text(0.1, 0.83, "Simulation Internal", verticalalignment='bottom', transform=ax.transAxes)
ax.text(0.04, 0.77, r'$ZH, H\rightarrow aa\rightarrow 4b$', verticalalignment='bottom', transform=ax.transAxes)
ax.text(0.04, 0.73, r'$m_{a} = 50\,$GeV', verticalalignment='bottom', transform=ax.transAxes)

plt.savefig('categories.png')
