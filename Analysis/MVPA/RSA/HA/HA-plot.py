# Import results from clf_within and plot them.
# Mostly copied code from that script, since it
# outputs the plots for each ROI.
# Here the masks are combined into one plot.

from collections import OrderedDict

import os
import argparse
import pandas as pd
import seaborn as sea
import matplotlib.pyplot as plt; plt.ion()
from matplotlib.patches import Patch
from matplotlib import rcParams

rcParams['axes.labelsize'] = 'x-large'
rcParams['xtick.labelsize'] = 'x-large'
rcParams['ytick.labelsize'] = 'large'


parser = argparse.ArgumentParser()
parser.add_argument('--mask',default='vtc',choices=['vtc','ips'],type=str)
parser.add_argument('--mode',default='ha',choices=['ha','ws'],type=str)
args = parser.parse_args()


MASK = args.mask
MODE = args.mode

home = os.path.expanduser('~')
res_dir  = os.path.join(home,'DBw/STUDY/RetroBlast/results')
if MASK == 'ips':
    res_dir = os.path.join(res_dir,MASK)

fn_template = os.path.join(res_dir,'within_del-{:s}-{:s}.csv')

if MODE == 'ha':
    export_base = 'fig3'
else:
    export_base = 'figS5'


# load data
df_acc = pd.read_csv(fn_template.format(MODE,'acc'))
df_null = pd.read_csv(fn_template.format(MODE,'null'))


# if within-subjects have to rearrange
# first get average accuracy per subj/cond (ie, across cv folds)
if MODE == 'ws':
    df_acc = df_acc.groupby(['subj','cond']
        ).aggregate({'accuracy':'mean'}
        ).reset_index()

COND_PALETTE = {'d1':        'lightgray',
                'd2-stay':   'tomato',
                'd2-switch': 'royalblue',
                'blast':     'gold'}
COND_XPOS = OrderedDict([('blast',     0.5),
                         ('d1',        0),
                         ('d2-stay',   1.0),
                         ('d2-switch', 1.2)])



# plot mean/sem across HA cv folds

fig, ax = plt.subplots(figsize=(6,7))

for cond, series in df_acc.groupby('cond').accuracy:
    y, yerr = series.mean(), series.sem()
    x = COND_XPOS[cond]
    col = COND_PALETTE[cond]
    ax.errorbar(x,y,yerr,c=col,marker='o',capsize=5,
        markersize=10,elinewidth=3,capthick=3,zorder=10)


# plot null distributions as violins
for cond, series in df_null.groupby('cond').accuracy:
    x = COND_XPOS[cond]
    # get null distribution
    if MODE == 'ha':
        null_distrn = series.values
    else:
        n_samples_persubj = series.size / df_acc.subj.nunique()
        # bootstrap that number of draws from all perms across subjs per cond
        null_distrn = pd.np.random.choice(series,size=n_samples_persubj,replace=True)
    viols = ax.violinplot([null_distrn],positions=[x],widths=[.125],showextrema=False)
    plt.setp(viols['bodies'],facecolor='whitesmoke',edgecolor='white')

# viols = ax.violinplot([vtc_distrn,ips_distrn],positions=[0,1],widths=[.125,.125],showextrema=False)
# plt.setp(viols['bodies'],facecolor='gray',edgecolor='white')

plt.axhline(.5,linestyle='--',linewidth=.5,color='black')
plt.xlim(-0.3,1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xticks(COND_XPOS.values())
ax.set_xticklabels(COND_XPOS.keys(),rotation=45)
plt.ylabel('Classifier Accuracy')
plt.xlabel('')
plt.tight_layout()


fname = os.path.join(res_dir,'{:s}.eps'.format(export_base))
plt.savefig(fname,dpi=300)
plt.close()



# # df_list = []
# # for metr in metrics:
# #     fn = fn_template.format(metr)
# #     df = pd.read_csv(fn)
# #     df.index = pd.Index(pd.np.repeat(metr,df.shape[0]),name='metric')
# #     df_list.append(df)
# # df = pd.concat(df_list).reset_index()


# # make background facetgrid with violinplots
# mask_order = masks
# metric_order = metrics
# event_order = ['delay1','delay2']
# ttype_order = ['stay','switch']
# ttype_palette = dict(stay='red',switch='royalblue')


# ##--- plotting both stay/switch trials in delay1 ---##

# g = sea.catplot(data=df[df.metric=='null'],
#         x='mask',y='accuracy',hue='ttype',col='event',
#         kind='violin',palette=ttype_palette,order=mask_order,
#         hue_order=ttype_order,col_order=event_order,
#         split=True,
#         linewidth=0,
#         saturation=1,
#         width=.6, # of violinplots
#         inner=None, # {box, quartile, point, stick, None}
#         cut=0, # "dont let density extend past extreme values"
#         aspect=.7,
#         # bw=1, # controls amount of smoothing
#     ).set_axis_labels('','Classifier accuracy'
#     ).set_xticklabels(['VTC','IPS']
#     ).set_titles('{col_name}')
# g._legend.remove()

# # set alpha of violins before drawing more
# for ax in g.axes.flat:
#     plt.setp(ax.collections,alpha=.2)
#     ax.axhline(.5,linestyle='--',linewidth=.5,color='black',zorder=0)

# # get mean and sem for all conditions
# plot_df = df.query("metric=='acc'"
#            ).groupby(['mask','event','ttype']
#            ).aggregate({'accuracy':['mean','sem']})

# # overlay test results
# for indx, row in plot_df.iterrows():
#     mask, event, ttype = indx
#     ax = g.axes.flat[event_order.index(event)]
#     cstr = ttype_palette[ttype]
#     xval = mask_order.index(mask)
#     xval += .13 * (-1 if ttype=='stay' else 1)
#     yval = row.accuracy['mean']
#     err  = row.accuracy['sem']
#     ax.errorbar(x=xval,y=yval,yerr=err,color=cstr,marker='o',capsize=4)

# # a few more aesthetics
# plt.ylim(.4,.68)
# ticks = pd.np.linspace(.4,.68,29)
# labels = [ '{:.02f}'.format(x)[1:] if x in pd.np.linspace(.4,.65,6) else '' for x in ticks ]
# plt.yticks(ticks,labels)

# # add legend
# # for ax in g.axes.flat:
# #     box = ax.get_position()
# #     ax.set_position([box.x0,box.y0,box.width*0.7,box.height])
# handles = [ Patch(facecolor=ttype_palette[x],edgecolor=ttype_palette[x],label=x) for x in ttype_order ]
# legend = plt.legend(handles=handles,title='Trial type',frameon=False,loc='upper left',bbox_to_anchor=(1,1))
# plt.setp(legend.get_title(),fontsize='x-large')

# # a few moreeeee aesthetics
# plt.xlim(-1.5,2.5)
# plt.tight_layout()

# # save
# fname = os.path.join(res_dir,'within_del-ha-fig-combined.pdf')
# plt.savefig(fname,dpi=300)
# plt.close()



# ##--- collapsing across stay/switch trials in delay1 ---##
# ##
# ## Here just make separate plots for delay1/2.
# ##

# ## delay 1
# fig, ax = plt.subplots(figsize=(4,6))

# # TODO: dont get random samples like below, something more formal

# # violin plots of random distributions
# vtc_distrn = df.query("metric=='null'"
#               ).query("event=='delay1'"
#               ).query("mask=='vtc'"
#               ).accuracy.sample(n=1000,replace=False
#               ).values
# ips_distrn = df.query("metric=='null'"
#               ).query("event=='delay1'"
#               ).query("mask=='ips'"
#               ).accuracy.sample(n=1000,replace=False
#               ).values
# viols = ax.violinplot([vtc_distrn,ips_distrn],positions=[0,1],widths=[.125,.125],showextrema=False)
# plt.setp(viols['bodies'],facecolor='gray',edgecolor='white')

# # Take random 9 values from stay/switch trials
# vtc_y, vtc_sem = df.query("metric=='acc'"
#                   ).query("event=='delay1'"
#                   ).query("mask=='vtc'"
#                   ).accuracy.sample(n=9,replace=False
#                   ).aggregate(['mean','sem'])
# ips_y, ips_sem = df.query("metric=='acc'"
#                   ).query("event=='delay1'"
#                   ).query("mask=='ips'"
#                   ).accuracy.sample(n=9,replace=False
#                   ).aggregate(['mean','sem'])

# ax.errorbar(0,vtc_y,yerr=vtc_sem,color='red',marker='o',capsize=5,
#             markersize=10,elinewidth=3,capthick=3)
# ax.errorbar(1,ips_y,yerr=ips_sem,color='red',marker='o',capsize=5,
#             markersize=10,elinewidth=3,capthick=3)


# plt.xlim(-.5,1.5)
# plt.ylim(.4,.68)
# ticks = pd.np.linspace(.4,.68,29)
# labels = [ '{:.02f}'.format(x)[1:] if x in pd.np.linspace(.4,.65,6) else '' for x in ticks ]
# plt.yticks(ticks,labels)
# plt.xticks([0,1],['VTC','IPS'])
# plt.ylabel('Classifier accuracy')
# ax.axhline(.5,linestyle='--',linewidth=.5,color='black',zorder=0)

# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.tight_layout()

# fname = os.path.join(res_dir,'within_del-ha-fig-combined_del1.eps')
# plt.savefig(fname,dpi=300)
# plt.close()


# ## delay 2
# fig, ax = plt.subplots(figsize=(4,6))



# # permutation distributions
# distrns = df.query("metric=='null'"
#            ).query("event=='delay2'")


# # get mean and sem for all conditions
# plot_df = df.query("metric=='acc'"
#            ).query("event=='delay2'"
#            ).groupby(['mask','ttype']
#            ).aggregate({'accuracy':['mean','sem']})

# # point plots of accuracy and distrns at same time
# for indx, row in plot_df.iterrows():
#     mask, ttype = indx
#     cstr = ttype_palette[ttype]
#     xval = mask_order.index(mask)
#     xval += .12 * (-1 if ttype=='stay' else 1)
#     yval = row.accuracy['mean']
#     err  = row.accuracy['sem']

#     dstn = distrns.loc[(distrns['mask']==mask)&(distrns['ttype']==ttype),'accuracy'].values
#     viol = ax.violinplot(dstn,positions=[xval],widths=.12,showextrema=False)
#     plt.setp(viol['bodies'],facecolor='gray',edgecolor='white')
#     ax.errorbar(x=xval,y=yval,yerr=err,color=cstr,marker='o',capsize=5,
#                 markersize=10,elinewidth=3,capthick=3)



# # a few more aesthetics
# plt.ylim(.4,.68)
# ticks = pd.np.linspace(.4,.68,29)
# labels = [ '{:.02f}'.format(x)[1:] if x in pd.np.linspace(.4,.65,6) else '' for x in ticks ]
# plt.yticks(ticks,labels)
# plt.xticks([0,1],['VTC','IPS'])
# plt.ylabel('Classifier accuracy')
# ax.axhline(.5,linestyle='--',linewidth=.5,color='black',zorder=0)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.xlim(-.5,1.5)
# plt.tight_layout()

# fname = os.path.join(res_dir,'within_del-ha-fig-combined_del2.eps')
# plt.savefig(fname,dpi=300)
# plt.close()

