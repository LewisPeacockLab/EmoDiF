import os
import glob
import tqdm
import argparse
import pandas as pd
import mvpa2.suite as mvpa2

from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import permutation_test_score
from sklearn.model_selection import LeaveOneGroupOut

from load_data import load_data


parser = argparse.ArgumentParser()
parser.add_argument('--mask',default='vtc',choices=['vtc','ips'],type=str)
parser.add_argument('--nperms',default=2,type=int)
args = parser.parse_args()

MASK = args.mask
N_PERMUTATIONS = args.nperms

home = os.path.expanduser('~')
res_dir  = os.path.join(home,'DBw/STUDY/RetroBlast/results')
if MASK == 'ips':
    res_dir = os.path.join(res_dir,MASK)

K_FEATS = 1000 if MASK=='vtc' else 500

featsel = SelectKBest(f_classif,k=K_FEATS)
clf = LogisticRegression(penalty='l2',multi_class='ovr',solver='liblinear')


#################
##  LOAD DATA  ##
#################

map_ds_dict, mem_ds_dict = load_data(MASK)

# preprocess
for d in [mem_ds_dict,map_ds_dict]:
    for ds in d.values():
        mvpa2.remove_invariant_features(ds)
        mvpa2.poly_detrend(ds,polyord=1,chunks_attr='chunks')
        mvpa2.zscore(ds,chunks_attr='chunks')



##############################################################
##  build and convert to common space using hyperalignment  ##
##############################################################

# select features based on localizer data
fsel_masks = [ featsel.fit(ds.samples,ds.targets).get_support() for ds in map_ds_dict.values() ]
# apply feature selection to all data (localizer and memory)
fs_mapds_list = [ ds[:,mask] for ds, mask in zip(map_ds_dict.values(),fsel_masks) ]
fs_memds_list = [ ds[:,mask] for ds, mask in zip(mem_ds_dict.values(),fsel_masks) ]

hyper = mvpa2.Hyperalignment()
hypmaps = hyper(fs_mapds_list) # returns list

# apply the hyperalignment maps to feature-selected test data
memds_hyper_list = [ ha.forward(ds) for ha, ds in zip(hypmaps,fs_memds_list)]

# stack all datasets and zscore, because now all in common space
ds_hyper = mvpa2.vstack(memds_hyper_list)
mvpa2.zscore(ds_hyper,chunks_attr='subj')



#######################################
##  use sklearn to perform decoding  ##
#######################################

# get the num of groups for leave-one-subj-out cross-validation
n_groups = len(pd.np.unique(ds_hyper.sa.subj))

cv_df_list = []   # holds cross-validation accuracies for each condition
null_df_list = [] # holds permutation distribution accuracies for each condition
pval_df_list = [] # holds mean accuracy and p-values for each condition

# choose the "target" sample attributes
# doesn't matter if ami is target on switch trials bc only 2 classes
ds_hyper.sa['targets'] = ds_hyper.sa['amiCat']

for cond in tqdm.tqdm(['d1','blast','d2-stay','d2-switch'],desc='HA-clf'):
    # extract relevant subset of dataset
    cond_ds = ds_hyper[ds_hyper.sa.condition==cond]

    # extract info from PyMVPA dataset for sklearn
    XXxx = cond_ds.samples
    YYyy = cond_ds.targets
    GGgg = cond_ds.sa.subj

    # run cross-val-score to get the array of accuracies for each fold
    # and permutation-test-score to get the shuffled distributions
    scores = cross_val_score(clf,XXxx,YYyy,GGgg,
                             cv=LeaveOneGroupOut())
    mu, perms, p = permutation_test_score(clf,XXxx,YYyy,GGgg,
                                          cv=LeaveOneGroupOut(),
                                          n_permutations=N_PERMUTATIONS,
                                          n_jobs=-1,verbose=2)

    # save into lists of dataframes
    cv_df = pd.DataFrame(scores,columns=['accuracy'],
                index=pd.MultiIndex.from_product([[MASK],[cond],range(n_groups)],
                names=['mask','cond','fold']))
    cv_df_list.append(cv_df)

    null_df = pd.DataFrame(perms,columns=['accuracy'],
                index=pd.MultiIndex.from_product([[MASK],[cond],range(N_PERMUTATIONS)],
                names=['mask','cond','fold']))
    null_df_list.append(null_df)

    pval_df = pd.DataFrame(pd.np.array([mu,p]).reshape(1,-1),columns=['mean_acc','pval'],
                index=pd.MultiIndex.from_product([[MASK],[cond]],
                names=['mask','cond']))
    pval_df_list.append(pval_df)


# save out results by stacking all lists into single dataframes and exporting
cv_out = pd.concat(cv_df_list)
null_out = pd.concat(null_df_list)
pval_out = pd.concat(pval_df_list)

fname = os.path.join(res_dir,'within_del-ha-acc.csv')
cv_out.to_csv(fname)
fname = os.path.join(res_dir,'within_del-ha-null.csv')
null_out.to_csv(fname)
fname = os.path.join(res_dir,'within_del-ha-pval.csv')
pval_out.to_csv(fname)