import os
import glob

SUBJECTS = [ 's2{:02d}'.format(i+1) for i in range(18) ]
SUBJECTS.remove('s211') # completed only 3 task runs
SUBJECTS.remove('s218') # did task all wrong


N_TRIAL_TRs = dict(map=10,early=10,stay=18,switch=18) # including ITI
N_ITI_TRs = dict(map=4,memory=3)
DELAY1_TRs = [ 4, 5]
BLAST_TRs  = [ 8, 9] # but ONLY on long trials
DELAY2_TRs = [12,13]

home = os.path.expanduser('~')
data_dir = os.path.join(home,'DBw/STUDY/RetroBlast/data/{:s}')
mask_template = os.path.join(data_dir,'ref/{:s}.nii.gz')
attr_template = os.path.join(data_dir,'behavior/{:s}_e01_data.csv')
run_templates = lambda subj,task: sorted(glob.glob(os.path.join(data_dir.format(subj),'func/{:s}*mc.nii.gz'.format(task))))


def load_data(mask):
    '''
    Load fMRI data and attributes from behavioral file
    into a single PyMVPA fMRI dataset for each subj.
    Export dictionaries of data for memory task
    and localizer task.
    Used for each classification analysis file.
    ARGS
    ----
    mask : 'vtc' or 'ips'
    
    '''
    from collections import OrderedDict
    from tqdm import tqdm
    import pandas as pd
    import mvpa2.suite as mvpa2


    mem_ds_dict = OrderedDict([])
    map_ds_dict = OrderedDict([])

    for subj in tqdm(SUBJECTS,desc='loading fmri data'):
        mask_fname = mask_template.format(subj,mask)

        # import memory data
        attr_fname = attr_template.format(subj,'memory')
        run_fnames = run_templates(subj,'memory')
        attr = pd.read_csv(attr_fname)
        # extracting attributes for tst data is tricky bc of diff trial lengths
        keep_cols = ['run','trial','amiCat','umiCat','trialType','acc','subj']
        attr_df = pd.concat( # same as above for each row/trial
            [ pd.DataFrame(data={ col: val for col, val in row[keep_cols].items() },
                           index=pd.Index(range(N_TRIAL_TRs[row['trialType']]),name='repetitions'))
              for _, row in attr.iterrows() ]
            ).rename(columns={'run':'chunks'}
            ).dropna( # some subjs didn't finish all runs
            ).reset_index()
        # # add some useful attributes: targets (cued item at each timepoint) and events
        # attr_df['targets'] = attr_df.apply(axis=1,func=lambda row: 
        #                 row['umiCat'] if 
        #                     (row['trialType']=='switch' and row['repetitions']>CUED_CUTOFF) 
        #                     else row['amiCat'])
        attr_df['events'] = attr_df.apply(axis=1,func=lambda row:
                                'delay1' if (row['repetitions'] in DELAY1_TRs
                                    ) else
                                'delay2' if (row['repetitions'] in DELAY2_TRs
                                    ) else
                                'blast'  if ((row['repetitions'] in BLAST_TRs
                                    ) and (row['trialType'] != 'early')
                                    ) else
                                'UNDECLARED')

        def define_cond(row):
            # Want to collapse trialType (early/stay/switch)
            # and events (delay1/delay2) into 4 "conditions"
            # (d1/blast/d2-stay/d2-switch)
            if row['events'] == 'delay1':
                cond = 'd1'
            elif row['events'] == 'blast':
                cond = 'blast'
            elif row['events'] == 'delay2':
                if row['trialType'] == 'stay':
                    cond = 'd2-stay'
                elif row['trialType'] == 'switch':
                    cond = 'd2-switch'
            else:
                cond = 'UNDECLARED'
            return cond
        attr_df['condition'] = attr_df.apply(define_cond,axis=1)

        # load in fMRI data
        ds = mvpa2.vstack([ mvpa2.fmri_dataset(samples=fn,mask=mask_fname) for fn in run_fnames ],a=True)
        # add attribute dataframe columns to fmri datasets
        for col, series in attr_df.items():
            ds.sa[col] = mvpa2.ArrayCollectable(name=col,doc=None,length=None,
                                                value=series.values if series.dtype!='O' else series.values.astype(str))
                                                      # change object dtypes to strings, for pymvpa
        mem_ds_dict[subj] = ds


        attr_fname = attr_template.format(subj,'mapping')
        run_fnames = run_templates(subj,'mapping')
        attr = pd.read_csv(attr_fname)
        ds = mvpa2.vstack([ mvpa2.fmri_dataset(samples=fn,mask=mask_fname) for fn in run_fnames ],a=True)
        # run_ds_list = []
        # for rnum, fn in enumerate(run_fnames):
        #     run_ds = mvpa2.fmri_dataset(samples=fn,mask=mask_fname)
        #     mvpa2.zscore(ds)
        #     for col, series in attr[attr.chunks==rnum].items():
        #         run_ds.sa[col] = mvpa2.ArrayCollectable(name=col,doc=None,length=None,value=series.values if series.dtype!='O' else series.values.astype(str))
        #     # # SHIFT LABELS bc of bold lag
        #     # run_ds.targets[3:] = run_ds.targets[:-3]
        #     # # SHIFT DATA
        #     # run_ds = run_ds[3:]
        #     run_ds_list.append(run_ds)
        # ds_list.append(mvpa2.vstack(run_ds_list,a=True))
        keep_cols = ['run','miniblock','category','subj']
        attr_df = pd.concat( # col, series dict gets unique values for each miniblock/trial
            [ pd.DataFrame(data={ col: series.unique()[0] for col, series in mb_df[keep_cols].items() },
                           index=pd.Index(range(N_TRIAL_TRs['map']),name='repetitions'))
              for _, mb_df in attr.groupby(['run','miniblock']) ]
            ).rename(columns={'run':'chunks','category':'targets','miniblock':'miniblocks'}
            ) # DO NOT reset index yet, bc use it to insert iti labels next line
        # place ITIs at end of each trial
        attr_df.loc[range(N_TRIAL_TRs['map']-N_ITI_TRs['map'],N_TRIAL_TRs['map']),'targets'] = 'iti'
        # now reset index
        attr_df.reset_index(inplace=True)
        for col, series in attr_df.items():
            ds.sa[col] = mvpa2.ArrayCollectable(name=col,doc=None,length=None,
                                                value=series.values if series.dtype!='O' else series.values.astype(str))
        map_ds_dict[subj] = ds


    return map_ds_dict, mem_ds_dict
