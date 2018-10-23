% Princeton MVPA Searchlight analysis tutorial workthrough
%
% Judy Chiu
% Oct. 2018

% Instructions
% https://github.com/PrincetonUniversity/princeton-mvpa-toolbox/wiki/TutorialSpheres

%% Set up subject structure

% Run function "init_sub" (within the tutorial_easy.m)
% make sure the data "working_set" and all Princeton Toolbox Subfolders are
% added to Matlab Path

[subj results]=tutorial_easy();

% remove previously established ANOVA models and results
subj = remove_group(subj,'pattern','epi_z_anova');
subj = remove_group(subj,'mask','epi_z_thresh0.05');
subj = remove_group(subj,'selector','runs_xval');



%%  Create cross-validation "selectors" --

%  which runs to include as training, test, and generalization set
%  Training data =1, test data= 2, searchlight validation data =3;

% previous selector  ( 1 by 1210, each number signify run in the exp for that TR)
runs=get_mat(subj,'selector','runs');

% Find out how many runs in this study
nRuns = max(runs);

% Find out how many total TR's
nTimepoints = length(runs);

% preallocate matrix that is (number of iterations of tranining, which is the # of turns,  by total TR's)
runs_xval_sl = ones(nRuns,nTimepoints);  % Default index is 1.

% Specify test (=2) indices
% (test runs, in this tutorial arbitrarilly, correspond to the nth iteration). hold 1st, 2nd, 3rd.... to Nth runs out.
for r=1:nRuns
    
    % find through each iteration where the positions of each run is
    cur_final_test_run = find(runs==r);
    
    % populate into runs_xval_sl matrix and set those to index 2 
    runs_xval_sl(r, cur_final_test_run) = 2;
    
end

% visualize which runs to leave out
imagesc(runs_xval_sl)
set(gca,'CLim',[1 3])
colorbar

% lable searchlight verification data (=3) 
% generalization runs, arbitrarily for this example, Total-nth iteration +1 
for r=1:nRuns
    
    % set searchight validation data as the current iteration -1 runs;
    cur_searchlight_gen_run=find(runs==nRuns-r+1);
    
    % populate into runs_xval_sl matrix and set those to index 3
    runs_xval_sl(r,cur_searchlight_gen_run)=3;
end

imagesc(runs_xval_sl)
colorbar

%%  Make the selector subfield in SUBJ
% using the selection matrices made above

% Feed in data row by row for FEATURE_SELECT.m and CROSS_VALIDATION.m
for r=1:nRuns
    
    cur_name = sprintf('runs_xval_sl_%i',r);  % name='runs_xval_sl_1,2,3,...etc...
    
    % initset_object = init_object, set_mat, set_objectfield methods all in
    % one go.
    subj = initset_object(subj, 'selector', cur_name, ...
        runs_xval_sl(r,:), ...
        'group_name', 'runs_xval_sl' );
end


%% Define Searchlight extent -- Adjacecy Matrix 
% creates a "adj_sphere" field within subj struct. Input for searchlight
% below.

subj.adj_sphere = create_adj_list(subj,'VT_category-selective', 'radius', 1); 
% specfies a radius of 1, for voxels within the mask VT_category-selective.


%% Feature Selection using searchlight spheres 
% 
% rather than using ANOVA for feature selection, use a searchlight
% adjacency matrix


%%%%% Set arguments for feature-selection related functions %%%%%%% 

% Specifies the classifier algorithm
class_args.train_funct_name = 'train_gnb';  % GNB -- Gaussian Naive Bayes; bp= back propagation
class_args.test_funct_name = 'test_gnb';

% Specifications for voxel selection
% Arguments for STATMAP_CLASSIFY.m (see below)
scratch.class_args = class_args;
scratch.perfmet_funct = 'perfmet_maxclass';   % chosen function to determine performance
scratch.perfmet_args = struct([]);

% Arguments for STATMAP_SEARCHLIGH.M  (see below, within feature selection)
statmap_srch_arg.adj_list = subj.adj_sphere;  % Adjacency matrix made above.
statmap_srch_arg.obj_funct = 'statmap_classify';  % specify function to use
statmap_srch_arg.scratch = scratch;  % specifies options for 'statmap_classify' above.


% run feature selection
subj = feature_select( ...
    subj, ...
    'epi_z', ... % data
    'conds', ... % binary regs (for GNB)
    'runs_xval_sl', ... % selector
    'statmap_funct','statmap_searchlight', ... % function
    'statmap_arg',statmap_srch_arg, ...
    'new_map_patname','epi_z_srch', ...
    'thresh',[]);


%%  Create masks from the searchligh feature selection results

% Ranked voxels from the searchlight
subj = create_sorted_mask( ...
    subj,'epi_z_srch', ...   % searchlight results patterns
    'epi_z_srch_200',200, ... % output mask, top #200 voxels
    'descending',true);       % sort by descending to get top/highest 200 voxels

%% Visualization 
% https://github.com/PrincetonUniversity/princeton-mvpa-toolbox/wiki/TutorialVisualization


% use in-house function to view the results

% load a whole-brain image as a mask (1,0)
subj=load_afni_mask(subj,'wholebrain','wholebrain+orig');  % Loads all non-zeros values as 1 (makes this a mask)


% downsample the brain mask so resolution matches our functional
% calculations. Do this in SHELL/bash





