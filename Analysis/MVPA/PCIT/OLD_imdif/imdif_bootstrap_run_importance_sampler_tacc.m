function [] = imdif_bootstrap_run_importance_sampler_tacc(bootstrap_num)
  
  %----------------------------------------------------------------------
  % This is the bootstrap wrapper script - this calls the
  % preprocessing_setup script and then bootstraps this the number of
  % bootstrap iterations (bootstrap_num). 
  % * bootstrap_id specifies name of bootstrap 
  % (e.g., '200_130303_0523_13_598')
  
  % THIS REPLACES IMDIF_RUN_IMPORTANCE_SAMPLER.M!!!!
  % tracy wang 2/25/16
  %
  %%----------------------------------------------------------------------
  
addpath /Applications/MATLAB_R2014a.app/toolbox/stats/stats %randsample was not working correctly - was using SPM12's version - so rather than replace SPM12's randsample, i'm recalling the toolbox for matlab.
data_dir = '~/Dropbox (LewPeaLab)/STUDY/IMDiF/Results/P-CIT';

%%% VARIABLES THAT CHANGE %%%%
analysis_id = '171027a';
group = 'full_N19';
%trial = 'F';
%time = '10to18sec';
trstart = '1';
trlens = '9';
run_num = 1;
target = 'target-nontarget'; %'preFR'
type = 'net';
condition = 'all';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Purpose
% 
% This scripts sets up the data matrix (number of samples x 6 columns) and
% the 'analysis_settings' struct with algorithm parameters, then it
% resamples with replacement data from 20 subjects and runs
% importance_sampler.m bootstrap_num (e.g. 500) number of times. 
% 
% Input
%
% --None - but modify the %%% variables that change section above %%%
% 
% Output
%
% --None
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Populating the analysis_settings struct with algorithm settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analysis_settings = struct(); % Creating a struct
analysis_settings.analysis_id = analysis_id; % analysis_id: specifies the target directory into which the output .mat will be located; if empty then the target directory is the timestamp of the form - YYYY-MM-DD-HH-MM
analysis_settings.target_dir = '~/Dropbox (LewPeaLab)/STUDY/IMDiF/Results/P-CIT/bootstrap_results';
%analysis_settings.target_dir = fullfile(analysis_settings.results_dir,analysis_id); % % analysis_settings.analysis_id = '090909asim'; % analysis_id: specifies the target directory into which the output .mat will be located; if empty then the target directory is the timestamp of the form - YYYY-MM-DD-HH-MM
analysis_settings.em_iterations = 20; % Number of expectation maximization iterations
analysis_settings.particles = 100000; % Number of particles to be used in the importance sampling algorithm
analysis_settings.curve_type = 'horz_indpnt';  % Name of the family of curves to be used. Refer to the family_of_curves.m file for more info
analysis_settings.distribution = 'normal'; % Name of the distribution (and the default canonical link function which maps the predictor variable to the dependent variable)
analysis_settings.dist_specific_params = struct(); % For normal distribution the additional parameter is sigma. We pass in sigma here.
analysis_settings.dist_specific_params.sigma = 1;

%analysis_settings.function_type = 'logreg'; % Name of the function that maps the independent variable to the dependent variable. For logistic regression it looks like ...
% f(z) = 1 / (1 + exp -(beta_0 + beta_1 x independent variable)), here f(z) gives you the probability of the dependent variable
analysis_settings.beta_0 = 0; % Initializing beta_0 for logistic regression function
analysis_settings.beta_1 = 1; % Initializing beta_1 for logistic regression function
analysis_settings.tau = 0.05; % Specifies the radius to sample curves in the curve space
analysis_settings.category = []; % Specifies if the analyses will need to run on a specific category. For instance '2' will cause the analyses to be run only on the second category; '-1' will run the analyses on all categories

analysis_settings.drop_outliers = 0; % note: default is 3 -- specifies how many std dev away from group mean will the independent variable outliers need to be dropped
analysis_settings.zscore_within_subjects = false; % if TRUE, the independednt variables will be zscored within each suibject

analysis_settings.data_matrix_columns = struct();
analysis_settings.data_matrix_columns.subject_id = 1;
analysis_settings.data_matrix_columns.trials = 2;
analysis_settings.data_matrix_columns.category = 3;
analysis_settings.data_matrix_columns.predictor_var = 4;
analysis_settings.data_matrix_columns.dependent_var = 5;
analysis_settings.data_matrix_columns.net_effect_clusters = 6;

analysis_settings.resolution = 4; % Denotes the resolution in which the data will be processed
analysis_settings.particle_chunks = 2; % Denotes the number of chunks you plan to partition the samples x particles matrix. An example chunk size will be 2 for a 3000 x 50,000 matrix

analysis_settings.repetition = 1;
analysis_settings.trs_of_interest = sprintf('from_TR%s_for_%sTRs',trstart,trlens);

analysis_settings.bootstrap = true; % indicates that this run is a bootstrap run
analysis_settings.bootstrap_num = bootstrap_num;
% reset bootstrap sample number before each boostrap iteration below. 
%analysis_settings.bootstrap_run = -1; % will need to specify a bootstrap sample number. This will need to be unique for each sample

analysis_settings.scramble = false; % indicates that this run is a scramble run
analysis_settings.scramble_run = -1; % will need to specify a scramble sample number. This will need to be unique for each sample
analysis_settings.scramble_style = -1; % choosing the appropriate scramble option from three options below
if analysis_settings.scramble_style > 0
	switch analysis_settings.scramble_style
	case 1, analysis_settings.scramble_style = 'within_subjects_within_categories';
	case 2, analysis_settings.scramble_style = 'within_subjects_across_categories';
	case 3, analysis_settings.scramble_style = 'across_subjects_across_categories';
	otherwise, error('Invalid scramble style!');
	end
end


analysis_settings.runs = run_num;
analysis_settings.subj_group = group;
analysis_settings.target = target;
analysis_settings.type = type;
analysis_settings.condition = condition;

%%%%%%%%%%%%%%%%%%%%%
% Reading in raw data
%%%%%%%%%%%%%%%%%%%%%
% % results_dir = fullfile(pwd, 'results');

start_dir = pwd;

  data_file = fullfile(data_dir, sprintf('/data/imdif_pcit_data_9TRS_%s_%s_%s.mat',condition,type,target));
  load(data_file);
  original_data = data;

for bs_iter = 1:bootstrap_num
    analysis_settings.bootstrap_run = bs_iter;
% for run = 1:run_num
%   analysis_settings.runs = run;

  %preprocess data into resample with replacement data 
  preprocessing_setup(original_data, analysis_settings)
  % input = all subjects data, output = resampled with replacement data.
  
  bs_data = data;
  imdif_importance_sampler(bs_data, analysis_settings);
  
% end
end

cd(start_dir);
  