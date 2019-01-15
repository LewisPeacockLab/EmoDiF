function [] = imdif_run_net_importance_sampler()


addpath /Applications/MATLAB_R2014a.app/toolbox/stats/stats %randsample was not working correctly - was using SPM12's version - so rather than replace SPM12's randsample, i'm recalling the toolbox for matlab.
data_dir = '~/Dropbox (LewPeaLab)/STUDY/IMDiF/Results/P-CIT';

%%% VARIABLES THAT CHANGE %%%%
analysis_id = '160226a';
group = 'full_N20';
%trial = 'F';
%time = '10to18sec';
trstart = '10';
trlens = '6';
run_num = 2; %can be tarsub, face or scene (in single quotes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [] = REPREF_RUN_IMPORTANCE_SAMPLER()
% 
% Purpose
% 
% This scripts sets up the data matrix (number of samples x 6 columns) and the 'analysis_settings' struct with algorithm parameters
% 
% Input
%
% --None
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
analysis_settings.target_dir = '~/Dropbox (LewPeaLab)/STUDY/IMDiF/Results/P-CIT/results';
%analysis_settings.target_dir = fullfile(analysis_settings.results_dir,analysis_id); % % analysis_settings.analysis_id = '090909asim'; % analysis_id: specifies the target directory into which the output .mat will be located; if empty then the target directory is the timestamp of the form - YYYY-MM-DD-HH-MM
analysis_settings.em_iterations = 20; % Number of expectation maximization iterations
analysis_settings.particles = 100000; % Number of particles to be used in the importance sampling algorithm
%analysis_settings.curve = 'horz_indpnt'; % Name of the family of curves to be used. Refer to the family_of_curves.m file for more info
analysis_settings.distribution = 'normal'; % Name of the distribution (and the default canonical link function which maps the predictor variable to the dependent variable)
%analysis_settings.function_type = 'logreg'; % Name of the function that maps the independent variable to the dependent variable. For logistic regression it looks like ...
% f(z) = 1 / (1 + exp -(beta_0 + beta_1 x independent variable)), here f(z) gives you the probability of the dependent variable
analysis_settings.beta_0 = 0; % Initializing beta_0 for logistic regression function
analysis_settings.beta_1 = 1; % Initializing beta_1 for logistic regression function
analysis_settings.tau = 0.05; % Specifies the radius to sample curves in the curve space
analysis_settings.category = []; % Specifies if the analyses will need to run on a specific category. For instance '2' will cause the analyses to be run only on the second category; '-1' will run the analyses on all categories
analysis_settings.drop_outliers = 0; % specifies how many std dev away from group mean will the independent variable outliers need to be dropped
analysis_settings.zscore_within_subjects = false; % if TRUE, the independednt variables will be zscored within each suibject
analysis_settings.resolution = 4; % Denotes the resolution in which the data will be processed
analysis_settings.particle_chunks = 2; % Denotes the number of chunks you plan to partition the samples x particles matrix. An example chunk size will be 2 for a 3000 x 50,000 matrix
analysis_settings.repetition = 1;
analysis_settings.trs_of_interest = sprintf('from_TR%s_for_%sTRs',trstart,trlens);
analysis_settings.runs = run_num;
analysis_settings.subj_group = group;

%%%%%%%%%%%%%%%%%%%%%
% Reading in raw data
%%%%%%%%%%%%%%%%%%%%%
% % results_dir = fullfile(pwd, 'results');


start_dir = pwd;


for run = 1:run_num
  analysis_settings.runs = run;
  %analysis_settings.analysis_id = sprintf('%s_%s_%s_%s_%s_%s_r%d',group,trial,time,mvpa,recog,datestamp,run);
  %analysis_settings.analysis_id = sprintf('%s_%s_%s_%s_%s_%s_r%d',group,trial,time,datestamp,run);
  data_file = fullfile(data_dir, sprintf('imdif_net_pcit_data_windowstart%s_windowlen%s.mat',trstart,trlens));
  load(data_file);
  raw_data = data;
  
  imdif_importance_sampler(raw_data, analysis_settings);
  
end

cd(start_dir);