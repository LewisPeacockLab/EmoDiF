function [] = imdif_batch_importance_sampler()

% [] = REPREF_BATCH_IMPORTANCE_SAMPLER()
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
% % analysis_settings.analysis_id = '090909asim'; % analysis_id: specifies the target directory into which the output .mat will be located; if empty then the target directory is the timestamp of the form - YYYY-MM-DD-HH-MM
analysis_settings.em_iterations = 10; % Number of expectation maximization iterations
analysis_settings.particles = 100000; % Number of particles to be used in the importance sampling algorithm
analysis_settings.curve = 'horz_indpnt'; % Name of the family of curves to be used. Refer to the family_of_curves.m file for more info
analysis_settings.function_type = 'logreg'; % Name of the function that maps the independent variable to the dependent variable. For logistic regression it looks like ...
% f(z) = 1 / (1 + exp -(beta_0 + beta_1 x independent variable)), here f(z) gives you the probability of the dependent variable
analysis_settings.beta_0 = 0; % Initializing beta_0 for logistic regression function
analysis_settings.beta_1 = 1; % Initializing beta_1 for logistic regression function
analysis_settings.tau = 0.05; % Specifies the radius to sample curves in the curve space
analysis_settings.category = -1; % Specifies if the analyses will need to run on a specific category. For instance '2' will cause the analyses to be run only on the second category; '-1' will run the analyses on all categories
analysis_settings.drop_outliers = 0; % specifies how many std dev away from group mean will the independent variable outliers need to be dropped
analysis_settings.zscore_within_subjects = false; % if TRUE, the independednt variables will be zscored within each suibject
analysis_settings.resolution = 4; % Denotes the resolution in which the data will be processed
analysis_settings.particle_chunks = 2; % Denotes the number of chunks you plan to partition the samples x particles matrix. An example chunk size will be 2 for a 3000 x 50,000 matrix
analysis_settings.repetition = 1;

%%%%%%%%%%%%%%%%%%%%%
% Reading in raw data
%%%%%%%%%%%%%%%%%%%%%
% % results_dir = fullfile(pwd, 'results');
data_dir = fullfile(pwd, 'data');

group = 'full_N21';
trial = 'switch';
time = '12to18sec';
mvpa = 'sceneminface';
recog = 'remember';

start_dir = pwd;
cd /jukebox/norman/jalewpea/fmri/repref/curvefit;

for run = 1:10
  
  analysis = sprintf('%s_%s_%s_%s_%s',group,trial,time,mvpa,recog);
  analysis_settings.analysis_id = sprintf('%s_%s_r%d',analysis,datestr(now,'yymmdd'),run);
  data_file = sprintf('%s/imdif_%s.mat',data_dir,analysis);
  
  load(data_file);
  raw_data = data;
  
  importance_sampler(raw_data, analysis_settings);
  
end

cd(start_dir);