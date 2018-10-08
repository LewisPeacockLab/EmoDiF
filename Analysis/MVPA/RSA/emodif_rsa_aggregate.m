function [results] = emodif_rsa_aggregate(type, preview_shift, TRstart, TRlength, maskName)
%emodif_rsa_aggregate(4, 2, 5, 3, 'scene-all_005')
%subj_list - mat file with subject numbers and folder analysis dates (for
%identification of correct analysis)
%type - intraphase, or preview_dfencode 
  %1) intraphase localizer - by category
  %2) intraphase decoding - by trial
  %3) intraphase preview - by trial
  %4) interphase preview - decoding - by trial, by decoding
  % 5) interphase Localizer- DF and Preview - by trial 

%This applies fisher's z transform, concatenates arrays along the 3rd dimension 
% and then averages along the 3rd dimension to create one array of zscores
% then we reverse teh fisher's z transform to report average Rs. 

%loading all the plots
%load subject list on astoria

load(sprintf('/Users/tw24955/EmoDif/Analysis/MVPA/subj_list_%s.mat', maskName));

%subj_list

for i = 1:length(subj_list);
subjNum = (subj_list{i,1});
subj_intraphase_folder = (subj_list{i,2});
subj_preview_dfencode_folder = (subj_list{i,3});
subj_preview_DF_local_folder = (subj_list{i,4});

args.base_dir = '/Users/tw24955/emodif_data';

args.subjID = sprintf('emodif_%s',num2str(subjNum));
args.maskName = maskName;
args.script_dir = pwd;
args.previewshift = preview_shift;

args.localizer.trialnum = 24; % 30 for subjects 1-3, but already coded in ...
%within phase rsa script. 24 should apply to everyone else. reducing to 24 happens at concatenation level. comment out 1-3 
%if not using these first 4 subjects.
args.preview.trialnum = 60;
args.DFencode.trialnum = 60;



args.subj_dir = sprintf('%s/%s', args.base_dir, args.subjID);


if type == 1 | type == 2 | type == 3
    
    args.data_dir = sprintf('%s/results/rsa_results/intraphase/%s/%d/%s',args.subj_dir,args.maskName, args.previewshift, subj_intraphase_folder);
    
elseif type == 4
    
    args.data_dir = sprintf('%s/results/rsa_results/preview_dfencode/%s/%d/%s',args.subj_dir,args.maskName, args.previewshift, subj_preview_dfencode_folder);
elseif type == 5
    args.data_dir = sprintf('%s/results/rsa_results/preview_DF_local/%s/%d/%s',args.subj_dir,args.maskName, args.previewshift, subj_preview_DF_local_folder);
end


results.bysubject.name = sprintf('%s/emodif_%d_TR%dto%d_rsa_results.mat',args.data_dir, subjNum, TRstart, ((TRstart+TRlength)-1));
results.bysubject.data(i) = load(results.bysubject.name);
results.bysubject.names{i} = results.bysubject.name;


end

if type == 1 | type == 2 | type == 3
    
args.output_dir = sprintf('%s/aggregate_results/RSA_intraphase_results/%s', args.base_dir, maskName);
args.outfname = sprintf('%s/rsa_intraphase_aggregate_%s', args.output_dir, date);
mkdir(args.output_dir);

%loaded all the data into results structure
preview_stack = [];

for x = 1:length(subj_list)
    preview_stack = cat(3, preview_stack, results.bysubject.data(x).rsa.results.smatrix.corr_matrix_match_previewz);
end

DFencode_stack = [];

for x = 1:length(subj_list)
    DFencode_stack = cat(3, DFencode_stack, results.bysubject.data(x).rsa.results.smatrix.corr_matrix_match_DFencodez);
end

Localizer_stack = [];

%this always requires 101, 102 and 103 being in the localizer stack

% s101l = Replace(results.bysubject.data(1).rsa.results.smatrix.corr_matrix_match_Localizerz,inf,NaN);
% s102l = Replace(results.bysubject.data(2).rsa.results.smatrix.corr_matrix_match_Localizerz,inf,NaN);
% s103l = Replace(results.bysubject.data(3).rsa.results.smatrix.corr_matrix_match_Localizerz,inf,NaN);
% 
% cat_previous_101 = cat(3,s101l(:,25:30), s101l(:,19:24));
% cat_previous_101m = nanmean(cat_previous_101,3);
% 
% cat_previous_102 = cat(3,s102l(:,25:30), s102l(:,19:24));
% cat_previous_102m = nanmean(cat_previous_102,3);
% 
% cat_previous_103 = cat(3,s103l(:,25:30), s103l(:,19:24));
% cat_previous_103m = nanmean(cat_previous_103,3);

results.bysubject.data(1).rsa.results.smatrix.corr_matrix_match_Localizer_bycatz= ...
    results.bysubject.data(1).rsa.results.smatrix.corr_matrix_match_Localizer_bycatz(1:24,1:24);

results.bysubject.data(2).rsa.results.smatrix.corr_matrix_match_Localizer_bycatz= ...
    results.bysubject.data(2).rsa.results.smatrix.corr_matrix_match_Localizer_bycatz(1:24,1:24);

results.bysubject.data(3).rsa.results.smatrix.corr_matrix_match_Localizer_bycatz= ...
    results.bysubject.data(4).rsa.results.smatrix.corr_matrix_match_Localizer_bycatz(1:24,1:24);


for x = 1:length(subj_list)
    Localizer_stack = cat(3, Localizer_stack, results.bysubject.data(x).rsa.results.smatrix.corr_matrix_match_Localizer_bycatz);
end

%replace infinity with NaN
preview_stack = Replace(preview_stack,inf,NaN);
DFencode_stack = Replace(DFencode_stack,inf,NaN);
Localizer_stack = Replace(Localizer_stack,inf,NaN);

%take the mean of the stack
preview_stack_mz = nanmean(preview_stack,3);
DFencode_stack_mz = nanmean(DFencode_stack,3);
Localizer_stack_mz = nanmean(Localizer_stack,3);

%set results file
results.groupmean.preview_rsa_z = preview_stack_mz;
results.groupmean.DFencode_rsa_z = DFencode_stack_mz;
results.groupmean.Localizer_rsa_z = Localizer_stack_mz;

%replace NaN with 1
results.groupmean.preview_rsa_z = Replace(preview_stack_mz,NaN, 1);
results.groupmean.DFencode_rsa_z = Replace(DFencode_stack_mz,NaN, 1);
results.groupmean.Localizer_rsa_z = Replace(Localizer_stack_mz,NaN, 1);

%fisher z transform everything back to pearson's R

% FISHERINV Equation 
%r = (exp(2*z)-1)/(exp(2*z)+1);
%FISHER Equation
%z = 0.5*log((1+r)/(1-r));

for x = 1:args.localizer.trialnum %trial number
    
    for y = 1:args.localizer.trialnum %trial number
        
        results.groupmean.Localizer_rsa(x,y) = (exp(2*results.groupmean.Localizer_rsa_z(x, y))-1)/(exp(2*results.groupmean.Localizer_rsa_z(x, y))+1);
    end
end

for x = 1:args.DFencode.trialnum %trial number
    
    for y = 1:args.DFencode.trialnum %trial number
        
        results.groupmean.DFencode_rsa(x,y) = (exp(2*results.groupmean.DFencode_rsa_z(x, y))-1)/(exp(2*results.groupmean.DFencode_rsa_z(x, y))+1);
    end
end

for x = 1:args.preview.trialnum %trial number
    
    for y = 1:args.preview.trialnum %trial number
        
        results.groupmean.preview_rsa(x,y) = (exp(2*results.groupmean.preview_rsa_z(x, y))-1)/(exp(2*results.groupmean.preview_rsa_z(x, y))+1);
        
    end
end
                                

cd(args.output_dir)
results.parameters = args;
save(args.outfname,'results')

%plot figures Z and R

    Localizerz_fig = figure;
    set(Localizerz_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.Localizer_rsa_z); colormap('parula'); colorbar; 
    
    
    ylabel(sprintf('Localizer zmean patterns averaged over TR %d: TR %d', TRstart, TRlength), 'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('Localizer zmean patterns averaged over TR %d: TR %d', TRstart, TRlength),'FontSize',15,'FontWeight','bold');
    saveas(Localizerz_fig, sprintf('emodif_average_run_Localizerz_TRstart%d_mean%d', TRstart, TRlength),'png');
    
    Localizer_fig = figure;
    set(Localizer_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.Localizer_rsa); colormap('parula'); colorbar;
    
    
    ylabel(sprintf('Localizer mean patterns averaged over TR %d: TR %d', TRstart, TRlength), 'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('Localizer mean patterns averaged over TR %d: TR %d', TRstart, TRlength),'FontSize',15,'FontWeight','bold');
    saveas(Localizer_fig, sprintf('emodif_average_run_Localizer_TRstart%d_mean%d', TRstart, TRlength),'png');
    
    %DFencode
    DFencodez_fig = figure;
    set(DFencodez_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.DFencode_rsa_z); colormap('parula'); colorbar; 
    
    
    ylabel(sprintf('DFencode zmean patterns averaged over TR %d: TR %d', TRstart, TRlength), 'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('DFencode zmean patterns averaged over TR %d: TR %d', TRstart, TRlength),'FontSize',15,'FontWeight','bold');
    saveas(DFencodez_fig, sprintf('emodif_average_run_DFencodez_TRstart%d_mean%d', TRstart, TRlength),'png');
    
    DFencode_fig = figure;
    set(DFencode_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.DFencode_rsa); colormap('parula'); colorbar;
    
    
    ylabel(sprintf('DFencode mean patterns averaged over TR %d: TR %d', TRstart, TRlength), 'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('DFencode mean patterns averaged over TR %d: TR %d', TRstart, TRlength),'FontSize',15,'FontWeight','bold');
    saveas(DFencode_fig, sprintf('emodif_average_run_DFencode_TRstart%d_mean%d', TRstart, TRlength),'png');

    %preview
    previewz_fig = figure;
    set(previewz_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.preview_rsa_z); colormap('parula'); colorbar; 
    
    
    ylabel(sprintf('preview zmean patterns averaged over TR %d: TR %d', TRstart, TRlength), 'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('preview zmean patterns averaged over TR %d: TR %d', TRstart, TRlength),'FontSize',15,'FontWeight','bold');
    saveas(previewz_fig, sprintf('emodif_average_run_previewz_TRstart%d_mean%d', TRstart, TRlength),'png');
    
    preview_fig = figure;
    set(preview_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.preview_rsa); colormap('parula'); colorbar;
    
    
    ylabel(sprintf('preview mean patterns averaged over TR %d: TR %d', TRstart, TRlength), 'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('preview mean patterns averaged over TR %d: TR %d', TRstart, TRlength),'FontSize',15,'FontWeight','bold');
    saveas(preview_fig, sprintf('emodif_average_run_preview_TRstart%d_mean%d', TRstart, TRlength),'png');
    
elseif type == 4
    
args.output_dir = sprintf('%s/aggregate_results/RSA_preview_dfencode_results/%s', args.base_dir, maskName);
args.outfname = sprintf('%s/rsa_preview_DFencode_aggregate_%s', args.output_dir, date);
mkdir(args.output_dir);
    
    preview_DFencode_stack = [];
    preview_DFencode_stack_F = [];
    preview_DFencode_stack_R = [];
    preview_DFencode_stack_mz = [];
    preview_DFencode_stack_F_mz = [];
    preview_DFencode_stack_R_mz = [];
    
    for x = 1:length(subj_list)
        preview_DFencode_stack = cat(3, preview_DFencode_stack, results.bysubject.data(x).rsa.results.smatrix.corr_matrix_match_full);
        preview_DFencode_stack_F = cat(3, preview_DFencode_stack_F, results.bysubject.data(x).rsa.results.smatrix.corr_matrix_match_F);
        preview_DFencode_stack_R = cat(3, preview_DFencode_stack_R, results.bysubject.data(x).rsa.results.smatrix.corr_matrix_match_R);
        preview_DFencode_stack_mz = cat(3, preview_DFencode_stack_mz , results.bysubject.data(x).rsa.results.smatrix.corr_matrix_match_fullz);
        preview_DFencode_stack_F_mz  = cat(3, preview_DFencode_stack_F_mz , results.bysubject.data(x).rsa.results.smatrix.corr_matrix_match_Fz);
        preview_DFencode_stack_R_mz  = cat(3, preview_DFencode_stack_R_mz , results.bysubject.data(x).rsa.results.smatrix.corr_matrix_match_Rz);
    end
    
    preview_DFencode_stack_m = mean(preview_DFencode_stack_mz,1);
    preview_DFencode_stack_m2 = mean(preview_DFencode_stack_m,2);
    preview_DFencode_stack_summary = reshape(preview_DFencode_stack_m2,24,1);
    
        preview_DFencode_stack_Fm = mean(preview_DFencode_stack_F_mz,1);
    preview_DFencode_stack_Fm2 = mean(preview_DFencode_stack_Fm,2);
    preview_DFencode_stack_Fsummary = reshape(preview_DFencode_stack_Fm2,24,1);

    
    preview_DFencode_stack_Rm = mean(preview_DFencode_stack_R_mz,1);
    preview_DFencode_stack_Rm2 = mean(preview_DFencode_stack_Rm,2);
    preview_DFencode_stack_Rsummary = reshape(preview_DFencode_stack_Rm2,24,1);
    
     for x = 1:length(subj_list)
        preview_DFencode_summary_1 = mean(results.bysubject.data(x).rsa.results.smatrix.corr_matrix_match_full,1);
        preview_DFencode_stack_F = cat(3, preview_DFencode_stack_F, results.bysubject.data(x).rsa.results.smatrix.corr_matrix_match_F);
        preview_DFencode_stack_R = cat(3, preview_DFencode_stack_R, results.bysubject.data(x).rsa.results.smatrix.corr_matrix_match_R);
        preview_DFencode_stack_mz = cat(3, preview_DFencode_stack_mz , results.bysubject.data(x).rsa.results.smatrix.corr_matrix_match_fullz);
        preview_DFencode_stack_F_mz  = cat(3, preview_DFencode_stack_F_mz , results.bysubject.data(x).rsa.results.smatrix.corr_matrix_match_Fz);
        preview_DFencode_stack_R_mz  = cat(3, preview_DFencode_stack_R_mz , results.bysubject.data(x).rsa.results.smatrix.corr_matrix_match_Rz);
    end
   

%take the mean of the z-scored stack
preview_DFencode_mz = nanmean(preview_DFencode_stack_mz,3);
preview_DFencode_F_mz  = nanmean(preview_DFencode_stack_F_mz ,3);
preview_DFencode_R_mz = nanmean(preview_DFencode_stack_R_mz,3);

%set results file
results.groupmean.preview_DFencode_mz= preview_DFencode_mz;
results.groupmean.preview_DFencode_F_mz  = preview_DFencode_F_mz;
results.groupmean.preview_DFencode_R_mz = preview_DFencode_R_mz;

%fisher z transform everything back to pearson's R

% FISHERINV Equation 
%r = (exp(2*z)-1)/(exp(2*z)+1);
%FISHER Equation
%z = 0.5*log((1+r)/(1-r));


for x = 1:args.DFencode.trialnum %trial number
    
    for y = 1:args.DFencode.trialnum %trial number
        
        results.groupmean.preview_DFencode(x,y) = (exp(2*results.groupmean.preview_DFencode_mz(x, y))-1)/(exp(2*results.groupmean.preview_DFencode_mz(x, y))+1);
    end
end

for x = 1:(args.DFencode.trialnum/2) %trial number
    
    for y = 1:(args.DFencode.trialnum/2) %trial number
        
        results.groupmean.preview_DFencode_F(x,y) = (exp(2*results.groupmean.preview_DFencode_F_mz(x, y))-1)/(exp(2*results.groupmean.preview_DFencode_F_mz(x, y))+1);
    end
end

for x = 1:(args.DFencode.trialnum/2) %trial number
    
    for y = 1:(args.DFencode.trialnum/2) %trial number
        
        results.groupmean.preview_DFencode_R(x,y) = (exp(2*preview_DFencode_stack_R_mz(x, y))-1)/(exp(2*preview_DFencode_stack_R_mz(x, y))+1);
    end
end




cd(args.output_dir)
results.parameters = args;
save(args.outfname,'results')

    %DFencode 
    
    preview_DFencode_fig = figure;
    set(preview_DFencode_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.preview_DFencode); colormap('parula'); colorbar; 
    
    
    ylabel(sprintf('preview_DFencode raw patterns averaged over TR %d: TR %d', TRstart, TRlength), 'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('preview_DFencode raw patterns averaged over TR %d: TR %d', TRstart, TRlength),'FontSize',15,'FontWeight','bold');
    saveas(preview_DFencode_fig, sprintf('emodif_average_run preview_DFencode_TRstart%d_mean%d', TRstart, TRlength),'png');
    
    preview_DFencodez_fig = figure;
    set(preview_DFencodez_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.preview_DFencode_mz); colormap('parula'); colorbar;
    
    
    ylabel(sprintf('preview_DFencode_mz mean patterns averaged over TR %d: TR %d', TRstart, TRlength), 'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('preview_DFencode_mz mean patterns averaged over TR %d: TR %d', TRstart, TRlength),'FontSize',15,'FontWeight','bold');
    saveas(preview_DFencodez_fig, sprintf('emodif_average_run_preview_DFencode_mz_TRstart%d_mean%d', TRstart, TRlength),'png');
    
    % remember DFencode
    
    preview_DFencode_R_fig = figure;
    set(preview_DFencode_R_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.preview_DFencode_R); colormap('parula'); colorbar; 
    
    
    ylabel(sprintf('preview_DFencode_R raw patterns averaged over TR %d: TR %d', TRstart, TRlength), 'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('preview_DFencode_R raw patterns averaged over TR %d: TR %d', TRstart, TRlength),'FontSize',15,'FontWeight','bold');
    saveas(preview_DFencode_R_fig, sprintf('emodif_average_run preview_DFencode_R_TRstart%d_mean%d', TRstart, TRlength),'png');
    
    preview_DFencode_Rz_fig = figure;
    set(preview_DFencode_Rz_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.preview_DFencode_R_mz); colormap('parula'); colorbar; 
    
    
    ylabel(sprintf('preview_DFencodez_R mean patterns averaged over TR %d: TR %d', TRstart, TRlength), 'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('preview_DFencodez_R mean patterns averaged over TR %d: TR %d', TRstart, TRlength),'FontSize',15,'FontWeight','bold');
    saveas(preview_DFencode_Rz_fig, sprintf('emodif_average_run preview_DFencode_Rmz_TRstart%d_mean%d', TRstart, TRlength),'png');
    
    %Forget
    
    preview_DFencode_F_fig = figure;
    set(preview_DFencode_F_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.preview_DFencode_F); colormap('parula'); colorbar; 
    
    
    ylabel(sprintf('preview_DFencode_F raw patterns averaged over TR %d: TR %d', TRstart, TRlength), 'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('preview_DFencode_F raw patterns averaged over TR %d: TR %d', TRstart, TRlength),'FontSize',15,'FontWeight','bold');
    saveas(preview_DFencode_F_fig, sprintf('emodif_average_run preview_DFencode_F_TRstart%d_mean%d', TRstart, TRlength),'png');
    
    preview_DFencodez_Fz_fig = figure;
    set(preview_DFencodez_Fz_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.preview_DFencode_F_mz); colormap('parula'); colorbar;
    
    
    ylabel(sprintf('preview_DFencode_Fz mean patterns averaged over TR %d: TR %d', TRstart, TRlength), 'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('preview_DFencode_Fz mean patterns averaged over TR %d: TR %d', TRstart, TRlength),'FontSize',15,'FontWeight','bold');
    saveas(preview_DFencodez_Fz_fig, sprintf('emodif_average_run_preview_DFencode_Fmz_TRstart%d_mean%d', TRstart, TRlength),'png');
    

cd(args.script_dir)

elseif type == 5
    
args.output_dir = sprintf('%s/aggregate_results/RSA_localizer_preview_dfencode_results/%s', args.base_dir, maskName);
args.outfname = sprintf('%s/rsa_preview_DF_local_aggregate_%s', args.output_dir, date);
mkdir(args.output_dir);
    

    local_preview_stack_mz = [];
    local_DFstack_mz  = [];
    local_DFstack_F_mz = [];
    local_DFstack_R_mz  = [];
    
    for x = 1:length(subj_list)
        local_preview_stack_mz = cat(3, local_preview_stack_mz, results.bysubject.data(x).rsa.results.smatrix.corr_matrix_LPz);
        local_DFstack_mz  = cat(3, local_DFstack_mz , results.bysubject.data(x).rsa.results.smatrix.corr_matrix_LDFz);
        local_DFstack_F_mz = cat(3, local_DFstack_F_mz, results.bysubject.data(x).rsa.results.smatrix.corr_matrix_LDF_Fz);
        local_DFstack_R_mz = cat(3, local_DFstack_R_mz , results.bysubject.data(x).rsa.results.smatrix.corr_matrix_LDF_Rz);

    end
    
   

%take the mean of the z-scored stack
local_preview_stack_mz = nanmean(local_preview_stack_mz,3);
local_DFstack_mz   = nanmean(local_DFstack_mz ,3);
local_DFstack_F_mz  = nanmean(local_DFstack_F_mz ,3);
 local_DFstack_R_mz = nanmean( local_DFstack_R_mz,3);
 
%set results file
results.groupmean.local_preview_stack_mz= local_preview_stack_mz;
results.groupmean.local_DFstack_mz  = local_DFstack_mz;
results.groupmean.local_DFstack_F_mz= local_DFstack_F_mz;
results.groupmean.local_DFstack_R_mz= local_DFstack_R_mz;
%fisher z transform everything back to pearson's R

% FISHERINV Equation 
%r = (exp(2*z)-1)/(exp(2*z)+1);
%FISHER Equation
%z = 0.5*log((1+r)/(1-r));


for x = 1:6 %trial number
    
    for y = 1:30 %trial number
        
        results.groupmean.local_preview_stack_m(x,y) = (exp(2*results.groupmean.local_preview_stack_mz(x, y))-1)/(exp(2*results.groupmean.local_preview_stack_mz(x, y))+1);
    end
end

for x = 1:6 %trial number
    
    for y = 1:30 %trial number
        
        results.groupmean.local_DFstack_m(x,y) = (exp(2*results.groupmean.local_DFstack_mz(x, y))-1)/(exp(2*results.groupmean.local_DFstack_mz(x, y))+1);
    end
end

for x = 1:6 %trial number
    
    for y = 1:30 %trial number
        
        results.groupmean.local_DFstack_F_m(x,y) = (exp(2*results.groupmean.local_DFstack_F_mz(x, y))-1)/(exp(2*results.groupmean.local_DFstack_F_mz(x, y))+1);
    end
end

for x = 1:6 %trial number
    
    for y = 1:30 %trial number
        
        results.groupmean.local_DFstack_R_m(x,y) = (exp(2*results.groupmean.local_DFstack_R_mz(x, y))-1)/(exp(2*results.groupmean.local_DFstack_R_mz(x, y))+1);
    end
end

results.groupmean.local_preview_stack_mz= local_preview_stack_mz;
results.groupmean.local_DFstack_mz  = local_DFstack_mz;
results.groupmean.local_DFstack_F_mz= local_DFstack_F_mz;
results.groupmean.local_DFstack_R_mz= local_DFstack_R_mz;



cd(args.output_dir)
results.parameters = args;
save(args.outfname,'results')

    %DFencode 
    
    local_preview_fig = figure;
    set(local_preview_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.local_preview_stack_m); colormap('parula'); colorbar; 
    
    
    ylabel(sprintf('localizer raw patterns averaged over TR %d: TR %d', TRstart, TRlength), 'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('preview raw patterns averaged over TR %d: TR %d', TRstart, TRlength),'FontSize',15,'FontWeight','bold');
    saveas(local_preview_fig, sprintf('emodif_Local_Preview_crossphase_avg_TRstart%d_mean%d', TRstart, TRlength),'png');
    
    local_DF_fig = figure;
    set(local_DF_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.local_DFstack_m); colormap('parula'); colorbar; 
    
    
    ylabel(sprintf('localizer raw patterns averaged over TR %d: TR %d', TRstart, TRlength), 'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('DFencode raw patterns averaged over TR %d: TR %d', TRstart, TRlength),'FontSize',15,'FontWeight','bold');
    saveas(local_DF_fig, sprintf('emodif_Local_DFencode_crossphase_avg_TRstart%d_mean%d', TRstart, TRlength),'png');
    
    
     local_DF_F_fig = figure;
    set(local_DF_F_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.local_DFstack_F_m); colormap('parula'); colorbar; 
    
    
    ylabel(sprintf('localizer raw patterns averaged over TR %d: TR %d', TRstart, TRlength), 'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('DFencodeF raw patterns averaged over TR %d: TR %d', TRstart, TRlength),'FontSize',15,'FontWeight','bold');
    saveas(local_DF_F_fig, sprintf('emodif_Local_DFencodeF_crossphase_avg_TRstart%d_mean%d', TRstart, TRlength),'png');
    
     local_DF_R_fig = figure;
    set(local_DF_R_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.local_DFstack_R_m); colormap('parula'); colorbar; 
    
    
    ylabel(sprintf('localizer raw patterns averaged over TR %d: TR %d', TRstart, TRlength), 'FontSize',15,'FontWeight','bold');
    xlabel(sprintf('DFencodeR raw patterns averaged over TR %d: TR %d', TRstart, TRlength),'FontSize',15,'FontWeight','bold');
    saveas(local_DF_R_fig, sprintf('emodif_Local_DFencodeR_crossphase_avg_TRstart%d_mean%d', TRstart, TRlength),'png');
    

cd(args.script_dir)



% %fisher's z transform
% 
% rsa_local_101replace = Replace(rsa_local_101.results.smatrix.corr_matrix_match_Localizer_bycatz,inf,NaN);
% 
% rsa_local_101.catword = cat(3,rsa_local_101replace(:,25:30),...
%     rsa_local_101replace(:,19:24));
% 
% rsa_local_101.catwordm = nanmean(rsa_local_101.catword,3);
% 
% rsa_local_101.replace_cat = vertcat(
% 
% stackedmean = cat(3,rsa_local_101.results.smatrix.corr_matrix_match_Localizer_bycatz(1:24,1:24),...
%     rsa_local_102.results.smatrix.corr_matrix_match_Localizer_bycatz(1:24,1:24),...
%     rsa_local_103.results.smatrix.corr_matrix_match_Localizer_bycatz(1:24,1:24),...
%     rsa_local_104.results.smatrix.corr_matrix_match_Localizer_bycatz(1:24,1:24));
% 
% stackedmeanreplace = Replace(stackedmean,inf,1);
% 
% rsa_mean = nanmean(stackedmeanreplace,3);
% 
%     Localizer_fig = figure;
%     set(Localizer_fig, 'Position', [0 0 1500 1500])
%     
%     subplot(1,1,1)
%     imagesc(rsa_mean); colormap('parula'); colorbar; 
%     
%     
%     ylabel('Localizer mean patterns averaged over TR 4: TR 3 ','FontSize',15,'FontWeight','bold');
%     xlabel('Localizer mean patterns averaged over TR 3: TR 4 ','FontSize',15,'FontWeight','bold');
%     saveas(Localizer_fig,'emodif_average_run_Localizer_TRstart4_mean3','png');
end