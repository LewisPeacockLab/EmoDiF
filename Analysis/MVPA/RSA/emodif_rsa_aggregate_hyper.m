function [results] = emodif_rsa_aggregate_hyper(preview_shift)
%hyperalignment from remy uses 200 voxels
% preview TRs 3-5
% DF encode TRs 2-4
% preloaded data - the question is whether the regressors 
args.base_dir = '/Users/tw24955/emodif_data';
args.preview.TR = 366;
args.DFencode.TR = 426;

args.subjID = sprintf('emodif_%s',num2str(subjNum));
args.script_dir = pwd;
args.previewshift = preview_shift;
args.preview.trialnum = 60;
args.DFencode.trialnum = 60;

args.subj_dir = sprintf('%s/%s', args.base_dir, args.subjID);
args.data_dir = sprintf('%s/results/rsa_results/preview_DF_preview/scene_200',args.subj_dir);

args.output_dir = sprintf('%s/aggregate_results/RSA_preview_dfencode_results/scene_200', args.base_dir);
args.outfname = sprintf('%s/rsa_preview_dfencode_scene_200_aggregate_%s', args.output_dir, date);

cd(args.data_dir);
results.bysubject.data(i) = load(sprintf('%s_TR2to4_Pagg_rsa_results.mat', args.subjID));
end

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

%this takes the mean for JUST the match
PDF_Fmatch = [];
PDF_Rmatch = [];
for x = 1:30
    R_match = results.groupmean.preview_DFencode_R_mz(x,x);
    F_match = results.groupmean.preview_DFencode_F_mz(x,x);
PDF_Fmatch = vertcat(PDF_Fmatch, F_match);
PDF_Rmatch = vertcat(PDF_Rmatch, R_match);
end


cd(args.output_dir)
results.parameters = args;
save(args.outfname,'results')

    %DFencode 
    
    preview_DFencode_fig = figure;
    set(preview_DFencode_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.preview_DFencode); colormap('parula'); colorbar; 
    
    
    ylabel('preview_DFencode raw patterns', 'FontSize',15,'FontWeight','bold');
    xlabel('preview_DFencode raw patterns','FontSize',15,'FontWeight','bold');
    saveas(preview_DFencode_fig, 'emodif_average_run preview_DFencode_TRstart%d_mean%d','png');
    
    preview_DFencodez_fig = figure;
    set(preview_DFencodez_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.preview_DFencode_mz); colormap('parula'); colorbar;
    
    
    ylabel('preview_DFencode mz patterns', 'FontSize',15,'FontWeight','bold');
    xlabel('preview_DFencode mz patterns','FontSize',15,'FontWeight','bold');
    saveas(preview_DFencodez_fig, 'emodif_average_run preview_DFencode_TRstart%d_mean%d','png');
    
    % remember DFencode
    
    preview_DFencode_R_fig = figure;
    set(preview_DFencode_R_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.preview_DFencode_R); colormap('parula'); colorbar; 
    
    
    ylabel('preview_DFencode R patterns', 'FontSize',15,'FontWeight','bold');
    xlabel('preview_DFencode R patterns','FontSize',15,'FontWeight','bold');
    saveas(preview_DFencode_R_fig, 'emodif_average_run preview_DFencode_TRstart%d_mean%d','png');
    
    preview_DFencode_Rz_fig = figure;
    set(preview_DFencode_Rz_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.preview_DFencode_R_mz); colormap('parula'); colorbar; 
    
    
    ylabel('preview_DFencode Rz patterns', 'FontSize',15,'FontWeight','bold');
    xlabel('preview_DFencode Rz patterns','FontSize',15,'FontWeight','bold');
    saveas(preview_DFencode_Rz_fig, 'emodif_average_run preview_DFencode_TRstart%d_mean%d','png');
 
    
    %Forget
    
    preview_DFencode_F_fig = figure;
    set(preview_DFencode_F_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.preview_DFencode_F); colormap('parula'); colorbar; 
    
    
    ylabel('preview_DFencode F patterns', 'FontSize',15,'FontWeight','bold');
    xlabel('preview_DFencode F patterns','FontSize',15,'FontWeight','bold');
    saveas(preview_DFencode_F_fig, 'emodif_average_run preview_DFencode_TRstart%d_mean%d','png');
    
    preview_DFencodez_Fz_fig = figure;
    set(preview_DFencodez_Fz_fig, 'Position', [0 0 1500 1500])
    
    subplot(1,1,1)
    imagesc(results.groupmean.preview_DFencode_F_mz); colormap('parula'); colorbar;
    
    
    ylabel('preview_DFencode Fz patterns', 'FontSize',15,'FontWeight','bold');
    xlabel('preview_DFencode Fz patterns','FontSize',15,'FontWeight','bold');
    saveas(preview_DFencodez_Fz_fig, 'emodif_average_run preview_DFencode_TRstart%d_mean%d','png');

cd(args.script_dir)
