function[] = imdif_analyze_outputs_logbayes_chunks(analysis_id, varargin)
  
  % REPREF_ANALYZE_OUTPUTS(...)
  %
  % Purpose:
  %
  % Analyse the expectation maximization plus importance sampling results
  %
  % Inputs:
  %
  % analysis_id: Valid analysis Id
  % varargin
  %       If 'regular' then it also expects the total number of runs; optional - original run analysis id
  %       If 'bootstrap' then it also expects the total number of bootstrap runs; optional - original run analysis id
  % 	  If 'scramble' then it also expects the total number of scramble runs; optional - original run analysis id
  %f
  % Outputs:
  %
  % Bunch of plots
  %
  % Example usage:
  %
  % analyze_outputs('my_analysis_id','regular', 10)
  %
  % analyze_outputs('my_analysis_id','regular_rescale_x', 10, x_min, x_max)
  % 
  % analyze_outputs('my_bootstrap_id', 'bootstrap', 100, 'my_original_id')
  % analyze_outputs('my_bootstrap_id', 'bootstrap', 100)
  % 
  % analyze_outputs('my_bootstrap_id', 'bootstrap_outliers_semipartial', 100, 'my_other_id') 
  % 
  % analyze_outputs('my_bootstrap_id', 'bootstrap_outliers_ptc', 100, 0.05)
  %     % 0.05 = cutoff fraction of lower PTC values
  %
  % analyze_outputs('my_bootstrap_id', 'bootstrap_outerloop', 100,'my_bootstrap_id2')
  %
  % analyze_outputs('my_scramble_id', 'scramble', 100, 0.750, 0.492) 
  %     % 0.750 = actual P(theory consistency) for non-scrambled analysis
  %     % 0.492 = actual beta1 for non-scrambled analysis
  % analyze_outputs('my_scramble_id', 'scramble', 100) % 
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % This is part of the P-CIT toolbox released under the BSD license.
  % Copyright (c) 2012, Princeton University
  % All rights reserved.
  %
  % Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
  %
  % Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
  % Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in %the documentation and/or other materials provided with the distribution.
  % Neither the name of the Princeton University nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
  % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT %NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  close all;
  
  % Checks if the correct number of arguments are passed in
  if nargin < 1, error('Missing input parameters'); end
  
  % Set paths
  % Setting the target directory
  root_dir = '~/Dropbox (LewPeaLab)/STUDY/IMDiF/Results/P-CIT/';
  results_dir = fullfile(root_dir, 'bootstrap_results');
  
  read_dir = fullfile(results_dir, analysis_id);
  write_dir = fullfile(results_dir, analysis_id);
  
  x_range = [0 1];
  
  if length(varargin) >= 3
    switch varargin{1}
      case 'regular_rescale_x'
        x_range(1) = varargin{3};
        x_range(2) = varargin{4};
      case {'bootstrap_outerloop', 'bootstrap_outliers_semipartial'}
        analysis_id2 = varargin{3};
        read_dir2 = fullfile(results_dir, analysis_id2);
        write_dir2 = fullfile(results_dir, analysis_id2);
      case 'bootstrap_outliers_ptc'
        cutoff_fraction = varargin{3};
      case 'scramble'
        actual_ptc = varargin{3};
        actual_beta1 = varargin{4};
        
        if ischar(actual_ptc); actual_ptc = str2double(actual_ptc); end    
        if ischar(actual_beta1); actual_beta1 = str2double(actual_beta1); end    
    end
  end  
  
  resolution = 2;
  credible_interval = 0.9; % credible interval a 0.9 implies that 90% of samples should lie in that interval
  n_bins = 10; % number of bins to plot histgrams
  image_format = 'eps';
  alpha = 0.05;
  visible_off = true; % control whether or not the generated figures pop up on screen or not before being saved to disk.
  
  % Check length of varargin to determine the type of analysis we're plotting, either one of the following: {regular, bootstrap, scrambled}
  if length(varargin) > 1 && length(varargin) <= 4
    
    nRuns = varargin{2};
    if ischar(nRuns); nRuns = str2double(nRuns); end    
    
    chunk_start = varargin{3};
    chunk_end = varargin{4};

    x = 0:(1/10^resolution):1;
    y_all_runs = NaN(nRuns, length(x));
    weight_all_runs = NaN(nRuns, 1);
    runs_valid = NaN(nRuns,1);
    
    min_pval = 0.000000000001;
    max_chi2 = chi2inv(1-min_pval,1);
    
    switch varargin{1}
      
      case {'regular', 'regular_rescale_x'}
        % Individually generate figures for each run, accumulating particle weights and y values for each run
                
        p_wgts = NaN(nRuns,1);
        beta1 = NaN(nRuns,1);
        fminunc = NaN(nRuns,1);
        lrt_beta1 = NaN(nRuns,1);
        chi2_beta1 = NaN(nRuns,1);
        
        for b = 1:nRuns
          fprintf('Regular run %d\n', b);
          [pcit_stats] = plot_figures(sprintf('%s_r%d', analysis_id, b), read_dir, write_dir, resolution,...
            x_range, credible_interval, n_bins, image_format, visible_off, true);
          
          runs_valid(b) = pcit_stats.valid;
          if runs_valid(b)
            y_all_runs(b,:) = pcit_stats.weighted_curve_struct.y_final;
            
            p_wgts(b) = pcit_stats.p_wgts;
            beta1(b) = pcit_stats.beta_1(end);
            fminunc(b) = pcit_stats.fminunc_fvals(end);
            lrt_beta1(b) = pcit_stats.likelihood_ratio_test_for_beta1;
            chi2_beta1(b) = pcit_stats.chi2_beta1;
            if lrt_beta1(b) < min_pval; chi2_beta1(b) = max_chi2; end
          end
          
        end
        
        bad_runs = (runs_valid ~= 1);
        y_all_runs(bad_runs,:) = [];
        
        p_wgts(bad_runs) = [];
        beta1(bad_runs) = [];
        fminunc(bad_runs) = [];
        lrt_beta1(bad_runs) = [];
        chi2_beta1(bad_runs) = [];
        
        % Print out the average stats!
        total_runs = sum(runs_valid);
        log_bayes_ptc = calc_log_bayes_PTC(mean(p_wgts)); 
        
        fprintf('\n\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n');
        fprintf('Regular runs stats (%d valid runs of %d)\n',total_runs,nRuns);
        fprintf('Analysis (%s)\n\n',analysis_id);
        fprintf('lgBayes\tPTC\tPTC>.583\tBeta1\tFminunc\tchi2_B1\tLRT_B1\n');
        fprintf('%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%.3f\t%.6f\n', log_bayes_ptc, mean(p_wgts), mean(p_wgts>0.5833), mean(beta1), mean(fminunc), mean(chi2_beta1), mean(lrt_beta1));
        fprintf('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n\n');

      case 'bootstrap'
        % Individually generate figures for each run, accumulating particle weights and y values for each run

        plot_boots = true;
        plot_boot_aggregate = true;
        plot_boxplots = true;
        plot_violinplots_distributionPlot = true;
        plot_violinplots_violin = true;
        
        p_wgts = NaN(nRuns,1);
        beta1 = NaN(nRuns,1);
        fminunc = NaN(nRuns,1);
        lrt_beta1 = NaN(nRuns,1);
        chi2_beta1 = NaN(nRuns,1);
        
        %         specific_runs = [108 121 130 137 175 190];
        
        for b = chunk_start:chunk_end
          %           myrun = specific_runs(b);
          myrun = b;
          fprintf('Bootstrap run %d of %d\n', myrun, nRuns);
          [pcit_stats] = plot_figures(sprintf('%s_b%d', analysis_id, myrun), read_dir, write_dir, resolution,...
            [], credible_interval, n_bins, image_format, visible_off, plot_boots);
          
          runs_valid(b) = pcit_stats.valid;       
          if runs_valid(b)
            y_all_runs(b, :) = pcit_stats.weighted_curve_struct.y_final;
                        
            p_wgts(b) = pcit_stats.p_wgts;
            beta1(b) = pcit_stats.beta_1(end);
            fminunc(b) = pcit_stats.fminunc_fvals(end);
            lrt_beta1(b) = pcit_stats.likelihood_ratio_test_for_beta1;
            chi2_beta1(b) = pcit_stats.chi2_beta1;
            if lrt_beta1(b) < min_pval; chi2_beta1(b) = max_chi2; end
          end
             
        end
        
        % Remove any bogus runs
        bad_runs = ~runs_valid;
        y_all_runs(bad_runs,:) = [];
        
        p_wgts(bad_runs) = [];
        beta1(bad_runs) = [];
        fminunc(bad_runs) = [];
        lrt_beta1(bad_runs) = [];
        chi2_beta1(bad_runs) = [];

        % Print out the bootstrap stats!
        total_boots = sum(runs_valid);
        log_bayes_ptc = calc_log_bayes_PTC(mean(p_wgts)); 
        log_bayes_ptc_dist = calc_log_bayes_PTC_distribution(p_wgts);
        
        % save distribution to disk
        lbptc_filename = sprintf('%s/%s_log_bayes_ptc_dist_%d_%d.mat', write_dir, analysis_id, chunk_start, chunk_end);
        save(lbptc_filename, 'total_boots', 'p_wgts', 'log_bayes_ptc', 'log_bayes_ptc_dist');
        
        fprintf('\n\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n');
        fprintf('Bootstrap stats (%d valid runs of %d)\n',total_boots,nRuns);
        fprintf('Analysis (%s)\n\n',analysis_id);
        fprintf('lgBayes\tPTC\tPTC>.5833\tBeta1\tFminunc\tchi-sq\tpb1\tpb1<.05\tlgBayes>0\n');
        fprintf('%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%.3f\t%.6f\t%.3f\t%.3f\n', ...
          log_bayes_ptc, mean(p_wgts), mean(p_wgts>0.5833), ...
          mean(beta1), mean(fminunc), mean(chi2_beta1),...
          mean(lrt_beta1), mean(lrt_beta1<.05), mean(log_bayes_ptc_dist>0));
        fprintf('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n\n');   
        
        % Print histogram  and box plot of log Bayes PTC values

        
        if plot_boxplots
          figure;
          
          h1 = subplot(1,2,1);
          hist(h1,log_bayes_ptc_dist);
          title(h1,'log_bayes_PTC histogram', 'interpreter','none');
          xlabel(h1,'log( PTC/(1-PTC)/(.583/.417) )'); ylabel(h1,'count');
          axis(h1,'square');
          
          h2 = subplot(1,2,2);
          boxplot(h2,log_bayes_ptc_dist,'symbol','r*');
          ylim(h2,[-3 9]);
          title(h2,'log_bayes_PTC boxplot', 'interpreter','none');
          xlabel(h2,'log( PTC/(1-PTC)/(.583/.417) )'); ylabel(h2,'count');
          axis(h2,'square');
          
          file_name = sprintf('%s/%s_logbayesptc_hist+box_%s', write_dir, analysis_id, datestamp);
          fprintf('Printing %s\n', file_name);
          print('-depsc','-painters','-r200',file_name);
        end
        
        if plot_violinplots_distributionPlot
          figure;
          
          h1 = subplot(1,3,1);
          hist(h1,log_bayes_ptc_dist);
          title(h1,'log_bayes_PTC histogram', 'interpreter','none');
          xlabel(h1,'log( PTC/(1-PTC)/(.583/.417) )'); ylabel(h1,'count');
          axis(h1,'square');
          
          h2 = subplot(1,3,2);          
          distributionPlot(h2, log_bayes_ptc_dist, 'showMM', 4);
          ylim(h2,[-2 10]);
          title(h2,'log_bayes_PTC violinplot', 'interpreter','none');
          xlabel(h2,'log( PTC/(1-PTC)/(.583/.417) )'); ylabel(h2,'count');
          axis(h2,'square');
          
          h3 = subplot(1,3,3);
          distributionPlot(h3, log_bayes_ptc_dist, 'showMM', 4, ...
            'distWidth', 0.95, 'addSpread', 1,  'histOpt', 1);
          ylim(h3,[-2 10]);
          title(h3,'log_bayes_PTC violinplot', 'interpreter','none');
          xlabel(h3,'log( PTC/(1-PTC)/(.583/.417) )'); ylabel(h2,'count');
          axis(h3,'square');
                              
          file_name = sprintf('%s/%s_logbayesptc_hist+violin_%s', write_dir, analysis_id, datestamp);
          fprintf('Printing %s\n', file_name);
          print('-depsc','-painters','-r200',file_name);
          
        end
        
        if plot_violinplots_violin
          h1 = figure;
          hist(log_bayes_ptc_dist);
          title('log_bayes_PTC histogram', 'interpreter','none');
          xlabel('log( PTC/(1-PTC)/(.583/.417) )'); ylabel('count');
          axis('square');
          
          [h2,L,MX,MED] = violin(log_bayes_ptc_dist,[]);
          ylim([-4 12]);
          title(sprintf('log_bayes_PTC violin (MX=%.3f, MED=%.3f', MX, MED), 'interpreter','none');
          xlabel('log( PTC/(1-PTC)/(.583/.417) )'); ylabel('count');
          axis('square');
          
          hist(log_bayes_ptc_dist);
          title('log_bayes_PTC histogram', 'interpreter','none');
          xlabel('log( PTC/(1-PTC)/(.583/.417) )'); ylabel('count');
          axis('square');
          
          file_name = sprintf('%s/%s_logbayesptc_hist_%s', write_dir, analysis_id, datestamp);
          fprintf('Printing %s\n', file_name);
          print(h1,'-depsc','-painters','-r200',file_name);
          
          file_name = sprintf('%s/%s_logbayesptc_violin_%s', write_dir, analysis_id, datestamp);
          fprintf('Printing %s\n', file_name);
          print(h2,'-depsc','-painters','-r200',file_name);
          
        end
        
        % Generate a bootstrap specific plot
        if plot_boot_aggregate
          plot_boot_results(analysis_id, write_dir, x, y_all_runs, p_wgts, alpha, resolution, credible_interval,...
            image_format, n_bins, varargin{:});
        end
        
      case 'bootstrap_outerloop'
        % Compare 2 different bootstrapped analysis IDs
        %
        % For each analysis:
        % - individually generate figures for each run
        % - accumulate particle weights and y values
        % 
        % Comparing each analysis:
        % - compare chi-sq likliehood test
        % - compare fvals
        % 
        
        plot_stuff = false;
        
        a1.p_wgts = NaN(nRuns,1);
        a1.beta1 = NaN(nRuns,1);
        a1.fminunc = NaN(nRuns,1);
        a1.lrt_beta1 = NaN(nRuns,1);
        
        a2 = a1;
        a1_vs_a2 = a1;
        
        boot_true_indices = 1:nRuns;
        
        for b = 1:nRuns
          fprintf('\nBootstrap run %d of %d\n', b, nRuns);
                  
          fprintf('Analysis 1 (%s)\n',analysis_id);
          [pcit_stats1] = plot_figures(sprintf('%s_r%d', analysis_id, b), read_dir, write_dir, resolution,...
            [], credible_interval, n_bins, image_format, visible_off, plot_stuff);
          
          fprintf('Analysis 2 (%s)\n',analysis_id2);
          [pcit_stats2] = plot_figures(sprintf('%s_r%d', analysis_id2, b), read_dir2, write_dir2, resolution,...
            [], credible_interval, n_bins, image_format, visible_off, plot_stuff);
         
          runs_valid(b) = pcit_stats1.valid && pcit_stats2.valid;
          
          if runs_valid(b)
            
            % A1 stats
            a1.p_wgts(b) = pcit_stats1.p_wgts;
            a1.beta1(b) = pcit_stats1.beta_1(end);
            a1.fminunc(b) = pcit_stats1.fminunc_fvals(end);
            a1.lrt_beta1(b) = pcit_stats1.likelihood_ratio_test_for_beta1;
            a1.chi2_beta1(b) = chi2inv(1-a1.lrt_beta1(b),1);
            if a1.lrt_beta1(b) < min_pval; a1.chi2_beta1(b) = max_chi2; end

            % A2 stats
            a2.p_wgts(b) = pcit_stats2.p_wgts;
            a2.beta1(b) = pcit_stats2.beta_1(end);
            a2.fminunc(b) = pcit_stats2.fminunc_fvals(end);
            a2.lrt_beta1(b) = pcit_stats2.likelihood_ratio_test_for_beta1;
            a2.chi2_beta1(b) = chi2inv(1-a2.lrt_beta1(b),1);
            if a2.lrt_beta1(b) < min_pval; a2.chi2_beta1(b) = max_chi2; end
            
            % Calculate if analysis 1 is BETTER than analysis 2
            a1_vs_a2.p_wgts(b) = a1.p_wgts(b) > a2.p_wgts(b);
            a1_vs_a2.beta1(b) = a1.beta1(b) > a2.beta1(b);
            a1_vs_a2.fminunc(b) = a1.fminunc(b) < a2.fminunc(b);
            a1_vs_a2.lrt_beta1(b) = a1.lrt_beta1(b) < a2.lrt_beta1(b);
            a1_vs_a2.chi2_beta1(b) = a1.chi2_beta1(b) > a2.chi2_beta1(b);
            
          end
        end
        
        % Remove any bogus runs
        bad_runs = ~runs_valid;
        boot_true_indices(bad_runs) = [];
        the_fields = {'p_wgts','beta1','fminunc','lrt_beta1','chi2_beta1'};
        for f = 1:length(the_fields)
          a1.(the_fields{f})(bad_runs) = [];
          a2.(the_fields{f})(bad_runs) = [];
          a1_vs_a2.(the_fields{f})(bad_runs) = [];
        end
        
        % get distribution of f-val and chi-sq val comparisons
        a1_vs_a2.dist.fminunc = a1.fminunc - a2.fminunc;
        a1_vs_a2.dist.chi2_beta1 = a1.chi2_beta1 - a2.chi2_beta1;
        hist_edges = [-inf -16:4:16 inf];
        
        [n,bins] = histc(a1_vs_a2.dist.fminunc,hist_edges);
        a1_vs_a2.hist.fminunc.n = n;
        a1_vs_a2.hist.fminunc.x = hist_edges;
        a1_vs_a2.hist.fminunc.bins = bins;
        
        [n,bins] = histc(a1_vs_a2.dist.chi2_beta1,hist_edges);
        a1_vs_a2.hist.chi2_beta1.n = n;
        a1_vs_a2.hist.chi2_beta1.x = hist_edges;
        a1_vs_a2.hist.chi2_beta1.bins = bins;
        
        % Print out the bootstrap comparision stats!
        total_boots = sum(runs_valid);
        a1.log_bayes_ptc = calc_log_bayes_PTC(mean(a1.p_wgts));
        a2.log_bayes_ptc = calc_log_bayes_PTC(mean(a2.p_wgts));
        
        fprintf('\n\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n');
        fprintf('Bootstrap outerloop comparision (%d valid boots of %d)\n',total_boots,nRuns);
        fprintf('Analysis 1 (%s) vs. Analysis 2 (%s)\n\n',analysis_id,analysis_id2);
        fprintf('\tlgBayes\tPTC\tPTC>.583\tBeta1\tFminunc\tchi-sq\tpb1\tpb1<.05\n');
        fprintf('A1\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%.3f\t%.6f\t%.3f\n', a1.log_bayes_ptc, mean(a1.p_wgts), mean(a1.p_wgts>0.5833), std(a1.p_wgts), mean(a1.beta1), mean(a1.fminunc), mean(a1.chi2_beta1), mean(a1.lrt_beta1), mean(a1.lrt_beta1<.05));
        fprintf('A2\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%.3f\t%.6f\t%.3f\n', a2.log_bayes_ptc, mean(a2.p_wgts), mean(a2.p_wgts>0.5833), std(a2.p_wgts), mean(a2.beta1), mean(a2.fminunc), mean(a2.chi2_beta1), mean(a2.lrt_beta1), mean(a2.lrt_beta1<.05));
        fprintf('A1>A2\t%.3f\t--\t%.3f\t%.3f\t--\t%.3f\t--\t--\n', sum(a1_vs_a2.p_wgts)/total_boots, sum(a1_vs_a2.beta1)/total_boots, sum(a1_vs_a2.chi2_beta1)/total_boots);
        fprintf('A1<A2\t--\t--\t--\t--\t%.2f\t--\t%.6f\t--\n', sum(a1_vs_a2.fminunc)/total_boots, sum(a1_vs_a2.lrt_beta1)/total_boots);
        
        fprintf('\n\nA1-A2:Fminunc\nbin\t');
        fprintf('%d\t',a1_vs_a2.hist.fminunc.x);
        fprintf('\ncumul%%\t');
        fprintf('%.3f\t',cumsum(a1_vs_a2.hist.fminunc.n)/total_boots);
        fprintf('\nA1-A2:chi2\nbin\t');
        fprintf('%d\t',a1_vs_a2.hist.chi2_beta1.x);
        fprintf('\ncumul%%\t');
        fprintf('%.3f\t',cumsum(a1_vs_a2.hist.chi2_beta1.n)/total_boots);
        fprintf('\n');
        fprintf('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n\n');  
        
        % Look at "worst offenders" and print out curves for each analysis
        worst_offenders_bincutoff = 3;
        worst_offenders_idx = a1_vs_a2.hist.chi2_beta1.bins < worst_offenders_bincutoff;
        worst_offenders_runs = find(worst_offenders_idx);
        worst_offenders_true_boot_indices = boot_true_indices(worst_offenders_runs);
                
        fprintf('\nHave %d "worst offenders" in %d boots (A1-A2 chisq < %d)\n', ...
          sum(worst_offenders_idx),total_boots,hist_edges(worst_offenders_bincutoff));
        fprintf('Worse Offender RUNS: %s\n',mat2str(worst_offenders_true_boot_indices));        
        fprintf('Worse Offender CHI2diff: %s\n\n',mat2str(a1_vs_a2.dist.chi2_beta1(worst_offenders_idx))); 
        
        for wo = 1:length(worst_offenders_true_boot_indices)
        
          plot_stuff = true;
          run_idx = worst_offenders_true_boot_indices(wo);
          fprintf('\n\nPrinting "worst offenders": run %d (# %d of %d offenders)\n', run_idx, wo, length(worst_offenders_runs));
                    
          fprintf('Analysis 1 (%s)\n',analysis_id);
          [pcit_stats1] = plot_figures(sprintf('%s_r%d', analysis_id, run_idx), read_dir, write_dir, resolution,...
            [], credible_interval, n_bins, image_format, visible_off, plot_stuff);
          
          fprintf('Analysis 2 (%s)\n',analysis_id2);
          [pcit_stats2] = plot_figures(sprintf('%s_r%d', analysis_id2, run_idx), read_dir2, write_dir2, resolution,...
            [], credible_interval, n_bins, image_format, visible_off, plot_stuff);
          
        end
         
                
        % %         % Generate a bootstrap specific plot
        % %         plot_boot_results(analysis_id, write_dir, x, y_all_runs, weight_all_runs, alpha, resolution, credible_interval,...
        % %           image_format, n_bins, varargin{:});
           
      case 'bootstrap_outliers_semipartial'
        % This case is specifically for the semipartial bootstraps in which
        % outliers are defined as semipartial runs with low beta1s/ high p(beta1s).
        %
        % For each analysis:
        % - accumulate particle weights and y values
        %
        % Comparing each analysis:
        % -
        %
        
        plot_boots = false;
        plot_boot_aggregate = true;
        
        p_wgts = NaN(nRuns,1);
        beta1 = NaN(nRuns,1);
        fminunc = NaN(nRuns,1);
        lrt_beta1 = NaN(nRuns,1);
        chi2_beta1 = NaN(nRuns,1);
        
        boot_true_indices = 1:nRuns;
        
        for b = 1:nRuns
          %           myrun = specific_runs(b);
          myrun = b;
          fprintf('Bootstrap run %d of %d\n', myrun, nRuns);
          [pcit_stats] = plot_figures(sprintf('%s_r%d', analysis_id, myrun), read_dir, write_dir, resolution,...
            [], credible_interval, n_bins, image_format, visible_off, plot_boots);
          
          runs_valid(b) = pcit_stats.valid;       
          if runs_valid(b)
            y_all_runs(b, :) = pcit_stats.weighted_curve_struct.y_final;
                        
            p_wgts(b) = pcit_stats.p_wgts;
            beta1(b) = pcit_stats.beta_1(end);
            fminunc(b) = pcit_stats.fminunc_fvals(end);
            lrt_beta1(b) = pcit_stats.likelihood_ratio_test_for_beta1;
            chi2_beta1(b) = pcit_stats.chi2_beta1;
            if lrt_beta1(b) < min_pval; chi2_beta1(b) = max_chi2; end
          end
             
        end
        
        % Remove any bogus runs
        bad_runs = ~runs_valid;
        boot_true_indices(bad_runs) = [];
        p_wgts(bad_runs) = [];
        beta1(bad_runs) = [];
        fminunc(bad_runs) = [];
        lrt_beta1(bad_runs) = [];
        chi2_beta1(bad_runs) = [];
        
        % Print out the bootstrap stats!
        total_boots = sum(runs_valid);
        fprintf('\n\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n');
        fprintf('Bootstrap stats (%d valid runs of %d)\n',total_boots,nRuns);
        fprintf('Analysis (%s)\n\n',analysis_id);
        fprintf('PTC\tPTC>.5\tBeta1\tFminunc\tchi-sq\tpb1\tpb1<.05\n');
        fprintf('%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%.3f\t%.6f\t%.3f\n', mean(p_wgts), mean(p_wgts>0.5), mean(beta1), mean(fminunc), mean(chi2_beta1), mean(lrt_beta1), mean(lrt_beta1<.05));
        fprintf('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n\n');   
        
        % Generate a bootstrap specific plot
        if plot_boot_aggregate
          plot_boot_results(analysis_id, write_dir, x, y_all_runs, p_wgts, alpha, resolution, credible_interval,...
            image_format, n_bins, varargin{:});
        end
        
        % Print histogram of beta1 values
        hist(beta1);
        title(sprintf('beta1 histogram: analysis_id=%s', analysis_id), 'interpreter','none');
        xlabel('beta1'); ylabel('count');
        file_name = sprintf('%s/%s_beta1_hist', write_dir, analysis_id);
        fprintf('Printing %s\n', file_name);
        print('-depsc','-painters','-r200',file_name);
        
        % Print histogram of p(beta1) values
        hist(lrt_beta1);
        title(sprintf('p(beta1) histogram: analysis_id=%s', analysis_id), 'interpreter','none');
        xlabel('p(beta1)'); ylabel('count');
        file_name = sprintf('%s/%s_p(beta1)_hist', write_dir, analysis_id);
        fprintf('Printing %s\n', file_name);
        print('-depsc','-painters','-r200',file_name);
                
        % Look at "worst offenders" and print out curves for each analysis
        worst_offenders_pcutoff = 0.5;
        worst_offenders_idx = lrt_beta1 > worst_offenders_pcutoff;
        worst_offenders_runs = find(worst_offenders_idx);
        worst_offenders_true_boot_indices = boot_true_indices(worst_offenders_runs);
        
        fprintf('\nHave %d "worst offenders" in %d boots [p(beta1) > %.2f]\n', ...
          sum(worst_offenders_idx), total_boots, worst_offenders_pcutoff);
        fprintf('Worse Offender RUNS: %s\n',mat2str(worst_offenders_true_boot_indices));
        fprintf('Worse Offender P(BETA1): %s\n\n',mat2str(lrt_beta1(worst_offenders_idx)));
        
        for wo = 1:length(worst_offenders_true_boot_indices)
        
          plot_stuff = true;
          run_idx = worst_offenders_true_boot_indices(wo);
          fprintf('\n\nPrinting "worst offenders": run %d (# %d of %d offenders)\n', run_idx, wo, length(worst_offenders_runs));
                    
          fprintf('Analysis 1 (%s)\n',analysis_id);
          [pcit_stats1] = plot_figures(sprintf('%s_r%d', analysis_id, run_idx), read_dir, write_dir, resolution,...
            [], credible_interval, n_bins, image_format, visible_off, plot_stuff);
          
          fprintf('Analysis 2 (%s)\n',analysis_id2);
          [pcit_stats2] = plot_figures(sprintf('%s_r%d', analysis_id2, run_idx), read_dir2, write_dir2, resolution,...
            [], credible_interval, n_bins, image_format, visible_off, plot_stuff);
          
        end 
        
      case 'bootstrap_outliers_ptc'
        % This case is for looking at the lowest 10% of PTC boots
        %
        
        plot_boots = false;
        plot_boot_aggregate = false;
        
        p_wgts = NaN(nRuns,1);
        beta1 = NaN(nRuns,1);
        fminunc = NaN(nRuns,1);
        lrt_beta1 = NaN(nRuns,1);
        chi2_beta1 = NaN(nRuns,1);
        
        boot_true_indices = 1:nRuns;
        
        for b = 1:nRuns
          myrun = b;
          fprintf('Bootstrap run %d of %d\n', myrun, nRuns);
          [pcit_stats] = plot_figures(sprintf('%s_r%d', analysis_id, myrun), read_dir, write_dir, resolution,...
            [], credible_interval, n_bins, image_format, visible_off, plot_boots);
          
          runs_valid(b) = pcit_stats.valid;       
          if runs_valid(b)
            y_all_runs(b, :) = pcit_stats.weighted_curve_struct.y_final;
                        
            p_wgts(b) = pcit_stats.p_wgts;
            beta1(b) = pcit_stats.beta_1(end);
            fminunc(b) = pcit_stats.fminunc_fvals(end);
            lrt_beta1(b) = pcit_stats.likelihood_ratio_test_for_beta1;
            chi2_beta1(b) = pcit_stats.chi2_beta1;
            if lrt_beta1(b) < min_pval; chi2_beta1(b) = max_chi2; end
          end
             
        end
        
        % Remove any bogus runs
        bad_runs = ~runs_valid;
        boot_true_indices(bad_runs) = [];
        p_wgts(bad_runs) = [];
        beta1(bad_runs) = [];
        fminunc(bad_runs) = [];
        lrt_beta1(bad_runs) = [];
        chi2_beta1(bad_runs) = [];
        
        % Print out the bootstrap stats!
        total_boots = sum(runs_valid);
        fprintf('\n\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n');
        fprintf('Bootstrap stats (%d valid runs of %d)\n',total_boots,nRuns);
        fprintf('Analysis (%s)\n\n',analysis_id);
        fprintf('PTC\tPTC>.5\tBeta1\tFminunc\tchi-sq\tpb1\tpb1<.05\n');
        fprintf('%.3f\t%.3f\t%.3f\t%.2f\t%.3f\t%.6f\t%.3f\n', mean(p_wgts), mean(p_wgts>0.5), mean(beta1), mean(fminunc), mean(chi2_beta1), mean(lrt_beta1), mean(lrt_beta1<.05));
        fprintf('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n\n');   
        
        % Generate a bootstrap specific plot
        if plot_boot_aggregate
          plot_boot_results(analysis_id, write_dir, x, y_all_runs, p_wgts, alpha, resolution, credible_interval,...
            image_format, n_bins, varargin{:});
        end
        
        % Print histogram of beta1 values
        hist(beta1);
        title(sprintf('beta1 histogram: analysis_id=%s', analysis_id), 'interpreter','none');
        xlabel('beta1'); ylabel('count');
        file_name = sprintf('%s/%s_beta1_hist', write_dir, analysis_id);
        fprintf('Printing %s\n', file_name);
        print('-depsc','-painters','-r200',file_name);
        
        % Print histogram of p(beta1) values
        hist(lrt_beta1);
        title(sprintf('p(beta1) histogram: analysis_id=%s', analysis_id), 'interpreter','none');
        xlabel('p(beta1)'); ylabel('count');
        file_name = sprintf('%s/%s_p(beta1)_hist', write_dir, analysis_id);
        fprintf('Printing %s\n', file_name);
        print('-depsc','-painters','-r200',file_name);
                
        % Look at "worst offenders" and print out curves for each analysis
        [sorted_ptc order_ptc] = sort(p_wgts);
        worst_offenders_fraction_cutoff = cutoff_fraction;
        worse_offenders_cutoff_maxidx = ceil(length(p_wgts) * worst_offenders_fraction_cutoff);
        worst_offenders_idx = 1:worse_offenders_cutoff_maxidx;
        
        worst_offenders_runs = order_ptc(worst_offenders_idx);
        worst_offenders_true_boot_indices = boot_true_indices(worst_offenders_runs);
        
        fprintf('\nHave %d "worst offenders" in %d boots [lower %d%% of PTC]\n', ...
          length(worst_offenders_idx), total_boots, 100*worst_offenders_fraction_cutoff);
        fprintf('Worse Offender RUNS: %s\n',mat2str(worst_offenders_true_boot_indices));
        fprintf('Worse Offender PTC: %s\n\n',mat2str(sorted_ptc(worst_offenders_idx)));
        
        for wo = 1:length(worst_offenders_true_boot_indices)
        
          plot_stuff = true;
          run_idx = worst_offenders_true_boot_indices(wo);
          fprintf('\n\nPrinting "worst offenders": run %d (# %d of %d offenders)\n', run_idx, wo, length(worst_offenders_runs));
                    
          fprintf('Analysis (%s)\n',analysis_id);
          [pcit_stats1] = plot_figures(sprintf('%s_r%d', analysis_id, run_idx), read_dir, write_dir, resolution,...
            [], credible_interval, n_bins, image_format, visible_off, plot_stuff);
                  
        end
        
      case 'scramble'
        
        plot_stuff = false;
        
        p_wgts = NaN(nRuns,1);
        beta1 = NaN(nRuns,1);
        fminunc = NaN(nRuns,1);
        lrt_beta1 = NaN(nRuns,1);
        chi2_beta1 = NaN(nRuns,1);
        
        for b = 1:nRuns
          fprintf('Scramble run %d of %d\n', b, nRuns);
          [pcit_stats] = plot_figures(sprintf('%s_r%d', analysis_id, b), read_dir, write_dir, resolution,...
            [], credible_interval, n_bins, image_format, visible_off, plot_stuff);
          
          runs_valid(b) = pcit_stats.valid;
          if runs_valid(b)            
            p_wgts(b) = pcit_stats.p_wgts;
            beta1(b) = pcit_stats.beta_1(end);
            fminunc(b) = pcit_stats.fminunc_fvals(end);
            lrt_beta1(b) = pcit_stats.likelihood_ratio_test_for_beta1;
            chi2_beta1(b) = pcit_stats.chi2_beta1;
            if lrt_beta1(b) < min_pval; chi2_beta1(b) = max_chi2; end
          end
          
        end
        
        % Remove any bogus runs
        bad_runs = ~runs_valid;
        
        p_wgts(bad_runs) = [];
        beta1(bad_runs) = [];
        fminunc(bad_runs) = [];
        lrt_beta1(bad_runs) = [];
        chi2_beta1(bad_runs) = [];
        
        % Print out the scramble stats!
        total_boots = sum(runs_valid);
        fprintf('\n\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n');
        fprintf('Scramble stats (%d valid runs of %d)\n',total_boots,nRuns);
        fprintf('Analysis (%s)\n\n',analysis_id);
        fprintf('logBayesPTC\tPTC\tPTC>%.3f\tBeta1\tB1>%.3f\tFminunc\tchi-sq\tpb1\tpb1<.05\tB1>%.3f_AND_PTC>%.3f\n', ...
          actual_ptc, actual_beta1, actual_beta1, actual_ptc);
        fprintf('%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%.3f\t%.6f\t%.3f\t%.3f\n', ...
          calc_log_bayes_PTC(mean(p_wgts)), mean(p_wgts), mean(p_wgts>actual_ptc), mean(beta1), mean(beta1>actual_beta1), ...
          mean(fminunc), mean(chi2_beta1), mean(lrt_beta1), mean(lrt_beta1<.05), mean((beta1>actual_beta1) & (p_wgts>actual_ptc)));
        fprintf('# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n\n');
                
        
        % %         % Individually generate figures for each run, accumulating particle weights and y values for each run
        % %         for s = 1:nRuns
        % %           fprintf('Scramble run %d\n', s);
        % %           weight_all_runs(s, :) = plot_figures(sprintf('%s_s%d', analysis_id, s), read_dir, write_dir, resolution,...
        % %             [], credible_interval, n_bins, image_format, visible_off, true);
        % %         end
        % %
        % Generate a scramble specific plot
        if plot_stuff
          plot_scram_results(analysis_id, write_dir, p_wgts, beta1, resolution, image_format, n_bins, varargin{:});
        end
        
      otherwise, error('Invalid analysis type!');
    end
  elseif isempty(varargin)
    plot_figures(analysis_id, read_dir, write_dir, resolution, [], credible_interval, n_bins, image_format, visible_off);
  else
    error('Invalid number of input arguments failed to trigger any of the analyses'' pipelines!');
  end
  
end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function[pcit_stats] = plot_figures(analysis_id, read_dir, write_dir, resolution, x_range, credible_interval, n_bins, image_format, visible_off, plot_stuff)
    
    % Load the target .mat file
    targ_file = sprintf('%s/%s_importance_sampler.mat', read_dir, analysis_id);
    if ~exist(targ_file,'file')
      pcit_stats.valid = false;
      return;
    else
      importance_sampler_mat = load(targ_file);
    end
    
    % Find the particles / curves indices that are theory consistent
    %%% JLP HACK %%% change this to 'horz_indpnt' for sanity check analysis
    % %     importance_sampler_mat.analysis_settings.curve_type = 'horz_indpnt';
    % %     fprintf('** CAUTION: forcing use of ''horz_indpnt'' criteria\n');
    th_con_idx = family_of_curves(importance_sampler_mat.analysis_settings.curve_type, 'count_particles', importance_sampler_mat.curve_params);
    
    % Using the above indices fetch the sum of the particle weights. This is what we call / refer throughout as p(theory consistent)
    % Here I am fetch the weights from the last em eiteration hence the 'end'
    p_wgts = sum(importance_sampler_mat.normalized_weights(end, th_con_idx));
    
    % Calculate the p-value that our beta1 is significantly different from 0
    likelihood_ratio_test_for_beta1 = importance_sampler_mat.likratiotest;
    
    % Calculate relative changes in fminunc over time
    fminunc_fvals = importance_sampler_mat.exp_max_fval;
    [max_fval, max_fval_iteration] = max(fminunc_fvals);
    [min_fval, min_fval_iteration] = min(fminunc_fvals);
    fminunc_val_change_dir = 'DEC';
    fminunc_relative_value_name = 'max';
    if max_fval_iteration > min_fval_iteration
      fminunc_val_change_dir = 'INC';
      fminunc_relative_value_name = 'min';
      fminunc_relative_change = ((max_fval - min_fval) / min_fval);
    else
      fminunc_relative_change = ((max_fval - min_fval) / max_fval);
    end
    
    % Get the weighted sum of the curves over all particles and the associated credible interval. This is the blue line and the grey envelope around the blue line
    weighted_curve_struct = common_to_all_curves(importance_sampler_mat.analysis_settings.curve_type, 'weighted_curve', importance_sampler_mat, credible_interval, resolution);
    x_final = weighted_curve_struct.x_final;
    y_final = weighted_curve_struct.y_final;
    interval = weighted_curve_struct.interval;
    
    if ~isempty(x_range)
      x_final = linspace(x_range(1),x_range(2),length(x_final));
    else
      x_range = [0 1];
    end
    
    beta_1 = importance_sampler_mat.hold_betas_per_iter(:, 2);
    fminunc_fvals = importance_sampler_mat.exp_max_fval;

    if isfield(importance_sampler_mat,'chi2val')	
        chi2val = importance_sampler_mat.chi2val;
    else	
    	chi2val = NaN;
    end	
    
    % Save results for return.
    pcit_stats.valid = true;
    pcit_stats.p_wgts = p_wgts;
    pcit_stats.likelihood_ratio_test_for_beta1 = likelihood_ratio_test_for_beta1;
    pcit_stats.weighted_curve_struct = weighted_curve_struct;
    pcit_stats.beta_1 = beta_1;
    pcit_stats.chi2_beta1 = chi2val;
    pcit_stats.fminunc_fvals = fminunc_fvals;
    
    if plot_stuff
      % Plot the curve
      figure();
      if visible_off
        set(0,'DefaultFigureVisible', 'off');
      end
      set(gcf, 'Position', [50, 900, 500, 500]);
      
      color = [107, 107, 107] ./ 255; transparency = 0.4;
      hhh = jbfill(x_final, interval(1, :), interval(2, :), color, color, 0, transparency); hold on;
      hAnnotation = get(hhh, 'Annotation');
      hLegendEntry = get(hAnnotation', 'LegendInformation');
      set(hLegendEntry, 'IconDisplayStyle', 'off');
      plot(x_final, y_final, 'b-', 'LineWidth', 2);
      
      ylabel('Change in Memory Strength', 'FontSize', 15, 'FontName', 'Helvetica');
      ylim([-1, 1]);
      xlabel('Activation', 'FontSize', 15, 'FontName', 'Helvetica');
      xlim(x_range);
      
      grid on; set(gca, 'Layer', 'top');
      title(sprintf('P(theory consistent) = %0.4f (%s)', p_wgts, importance_sampler_mat.analysis_settings.curve_type),'Interpreter','none');
      file_name = sprintf('%s/%s_weighted_curve', write_dir, analysis_id);
      print(gcf, sprintf('-d%s', image_format), '-painters', file_name);
      disp(sprintf('1. Recovered curve plot is saved as %s.%s', file_name, image_format));
      
      figure();
      if visible_off
        set(0,'DefaultFigureVisible', 'off');
      end
      set(gcf, 'Position', [50, 900, 1000, 1000]);
      
      subplot(2, 2, 1);
      plot(0:importance_sampler_mat.analysis_settings.em_iterations, beta_1, 'bo-', 'MarkerFaceColor', 'b');
      
      ylabel('beta 1');
      ylim([-0.2, 2.2]);
      xlabel('EM iterations');
      
      title({'Beta 1 over em iterations' sprintf('P(Beta 1 = 0) = %0.6f',likelihood_ratio_test_for_beta1 )});
      grid on; set(gca, 'Layer', 'top');
      
      subplot(2, 2, 2);
      hist(importance_sampler_mat.normalized_weights(end, :), n_bins); hold on;
      h1 = plot(max(importance_sampler_mat.normalized_weights(end, :)), 1, 'ro', 'MarkerFaceColor', 'r');
      
      ylabel('count', 'FontSize', 12, 'FontName', 'Helvetica');
      xlabel('posterior weights' , 'FontSize', 12, 'FontName', 'Helvetica');
      
      title({ sprintf('Distribution of posterior weights (%d iteration)', importance_sampler_mat.analysis_settings.em_iterations)  sprintf('sample max weight=%0.4f', max(importance_sampler_mat.normalized_weights(end, :))) });
      
      grid on; set(gca, 'Layer', 'top');
      legend([h1], 'Max weight');
      
      subplot(2, 2, [3 4]);
      plot(1:importance_sampler_mat.analysis_settings.em_iterations, fminunc_fvals, 'bo-', 'MarkerFaceColor', 'b');
      
      ylabel('fminunc fval', 'FontSize', 12, 'FontName', 'Helvetica');
      xlabel('EM iterations', 'FontSize', 12, 'FontName', 'Helvetica');
      
      title({'fminunc fval over em iterations' sprintf('fminunc %s over iterations by %0.4f of %s fval',fminunc_val_change_dir,fminunc_relative_change,fminunc_relative_value_name)});
      grid on; set(gca, 'Layer', 'top');
      
      file_name = sprintf('%s/%s_report_plot', write_dir, analysis_id);
      savesamesize(gcf, 'file', file_name, 'format', image_format, 'renderer', 'painters');
      disp(sprintf('2. Toolbox report plot is saved as %s.%s', file_name, image_format));
      
    end
                
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [] = plot_boot_results(analysis_id, write_dir, x, y_all_runs, weight_all_runs, alpha, resolution, credible_interval, image_format, n_bins, varargin)
    
    %     root_dir = pwd;
    %     total_bootstrap_runs = varargin{2};
    
    root_dir = '/jukebox/norman/jalewpea/fmri/repref/p-cit';
    
    % allow for some failed boots
    total_bootstrap_runs = size(y_all_runs,1);
    
    sample_idx(1) = floor(total_bootstrap_runs * (alpha / 2));
    if sample_idx(1) == 0, sample_idx(1) = ceil(total_bootstrap_runs * (alpha / 2)); end
    sample_idx(2) = ceil(total_bootstrap_runs * (alpha / 2));
    sample_idx(3) = floor(total_bootstrap_runs * (1 - (alpha / 2)));
    if sample_idx(3) == 0, ceil(total_bootstrap_runs * (1 - (alpha / 2))); end
    sample_idx(4) = ceil(total_bootstrap_runs * (1 - (alpha / 2)));
    
    sorted_y_all_runs = sort(y_all_runs);
    envelope_bounds = NaN(2, length(x));
    envelope_bounds(1, :) = (sorted_y_all_runs(sample_idx(1), :) + sorted_y_all_runs(sample_idx(2), :)) / 2;
    envelope_bounds(2, :) = (sorted_y_all_runs(sample_idx(3), :) + sorted_y_all_runs(sample_idx(4), :)) / 2;
    
    legend_str = '';
    if length(varargin) == 3 && ~strcmp(varargin{1},'bootstrap_outerloop')
      results_dir = fullfile(root_dir, 'results');
      read_dir = fullfile(results_dir, varargin{3});
      % Load the original .mat file
      importance_sampler_mat = load(sprintf('%s/%s_r1_importance_sampler.mat', read_dir, varargin{3}));
      % Find the particles / curves indices that are theory consistent
      th_con_idx = family_of_curves(importance_sampler_mat.analysis_settings.curve_type, 'count_particles', importance_sampler_mat.curve_params);
      % Using the above indices fetch the sum of the particle weights. This is what we call / refer throughout as p(theory consistent)
      % Here I am fetch the weights from the last em eiteration hence the 'end'
      weight_original = sum(importance_sampler_mat.normalized_weights(end, th_con_idx));
      
      % Get the weighted sum of the curves over all particles and the associated credible interval. This is the blue line and the grey envelope around the blue line
      weighted_curve_struct = common_to_all_curves(importance_sampler_mat.analysis_settings.curve_type, 'weighted_curve', importance_sampler_mat, credible_interval, resolution);
      y_original = weighted_curve_struct.y_final;
      legend_str = 'Original recovered curve';
    else
      weight_original = mean(weight_all_runs);
      y_original = mean(y_all_runs);
      legend_str = 'Mean(bootstrap runs)';
    end
    proportion = sum(weight_all_runs > 0.5) / length(weight_all_runs);
    proportion_str = sprintf('proportion(runs)>0.5=%0.2f', proportion);
    
    % Plot the curve
    figure(); set(gcf, 'Position', [50, 900, 1200, 500]);
    
    subplot(1, 2, 1);
    color = [107, 107, 107] ./ 255; transparency = 0.4;
    hhh = jbfill(x, envelope_bounds(1, :), envelope_bounds(2, :), color, color, 0, transparency); hold on;
    % hAnnotation = get(hhh, 'Annotation');
    % hLegendEntry = get(hAnnotation', 'LegendInformation');
    % set(hLegendEntry, 'IconDisplayStyle', 'off');
    h2 = plot(x, y_original, 'b-', 'LineWidth', 2);
    
    ylabel('Change in Memory Strength', 'FontSize', 12, 'FontName', 'Helvetica');
    ylim([-1, 1]);
    xlabel('Activation', 'FontSize', 12, 'FontName', 'Helvetica');
    xlim([0, 1]);
    
    grid on; set(gca, 'Layer', 'top');
    title(sprintf('%d bootstrap runs, recovered curves', total_bootstrap_runs));
    legend([h2(1), hhh(1)], legend_str, sprintf('%%95 confidence interval'), 'Orientation', 'Horizontal');
    
    sorted_weights = sort(weight_all_runs);
    wgt_lb = (sorted_weights(sample_idx(1)) + sorted_weights(sample_idx(2))) / 2;
    wgt_ub = (sorted_weights(sample_idx(3)) + sorted_weights(sample_idx(4))) / 2;
    
    subplot(1, 2, 2);
    numbers = hist(weight_all_runs, n_bins);
    hist(weight_all_runs, n_bins); hold on;
    h = findobj(gca, 'Type', 'patch'); set(h, 'FaceColor', [255,193,193] ./ 255, 'EdgeColor', 'w');
    hAnnotation = get(h, 'Annotation');
    hLegendEntry = get(hAnnotation', 'LegendInformation');
    set(hLegendEntry, 'IconDisplayStyle', 'off');
    
    h1 = plot(weight_original, 0:3, 'b*');
    h2 = plot([repmat(wgt_lb, 1, 4), repmat(wgt_ub, 1, 4)], [0:3, 0:3], 'b*', 'Color', color, 'Marker', '*');
    
    ylabel('Count', 'FontSize', 12, 'FontName', 'Helvetica');
    ylim([0, max(numbers)+5]);
    xlabel('Distribution of posterior weights', 'FontSize', 12, 'FontName', 'Helvetica');
    xlim([-0.02, 1.02]);
    title(sprintf('%d bootstrap runs, posterior weights\n%s', total_bootstrap_runs, proportion_str));
    legend([h1(1), h2(1)], legend_str, sprintf('%%95 confidence interval'), 'Location', 'NorthWest', 'Orientation', 'Vertical');
    
    file_name = sprintf('%s/%s_bootstrap_results_%s', write_dir, analysis_id, datestamp);
    savesamesize(gcf, 'file', file_name, 'format', image_format, 'renderer', 'painters');
    disp(sprintf('Bootstrap results plot is saved as %s.%s', file_name, image_format));
    
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function [] = plot_scram_results(analysis_id, write_dir, weight_all_runs, beta1_all_runs, resolution, image_format, n_bins, varargin)
    
    %     root_dir = pwd;
    root_dir = '/jukebox/norman/jalewpea/fmri/repref/p-cit';
    total_scramble_runs = varargin{2};
    
    legend_str = '';
    if length(varargin) == 3
      results_dir = fullfile(root_dir, 'results');
      read_dir = fullfile(results_dir, varargin{3});
      % Load the original .mat file
      importance_sampler_mat = load(sprintf('%s/%s_importance_sampler.mat', read_dir, varargin{3}));
      % Find the particles / curves indices that are theory consistent
      th_con_idx = family_of_curves(importance_sampler_mat.analysis_settings.curve_type, 'count_particles', importance_sampler_mat.curve_params);
      % Using the above indices fetch the sum of the particle weights. This is what we call / refer throughout as p(theory consistent)
      % Here I am fetch the weights from the last em eiteration hence the 'end'
      weight_original = sum(importance_sampler_mat.normalized_weights(end, th_con_idx));
      legend_str = 'Original recovered curve';
    else
      weight_original = mean(weight_all_runs);
      legend_str = 'Mean(scramble runs)';
    end
    pval = sum(weight_all_runs >= weight_original) / (total_scramble_runs + 1);
    
    % Plot the curve
    figure(); set(gcf, 'Position', [50, 900, 600, 600]);
    
    numbers = hist(weight_all_runs, n_bins);
    hist(weight_all_runs, n_bins); hold on;
    h = findobj(gca, 'Type', 'patch'); set(h, 'FaceColor', [255,193,193] ./ 255, 'EdgeColor', 'w');
    hAnnotation = get(h, 'Annotation');
    hLegendEntry = get(hAnnotation', 'LegendInformation');
    set(hLegendEntry, 'IconDisplayStyle', 'off');
    
    plot(weight_original, 0:3, 'b*');
    
    ylabel('Count', 'FontSize', 15, 'FontName', 'Helvetica');
    ylim([0, max(numbers)+5]);
    xlabel('Distribution of posterior weights', 'FontSize', 15, 'FontName', 'Helvetica');
    xlim([-0.02, 1.02]);
    title(sprintf('%d scramble samples posterior weights\npval=%0.4f', total_scramble_runs, pval));
    legend(legend_str, 'Orientation', 'Horizontal');
    
    file_name = sprintf('%s/%s_scramble_results', write_dir, analysis_id);
    savesamesize(gcf, 'file', file_name, 'format', image_format, 'renderer', 'painters');
    disp(sprintf('Scramble results plot is saved as %s.%s', file_name, image_format));
    
    % Print histogram of beta1 values
    figure;
    hist(beta1_all_runs);
    title(sprintf('beta1 histogram: analysis_id=%s', analysis_id), 'interpreter','none');
    xlabel('beta1'); ylabel('count');
    file_name = sprintf('%s/%s_beta1_hist', write_dir, analysis_id);
    fprintf('Printing %s\n', file_name);
    print('-depsc','-painters','-r200',file_name);
    
  end
  
  
  function [lb_ptc] = calc_log_bayes_PTC(meanPTC)
    
    % vert_rel ratio of theory consistent / theory inconsistent
    vert_rel_ratio = .583 / .417;
    ptc_ratio = meanPTC / (1 - meanPTC);
    
    bayes_ptc = ptc_ratio / vert_rel_ratio;
    lb_ptc =  log(bayes_ptc);
    
  end
  
  function [lb_ptc_dist] = calc_log_bayes_PTC_distribution(PTCs)
    
    % vert_rel ratio of theory consistent / theory inconsistent
    vert_rel_ratio = .583 / .417;
    ptc_ratio = PTCs ./ (1 - PTCs);
    
    bayes_ptc = ptc_ratio / vert_rel_ratio;
    lb_ptc_dist =  log(bayes_ptc);
    
  end
  
