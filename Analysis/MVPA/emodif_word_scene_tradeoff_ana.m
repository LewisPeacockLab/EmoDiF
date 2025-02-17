function emodif_word_scene_tradeoff_ana()

%to where results directory is
load('emodif_aggregate_results_complete.mat')

%load data
R_data_word = [];
F_data_word = [];
R_data_scene = [];
F_data_scene = [];

%get all subjects all trials

for s = 1:24
    r_data_w = results.bysubject.data(s).results.remember.all.acts.word;
    f_data_w = results.bysubject.data(s).results.forget.all.acts.word;    
    r_data_s = results.bysubject.data(s).results.remember.all.acts.scene;
    f_data_s = results.bysubject.data(s).results.forget.all.acts.scene;
    
    for t = 1:length(r_data_w)
    r_data_wpost(t) = mean(r_data_w(t,5:7));
    f_data_wpost(t) = mean(f_data_w(t,5:7));
    r_data_spost(t) = mean(r_data_s(t,5:7));
    f_data_spost(t) = mean(f_data_s(t,5:7));

    end
    
    all_data_wpost = horzcat(r_data_wpost,f_data_wpost);
    all_data_spost = horzcat(r_data_spost,f_data_spost);

%     R_data_word = vertcat(r_data_w, R_data_word);
%     F_data_word = vertcat(f_data_w, F_data_word);
%     R_data_scene = vertcat(r_data_s, R_data_scene);
%     F_data_scene = vertcat(f_data_s, F_data_scene);
    


    %compute correlation for every trial between word and scene

    
    
    %correlate trials word and scene


        [R_correl, R_corrP] = corr(r_data_wpost',r_data_spost');
        [F_correl, F_corrP] = corr(f_data_wpost',f_data_spost');
        [all_correl, all_corrP] = corr(all_data_wpost',all_data_spost');
        
        results_tradeoff.Rcorrel(s,1) = R_correl;
        results_tradeoff.Fcorrel(s,1)  = F_correl;
        results_tradeoff.allcorrel(s,1)  = all_correl;
        results_tradeoff.RcorrelP(s,1) = R_corrP;
        results_tradeoff.FcorrelP(s,1)  = F_corrP;
        results_tradeoff.allcorrelP(s,1)  = all_corrP;
end
        
%Fisher's transform 
        
%FISHER Equation
%z = 0.5*log((1+r)/(1-r));

for x = 1:length(results_tradeoff.Rcorrel)

        results_tradeoff.Rcorrel(x,2) = 0.5*log((1+results_tradeoff.Rcorrel(x,1))/(1-results_tradeoff.Rcorrel(x,1)));
        results_tradeoff.Fcorrel(x,2) = 0.5*log((1+results_tradeoff.Fcorrel(x,1))/(1-results_tradeoff.Fcorrel(x,1)));
        results_tradeoff.allcorrel(x,2) = 0.5*log((1+results_tradeoff.allcorrel(x,1))/(1-results_tradeoff.allcorrel(x,1)));
end
        
        results_tradeoff.R_correl_mean_z = mean(results_tradeoff.Rcorrel(:,2));
        results_tradeoff.F_correl_mean_z = mean(results_tradeoff.Fcorrel(:,2));
        results_tradeoff.all_correl_mean_z = mean(results_tradeoff.allcorrel(:,2));
        
        %ttest for correlation
        
        chance = zeros(24,1);
        [results_tradeoff.all_ttest_chance,results_tradeoff.all_ttest_chanceP, results_tradeoff.all_ttest_chanceCI, ...
            results_tradeoff.all_ttest_chanceSTATS] = ttest(results_tradeoff.allcorrel(:,2),chance);
        
        [results_tradeoff.R_ttest_chance,results_tradeoff.R_ttest_chanceP, results_tradeoff.R_ttest_chanceCI, ...
            results_tradeoff.R_ttest_chanceSTATS] = ttest(results_tradeoff.Rcorrel(:,2),chance);
        
        [results_tradeoff.F_ttest_chance,results_tradeoff.F_ttest_chanceP, results_tradeoff.F_ttest_chanceCI, ...
            results_tradeoff.aF_ttest_chanceSTATS] = ttest(results_tradeoff.Fcorrel(:,2),chance);
        
        [results_tradeoff.RF_ttest,results_tradeoff.RF_ttestP, results_tradeoff.RF_ttestCI, ...
            results_tradeoff.RF_ttestSTATS] = ttest(results_tradeoff.Fcorrel(:,2),results_tradeoff.Rcorrel(:,2));
        
       
    
outfilename = 'emodif_tradeoff_ana.mat';

save(outfilename,'results_tradeoff');
    
    
end




