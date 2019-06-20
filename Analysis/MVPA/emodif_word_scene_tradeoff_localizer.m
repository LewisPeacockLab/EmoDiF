function emodif_word_scene_tradeoff_localizer()

%to where results directory is
%load data
load('emodif_aggregate_localizer_corr.mat')

%aggregate corr data
%also have fisherz as second value
%FISHER Equation
%z = 0.5*log((1+r)/(1-r));


for s = 1:24
    corr.all_ws(s,1) = results.bysubject.data(s).localizer.corr.corr_all_ws;
    corr.all_ws(s,2) =  0.5*log((1+corr.all_ws(s,1))/(1-corr.all_ws(s,1)));
    corr.all_wf(s,1) = results.bysubject.data(s).localizer.corr.corr_all_wf;
    corr.all_wf(s,2) =  0.5*log((1+corr.all_wf(s,1))/(1-corr.all_wf(s,1)));
    
    
    corr.face_ws(s,1) = results.bysubject.data(s).localizer.corr.corr_face_ws;
    corr.face_ws(s,2) =  0.5*log((1+corr.face_ws(s,1))/(1-corr.face_ws(s,1)));
    corr.face_wf(s,1) = results.bysubject.data(s).localizer.corr.corr_face_wf;
    corr.face_wf(s,2) =  0.5*log((1+corr.face_wf(s,1))/(1-corr.face_wf(s,1)));
    
    corr.scene_ws(s,1) = results.bysubject.data(s).localizer.corr.corr_scene_ws;
    corr.scene_ws(s,2) =  0.5*log((1+corr.scene_ws(s,1))/(1-corr.scene_ws(s,1)));
    corr.scene_wf(s,1) = results.bysubject.data(s).localizer.corr.corr_scene_wf;
    corr.scene_wf(s,2) =  0.5*log((1+corr.scene_wf(s,1))/(1-corr.scene_wf(s,1)));
    
    corr.word_ws(s,1) = results.bysubject.data(s).localizer.corr.corr_word_ws;
    corr.word_ws(s,2) =  0.5*log((1+corr.word_ws(s,1))/(1-corr.word_ws(s,1)));
    corr.word_wf(s,1) = results.bysubject.data(s).localizer.corr.corr_word_wf;
    corr.word_wf(s,2) =  0.5*log((1+corr.word_wf(s,1))/(1-corr.word_wf(s,1)));
    
end

%means
results.corr = corr;
results.corr_m.all_ws = mean(corr.all_ws(:,2));
results.corr_m.all_wf = mean(corr.all_wf(:,2));
results.corr_m.face_ws = mean(corr.face_ws(:,2));
results.corr_m.face_wf = mean(corr.face_wf(:,2));
results.corr_m.scene_ws = mean(corr.scene_ws(:,2));
results.corr_m.scene_wf = mean(corr.scene_wf(:,2));
results.corr_m.word_ws = mean(corr.word_ws(:,2));
results.corr_m.word_wf = mean(corr.word_wf(:,2));

%ttests

chance = zeros(24,1);

[results.corr_ttest.all_ws,results.corr_ttest.all_wsP, results.corr_ttest.all_wsCI, ...
    results.corr_ttest.all_wsSTATS] = ttest(corr.all_ws(:,2),chance);

[results.corr_ttest.all_wf,results.corr_ttest.all_wfP, results.corr_ttest.all_wfCI, ...
    results.corr_ttest.all_wfSTATS] = ttest(corr.all_wf(:,2),chance);



[results.corr_ttest.face_ws,results.corr_ttest.face_wsP, results.corr_ttest.face_wsCI, ...
    results.corr_ttest.face_wsSTATS] = ttest(corr.face_ws(:,2),chance);

[results.corr_ttest.face_wf,results.corr_ttest.face_wfP, results.corr_ttest.face_wfCI, ...
    results.corr_ttest.face_wfSTATS] = ttest(corr.face_wf(:,2),chance);



[results.corr_ttest.scene_ws,results.corr_ttest.scene_wsP, results.corr_ttest.scene_wsCI, ...
    results.corr_ttest.scene_wsSTATS] = ttest(corr.scene_ws(:,2),chance);

[results.corr_ttest.scene_wf,results.corr_ttest.scene_wfP, results.corr_ttest.scene_wfCI, ...
    results.corr_ttest.scene_wfSTATS] = ttest(corr.scene_wf(:,2),chance);



[results.corr_ttest.word_ws,results.corr_ttest.word_wsP, results.corr_ttest.word_wsCI, ...
    results.corr_ttest.word_wsSTATS] = ttest(corr.word_ws(:,2),chance);

[results.corr_ttest.word_wf,results.corr_ttest.word_wfP, results.corr_ttest.word_wfCI, ...
    results.corr_ttest.word_wfSTATS] = ttest(corr.word_wf(:,2),chance);

outfilename = 'emodif_tradeoff_localizer.mat';

save(outfilename,'results');



