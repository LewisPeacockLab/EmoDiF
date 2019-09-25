
function MVPA_create_evidence_curves()

%create evidence curves.  this is to create figure 3 and also supplementary
%figures for word out analysis and pilot out analysis. 

%note for figure 3, neural separation SEMS will be used.

%load complete data

load('emodif_aggregate_results_complete');

emodif_fsowr.remember_scene = [];
emodif_fsowr.remember_word = [];
emodif_fsowr.remember_face = [];

emodif_fsowr.forget_scene = [];
emodif_fsowr.forget_word = [];
emodif_fsowr.forget_face = [];

for s = 1:24
emodif_fsowr.remember_wordo = results.bysubject.data(s).results.remember.all.mean.word(:)';
emodif_fsowr.remember_wordb = emodif_fsowr.remember_wordo(1,:);
emodif_fsowr.remember_wordbb = emodif_fsowr.remember_wordb(1,:)-emodif_fsowr.remember_wordb(1,1);
emodif_fsowr.remember_word = vertcat(emodif_fsowr.remember_word, emodif_fsowr.remember_wordbb);

emodif_fsowr.remember_sceneo = results.bysubject.data(s).results.remember.all.mean.scene(:)';
emodif_fsowr.remember_sceneb = emodif_fsowr.remember_sceneo(1,:);
emodif_fsowr.remember_scenebb = emodif_fsowr.remember_sceneb(1,:)-emodif_fsowr.remember_sceneb(1,1);
emodif_fsowr.remember_scene = vertcat(emodif_fsowr.remember_scene, emodif_fsowr.remember_scenebb);


emodif_fsowr.remember_faceo = results.bysubject.data(s).results.remember.all.mean.face(:)';
emodif_fsowr.remember_faceb = emodif_fsowr.remember_faceo(1,:);
emodif_fsowr.remember_facebb = emodif_fsowr.remember_faceb(1,:)-emodif_fsowr.remember_faceb(1,1);
emodif_fsowr.remember_face = vertcat(emodif_fsowr.remember_face, emodif_fsowr.remember_facebb);

emodif_fsowr.forget_wordo = results.bysubject.data(s).results.forget.all.mean.word(:)';
emodif_fsowr.forget_wordb = emodif_fsowr.forget_wordo(1,:);
emodif_fsowr.forget_wordbb = emodif_fsowr.forget_wordb(1,:)-emodif_fsowr.forget_wordb(1,1);
emodif_fsowr.forget_word = vertcat(emodif_fsowr.forget_word, emodif_fsowr.forget_wordbb);

emodif_fsowr.forget_sceneo = results.bysubject.data(s).results.forget.all.mean.scene(:)';
emodif_fsowr.forget_sceneb = emodif_fsowr.forget_sceneo(1,:);
emodif_fsowr.forget_scenebb = emodif_fsowr.forget_sceneb(1,:)-emodif_fsowr.forget_sceneb(1,1);
emodif_fsowr.forget_scene = vertcat(emodif_fsowr.forget_scene, emodif_fsowr.forget_scenebb);


emodif_fsowr.forget_faceo = results.bysubject.data(s).results.forget.all.mean.face(:)';
emodif_fsowr.forget_faceb = emodif_fsowr.forget_faceo(1,:);
emodif_fsowr.forget_facebb = emodif_fsowr.forget_faceb(1,:)-emodif_fsowr.forget_faceb(1,1);
emodif_fsowr.forget_face = vertcat(emodif_fsowr.forget_face, emodif_fsowr.forget_facebb);

end

%mean

emodif_fsowr.remember_wordm = mean(emodif_fsowr.remember_word);
emodif_fsowr.remember_scenem = mean(emodif_fsowr.remember_scene);
emodif_fsowr.remember_facem = mean(emodif_fsowr.remember_face);

emodif_fsowr.forget_wordm = mean(emodif_fsowr.forget_word);
emodif_fsowr.forget_scenem = mean(emodif_fsowr.forget_scene);
emodif_fsowr.forget_facem = mean(emodif_fsowr.forget_face);

%SEM of the difference

%create difference

emodif_fsowr.RFdiff_word = emodif_fsowr.remember_word - emodif_fsowr.forget_word;
emodif_fsowr.RFdiff_scene = emodif_fsowr.remember_scene - emodif_fsowr.forget_scene;
emodif_fsowr.RFdiff_face = emodif_fsowr.remember_face - emodif_fsowr.forget_face;

%create SD

emodif_fsowr.RFdiff_wordSD = std(emodif_fsowr.RFdiff_word);
emodif_fsowr.RFdiff_sceneSD = std(emodif_fsowr.RFdiff_scene);
emodif_fsowr.RFdiff_faceSD = std(emodif_fsowr.RFdiff_face);

%create SE

emodif_fsowr.RFdiff_wordSE = emodif_fsowr.RFdiff_wordSD/sqrt(length(emodif_fsowr.RFdiff_word ));
emodif_fsowr.RFdiff_sceneSE = emodif_fsowr.RFdiff_sceneSD/sqrt(length(emodif_fsowr.RFdiff_scene));
emodif_fsowr.RFdiff_faceSE = emodif_fsowr.RFdiff_faceSD/sqrt(length(emodif_fsowr.RFdiff_face ));

%complete - pilot data is just this without the first 4 subjects

load('emodif_aggregate_results_complete');

emodif_fsowr20.remember_scene = [];
emodif_fsowr20.remember_word = [];
emodif_fsowr20.remember_face = [];

emodif_fsowr20.forget_scene = [];
emodif_fsowr20.forget_word = [];
emodif_fsowr20.forget_face = [];

for s = 5:24

emodif_fsowr20.remember_wordo = results.bysubject.data(s).results.remember.all.mean.word(:)';
emodif_fsowr20.remember_wordb = emodif_fsowr20.remember_wordo(1,:);
emodif_fsowr20.remember_wordbb = emodif_fsowr20.remember_wordb(1,:)-emodif_fsowr20.remember_wordb(1,1);
emodif_fsowr20.remember_word = vertcat(emodif_fsowr20.remember_word, emodif_fsowr20.remember_wordbb);

emodif_fsowr20.remember_sceneo = results.bysubject.data(s).results.remember.all.mean.scene(:)';
emodif_fsowr20.remember_sceneb = emodif_fsowr20.remember_sceneo(1,:);
emodif_fsowr20.remember_scenebb = emodif_fsowr20.remember_sceneb(1,:)-emodif_fsowr20.remember_sceneb(1,1);
emodif_fsowr20.remember_scene = vertcat(emodif_fsowr20.remember_scene, emodif_fsowr20.remember_scenebb);


emodif_fsowr20.remember_faceo = results.bysubject.data(s).results.remember.all.mean.face(:)';
emodif_fsowr20.remember_faceb = emodif_fsowr20.remember_faceo(1,:);
emodif_fsowr20.remember_facebb = emodif_fsowr20.remember_faceb(1,:)-emodif_fsowr20.remember_faceb(1,1);
emodif_fsowr20.remember_face = vertcat(emodif_fsowr20.remember_face, emodif_fsowr20.remember_facebb);

emodif_fsowr20.forget_wordo = results.bysubject.data(s).results.forget.all.mean.word(:)';
emodif_fsowr20.forget_wordb = emodif_fsowr20.forget_wordo(1,:);
emodif_fsowr20.forget_wordbb = emodif_fsowr20.forget_wordb(1,:)-emodif_fsowr20.forget_wordb(1,1);
emodif_fsowr20.forget_word = vertcat(emodif_fsowr20.forget_word, emodif_fsowr20.forget_wordbb);

emodif_fsowr20.forget_sceneo = results.bysubject.data(s).results.forget.all.mean.scene(:)';
emodif_fsowr20.forget_sceneb = emodif_fsowr20.forget_sceneo(1,:);
emodif_fsowr20.forget_scenebb = emodif_fsowr20.forget_sceneb(1,:)-emodif_fsowr20.forget_sceneb(1,1);
emodif_fsowr20.forget_scene = vertcat(emodif_fsowr20.forget_scene, emodif_fsowr20.forget_scenebb);


emodif_fsowr20.forget_faceo = results.bysubject.data(s).results.forget.all.mean.face(:)';
emodif_fsowr20.forget_faceb = emodif_fsowr20.forget_faceo(1,:);
emodif_fsowr20.forget_facebb = emodif_fsowr20.forget_faceb(1,:)-emodif_fsowr20.forget_faceb(1,1);
emodif_fsowr20.forget_face = vertcat(emodif_fsowr20.forget_face, emodif_fsowr20.forget_facebb);

end

%mean

emodif_fsowr20.remember_wordm = mean(emodif_fsowr20.remember_word);
emodif_fsowr20.remember_scenem = mean(emodif_fsowr20.remember_scene);
emodif_fsowr20.remember_facem = mean(emodif_fsowr20.remember_face);

emodif_fsowr20.forget_wordm = mean(emodif_fsowr20.forget_word);
emodif_fsowr20.forget_scenem = mean(emodif_fsowr20.forget_scene);
emodif_fsowr20.forget_facem = mean(emodif_fsowr20.forget_face);

%SEM of the difference

%create difference

emodif_fsowr20.RFdiff_word = emodif_fsowr20.remember_word - emodif_fsowr20.forget_word;
emodif_fsowr20.RFdiff_scene = emodif_fsowr20.remember_scene - emodif_fsowr20.forget_scene;
emodif_fsowr20.RFdiff_face = emodif_fsowr20.remember_face - emodif_fsowr20.forget_face;

%create SD

emodif_fsowr20.RFdiff_wordSD = std(emodif_fsowr20.RFdiff_word);
emodif_fsowr20.RFdiff_sceneSD = std(emodif_fsowr20.RFdiff_scene);
emodif_fsowr20.RFdiff_faceSD = std(emodif_fsowr20.RFdiff_face);

%create SE

emodif_fsowr20.RFdiff_wordSE = emodif_fsowr20.RFdiff_wordSD/sqrt(length(emodif_fsowr20.RFdiff_word ));
emodif_fsowr20.RFdiff_sceneSE = emodif_fsowr20.RFdiff_sceneSD/sqrt(length(emodif_fsowr20.RFdiff_scene));
emodif_fsowr20.RFdiff_faceSE = emodif_fsowr20.RFdiff_faceSD/sqrt(length(emodif_fsowr20.RFdiff_face ));


%for fsor data

load('emodif_aggregate_results_fsor_complete');

emodif_fsor.remember_scene = [];
emodif_fsor.remember_face = [];

emodif_fsor.forget_scene = [];
emodif_fsor.forget_face = [];

for s = 1:24
for s = 1:24

emodif_fsor.remember_sceneo = results.bysubject.data(s).results.remember.all.mean.scene(:)';
emodif_fsor.remember_sceneb = emodif_fsor.remember_sceneo(1,:);
emodif_fsor.remember_scenebb = emodif_fsor.remember_sceneb(1,:)-emodif_fsor.remember_sceneb(1,1);
emodif_fsor.remember_scene = vertcat(emodif_fsor.remember_scene, emodif_fsor.remember_scenebb);


emodif_fsor.remember_faceo = results.bysubject.data(s).results.remember.all.mean.face(:)';
emodif_fsor.remember_faceb = emodif_fsor.remember_faceo(1,:);
emodif_fsor.remember_facebb = emodif_fsor.remember_faceb(1,:)-emodif_fsor.remember_faceb(1,1);
emodif_fsor.remember_face = vertcat(emodif_fsor.remember_face, emodif_fsor.remember_facebb);


emodif_fsor.forget_sceneo = results.bysubject.data(s).results.forget.all.mean.scene(:)';
emodif_fsor.forget_sceneb = emodif_fsor.forget_sceneo(1,:);
emodif_fsor.forget_scenebb = emodif_fsor.forget_sceneb(1,:)-emodif_fsor.forget_sceneb(1,1);
emodif_fsor.forget_scene = vertcat(emodif_fsor.forget_scene, emodif_fsor.forget_scenebb);


emodif_fsor.forget_faceo = results.bysubject.data(s).results.forget.all.mean.face(:)';
emodif_fsor.forget_faceb = emodif_fsor.forget_faceo(1,:);
emodif_fsor.forget_facebb = emodif_fsor.forget_faceb(1,:)-emodif_fsor.forget_faceb(1,1);
emodif_fsor.forget_face = vertcat(emodif_fsor.forget_face, emodif_fsor.forget_facebb);

end

%mean

emodif_fsor.remember_scenem = mean(emodif_fsor.remember_scene);
emodif_fsor.remember_facem = mean(emodif_fsor.remember_face);

emodif_fsor.forget_scenem = mean(emodif_fsor.forget_scene);
emodif_fsor.forget_facem = mean(emodif_fsor.forget_face);


%create difference

emodif_fsor.RFdiff_scene = emodif_fsor.remember_scene - emodif_fsor.forget_scene;
emodif_fsor.RFdiff_face = emodif_fsor.remember_face - emodif_fsor.forget_face;

%create SD

emodif_fsor.RFdiff_sceneSD = std(emodif_fsor.RFdiff_scene);
emodif_fsor.RFdiff_faceSD = std(emodif_fsor.RFdiff_face);

%create SE

emodif_fsor.RFdiff_sceneSE = emodif_fsor.RFdiff_sceneSD/sqrt(length(emodif_fsor.RFdiff_scene));
emodif_fsor.RFdiff_faceSE = emodif_fsor.RFdiff_faceSD/sqrt(length(emodif_fsor.RFdiff_face ));


evidence_curves.fsor = emodif_fsor;
evidence_curves.fsowr = emodif_fsowr;
evidence_curves.fsowr20 = emodif_fsowr20;

outfilename = 'evidence_curve_data.mat';
save(outfilename, 'evidence_curves');
end







