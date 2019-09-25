%emodif shaded error bars

x=1:7;
shadedErrorBar(x,evidence_curves.fsor.remember_scenem,evidence_curves.fsor.RFdiff_sceneSE,'lineprops','o');
hold on
shadedErrorBar(x,evidence_curves.fsor.forget_scenem,evidence_curves.fsor.RFdiff_sceneSE,'lineprops','y');

x=1:7;
shadedErrorBar(x,evidence_curves.fsowr.remember_scenem,evidence_curves.fsowr.RFdiff_sceneSE,'lineprops','r');
hold on
shadedErrorBar(x,evidence_curves.fsowr.forget_scenem,evidence_curves.fsowr.RFdiff_sceneSE,'lineprops','o');

x=1:7;
shadedErrorBar(x,evidence_curves.fsowr.remember_wordm,evidence_curves.fsowr.RFdiff_wordSE,'lineprops','r');
hold on
shadedErrorBar(x,evidence_curves.fsowr.forget_wordm,evidence_curves.fsowr.RFdiff_wordSE,'lineprops','o');

x=1:7;
shadedErrorBar(x,evidence_curves.fsowr.remember_facem,evidence_curves.fsowr.RFdiff_faceSE,'lineprops','r');
hold on
shadedErrorBar(x,evidence_curves.fsowr.forget_facem,evidence_curves.fsowr.RFdiff_faceSE,'lineprops','o');

x=1:7;
shadedErrorBar(x,evidence_curves.fsowr20.remember_scenem,evidence_curves.fsowr20.RFdiff_sceneSE,'lineprops','r');
hold on
shadedErrorBar(x,evidence_curves.fsowr20.forget_scenem,evidence_curves.fsowr20.RFdiff_sceneSE,'lineprops','o');

x=1:7;
shadedErrorBar(x,evidence_curves.fsowr20.remember_wordm,evidence_curves.fsowr20.RFdiff_wordSE,'lineprops','r');
hold on
shadedErrorBar(x,evidence_curves.fsowr20.forget_wordm,evidence_curves.fsowr20.RFdiff_wordSE,'lineprops','o');

x=1:7;
shadedErrorBar(x,evidence_curves.fsowr20.remember_facem,evidence_curves.fsowr20.RFdiff_faceSE,'lineprops','r');
hold on
shadedErrorBar(x,evidence_curves.fsowr20.forget_facem,evidence_curves.fsowr20.RFdiff_faceSE,'lineprops','o');