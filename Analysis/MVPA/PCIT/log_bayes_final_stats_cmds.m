load(log_bayes_R_RF_F_violin.mat)
log_bayes_R_RF_F_violin_mean = mean(log_bayes_R_RF_F_violin,1);
SD_exclusion = (log_bayes_R_RF_F_violin_mean*3)+log_bayes_R_RF_F_violin_mean;
SD_exclusion_ID_R = log_bayes_R_RF_F_violin(:,1) > SD_exclusion(:,1);
SD_exclusion_ID_RF = log_bayes_R_RF_F_violin(:,2) > SD_exclusion(:,2);
SD_exclusion_ID_F = log_bayes_R_RF_F_violin(:,3) > SD_exclusion(:,3);
SD_exclusion_ID = horzcat(SD_exclusion_ID_R,SD_exclusion_ID_RF, SD_exclusion_ID_F);
log_bayes_R_RF_F_violin_clean = log_bayes_R_RF_F_violin;
log_bayes_R_RF_F_violin_clean(SD_exclusion_ID) = NaN;
violin(log_bayes_R_RF_F_violin_clean, 'xlabel', {'Remember', 'Remember+Forget', 'Forget'}, 'facecolor',[.7 .7 .7;.2 0 0;.9 .9 0],'mc','k','medc','')

%for jsut remember and forget

log_bayes_R_F_violin_clean = horzcat(log_bayes_R_RF_F_violin_clean(:,1),log_bayes_R_RF_F_violin_clean(:,3));
violin(log_bayes_R_F_violin_clean, 'xlabel', {'Remember', 'Forget'},'mc','k','medc','')




