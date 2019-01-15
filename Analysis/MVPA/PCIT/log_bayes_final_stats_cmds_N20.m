
log_bayes_N20_F_R_mean = mean(log_bayes_N20_F_R,1);
SD_exclusion = (log_bayes_N20_F_R_mean*3)+log_bayes_N20_F_R_mean;
SD_exclusion_ID_F = log_bayes_N20_F_R(:,1) > SD_exclusion(:,1);
SD_exclusion_ID_R = log_bayes_N20_F_R(:,2) > SD_exclusion(:,2);

SD_exclusion_ID = horzcat(SD_exclusion_ID_F, SD_exclusion_ID_R);
log_bayes_N20_F_R_violin_clean = log_bayes_N20_F_R;
log_bayes_N20_F_R_violin_clean(SD_exclusion_ID) = NaN;
violin(log_bayes_N20_F_R_violin_clean, 'xlabel', {'Forget','Remember'},'mc','k','medc','')

%for jsut remember and forget

log_bayes_R_F_violin_clean = horzcat(log_bayes_N20_F_R_violin_clean(:,1),log_bayes_N20_F_R_violin_clean(:,3));
violin(log_bayes_R_F_violin_clean, 'xlabel', {'Remember', 'Forget'},'mc','k','medc','')




