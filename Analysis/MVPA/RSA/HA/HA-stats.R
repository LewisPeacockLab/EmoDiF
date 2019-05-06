######                        ######
######  Hyperalignment stats  ######
######                        ######


suppressPackageStartupMessages(library(nlme))
suppressPackageStartupMessages(library(multcomp))
suppressPackageStartupMessages(library(pwr))
suppressPackageStartupMessages(library(effsize))
suppressPackageStartupMessages(library(optparse))


option_list <- list( 
    make_option("--mask",default="vtc"),
    make_option("--ws",default=FALSE,action="store_true")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

MASK <- opt$mask
if (opt$ws) {
    HAorWS <- "ws"
} else {
    HAorWS <- "ha"
}


# hyperalignment results are in 2 files,
# on holding test accuracies and one hold permutation distributions
resdir <- '~/DBw/STUDY/RetroBlast/results'
if (MASK == "ips") {
    resdir <- sprintf('%s/%s',resdir,MASK)
}
fname <- sprintf('%s/fig2a.csv',resdir)
df <- read.csv(fname)

acc_fname <- sprintf('%s/within_del-%s-acc.csv',resdir,HAorWS)
df_acc <- read.csv(acc_fname)
null_fname <- sprintf('%s/within_del-%s-null.csv',resdir,HAorWS)
df_null <- read.csv(null_fname)


# if within-subjects rather than hyperalignment,
# then need to collapse and get mean acc per subject
if (opt$ws) {
    df_acc <- aggregate(accuracy~subj+cond, df_acc, mean)
    # rename subj to fold so df looks like that from hyperalignment
    names(df_acc)[names(df_acc)=="subj"] <- "fold"
}

##
## compare accuracies across conditions (d1/stay/switch)
##

anova_df <- df_acc[df_acc$cond!="blast",]
anova_df$cond <- factor(anova_df$cond)

lme_model <- lme(accuracy~cond, random=~1|fold, data=anova_df)
# print(summary(lme_model))
aov_model <- anova(lme_model)
cat('\n**** HA-CLF-ACC ANOVA RESULTS ****\n')
print(aov_model[2,,])

# power and effect size (partial eta squared) for anova
aov.fit <- aov(accuracy~cond+Error(fold), data=anova_df)
within_effects <- summary(aov.fit$"Within")
ss_cond <- within_effects[[1]]["cond","Sum Sq"]
ss_residuals <- within_effects[[1]]["Residuals","Sum Sq"]
part_etasq <- ss_cond / (ss_cond+ss_residuals)
cat(sprintf('\tANOVA partial eta squared = %f\n',part_etasq))


posthoc <- glht(lme_model,linfct=mcp(cond="Tukey"))
summ_posthoc <- summary(posthoc)
cat('\n**** HA-CLF-ACC 2WAY RESULTS ****\n')
print(summ_posthoc)

nsamples <- nrow(anova_df) / 3
d1 <- subset(anova_df,cond=="d1",accuracy,drop=TRUE)
stay <- subset(anova_df,cond=="d2-stay",accuracy,drop=TRUE)
switch <- subset(anova_df,cond=="d2-switch",accuracy,drop=TRUE)
cohd <- cohen.d(d1,stay,paired=TRUE)
pow <- pwr.t.test(n=nsamples,d=cohd$estimate,sig.level=0.05,type="paired")
cat(sprintf('\td1 vs stay: d=%f, pow=%f\n',cohd$estimate,pow$power))
cohd <- cohen.d(d1,switch,paired=TRUE)
pow <- pwr.t.test(n=nsamples,d=cohd$estimate,sig.level=0.05,type="paired")
cat(sprintf('\td1 vs switch: d=%f, pow=%f\n',cohd$estimate,pow$power))
cohd <- cohen.d(stay,switch,paired=TRUE)
pow <- pwr.t.test(n=nsamples,d=cohd$estimate,sig.level=0.05,type="paired")
cat(sprintf('\tstay vs switch: d=%f, pow=%f\n',cohd$estimate,pow$power))



##
## descriptive stats and permutation tests
## of clf accuracies across conditions (d1/stay/switch)
##

cat('\n**** HA-CLF-ACC PERMUTATION TESTS ****\n')

for (c in c("d1","blast","d2-stay","d2-switch")){
    vals <- subset(df_acc,cond==c,accuracy,drop=TRUE)
    vmean <- mean(vals)
    vsd <- sd(vals)

    if (opt$ws) {
        # get the permutation distrns of all subjs
        # and resample from them to get same count as hyperalignment
        null_vals <- subset(df_null,cond==c,accuracy,drop=TRUE)
        nsamples <- length(null_vals) / length(unique(df_acc$fold)) # fold = subj
        null_distrn <- sample(null_vals,nsamples,replace=TRUE)
    } else {
        null_distrn <- subset(df_null,cond==c,accuracy,drop=TRUE)
    }

    # get 95% CI of null distrn
    loCI <- quantile(null_distrn,.05)
    hiCI <- quantile(null_distrn,.95)

    # get proportion of null distribution higher than mean acc
    null_distrn <- append(null_distrn,1) # prevent p=0
    pval <- mean(null_distrn>=vmean)

    cat(sprintf('\t%s, M=%f, SD=%f, p=%f, nullCI=(%f,%f)\n',
        c,vmean,vsd,pval,loCI,hiCI))
}