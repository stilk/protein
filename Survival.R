#install.packages(c("survival", "survminer"))
library("survival")
library("survminer")
library("cowplot")
library("gridExtra")
library("ggplot2")
library("grid")
#library("scater")
OutDir = "/labs/ccurtis2/tilk/scripts/protein/Figures/"


GetForestPlot = function(df, Dataset) {
    if (Dataset == 'MMRF') {
        df$age = as.numeric(as.character(df$age_at_index))
        df$MutLoad = as.numeric(as.character(df$MutLoad))
        df$LogMuts = log10(as.numeric(as.character(df$MutLoad)) + 1)
        res.cox <- coxph(Surv(as.numeric(as.character(days_to_last_follow_up)), as.numeric(as.character(survival))) ~ sex + age + LogMuts + strata(iss_stage), data =  df)
    } else if (Dataset == 'TCGA') {
        df$age = as.numeric(as.character(df$age_at_initial_pathologic_diagnosis))
        df$MutLoad = as.numeric(as.character(df$MutLoad))
        df$LogMuts = log10(as.numeric(as.character(df$MutLoad)) + 1)
        df$OverExpressedTranscript = as.numeric(as.character(df$OverExpressedTranscript))
        print(head(df))
        res.cox <- coxph(Surv(as.numeric(as.character(OS.time)), as.numeric(as.character(OS))) ~ ProteasomeGroup * MutLoadGroup , data =  df)
    }
       print(summary(res.cox))
    PlotOut = ggforest(res.cox, data=df)
    print('foo')
    ggsave(paste0(OutDir, Dataset, "_ForestPlot.pdf"))
    #print(paste0('Plot saved to: ', OutDir))
}

PlotSurvivalCurve = function(df, Dataset) {
    if (Dataset == 'MMRF') {
        df$survival = as.numeric(as.character(df$survival))
        df$survival_time = as.numeric(as.character(df$days_to_last_follow_up)) 
        df$age = as.numeric(as.character(df$age_at_index))
        df$LogMuts = log10(as.numeric(as.character(df$MutLoad)))
        res.cox <- survfit(Surv(survival_time, survival) ~ MutLoadGroup, data = df)
        PlotSurvival = ggsurvplot(res.cox, conf.int = TRUE, pval=TRUE, data=df)
    } else if (Dataset == 'TCGA') {
        #df= subset(df, df$StrataGroup != '')
        #df = df[df$StrataGroup %in% c('LowPro_LowTMB','HighPro_LowTMB'),]
        #print(unique(df$StrataGroup))
        df = df[df$StrataGroup %in% c('LowPro_HighTMB','HighPro_HighTMB'),]
        print(survdiff(Surv(as.numeric(as.character(PFI.time)), as.numeric(as.character(PFI))) ~  StrataGroup, data = df))
        #df = df[df$StrataGroup %in% c('HighPro_MidTMB', 'LowPro_MidTMB'),]
        res.cox <- survfit(Surv(as.numeric(as.character(PFI.time)), as.numeric(as.character(PFI))) ~ StrataGroup, data = df)
        PlotSurvival = ggsurvplot(res.cox, title = unique(df$type_x), conf.int = TRUE, pval=TRUE,data=df)
        return(PlotSurvival)
    }
     #PlotTable = ggsurvplot(res.cox, risk.table.y.text.col = TRUE,risk.table = TRUE, data=df)
    #PlotOut = plot_grid(PlotSurvival, PlotTable, ncol = 1, rel_heights = c(0.66, 0.33))
    #PlotSurvival
    ggsave(paste0(OutDir, Dataset, "_KaplanMeierCurve.pdf"), width=8, height=6, units="in")

}


PlotSurvivalCurveByCancerTypeFacet = function(df) {
    df = df[df$StrataGroup %in% c('LowPro_HighTMB','HighPro_HighTMB'),]
    df$PFI.time = as.numeric(as.character(df$PFI.time))
    print(survdiff(Surv(as.numeric(as.character(PFI.time)), as.numeric(as.character(PFI))) ~  StrataGroup, data = df))
    res.cox <- survfit(Surv(as.numeric(as.character(PFI.time)), as.numeric(as.character(PFI))) ~ StrataGroup, data = df)
    PlotSurvival = ggsurvplot_facet(res.cox, facet.by="type_x", conf.int = TRUE, pval=TRUE,data=df)
    ggsave(paste0(OutDir, "TCGA_FacetedByCancerTypeKaplanMeierCurve.pdf"), width=12, height=12, units="in")
}
