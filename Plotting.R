# This script plots all of the figures in the manuscript using ggplot

PlotDir=paste0(getwd(),'/Figures/')
source(paste0(getwd(), '/DrugAnnotationMOA.R'))
library(scales)
library(cowplot)
library(dplyr)
library(ggpubr)
library(boot)
library(tidyr)
library(ggridges)


BootAS = function(df, NumBoot=100) {
    out = data.frame()
    for (i in 1:NumBoot) {
        rep = data.frame(df %>% group_by(Bin, AS, Threshold) %>% sample_frac(1, replace=TRUE) %>% 
                    summarize(retained = sum(True), not_retained = sum(False)) %>% 
                    group_by(Bin,AS, Threshold) %>% 
                    summarize(fraction = retained/(retained + not_retained)))
        rep$num = i
        out = rbind(out, rep)
    }
    return(out)
}

VisualizeAllAS = function(Threshold='0.8') {
    Dir =paste0(getwd(), '/Data/AS_Tables/')
    AS_Types = c('RI','ES','AD','AP','AT','AA','ME')
    Combined = data.frame()
    df = data.frame()
    for (AS in AS_Types) {
        tmp = read.table(paste0(Dir, 'TCGA_',AS,'_Counts_ThresholdByPSI_', Threshold), sep=',', header=TRUE)
        tmp$AS = AS
        df = rbind(tmp, df)
    }
    df$Bin = cut(df$MutLoad,  breaks=c(0,50,100,500,1000,50000), 
                        labels=c('0-50','50-100','100-500','500-1000','1000->10,000'))
    df = na.omit(df); df$Threshold=Threshold
    df$FractionNMD= df$True/(df$False + df$True)
    booted= BootAS(df)
    booted_CI = data.frame(booted %>% group_by(Bin, AS) %>% summarise(mean = mean(fraction, na.rm = TRUE), sd = sd(fraction, na.rm = TRUE), n= n()) %>%
                mutate(se = sd / sqrt(n), lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se, upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se))
    booted_CI$Bin2 = factor(booted_CI$Bin, levels=unique(booted_CI$Bin))
    PlotOut = ggplot(booted_CI, aes(x=Bin2, y=mean, group=AS, color=AS)) + 
        geom_point() + theme_minimal() + geom_line() + #ylim(0.08,0.1) +
        geom_errorbar(aes(ymin=upper.ci, ymax=lower.ci), width=.2) +
        labs(y=expression(atop(underline("Under-expressed Transcripts With AS Event"), 
             paste("All Transcripts With AS Event"))), x= 'Number of Protein Coding Mutations', color='')
    ggsave(paste0(PlotDir, 'AS_RI_AllASType_TCGA.pdf' ), width=6, height=5, units='in') 

}


VisualizeAllASThresholds = function(AS_Type='RI') { 
    Dir = paste0(getwd(), '/Data/AS_Tables/')
    Thresholds = c('0.5','0.6','0.7','0.8','0.9')
    Combined = data.frame()
    for (T in Thresholds) {
        df = read.table(paste0(Dir, 'TCGA_',AS_Type,'_Counts_ThresholdByPSI_', T), sep=',', header=TRUE)
        df$Threshold = T
        Combined = rbind(Combined, df)
    }
    Combined$Bin = cut(Combined$MutLoad,  breaks=c(0,50,100,500,1000,50000), 
                        labels=c('0-50','50-100','100-500','500-1000','1000->10,000'))
    Combined = na.omit(Combined); Combined$AS='RI'
    Combined$FractionNMD= Combined$True/(Combined$False + Combined$True)
    booted= BootAS(Combined)
    booted_CI = data.frame(booted %>% group_by(Bin, Threshold) %>% summarise(mean = mean(fraction, na.rm = TRUE), sd = sd(fraction, na.rm = TRUE), n= n()) %>%
                mutate(se = sd / sqrt(n), lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se, upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se))
    booted_CI$Bin2 = factor(booted_CI$Bin, levels=unique(booted_CI$Bin))
    PlotOut = ggplot(booted_CI, aes(x=Bin2, y=mean, group=Threshold, color=Threshold)) + 
        geom_point() + theme_minimal() + geom_line() + #ylim(0.08,0.1) +
        geom_errorbar(aes(ymin=upper.ci, ymax=lower.ci), width=.2) +
        labs(y=expression(atop(underline("Under-expressed Transcripts With Intron Retention"), 
             paste("All Transcripts With Intron Retention"))), x= 'Number of Protein Coding Mutations')
    ggsave(paste0(PlotDir, 'AS_RI_AllThresholds_TCGA.pdf' ), width=6, height=5, units='in') 
}



VisualizeAS = function(df) {
    df$Bin = cut(df$MutLoad, breaks=c(0,50,100,500,1000,50000), 
                    labels=c('0-50','50-100','100-500','500-1000','1000->10,000'))
    df = na.omit(df); df$AS = 'RI'; df$Threshold = 0.8
    df$FractionNMD= df$True/(df$False + df$True)
    foo = df %>%  group_by(Bin) %>% sample_frac(1, replace=TRUE) %>% 
                    summarize(retained = sum(True), not_retained = sum(False)) %>% 
                    group_by(Bin) %>% 
                    summarize(fraction = retained/(retained + not_retained))
    booted= BootAS(df)

    booted_CI = data.frame(booted %>% group_by(Bin) %>% summarise(mean = mean(fraction, na.rm = TRUE), sd = sd(fraction, na.rm = TRUE), n= n()) %>%
                mutate(se = sd / sqrt(n), lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se, upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se))
    booted_CI$Bin2 = factor(booted_CI$Bin, levels=booted_CI$Bin )
    PlotOut = ggplot(booted_CI, aes(x=Bin2, y=mean, group=1)) + 
        geom_point() + theme_minimal() + geom_line() + #ylim(0.08,0.1) +
        geom_ribbon(aes(ymin=upper.ci, ymax=lower.ci), alpha=.7, position=position_dodge(.9)) +
        labs(y=expression(atop(underline("Under-expressed Transcripts With Intron Retention"), paste("All Transcripts With Intron Retention"))),
             x= 'Number of Protein Coding Mutations')
    if (df$Filtered == 'Filtered') {
         ggsave(paste0(PlotDir, 'AS_RI_PerTumor_PSI_TCGA.pdf' ), width=4.5, height=4.5, units='in') 
    } else {
         ggsave(paste0(PlotDir, 'AS_RI_PerTumor_PSI_TCGA_NoFilterForeQTL.pdf' ), width=5, height=4.5, units='in') 
    }

}



PlotCircularCORUMTCGA = function(df) { 

    df$negative_log10_of_adjusted_p_value = -log10(df$p_value)
    df= subset(df, df$source == 'CORUM')
    # Match CORUM terms to broader categories
    df$grouping = ''
    df[grep('ribo',df[,1]),]$grouping = 'Translation' 
    df[grep('roteasome|PA',df[,1]),]$grouping = 'Protein Degradation' 
    df[grep('CCT',df[,1]),]$grouping = 'Chaperone'
    df[grep('synthesome|BRCA|FA|BLM|MCM|RC|BRAFT|DNA-PK',df[,1]),]$grouping = 'DNA Replication/Repair'
    df[grep('snRNP|plice|SMN|Rnase|Sm|Exosome',df[,1]),]$grouping = 'Transcription'
    df[grep('RC|CEN|Nup',df[,1]),]$grouping = 'DNA Replication/Repair'
    # Add colors
    colors = data.frame(
        grouping = c('Protein Degradation','Translation','Transcription','Chaperone','DNA Replication/Repair', 'Immune','Other' ),
        colors = as.character(c('#F7B530','#3EA612','#A3E189','#F73036','#82CEF5','#3366C7','#F97FA2'))
    )
    ColorsForGrouping=colors$grouping
    df = merge(df, colors, by='grouping')
    df$term_name = gsub("\\(.*","",df$term_name) # remove parenthesis for visualization and repeat categories
    df= df[!duplicated(df$term_name),] # Remove duplicate enrichment categoris
    df$term_name = gsub("ribosome,","ribosome,\n",df$term_name) # add line break in term names for visualization
    df$term_name = gsub("20S-PA28","20S-PA28\n",df$term_name) # add line break in term names for visualization
    df$term_name = gsub("of cell cycle","",df$term_name) # add line break in term names for visualization
    GapForLegend = 1
    GroupsForOrdering = colors$grouping
    df = df %>% 
        arrange(grouping, -negative_log10_of_adjusted_p_value) %>%
        arrange(factor(grouping, levels = ColorsForGrouping)) %>%
        mutate(term_name = factor(term_name, levels = unique(term_name) )) %>%
        add_row(term_name = rep(NA,GapForLegend), negative_log10_of_adjusted_p_value = rep(NA,GapForLegend)) %>%
        mutate(ID = 1:(nrow(df) + GapForLegend)) %>%
        mutate(angle = 90-360*(ID-0.5)/(nrow(df) + 1)) %>%
        mutate(hjust = ifelse(angle < -90, 1, 0)) %>%
        mutate(angle = ifelse(angle < -90, angle+180, angle)) %>%  
        mutate(ID = factor(ID, levels = unique(ID)))

    PlotOut = ggplot(df, aes(x=ID, y=negative_log10_of_adjusted_p_value, fill=colors)) +
        geom_segment(x = 70, y = 10, xend = 1, yend = 10, colour = "grey", alpha=0.05, size=0.3 ) + # this set the scale line for y = 10
        geom_segment(x = 70, y = 5, xend = 1, yend = 5, colour = "grey", alpha=0.05, size=0.3 ) + # this set the scale line for y = 5
        geom_bar(stat='identity') +
        ylim(0,20)+
        labs(x='', y='Negative Log10 of Adjusted P-Value') + 
        #scale_fill_brewer(palette="Paired") +
        scale_fill_identity(guide = "legend", 
                            labels = unique(df$grouping),
                            breaks = unique(df$colors))  +
        theme_minimal() +
        theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            legend.position='none',
            plot.margin = unit(rep(-1,4), "cm") ) +  # Adjust the margin to make in sort labels are not truncated!
        #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_blank())+
        #theme(strip.text.y = element_blank(), legend.position='bottom') +
        #guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
        geom_text(aes(label = term_name, y = ifelse(negative_log10_of_adjusted_p_value<4,5, negative_log10_of_adjusted_p_value+ GapForLegend), 
            hjust = hjust, angle = angle), size = 8)+
          annotate("text", x = rep(nrow(df),2), y = c(5,10), label = c("5","10") , color="grey", size=8 , alpha=1.6, angle=0, fontface="bold", hjust=0.25, vjust=2) + # Add annotation
        coord_polar() 
        print(head(df))
    ggsave(paste0(PlotDir, 'CircularRegressionGeneSetEnrichmentGlobalCORUMTCGA.pdf' ), width=15, height=15, units='in') 

}

PlotCircularKeggTCGA = function(df) { 

    df = subset(df, df$source == 'KEGG')
    df = df[!grepl("disease", df$term_name),]
    df = df[!df$term_name %in% c('Amyotrophic lateral sclerosis','Epstein-Barr virus infection','Human T-cell leukemia virus 1 infection',
                           'Viral myocarditis','Influenza A','Spinocerebellar ataxia','Legionellosis'),]
    df$grouping = ''
    df$term_name = gsub('Protein processing in endoplasmic reticulum','Protein processing in ER',df$term_name)
    df$term_name = gsub('Antigen processing and presentation','Antigen processing\nand presentation',df$term_name)
    df[grep('Ribo|ER|amino',df[,1]),]$grouping = 'Translation' 
    df[grep('Proteasome|Phagosome',df[,1]),]$grouping = 'Protein Degradation' 
    df[grep('RNA|Splic|mRNA',df[,1]),]$grouping = 'Transcription' 
    df[grep('repair|recombination|Cell|replication|sugar|backbone',df[,1]),]$grouping = 'DNA Replication/Repair'
    df[grep('Apop|senesc|Amyo|Eps|Hum|myoc|Inf|Spin|Leg|cofac|Carb|disease',df[,1]),]$grouping = 'Other'
    df[grep('Antigen|killer',df[,1]),]$grouping = 'Immune'
    colors = data.frame(
        grouping = c('Protein Degradation','Translation','Transcription','Chaperone','DNA Replication/Repair', 'Immune','Other' ),
        colors = as.character(c('#F7B530','#3EA612','#A3E189','#F73036','#82CEF5','#3366C7','#F97FA2'))
    )
    ColorsForGrouping=colors$grouping
    df = merge(df, colors, by='grouping', all.y=TRUE)
    print(df)
    df$negative_log10_of_adjusted_p_value = -log10(df$p_value)
    rankingDF = df %>% group_by(grouping) %>% summarise(max = max(negative_log10_of_adjusted_p_value))
    rankingDF = rankingDF[order(rankingDF$max),]
    rankingDF$rank = nrow(rankingDF):1
    df = merge(df, rankingDF, by='grouping')
    df = df[order(df$rank, -df$negative_log10_of_adjusted_p_value), ]
    df$ID = 1:nrow(df)
    df$ID2 = factor(df$ID, levels=df$ID)
    df$grouping = factor(df$grouping, levels=unique(df$grouping))
    df$term_name = gsub("\\(.*","",df$term_name) # remove parenthesis for visualization and repeat categories
    df= df[!duplicated(df$term_name),] # Remove duplicate hits
    df = df[order(df$ID),]
    GapForLegend = 1
   
    df = df %>% 
        arrange(grouping, -negative_log10_of_adjusted_p_value) %>%
        arrange(factor(grouping, levels = ColorsForGrouping)) %>%
        mutate(individual = factor(term_name, levels = unique(term_name))) %>%
        add_row(term_name = rep(NA,GapForLegend), negative_log10_of_adjusted_p_value = rep(NA,GapForLegend)) %>%
        mutate(ID = 1:(nrow(df) + GapForLegend)) %>%
        mutate(angle = 90-360*(ID-0.25)/(nrow(df) + 1)) %>%
        mutate(hjust = ifelse(angle < -90, 1, 0)) %>%
        mutate(angle = ifelse(angle < -90, angle+180, angle)) %>%
        mutate(ID = factor(ID, levels = unique(ID)))
  
    
    PlotOut = ggplot(df, aes(x=ID, y=negative_log10_of_adjusted_p_value, fill=colors)) +
        geom_segment(x = 70, y = 20, xend = 1, yend = 20, colour = "grey", alpha=0.05, size=0.3 ) + # this set the scale line for y = 20
        geom_segment(x = 70, y = 15, xend = 1, yend = 15, colour = "grey", alpha=0.05, size=0.3 ) + # this set the scale line for y = 15
        geom_segment(x = 70, y = 10, xend = 1, yend = 10, colour = "grey", alpha=0.05, size=0.3 ) + # this set the scale line for y = 10
        geom_segment(x = 70, y = 5, xend = 1, yend = 5, colour = "grey", alpha=0.05, size=0.3 ) + # this set the scale line for y = 5
        geom_bar(stat='identity') +
        labs(x='', y='Negative Log10 of Adjusted P-Value') + 
        #scale_fill_brewer(palette="Paired") +
        scale_fill_identity(guide = "legend", labels = unique(df$grouping),breaks = unique(df$colors))  +
        theme_minimal() +
        theme(axis.text = element_blank(),axis.title = element_blank(), panel.grid = element_blank(),
            legend.text=element_text(size=17),
            legend.position='bottom', plot.margin = unit(rep(-1,4), "cm"), legend.title = element_blank()) +  
        guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
        geom_text(aes(label = term_name, y = ifelse(negative_log10_of_adjusted_p_value< 6,5, negative_log10_of_adjusted_p_value+ GapForLegend), 
            hjust = hjust, angle = angle), size = 8)+ 
            annotate("text", x = rep(nrow(df),4), y = c(5,10,15,20), label = c("5","10","15","20") , 
            color="grey", size=8 , alpha=1.6, angle=0, fontface="bold", hjust=0.25, vjust=2) +
        coord_polar() 
    ggsave(paste0(PlotDir, 'CircularRegressionGeneSetEnrichmentGlobalKEGGTCGA.pdf' ), width=15, height=18, units='in') 

}



PlotRegCoefPerGroup = function(df) { 
   
    df$Dataset = gsub('TCGA', 'TCGA (Human Tumors)', df$Dataset); df$Dataset = gsub('CCLE', 'CCLE (Cancer Cell Lines)', df$Dataset)
    df$Dataset2 = factor(df$Dataset, levels=c('TCGA (Human Tumors)','CCLE (Cancer Cell Lines)'))
    quantile = subset(df, df$subgroup == 'Quantile')
    quantile$Label = paste0(as.numeric(as.character(quantile$Group))*100, '%')
    df = subset(df, df$subgroup != 'Quantile')
    df$AdjPval = p.adjust(df$pval, method= 'fdr')
    df$NegLog10Pval = -log10(as.numeric(as.character(df$AdjPval)))
    df$SortedLevel = factor(df$Group, levels=c('Mitochondrial Chaperones', 'ER Chaperones','Small HS','HSP 100','HSP 90',
                'HSP 70','HSP 60', 'HSP 40','20S Core','20S Catalytic Core','Immunoproteasome Core','Immunoproteasome Catalytic Core',
                 '11S Regulatory Particle', '19S Regulatory Particle','Mitochondrial Ribosomes', 'Cytoplasmic Ribosomes'))
        
    print(head(quantile))
    Estimate = ggplot(data=df, aes(y = SortedLevel, x=as.numeric(as.character(Estimate)), color=subgroup)) +  
                geom_vline(data = quantile, aes(xintercept = as.numeric(as.character(Estimate))), alpha=0.6, color='grey', size=1) +
                geom_boxplot() +
                theme_minimal() + labs(x='Effect Size (Beta Coefficient)', y='') + theme(legend.position='none') +
                facet_wrap(~Dataset2, scales='free_x') +
                scale_color_manual(values=c('#8856a7','#F7B530','#3EA612')) + 
                #geom_vline(data=filter(quantile, Group == '0.05'), aes(xintercept=Estimate))+
                #geom_vline(xintercept=0, linetype='dashed', col = 'black') + 
                geom_text(data = quantile, aes(xintercept = as.numeric(as.character(Estimate)), 
                    label=Label, y=12), colour="black",  vjust="inward", size=3) +
                theme(strip.text = element_text(face="bold", size=12),  panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
                theme(panel.spacing = unit(2, "lines"), legend.title=element_blank())
                
    Rank = ggplot(data=df, aes(y = SortedLevel, x=as.numeric(as.character(NegLog10Pval)), color=subgroup)) + 
                geom_vline(xintercept= -log10(0.05), alpha=0.6, size=1, col = 'grey') + facet_wrap(~Dataset2, scales='free_x') +
                geom_boxplot() +
                theme_minimal() + labs(x='Negative Log10 of Adjusted P-Value', y='')  + theme(legend.position='bottom') +
                facet_wrap(~Dataset2, scales='free_x') +
                theme(legend.position='bottom', legend.title=element_blank()) +
                scale_color_manual(values=c('#8856a7','#F7B530','#3EA612')) + 
                theme(strip.text = element_text(face="bold", size=12))

    Combined = plot_grid(Estimate, Rank, rel_heights=c(0.88, 1), ncol = 1) 
    ggsave(paste0(PlotDir, 'RegCoefPerGroups_TCGAandCCCLE.pdf' ), width=8, height=8, units='in')


}



PlotDeltaPSI = function(df) {
    df$term_name  = gsub(', ATP synthesis by chemiosmotic coupling, and heat production by uncoupling proteins.', "", df$term_name)  
    df$negative_log10_of_adjusted_p_value = -log10(df$p_value)
    df$Group = gsub('PosPSI','Less Intron Retention\n(In High vs. Low)', df$Group)
    df$Group = gsub('NegPSI','More Intron Retention\n(In High vs. Low)', df$Group)
    rankingDF = df %>% group_by(Group) %>% summarise(max = max(negative_log10_of_adjusted_p_value))
    rankingDF = rankingDF[order(rankingDF$max),]
    rankingDF$rank = nrow(rankingDF):1
    df = merge(df, rankingDF, by='Group')
    df = df[order(df$rank, -df$negative_log10_of_adjusted_p_value), ]
    df$ID = 1:nrow(df)
    df$ID2 = factor(df$ID, levels=df$ID)
    df$grouping = factor(df$Group, levels=unique(df$Group))
    df$term_name = gsub("\\(.*","",df$term_name) # remove parenthesis for visualization and repeat categories
    df= df[!duplicated(df$term_name),] # Remove duplicate hits
    Out = ggplot(df, aes(x=reorder(term_name,-ID), y=negative_log10_of_adjusted_p_value, fill=source)) +
        geom_col() +
        coord_flip() +
        facet_grid(rows=vars(Group), scales='free_y', space='free', switch='y') +
        labs(x='', y='Negative Log10 of Adjusted P-Value') + 
        scale_fill_manual(values = c('red','#3279D8')) +
        theme_minimal() + theme(legend.position='bottom', legend.title = element_blank(), legend.title.align=0.5) 
    ggsave(paste0(PlotDir, 'AS_Delta_PSI_IntonRetention_TCGA.pdf' ), width=6.5, height=4.5, units='in')
}



PlotWithinCancerTypeExpressionAcrossGroups = function(df, dataset) { 
    print(head(df))
    df$AdjPval = p.adjust(df$Pr...t.., method= 'fdr')
    df$NegLog10Pval = -log10(as.numeric(as.character(df$AdjPval)))
    df$SortedLevel = factor(df$Group, levels=c('Mitochondrial Chaperones', 'ER Chaperones','Small HS','HSP 100','HSP 90',
                'HSP 70','HSP 60', 'HSP 40','20S Core','20S Catalytic Core','Immunoproteasome Core','Immunoproteasome Catalytic Core',
                 '11S Regulatory Particle', '19S Regulatory Particle','Mitochondrial Ribosomes', 'Cytoplasmic Ribosomes'))
    Estimate = ggplot(data=df, aes(x = CancerType, y=as.numeric(as.character(Estimate)), color=SortedLevel)) +  geom_boxplot() +
                facet_wrap(~SortedLevel, ncol=5) +
                 theme_minimal() + labs(y='Effect Size (Beta Coefficient)', x='') +  theme(legend.position='none') +
                 geom_hline(yintercept=0, linetype='dashed', col = 'black') + #facet_wrap(~SortedLevel, scales='free_x', ncol=12) +
                 theme(strip.text = element_text(face="bold", size=8)) + ggtitle(dataset) +
                 theme(plot.title = element_text(hjust = 0.5)) + # Center title
                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7))
    ggsave(paste0(PlotDir, 'RegCoefPerGroupsWithinCancerType_', dataset,'.pdf' ), width=10, height=6, units='in')

}




PlotJacknifedExpressionAcrossGroups = function(df, dataset) { 
    print(head(df))
    df$AdjPval = p.adjust(df$Pr...t.., method= 'fdr')
    df$NegLog10Pval = -log10(as.numeric(as.character(df$AdjPval)))
    df$SortedLevel = factor(df$Group, levels=c('Mitochondrial Chaperones', 'ER Chaperones','Small HS','HSP 100','HSP 90',
                'HSP 70','HSP 60', 'HSP 40','20S Core','20S Catalytic Core','Immunoproteasome Core','Immunoproteasome Catalytic Core',
                 '11S Regulatory Particle', '19S Regulatory Particle','Mitochondrial Ribosomes', 'Cytoplasmic Ribosomes'))
    Estimate = ggplot(data=df, aes(x = CancerTypeRemoved, y=as.numeric(as.character(Estimate)), color=SortedLevel)) +  geom_boxplot() +
                facet_wrap(~SortedLevel, ncol=5) +
                 theme_minimal() + labs(y='Effect Size (Beta Coefficient)', x='') +  theme(legend.position='none') +
                 geom_hline(yintercept=0, linetype='dashed', col = 'black') + #facet_wrap(~SortedLevel, scales='free_x', ncol=12) +
                 theme(strip.text = element_text(face="bold", size=8)) + ggtitle(dataset) +
                 theme(plot.title = element_text(hjust = 0.5)) + # Center title
                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7))
    ggsave(paste0(PlotDir, 'RegCoefPerGroupsJacknifed_', dataset,'.pdf' ), width=10, height=6, units='in')

}

PlotJacknifedDrugsAcrossGroups = function(df) { 
    df$AdjPval = p.adjust(df$Pr...t.., method= 'fdr')
    df$Sig = ifelse(df$Pr...t.. < 0.05, '*','')
    df$Sig = ifelse(df$Pr...t.. < 0.005, '**', df$Sig)
    df$Sig = ifelse(df$Pr...t.. < 0.0005, '***', df$Sig)
    df$Dummy = '' #'Protein Coding\nMutations'
    df = subset(df, (df$subgroup != 'Protein Synthesis Inhibitor') & (df$subgroup != 'RNA Synthesis Inhibitor'))
    df$subgroup2 = gsub(' Inhibitor', '\nInhibitor', df$subgroup)
    df$subgroup2 = gsub('cific Pro', 'cific\nPro', df$subgroup)
    print(head(df))
    ggplot(df, aes(y=CancerTypeRemoved, x = name ,fill=as.numeric(as.character(Estimate)))) + geom_tile()+
        #scale_fill_viridis(discrete=FALSE) +
       # scale_fill_gradientn(colours=c("red","white","grey"), limits=c(1,-1)) +
        scale_fill_gradientn(colors=c("red","white","blue"), 
           values=rescale(c(low=-0.5, mid=0, high=0.1)), breaks=c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1),limits=c(-0.5,0.1))+
        geom_text(aes(label = Sig, color=Sig), show.legend = FALSE, size= 7, vjust = 0.85) +
        theme_minimal() +
        #coord_flip() +
        #facet_wrap(~subgroup, scales='free', ncol=5) +
        #facet_wrap(~Group, scales='free_x', space='free') +
        facet_grid(cols=vars(subgroup2), scales='free_x', space='free_x') +
        scale_colour_manual(values=c("black","black",'black','black')) +
        labs(y='',x='Drug Name', fill='Beta\nCoefficient') + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size=12)) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) + ggtitle('PRISM (Drug Screen)') 
    ggsave(paste0(PlotDir, 'CCLEDrugJacknifed.pdf' ), width=12, height=6, units='in')
}
    

VisualizeBetaAfterShuffle = function(df) {
    df$NegLog10Pval = -log10(as.numeric(as.character(df$Pr...t..)))
    Beta = ggplot(df, aes(x=Estimate, color=Group, fill=Group)) +
        geom_histogram() + theme_minimal() + labs(x='Beta Coefficient (TMB vs Expression)', y='Counts')
    Pval = ggplot(df, aes(x=NegLog10Pval, color=Group, fill=Group)) +
        geom_histogram() + theme_minimal() + labs(x='Negative Log10 P-Value', y='Counts')
    Combined = plot_grid(Beta, Pval, rel_heights=c(0.88, 1), ncol = 1) 
    ggsave(paste0(PlotDir, 'TMBShuffling_BetaDist_Expression_TCGA.pdf' ), width=4, height=4, units='in')
}

ComparePValueDistributions = function(df) {
    original = df[df$Group == 'Observed',]
    null = df[df$Group == 'Null',]
    CalcValues = function(x) { # used to calculate adjusted p values; how many values are greater than estimate under null?
        if (x > 0) {
            return(nrow(null[null$Estimate > x,])/nrow(null))
        } else {
            return(nrow(null[null$Estimate < x,])/nrow(null))
        }
    }
    original$adj_pval = sapply(original$Estimate,CalcValues)
    plotting = rbind(
        cbind(setNames(original[c('GeneName','adj_pval')], c('Gene','Pval')), data.frame(Group=c('Adjusted'))),
        cbind(setNames(original[c('GeneName','Pr...t..')], c('Gene','Pval')), data.frame(Group=c('Original')))
    )
    plotting$NegLog10Pval = -log10(plotting$Pval)
    Pval = ggplot(plotting, aes(x=NegLog10Pval, color=Group, fill=Group)) +
        geom_histogram() + theme_minimal() + labs(x='Negative Log10 P-Value', y='Counts') +
        scale_y_log10() + facet_wrap(~Group, scales="free_x")
    ggsave(paste0(PlotDir, 'TMBShuffling_AdjustedPval_TCGA.pdf' ), width=4, height=4, units='in')

}

PlotVarExprByCancerType = function() {
    df = read.table(paste0(getwd(),'/Data/Regression/TCGAVarianceInGeneExpression'), sep=',', header=T)
    grps = df[df['Group'] != '',]
    grps$Group = gsub(" ","\n", grps$Group)

    ggplot(grps, aes(y=type, x=Value, color=Group, fill=Group)) + 
        geom_boxplot() + theme_minimal() +  theme(legend.position="none") +
        facet_grid(cols=vars(Group), scales='free_x') + labs(y='',x='Variance In Gene Expression') +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7)) +
        theme(strip.text.y = element_text(size = 7)) + 
        theme(strip.text.x = element_text(size = 7))
    
    ggsave(paste0(PlotDir, 'VarianceInExpAcrossCancerType_TCGA.pdf' ), width=8, height=6, units='in')

}


# PlotAllCorumNullvsObs = function(df) {
#     original = df[df$Group == 'Observed',]
#     null = df[df$Group == 'Null',]
#     CalcValues = function(x) { # used to calculate adjusted p values; how many values are greater than estimate under null?
#         if (x > 0) {
#             return(nrow(null[null$Estimate > x,])/nrow(null))
#         } else {
#             return(nrow(null[null$Estimate < x,])/nrow(null))
#         }
#     }
#     original$adj_pval = sapply(original$Estimate,CalcValues)


# }


PlotCorrelationCoefficientsBetweenExpAndProt = function(df) {

    df$Exp_Estimate = as.numeric(as.character(df$Exp_Estimate))
    df$Protein_Estimate = as.numeric(as.character(df$Protein_Estimate))
    # print(head(df))
    # cors = df %>% group_by(Group, subgroup) %>% summarise(
    #     cor = cor.test(Protein_Estimate, Exp_Estimate)$estimate,
    #     pval = cor.test(Protein_Estimate, Exp_Estimate)$p.value)

 
    cors = df %>% group_by(Group, subgroup) %>% summarise(
        abs = abs(Protein_Estimate - Exp_Estimate))

    print(head(cors))

    PlotEstimate = ggplot(cors, aes(x=Group, y=abs, color=subgroup, fill=subgroup)) + 
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_blank()) +
    facet_wrap(~subgroup, scales='free_x') + labs(x='', y="Difference Between Protein and Expression ") 
    #scale_y_continuous(breaks = seq(-1, 1, by = 0.2))
    
    #PlotPval = ggplot(cors, aes(x=Group, y=pval)) + geom_point()
    #Combined = plot_grid(PlotEstimate, PlotPval, rel_heights=c(0.88, 1), ncol = 1) 

    ggsave(paste0(PlotDir,'CCLE_Protein_vs_Exp_Correlation_BoxPlot.pdf' ), width=8, height=5, units='in')

}



PlotAdjustedPValGeneSetEnrich = function(df) {
    
    original = df[df$Group == 'Observed',]
    null = df[df$Group == 'Null',]
    CalcValues = function(x) { # used to calculate adjusted p values; how many values are greater than estimate under null?
        if (x > 0) {
            return(nrow(null[null$Estimate > x,])/nrow(null))
        } else {
            return(nrow(null[null$Estimate < x,])/nrow(null))
        }
    }
    original$adj_pval = sapply(original$Estimate,CalcValues)
    library('gprofiler2')
    result = data.frame(gost(original[original$adj_pval  < 0.05,]$GeneName)$result)
    df= subset(result, result$source == 'CORUM')[c('source','term_name','p_value')]
    # Match CORUM terms to broader categories
    df$grouping = ''
    #df[grep('ribo',df[,1]),]$grouping = 'Translation' 
    df[grep('roteasome|PA',df[,2]),]$grouping = 'Protein Degradation' 
    df[grep('CCT',df[,2]),]$grouping = 'Chaperone'
    df[grep('synthesome|BRCA|FA|BLM|MCM|RC|BRAFT|DNA-PK|RFC2|PCNA|GINS|MCC',df[,2]),]$grouping = 'DNA Replication/Repair'
    #df[grep('snRNP|plice|SMN|Rnase|Sm|Exosome',df[,2]),]$grouping = 'Transcription'
    df[grep('Conden|kinetochore',df[,2]),]$grouping = 'Chromosome Segregation'
    df[grep('RC|CEN|Nup',df[,2]),]$grouping = 'DNA Replication/Repair'
    # Add colors
    colors = data.frame(
        grouping = c('Protein Degradation','Translation','Transcription','Chaperone','DNA Replication/Repair', 'Immune','Other' ),
        colors = as.character(c('#F7B530','#3EA612','#A3E189','#F73036','#82CEF5','#3366C7','#F97FA2'))
    )
    ColorsForGrouping=colors$grouping
    df = merge(df, colors, by='grouping')
    df$negative_log10_of_adjusted_p_value = -log10(df$p_value)
    rankingDF = df %>% group_by(grouping) %>% summarise(max = max(negative_log10_of_adjusted_p_value))
    rankingDF = rankingDF[order(rankingDF$max),]
    rankingDF$rank = nrow(rankingDF):1
    df = merge(df, rankingDF, by='grouping')
    df = df[order(df$rank, -df$negative_log10_of_adjusted_p_value), ]
    df$ID = 1:nrow(df)
    df$ID2 = factor(df$ID, levels=df$ID)
    df$grouping = factor(df$grouping, levels=unique(df$grouping))
    df$term_name = gsub("\\(.*","",df$term_name) # remove parenthesis for visualization and repeat categories
    df= df[!duplicated(df$term_name),] # Remove duplicate hits
    print(head(df))
    PlotOut = ggplot(df, aes(x=reorder(term_name,-ID), y=negative_log10_of_adjusted_p_value, fill=colors)) +
        geom_col() +
        coord_flip() +
        facet_grid(rows=vars(grouping), scales='free_y', space='free', switch='y') +
        labs(x='', y='Negative Log10 of Adjusted P-Value') + 
        #scale_fill_brewer(palette="Paired") +
        scale_fill_identity(guide = "legend", 
                            labels = unique(df$grouping),
                            breaks = unique(df$colors))  +
        # scale_fill_manual(values = df$colors, labels=df$grouping) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.title = element_blank())+
        theme(strip.text.y = element_blank(), legend.position='bottom') +
        guides(fill=guide_legend(nrow=2,byrow=TRUE))
    ggsave(paste0(PlotDir, 'TCGA_AdjustedPVal_CORUM_GSE.pdf' ), width=7, height=5, units='in') 
}


ProteinVsExpressionByIndividualGene = function(df, Dataset) {

    print(head(df))
    quantile = subset(df, df$subgroup == 'Quantile')
    quantile$Label = paste0(as.numeric(as.character(quantile$Group))*100, '%')
    quantile = subset(quantile, (quantile$Label != '0%') & (quantile$Label != '100%'))
    df = subset(df, df$subgroup != 'Quantile')
    df$NegLog10Pval = -log10(as.numeric(as.character(df$AdjPval)))
    df$SortedLevel = factor(df$Group, levels=c('Mitochondrial Chaperones', 'ER Chaperones','Small HS','HSP 100','HSP 90',
                'HSP 70','HSP 60', 'HSP 40','20S Core','20S Catalytic Core','Immunoproteasome Core','Immunoproteasome Catalytic Core',
                 '11S Regulatory Particle', '19S Regulatory Particle','Mitochondrial Ribosomes', 'Cytoplasmic Ribosomes'))
        
    print(head(quantile))
    Estimate = ggplot(data=df, aes(y = SortedLevel, x=as.numeric(as.character(Estimate)), color=subgroup)) +  
                geom_vline(data = quantile, aes(xintercept = as.numeric(as.character(Estimate))), alpha=0.6, color='grey', size=1) +
                geom_boxplot() +
                theme_minimal() + labs(x='Effect Size (Beta Coefficient)', y='') + theme(legend.position='none') +
                facet_wrap(~DataType, scales='free_x') +
                scale_color_manual(values=c('#8856a7','#F7B530','#3EA612')) + 
                #geom_vline(data=filter(quantile, Group == '0.05'), aes(xintercept=Estimate))+
                #geom_vline(xintercept=0, linetype='dashed', col = 'black') + 
                geom_text(data = quantile, aes(xintercept = as.numeric(as.character(Estimate)), 
                    label=Label, y=12), colour="black",  vjust="inward", size=3) +
                theme(strip.text = element_text(face="bold", size=12),  panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
                theme(panel.spacing = unit(2, "lines"), legend.title=element_blank())
                
    Rank = ggplot(data=df, aes(y = SortedLevel, x=as.numeric(as.character(NegLog10Pval)), color=subgroup)) + 
                geom_vline(xintercept= -log10(0.05), alpha=0.6, size=1, col = 'grey') + 
                geom_boxplot() +
                theme_minimal() + labs(x='Negative Log10 of Adjusted P-Value', y='')  + theme(legend.position='bottom') +
                facet_wrap(~DataType, scales='free_x') +
                theme(legend.position='bottom', legend.title=element_blank()) +
                scale_color_manual(values=c('#8856a7','#F7B530','#3EA612')) + 
                theme(strip.text = element_text(face="bold", size=12))


    Combined = plot_grid(Estimate, Rank, rel_heights=c(0.88, 1), ncol = 1) 
    ggsave(paste0(PlotDir, Dataset,'_Protein_vs_Exp_BoxPlot.pdf' ), width=8, height=8, units='in')


}

PlotProteinAcrossCancerTypesInGroups = function(df) { 
    df = na.omit(df)
    df = subset(df, df$type =='All')
    df$AdjPval = p.adjust(df$Pr...t.., method= 'fdr')
    df$NumSamples = as.numeric(as.character(df$NumSamples))
    df = subset(df, df$NumSamples > 10)
    df$NegLog10Pval = -log10(as.numeric(as.character(df$AdjPval)))
    df$SortedLevel = factor(df$Group, levels=c('Mitochondrial Chaperones', 'ER Chaperones','Small HS','HSP 100','HSP 90',
                'HSP 70','HSP 60', 'HSP 40','20S Core','20S Catalytic Core','Immunoproteasome Core','Immunoproteasome Catalytic Core',
                 '11S Regulatory Particle', '19S Regulatory Particle','Mitochondrial Ribosomes', 'Cytoplasmic Ribosomes'))
    Estimate = ggplot(data=df, aes(x = type, y=as.numeric(as.character(Estimate)), color=SortedLevel)) +  geom_boxplot() +
                facet_wrap(~SortedLevel, ncol=5, scales='free_y') +
                 theme_minimal() + labs(y='Effect Size (Beta Coefficient)', x='') +  theme(legend.position='none') +
                 geom_hline(yintercept=0, linetype='dashed', col = 'black') + #facet_wrap(~SortedLevel, scales='free_x', ncol=12) +
                 theme(strip.text = element_text(face="bold", size=8)) + 
                 theme(plot.title = element_text(hjust = 0.5)) + # Center title
                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7))
    ggsave(paste0(PlotDir, 'CCLE_ProteinGrouped_WithinType.pdf' ), width=12, height=8, units='in')
}


PlotCancerTypeInComplexesAcrossGroups = function(df) { 


    # data <- scale(t(data))
    # ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
    # ord

    print(head(df))
    df$NumSamples = as.numeric(as.character(df$NumSamples))
    df$MutLoad = as.numeric(as.character(df$MutLoad))
    #df = subset(df, df$NumSamples > 250)
    #df = subset(df, df$MutLoad > 25)
    df$AdjPval = p.adjust(df$Pr...t.., method= 'fdr')
    df$Sig = ifelse(df$Pr...t.. < 0.05, '*','')
    df$Sig = ifelse(df$Pr...t.. < 0.005, '**', df$Sig)
    df$Sig = ifelse(df$Pr...t.. < 0.0005, '***', df$Sig)
    df$Dummy = '' #'Protein Coding\nMutations'
    df$Group = gsub('mic Rib', 'mic\nRib', df$Group); df$Group = gsub('ial Rib', 'ial\nRib', df$Group); df$Group = gsub('ial Ch', 'ial\nCh', df$Group)
    df$Group = gsub('ory Pa', 'ory\nPa', df$Group)
    df$subgroup = gsub('Ribosomes','Translation',df$subgroup)
    df$Group2 = factor(df$Group, levels= c('HSP 100','HSP 90','HSP 70','HSP 60','Mitochondrial\nChaperones','ER Chaperones','Small HS',
                '19S Regulatory\nParticle','20S Core','Mitochondrial\nRibosomes','Cytoplasmic\nRibosomes'))
    ggplot(df, aes(y=type, x = Group2 ,fill=Estimate)) + geom_tile()+
        #scale_fill_viridis(discrete=FALSE) +
        #scale_fill_gradientn(colours=c("blue","black","red")) +
        scale_fill_gradientn(colors=c("red","white","blue"), 
            values=rescale(c(min(df$Estimate),0, max(df$Estimate))),limits=c(min(df$Estimate),max(df$Estimate)))+
        geom_text(aes(label = Sig, color=Sig), show.legend = FALSE, size= 8, vjust = 0.85) +
        theme_minimal() +
        #coord_flip() +
        #facet_wrap(~subgroup, scales='free', ncol=5, space='free_x') +
        #facet_wrap(~Group, scales='free_x', space='free') +
        facet_grid(cols=vars(subgroup), scales='free_x', space='free_x') +
        scale_colour_manual(values=c("black","black",'black','black')) +
        labs(y='',x='', fill='Beta\nCoefficient') + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(text = element_text(size=18))
    ggsave(paste0(PlotDir, 'TCGA_ExpressionGrouped_WithinType.pdf' ), width=12, height=8, units='in')
}


PlotJacknifedshRNAAcrossGroups = function(df) { 

    df$AdjPval = p.adjust(df$Pr...t.., method= 'fdr')
    df$Sig = ifelse(df$Pr...t.. < 0.05, '*','')
    df$Sig = ifelse(df$Pr...t.. < 0.005, '**', df$Sig)
    df$Sig = ifelse(df$Pr...t.. < 0.0005, '***', df$Sig)
    df$Dummy = '' #'Protein Coding\nMutations'
    df$Group = gsub('mic Rib', 'mic\nRib', df$Group); df$Group = gsub('ial Rib', 'ial\nRib', df$Group); df$Group = gsub('ial Ch', 'ial\nCh', df$Group)
    df$Group = gsub('ory Pa', 'ory\nPa', df$Group)
    df$subgroup = gsub('Ribosomes','Translation',df$subgroup)
    df$Group2 = factor(df$Group, levels= c('HSP 100','HSP 90','HSP 70','HSP 60','Mitochondrial\nChaperones','ER Chaperones','Small HS',
                '19S Regulatory\nParticle','20S Core','Mitochondrial\nRibosomes','Cytoplasmic\nRibosomes'))
    ggplot(df, aes(y=CancerTypeRemoved, x = Group2 ,fill=Estimate)) + geom_tile()+
        #scale_fill_viridis(discrete=FALSE) +
        #scale_fill_gradientn(colours=c("blue","black","red")) +
        scale_fill_gradientn(colors=c("red","white","blue"), 
            values=rescale(c(min(df$Estimate),0, max(df$Estimate))),limits=c(min(df$Estimate),max(df$Estimate)))+
        geom_text(aes(label = Sig, color=Sig), show.legend = FALSE, size= 8, vjust = 0.85) +
        theme_minimal() +
        #coord_flip() +
        #facet_wrap(~subgroup, scales='free', ncol=5, space='free_x') +
        #facet_wrap(~Group, scales='free_x', space='free') +
        facet_grid(cols=vars(subgroup), scales='free_x', space='free_x') +
        scale_colour_manual(values=c("black","black",'black','black')) +
        labs(y='',x='', fill='Beta\nCoefficient') + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(text = element_text(size=18))+
        ggtitle('Achilles (shRNA Screen)')
    ggsave(paste0(PlotDir, 'CCLEshRNAJacknifed.pdf' ), width=12, height=8, units='in')
}



PlotFractionOfSigDrugsForAll = function(df) { 
    df$log10pval= log10(as.numeric(as.character(df$pVal)))
    df$SigGroups = ifelse(df$pVal < 0.05, 'Significant', 'NotSignificant')
    df$EstimateGroups = ifelse(df$Estimate > 0, 'Increase In\nViability\nWith Load','Decrease In\nViability\nWith Load')
    #print(head(df))
    PlotOfCounts = ggplot(df, aes(x = EstimateGroups, color=SigGroups, fill=SigGroups)) + geom_bar() +
        theme_minimal() + labs(y='Number of Drugs', x='') + scale_colour_manual(values=c('grey','black')) +
        scale_fill_manual(values=c('grey','black')) +
        theme(legend.position='bottom', legend.title=element_blank())# +
       # guides(fill=guide_legend(nrow=2)) + guides(colour = guide_legend(nrow = 2))
    ggsave(paste0(PlotDir, 'RegCoefAllDrugs_CCLE.pdf' ), width=3, height=5, units='in')
}

PlotBootstrappedNegativelyAssociatedDrugs = function() { 
    library(viridis)
    df = GetAnnotatedMOA()
    df$log10pval= log10(as.numeric(as.character(df$pVal)))
    df$SigGroups = ifelse(df$pVal < 0.05, 'Significant', 'NotSignificant')
    df$EstimateGroups = ifelse(df$Estimate > 0, 'Increase In Viability With Load','Decrease In Viability With Load')
    CountsPerDrug = data.frame(table(as.character(df$PlottingGroup)))
    df = merge(df, CountsPerDrug, by.x='PlottingGroup', by.y='Var1') 
    df = subset(df, (df$Estimate < 0))
    # Bootstrap 50 drugs from each broad group 100 times; How often is the group significant? 
    set.seed(123)
    AllBootstraps = data.frame() # Empty df where to append results
    for (i in 1:100) {
        # Sample within broad groups
        SampledByGroup = df %>% group_by(PlottingGroup) %>% sample_n(size = 50, replace=TRUE) %>% # Sample 50 drugs
                group_by(PlottingGroup, SigGroups) %>% tally() %>% # Count how many times each drug is sig assoc w/ load
                spread(SigGroups, n) %>% replace(is.na(.), 0) %>% # Long to wide transformation of significant groups
                mutate(FracSig = Significant/(NotSignificant + Significant)) # Calculate fraction significant drugs in group
        # Sample from all drugs (aka not by groups)
        SampledRandom = df %>% sample_n(size = 50, replace=TRUE) %>% 
                group_by(SigGroups) %>% tally() %>% spread(SigGroups, n) %>% replace(is.na(.), 0) %>% 
                mutate(FracSig = Significant/(NotSignificant + Significant)) 
        SampledRandom = data.frame('PlottingGroup' = c('Any Category'), SampledRandom)
        SampledByGroup = rbind(data.frame(SampledByGroup), SampledRandom)
        AllBootstraps = rbind(AllBootstraps, data.frame(SampledByGroup))
    }
    # Order results by largest effect size
    Rank = data.frame(AllBootstraps %>% group_by(PlottingGroup) %>% summarise(mean = mean(FracSig)))
    Rank$FracSigRank = rank(Rank$mean)
    Rank = Rank[order(Rank$FracSigRank),]
    AllBootstraps$PlottingGroup2 = factor(AllBootstraps$PlottingGroup, levels=as.character(Rank$PlottingGroup))
    GroupsWithZero = as.character(data.frame(AllBootstraps %>% group_by(PlottingGroup) %>% summarise(med=median(FracSig)) %>% filter(med == 0))$PlottingGroup)
    AllBootstraps = AllBootstraps[!(AllBootstraps$PlottingGroup %in% GroupsWithZero),] # Remove categories with 0 sig assoc with load
    # Box plot of all boostraps
    PlotOfCounts = ggplot(AllBootstraps, aes(y = PlottingGroup2, x=FracSig, color=PlottingGroup2, fill=PlottingGroup2)) + 
        geom_boxplot() + scale_color_viridis(discrete=TRUE) +  scale_fill_viridis(discrete=TRUE) +
        theme_minimal() + labs(y='', x='Fraction of Significant Drugs In Category') + 
        theme(legend.position='none', legend.title=element_blank()) +
        geom_vline(xintercept=mean(subset(AllBootstraps, AllBootstraps$PlottingGroup == 'Any Category')$FracSig), linetype="dashed")
    ggsave(paste0(PlotDir, 'RegCoefBootstrappedByDrugGroupSize_CCLE.pdf' ), width=6, height=5, units='in')

}



PlotMultiCollinearityOfSNVsAndCNVs = function(df) { 
    library('ggcorrplot')
    print(head(df))
    df$CNA = df$CNV; df$CNV = NULL # change cnvs = cnas to be consistent with text
    ccle = subset(df, df$Dataset == 'CCLE'); row.names(ccle) = ccle$Barcode; ccle$Dataset=NULL; ccle$Barcode=NULL
    tcga = subset(df, df$Dataset == 'TCGA'); row.names(tcga) = tcga$Barcode; tcga$Dataset=NULL; tcga$Barcode=NULL
    ccle.cor = cor(ccle, method = c("pearson")); tcga.cor = cor(tcga, method = c("pearson"))
    ccle_plot = ggcorrplot(ccle.cor, lab = TRUE) + ggtitle('CCLE') + theme(plot.title = element_text(hjust = 0.5,face="bold", size=12))
    tcga_plot = ggcorrplot(tcga.cor, lab = TRUE) + ggtitle('TCGA') + theme(plot.title = element_text(hjust = 0.5,face="bold", size=12))
    plot_grid(tcga_plot, ccle_plot, labels = c("A", "B"), ncol=1)
    ggsave(paste0(PlotDir, 'CCLE_TCGA_Multicollinearity.pdf' ), width=7, height=11, units='in')
}


PlotGLMMRegressionCoefficientsByAge = function(df) { 
    print(head(df))
    ggplot(df, aes(x=AgeBin, y=Estimate, fill=Group)) + 
        geom_boxplot() + facet_wrap(~Group) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        labs(x='Age Groups', y='Association Between Mutational Load and Expression (Beta Coefficient)')
    ggsave(paste0(PlotDir, 'TCGA_GLMM_RegressionByAge.pdf' ), width=11, height=8, units='in')

}



PlotNumASEventsFiltered = function(df) { 
    df$Group = gsub('_',' ', df$Group)
    df = melt(df)
    print(head(df))
    PlotOut = ggplot(df, aes(x=Group, y=value , fill=variable)) + geom_bar(stat='identity', position = "dodge") + theme_minimal() +
        labs(x='', y='Number of Alternative Splicing Events') +  theme(legend.title = element_blank()) 
        #scale_y_log10("y", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))
    ggsave(paste0(PlotDir, 'TCGA_CountsFilteredASEvents.pdf' ), width=5, height=5, units='in')


}