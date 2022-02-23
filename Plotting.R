PlotDir='/labs/ccurtis2/tilk/scripts/protein/Figures/'
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
    Dir = '/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/'
    AS_Types = c('RI','ES','AD','AP','AT','AA','ME')
    Combined = data.frame()
    df = data.frame()
    for (AS in AS_Types) {
        tmp = read.table(paste0(Dir, 'TCGA_',AS,'_Counts_ThresholdByPSI_', Threshold), sep=',', header=TRUE)
        tmp$AS = AS
        df = rbind(tmp, df)
    }
    df$Bin = cut(df$MutLoad, breaks=c(0,10,100,1000,50000), labels=c('0-10','10-100','100-1000','1000->10,000'))
    df = na.omit(df); df$Threshold=Threshold
    df$FractionNMD= df$True/(df$False + df$True)
    booted= BootAS(df)
    booted_CI = data.frame(booted %>% group_by(Bin, AS) %>% summarise(mean = mean(fraction, na.rm = TRUE), sd = sd(fraction, na.rm = TRUE), n= n()) %>%
                mutate(se = sd / sqrt(n), lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se, upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se))
    PlotOut = ggplot(booted_CI, aes(x=as.character(Bin), y=mean, group=AS, color=AS)) + 
        geom_point() + theme_minimal() + geom_line() + #ylim(0.08,0.1) +
        geom_errorbar(aes(ymin=upper.ci, ymax=lower.ci), width=.2) +
        labs(y=expression(atop(underline("Under-expressed Transcripts With AS Event"), 
             paste("All Transcripts With AS Event"))), x= 'Number of Protein Coding Mutations', color='')
    ggsave(paste0(PlotDir, 'AS_RI_AllASType_TCGA.pdf' ), width=5, height=5, units='in') 

}






VisualizeAllASThresholds = function(AS_Type='RI') {
    Dir = '/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/'
    Thresholds = c('0.5','0.6','0.7','0.8','0.9')
    Combined = data.frame()
    for (T in Thresholds) {
        df = read.table(paste0(Dir, 'TCGA_',AS_Type,'_Counts_ThresholdByPSI_', T), sep=',', header=TRUE)
        df$Threshold = T
        Combined = rbind(Combined, df)
    }
    Combined$Bin = cut(Combined$MutLoad, breaks=c(0,10,100,1000,50000), labels=c('0-10','10-100','100-1000','1000->10,000'))
    Combined = na.omit(Combined); Combined$AS='RI'
    Combined$FractionNMD= Combined$True/(Combined$False + Combined$True)
    booted= BootAS(Combined)
    booted_CI = data.frame(booted %>% group_by(Bin, Threshold) %>% summarise(mean = mean(fraction, na.rm = TRUE), sd = sd(fraction, na.rm = TRUE), n= n()) %>%
                mutate(se = sd / sqrt(n), lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se, upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se))
    PlotOut = ggplot(booted_CI, aes(x=as.character(Bin), y=mean, group=Threshold, color=Threshold)) + 
        geom_point() + theme_minimal() + geom_line() + #ylim(0.08,0.1) +
        geom_errorbar(aes(ymin=upper.ci, ymax=lower.ci), width=.2) +
        labs(y=expression(atop(underline("Under-expressed Transcripts With Intron Retention"), 
             paste("All Transcripts With Intron Retention"))), x= 'Number of Protein Coding Mutations')
    ggsave(paste0(PlotDir, 'AS_RI_AllThresholds_TCGA.pdf' ), width=5, height=5, units='in') 
}



VisualizeAS = function() {
    Dir = '/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/'
    df = read.table('/labs/ccurtis2/tilk/scripts/protein/Data/AS_Tables/TCGA_RI_Counts_ThresholdByPSI_0.8' ,sep=',', header=T)
    df$Bin = cut(df$MutLoad, breaks=c(0,10,100,1000,50000), labels=c('0 - 10','10 - 100','100 - 1000','1000 - >10,000'))
    df = na.omit(df)
    df$FractionNMD= df$True/(df$False + df$True)
    foo = df %>%  group_by(Bin) %>% sample_frac(1, replace=TRUE) %>% 
                    summarize(retained = sum(True), not_retained = sum(False)) %>% 
                    group_by(Bin) %>% 
                    summarize(fraction = retained/(retained + not_retained))
    booted= BootAS(df)

    booted_CI = data.frame(booted %>% group_by(Bin) %>% summarise(mean = mean(fraction, na.rm = TRUE), sd = sd(fraction, na.rm = TRUE), n= n()) %>%
                mutate(se = sd / sqrt(n), lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se, upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se))
    PlotOut = ggplot(booted_CI, aes(x=as.character(Bin), y=mean, group=1)) + 
        geom_point() + theme_minimal() + geom_line() + #ylim(0.08,0.1) +
        geom_ribbon(aes(ymin=upper.ci, ymax=lower.ci), alpha=.7, position=position_dodge(.9)) +
        labs(y=expression(atop(underline("Under-expressed Transcripts With Intron Retention"), paste("All Transcripts With Intron Retention"))),
             x= 'Number of Protein Coding Mutations')
    ggsave(paste0(PlotDir, 'AS_RI_PerTumor_PSI_TCGA.pdf' ), width=4.5, height=4.5, units='in') 

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
    df= df[!duplicated(df$term_name),] # Remove duplicate hits
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
        geom_segment(x = 70, y = 10, xend = 1, yend = 10, colour = "grey", alpha=0.03, size=0.3 ) + # this set the scale line for y = 10
        geom_segment(x = 70, y = 5, xend = 1, yend = 5, colour = "grey", alpha=0.03, size=0.3 ) + # this set the scale line for y = 5
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
            hjust = hjust, angle = angle), size = 5)+
          annotate("text", x = rep(nrow(df),2), y = c(5,10), label = c("5","10") , color="grey", size=5 , alpha=1.4, angle=0, fontface="bold", hjust=0.25, vjust=2) + # Add annotation
        coord_polar() 
        print(head(df))
    ggsave(paste0(PlotDir, 'CircularRegressionGeneSetEnrichmentGlobalCORUMTCGA.pdf' ), width=12, height=12, units='in') 

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
        geom_segment(x = 70, y = 20, xend = 1, yend = 20, colour = "grey", alpha=0.03, size=0.3 ) + # this set the scale line for y = 20
        geom_segment(x = 70, y = 15, xend = 1, yend = 15, colour = "grey", alpha=0.03, size=0.3 ) + # this set the scale line for y = 15
        geom_segment(x = 70, y = 10, xend = 1, yend = 10, colour = "grey", alpha=0.03, size=0.3 ) + # this set the scale line for y = 10
        geom_segment(x = 70, y = 5, xend = 1, yend = 5, colour = "grey", alpha=0.03, size=0.3 ) + # this set the scale line for y = 5
        geom_bar(stat='identity') +
        labs(x='', y='Negative Log10 of Adjusted P-Value') + 
        #scale_fill_brewer(palette="Paired") +
        scale_fill_identity(guide = "legend", labels = unique(df$grouping),breaks = unique(df$colors))  +
        theme_minimal() +
        theme(axis.text = element_blank(),axis.title = element_blank(), panel.grid = element_blank(),
            legend.text=element_text(size=14),
            legend.position='bottom', plot.margin = unit(rep(-1,4), "cm"), legend.title = element_blank()) +  
        guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
        geom_text(aes(label = term_name, y = ifelse(negative_log10_of_adjusted_p_value< 6,5, negative_log10_of_adjusted_p_value+ GapForLegend), 
            hjust = hjust, angle = angle), size = 5)+ 
            annotate("text", x = rep(nrow(df),4), y = c(5,10,15,20), label = c("5","10","15","20") , 
            color="grey", size=5 , alpha=1, angle=0, fontface="bold", hjust=0.25, vjust=2) +
        coord_polar() 
    ggsave(paste0(PlotDir, 'CircularRegressionGeneSetEnrichmentGlobalKEGGTCGA.pdf' ), width=12, height=14, units='in') 

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
    Estimate = ggplot(data=df, aes(y = SortedLevel, x=as.numeric(as.character(Estimate)), color=subgroup)) +  geom_boxplot() +
                theme_minimal() + labs(x='Effect Size (Beta Coefficient)', y='') + theme(legend.position='bottom') +
                facet_wrap(~Dataset2, scales='free_x') +
                scale_color_manual(values=c('#8856a7','#F7B530','#3EA612')) + 
                #geom_vline(data=filter(quantile, Group == '0.05'), aes(xintercept=Estimate))+
                #geom_vline(xintercept=0, linetype='dashed', col = 'black') + 
                geom_vline(data = quantile, aes(xintercept = as.numeric(as.character(Estimate))), alpha=0.3, color='grey', size=1) +
                geom_text(data = quantile, aes(xintercept = as.numeric(as.character(Estimate)), 
                    label=Label, y=12), colour="grey",  vjust="inward", size=3) +
                theme(strip.text = element_text(face="bold", size=12),  panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
                theme(panel.spacing = unit(2, "lines"), legend.title=element_blank())
                
                #theme(strip.text = element_blank())
    # Rank = ggplot(data=df, aes(y = SortedLevel, x=as.numeric(as.character(NegLog10Pval)), color=subgroup)) +  geom_boxplot() +
    #             theme_minimal() + labs(x='Negative Log10 of Adjusted P-Value', y='') + 
    #             theme(legend.position='bottom', legend.title=element_blank()) +
    #             #theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    #             geom_vline(xintercept= -log10(0.05), alpha=0.5, size=1, col = 'grey') + facet_wrap(~Dataset2, scales='free_x') +
    #             theme(strip.text = element_text(face="bold", size=12))
    # Combined = plot_grid(Estimate, Rank, rel_heights=c(0.85, 1), ncol = 1) 
    #plot_grid(Title, Combined, ncol=1, rel_heights=c(0.1, 1)) 
    ggsave(paste0(PlotDir, 'RegCoefPerGroups_TCGAandCCCLE.pdf' ), width=8, height=4, units='in')


}


PlotDeltaPSICircular = function(df) {

    df$negative_log10_of_adjusted_p_value = -log10(df$p_value)
    df$Group = gsub('PosPSI','More IR\n(High vs. Low)', df$Group)
    df$Group = gsub('NegPSI','Less IR\n(High vs. Low)', df$Group)
    df$term_name = gsub('subunit, cytoplasmic','subunit,\ncytoplasmic',df$term_name)
    df$term_name = gsub('Capped ','Capped\n',df$term_name)
    df$term_name = gsub('Containing ','Containing\n',df$term_name)
    df$term_name = gsub(' Termination','\nTermination',df$term_name)
    df$term_name = gsub('Splicing -','Splicing -\n',df$term_name)
    to_add = setNames(data.frame( matrix(NA, 4*2, ncol(df))), colnames(df)) # Add empty space between groups
    to_add$Group = rep(unique(df$Group), each=4)
    df = rbind(df, to_add) %>% arrange(Group)
    df$term_name = gsub("\\(.*","",df$term_name) # remove parenthesis for visualization and repeat categories
    df= df[!duplicated(df$term_name),] # Remove duplicate hits
    GapForLegend=1

    # Add angles for lables in bar plot and sort by value
    df = df %>% 
        arrange(Group, source, -negative_log10_of_adjusted_p_value) %>%
        #arrange(factor(source, levels = unique(source))) %>%
        mutate(individual = factor(term_name, levels = unique(term_name))) %>%
        add_row(term_name = rep(NA,GapForLegend), negative_log10_of_adjusted_p_value = rep(NA,GapForLegend)) %>%
        mutate(ID = 1:(nrow(df) + GapForLegend)) %>%
        mutate(angle = 90-360*(ID-0.25)/(nrow(df) + 1)) %>%
        mutate(hjust = ifelse(angle < -90, 1, 0)) %>%
        mutate(angle = ifelse(angle < -90, angle+180, angle)) %>%
        mutate(ID = factor(ID, levels = unique(ID)))
  
  
        # prepare a data frame for lines in each group
        base_data = na.omit(df) %>% 
        group_by(Group) %>% 
        summarize(start=min(as.numeric(ID)), end=max(as.numeric(ID)) + 0.55) %>% 
        rowwise() %>% 
        mutate(title=mean(c(start, end)))
        
  
    
    PlotOut = ggplot(df, aes(x=ID, y=negative_log10_of_adjusted_p_value, fill=source)) +
        geom_segment(x = 70, y = 8, xend = 1, yend = 8, colour = "grey", alpha=0.03, size=0.3 ) + # this set the scale line for y = 20
        geom_segment(x = 70, y = 6, xend = 1, yend = 6, colour = "grey", alpha=0.03, size=0.3 ) + # this set the scale line for y = 15
        geom_segment(x = 70, y = 4, xend = 1, yend = 4, colour = "grey", alpha=0.03, size=0.3 ) + # this set the scale line for y = 10
        geom_segment(x = 70, y = 2, xend = 1, yend = 2, colour = "grey", alpha=0.03, size=0.3 ) + # this set the scale line for y = 5
        geom_bar(stat='identity') + 
        labs(x='', y='Negative Log10 of Adjusted P-Value') + 
        scale_fill_identity(guide = "legend", labels = unique(df$grouping),breaks = unique(df$colors))  +
        theme_minimal() + ylim(-4,10) + 
        theme(axis.text = element_blank(),axis.title = element_blank(), panel.grid = element_blank(),
            legend.text=element_text(size=14),
            legend.position='bottom', plot.margin = unit(rep(-1,4), "cm"), legend.title = element_blank()) +  
        guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
        geom_text(aes(label = term_name, y = ifelse(negative_log10_of_adjusted_p_value< 1,6, negative_log10_of_adjusted_p_value + GapForLegend), 
            hjust = hjust, angle = angle), size = 5)+ 
            annotate("text", x = rep(nrow(df),4), y = c(2,4,6,8), label = c("5","10","15","20") , 
            color="grey", size=8 , alpha=1, angle=0, fontface="bold", hjust=0.5, vjust=2) +
        coord_polar() + scale_fill_manual("legend", values = c('#F59F1A','#3279D8')) +
          # Add base line information
        geom_segment(data=base_data, aes(x = start, y = -0.5, xend = end, yend = -0.5), colour = "black", alpha=0.8, 
            size=0.6 , inherit.aes = FALSE )  + 
        geom_text(data=base_data, aes(x = title, y = -2, label=Group), hjust=c(0.6,0.5), vjust=c(0.3,0.5), colour = "black", 
            alpha=0.8, angle=60,size=5, inherit.aes = FALSE)
 
    
    ggsave(paste0(PlotDir, 'Circular_AS_Delta_PSI_IntonRetention_TCGA.pdf' ), width=12, height=20, units='in')
}


PlotDeltaPSI = function(df) {
    print(df)
    df$negative_log10_of_adjusted_p_value = -log10(df$p_value)
    df$Group = gsub('PosPSI','More Intron Retention\n(In High vs. Low)', df$Group)
    df$Group = gsub('NegPSI','Less Intron Retention\n(In High vs. Low)', df$Group)
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
    ggsave(paste0(PlotDir, 'AS_Delta_PSI_IntonRetention_TCGA.pdf' ), width=6, height=4.5, units='in')
}

PlotProtein = function(df) {
    df$Value = as.numeric(as.character(df$Value))
    df = df[df$Group %in% c('Cytoplasmic Ribosomes','Mitochondrial Ribosomes'),]
    ggviolin(df, x = "MutGroup", y = "Value", xlab = "", ylab = "Normalized Protein Expression",
        add = "boxplot", add.params = list(fill = "white"),
        fill = "Group", palette = "jco", facet.by = "Group", ylim = c(-2.5,2.5)) + 
        stat_compare_means(method = "t.test", label =  "p.signif", label.x = 1.5)
    # ggplot(data=df, aes(x = MutGroup, y=Value, color=Group)) +  geom_boxplot() + 
    # theme_minimal() + labs(y='Normalized Protein Expression',x='') + facet_wrsap(~Group)
    ggsave(paste0(PlotDir, 'Protein_Expression_CCLE.pdf' ), width=5, height=5, units='in')

}



PlotJacknifedExpressionAcrossGroups = function(df, dataset) {
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
    ggplot(df, aes(y=CancerTypeRemoved, x = name ,fill=Estimate)) + geom_tile()+
        #scale_fill_viridis(discrete=FALSE) +
        #scale_fill_gradientn(colours=c("blue","black","red")) +
        scale_fill_gradientn(colors=c("red","white","grey"), 
            values=rescale(c(min(df$Estimate),0, max(df$Estimate))),limits=c(min(df$Estimate),max(df$Estimate)))+
        geom_text(aes(label = Sig, color=Sig), show.legend = FALSE, size= 8) +
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


PlotJacknifedshRNAAcrossGroups = function(df) {
    df$AdjPval = p.adjust(df$Pr...t.., method= 'fdr')
    df$Sig = ifelse(df$Pr...t.. < 0.05, '*','')
    df$Sig = ifelse(df$Pr...t.. < 0.005, '**', df$Sig)
    df$Sig = ifelse(df$Pr...t.. < 0.0005, '***', df$Sig)
    df$Dummy = '' #'Protein Coding\nMutations'
    df$Group = gsub('mic Rib', 'mic\nRib', df$Group); df$Group = gsub('ial Rib', 'ial\nRib', df$Group); df$Group = gsub('ial Ch', 'ial\nCh', df$Group)
    df$Group = gsub('ory Pa', 'ory\nPa', df$Group)
    ggplot(df, aes(y=CancerTypeRemoved, x = Group ,fill=Estimate)) + geom_tile()+
        #scale_fill_viridis(discrete=FALSE) +
        #scale_fill_gradientn(colours=c("blue","black","red")) +
        scale_fill_gradientn(colors=c("red","grey","blue"), 
            values=rescale(c(min(df$Estimate),0, max(df$Estimate))),limits=c(min(df$Estimate),max(df$Estimate)))+
        geom_text(aes(label = Sig, color=Sig), show.legend = FALSE, size= 8) +
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


GetAnnotatedMOA = function() {
    AllDrugs = read.table('/labs/ccurtis2/tilk/scripts/protein/Data/Regression/AllDrugsOLSRegressionEstimatesKsKaCCLE', sep=',', header=T, fill=T)
    AllDrugs = subset(AllDrugs, AllDrugs$Coefficient == 'LogScore')
    AllDrugs = subset(AllDrugs, AllDrugs$Estimate < 0)
    # Replace specific drugs with moa categories into more specific groups
    AllDrugs$subgroup[which(AllDrugs$name == "aloe-emodin")] = 'Pro-Apoptosis'
    AllDrugs$subgroup[which(AllDrugs$name == "lenalidomide")] = 'ubiquitin specific protease inhibitor'
    AllDrugs$subgroup[which(AllDrugs$name == "10-deacetylbaccatin")] = 'microtubule inhibitor'
    AllDrugs$subgroup[which(AllDrugs$name == "mitoflaxone")] = 'angiogenesis inhibitor'
    AllDrugs$subgroup[which(AllDrugs$name == "fenspiride")] = 'Inflammatory/Immune'
    AllDrugs$subgroup[which(AllDrugs$name == "ibrutinib")] = 'Pro-Apoptosis'
    AllDrugs$subgroup[which(AllDrugs$name == "CNX-774")] = 'Pro-Apoptosis'  
    AllDrugs$subgroup[which(AllDrugs$name == "acalabrutinib")] = 'Pro-Apoptosis'
    AllDrugs$subgroup[which(AllDrugs$name == "ornithine")] = 'Nitrogen Metabolism Inhibitor'
    AllDrugs$subgroup[which(AllDrugs$name == "heptaminol")] = 'Angiogenesis Activator'
    PlottingGroup = rbind(
        data.frame(PlottingGroup = 'Growth Factor Inhibitors', NewGroup = c('dipeptidyl peptidase inhibitor', 
                                        'EGFR inhibitor',"FGFR inhibitor, VEGFR inhibitor","PDGFR tyrosine kinase receptor inhibitor",'CDK inhibitor, growth factor receptor inhibitor',
                                        "retinoid receptor agonist", 'VEGFR inhibitor','Abl kinase inhibitor, src inhibitor, VEGFR inhibitor','EGFR inhibitor, VEGFR inhibitor',
                                        'EGFR inhibitor, protein tyrosine kinase inhibitor','FLT3 inhibitor, growth factor receptor inhibitor, JAK inhibitor',
                                        'RAF inhibitor, VEGFR inhibitor','tropomyosin receptor kinase inhibitor','PDGFR tyrosine kinase receptor inhibitor, VEGFR inhibitor',
                                        'hepatocyte growth factor receptor inhibitor','insulin growth factor receptor inhibitor','MER tyrosine kinase inhibitor',
                                        'EGFR inhibitor, FGFR inhibitor, FLT3 inhibitor, PDGFR tyrosine kinase receptor inhibitor, VEGFR inhibitor','GK0582 inhibitor',
                                        'FGFR inhibitor, KIT inhibitor, PDGFR tyrosine kinase receptor inhibitor, RAF inhibitor, RET tyrosine kinase inhibitor, VEGFR inhibitor',
                                        'Aurora kinase inhibitor, growth factor receptor inhibitor','bone morphogenic protein inhibitor','growth factor receptor inhibitor')),
        data.frame(PlottingGroup = 'Transport Inhibitor', NewGroup = c('','MRP inhibitor, P glycoprotein inhibitor')),
        data.frame(PlottingGroup = 'Pro-Apoptosis' , NewGroup = c('Pro-Apoptosis',"AKT inhibitor","ALK tyrosine kinase receptor inhibitor","Aurora kinase inhibitor",'CDK inhibitor',
                                        'JNK inhibitor','mTOR inhibitor, PI3K inhibitor','NFkB pathway inhibitor','other antibiotic','phosphodiesterase inhibitor',
                                        'IKK inhibitor, NFkB pathway inhibitor','IKK inhibitor','K-ras inhibitor','KIT inhibitor, VEGFR inhibitor',
                                        'KIT inhibitor, PDGFR tyrosine kinase receptor inhibitor, VEGFR inhibitor','MAP kinase inhibitor','MAP kinase inhibitor, MEK inhibitor',
                                        'MAP kinase phosphatase inhibitor','MET inhibitor','survivin inhibitor','src inhibitor','src inhibitor, syk inhibitor',
                                        'PI3K inhibitor, mTOR inhibitor','mTOR inhibitor','monopolar spindle 1 kinase inhibitor','NAMPT inhibitor',
                                        'protein tyrosine kinase inhibitor, tyrosine kinase inhibitor','tumor necrosis factor production inhibitor','MRP inhibitor, P glycoprotein inhibitor',
                                        'SIRT inhibitor','WEE1 kinase inhibitor','WDR5/MLL interaction inhibitor','vitamin D receptor agonist','phosphoinositide dependent kinase inhibitor',
                                        'protein kinase inhibitor','PKC inhibitor, protein tyrosine kinase inhibitor','MDM inhibitor','serine arginine protein kinase inhibitor',
                                        'tumor necrosis factor receptor antagonist','tyrosine phosphatase inhibitor','Wnt pathway inhibitor','TP53 inhibitor',
                                        'maternal embryonic leucine zipper kinase inhibitor','MEK inhibitor','PARP inhibitor','PKC activator','Ras GTPase inhibitor',
                                        'PI3K inhibitor','PKC inhibitor','PLK inhibitor','protein tyrosine kinase inhibitor','RAF inhibitor','Bcr-Abl kinase inhibitor',
                                        'tyrosine kinase inhibitor','XIAP inhibitor',"ALK tyrosine kinase receptor inhibitor, insulin growth factor receptor inhibitor",
                                        'AKT inhibitor, mTOR inhibitor','AKT inhibitor, MAP kinase inhibitor','AKT inhibitor, pyruvate dehydrogenase inhibitor',
                                        'breast cancer resistance protein inhibitor','caspase activator','carnitine palmitoyltransferase inhibitor',
                                        'thymidine phosphorylase inhibitor','focal adhesion kinase inhibitor','exportin antagonist',
                                        'p21 activated kinase inhibitor','nitric oxide production inhibitor','gap junction modulator, nitric oxide production inhibitor',
                                        'smoothened receptor antagonist','stearoyl-CoA desaturase inhibitor','ornithine decarboxylase inhibitor','mediator release inhibitor, SYK inhibitor',
                                        'DNA protein kinase inhibitor, PI3K inhibitor','FLT3 inhibitor, JAK inhibitor','farnesyltransferase inhibitor',
                                        'FLT3 inhibitor, KIT inhibitor, PDGFR tyrosine kinase receptor inhibitor, RAD51 inhibitor, RET tyrosine kinase inhibitor',
                                        'CDK inhibitor, glycogen synthase kinase inhibitor','CDK inhibitor, FLT3 inhibitor, JAK inhibitor','cFMS kinase inhibitor',
                                        'ALK tyrosine kinase receptor inhibitor, insulin growth factor receptor inhibitor','HDAC inhibitor','ARF inhibitor','BMI-1 inhibitor',
                                        'Aurora kinase inhibitor, JAK inhibitor','BCL inhibitor','BCL inhibitor, MCL1 inhibitor','Bcl-XL downregulator','beta-catenin inhibitor',
                                        'apoptosis stimulant, ribonucleotide reductase inhibitor','apoptosis stimulant','ribonucleotide reductase inhibitor','c-Myc inhibitor',
                                        'DNA damage inducer','CDC inhibitor','ephrin inhibitor','hedgehog pathway inhibitor','HCV inhibitor, XIAP inhibitor',
                                        'telomerase inhibitor','rho associated kinase inhibitor','porcupine inhibitor','phosphoinositide dependent kinase inhibitor',
                                        'Bcr-Abl kinase inhibitor, LYN tyrosine kinase inhibitor','Bcr-Abl kinase inhibitor, FLT3 inhibitor, PDGFR tyrosine kinase receptor inhibitor',
                                        'adenosine deaminase inhibitor, ribonucleotide reductase inhibitor','apoptosis stimulant, NFkB pathway inhibitor','ATM kinase inhibitor',
                                        'casein kinase inhibitor, cell proliferation inhibitor','casein kinase inhibitor, mTOR inhibitor, PI3K inhibitor','casein kinase inhibitor')),
        data.frame(PlottingGroup = 'Anti-Apoptosis' , NewGroup = c('protein tyrosine kinase activator','nitric oxide stimulant','hepatocyte function enhancer',
                                        'src activator','lysophosphatidic acid receptor antagonist','protein kinase activator','p53 activator')),   
        data.frame(PlottingGroup = 'Energy Metabolism Inhibitor' , NewGroup = c("ATPase inhibitor",'ATP synthase inhibitor','AMPK inhibitor','ATP channel blocker',
                                        'ATP synthase inhibitor, ATPase inhibitor','ATP-sensitive potassium channel antagonist','ATPase inhibitor, gastrin inhibitor',
                                        'isocitrate dehydrogenase inhibitor','P glycoprotein inhibitor','pyruvate dehydrogenase kinase inhibitor','NTPDase inhibitor',
                                        'pyruvate dehydrogenase inhibitor',
                                        'ATPase inhibitor, TRPV agonist','ATP channel blocker, insulin secretagogue','autotaxin inhibitor','carbonic anhydrase inhibitor',
                                        'cardiac myosin activator')),
        data.frame(PlottingGroup = 'Inflammatory/Immune' , NewGroup = c('analgesic agent','CC chemokine receptor antagonist','cyclooxygenase inhibitor',
                                        'glucocorticoid receptor agonist','histamine receptor antagonist','histamine receptor agonist','nitric oxide synthase inhibitor','local anesthetic',
                                        'p38 MAPK inhibitor','SYK inhibitor','anti-inflammatory agent, glucocorticoid receptor agonist','antihistamine',
                                        'anti-inflammatory agent','antiinflammatory agent, cyclooxygenase inhibitor','antiinflammatory agent','NOD1 inhibitor',
                                        'macrophage migration inhibiting factor inhibitor','lipoxygenase inhibitor','TGF beta receptor inhibitor',
                                        'TGF beta receptor inhibitor, p38 MAPK inhibitor','MALT1 inhibitor (JH)','mediator release inhibitor','melatonin receptor antagonist',
                                        'metalloproteinase inhibitor','sphingosine 1 phosphate receptor agonist','melanocortin receptor agonist','neprilysin inhibitor',
                                        'matrix metalloprotease inhibitor, tumor necrosis factor production inhibitor',' ediator release inhibitor, SYK inhibitor',
                                        'cyclooxygenase inhibitor, glutathione peroxidase agonist, H+/K+-ATPase inhibitor, nitric oxide synthase inhibitor',
                                        'cyclooxygenase inhibitor, histone acetyltransferase inhibitor, lipoxygenase inhibitor, NFkB pathway inhibitor',
                                        'glycogen synthase kinase inhibitor, lipoxygenase inhibitor','histamine receptor modulator','integrin antagonist',
                                        'cyclooxygenase inhibitor, lipoxygenase inhibitor','corticosteroid agonist','cyclooxygenase inhibitor, prostanoid receptor antagonist',
                                        'cyclooxygenase inhibitor, prostaglandin inhibitor','cyclooxygenase inhibitor, platelet aggregation inhibitor',
                                        'neurotrophic agent','steroid','syk inhibitor',
                                        'STAT inhibitor','prostanoid receptor agonist','prostanoid receptor antagonist','tachykinin antagonist','antitussive',
                                        'cyclooxygenase inhibitor, NFkB pathway inhibitor','cyclooxygenase inhibitor, prostanoid receptor agonist','toll-like receptor antagonist',
                                        'inosine monophosphate dehydrogenase inhibitor','interferon inducer','integrin antagonist','interleukin receptor antagonist, STAT inhibitor',
                                        'interleukin synthesis inhibitor','interleukin inhibitor','JAK inhibitor','JAK inhibitor, STAT inhibitor','JAK inhibitor, syk inhibitor',
                                        'immunosuppressant, inosine monophosphate dehydrogenase inhibitor','immunosuppressant, protein synthesis inhibitor, purine antagonist',
                                        'glutamate inhibitor, macrophage migration inhibiting factor inhibitor','histamine receptor antagonist, platelet activating factor receptor antagonist',
                                        'antioxidant, capillary stabilizing agent, nitric oxide scavenger','antihistamine','calcineurin inhibitor','toll-like receptor agonist',
                                        'leukotriene receptor antagonist, phosphodiesterase inhibitor','leukotriene synthesis inhibitor','NOD1 inhibitor',
                                        'indoleamine 2,3-dioxygenase inhibitor','ICAM1 expression inhibitor','ICAM1 antagonist','leukotriene receptor antagonist',
                                        'CCR antagonist','cytokine production inhibitor','cytokine production inhibitor, NFkB pathway inhibitor, TP53 activator',
                                        'dehydrogenase inhibitor, inositol monophosphatase inhibitor','folate receptor ligand','glucocorticoid receptor agonist, immunosuppressant',
                                        'glucocorticoid receptor antagonist, progesterone receptor antagonist','immunostimulant','immunosuppressant',
                                        'aryl hydrocarbon receptor agonist, indoleamine 2,3-dioxygenase inhibitor')),
        data.frame(PlottingGroup = 'Protein Synthesis Inhibitor' , NewGroup = c('antimalarial agent','bacterial 50S ribosomal subunit inhibitor','protein synthesis inhibitor',
                                        'bacterial 30S ribosomal subunit inhibitor','bacterial 30S ribosomal subunit inhibitor, bacterial 50S ribosomal subunit inhibitor',
                                        'FOXM1 inhibitor, protein synthesis inhibitor','ribosomal protein inhibitor','cytochrome P450 inhibitor, protein synthesis inhibitor',
                                        'cytochrome P450 inhibitor, imidazoline receptor ligand',
                                        'cytochrome P450 inhibitor','eukaryotic translation elongation factor 2 inhibitor','eukaryotic translation initiation factor inhibitor')),
        data.frame(PlottingGroup = 'Protein Synthesis Activator' , NewGroup = c('PERK inhibitor')),
        data.frame(PlottingGroup = 'DNA Replication Inhibitor' , NewGroup = c('antiprotozoal agent','ATR kinase inhibitor','bacterial DNA gyrase inhibitor','CHK inhibitor',
                                        'DNA inhibitor','DNA polymerase inhibitor','topoisomerase inhibitor','bacterial DNA inhibitor','CDK inhibitor, cell cycle inhibitor, MCL1 inhibitor',
                                        'anticancer agent, aryl hydrocarbon receptor antagonist','antiviral','bacterial antifolate','CDK inhibitor, CHK inhibitor, PKC inhibitor',
                                        'purine antagonist',
                                        'DNA alkylating agent, DNA synthesis inhibitor','DNA inhibitor, topoisomerase inhibitor','DNA directed DNA polymerase inhibitor, thymidylate synthase inhibitor',
                                        'DNA alkylating agent, DNA inhibitor','dihydropteroate synthetase inhibitor','DNA directed DNA polymerase inhibitor','dehydrogenase inhibitor',
                                        'DNA alkylating agent','IGF-1 inhibitor','mitotic kinase inhibitor','PABA antagonist','Pim kinase inhibitor','thymidylate synthase inhibitor',
                                        'DNA replication inhibitor, STAT inhibitor','DNA synthesis inhibitor','DNA dependent protein kinase inhibitor','DNA synthesis inhibitor, RNA synthesis inhibitor',
                                        'cell cycle inhibitor, PLK inhibitor','chelating agent, topoisomerase inhibitor','chelating agent','cyclin D inhibitor','dihydrofolate reductase inhibitor',
                                        'dihydropteroate synthase inhibitor','dihydrofolate reductase inhibitor, thymidylate synthase inhibitor','dihydroorotate dehydrogenase inhibitor')),
        data.frame(PlottingGroup = 'Cytoskeleton Organization Inhibitor' , NewGroup = c('bacterial cell wall synthesis inhibitor','fungal lanosterol demethylase inhibitor',
                                        'gap junction modulator','SHIP2 phosphatase inhibitor','sterol demethylase inhibitor','sterol methyltransferase inhibitor','utrophin enhancer',
                                        'kinesin-like spindle protein inhibitor','membrane integrity inhibitor','microtubule inhibitor','beta lactamase inhibitor','DNA intercalating agent',
                                        'microtubule stabilizing agent, tubulin polymerization inhibitor','tubulin polymerization inhibitor','DNA synthesis inhibitor, microtubule inhibitor',
                                        'microtubule stabilizing agent','microtubule inhibitor, tubulin polymerization inhibitor','src inhibitor, tubulin polymerization inhibitor',
                                        'DPRE1 inhibitor','kinesin inhibitor, kinesin-like spindle protein inhibitor','kinesin inhibitor','tubulin polymerization inhibitor, VE-cadherin antagonist',
                                        'bradykinin receptor antagonist, tubulin polymerization inhibitor','penicillin binding protein inhibitor')),
        data.frame(PlottingGroup = 'Transcription Inhibitor' , NewGroup = c('bromodomain inhibitor','RNA polymerase inhibitor','RNA synthesis inhibitor','RNA synthesis inhibitor, topoisomerase inhibitor')),
        data.frame(PlottingGroup = 'Ion Channel Regulation', NewGroup =  c('calcium channel blocker','potassium channel blocker','potassium channel activator','TRPV agonist','TRPV antagonist',
                                        'sodium channel blocker','antiarrhythmic medication','antiarrhythmic','anticestodal agent','calcium sensitizer','calmodulin antagonist',
                                        'calcium receptor antagonist','calcium channel modulator','calcium channel activator','calcium sensitizer, phosphodiesterase inhibitor',
                                        'cardiac glycoside','cationic surfactant','chloride channel blocker','chloride channel antagonist','electrolyte reabsorption inhibitor',
                                        'voltage-gated sodium channel modulator','voltage-gated sodium channel blocker','voltage-gated calcium channel ligand','bacterial permeability inducer',
                                        'succinimide antiepileptic','slow afterhyperpolarization channel blocker','TRPA1 channel blocker','CFTR channel agonist',
                                        'T-type calcium channel binder','T-type calcium channel blocker','sodium/glucose cotransporter inhibitor','sodium/potassium/chloride transporter inhibitor',
                                        'sodium channel blocker, T-type calcium channel blocker',' slow afterhyperpolarization channel blocker','phosphate antagonist',
                                        'potassium channel blocker, sodium channel blocker','N-type calcium channel blocker','Na/K-ATPase inhibitor','mineralocorticoid receptor agonist',
                                        'mineralocorticoid receptor antagonist','potassium-competitive acid antagonist','serum/glucocorticoid regulated kinase inhibitor',
                                        'phosphodiesterase inhibitor, sodium channel blocker','nitric oxide donor, potassium channel activator','nitric oxide donor',
                                        'electrolyte reabsorption inhibitor, thromboxane receptor antagonist','electrolyte reabsorption inhibitor, thromboxane receptor antagonist',
                                        'collapsin response mediator protein inhibitor, collapsin response mediator protein stimulant, sodium channel blocker','Kir6 channel (KATP) activator',
                                        'HCN channel antagonist, potassium channel blocker, sodium channel blocker','HCN channel blocker','hydantoin antiepileptic','thiazide diuretic',
                                        'L-type calcium channel activator','L-type calcium channel blocker','KATP activator, Kir6 channel (KATP) activator, vasodilator','vasoconstrictor',
                                        'glutamate inhibitor','glutamate receptor antagonist, calcium channel blocker','glutamate receptor modulator','ionotropic receptor IR40a activator',
                                        'glutamate receptor negative allosteric modulator',' glutamate receptor positive allosteric modulator','intermediate conductance potassium channel blocker',
                                        'chloride channel agonist, benzodiazepine receptor agonist','chloride reabsorption inhibitor','electron acceptor','membrane permeability inhibitor',
                                        'chloride channel activator','chloride channel blocker, gap junction modulator, glutamate inhibitor')),
        data.frame(PlottingGroup = 'Protein Degradation Inhibitors', NewGroup =  c("HCV inhibitor",'HIV protease inhibitor',"proteasome inhibitor",'ubiquitin specific protease inhibitor',
                                        'NFkB pathway inhibitor, proteasome inhibitor','ubiquitin C-terminal hydrolase inhibitor','ubiquitin-conjugating enzyme inhibitor',
                                        'deubiquitinase inhibitor','cysteine protease inhibitor, calpain inhibitor','DUB inhibitor','nedd activating enzyme inhibitor')),
        data.frame(PlottingGroup = 'Chaperone Inhibitors', NewGroup =  c("HSP inhibitor",'HSP antagonist')),
        data.frame(PlottingGroup = 'Chaperone Activators', NewGroup =  c("HSP inducer")),
        data.frame(PlottingGroup = 'Oxidative Stress Activator', NewGroup =  c("gamma glutamyltransferase Inhibitors",'glucose dependent insulinotropic receptor agonist',
                                        'mitochondrial electron transport inhibitor','oxidative stress inducer','MTH1 inhibitor','superoxide dismutase inhibitor')),
        data.frame(PlottingGroup = 'Lipid Metabolism Inhibitor', NewGroup = c('3-ketoacyl CoA thiolase inhibitor','ACAT inhibitor','atherogenesis inhibitor','cholesteryl ester transfer protein inhibitor',
                                        'ACAT inhibitor, sterol regulatory element binding protein (SREBP) inhibitor','adipose triglyceride lipase inhibitor','cholesterol inhibitor',
                                        'lipase inhibitor','lipid peroxidase inhibitor','microsomal trigylceride transfer protein inhibitor','phosphatidylglycerophosphatase inhibitor',
                                        'fatty acid synthase inhibitor','11-beta hydroxysteroid dehydrogenase inhibitor, FXR agonist','FXR agonist','phospholipase inhibitor',
                                        'diacylglycerol kinase inhibitor, protein kinase inhibitor','PPAR receptor antagonist',
                                        'PKA inhibitor','PPAR receptor agonist','secretory phospholipase inhibitor','oxidosqualene cyclase inhibitor','omega 3 fatty acid stimulant',
                                        'cholesterol inhibitor, Niemann-Pick C1-like 1 protein antagonist','diacylglycerol O acyltransferase inhibitor')),
        data.frame(PlottingGroup = 'Lipid Metabolism Activator', NewGroup = c('apolipoprotein expression enhancer','constitutive androstane receptor (CAR) agonist','FABI inhibitor',
                                        'lipoprotein lipase activator','LDL receptor activator','monoacylglucerol lipase inhibitor','mycolic synthesis inhibitor',
                                        'phospholipase activator',
                                        'free fatty acid receptor agonist','fungal squalene epoxidase inhibitor','gamma secretase inhibitor','gamma secretase modulator','coenzyme A precursor')),
        data.frame(PlottingGroup = 'Growth Hormone Inhibitor', NewGroup = c('5 alpha reductase inhibitor','androgen receptor antagonist', "progesterone receptor agonist",'adenylyl cyclase inhibitor',
                                        'androgen biosynthesis inhibitor','androgen receptor agonist','androgen receptor agonist, estrogen receptor agonist',
                                        'aromatase inhibitor, TRPV antagonist','estrogen receptor agonist','estrogen receptor agonist, estrogenic hormone',
                                        'somatostatin receptor agonist',
                                        'thyrotropin releasing hormone receptor agonist','thyroid hormone inhibitor','selective estrogen receptor destabilizer',
                                        'progesterone receptor antagonist','prostaglandin inhibitor','prostaglandin receptor agonist','prostaglandin synthesis inhibitor',
                                        'estrogen receptor antagonist, selective estrogen receptor modulator (SERM)','estrogen receptor antagonist','estrogen-related receptor agonist',
                                        'estrogenic hormone','estrogen receptor antagonist, progesterone receptor agonist','GLP receptor agonist','progestogen hormone',
                                        'gonadotropin releasing factor hormone receptor antagonist','steroid sulfatase inhibitor','selective estrogen receptor modulator (SERM)',
                                        'growth hormone secretagogue receptor agonist','gonadotropin releasing factor hormone receptor agonist','gonadotropin inhibitor',
                                        'estrogen receptor agonist, glucocorticoid receptor antagonist, progesterone receptor agonist, progesterone receptor antagonist',
                                        'androgen receptor agonist, estrogen receptor agonist, progesterone receptor agonist','aromatase inhibitor','calcitonin antagonist')),
        data.frame(PlottingGroup = 'Growth Hormone Activator', NewGroup = c('androgen receptor enhancer',' androgen receptor modulator','thyroid hormone stimulant')),
        data.frame(PlottingGroup = 'Angiogenesis Inhibitor', NewGroup = c('angiogenesis inhibitor','angiogenesis inhibitor, tumor necrosis factor production inhibitor',
                                        'endothelin receptor antagonist','hypoxia inducible factor inhibitor','matrix metalloprotease inhibitor','P selectin inhibitor',
                                        'platelet aggregation inhibitor','platelet activating factor receptor antagonist','thromboxane receptor antagonist',
                                        'vasopressin receptor agonist','vasopressin receptor antagonist',
                                        'thrombopoietin receptor agonist','vitamin K antagonist','coagulation factor inhibitor','plasminogen activator inhibitor',
                                        'thromboxane receptor antagonist, thromboxane synthase inhibitor','thromboxane receptor antagonist','thrombin inhibitor',
                                        'angiotensin antagonist','angiotensin converting enzyme inhibitor','angiotensin receptor agonist','angiotensin receptor antagonist',
                                        'bradykinin receptor antagonist','phosphodiesterase inhibitor, platelet aggregation inhibitor','thromboxane synthase inhibitor')),
        data.frame(PlottingGroup = 'Angiogenesis Activator', NewGroup = c('renin inhibitor','haemostatic agent')),                         
        data.frame(PlottingGroup = 'Neurotransmitter Inhibitor', NewGroup = c('acetylcholine receptor agonist','acetylcholine receptor agonist, benzodiazepine receptor agonist',
                                        'acetylcholine receptor allosteric modulator', 'dopamine receptor agonist','dopamine receptor antagonist','cholinergic receptor antagonist',
                                        'cholinesterase inhibitor','cholinesterase reactivator','cholinergic receptor agonist','dopamine receptor agonist, serotonin receptor antagonist',
                                        'dopamine receptor ligand','dopamine beta hydroxylase inhibitor','dopamine uptake inhibitor','dopamine reuptake inhibitor, serotonin reuptake inhibitor',
                                        'dopamine reuptake inhibitor','dopamine receptor partial agonist','dopamine reuptake inhibitor, selective serotonin reuptake inhibitor (SSRI)',
                                        'sigma receptor agonist','sigma receptor agonist, sigma receptor antagonist','sigma receptor antagonist','sigma receptor ligand',
                                        'nicotinic receptor agonist','norepinephrine reuptake inhibitor, serotonin-norepinephrine reuptake inhibitor (SNRI)',
                                        'norepinephrine reuptake inhibitor, serotonin-norepinephrine reuptake inhibitor (SNRI), tricyclic antidepressant',
                                        'neuromuscular blocker','opioid receptor agonist','opioid receptor agonist, opioid receptor antagonist','opioid receptor antagonist',
                                        'opioid receptor modulator','purinergic receptor antagonist','tricyclic antidepressant','sedative',
                                        'norepinephrine transporter inhibitor','norepinephrine reuptake inhibitor','trace amine associated receptor agonist',
                                        'norepinephrine reputake inhibitor, opioid receptor agonist, serotonin reuptake inhibitor','norepinephrine reputake inhibitor, tricyclic antidepressant',
                                        'norepinephrine reputake inhibitor','noradrenaline uptake inhibitor','nootropic agent','tyrosine hydroxylase inhibitor',
                                        'vesicular monoamine transporter inhibitor','neuropeptide receptor antagonist','tryptophan hydroxylase inhibitor','antispasmodic',
                                        'norepinephrine inhibitor, norepinephrine reuptake inhibitor, serotonin receptor antagonist, serotonin-norepinephrine reuptake inhibitor (SNRI)',
                                        'melatonin receptor agonist','serotonin transporter (SERT) inhibitor','serotonin reuptake inhibitor','serotonin receptor agonist, serotonin receptor antagonist',
                                        'dopamine-norepinephrine reuptake inhibitor','excitatory amino acid transporter inhibitor','GABA receptor antagonist, GABA receptor modulator',
                                        'GABA benzodiazepine site receptor inverse agonist','GABA benzodiazepine site receptor agonist','GABA aminotransferase inhibitor','GABA uptake inhibitor',
                                        'GABA uptake inhibitor, GAT inhibitor','GABA receptor modulator','GABA receptor antagonist','GABA receptor negative allosteric modulator',
                                        'GABA receptor antagonist, serotonin receptor antagonist','GABA transaminase inhibitor','gamma hydroxybutyric acid ligand','imidazoline receptor antagonist',
                                        'Glycine transporter 1 inhibitor','glycine transporter inhibitor, GlyT-1 inhibitor','GAT inhibitor','imidazoline receptor agonist',
                                        'neurotransmitter agonist','orexin receptor antagonist','oxytocin receptor antagonist','serotonin-norepinephrine reuptake inhibitor (SNRI)',
                                        'kainate receptor agonist','serotonin receptor agonist, dopamine receptor agonist','selective serotonin reuptake inhibitor (SSRI)',
                                        'dopamine receptor antagonist, serotonin receptor antagonist',"acetylcholine receptor antagonist","adenosine receptor antagonist","adenosine receptor agonist",
                                        "adrenergic receptor agonist","adrenergic receptor antagonist",'benzodiazepine receptor agonist',"glutamate receptor antagonist",'cannabinoid receptor antagonist',
                                        'serotonin receptor agonist','serotonin receptor antagonist','acetylcholinesterase inhibitor','acetylcholinesterase inhibitor, beta amyloid synthesis inhibitor',
                                        'acetylcholinesterase inhibitor, microtubule inhibitor','acetylcholinesterase inhibitor, monoamine oxidase inhibitor',' barbiturate antiepileptic, GABA receptor modulator',
                                        "adenosine receptor antagonist, hemoglobin antagonist",'adenosine receptor antagonist, phosphodiesterase inhibitor','anxiolytic','barbiturate antiepileptic, GABA receptor modulator',
                                        'adenosine reuptake inhibitor, phosphodiesterase inhibitor','adrenergic receptor agonist, dopamine receptor agonist','anticonvulsant','cannabinoid receptor agonist',
                                        'adrenergic receptor agonist, serotonin receptor agonist','benzodiazepine receptor antagonist','benzodiazepine receptor ligand','barbiturate antiepileptic, GABA receptor modulator',
                                        'adrenergic receptor antagonist, glutamate receptor antagonist','adrenergic receptor ligand','adrenergic inhibitor','anthelmintic agent',' cannabinoid receptor inverse agonist',
                                        'cannabinoid receptor inverse agonist','cannabinoid receptor agonist, glucose dependent insulinotropic receptor agonist, potassium channel blocker, PPAR receptor agonist',
                                        'beta amyloid protein neurotoxicity inhibitor','FAAH inhibitor','FAAH inhibitor, FAAH reuptake inhibitor','glutamate receptor agonist','imidazoline receptor ligand',
                                        'cannabinoid receptor modulator','catechol O methyltransferase inhibitor')),
        data.frame(PlottingGroup = 'Neurotransmitter Enhancer', NewGroup = c('acetylcholine release enhancer','androgen receptor modulator','barbiturate antiepileptic, GABA receptor modulator',
                                        'enkephalinase inhibitor','glutamate receptor positive allosteric modulator','monoamine oxidase inhibitor')),
        data.frame(PlottingGroup = 'Oxidative Stress Inhibitors', NewGroup = c('antioxidant', "nuclear factor erythroid derived, like (NRF2) activator",'11-beta hydroxysteroid dehydrogenase inhibitor',
                                        'free radical scavenger','flavanone glycoside','glutathione transferase stimulant','NADPH oxidase inhibitor','NADPH inhibitor','melanin inhibitor',
                                        'RAGE receptor antagonist','NAD precursor, vitamin B','xanthine oxidase inhibitor')),                          
        data.frame(PlottingGroup = 'Other', NewGroup = c('antacid','antiinfective drug','levodropropizine','bone resorption inhibitor','contrast agent','expectorant','laxative','topical sunscreen agent',
                                        'topical anesthetic','radiopaque medium','psychoactive drug','antiseptic','mucolytic agent','mucus protecting agent','antiseptic','choleretic agent','pharmacological chaperone',
                                        'panipenem uptake inhibitor','pregnane X receptor agonist','tyrosinase inhibitor','protein synthesis stimulant',
                                        'solute carrier family member inhibitor','transthyretin amyloid inhibitor','anti-pneumocystis agent',
                                        'gene expression stimulant','hypercalcaemic agent','motilin receptor agonist','muscle relaxant','myorelaxant','nematocide','osmosis stimulant', 'other antifungal',
                                        'cosmetic','cosmetic moisturizer','coloring agent','caspase inhibitor','contraceptive agent','excipient','enzyme inducer','growth factor receptor activator')),
        data.frame(PlottingGroup = 'Histone/Methylation Inhibitor', NewGroup = c('adenosylhomocysteinase inhibitor','DNA methyltransferase inhibitor','histone lysine methyltransferase inhibitor',
                                        'L3MBTL antagonist')),
        data.frame(PlottingGroup = 'Alcohol Metabolism Inhibitor', NewGroup = c('alcohol dehydrogenase inhibitor')),
        data.frame(PlottingGroup = 'Sugar Metabolism Inhibitor', NewGroup = c('aldose reductase inhibitor','alpha mannosidase inhibitor','fungal 1,3-beta-D-glucan synthase inhibitor',
                                        'glucokinase inhibitor','glucosidase inhibitor','glycolysis inhibitor','insulin sensitizer, PPAR receptor agonist','insulin sensitizer',
                                        'thiamine uptake blocker','neprilysin','hypercalcaemic','sulfonylurea',
                                        'insulin secretagogue','phosphofructokinase inhibitor','niacinamide phosphoribosyltransferase inhibitor','neuraminidase inhibitor')),
        data.frame(PlottingGroup = 'Sugar Metabolism Activator', NewGroup = c('glucokinase activator','glucagon receptor antagonist','glycogen synthase kinase inhibitor',
                                        'phosphoenolpyruvate carboxylase activator, serine/threonine protein phosphatase activator')),
        
        data.frame(PlottingGroup = 'G-Coupled Receptor Regulator', NewGroup = c('aryl hydrocarbon receptor antagonist','aryl hydrocarbon receptor agonist, indoleamine 2,3-dioxygenase inhibitor',
                                        'CCK receptor agonist','CCK receptor antagonist','G protein-coupled receptor agonist','guanylate cyclase stimulant','guanylate cyclase activator',
                                        'IP1 prostacyclin receptor agonist','neurokinin receptor antagonist','protease-activated receptor inhibitor','nociceptin/orphanin FQ receptor antagonist')),
        data.frame(PlottingGroup = 'Protein Metabolism Inhibitor', NewGroup =  c("aspartic protease inhibitor",'beta-secretase inhibitor','cathepsin inhibitor','collagenase inhibitor',
                                        'collagenase inhibitor, metalloproteinase inhibitor','collagenase inhibitor','glycosylation inhibitor','Glycosyl transferase inhibitor',
                                        'Procollagen C-Endopeptidase Inhibitors','vitamin B','protease inhibitor','serine protease inhibitor',
                                        'hydroxyphenylpyruvate dioxygenase inhibitor','peptidase inhibitor','methylmalonyl CoA mutase stimulant, vitamin B','mannosidase inhibitor',
                                        'HCV inhibitor, protease inhibitor','elastase inhibitor','HMGCR inhibitor')),
        data.frame(PlottingGroup = 'Nitrogen Metabolism Inhibitor', NewGroup =  c('carbamoyl phosphate synthase activator','urate transporter inhibitor',
                                        'urease inhibitor','uricase inhibitor','uricosuric blocker')),
        data.frame(PlottingGroup = 'Viral Replication Inhibitor', NewGroup = c('nucleoside reverse transcriptase inhibitor','non-nucleoside reverse transcriptase inhibitor',
                                        'HIV attachment inhibitor','HIV integrase inhibitor','RSV fusion inhibitor'))
    )
    return(merge(AllDrugs, PlottingGroup, by.x='subgroup', by.y='NewGroup', all.x=T))
    

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


### Heat - map version !
# PlotModelDiagnosticsCCLEDrug = function(df) {
#     # Fitted values vs residuals for full model and all explanatory variables
#     df = subset(df, (df$subgroup != 'RNA Synthesis Inhibitor') &  (df$subgroup != 'Protein Synthesis Inhibitor'))
#     df$Diagnostic = gsub('MutLoad','Mutational Load',df$Diagnostic);df$Diagnostic = gsub('FullModel','Full Model',df$Diagnostic);
#     df$subgroup = gsub(' ','\n', df$subgroup)
#     print(unique(df$subgroup))
#     df$PearsonsR = as.numeric(as.character(df$PearsonsR))
#     Plot = ggplot(df, aes(y=name, x = Diagnostic ,fill=PearsonsR)) + geom_tile()+
#         #scale_fill_gradientn(colors=c("red","grey","blue"), 
#         #    values=rescale(c(min(df$PearsonsR),0, max(df$PearsonsR))),limits=c(min(df$PearsonsR),max(df$PearsonsR)))+
#         theme_minimal() +
#         facet_grid(rows=vars(subgroup), scales='free_y', space='free_y') +
#         labs(y='',x='', fill="Pearson's R") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
#         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(text = element_text(size=18))+
#         ggtitle('PRISM (Drug Screen)')
#     ggsave(paste0(PlotDir, 'ModelDiagnosticCCLEDrug.pdf'), width=8, height=12, units='in')
# }


PlotModelDiagnosticsCCLEDrug = function(df) {
    # Fitted values vs residuals for full model and all explanatory variables
    df = subset(df, (df$subgroup != 'RNA Synthesis Inhibitor') &  (df$subgroup != 'Protein Synthesis Inhibitor'))
    df$Diagnostic = gsub('MutLoad','Mutational Load',df$Diagnostic);df$Diagnostic = gsub('FullModel','Full Model',df$Diagnostic);
    df$subgroup = gsub(' ','\n', df$subgroup)
    print(unique(df$subgroup))
    df$PearsonsR = as.numeric(as.character(df$PearsonsR))
    Plot = ggplot(df, aes(x=PearsonsR, y = Diagnostic, color=subgroup)) + geom_point(size=4)+
        #scale_fill_gradientn(colors=c("red","grey","blue"), 
        #    values=rescale(c(min(df$PearsonsR),0, max(df$PearsonsR))),limits=c(min(df$PearsonsR),max(df$PearsonsR)))+
        theme_minimal() + guides(color = guide_legend(ncol = 1))  +
        #facet_grid(rows=vars(subgroup), scales='free_y', space='free_y') +
        labs(y='Model Diagnostic',x="Pearson's R", color='') + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(text = element_text(size=13))+
        ggtitle('PRISM (Drug Screen)') +theme(legend.position='bottom') 
    ggsave(paste0(PlotDir, 'ModelDiagnosticCCLEDrug.pdf'), width=4, height=7, units='in')
}

PlotModelDiagnosticsCCLEshRNA = function(df) {
    df$Diagnostic = gsub('MutLoad','Mutational Load',df$Diagnostic);df$Diagnostic = gsub('FullModel','Full Model',df$Diagnostic);
    #df$Group= gsub(' ','\n', df$Group)
    print(unique(df$subgroup))
    df$PearsonsR = as.numeric(as.character(df$PearsonsR))
    Plot = ggplot(df, aes(x=PearsonsR, y = Diagnostic, color=Group)) +
        #scale_fill_gradientn(colors=c("red","grey","blue"), 
        #    values=rescale(c(min(df$PearsonsR),0, max(df$PearsonsR))),limits=c(min(df$PearsonsR),max(df$PearsonsR)))+
        theme_minimal() + geom_point(size=4)+  theme(legend.position='bottom') +
        #facet_grid(rows=vars(subgroup), scales='free_y', space='free_y') +
        labs(y='Model Diagnostic',x="Pearson's R", color='') + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(text = element_text(size=13))+
        ggtitle('Achilles (shRNA Screen)') + guides(color = guide_legend(ncol = 1)) 
    ggsave(paste0(PlotDir, 'ModelDiagnosticCCLEshRNA.pdf'), width=4, height=7, units='in')
}

PlotModelDiagnosticsCCLEExpression = function(df) {

    df$Diagnostic = gsub('MutLoad','Mutational Load',df$Diagnostic);df$Diagnostic = gsub('FullModel','Full Model',df$Diagnostic);
    #df$Group= gsub(' ','\n', df$Group)
    print(unique(df$subgroup))
    df$PearsonsR = as.numeric(as.character(df$PearsonsR))
    Plot = ggplot(df, aes(x=PearsonsR, y = Diagnostic , color=Group)) + geom_boxplot() +
        #scale_fill_gradientn(colors=c("red","grey","blue"), 
        #    values=rescale(c(min(df$PearsonsR),0, max(df$PearsonsR))),limits=c(min(df$PearsonsR),max(df$PearsonsR)))+
        theme_minimal() + theme(legend.position='bottom') +
        #facet_grid(rows=vars(subgroup), scales='free_y', space='free_y') +
        labs(x="Pearson's R",y='Model Diagnostic', color="") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(text = element_text(size=14))+
        ggtitle('Gene Expression') + guides(color = guide_legend(ncol = 1)) 
    ggsave(paste0(PlotDir, 'ModelDiagnosticCCLEExpression.pdf'), width=4, height=7, units='in')

}


PlotLinearityOfExpVars = function(df, Dataset) {
    # Linearity of explanatory variables - purity, exp, load
    print(head(df))
    df= na.omit(df)
    library(MASS)
    Linear = ggplot(df, aes(x = value)) + 
        geom_histogram(aes(y=..density..),colour = "black", fill = "white", bins=20) +
        geom_density(aes(color = 'emprical'), size = 1)+
        stat_function(fun = ~ dnorm(.value,
                                MASS::fitdistr(.value, 'normal')$estimate[[1]],
                                MASS::fitdistr(.value, 'normal')$estimate[[2]]),
                  aes(color = 'Normal'))
    if (Dataset == 'TCGA') {
        Linear = Linear + facet_wrap(~variable, scale='free') + ggtitle(Dataset)
    }

    ggsave(paste0(PlotDir, 'LinearityOfLoadAndPurity_', Dataset, '.pdf' ), width=6, height=5, units='in')
}


PlotMultiCollinearityOfSNVsAndCNVs = function(df) {
    library('ggcorrplot')
    print(head(df))
    ccle = subset(df, df$Dataset == 'CCLE'); row.names(ccle) = ccle$Barcode; ccle$Dataset=NULL; ccle$Barcode=NULL
    tcga = subset(df, df$Dataset == 'TCGA'); row.names(tcga) = tcga$Barcode; tcga$Dataset=NULL; tcga$Barcode=NULL
    ccle.cor = cor(ccle, method = c("pearson")); tcga.cor = cor(tcga, method = c("pearson"))
    ccle_plot = ggcorrplot(ccle.cor, lab = TRUE) + ggtitle('CCLE') + theme(plot.title = element_text(hjust = 0.5,face="bold", size=12))
    tcga_plot = ggcorrplot(tcga.cor, lab = TRUE) + ggtitle('TCGA') + theme(plot.title = element_text(hjust = 0.5,face="bold", size=12))
    plot_grid(tcga_plot, ccle_plot, labels = c("A", "B"), ncol=1)
    ggsave(paste0(PlotDir, 'CCLE_TCGA_Multicollinearity.pdf' ), width=5, height=7, units='in')
}