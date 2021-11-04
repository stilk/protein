
# This script does differential gene expression analysis on high vs low mutational burden tumors
# as an agnostic way of trying to understand how high mutational burden tumors overcome their
# mutational load. The goal is to find genes that are highly expressed in high mutational burden
# tumors, but weakly expressed in low mutational burden tumors.
#

library(DESeq2)
library(dplyr)


HT_Seq_Dir="/labs/ccurtis2/tilk/scripts/protein/Data/Raw/Expression/TCGA/HTSeq-Concatenated/"
DGE_Input="/labs/ccurtis2/tilk/scripts/protein/Data/DGE/Input/" 
DGE_Output="/labs/ccurtis2/tilk/scripts/protein/Data/DGE/Output/"

#################
### Functions ###
#################

# Returns a dataframe that contains the cancer type and barcode of low and high mutational burden tumors
GetAllSamples = function() { 
    samples = read.table(paste0(getwd(), '/Data/Raw/Mutations/TCGA/CancerTypesAndMutLoadForDGE'), sep=',', header=TRUE)
    return(samples)
}

# Returns a data frame of barcodes for each condition type (low/high)
GetBarcodesOfInterest = function(Condition) {
  Samples = GetAllSamples()
  if (Condition == 'PanCancer') {
      out = rbind( data.frame(Condition='high', Samples[Samples[c('MutLoad')] > 1000,]),
             data.frame(Condition='low', Samples[Samples[c('MutLoad')] < 10,]))
  } else if (Condition == 'PanCancerMatched') {
    Samples$Low = (Samples$MutLoad <= 50)
    Samples$High = (Samples$MutLoad >= 1000)
    Samples = subset(Samples, (Samples$Low == TRUE) | (Samples$High == TRUE))
    Samples$Condition = ifelse(Samples$Low == 'TRUE','Low', 'High')
    NumSamples = as.data.frame.matrix(table(out[c('type','Condition')]))
    NumSamples = cbind(data.frame(MaxSamplesPerBothGroup = apply(NumSamples, 1, FUN = min, na.rm = TRUE)),
                            data.frame(type=row.names(data.frame(apply(NumSamples, 1, FUN = min, na.rm = TRUE)))))
    Samples = merge(Samples, NumSamples, by='type', all.x=TRUE)
    Samples = subset(Samples, Samples$MaxSamplesPerBothGroup != 0 )
    # Select maximum # of barcodes for each cancer type that match conditions for low and high
    out = data.frame(
            rbind(Samples %>% filter(Condition == "High") %>% group_by(type, MaxSamplesPerBothGroup) %>% top_n(MaxSamplesPerBothGroup, MutLoad) %>% arrange(MutLoad, .by_group = TRUE),
            Samples %>% filter(Condition == "Low") %>% group_by(type, MaxSamplesPerBothGroup) %>% top_n(-MaxSamplesPerBothGroup, MutLoad) %>% arrange(MutLoad, .by_group = TRUE)))
    # For ties of mut load sample one randomly that meets condition
    out = data.frame(out %>% group_by(Condition, type, MaxSamplesPerBothGroup) %>% filter(1:n() <= MaxSamplesPerBothGroup)) 
  }
  return(out)
}


# Writes individual barcodes that match conditions of interest into separate dirs to run DGE.
# @Barcodes = barcodes to split raw htseq counts by 
# @Directory = directory where each barcode is placed into
SeparateBarcodesIntoDirectoriesOfInterest = function(Barcodes, Directory) {
    for (FileName in list.files(HT_Seq_Dir)) { # Loop through all raw HT-seq files concatenated by cancer type
        print(paste0('Reading in: ', FileName))
        htseq = read.table(paste0(HT_Seq_Dir, FileName), sep='\t', header=TRUE, row.names="ENSG")
        names(htseq) = gsub('\\.','-',substr(colnames(htseq), 1, 15)) # Set barcode col names as the same as in conditions
        row.names(htseq) = do.call(rbind,strsplit(row.names(htseq), '\\.'))[,1] # Remove ENSG version # (the . after the ENSG ID)
        SamplesOfInterest = htseq[colnames(htseq) %in% Barcodes] # Dataframe of all Barcodes that match 
        if (ncol(SamplesOfInterest) == 0) {
           next # Skip cancer type if no samples match
        }
        dir.create(Directory, showWarnings = FALSE) # Create dir if doesnt exist
        for (i in 1:ncol(SamplesOfInterest)) { # Loop through each barcode that matches and write to file
          write.table(tail(SamplesOfInterest[i], -5), paste0(Directory, colnames(SamplesOfInterest[i])), 
              col.names=FALSE, row.names=TRUE, quote=FALSE)
        }
    }
}


# Returns a dataframe for DGE input that contains:
# 1) sampleName = file name for the sample
# 2) fileName = file paths to raw HT-Seq counts for every sample
# 3) condition = either low or high mutational burden tumor
GetSampleTable = function(Directory) {
    SampleTable = data.frame(FileName=Sys.glob(file.path(paste0(Directory,"/*"), "*")))
    SampleTable$SampleName = do.call(rbind,strsplit(as.character(SampleTable$FileName), "/"))[,12]
    SampleTable$Condition = do.call(rbind,strsplit(as.character(SampleTable$FileName), "/"))[,11]
    SampleTable$SampleName = paste0("/", SampleTable$Condition, "/", SampleTable$SampleName)
    return(SampleTable)
}

# To run DEseq2, 2 parameters are required:
# @SampleTable = a dataframe of columns containing 
#               `SampleName` :  file name for the sample
#               `FileName` : file paths to raw HT-Seq counts for every sample
#               `Condition` : either low or high mutational burden tumor
# @Directory = a directory of expression counts for each individual sample
DoDGE <- function(SampleTable, Directory) {
    ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = SampleTable, directory = Directory, design= ~ Condition)
    dds = DESeq(ddsHTSeq)
    results = results(dds)
    return(results)
}


GetDGEResults = function(Condition) {
    Directory = paste0(DGE_Input, Condition)
    df = GetBarcodesOfInterest(Condition)
    #SeparateBarcodesIntoDirectoriesOfInterest(df[df$Condition == 'high',]$Barcode, paste0(Directory, '/high/'))
    #SeparateBarcodesIntoDirectoriesOfInterest(df[df$Condition == 'low',]$Barcode, paste0(Directory, '/low/'))
    SampleTable = GetSampleTable(Directory)
    Results = DoDGE(SampleTable, Directory)
    return(Results)
}

# write.table(GetDGEResults('PanCancerMatched'), paste0(DGE_Output, 'TCGA_PanCancerMatched_DGE.txt'), quote=FALSE)
write.table(GetDGEResults('PanCancer'), paste0(DGE_Output, 'TCGA_PanCancer_DGE.txt'), quote=FALSE)