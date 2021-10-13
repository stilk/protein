# Download un-normalized gene expression quantification across high vs low mutational burden patients
# to perform differential gene expression analysis with DESeq2.

#################
### Libraries ###
#################

library(GenomicDataCommons)
library(magrittr)
library(SummarizedExperiment) # Load required packages to download TCGA expression data 
library(TCGAbiolinks) # Website as a reference to get more info on TCGAbiolinks: https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html

###################
### Directories ###
###################


HTSeqCountsDir = "/labs/ccurtis2/tilk/scripts/protein/Data/Raw/Expression/TCGA/HTSeq/" # Directory of raw HT seq counts that are downloaded
OutDir = "/labs/ccurtis2/tilk/scripts/protein/Data/Raw/Expression/TCGA/HTSeq-Concatenated/" # Final direcotry of concatenated expr profiles
#################
### Functions ###
#################


GetAllCancerTypes <- function() {
    # Returns a dataframe that contains the cancer type and barcode of low and high mutational burden tumors
    samples = read.table(paste0(getwd(), '/Data/Raw/Mutations/TCGA/CancerTypesAndMutLoadForDGE'), sep=',', header=T)
    return(unique(samples$type))
}


TranslateTCGAUUIDs <- function(file_ids, legacy = FALSE) {
    # Takes an input of TCGA UUIDs and returns a dataframe that maps them to TCGA Barcodes
    info = files(legacy = legacy) %>%
        filter( ~ file_id %in% file_ids) %>%
        select('cases.samples.submitter_id') %>%
        results_all()
    # The mess of code below is to extract TCGA barcodes
    # id_list will contain a list (one item for each file_id)
    # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
    id_list = lapply(info$cases,function(a) {
        a[[1]][[1]][[1]]})
    # so we can later expand to a data.frame of the right size
    barcodes_per_file = sapply(id_list,length)
    # And build the data.frame
    return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                      submitter_id = unlist(id_list), stringsAsFactors=FALSE))
    }

ConcatenateSamples <- function(UUID_Files) {
    # Takes in a path of HT-Seq files downloaded from GDC for a particular cancer type
    # Downloaded files are deposited according to UUIDs
    # This function merges all downloaded UUID files by gene ID 
    # and adds the TGGA barcode ID associated with the sample as a column.
    UUID_To_Barcode = TranslateTCGAUUIDs(UUID_Files) ### Get dataframe of UUIDs to TCGA barcodes
    output = data.frame()
    for (UUID in UUID_Files) {
        fileName = Sys.glob(file.path(UUID, "*"))
        sampleCounts = read.table(fileName, header=F, sep='\t')
        TGCA_Barcode = subset(UUID_To_Barcode, UUID_To_Barcode$file_id == UUID)[c('submitter_id')]
        names(sampleCounts) = c('ENSG',TGCA_Barcode)
        if (length(output) == 0) {
          output = sampleCounts
        } else {
          output = merge(output, sampleCounts, by=c('ENSG'))
        }
    }
    return(output)
}

GetHTSeqCountsForPatients <- function() {
    # Writes to file HT-seq counts for each cancer type
    for (CancerType in GetAllCancerTypes()) { ### loop through all the cancer types
        SampleCancerType = as.vector(subset(samples, samples$type == CancerType)[c('Barcode')]) ### get the TCGA barcodes
        ### query GDC for HT-Seq files that exist for a particular cancer type
        tryCatch(
          {
            query.exp.proj.gene = GDCquery(project = paste0("TCGA-",CancerType),
                                                    legacy = FALSE,
                                                    data.category = "Transcriptome Profiling",
                                                    data.type = "Gene Expression Quantification",
                                                    workflow.type = "HTSeq - Counts"
                                                )
            GDCdownload(query.exp.proj.gene, directory = HTSeqCountsDir) ### download the raw HT-Seq count data
            setwd(paste0(HTSeqCountsDir, "TCGA-", CancerType, "/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/"))
            UUID_Files = gsub("\\./*", "", list.dirs(recursive=FALSE)) ### Get the list of UUIDs that were downloaded 
            MergedCountsByCancerType = ConcatenateSamples(UUID_Files)
            write.table(MergedCountsByCancerType, paste0(OutDir, CancerType), quote = FALSE, sep = "\t")
              
          }, error=function(error_message) {
                message(error_message)
                message('Skipping this cancer type.')
          }
        )
    }
}

