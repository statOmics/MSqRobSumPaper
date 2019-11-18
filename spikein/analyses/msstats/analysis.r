library(tidyverse)
library(msqrobsum)
## library(MSnbase)
library(MSstats)
source("../functions.r")

fasta_path = '../../data/fasta/'
mq_path = '../../data/maxquant/'
output_path = '../output/'
results_path = '../results/'

##############################
## default msstats analysis ##
##############################
## Read in MaxQuant files
proteinGroups <- read.table(str_glue("{mq_path}proteinGroups.txt"), sep="\t", header=TRUE)
infile <- read.table(str_glue("{mq_path}evidence.txt"), sep="\t", header=TRUE)
## Read in annotation including condition and biological replicates per run.
## Users should make this annotation file. It is not in the output from MaxQuant.
annot = data_frame(Run = paste0(rep(letters[1:5],each = 4), rep(1:4,5)), Condition = rep(letters[1:5],each = 4)
             , BioReplicate = 1, IsotopeLabelType = 'L')

annot = transmute(infile,Raw.file, Run = Experiment) %>% distinct %>% left_join(annot)

msstat_input <- MaxQtoMSstatsFormat(evidence=infile,
                           annotation=annot,
                           proteinGroups=proteinGroups)

dproc_default = dataProcess(msstat_input)

form_pep = expression ~ (1|condition) + (1|sample) + (1|feature)
contrasts_msstat = t(msqrobsum:::make_simple_contrast(msstat_input,'Condition'))
row.names(contrasts_msstat) = str_replace_all(row.names(contrasts_msstat),'Condition', '')
out = groupComparison(contrast.matrix=contrasts_msstat, data = dproc_default)

## out %>% write_rds(paste0(output_path, 'msstats')) ## Huge file
out %>% msstats2res(fasta_path) %>% write_rds(paste0(results_path, 'msstats'))



####################################
## msstats analysis no imputation ##
####################################
## Read in MaxQuant files
proteinGroups <- read.table(str_glue("{mq_path}proteinGroups.txt"), sep="\t", header=TRUE)
infile <- read.table(str_glue("{mq_path}evidence.txt"), sep="\t", header=TRUE)
## Read in annotation including condition and biological replicates per run.
## Users should make this annotation file. It is not in the output from MaxQuant.
annot = data_frame(Run = paste0(rep(letters[1:5],each = 4), rep(1:4,5)), Condition = rep(letters[1:5],each = 4)
             , BioReplicate = 1, IsotopeLabelType = 'L')

annot = transmute(infile,Raw.file, Run = Experiment) %>% distinct %>% left_join(annot)

msstat_input <- MaxQtoMSstatsFormat(evidence=infile,
                           annotation=annot,
                           proteinGroups=proteinGroups)

dproc_default = dataProcess(msstat_input
                          , MBimpute = FALSE
                          , censoredInt = NULL
                            )

form_pep = expression ~ (1|condition) + (1|sample) + (1|feature)
contrasts_msstat = t(msqrobsum:::make_simple_contrast(msstat_input,'Condition'))
row.names(contrasts_msstat) = str_replace_all(row.names(contrasts_msstat),'Condition', '')
out = groupComparison(contrast.matrix=contrasts_msstat, data = dproc_default)

## out %>% write_rds(paste0(output_path, 'msstats')) ## Huge file
out %>% msstats2res(fasta_path) %>% write_rds(paste0(results_path, 'msstats_noImp'))

####################################
## msstats analysis no norm ##
####################################
## Read in MaxQuant files
proteinGroups <- read.table(str_glue("{mq_path}proteinGroups.txt"), sep="\t", header=TRUE)
infile <- read.table(str_glue("{mq_path}evidence.txt"), sep="\t", header=TRUE)
## Read in annotation including condition and biological replicates per run.
## Users should make this annotation file. It is not in the output from MaxQuant.
annot = data_frame(Run = paste0(rep(letters[1:5],each = 4), rep(1:4,5)), Condition = rep(letters[1:5],each = 4)
             , BioReplicate = 1, IsotopeLabelType = 'L')

annot = transmute(infile,Raw.file, Run = Experiment) %>% distinct %>% left_join(annot)

msstat_input <- MaxQtoMSstatsFormat(evidence=infile,
                           annotation=annot,
                           proteinGroups=proteinGroups)

dproc_default = dataProcess(msstat_input
                            , normalization = FALSE
                          ## , MBimpute = FALSE
                          ## , censoredInt = Null
                            )

form_pep = expression ~ (1|condition) + (1|sample) + (1|feature)
contrasts_msstat = t(msqrobsum:::make_simple_contrast(msstat_input,'Condition'))
row.names(contrasts_msstat) = str_replace_all(row.names(contrasts_msstat),'Condition', '')
out = groupComparison(contrast.matrix=contrasts_msstat, data = dproc_default)

## out %>% write_rds(paste0(output_path, 'msstats')) ## Huge file
out %>% msstats2res(fasta_path) %>% write_rds(paste0(results_path, 'msstats_noNorm'))

############################################
## msstats analysis no imputation no norm ##
############################################
## Read in MaxQuant files
proteinGroups <- read.table(str_glue("{mq_path}proteinGroups.txt"), sep="\t", header=TRUE)
infile <- read.table(str_glue("{mq_path}evidence.txt"), sep="\t", header=TRUE)
## Read in annotation including condition and biological replicates per run.
## Users should make this annotation file. It is not in the output from MaxQuant.
annot = data_frame(Run = paste0(rep(letters[1:5],each = 4), rep(1:4,5)), Condition = rep(letters[1:5],each = 4)
             , BioReplicate = 1, IsotopeLabelType = 'L')

annot = transmute(infile,Raw.file, Run = Experiment) %>% distinct %>% left_join(annot)

msstat_input <- MaxQtoMSstatsFormat(evidence=infile,
                           annotation=annot,
                           proteinGroups=proteinGroups)

dproc_default = dataProcess(msstat_input
                          , normalization = FALSE
                          , MBimpute = FALSE
                          , censoredInt = NULL
                            )

form_pep = expression ~ (1|condition) + (1|sample) + (1|feature)
contrasts_msstat = t(msqrobsum:::make_simple_contrast(msstat_input,'Condition'))
row.names(contrasts_msstat) = str_replace_all(row.names(contrasts_msstat),'Condition', '')
out = groupComparison(contrast.matrix=contrasts_msstat, data = dproc_default)

## out %>% write_rds(paste0(output_path, 'msstats')) ## Huge file
ut %>% msstats2res(fasta_path) %>% write_rds(paste0(results_path, 'msstats_noImp_noNorm'))
