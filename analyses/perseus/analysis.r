library(tidyverse)
library(msqrobsum)
library(DEP)
library(MSnbase)
source("../functions.r")

fasta_path = '../../data/fasta/'
mq_path = '../../data/maxquant/'
output_path = '../output/'
results_path = '../results/'

#####################
## default Perseus ##
#####################
read_tsv('./Perseus_Perseuspreprocess.txt', comment = '#') %>%
    rename(protein = `Protein IDs`) %>%
    perseus2res(fasta_path) %>%
    write_rds(paste0(results_path, 'perseus'))

####################################
## Perseus with our preprocessing ##
####################################
## NOTE THere are some NA log FC with a high qvalue/pvalue
## this happens if one of the 2 conditions in the contrast is completly missing (so no comparision is possible)
## dou our preprocessing on maxlfq values and save to file to use as perseus input
read_maxlfq(mq_path, fasta_path) %>%
    preprocess(norm = 'vsn',logtrans = FALSE) %>%
    exprs %>% as.data.frame %>%
    rownames_to_column('protein') %>% as_data_frame %>%
    write_tsv('./perseus_input.tsv')

read_tsv('./Perseus_mypreprocess_vsn.txt',comment = '#') %>%
    perseus2res(fasta_path) %>%
    write_rds(paste0(results_path, 'perseus_my'))
