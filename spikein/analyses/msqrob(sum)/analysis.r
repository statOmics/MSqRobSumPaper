library(tidyverse)
library(msqrobsum)
library(DEP)
library(MSnbase)
source("../functions.r")

fasta_path = '../../data/fasta/'
mq_path = '../../data/maxquant/'
output_path = '../output/'
results_path = '../results/'

################################
## Import and preprocess data ##
################################
set = read_mq(mq_path,fasta_path) %>% preprocess(norm = 'vsn',logtrans = FALSE)

############
## MSqRob ##
############
form =  c(expression ~ (1|condition) + (1|sample) + (1|feature),expression ~ (1|condition))
out = msqrobsum(set, formulas = form, group_vars = c('protein', 'ecoli', 'human')
              , contrasts = 'condition', mode = 'msqrob')

out %>% write_rds(paste0(output_path, 'msqrob'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrob'))

###############
## MSqRobSum ##
###############
form =  expression ~ (1|condition)
out = msqrobsum(set, formulas = form, group_vars = c('protein', 'ecoli', 'human')
              , contrasts = 'condition')

out %>% write_rds(paste0(output_path, 'msqrobsum'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrobsum'))

###########################
## MSqRob without Q9BZJ0 ##
###########################
set2 = set[fData(set)$protein != 'Q9BZJ0',]
form =  c(expression ~ (1|condition) + (1|sample) + (1|feature),expression ~ (1|condition))
out = msqrobsum(set2, formulas = form, group_vars = c('protein', 'ecoli', 'human')
              , contrasts = 'condition', mode = 'msqrob')

out %>% write_rds(paste0(output_path, 'msqrob_wo_Q9BZJ0'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrob_wo_Q9BZJ0'))

##############################
## MSqRobSum without Q9BZJ0 ##
##############################
form =  expression ~ (1|condition)
out = msqrobsum(set2, formulas = form, group_vars = c('protein', 'ecoli', 'human')
              , contrasts = 'condition')

out %>% write_rds(paste0(output_path, 'msqrobsum_wo_Q9BZJ0'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrobsum_wo_Q9BZJ0'))

######################################
## MSqRobSum with DEP preprocessing ##
######################################
protset_DEP = DEP_preprocess(mq_path,fasta_path)
form =  expression ~ (1|condition)
out = msqrobsum(protset_DEP ,formulas = form, group_vars = c('protein', 'ecoli', 'human')
              , contrasts = 'condition')

out %>% write_rds(paste0(output_path, 'msqrobsum_DEP'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrobsum_DEP'))

#####################################
## MSqRobSum with mixed imputation ##
#####################################
## two sample condition to false because missing values will be imputed
## Use MSnbase robust summmarization implementation because I wnat to compare imputation
##before and after summarization
form =  expression ~ (1|condition)
## Impute on peptide level
protset_dep_rr = read_mq(mq_path, fasta_path) %>%
    preprocess(logtrans = FALSE,norm = 'vsn', two_sample_condition = FALSE) %>%
    DEP_impute %>%
    combineFeatures(., fun="robust", groupBy = fData(.)$protein,cv = FALSE)

out = msqrobsum(formulas = form, protset_dep_rr , group_vars = c('protein', 'ecoli', 'human')
              , contrasts = 'condition', mode = 'msqrob')

out %>% write_rds(paste0(output_path, 'msqrobsum_miximp_pep'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrobsum_miximp_pep'))

## Impute on protein level
protset_rr_dep = read_mq(mq_path, fasta_path) %>%
    preprocess(logtrans = FALSE,norm = 'vsn', two_sample_condition = FALSE) %>%
    combineFeatures(., fun="robust", groupBy = fData(.)$protein,cv = FALSE) %>%
    DEP_impute

out = msqrobsum(formulas = form, protset_rr_dep , group_vars = c('protein', 'ecoli', 'human')
              , contrasts = 'condition', mode = 'msqrob')

out %>% write_rds(paste0(output_path, 'msqrobsum_miximp_prot'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrobsum_miximp_prot'))

##################################
## MSqRobSum with MaxLFQ values ##
##################################
form =  expression ~ (1|condition)
protset_maxlfq = read_maxlfq(mq_path, fasta_path) %>%
    preprocess(norm = 'vsn',logtrans = FALSE)

out = msqrobsum(formulas = form, protset_maxlfq , group_vars = c('protein', 'ecoli', 'human')
              , contrasts = 'condition', mode = 'msqrob')

out %>% write_rds(paste0(output_path, 'msqrobsum_maxlfq'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrobsum_maxlfq'))


#######################
## MSqRob on pairwise ##
#######################
set = read_mq(mq_path,fasta_path)
ids = combn(unique(pData(set)$condition)[1:4],2,simplify = FALSE)
form =  c(expression ~ (1|condition) + (1|sample) + (1|feature),expression ~ (1|condition))

a = Sys.time()
out =map_dfr(ids,~{
  print(.)
  subset = set[,pData(set)$condition %in% .]  %>% preprocess(norm = 'vsn',logtrans = FALSE)
  msqrobsum(subset, formulas = form, group_vars = c('protein', 'ecoli', 'human')
          , contrasts = 'condition', mode = 'msqrob',parallel_args = list(strategy = "sequential"))
})
b = Sys.time()
b-a
## > > Time difference of 1.144766 hours

## select(out,contrasts) %>% unnest %>% count(contrast)
out %>% write_rds(paste0(output_path, 'msqrob_pairwise2'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrob_pairwise2'))

########################
## MSqRobSum pairwise ##
########################
set = read_mq(mq_path,fasta_path)
ids = combn(unique(pData(set)$condition)[1:4],2,simplify = FALSE)

form =  expression ~ (1|condition)
a = Sys.time()
out =map_dfr(ids,~{
  print(.)
  subset = set[,pData(set)$condition %in% .]  %>% preprocess(norm = 'vsn',logtrans = FALSE)
  msqrobsum(subset, formulas = form, group_vars = c('protein', 'ecoli', 'human')
          , contrasts = 'condition',parallel_args = list(strategy = "sequential"))
})
b = Sys.time()
b-a
## > > Time difference of 30.77726 mins

out %>% write_rds(paste0(output_path, 'msqrobsum_pairwise2'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrobsum_pairwise2'))
