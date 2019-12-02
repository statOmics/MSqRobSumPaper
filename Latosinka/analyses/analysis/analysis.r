library(tidyverse)
library(msqrobsum)
library(DEP)
library(MSnbase)
source("../functions.r")

mq_path = '../../data/maxquant/'
output_path = '../output/'
results_path = '../results/'

################################
## Import and preprocess data ##
################################
## txt_path = mq_path

set = read_mq(mq_path) %>% preprocess(norm = 'vsn',logtrans = FALSE)
############
## MSqRob ##
############
form =  c(expression ~ (1|condition) + (1|sample) + (1|feature),expression ~ (1|condition))
out = msqrobsum(set, formulas = form, group_vars = c('protein')
              , contrasts = 'condition', mode = 'msqrob',keep_model = TRUE)

out %>% write_rds(paste0(output_path, 'msqrob'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrob'))

## future::nbrOfWorkers()

## system.time({out = msqrobsum(set, formulas = form, group_vars = c('protein'), contrasts = 'condition', mode = 'msqrob', parallel_args = list(strategy = "sequential"))}) #326
## system.time({out = msqrobsum(set, formulas = form, group_vars = c('protein'), contrasts = 'condition', mode = 'msqrob', parallel_args = list(strategy = "multisession"))}) # 1042.432
## system.time({out = msqrobsum(set, formulas = form, group_vars = c('protein'), contrasts = 'condition', mode = 'msqrob', parallel_args = list(strategy = "multisession", workers = 8))})

## system.time({out = msqrobsum(set, formulas = form, group_vars = c('protein'), contrasts = 'condition', mode = 'msqrob', parallel_args = list(strategy = "multiprocess", workers = 4))})# 1042.432

## packageurl <- 'https://cran.r-project.org/src/contrib/Archive/future/future_1.11.1.1.tar.gz'
## install.packages(packageurl, repos=NULL, type="source")
## packageurl <- 'https://cran.r-project.org/src/contrib/Archive/future.apply/future.apply_1.1.0.tar.gz'
## install.packages(packageurl, repos=NULL, type="source")
## system.time({out = msqrobsum(set, formulas = form, group_vars = c('protein'), contrasts = 'condition', mode = 'msqrob', parallel_args = list(strategy = "multisession"))}) # 855.4
## system.time({out = msqrobsum(set, formulas = form, group_vars = c('protein'), contrasts = 'condition', mode = 'msqrob', parallel_args = list(strategy = "multisession", workers = 8L))})
###############
## MSqRobSum ##
###############
form =  expression ~ (1|condition)
out = msqrobsum(set, formulas = form, group_vars = c('protein')
              , contrasts = 'condition',keep_model = TRUE)

out %>% write_rds(paste0(output_path, 'msqrobsum'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrobsum'))

#############
## Perseus ##
#############
read_tsv('./perseus_S0.txt', comment = '#') %>%
  mutate(protein = as.character(row_number())) %>%
  perseus2res() %>%
  ## mutate(qvalue = p.adjust(pvalue,'BH')) %>%
  write_rds(paste0(results_path, 'perseus'))


## a = read_tsv('./perseus_S0.txt', comment = '#')
## c = a %>%
##   transmute(protein = as.character(row_number())
##           , feature = protein
##           , contrast = 'pTa-pT2+'
##           , sample1 = 'pTa'
##           , sample2 ='pT2+' 
##           , logFC = as.numeric(`N: Welch's T-test Difference pTa_pT2+`)
##           , t = as.numeric(`N: Welch's T-test Test statistic pTa_pT2+`)
##           , pvalue = as.numeric(`N: Welch's T-test p-value pTa_pT2+`)
##           , qvalue = as.numeric(`N: Welch's T-test q-value pTa_pT2+`)
##             )

###############################
## DEP with mixed imputation ##
###############################
protset_DEP = DEP_preprocess(mq_path)
fd = mutate(fData(protset_DEP)
          , name = featureNames(protset_DEP)
          , ID = name)
pd = mutate(pData(protset_DEP)
          , label = paste0(condition,'_',replicate))
e = 2^exprs(protset_DEP) %>% as_data_frame %>% bind_cols(fd,.)
colid = 6:25
depin = make_se(e, colid, pd)

diff <- test_diff(depin, 'all')
dep <- add_rejections(diff, .05, 1)
out <- get_results(dep)

out %>% write_rds(paste0(output_path, 'dep_miximp'))
out %>% dep2res(fasta_path) %>% write_rds(paste0(results_path, 'dep_miximp'))

#####################
## MSqRob quantile ##
#####################
set = read_mq(mq_path) %>% preprocess(norm = 'quantiles',logtrans = FALSE)
form =  c(expression ~ (1|condition) + (1|sample) + (1|feature),expression ~ (1|condition))
out = msqrobsum(set, formulas = form, group_vars = c('protein')
              , contrasts = 'condition', mode = 'msqrob',keep_model = TRUE)

out %>% write_rds(paste0(output_path, 'msqrob_quant'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrob_quant'))

########################
## MSqRobSum quantile ##
########################
set = read_mq(mq_path) %>% preprocess(norm = 'quantiles',logtrans = FALSE)
form =  expression ~ (1|condition)
out = msqrobsum(set, formulas = form, group_vars = c('protein')
              , contrasts = 'condition',keep_model = TRUE)

out %>% write_rds(paste0(output_path, 'msqrobsum_quant'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrobsum_quant'))
