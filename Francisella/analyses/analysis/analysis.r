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
set = read_mq(mq_path) %>% preprocess(norm = 'vsn',logtrans = FALSE)

############
## MSqRob ##
############
form =  c(expression ~ (1|condition) +(1|biorep) + (1|sample) + (1|feature),
          expression ~ (1|condition) +(1|biorep))

system.time({
out = msqrobsum(set, formulas = form, group_vars = c('protein')
              , contrasts = 'condition', mode = 'msqrob')
}) # 592.037

out %>% write_rds(paste0(output_path, 'msqrob'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrob'))

out = read_rds(paste0(output_path, 'msqrob'))
glimpse(out)
count(out,formula)

############
## MSqRob quant norm ##
############
set = read_mq(mq_path) %>% preprocess(norm = 'quantiles',logtrans = TRUE)

set = read_mq(mq_path) %>% preprocess(norm = 'vsn',logtrans = FALSE)
set = read_mq(mq_path) %>% preprocess(norm = FALSE,logtrans = TRUE)
MSnSet2df(set) %>% ggplot + geom_density(aes(expression,group = sample,colour = condition))
MSnSet2df(set) %>% ggplot + geom_density(aes(expression,group = sample,colour = biorep))

hist(exprs(set))

form =  c(expression ~ (1|condition) +(1|biorep) + (1|sample) + (1|feature),
          expression ~ (1|condition) +(1|biorep))
system.time({
  out = msqrobsum(set, formulas = form, group_vars = c('protein')
                , contrasts = 'condition', mode = 'msqrob')
}) # 592.037

out %>% write_rds(paste0(output_path, 'msqrob_quant'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrob_quant'))

###############
## MSqRobSum ##
###############
form =  c(expression ~ (1|condition) +(1|biorep))
system.time({
out = msqrobsum(set, formulas = form, group_vars = c('protein')
              , contrasts = 'condition')
}) #515
out %>% write_rds(paste0(output_path, 'msqrobsum'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrobsum'))


###############
## MSqRobSum quant norm ##
###############
set = read_mq(mq_path) %>% preprocess(norm = 'quantiles',logtrans = TRUE)
hist(exprs(set))
form =  c(expression ~ (1|condition) +(1|biorep))
system.time({
  out = msqrobsum(set, formulas = form, group_vars = c('protein')
                , contrasts = 'condition')
}) #515
out %>% write_rds(paste0(output_path, 'msqrobsum_quant'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrobsum_quant'))


############################
## ## MSqRob no shrinkage ##
############################

form =  c(expression ~ condition +(1|biorep) + (1|sample) + (1|feature),
          expression ~ condition +(1|biorep))
system.time({
  out = msqrobsum(set, formulas = form, group_vars = c('protein')
                , mode = 'msqrob',keep_model = TRUE)
}) # 592.037

out = filter(out, !is.na(formula)) %>% 
  mutate(
    contrast = 'conditionWT-conditionKO'
  , coefsum = map(model,~{summary(.)$coefficients[2,]})
  , logFC = coefsum[[1]]['Estimate']
  , sigma_contrast = map_dbl(model,~{sqrt(vcov(.)[2,2])})
  ## , sigma_contrast = coefsum[[1]]['Std. Error']/sigma
  , se = sigma_contrast * sigma_post, t = logFC/se
  , pvalue = pt(-abs(t), df_post) * 2
  , qvalue = p.adjust(pvalue,'BH')
  ) %>% group_by(protein) %>%
  mutate(contrasts = list(tibble( contrast, logFC, sigma_contrast,se,pvalue,qvalue))) %>%
  bind_rows(filter(out, is.na(formula))) %>% ungroup

out %>% write_rds(paste0(output_path, 'msqrob_noshrink'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrob_noshrink'))


############################
## MSqRobSum no shrinkage ##
############################
form =  c(expression ~ condition +(1|biorep))
system.time({
  out = msqrobsum(set, formulas = form, group_vars = c('protein'), keep_model = TRUE)
}) #515

out =
  filter(out, !is.na(formula)) %>%
  mutate(
    contrast = 'conditionWT-conditionKO'
  , coefsum = map(model,~{summary(.)$coefficients[2,]})
  , logFC = coefsum[[1]]['Estimate']
  , sigma_contrast = map_dbl(model,~{sqrt(vcov(.)[2,2])})
    ## , sigma_contrast = coefsum[[1]]['Std. Error']/sigma
  , se = sigma_contrast * sigma_post, t = logFC/se
  , pvalue = pt(-abs(t), df_post) * 2
  , qvalue = p.adjust(pvalue,'BH')
  ) %>% group_by(protein) %>%
  mutate(contrasts = list(tibble( contrast, logFC, sigma_contrast,se,pvalue,qvalue))) %>%
  bind_rows(filter(out, is.na(formula))) %>% ungroup

out %>% write_rds(paste0(output_path, 'msqrobsum_noshrinkage'))
out %>% rmm2res %>% write_rds(paste0(results_path, 'msqrobsum_noshrinkage'))

##############
## original ##
##############
readxl::read_xlsx('../../data/mcp.M115.055897-2.xlsx',sheet = 2) %>%
  transmute(protein = accession
          , feature = protein
          , contrast = 'WT-KO'
          , sample1 = 'WT'
          , sample2 = 'KO'
          , logFC = as.numeric(`log2 FC estimate`)
          , t = NA
          , pvalue = as.numeric(pval)
          , qvalue = p.adjust(pvalue, 'BH')
            ) %>%
  as_tibble  %>%
  write_rds(paste0(results_path, 'original'))


readxl::read_xlsx('../../data/mcp.M115.055897-2.xlsx',sheet = 1) %>%
  transmute(protein = accession
          , feature = protein
          , contrast = 'WT-KO'
          , sample1 = 'WT'
          , sample2 = 'KO'
          , logFC = as.numeric(`log2 FC estimate`)
          , t = NA
          , pvalue = as.numeric(pval)
          , qvalue = p.adjust(pvalue, 'BH')
            ) %>%
  as_tibble  %>%
  write_rds(paste0(results_path, 'original_ludger'))

