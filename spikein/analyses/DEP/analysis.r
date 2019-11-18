library(tidyverse)
library(msqrobsum)
library(DEP)
library(MSnbase)
source("../functions.r")

fasta_path = '../../data/fasta/'
mq_path = '../../data/maxquant/'
output_path = '../output/'
results_path = '../results/'

###############################
## DEP with mixed imputation ##
###############################
protset_DEP = DEP_preprocess(mq_path,fasta_path)
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


########################################
## DEP with mixed imputation pairwise ##
########################################
## protset_DEP = DEP_preprocess(mq_path,fasta_path)
set_DEP = DEP_read(mq_path)
ids = combn(unique(set_DEP@colData$condition),2,simplify = FALSE)
## protset_DEP = set_DEP[,set_DEP@colData$condition %in% i] %>% DEP_process(fasta_path)

a = Sys.time()
out =map(ids,~{
  print(.)
  protset_DEP= set_DEP[,set_DEP@colData$condition %in% .] %>% DEP_process(fasta_path)

  fd = mutate(fData(protset_DEP)
            , name = featureNames(protset_DEP)
            , ID = name)
  pd = mutate(pData(protset_DEP)
            , label = paste0(condition,'_',replicate))
  e = 2^exprs(protset_DEP) %>% as_data_frame %>% bind_cols(fd,.)
  colid = 6:13
  depin = make_se(e, colid, pd)

  diff <- test_diff(depin, 'all')
  dep <- add_rejections(diff, .05, 1)
  get_results(dep)
})
b = Sys.time()
b-a
## > > Time difference of 17.30088 secs

out %>% write_rds(paste0(output_path, 'dep_miximp_pairwise'))
out %>% map_dfr(~dep2res(.,fasta_path)) %>% write_rds(paste0(results_path, 'dep_miximp_pairwise'))
