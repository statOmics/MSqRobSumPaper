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
