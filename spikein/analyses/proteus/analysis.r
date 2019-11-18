library(tidyverse)
library(proteus)
library(MSnbase)
source("../functions.r")

fasta_path = '../../data/fasta/'
mq_path = '../../data/maxquant/'
output_path = '../output/'
results_path = '../results/'

####################################
## Proteus with vsn normalization ##
####################################
meta = data_frame(condition = rep(letters[1:5], each = 4)
                , replicate = rep(1:4,5)
                , sample = paste0(condition,replicate)
                , experiment = sample
                , measure = 'Intensity')

evi <- readEvidenceFile(paste0(mq_path,"evidence.txt"))

pepdat <- makePeptideTable(evi, meta)
prodat <- makeProteinTable(pepdat)
## use vsn normalization
prodat.vsn <- normalizeData(prodat,limma::normalizeVSN)
write_rds(prodat.vsn,'proteus_summarized_data.rds')
contrasts = list(c('b','a'), c('c','a'), c('d','a'), c('e','a')
               , c('c','b'), c('d','b'), c('e','b')
               , c('d','c'), c('e','c'), c('e','d'))
out = map_dfr(contrasts, ~{
    limmaDE(prodat.vsn,conditions = .x,transform.fun = identity) %>% mutate(sample1 = .x[1], sample2 = .x[2])
})

proteus2res = function(df,fasta_path){
  as_data_frame(res) %>%
    transmute(protein,feature = protein, sample1, sample2
            , logFC,t, pvalue = P.Value,qvalue = adj.P.Val
            , contrast = paste0(sample1,'-',sample2)) %>%
    add_id(fasta_path = fasta_path)
}

out %>% proteus2res(fasta_path) %>% write_rds(paste0(results_path, 'proteus'))
