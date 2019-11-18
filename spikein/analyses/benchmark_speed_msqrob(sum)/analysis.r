library(tidyverse)
library(msqrobsum)
library(MSnbase)
library(bench)
source("../functions.r")

fasta_path = '../../data/fasta/'
mq_path = '../../data/maxquant/'

set = read_mq(mq_path,fasta_path) %>% preprocess(norm = 'vsn',logtrans = FALSE)

form_prot =  expression ~ (1|condition)
form_pep =  c(expression ~ (1|condition) + (1|sample) + (1|feature),expression ~ (1|condition))

results =
    mark(
        rmm_rr_parallel = {suppressMessages(msqrobsum(formulas = form_prot, set, group_vars = c('protein', 'ecoli', 'human')
                                  , contrasts = 'condition'));1}
        , rmm_rr = {suppressMessages(msqrobsum(formulas = form_prot, set, group_vars = c('protein', 'ecoli', 'human')
                         , contrasts = 'condition', parallel_strategy = 'sequential'));1}
        , rmm_parallel= {suppressMessages(msqrobsum(formulas = form_pep, set, group_vars = c('protein', 'ecoli', 'human')
                              , contrasts = 'condition', mode = 'msqrob'));1}
        , rmm = {suppressMessages(msqrobsum(formulas = form_pep, set, group_vars = c('protein', 'ecoli', 'human')
                      , contrasts = 'condition', mode = 'msqrob', parallel_strategy = 'sequential'));1}
      , iterations= 20L
      , filter_gc = FALSE
      , check =FALSE)

results = read_rds('results.rds')

p = unnest(results) %>%
    mutate(method = as.character(expression), minutes = as.numeric(time)/60
         , method = recode_factor(method
                                , rmm = 'MSqRob'
                                , rmm_parallel = 'MSqRob\nparallel'
                                , rmm_rr = 'MSqRobSum'
                                , rmm_rr_parallel = 'MSqRobSum\nparallel'
                                  )) %>%
    ggplot + geom_boxplot(aes(method, minutes)) + scale_y_continuous(breaks = c(1:10)) +
    theme(axis.title.y=element_blank()) + coord_flip()
ggsave('../../figures/parallel_msqrob_msqrobsum_comparison_times.png',p, width = 8,height = 4);p
