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

#################################################33
library(MSqRob)
set2 = set
fData(set2)$protein = as.factor(fData(set2)$protein)
fData(set2)$peptide = rownames(fData(set2))
pData(set2)$run = rownames(pData(set2))
glimpse(fData(set2))
fData(set2) = fData(set2)[,c(1,4)]
glimpse(pData(set2))
pData(set2) = pData(set2)[,c(1,3)]

system.time({
  L <- makeContrast(contrasts= c(
                      "conditionb-conditiona", "conditionc-conditiona", "conditiond-conditiona"
                    , "conditione-conditiona", "conditionc-conditionb", "conditiond-conditionb"
                    , "conditione-conditionb", "conditiond-conditionc", "conditione-conditionc"
                    , "conditione-conditiond"),
                    levels=c("conditiona","conditionb", "conditionc", "conditiond", "conditione"))
  proteinsCPTAC <- MSnSet2protdata(set2, accession="protein")
  modelCPTAC_RR <- fit.model(protdata=proteinsCPTAC, response="quant_value", fixed=c("condition"), random=c("peptide","run"), shrinkage.fixed=NULL, weights="Huber", squeezeVar=TRUE)
  result_old <- test.protLMcontrast(modelCPTAC_RR, L)
  result_old2 <- prot.p.adjust(result_old, method="fdr")
}) ##624

system.time({rmm = msqrobsum(formulas = form_pep, set, group_vars = c('protein', 'ecoli', 'human')
              , contrasts = 'condition', mode = 'msqrob', parallel_args = list(strategy = 'sequential'))}) ## 701

########
form_pep =  c(expression ~ (1|condition) + (1|sample) + (1|feature),expression ~ (1|condition))
rmm = msqrobsum(formulas = form_pep, set, group_vars = c('protein', 'ecoli', 'human')
                      , contrasts = 'condition', mode = 'msqrob')

########
res_old = imap_dfr(result_old,
                   ~{rownames_to_column(.x, 'protein') %>%
                       mutate(contrast = str_replace_all(.y, 'condition',''))}
                   ) %>%
  ## rownames_to_column('protein') %>%
  as_tibble %>%
  transmute(protein, contrast, logFC = estimate, se, t = Tval, pvalue = pval ,qvalue = qval,method = 'old') %>%
  separate(contrast, c('sample1','sample2'),remove = FALSE) %>%
  filter(!is.na(logFC))
res_old = MSnSet2df(set) %>% distinct(protein,ecoli,human) %>% right_join(res_old)

res_new = rmm %>% rmm2res %>% mutate(method = 'new')
res = bind_rows(res_old, res_new)
########################

calculate_fc = function(condition_concentration){
  data_frame(contrast = combn(condition_concentration$condition ,2,
                              function(x) paste(sort(x,decreasing=TRUE),collapse="-"))) %>%
    unnest(condition = str_split(contrast,'-')) %>%
    left_join(condition_concentration) %>%
    group_by(contrast) %>%
    summarise(conditions = list(condition)
            , concentrations = list(concentration)
            , conc_diff = diff(-concentration)
            , logFC_real = diff(-log2(concentration))
            , FC_real = 2^logFC_real
              )
}

preprocess_results = function(results){
  real_fc =
    tribble(~condition, ~concentration
          , 'a', 3
          , 'b', 4.5
          , 'c', 6
          , 'd', 7.5
          , 'e', 9) %>%
    calculate_fc

  ## negagative ecoli FC are also negatives
  results = results %>%
    ## Limma sometimes returns qvalue == NA , remove these
    ## perseus sometimes returns infinite foldchanges , remove these
    ## remove all human + ecoli assignments
    filter(!is.na(qvalue),!is.infinite(logFC),!is.na(logFC)
         , !(human & ecoli)
           ) %>%
    mutate(method = as.factor(method)
         , TP = ifelse(!ecoli | (logFC<0), FALSE,TRUE)
         , contrast = str_glue('{sample1}-{sample2}')
         , species = ifelse(ecoli,'E. coli','Human')) %>%
    group_by(contrast) %>% mutate(TP_total = length(unique(protein[ecoli]))) %>%
    group_by(method, contrast) %>% arrange(pvalue) %>%
    mutate(FDP = cummean(!TP)
         , FPR = cumsum(!TP)/sum(!TP)
         , TPR  = cumsum(TP)/TP_total
         , TP_hits = cumsum(TP)
           ) %>% ungroup %>%
    left_join(real_fc)
}

##########
results = preprocess_results(res)

make_msqrob_plot = function(results,max_TPR = .8){
  max_FDP = .1
  pd = filter(results
            , method %in% c('new', 'old')) %>%
    ## mutate(method = recode_factor(method
    ##                             , 'msqrob' = 'MSqRob', 'msqrobsum' = 'MSqRobSum'
    ##                               )) %>%
    mutate(contrast = factor(contrast,
                             levels = c("b-a", "c-a", "d-a", "e-a", "c-b", "d-b", "d-c", "e-b", "e-c", "e-d")))

  pd = group_by(pd,contrast,FC_real) %>% summarise_at(vars(conditions,concentrations),~{list(.x[[1]])}) %>%
    ungroup %>% transmute(contrast,contrast_info = factor(
                                     str_glue('{contrast}
Concentration (% wt/wt):\t{map_chr(conditions,~(.x[1]))}: {map_dbl(concentrations,~{(.x[1]) %>% round(2)})}     {map_chr(conditions,~last(.x))}: {map_dbl(concentrations,~{last(.x) %>% round(2)})}
Fold Change:\t\t\t\t{round(FC_real,2)}'))) %>%
left_join(pd,.)
  pd$contrast_info = factor(pd$contrast_info, levels = levels(pd$contrast_info)[c(1,2,4,3,5,6,7:10)])

  q = filter(pd,contrast %in% c('b-a','c-a','d-a','c-b','d-b','d-c')) %>%
      mutate('1%' = .01,'5%' = .05) %>% gather(fdr, fdr_value, '1%', '5%') %>%
      group_by(method, contrast,fdr) %>%
      filter(qvalue == max(qvalue[qvalue < fdr_value])) %>% arrange(pvalue) %>%
      dplyr::slice(1)

p = pd %>%
    filter(contrast %in% c('b-a','c-a','d-a','c-b','d-b','d-c')) %>%
    ggplot(aes(FDP,TPR,colour = method))+ geom_path() +
    geom_vline(xintercept = .05, colour = 'grey') +
    geom_vline(xintercept = .01, colour = 'grey') +
    xlim(0,max_FDP) +
    ylim(0,max_TPR) +
      geom_point(aes(FDP,TPR,colour = method,shape = fdr),data = q,width = .002,fill = 'white',
                 height = .002,size = 2)+
    scale_shape_manual(name = 'FDR',values = c('1%' = 21, '5%' = 24)) +
    facet_wrap(~contrast_info,ncol = 3) +
      theme(strip.text = element_text(hjust = 0)) +
    xlab('False Discovery Proportion') +
    ylab('True Positive Rate')
  p
}
p = make_msqrob_plot(results,1);p

