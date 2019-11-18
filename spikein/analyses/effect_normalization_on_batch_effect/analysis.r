library(tidyverse)
library(msqrobsum)
library(MSnbase)
source("../functions.r")

fasta_path = '../../data/fasta/'
mq_path = '../../data/maxquant/'
## output_path = '../output/'
## results_path = '../results/'

################################
## Import and preprocess data ##
################################
set = read_mq(mq_path,fasta_path) %>% preprocess(norm = 'vsn',logtrans = FALSE)
set_nn = read_mq(mq_path,fasta_path) %>% preprocess(norm = FALSE)

###############
## MSqRobSum ##
###############

form_prot =  expression ~ (1|condition)
system.time({
    rmm_rr = msqrobsum(formulas = form_prot, set, group_vars = c('protein', 'ecoli', 'human')
                     , contrasts = 'condition')
})
system.time({
    rmm_nn_rr = msqrobsum(formulas = form_prot, set_nn, group_vars = c('protein', 'ecoli', 'human')
                     , contrasts = 'condition')
})

form_prot_batch =  c(expression ~ (1|condition) + (1|batch), expression ~ (1|condition))
system.time({
    rmm_rr_batch = msqrobsum(formulas = form_prot_batch, set, group_vars = c('protein', 'ecoli', 'human')
                     , contrasts = 'condition')
})
system.time({
    rmm_nn_rr_batch = msqrobsum(formulas = form_prot_batch, set_nn, group_vars = c('protein', 'ecoli', 'human')
                        , contrasts = 'condition')
})

ml = list(rmm_rr =rmm_rr, rmm_nn_rr = rmm_nn_rr,
          rmm_nn_rr_batch = rmm_nn_rr_batch, rmm_rr_batch = rmm_rr_batch)
res = imap_dfr(ml, ~{rmm2res(.x) %>% mutate(method = .y)})

write_rds(res,'results_rds')

#############################################################################################################
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
           ## , FDR.05 = -log(first(pvalue[qvalue == max(qvalue[qvalue < .05])]))
           ) %>% ungroup %>%
    left_join(real_fc)
}


res = read_rds('results_rds')
res = list.files('../results', full.names = TRUE) %>%
    map_dfr(~{print(.x);read_rds(.x) %>% mutate(method = str_replace(.x,'.+\\/([^\\/]+)','\\1'))}) %>%
    filter(!method %in% unique(res$method)) %>% bind_rows(res)

results = preprocess_results(res)

d = list('vsn normalization' = set, 'no normalization' = set_nn) %>%
    imap_dfr(~{limma::plotMDS(exprs(.x), top = Inf)[[3]] %>% cbind(pData(.x)) %>%
                   mutate(method = .y,sample = sampleNames(.x))})
pmds = ggplot(d) + geom_label(aes(`1`,`2`, label = sample,colour = condition)) +
    facet_grid(~method) + guides(colour = FALSE) +
    xlab('Leading logFC dim 1') +
    ylab('Leading logFC dim 2')

make_msqrob_plot = function(results,max_TPR = .8){
  max_FDP = .1
  pd = results %>%
      filter(method %in% c('rmm_rr', 'rmm_rr_batch','rmm_nn_rr', 'rmm_nn_rr_batch')) %>%
    mutate(method = recode_factor(method
                                , 'rmm_nn_rr' = 'MSqRobSum\nno norm.'
                                , 'rmm_nn_rr_batch' = 'MSqRobSum\nno norm. + batch'
                                , 'rmm_rr' = 'MSqRobSum\nVSN norm.'
                                , 'rmm_rr_batch' = 'MSqRobSum\nVSN norm. + batch'
                                  )) %>%
    mutate(contrast = factor(contrast,
                             levels = c("b-a", "c-a", "d-a", "e-a", "c-b", "d-b", "d-c", "e-b", "e-c", "e-d")))

  pd = group_by(pd,contrast,FC_real) %>% summarise_at(vars(conditions,concentrations),~{list(first(.x))}) %>%
    ungroup %>% transmute(contrast,contrast_info = factor(
                                     str_glue('{contrast}
Concentration (% wt/wt):\t{map_chr(conditions,~first(.x))}: {map_dbl(concentrations,~{first(.x) %>% round(2)})}     {map_chr(conditions,~last(.x))}: {map_dbl(concentrations,~{last(.x) %>% round(2)})}
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
p = make_msqrob_plot(results);p

pall = cowplot::plot_grid(pmds,p,rel_heights = c(3,3), labels = c('A','B'),ncol=1)

ggsave('../../figures/vsn_vs_no_norm_batch_effect_ABCD.png', pall, width = 10.5, height = 10);pall
