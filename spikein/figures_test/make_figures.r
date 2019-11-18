library(tidyverse)

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
glimpse(results)
##########
res = list.files('../analyses/results', full.names = TRUE) %>%
    map_dfr(~{print(.x);read_rds(.x) %>% mutate(method = str_replace(.x,'.+\\/([^\\/]+)','\\1'))})

results = preprocess_results(res)

results_filtered = group_by(res,protein, contrast) %>% filter(n() == nlevels(results$method)) %>%
  group_by(method, contrast) %>% mutate(qvalue = p.adjust(pvalue, method = 'BH')) %>% ungroup
results_filtered = preprocess_results(results_filtered)

## Introduction plots
########################
## abcd
##########
make_intro_plot = function(results,max_TPR = .7){
  max_FDP = .1
  pd = filter(results
            , method %in% c('msqrob','perseus'
                           ,'msstats','dep_miximp','proteus')) %>%
    mutate(method = recode_factor(method
                                , 'msqrob' = 'MSqRob'
                                , 'perseus' = 'Perseus\n(maxLFQ)'
                                , 'proteus' = 'proteus\n(high-flyers)'
                                , 'dep_miximp' = 'DEPmiximp\n(maxLFQ)'
                                , 'msstats' = 'MSstats\n(median polish)'
                                  )
         , contrast = factor(contrast,
                             levels = c("b-a", "c-a", "d-a", "e-a", "c-b", "d-b", "d-c", "e-b", "e-c", "e-d")))

  pd = group_by(pd,contrast,FC_real) %>% summarise_at(vars(conditions,concentrations),~{list(first(.x))}) %>%
    ungroup %>% transmute(contrast,contrast_info = factor(
                                     str_glue('{contrast}
Concentration (% wt/wt): {map_chr(conditions,~first(.x))} = {map_dbl(concentrations,~{first(.x) %>% round(2)})}%, {map_chr(conditions,~last(.x))} = {map_dbl(concentrations,~{last(.x) %>% round(2)})}%
Fold Change ({map_chr(conditions,~first(.x))}/{map_chr(conditions,~last(.x))}) = {round(FC_real,2)}'))) %>%
left_join(pd,.)
  pd$contrast_info = factor(pd$contrast_info, levels = levels(pd$contrast_info)[c(1,2,4,3,5,6,7:10)])

  p1 = ungroup(pd) %>%
    ## filter(contrast %in% c('b-a','c-a','d-a','c-b','d-b','d-c')) %>%
    filter(contrast %in% c('b-a')) %>%
    ggplot() + geom_hline(yintercept = 0) +
    geom_hline(aes(yintercept = logFC_real),colour = 'grey',size = 1.5) +
    geom_boxplot(aes(method,logFC,colour = method), outlier.size = 1) +
    facet_wrap(contrast~species,ncol = 2) +
    ylim(-3,3) +
    ## xlab('Method')  +
    xlab(NULL)  +
    ## guides(colour = FALSE) +
    ylab('Log Fold Change') +theme_bw() +
    theme(legend.position = 'bottom') +
    guides(colour = guide_legend(nrow = 1))

  p2 = pd %>%
    filter(contrast %in% c('b-a','c-b','d-c')) %>%
    ggplot(aes(FDP,TPR,colour = method))+ geom_path() +
    facet_wrap(~contrast_info,ncol = 3) +
    xlim(0,max_FDP) +
    ylim(0,max_TPR) +
    xlab('False Discovery Proportion') +
    ylab('True Positive Rate') +
    guides(colour = 'none') + theme_bw() +
    theme(strip.text = element_text(hjust = 0))
  p = cowplot::plot_grid(p2,p1,rel_heights = c(5,6), labels = c('A','B'),ncol=1)
  ## p = cowplot::plot_grid(p2,p1,rel_widths = c(11,9), labels = c('A','B'),ncol=2)
  p
}

p = make_intro_plot(results);p
ggsave('introduction_ABCD_3.png', p, width = 8, height = 6);p
p = make_intro_plot(results_filtered, max_TPR=1)
ggsave('introduction_ABCD_common.png', p, width = 10.5, height = 7.5);p

#### all
##########
max_TPR = .8
max_FDP = .1
pd = filter(results
          , method %in% c('msqrob','perseus' ,'msstats','dep_miximp','proteus')) %>%
    mutate(method = recode_factor(method
                                , 'msqrob' = 'MSqRob'
                                , 'perseus' = 'Perseus\n(maxLFQ)'
                                , 'proteus' = 'proteus\n(high-flyers)'
                                , 'dep_miximp' = 'DEPmiximp\n(maxLFQ)'
                                , 'msstats' = 'MSstats\n(median polish)'
                                  ))

p1 = ungroup(pd) %>%
  ggplot() + geom_hline(yintercept = 0) +
  geom_hline(aes(yintercept = logFC_real),colour = 'grey',size = 1.5) +
  geom_boxplot(aes(method,logFC,colour = method), outlier.size = 1) +
  facet_grid(sample1 + species~ sample2) +
  xlab('Method') +
  ylab('Log Fold Change') +
  ylim(-2,3) +
 guides(colour = FALSE)
p1
ggsave('supl_introduction_logfc_all.png', p1, width = 12, height = 15)

p2 = pd %>%
  ggplot(aes(FDP,TPR,colour = method))+ geom_path() +
  xlim(0,max_FDP) +
  ylim(0,.85) +
  facet_grid(sample1 ~ sample2) +
  xlab('False Discovery Proportion') +
  ylab('True Positive Rate')

p2
ggsave('supl_introduction_FDR_FDP_all.png', p2, width = 12, height = 12)

max_TPR = .8
max_FDP = .1
pd = filter(results
          , method %in% c('msqrob','perseus' ,'msstats','dep_miximp','proteus','msqrobsum')) %>%
    mutate(method = recode_factor(method
                                , 'msqrob' = 'MSqRob'
                                , 'msqrobsum' = 'MSqRobSum\nrobust regression'
                                , 'perseus' = 'Perseus\n(maxLFQ)'
                                , 'proteus' = 'proteus\n(high-flyers)'
                                , 'dep_miximp' = 'DEPmiximp\n(maxLFQ)'
                                , 'msstats' = 'MSstats\n(median polish)'
                                  ))

p1 = ungroup(pd) %>%
  ggplot() + geom_hline(yintercept = 0) +
  geom_hline(aes(yintercept = logFC_real),colour = 'grey',size = 1.5) +
  geom_boxplot(aes(method,logFC,colour = method), outlier.size = 1) +
  facet_grid(sample1 + species~ sample2) +
  xlab('Method') +
  ylab('Log Fold Change') +
  ylim(-2,3) +
  guides(colour = FALSE) + theme(axis.text.x = element_text(angle = 45,hjust = 1))
p1
ggsave('supl_introduction_logfc_all_msqrobsum.png', p1, width = 12, height = 15)

#################################
## msqrob vs msqrobSum comparison
#################################
make_msqrob_plot = function(results,max_TPR = .7){
  max_FDP = .1
  pd = filter(results
            , method %in% c('msqrob', 'msqrobsum')) %>%
    mutate(method = recode_factor(method
                                , 'msqrob' = 'MSqRob', 'msqrobsum' = 'MSqRobSum'
                                  )) %>%
    mutate(contrast = factor(contrast,
                             levels = c("b-a", "c-a", "d-a", "e-a", "c-b", "d-b", "d-c", "e-b", "e-c", "e-d")))

  pd = group_by(pd,contrast,FC_real) %>% summarise_at(vars(conditions,concentrations),~{list(first(.x))}) %>%
    ungroup %>% transmute(contrast,contrast_info = factor(
                                     str_glue('{contrast}
Concentration (% wt/wt): {map_chr(conditions,~first(.x))} = {map_dbl(concentrations,~{first(.x) %>% round(2)})}%, {map_chr(conditions,~last(.x))} = {map_dbl(concentrations,~{last(.x) %>% round(2)})}%
Fold Change ({map_chr(conditions,~first(.x))}/{map_chr(conditions,~last(.x))}) = {round(FC_real,2)}'))) %>%
left_join(pd,.)
  pd$contrast_info = factor(pd$contrast_info, levels = levels(pd$contrast_info)[c(1,2,4,3,5,6,7:10)])

  ## q = filter(pd,contrast %in% c('b-a','c-a','d-a','c-b','d-b','d-c')) %>%
  q = pd %>%
    filter(contrast %in% c('b-a','c-b','d-c')) %>%
      mutate('1%' = .01,'5%' = .05) %>% gather(fdr, fdr_value, '1%', '5%') %>%
      group_by(method, contrast,fdr) %>%
      filter(qvalue == max(qvalue[qvalue < fdr_value])) %>% arrange(pvalue) %>%
      dplyr::slice(1)

p = pd %>%
  filter(contrast %in% c('b-a','c-b','d-c')) %>%
    ggplot(aes(FDP,TPR,colour = method))+ geom_path() +
    geom_vline(xintercept = .05, colour = 'grey') +
    geom_vline(xintercept = .01, colour = 'grey') +
    xlim(0,max_FDP) +
    ylim(0,max_TPR) +
      geom_point(aes(FDP,TPR,colour = method,shape = fdr),data = q,width = .002,fill = 'white',
                 height = .002,size = 2)+
    scale_shape_manual(name = 'FDR',values = c('1%' = 21, '5%' = 24)) +
  facet_wrap(~contrast_info,ncol = 3) +theme_bw() +
      theme(strip.text = element_text(hjust = 0)) +
    xlab('False Discovery Proportion') +
    ylab('True Positive Rate')
  p +theme(legend.position = 'bottom') +
    guides(colour = guide_legend(nrow = 1,title.position = 'top'),
           shape = guide_legend(nrow = 1,title.position = 'top'))
}
p = make_msqrob_plot(results);p
## ggsave('msqrob_ABCD.png', p, width = 10.5, height = 5);p
ggsave('msqrob_ABCD_3.png', p, width = 8, height = 4);p
##
## ALL
######
max_FDP = .1
pd = filter(results
          , method %in% c('msqrob', 'msqrobsum')) %>%
    mutate(method = recode_factor(method
                                , 'msqrob' = 'MSqRob', 'msqrobsum' = 'MSqRobSum'
                                  ))
q = filter(pd,contrast %in% c('b-a','c-a','d-a','c-b','d-b','d-c')) %>%
    mutate('1%' = .01,'5%' = .05) %>% gather(fdr, fdr_value, '1%', '5%') %>%
    group_by(method, contrast,fdr) %>%
    filter(qvalue == max(qvalue[qvalue < fdr_value])) %>% arrange(pvalue) %>%
    dplyr::slice(1)

p = pd %>%
  ggplot(aes(FDP,TPR,colour = method))+ geom_path() +
  geom_vline(xintercept = .05, colour = 'grey') +
  geom_vline(xintercept = .01, colour = 'grey') +
  xlim(0,max_FDP) +
  ylim(0,max_TPR) +
    geom_point(aes(FDP,TPR,colour = method,shape = fdr),data = q,width = .002,fill = 'white',
               height = .002,size = 2)+
    scale_shape_manual(name = 'FDR',values = c('1%' = 21, '5%' = 24)) +
  facet_grid(sample1 ~ sample2) +
  xlab('False Discovery Proportion') +
  ylab('True Positive Rate')
p

ggsave('supl_msqrob_all.png', p2, width = 12, height = 12)

#####################
## improvement plots
#####################
make_improvment_plot = function(results,max_TPR = .65){
    max_FDP = .1

    pd = filter(results
              , method %in% c(
                                'msqrobsum'
                              , 'perseus'
                              , 'perseus_my'
                              , 'msqrobsum_maxlfq'
                              , 'msqrobsum_DEP'
                            )) %>%
        mutate(method = recode_factor(method
                                    , 'perseus' = 'perseus default \nmaxLFQ + t-tests\n'
                                    , 'perseus_my' = 'perseus vsn\nmaxLFQ + vsn + t-tests\n'
                                    , 'msqrobsum_maxlfq' = 'MSqRobSum maxLFQ \nmaxLFQ + vsn + MSqRob\n'
                                    , 'msqrobsum_DEP' = 'MSqRobSum DEP\nmaxLFQ + vsn + imp. + MSqRob\n'
                                    , 'msqrobsum' = 'MSqRobSum default\nvsn + rob. sum.  + MSqRob\n'
                                      )) %>%
        mutate(contrast = factor(contrast,
                                 levels = c("b-a", "c-a", "d-a", "e-a", "c-b", "d-b", "d-c", "e-b", "e-c", "e-d")))

    pd = group_by(pd,contrast,FC_real) %>% summarise_at(vars(conditions,concentrations),~{list(first(.x))}) %>%
      ungroup %>% transmute(contrast,contrast_info = factor(
                                       str_glue('{contrast}
Concentration (% wt/wt): {map_chr(conditions,~first(.x))} = {map_dbl(concentrations,~{first(.x) %>% round(2)})}%, {map_chr(conditions,~last(.x))} = {map_dbl(concentrations,~{last(.x) %>% round(2)})}%
Fold Change ({map_chr(conditions,~first(.x))}/{map_chr(conditions,~last(.x))}) = {round(FC_real,2)}'))) %>%
left_join(pd,.)
    pd$contrast_info = factor(pd$contrast_info, levels = levels(pd$contrast_info)[c(1,2,4,3,5,6,7:10)])

    ## q = filter(pd,contrast %in% c('b-a','c-a','d-a','c-b','d-b','d-c')) %>%
      q = pd %>%
      filter(contrast %in% c('b-a','c-b','d-c')) %>%
        mutate('1%' = .01,'5%' = .05) %>% gather(fdr, fdr_value, '1%', '5%') %>%
        group_by(method, contrast,fdr) %>%
        filter(qvalue == max(qvalue[qvalue < fdr_value])) %>% arrange(pvalue) %>%
        dplyr::slice(1)

    p = pd %>%
      filter(contrast %in% c('b-a','c-b','d-c')) %>%
        ggplot(aes(FDP,TPR,colour = method))+ geom_path() +
        geom_vline(xintercept = .05, colour = 'grey') +
        geom_vline(xintercept = .01, colour = 'grey') +
        xlim(0,max_FDP) +
        ylim(0,max_TPR) +
        geom_point(aes(FDP,TPR,colour = method,shape = fdr),data = q,width = .002,fill = 'white',
                   height = .002,size = 2)+
        scale_shape_manual(name = 'FDR',values = c('1%' = 21, '5%' = 24)) +
      facet_wrap(~contrast_info,ncol = 3) + theme_bw() +
        theme(strip.text = element_text(hjust = 0)) +
        xlab('False Discovery Proportion') +
        ylab('True Positive Rate')
    p + theme(legend.position = 'bottom') +
      guides(colour = guide_legend(nrow = 2,title.position = 'top'),
        shape = guide_legend(nrow = 2,title.position = 'top'))
}

p = make_improvment_plot(results)
ggsave('improvement_ABCD_3.png', p, width = 8, height = 4.5);p

p = make_improvment_plot(results_filtered,max_TPR = 1)
ggsave('improvement_ABCD_common.png', p, width = 10.5, height = 5);p

###############################################
## Suplementary figure, failing summarisations
###############################################
## Select good example proteins where high flyers and maxlfq fail
tr = .01
aa = filter(res, method %in% c('msqrobsum', 'proteus', 'msqrobsum_maxlfq')
          , contrast %in% c('b-a', 'c-b', 'd-c', 'e-d')) %>%
    mutate(method = recode_factor(method
                                , 'proteus' = 'proteus'
                                , 'msqrobsum' = 'rr'
                                , 'msqrobsum_maxlfq' = 'maxLFQ'
                                  )) %>%
    mutate(contrast = factor(contrast, levels = c("b-a", "c-b", "d-c", "e-d"))) %>%
    select(protein,contrast,method,qvalue,human) %>% group_by(protein) %>%
    filter(n() == 12) %>% spread(method, qvalue) %>% select(- contrast) %>%
    group_by(protein,human) %>%
    summarise(rr_min = min(rr), proteus_min = min(proteus), maxLFQ_min = min(maxLFQ)
            , rr_max = max(rr), proteus_max = max(proteus), maxLFQ_max = max(maxLFQ)
            , rr_mean = mean(rr), proteus_mean = mean(proteus), maxLFQ_mean = mean(maxLFQ)
            , rr_median = median(rr), proteus_median = median(proteus), maxLFQ_median = median(maxLFQ)
            , rr = sum(rr < tr), proteus = sum(proteus < tr), maxLFQ = sum(maxLFQ < tr)) %>%
    ungroup

aaa = filter(aa,!human, rr==4,proteus == 0, maxLFQ == 0) %>% arrange(-proteus_mean);id1 = aaa$protein[1];id1
aaa = filter(aa,human, rr==0) %>% arrange(proteus_median);id2 = aaa$protein[1];id2
aaa = filter(aa,human, rr==0) %>% arrange(maxLFQ_median);id3 = aaa$protein[1];id3
ids = c(id1, id2, id3)
######################
## get intensity data
df1= read_rds('../analyses/output/msqrob') %>% select(protein, data,human) %>%
    filter(protein %in% ids) %>% mutate(sum = 'no summarization') %>% unnest
df2 = read_rds('../analyses/proteus/proteus_summarized_data.rds')$tab %>% 
        as_data_frame(rownames = 'protein') %>%
        gather(sample,expression,-protein) %>%
        mutate(condition = str_extract(sample, '^.')) %>%
        filter(protein  %in% ids) %>% mutate(sum = 'high flyer summarization')
df3 = read_rds('../analyses/output/msqrobsum_maxlfq') %>% select(protein, data) %>%
    filter(protein %in% ids) %>% mutate(sum = 'maxLFQ summarization') %>% unnest %>% select(-feature)
df4 = read_rds('../analyses/output/msqrobsum') %>% select(protein, data_summarized) %>%
    filter(protein %in% ids) %>% mutate(sum = 'robust regression summarization') %>% unnest

species = ifelse(df4$human[1],'human', 'E. coli')

df = bind_rows(df1, df2, df3,df4) %>% group_by(protein) %>%
    mutate(sum = recode_factor(as.factor(sum)
                              ,'no summarization' = '1\tno summarization'
                             , 'high flyer summarization' = '2\thigh flyers summarization'
                             , 'maxLFQ summarization' = '3\tmaxLFQ summarization'
                             , 'robust regression summarization' =  '4\trobust regression summarization'
                               )
         , title = str_glue("{ifelse(protein == id1,'A\t',ifelse(protein == id2,'B\t','C\t'))}{ifelse(any(human,na.rm = TRUE),'human', 'E. coli')} protein: {protein}") ) %>%
    mutate(feature = map_chr(feature,~{paste(rev(str_split(.x,'')[[1]]),collapse = '')}))

  p = ggplot(df) + geom_boxplot(aes(sample,expression,group = condition),outlier.size = -1) +
    geom_point(aes(sample,expression,colour = as.factor(feature),group = condition)) + guides(colour = 'none') +
    theme(strip.text = element_text(hjust = 0)) +
      facet_grid(title~sum,scales = 'free_y') +
      theme(axis.text.x = element_text(angle = 55, hjust = 1)
          , strip.text = element_text(size = 11)) + ylab('log2 intensity')
  p

ggsave(str_glue('protein_comparison_boxplots.png'),p,width=12,height=8)

###################################################
## improvement plots rmm DEP miximp vs no miximp ##
###################################################
make_dep_plot = function(results,max_TPR = .8){
  max_FDP = .1

pd = filter(results
          , method %in% c(
                          'msqrobsum'
                        , 'msqrobsum_miximp_pep'
                        , 'msqrobsum_miximp_prot'
                        )) %>%
  mutate(method = recode_factor(method
                              , 'msqrobsum' = 'MSqRobSum'
                              , 'msqrobsum_miximp_pep' = 'MSqRobSum with\nmix. imp. peptide level'
                              , 'msqrobsum_miximp_prot' = 'MSqRobSum with\nmix. imp. protein level'
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

  q.01 = group_by(pd, method, contrast) %>%
    filter(contrast %in% c('b-a','c-a','d-a','c-b','d-b','d-c')) %>%
    filter(qvalue == max(qvalue[qvalue < .01])) %>% arrange(pvalue) %>%
    dplyr::slice(1)

  q.05 = group_by(pd, method, contrast) %>%
    filter(contrast %in% c('b-a','c-a','d-a','c-b','d-b','d-c')) %>%
    filter(qvalue == max(qvalue[qvalue < .05])) %>% arrange(pvalue) %>%
    dplyr::slice(1)

  p = pd %>%
    filter(contrast %in% c('b-a','c-a','d-a','c-b','d-b','d-c')) %>%
    ggplot(aes(FDP,TPR,colour = method))+ geom_path() +
    geom_vline(xintercept = .05, colour = 'grey') +
    geom_vline(xintercept = .01, colour = 'grey') +
    xlim(0,max_FDP) +
    ylim(0,max_TPR) +
    geom_point(aes(FDP,TPR,colour = method),data = q.05,width = .002, shape = 24,fill = 'white',
               height = .002,size = 2)+
    geom_point(aes(FDP,TPR,colour = method),data = q.01,width = .002, shape = 21,fill = 'white',
               height = .002,size = 2) +
    facet_wrap(~contrast_info,ncol = 3) +
    theme(strip.text = element_text(hjust = 0)) +
    xlab('False Discovery Proportion') +
    ylab('True Positive Rate')
  p
}

p = make_dep_plot(results)
ggsave('dep_improvement_ABCD.png', p, width = 10, height = 5);p

res2 = filter(res, method %in% c('msqrobsum', 'msqrobsum_miximp_pep', 'msqrobsum_miximp_prot'))
results2_filtered = group_by(res2,protein, contrast) %>% filter(n() == 3) %>%
    group_by(method, contrast) %>% mutate(qvalue = p.adjust(pvalue, method = 'BH')) %>% ungroup
results2_filtered = preprocess_results(results2_filtered)

p = make_dep_plot(results_filtered,max_TPR = 1)
ggsave('dep_improvement_ABCD_common.png', p, width = 10, height = 5);p

########################
## overspiking effect ##
########################
pd = filter(results
          , method %in% c('msqrob'
                         ,'msqrobsum')) %>%
    mutate(method = recode_factor(method
                                , 'msqrob' = 'MSqRob'
                                , 'msqrobsum' = 'MSqRobSum'
                                  ))

q.05 = group_by(pd, method, contrast) %>%
  filter(qvalue == max(qvalue[qvalue < .05])) %>% arrange(pvalue) %>%
  dplyr::slice(1) %>% mutate(FDR = '5%', treshold = .05)
q.01 = group_by(pd, method, contrast) %>%
  filter(qvalue == max(qvalue[qvalue < .01])) %>% arrange(pvalue) %>%
  dplyr::slice(1) %>% mutate(FDR = '1%', treshold = .01)
q = bind_rows(q.05,q.01)
p1 = ggplot(q) +
  geom_hline(aes(yintercept = treshold), colour = 'grey') +
  geom_hline(yintercept = .0, colour = 'black') +
  geom_label(aes(conc_diff,FDP,label = contrast), size = 2.5,label.size=.2) +
  xlab('Difference in E.coli spike-in concentration (wt/wt%)') +
  ylab('False Discovery Proportion') +
  scale_x_continuous(breaks = unique(pd$conc_diff)) +
    facet_grid(method~FDR)
p1

p2 = filter(pd,!ecoli) %>%
  group_by(method,contrast,conc_diff) %>% summarise(logFC = mean(logFC)) %>%
  ggplot +
  geom_hline(yintercept = .0, colour = 'black') +
  geom_label(aes(conc_diff,logFC,label = contrast), size = 2.5,label.size=.2) +
  xlab('Difference in E.coli spike-in concentration (wt/wt%)') +
  ylab('Log Fold Change in human proteins') +
  scale_x_continuous(breaks = unique(pd$conc_diff)) +
  facet_grid(method~.)
## facet_wrap(~method)
p2

p = cowplot::plot_grid(p2,p1, labels = c('A','B'),ncol=2,rel_widths = c(4,6));p
ggsave('overspiking_rmm_vs_rmm_rr.png', p, width = 10, height = 8);p

###############################################
## DEP performance drop due to fdrtool vs BH ##
###############################################
## DEP pvalues are limma pvalues so ordening remains the same, also the performance curve

make_msqrob_plot = function(results,max_TPR = .8){
  max_FDP = .1
  pd = filter(results
            , method %in% c('rmm_lib','DEP_miximp','DEP_miximp_BH')) %>%
      mutate(method = recode_factor(method
                                  , 'rmm_lib' = 'MSqRob\n(BH)'
                                  , 'DEP_miximp' = 'DEPmiximp\n(emp. FDR)'
                                  , 'DEP_miximp_BH' = 'DEPmiximp\n(BH)'
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


  q = filter(pd,contrast %in% c('b-a','c-a','d-a','c-b','d-b' ,'d-c')) %>%
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


ggsave('DEP_BH_ABCD.png', p, width = 10.5, height = 7.5);p

###################################################
## msqrob vs msqrobSum comparison without Q9BZJ0 ##
###################################################
make_msqrob_plot = function(results,max_TPR = .8){
  max_FDP = .1
  pd = filter(results
            , method %in% c('msqrob', 'msqrobsum'
                          , 'msqrob_wo_Q9BZJ0', 'msqrobsum_wo_Q9BZJ0')) %>%
    mutate(method = recode_factor(method
                                , 'msqrob' = 'MSqRob', 'msqrobsum' = 'MSqRobSum'
                                , 'msqrob_wo_Q9BZJ0' = 'MSqRob\nwithout Q9BZJ0'
                                , 'msqrobsum_wo_Q9BZJ0' = 'MSqRobSum\nwithout Q9BZJ0'
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
ggsave('msqrob_ABCD_wo_Q9BZJ0.png', p, width = 10.5, height = 5);p

#########################################
## Plot intensties from protein Q9BZJ0 ##
#########################################
## get intensity data
ids = 'Q9BZJ0'
df1= read_rds('../analyses/output/msqrob') %>% select(protein, data,human) %>%
    filter(protein %in% ids) %>% mutate(sum = 'no summarization') %>% unnest
df2 = read_rds('../analyses/proteus/proteus_summarized_data.rds')$tab %>% 
        as_data_frame(rownames = 'protein') %>%
        gather(sample,expression,-protein) %>%
        mutate(condition = str_extract(sample, '^.')) %>%
        filter(protein %in% ids) %>% mutate(sum = 'high flyer summarization')
df3 = read_rds('../analyses/output/msqrobsum_maxlfq') %>% select(protein, data) %>%
    filter(protein %in% ids) %>% mutate(sum = 'maxLFQ summarization') %>% unnest %>% select(-feature)
df4 = read_rds('../analyses/output/msqrobsum') %>% select(protein, data_summarized) %>%
    filter(protein %in% ids) %>% mutate(sum = 'robust regression summarization') %>% unnest

species = ifelse(df4$human[1],'human', 'E. coli')

df = bind_rows(df1, df2, df3,df4) %>% group_by(protein) %>%
    mutate(sum = recode_factor(as.factor(sum)
                              ,'no summarization' = 'no summarization'
                             , 'high flyer summarization' = 'high flyers summarization'
                             , 'maxLFQ summarization' = 'maxLFQ summarization'
                             , 'robust regression summarization' =  'robust regression summarization'
                               )
         , title = str_glue("{ifelse(protein == id1,'A\t',ifelse(protein == id2,'B\t','C\t'))}{ifelse(any(human,na.rm = TRUE),'human', 'E. coli')} protein: {protein}") ) %>%
    mutate(feature = map_chr(feature,~{paste(rev(str_split(.x,'')[[1]]),collapse = '')}))

  p = ggplot(df) + geom_boxplot(aes(sample,expression,group = condition),outlier.size = -1) +
    geom_point(aes(sample,expression,colour = as.factor(feature),group = condition)) + guides(colour = 'none') +
    theme(strip.text = element_text(hjust = 0)) +
      facet_wrap(~sum,2) +
      theme(axis.text.x = element_text(angle = 55, hjust = 1)
          , strip.text = element_text(size = 11)) + ylab('log2 intensity')
  p

ggsave(str_glue('protein_Q9BZJ0_boxplots.png'),p,width=10,height=7)

###############################################################
## Effect of vsn on normalization and batch parameter effect ##
###############################################################

#################################
## msqrob vs msqrobSum vs pairwise comparison
#################################
make_msqrob_plot = function(results,max_TPR = .8){
  max_FDP = .1
  pd = filter(results
            ## , method %in% c('msqrob', 'msqrobsum','msqrob_pairwise', 'msqrobsum_pairwise')) %>%
            , method %in% c('msqrob', 'msqrobsum','msqrob_pairwise2',
                            'msqrobsum_pairwise2','dep_miximp','dep_miximp_pairwise')) %>%
    mutate(method = recode_factor(method
                                , 'msqrob' = 'MSqRob', 'msqrobsum' = 'MSqRobSum'
                                , 'msqrob_pairwise2' = 'MSqRob pairwise', 'msqrobsum_pairwise2' = 'MSqRobSum pairwise'
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
ggsave('msqrob_vs_pairwise_ABCD.png', p, width = 10.5, height = 5);p

#################################
## msqrob vs msqrobSum vs pairwise comparison
#################################
make_msqrob_plot = function(results,max_TPR = .8){
  max_FDP = .1
  pd = filter(results
            ## , method %in% c('msqrob', 'msqrobsum','msqrob_pairwise', 'msqrobsum_pairwise')) %>%
            , method %in% c(
                            'msqrob', 'msqrobsum',
                            'msqrob_pairwise', 'msqrobsum_pairwise',
                            'msqrob_pairwise2', 'msqrobsum_pairwise2'
                            ## 'msqrobsum', 'msqrobsumnew','msqrobsum_pairwise2'
                            ## 'dep_miximp','dep_miximp_pairwise'
                          )) %>%
    mutate(method = recode_factor(method
                                ## , 'msqrob' = 'MSqRob', 'msqrobsum' = 'MSqRobSum'
                                , 'msqrob_pairwise' = 'MSqRob pairwise', 'msqrobsum_pairwise' = 'MSqRobSum pairwise'
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


a = filter(res, method %in% c('msqrobsum', 'msqrobsum_pairwise2'), contrast == 'b-a')

count(a,method)
aa = select(a, protein,method,pvalue) %>% spread(method,pvalue)
plot(aa$msqrobsum,aa$msqrobsum_pairwise2)
group_by(a,protein)

head(res)
filter(res, contrast == 'b-a', !is.na(qvalue)) %>% count( method)



####################################
## msqrob number of proteins plot ##
####################################

count(res,method)

r = filter(results, method %in% c('msqrob', 'msqrobsum','msqrob_pairwise2', 'msqrobsum_pairwise2'),
           contrast == 'b-a')
r = group_by(r, method) %>% mutate(n_proteins = row_number(qvalue))

r = arrange(r,qvalue)
a = ggplot(r) + geom_line(aes(qvalue,n_proteins,colour = method)) + xlim(0, .1) + ylim(0,500)
b = ggplot(r) + geom_line(aes(qvalue,TPR,colour = method)) + xlim(0, .1) + ylim(0,1)
cowplot::plot_grid(a,b)

r = arrange(r,TPR)
a = ggplot(r) + geom_path(aes(FDP,n_proteins,colour = method)) + xlim(0, .1) + ylim(0,500)
b = ggplot(r) + geom_path(aes(FDP,TPR,colour = method)) + xlim(0, .1) + ylim(0,1)
cowplot::plot_grid(a,b)


###########################################################################
## improvement plots msstats with or without normalization or imputation ##
###########################################################################
make_dep_plot = function(results,max_TPR = .8){
  max_FDP = .1

pd = filter(results
          , method %in% c(
                          'msqrobsum'
                        , 'msstats'
                        , 'msstats_noImp'
                        , 'msstats_noNorm'
                        , 'msstats_noImp_noNorm'
                        )) %>%
  mutate(method = recode_factor(method
                              , 'msqrobsum' = 'MSqRobSum'
                              , 'msstats' = 'MSstats'
                              , 'msstats_noImp' = 'MSstats\nwithout imputation'
                              , 'msstats_noNorm' = 'MSstats\nwithout normalization'
                              , 'msstats_noImp_noNorm' = 'MSstats\nwithout imputation and normalization'
                              ## , 'msqrobsum_miximp_pep' = 'MSqRobSum with\nmix. imp. peptide level'
                              ## , 'msqrobsum_miximp_prot' = 'MSqRobSum with\nmix. imp. protein level'
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

  ## q.01 = group_by(pd, method, contrast) %>%
  ##   filter(contrast %in% c('b-a','c-a','d-a','c-b','d-b','d-c')) %>%
  ##   filter(qvalue == max(qvalue[qvalue < .01])) %>% arrange(pvalue) %>%
  ##   dplyr::slice(1)

  ## q.05 = group_by(pd, method, contrast) %>%
  ##   filter(contrast %in% c('b-a','c-a','d-a','c-b','d-b','d-c')) %>%
  ##   filter(qvalue == max(qvalue[qvalue < .05])) %>% arrange(pvalue) %>%
  ##   dplyr::slice(1)

  q = pd %>%
    filter(contrast %in% c('b-a','c-b','d-c')) %>%
    mutate('1%' = .01,'5%' = .05) %>% gather(fdr, fdr_value, '1%', '5%') %>%
    group_by(method, contrast,fdr) %>%
    filter(qvalue == max(qvalue[qvalue < fdr_value])) %>% arrange(pvalue) %>%
    dplyr::slice(1)

  p = pd %>%
    filter(contrast %in% c('b-a','c-b','d-c')) %>%
    ggplot(aes(FDP,TPR,colour = method))+ geom_path() +
    geom_vline(xintercept = .05, colour = 'grey') +
    geom_vline(xintercept = .01, colour = 'grey') +
    xlim(0,max_FDP) +
    ylim(0,max_TPR) +
    geom_point(aes(FDP,TPR,colour = method,shape = fdr),data = q,width = .002,fill = 'white',
               height = .002,size = 2)+
    scale_shape_manual(name = 'FDR',values = c('1%' = 21, '5%' = 24)) +
    facet_wrap(~contrast_info,ncol = 3) +theme_bw() +
    theme(strip.text = element_text(hjust = 0)) +
    xlab('False Discovery Proportion') +
    ylab('True Positive Rate')
  p +theme(legend.position = 'bottom') +
    guides(colour = guide_legend(nrow = 1,title.position = 'top'),
           shape = guide_legend(nrow = 1,title.position = 'top'))
  ## p = pd %>%
  ##   filter(contrast %in% c('b-a','c-a','d-a','c-b','d-b','d-c')) %>%
  ##   ggplot(aes(FDP,TPR,colour = method))+ geom_path() +
  ##   geom_vline(xintercept = .05, colour = 'grey') +
  ##   geom_vline(xintercept = .01, colour = 'grey') +
  ##   xlim(0,max_FDP) +
  ##   ylim(0,max_TPR) +
  ##   geom_point(aes(FDP,TPR,colour = method),data = q.05,width = .002, shape = 24,fill = 'white',
  ##              height = .002,size = 2)+
  ##   geom_point(aes(FDP,TPR,colour = method),data = q.01,width = .002, shape = 21,fill = 'white',
  ##              height = .002,size = 2) +
  ##   facet_wrap(~contrast_info,ncol = 3) +
  ##   theme(strip.text = element_text(hjust = 0)) +
  ##   xlab('False Discovery Proportion') +
  ##   ylab('True Positive Rate')
  ## p
}

p = make_dep_plot(results,.65)
ggsave('msstats_norm_imp_improvement_ABCD.png', p, width = 10, height = 5);p


#### density plots to show effect of normalization on expresson
#######################################################
library(msqrobsum)
source("../analyses/functions.r")
mq_path = '../data/maxquant/'
fasta_path = '../data/fasta/'
a = read_mq(mq_path, fasta_path) %>% preprocess(norm = 'vsn',logtrans = FALSE)
b = read_mq(mq_path, fasta_path) %>% preprocess(norm = 'quantiles',logtrans = TRUE)
c = read_mq(mq_path, fasta_path) %>% preprocess(norm = FALSE, logtrans = TRUE)

d = list('vsn normalization' = a, 'no normalization' = c, 'quantile normalization' = b)

dt = d %>% imap_dfr(~{limma::plotMDS(exprs(.x), top = Inf)[[3]] %>% cbind(pData(.x)) %>%
               mutate(method = .y,sample = sampleNames(.x))})
pmds = ggplot(dt) + geom_label(aes(`1`,`2`, label = sample,colour = condition)) +
  facet_grid(method~.) + guides(colour = FALSE) +
  xlab('Leading logFC dim 1') +
  ylab('Leading logFC dim 2')
pmds

dt = d %>% imap_dfr(~{MSnSet2df(.x) %>% mutate(method = .y)})
pdens = ggplot(dt) +  geom_density(aes(expression,group = sample,colour = condition),linetype = 3) +
  ## theme(legend.position = 'top') +
  facet_grid(method~.)
pdens
p = cowplot::plot_grid(pmds,pdens,ncol = 2,rel_widths = c(1,1.5))
ggsave('spikein_normalization_density_mds_plots.png', p, width = 10, height = 8 );p
