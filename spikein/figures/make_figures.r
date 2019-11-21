library(tidyverse)

#################################
## functions for preprocessing ##
#################################
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

##################
## main Figures ##
##################
## abc
######
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
  p = cowplot::plot_grid(p2,p1,rel_heights = c(3.5,3), labels = c('A','B'),ncol=1)
  ## p = cowplot::plot_grid(p2,p1,rel_widths = c(11,9), labels = c('A','B'),ncol=2)
  p
}
p = make_intro_plot(results);p
ggsave('spikein_introduction_ABC.png', p, width = 8, height = 6.5);p
postscript('spikein_introduction_ABC.eps', width = 8, height = 6.5);p;dev.off()
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
ggsave('spikein_msqrob_ABC.png', p, width = 8, height = 4);p
postscript('spikein_msqrob_ABC.eps', width = 8, height = 4);p;dev.off()

##
## ALL
######
## max_FDP = .1
## pd = filter(results
##           , method %in% c('msqrob', 'msqrobsum')) %>%
##     mutate(method = recode_factor(method
##                                 , 'msqrob' = 'MSqRob', 'msqrobsum' = 'MSqRobSum'
##                                   ))
## q = filter(pd,contrast %in% c('b-a','c-a','d-a','c-b','d-b','d-c')) %>%
##     mutate('1%' = .01,'5%' = .05) %>% gather(fdr, fdr_value, '1%', '5%') %>%
##     group_by(method, contrast,fdr) %>%
##     filter(qvalue == max(qvalue[qvalue < fdr_value])) %>% arrange(pvalue) %>%
##     dplyr::slice(1)

## p = pd %>%
##   ggplot(aes(FDP,TPR,colour = method))+ geom_path() +
##   geom_vline(xintercept = .05, colour = 'grey') +
##   geom_vline(xintercept = .01, colour = 'grey') +
##   xlim(0,max_FDP) +
##   ylim(0,max_TPR) +
##     geom_point(aes(FDP,TPR,colour = method,shape = fdr),data = q,width = .002,fill = 'white',
##                height = .002,size = 2)+
##     scale_shape_manual(name = 'FDR',values = c('1%' = 21, '5%' = 24)) +
##   facet_grid(sample1 ~ sample2) +
##   xlab('False Discovery Proportion') +
##   ylab('True Positive Rate')
## p

## ggsave('supl_msqrob_all.png', p2, width = 12, height = 12)

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
ggsave('spikein_improvement_ABC.png', p, width = 8, height = 4.5);p
postscript('spikein_improvement_ABC.eps', width = 8, height = 4.5);p;dev.off()

#########################
## Suplementary Figures ##
#########################
###############################################
## failing summarisations
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

## select proteins to plot
aaa = filter(aa,!human, rr==4,proteus == 0, maxLFQ == 0) %>% arrange(-proteus_mean);id1 = aaa$protein[1];id1
aaa = filter(aa,human, rr==0) %>% arrange(proteus_median);id2 = aaa$protein[1];id2
aaa = filter(aa,human, rr==0) %>% arrange(maxLFQ_median);id3 = aaa$protein[1];id3
ids = c(id1, id2, id3)
## get intensity data of selected proteins
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
  p = p + theme_bw()
p
ggsave(str_glue('spikein_supl_protein_comparison_boxplots.png'),p,width=12,height=8)

###########################################################################
## improvement plots msstats with or without normalization or imputation ##
###########################################################################
make_msstats_plot = function(results,max_TPR = .8){
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
    xlab('False Discovery Proportion') +
    ylab('True Positive Rate')
  p + theme_bw() +theme(strip.text = element_text(hjust = 0),legend.position = 'bottom') +
    guides(colour = guide_legend(nrow = 1,title.position = 'top'),
           shape = guide_legend(nrow = 1,title.position = 'top'))
}

p = make_msstats_plot(results,.65)
ggsave('spikein_supl_msstats_norm_imp_improvement_ABC.png', p, width = 10, height = 5);p

#######################################
## abcd comparisons shared and unshared
#######################################

make_introsup_plot = function(results,max_TPR = .7){
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

  p1 = pd %>%
    filter(contrast %in% c('c-a','d-a','d-b')) %>%
    ggplot(aes(FDP,TPR,colour = method))+ geom_path() +
    facet_wrap(~contrast_info,ncol = 3) +
    xlim(0,max_FDP) +
    ylim(0,max_TPR) +
    xlab('False Discovery Proportion') +
    ylab('True Positive Rate') +
    guides(colour = 'none') + theme_bw() +
    theme(strip.text = element_text(hjust = 0))

  ## Filter out common proteins
  pd2 = group_by(pd,protein, contrast) %>% filter(n() == 5) %>%
    group_by(method, contrast) %>% mutate(qvalue = p.adjust(pvalue, method = 'BH')) %>%
    arrange(pvalue) %>%
    mutate(FDP = cummean(!TP)
         , FPR = cumsum(!TP)/sum(!TP)
         , TPR  = cumsum(TP)/TP_total
         , TP_hits = cumsum(TP)) %>% ungroup

  p2 = pd2 %>%
    filter(contrast %in% c('b-a','c-b','d-c','c-a','d-a','d-b')) %>%
    ggplot(aes(FDP,TPR,colour = method))+ geom_path() +
    facet_wrap(~contrast_info,ncol = 3) +
    xlim(0,max_FDP) +
    ylim(0,max_TPR) +
    xlab('False Discovery Proportion') +
    ylab('True Positive Rate') +
    guides(colour = 'none') + theme_bw() +
    theme(strip.text = element_text(hjust = 0),
          legend.position = 'bottom') +
    guides(colour = guide_legend(nrow = 1))


  p = cowplot::plot_grid(p1,p2,rel_heights = c(3.5,7.5), labels = c('A','B'),ncol=1)
  p
}

p = make_introsup_plot(results, max_TPR=.77);p
ggsave('spikein_supl_introduction_ABCD_common.png', p, width = 8.15, height = 8);p

####################
#### logFC of all ##
####################
max_TPR = .8
max_FDP = .1
pd = filter(results
          , method %in% c('msqrob','perseus' ,'msstats','dep_miximp','proteus')) %>%
    mutate(method = recode_factor(method
                                , 'msqrob' = 'MSqRob'
                                , 'msqrobsum' = 'MSqRobSum'
                                , 'perseus' = 'Perseus\n(maxLFQ)'
                                , 'proteus' = 'proteus\n(high-flyers)'
                                , 'dep_miximp' = 'DEPmiximp\n(maxLFQ)'
                                , 'msstats' = 'MSstats\n(median polish)'
                                  ))

p1 = ungroup(pd) %>%
  ggplot() + geom_hline(yintercept = 0) +
  geom_hline(aes(yintercept = logFC_real),colour = 'grey',size = 1.5) +
  geom_boxplot(aes(method,logFC,colour = method), outlier.size = 1) +
  ## facet_grid(sample1 + species~ sample2) +
  facet_wrap(species~ contrast) +
  xlab('Method') +
  ylab('Log Fold Change') +
  ylim(-2,3) +
 guides(colour = FALSE) + scale_x_discrete(breaks = NULL) + xlab(NULL)
p1 = p1 + theme_bw() + theme(legend.position = 'bottom') +
  guides(colour = guide_legend(nrow = 1))
p1
ggsave('spikein_supl_introduction_logfc_all.png', p1, width = 11, height = 10)

#################################
## msqrob vs msqrobSum comparison
#################################
make_msqrobsupl_plot = function(results,max_TPR = .7){
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

  q = pd %>%
    filter(contrast %in% c('c-a','d-a','d-b')) %>%
      mutate('1%' = .01,'5%' = .05) %>% gather(fdr, fdr_value, '1%', '5%') %>%
      group_by(method, contrast,fdr) %>%
      filter(qvalue == max(qvalue[qvalue < fdr_value])) %>% arrange(pvalue) %>%
      dplyr::slice(1)

p = pd %>%
  filter(contrast %in% c('c-a','d-a','d-b')) %>%
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
p = make_msqrobsupl_plot(results,max_TPR = .8);p
ggsave('spikein_supl_msqrob_ABC.png', p, width = 8, height = 4);p

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
                               ))

  p = ggplot(df) + geom_boxplot(aes(sample,expression,group = condition),outlier.size = -1) +
    geom_point(aes(sample,expression,colour = as.factor(feature),group = condition)) + guides(colour = 'none') +
    theme(strip.text = element_text(hjust = 0)) +
      facet_wrap(~sum,2) +
      theme(axis.text.x = element_text(angle = 55, hjust = 1)
          , strip.text = element_text(size = 11)) + ylab('log2 intensity') + theme_bw()
  p

ggsave(str_glue('spikein_supl_protein_Q9BZJ0_boxplots.png'),p,width=10,height=7)

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
Concentration (% wt/wt): {map_chr(conditions,~first(.x))} = {map_dbl(concentrations,~{first(.x) %>% round(2)})}%, {map_chr(conditions,~last(.x))} = {map_dbl(concentrations,~{last(.x) %>% round(2)})}%
Fold Change ({map_chr(conditions,~first(.x))}/{map_chr(conditions,~last(.x))}) = {round(FC_real,2)}'))) %>%
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
  p +theme_bw() + theme(strip.text = element_text(hjust = 0),legend.position = 'bottom') +
    guides(colour = guide_legend(nrow = 1,title.position = 'top'),
           shape = guide_legend(nrow = 1,title.position = 'top'))
}
p = make_msqrob_plot(results);p
ggsave('spikein_supl_msqrob_ABCD_wo_Q9BZJ0.png', p, width = 8.2, height = 6);p

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
    facet_grid(method~FDR) + theme_bw()
p1

p2 = filter(pd,!ecoli) %>%
  group_by(method,contrast,conc_diff) %>% summarise(logFC = mean(logFC)) %>%
  ggplot +
  geom_hline(yintercept = .0, colour = 'black') +
  geom_label(aes(conc_diff,logFC,label = contrast), size = 2.5,label.size=.2) +
  xlab('Difference in E.coli spike-in concentration (wt/wt%)') +
  ylab('Log Fold Change in human proteins') +
  scale_x_continuous(breaks = unique(pd$conc_diff)) +
  facet_grid(method~.) + theme_bw()
## facet_wrap(~method)
p2

p = cowplot::plot_grid(p2,p1, labels = c('A','B'),ncol=2,rel_widths = c(4,6));p
ggsave('spikein_supl_overspiking_rmm_vs_rmm_rr.png', p, width = 10, height = 8);p


###############################################################
## Effect of vsn on normalization and batch parameter effect ##
################V###############################################

#####################
## rest of improvement plots
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
        filter(contrast %in% c('c-a','d-a','d-b')) %>%
        mutate('1%' = .01,'5%' = .05) %>% gather(fdr, fdr_value, '1%', '5%') %>%
        group_by(method, contrast,fdr) %>%
        filter(qvalue == max(qvalue[qvalue < fdr_value])) %>% arrange(pvalue) %>%
        dplyr::slice(1)

    p = pd %>%
      filter(contrast %in% c('c-a','d-a','d-b')) %>%
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

p = make_improvment_plot(results,.75)
ggsave('spikein_sup_improvement_ABC.png', p, width = 8, height = 4.5);p

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
Concentration (% wt/wt): {map_chr(conditions,~first(.x))} = {map_dbl(concentrations,~{first(.x) %>% round(2)})}%, {map_chr(conditions,~last(.x))} = {map_dbl(concentrations,~{last(.x) %>% round(2)})}%
Fold Change ({map_chr(conditions,~first(.x))}/{map_chr(conditions,~last(.x))}) = {round(FC_real,2)}'))) %>%
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
  p +theme_bw() + theme(strip.text = element_text(hjust = 0),legend.position = 'bottom') +
    guides(colour = guide_legend(nrow = 1,title.position = 'top'),
           shape = guide_legend(nrow = 1,title.position = 'top'))
}

p = make_dep_plot(results,max_TPR = .9);p
ggsave('spikein_supl_dep_improvement_ABCD.png', p, width = 8.2, height = 6);p

