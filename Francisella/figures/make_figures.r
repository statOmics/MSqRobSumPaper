library(tidyverse)
##########
res = list.files('../analyses/results', full.names = TRUE) %>%
    map_dfr(~{print(.x);read_rds(.x) %>% mutate(method = str_replace(.x,'.+\\/([^\\/]+)','\\1'))})

res = filter(res, !is.na(qvalue)) %>% group_by(method) %>%
  arrange(pvalue) %>% mutate(n_proteins = row_number(pvalue))

########################################
## performance with to or z statistic ##
########################################
r = filter(res, method %in% c('msqrobsum_quant','msqrob_quant'))
r  = group_by(r, method) %>%
  mutate(pvalue = pnorm(-abs(t)) * 2) %>% arrange(pvalue) %>% 
  mutate(qvalue = p.adjust(pvalue, 'BH')
       , n_proteins = row_number(pvalue)) %>% 
  ungroup %>% mutate(method = str_glue('{method}_z')) %>%
  bind_rows(r)
r = filter(res,method %in% c('original')) %>% bind_rows(r)

p = ungroup(r) %>% mutate(method = recode_factor(method,'original' = 'Ramond et al.'
                                            , 'msqrob_quant' = 'MSqRob (t-test)'
                                            , 'msqrobsum_quant' = 'MSqRobSum (t-test)'
                                            , 'msqrob_quant_z' = 'MSqRob (z-test)'
                                            , 'msqrobsum_quant_z' = 'MSqRobSum (z-test)'
                                              )) %>%
  ggplot + geom_line(aes(qvalue, n_proteins, colour = method)) +
  xlim(0, .1) + ylim(0,330) + xlab('False Discovery Rate') + ylab('Number of proteins returned') +theme_bw()
p
ggsave('Francisella_performance_z_vs_t_test.png', p, width = 5, height = 3.5);p


p = ungroup(r) %>% filter(method %in% c('original', 'msqrob_quant', 'msqrobsum_quant')) %>%
  mutate(method = recode_factor(method,'original' = 'Ramond et al.'
                              , 'msqrob_quant' = 'MSqRob'
                              , 'msqrobsum_quant' = 'MSqRobSum'
                                                 )) %>%
  ggplot + geom_line(aes(qvalue, n_proteins, colour = method)) +
  xlim(0, .1) + ylim(0,330) + xlab('False Discovery Rate') + ylab('Number of proteins returned') +theme_bw()
p
ggsave('Francisella_performance.png', p, width = 5, height = 3.5);p
postscript('Francisella_performance.eps', width = 5, height = 3.5);p;dev.off()


## make plots to show why the proteins that are found by orignal are bad
########################################################################
pd = filter(res, method %in% c('msqrob_quant', 'msqrobsum_quant', 'original')
            , qvalue < .01) %>%
  mutate(logFC = ifelse(method == 'original', -logFC, logFC),
         upregulated = logFC > 0)
ungroup(pd) %>% mutate(method = recode_factor(method,'original' = 'Ramond et al.'
                                              , 'msqrob_quant' = 'MSqRob'
                                            , 'msqrobsum_quant' = 'MSqRobSum'
                                              )) %>%
  ggplot + geom_histogram(aes(logFC),bins = 30) + facet_grid(~method) + theme_bw()
ggsave('Francisella_foldchange_histogram.png', width = 6, height = 3)

#### density plots to show effect of normalization on expresson
#######################################################
library(msqrobsum)
source("../analyses/functions.r")
mq_path = '../data/maxquant/'
a = read_mq(mq_path) %>% preprocess(norm = 'vsn',logtrans = FALSE)
b = read_mq(mq_path) %>% preprocess(norm = 'quantiles',logtrans = TRUE)
c = read_mq(mq_path) %>% preprocess(norm = FALSE, logtrans = TRUE)

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
ggsave('normalization_density_mds_plots.png', p, width = 10, height = 8 );p
