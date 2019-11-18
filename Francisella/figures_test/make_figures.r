library(tidyverse)

##########
res = list.files('../analyses/results', full.names = TRUE) %>%
    map_dfr(~{print(.x);read_rds(.x) %>% mutate(method = str_replace(.x,'.+\\/([^\\/]+)','\\1'))})

res = filter(res, !is.na(qvalue)) %>% group_by(method) %>%
  arrange(pvalue) %>% mutate(n_proteins = row_number(pvalue))
## res = mutate(res,pvalue = ifelse(is.na(pvalue),0,pvalue)) %>% group_by(method) %>%
##   arrange(pvalue) %>% mutate(n_proteins = row_number(pvalue))

## results = preprocess_results(res)
## results_filtered = group_by(res,protein, contrast) %>% filter(n() == nlevels(results$method)) %>%
##   group_by(method, contrast) %>% mutate(qvalue = p.adjust(pvalue, method = 'BH')) %>% ungroup
## results_filtered = preprocess_results(results_filtered)

###################
group_by(res,method) %>% count(qvalue < .01)
count(res,method)

filter(res,method == 'original') %>% summary
filter(res,method == 'original',is.na(t)) %>% View

filter(res) %>% ggplot + geom_line(aes(qvalue, n_proteins, colour = method))
p = filter(res, method %in% c('msqrob','msqrob_quant','msqrobsum','msqrobsum_quant','original')) %>%
  ggplot + geom_line(aes(qvalue, n_proteins, colour = method)) +
  xlim(0, .1) + ylim(0,330) + xlab('q value') + ylab('Number of proteins')
p
## ggsave('performance.png', p, width = 5, height = 4);p

p = filter(res) %>% ggplot + geom_line(aes(qvalue, n_proteins, colour = method))
p = filter(res, method %in% c('msqrob','msqrob_noshrink')) %>%
  ggplot + geom_line(aes(qvalue, n_proteins, colour = method)) +
  xlim(0, .1) + ylim(0,330) + xlab('q value') + ylab('Number of proteins')
p
ggsave('performance_noshrink.png', p, width = 5, height = 4);p


filter(res, method %in% c('msqrob','original')) %>% select(method,protein,qvalue) %>% spread(method,qvalue) %>%
  ggplot + geom_point(aes(msqrob,original))
select(res)

filter_func = function(signif_method = 'msqrob', other_method ='original'){
  filter(res, method %in% c(signif_method,other_method)) %>%
    select(method,protein,qvalue) %>% spread(method,qvalue) %>%
    filter_at(signif_method, ~{. < .01}) %>% arrange_at(other_method,desc)
}

filter_func(c('msqrob','original'))
filter_func(c('original', 'msqrob'))

filter_func(c('msqrob','msqrobsum'))
filter_func(c('msqrobsum', 'msqrob'))

intensities = read_rds('../analyses/output/msqrob') %>% select(protein,data) %>% unnest
makeplot = function(signif_method = 'msqrob', other_method ='original'){
  id = filter_func(signif_method,other_method) %>% head(10)
  id = mutate(id, title =  str_glue('{protein}   {other_method} qvalue: {round(unlist(id[,other_method]),3)*100}%'))
  ## pd = left_join(id,intensities) %>% ggplot +
  ##   geom_point(aes(sample,expression,colour = condition)) + facet_wrap(~protein,2)
  pd = left_join(id,intensities) %>% mutate(biorep = str_glue('{condition}{str_extract(biorep,".$")}'))
  pd2 = pd %>% group_by(title,condition,feature,biorep) %>% summarise(expression = median(expression))
  ggplot(pd) + geom_point(aes(biorep,expression,colour = feature))  +
    geom_line(aes(biorep,expression,group = feature,colour = feature),pd2, linetype = 2,alpha = .5) + facet_wrap(~title,ncol = 2) + guides(colour=FALSE) + xlab('') +
    scale_x_discrete(breaks = unique(pd$biorep),labels = str_extract(unique(pd$biorep),'^..'))
}
p + scale_x_discrete(breaks = unique(pd$biorep),labels = str_extract(unique(pd$biorep),'^..'))


makeplot('msqrob','original') 
makeplot('original','msqrob')
makeplot('msqrobsum','original')
makeplot('original','msqrobsum')
makeplot('msqrob','msqrobsum')
makeplot('msqrobsum','msqrob')

saveplot =function(signif_method = 'msqrob', other_method ='original'){
  p = makeplot(signif_method,other_method)
  ggsave(str_glue('{signif_method}_significant_{other_method}_nonsignificant_FDR1.pdf'),p,width = 6,height = 8)
}

saveplot('msqrob','original') 
saveplot('original','msqrob')
saveplot('msqrobsum','original')
saveplot('original','msqrobsum')
saveplot('msqrob','msqrobsum')
saveplot('msqrobsum','msqrob')

count(res,method)

pd = filter(res, method %in% c('msqrob_quant', 'original_ludger')) %>% select(protein,method, qvalue) %>% spread(method,qvalue)

ggplot(pd) + geom_point(aes(msqrob_quant,original_ludger))
ggplot(pd) + geom_point(aes(msqrob_quant,original_ludger)) + xlim(0,.1) + ylim(0,.1)

dev.new()

#################################
## check pvalue on z statistic ##
#################################
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
ggsave('Francisella_performance_z_vs_t_test.png', p, width = 6, height = 4);p


p = ungroup(r) %>% filter(method %in% c('original', 'msqrob_quant', 'msqrobsum_quant')) %>%
  mutate(method = recode_factor(method,'original' = 'Ramond et al.'
                              , 'msqrob_quant' = 'MSqRob'
                              , 'msqrobsum_quant' = 'MSqRobSum'
                                                 )) %>%
  ggplot + geom_line(aes(qvalue, n_proteins, colour = method)) +
  xlim(0, .1) + ylim(0,330) + xlab('False Discovery Rate') + ylab('Number of proteins returned') +theme_bw()
p
ggsave('Francisella_performance.png', p, width = 6, height = 4);p

#### density plots to show why vsn is worse then msqrob
#######################################################
library(msqrobsum)
source("../analyses/functions.r")
mq_path = '../data/maxquant/'
a = set = read_mq(mq_path) %>% preprocess(norm = 'vsn',logtrans = FALSE) %>%
  MSnSet2df %>% mutate(normalization = 'VSN')
b = read_mq(mq_path) %>% preprocess(norm = 'quantiles',logtrans = TRUE) %>%
  MSnSet2df %>% mutate(normalization = 'quantile')
c = read_mq(mq_path) %>% preprocess(norm = FALSE, logtrans = TRUE) %>%
  MSnSet2df %>% mutate(normalization = 'none')

d = bind_rows(a,b,c)

p = ggplot(d) +  geom_density(aes(expression,group = sample,colour = condition)) + facet_grid(normalization~.)
ggsave('normalization_densityplots.png', p, width = 6, height = 4 );p


glimpse(res)
count(res, method)

a = filter(res, method == 'msqrob')
boxplot(a$logFC)
abline(h = 0)
b = filter(res, method == 'msqrob',qvalue < .01)
boxplot(b$logFC)
abline(h = 0)

par(mfrow = c(1,3))
a = filter(res, method == 'msqrob_quant')
boxplot(a$logFC)
abline(h = 0)
b = filter(res, method == 'msqrob_quant',qvalue < .01)
boxplot(b$logFC)
ggplot(b) + geom_violin(aes(contrast,logFC))
vi
abline(h = 0)
filter(b, protein == 'gi|118496676') %>% glimpse

b = filter(res, method == 'original',qvalue < .01)
ggplot(b) + geom_violin(aes(contrast,-logFC))
boxplot(b$logFC)
abline(h = 0)
filter(b, protein == 'gi|118496676') %>% glimpse

b = filter(res, method == 'original_ludger',qvalue < .01)
ggplot(b) + geom_violin(aes(contrast,logFC))
boxplot(b$logFC)
abline(h = 0)

b = filter(res, method == 'msqrobsum_quant',qvalue < .01)
ggplot(b) + geom_violin(aes(contrast,logFC))
glimpse(b)
filter(b, protein == 'gi|118496676') %>% glimpse

o = read_rds('/home/st/MSqRobSumPaper/Francisella/analyses/output/msqrobsum_quant')

p = filter(o, protein == 'gi|118496676')$data[[1]]

ggplot(p) + geom_point(aes(condition,expression))

plot(p$condition,p$expression)

%>% glimpse
filter(o)

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
## ggplot(pd) + geom_violin(aes(method, logFC),scale = 'count')
## ggsave('Francisella_foldchange_violin.png')
## ggplot(pd) + geom_boxplot(aes(upregulated,logFC),varwidth = TRUE) + facet_grid(~method)
## ggsave('Francisella_foldchange_boxplot.png')
## ggplot(pd) + geom_boxplot(aes(method,logFC, colour = upregulated),varwidth = TRUE)

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
