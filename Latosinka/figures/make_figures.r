library(tidyverse)
##########
res = list.files('../analyses/results', full.names = TRUE) %>%
    map_dfr(~{print(.x);read_rds(.x) %>% mutate(method = str_replace(.x,'.+\\/([^\\/]+)','\\1'))})

## remove NA pvalues
res = filter(res,!is.na(pvalue)) %>% group_by(method) %>%
  arrange(pvalue) %>% mutate(n_proteins = row_number(pvalue))

#######################################
## performance plot on same proteins ##
#######################################
dt = filter(res, method %in% c('msqrob', 'msqrobsum', 'dep_miximp', 'msstats'))
res2 = dt %>%
  group_by(protein) %>% filter(n() == 4) %>%
  group_by(method) %>% mutate(qvalue = p.adjust(pvalue, 'BH')) %>% mutate(n_proteins = row_number(pvalue))
pd = bind_rows(mutate(dt,dataset = 'all proteins')
             , mutate(res2,dataset = 'common proteins'))
pd = ungroup(pd) %>%
  mutate(method = recode_factor(method
                              , 'msqrob' = 'MSqRob'
                              , 'msqrobsum' = 'MSqRobSum'
                              , 'dep_miximp' = 'DEP'
                              , 'msstats' = 'MSstats'
                                ))

p = pd %>% ggplot + geom_point(aes(qvalue, n_proteins, colour = method))+ geom_line(aes(qvalue, n_proteins, colour = method)) + xlim(0, .5) + ylim(0,100) + facet_grid(~dataset) + geom_vline(xintercept = c(.01,.05),colour = 'grey')
p
ggsave('latosinka_comparisons_commonVSall_proteins.pdf',p,width = 12,height = 8)

p = filter(pd,dataset != 'common proteins') %>% ggplot +
  geom_point(aes(qvalue, n_proteins, colour = method), size = .5)+ 
  geom_line(aes(qvalue, n_proteins, colour = method), size = .3) + xlim(0, .205) + ylim(0,100) +
  geom_vline(xintercept = c(.01,.05),colour = 'grey') +
  xlab('False Discovery Rate') +
  ylab('Number of proteins returned') + theme_bw()
p
ggsave('latosinka_comparisons_all_proteins.pdf',p,width = 5,height = 3.5)
postscript('latosinka_comparisons_all_proteins.eps', width = 5, height = 3.5);p;dev.off()


dp = filter(pd,dataset == 'common proteins') %>% ggplot +
  geom_point(aes(qvalue, n_proteins, colour = method), size = .5)+
  geom_line(aes(qvalue, n_proteins, colour = method), size = .3) + xlim(0, .205) + ylim(0,100) +
  geom_vline(xintercept = c(.01,.05),colour = 'grey') +
  xlab('False Discovery Rate') +
  ylab('Number of proteins returned') + theme_bw()
p
ggsave('latosinka_comparisons_common_proteins.pdf',p,width = 5,height = 3.5)

dt = filter(res, method != 'perseus')
res2 = dt %>% group_by(protein) %>% filter(n() == 6) %>%
  group_by(method) %>% mutate(qvalue = p.adjust(pvalue, 'BH')) %>% mutate(n_proteins = row_number(pvalue))
pd = bind_rows(mutate(dt, dataset = 'all proteins'), mutate(res2,dataset = 'common proteins'))
pd = ungroup(pd) %>%
  mutate(method = recode_factor(method
                              , 'msqrob' = 'MSqRob'
                              , 'msqrobsum' = 'MSqRobSum'
                              , 'dep_miximp' = 'DEP'
                              , 'msstats' = 'MSstats'
                                ))

p = pd %>% ggplot + geom_point(aes(qvalue, n_proteins, colour = method))+ geom_line(aes(qvalue, n_proteins, colour = method)) + xlim(0, .5) + ylim(0,100) + facet_grid(~dataset) + geom_vline(xintercept = c(.01,.05),colour = 'grey')
p
ggsave('latosinka_comparisons_commonVSall_proteins_witquant.pdf',p,width = 12,height = 8)

###########################################################################
## Check for differencess between msqrob and msqrobsum modellen on ttest ##
###########################################################################
out1 = read_rds('../analyses/output/msqrob') %>% mutate(method = 'msqrob')
out2 = read_rds('../analyses/output/msqrobsum') %>% mutate(method = 'msqrobsum')

r = bind_rows(out1,out2) %>% inner_join(res, by = c("protein", "method"))
r = group_by(r, method) %>% mutate(pvalue_z = pnorm(-abs(t)) * 2
                                 , qvalue_z = p.adjust(pvalue_z, 'BH'))
r = group_by(r, protein) %>%
  mutate(df_post_msqrob = df_post[method == 'msqrob']) %>%
  group_by(method) %>%
  mutate(pvalue_dfmsqrob = pt(-abs(t),df_post_msqrob) * 2
       , qvalue_dfmsqrob = p.adjust(pvalue_dfmsqrob, 'BH'))

id = group_by(r, protein) %>%
  summarise(msqrob = qvalue[method == 'msqrob'], msqrobsum = qvalue[method == 'msqrobsum']
           ,msqrobsum_z = qvalue_z[method == 'msqrobsum']
           ,msqrobsum_dfmsqrob = qvalue_dfmsqrob[method == 'msqrobsum']
            ) %>%
  mutate(
   msqrob_signif = msqrob < .05
  , msqrobsum_signif = msqrobsum < .05
 , msqrobsum_signif_z = msqrobsum_z < .05
 , msqrobsum_signif_dfmsqrob = msqrobsum_dfmsqrob < .05
  )

############
## log FC ##
############
pd1 = left_join(id,r) %>%
  transmute(protein,method,var = abs(logFC),msqrob_signif,msqrobsum_signif,msqrobsum_signif_z
           ,msqrobsum_signif_dfmsqrob
            ,statistic = 'absolute log fold changes') %>%
  spread(method, var)

########
## SE ##
########
pd2 = left_join(id,r) %>%
  transmute(protein,method,var = se,msqrob_signif,msqrobsum_signif,msqrobsum_signif_z
           ,msqrobsum_signif_dfmsqrob
           ,statistic = 'standard errors') %>%
  spread(method, var)

########
## DF ##
########
pd3 = left_join(id,r) %>%
  transmute(protein,method,var = df_post,msqrob_signif,msqrobsum_signif,msqrobsum_signif_z
           ,msqrobsum_signif_dfmsqrob
           ,statistic = 'degrees of freedom') %>%
  spread(method, var)

#############
## t-value ##
#############
pd4 = left_join(id,r) %>%
  transmute(protein,method,var = abs(t),msqrob_signif,msqrobsum_signif,msqrobsum_signif_z
           ,msqrobsum_signif_dfmsqrob
           ,statistic = 'absolute t-values') %>%
  spread(method, var)

#############
id = filter(out1, !is.na(formula)) %>%
  transmute(protein
          , theta_sample = map_dbl(model,~{a = lme4::getME(.,'theta')['sample.(Intercept)']})
          , small_theta = theta_sample < .01
          , small_theta = ifelse(is.na(small_theta),FALSE,small_theta)
          , small_theta = ifelse(small_theta, 'Yes', 'No')
            )
pd = bind_rows(pd1, pd2, pd3, pd4) %>% left_join(id)
pd = pd  %>%
  filter(!((statistic == 'standard errors') & (msqrob > 1.5))) %>%
  filter(!((statistic == 'degrees of freedom') & (msqrob > 300))) %>%
  mutate(msqrobsum_significance = ifelse(msqrobsum_signif, 'Yes'
                                       , ifelse(msqrobsum_signif_dfmsqrob,'Yes\nwhen using MSqRob df','No'))
       , statistic = ifelse(statistic == 'absolute log fold changes', '(A) absolute log fold changes',
                     ifelse(statistic == 'absolute t-values', '(B) absolute t-values',
                     ifelse(statistic == 'degrees of freedom', '(C) degrees of freedom',
                     ifelse(statistic == 'standard errors', '(D) standard errors',NA
                            )))))
p = filter(pd,!msqrob_signif) %>% ggplot + geom_point(aes(msqrob,msqrobsum), colour = 'grey') +
  geom_abline(slope = 1) +
  facet_wrap(~statistic,nrow = 2,scales = 'free')
  p = p + geom_point(aes(msqrob,msqrobsum
                       , colour = small_theta
                       , shape = msqrobsum_significance
                       )
                 , filter(pd,msqrob_signif)) +
  scale_shape_manual(values = c(1,8,16))
p = p + theme_bw() +
  labs(colour = 'theta sample > 0.01', shape = 'Protein significant with\nMSqRobSum'
     , x = 'MSqRob', y = 'MSqRobSum')
ggsave('Latosinka_msqrob_vs_msqrobsum_teststatistics_ttest_dfmsqrobttestInfo_smallTheta001_alldataGrey.png',p, width = 9, height = 7);p

#### density plots to show effect of normalization on expression values
#######################################################################
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
  facet_grid(~method) + guides(colour = FALSE) +
  xlab('Leading logFC dim 1') +
  ylab('Leading logFC dim 2') + theme_bw()
pmds

dt = d %>% imap_dfr(~{MSnSet2df(.x) %>% mutate(method = .y)})
pdens = ggplot(dt) +  geom_density(aes(expression,group = sample,colour = condition)) +
  theme_bw() +
  theme(legend.position = 'top') +
  facet_grid(~method)
pdens
p = cowplot::plot_grid(pdens,pmds,ncol = 1)
ggsave('latosinka_normalization_density_mds_plots.png', p, width = 10, height = 8 );p

