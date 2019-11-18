library(tidyverse)

##########
res = list.files('../analyses/results', full.names = TRUE) %>%
    map_dfr(~{print(.x);read_rds(.x) %>% mutate(method = str_replace(.x,'.+\\/([^\\/]+)','\\1'))})

## res = group_by(res, method) %>% arrange(pvalue) %>% mutate(n_proteins = row_number(pvalue))
## res = mutate(res,pvalue = ifelse(is.na(pvalue) & (qvalue == 0),0,pvalue)) %>% group_by(method) %>%
##   arrange(pvalue) %>% mutate(n_proteins = row_number(pvalue))

## remove NA pvalues
res = filter(res,!is.na(pvalue)) %>% group_by(method) %>%
  arrange(pvalue) %>% mutate(n_proteins = row_number(pvalue))

## res = group_by(res,method) %>% mutate(qval2 = p.adjust(pvalue,'BH'))
## a = filter(res, method == 'msstats')
## plot(a$qvalue,a$qval2-a$qvalue)

## results = preprocess_results(res)
## results_filtered = group_by(res,protein, contrast) %>% filter(n() == nlevels(results$method)) %>%
##   group_by(method, contrast) %>% mutate(qvalue = p.adjust(pvalue, method = 'BH')) %>% ungroup
## results_filtered = preprocess_results(results_filtered)

###################
count(res, qvalue < .01)
count(res,method)
tst = filter(res, method == 'msstats', qvalue < .01)
tst = filter(res, method == 'dep_miximp', is.na(qvalue))
tst %>% View
head(tst$protein)

filter(res,method == 'perseus') %>% summary

filter(res) %>% ggplot + geom_line(aes(qvalue, n_proteins, colour = method))
## res %>% ggplot + geom_line(aes(qvalue, n_proteins, colour = method)) + xlim(0, .2) + ylim(0,70)
p = res %>% ggplot + geom_line(aes(qvalue, n_proteins, colour = method)) + xlim(0, .4) + ylim(0,100)
ggsave('latsinka_comparisons.pdf',p)

select(res,protein, method,qvalue) %>% spread(method,qvalue) %>% ggplot + geom_point(aes(msqrob,msqrobsum))
filter(res, pvalue < .1) %>% select(protein, method,pvalue) %>% spread(method,pvalue) %>%
  ggplot + geom_point(aes(msqrob,msqrobsum))

r = filter(res, method %in% c('msqrob','msqrobsum'))
r = ungroup(r)
id = bind_rows(
  filter(r, method == 'msqrobsum') %>% transmute(method = 'msqrob', protein,pass_other_method = qvalue < .01)
 , filter(r, method == 'msqrob') %>% transmute(method = 'msqrobsum', protein,pass_other_method = qvalue < .01)
)
pd = left_join(r, id)
pass = group_by(r, method,.drop = TRUE) %>% summarise(cutoff = max(pvalue[qvalue < .01]))
p = ggplot(pd) + geom_point(aes(logFC, -log(pvalue),colour = pass_other_method)) + facet_grid(~method) + geom_hline(aes(yintercept = -log(cutoff)),data = pass)

ggsave('latsinka_vulcano_msqrobsum.pdf',p,width = 10,height = 6)

  ggplot + geom_point(aes(msqrob,msqrobsum))
filter(res, pvalue < .1) %>% select(protein, method,pvalue) %>% spread(method,pvalue) %>%

  filter(res,method == 'perseus') %>% arrange(qvalue) %>% View

#######################################
## performance plot on same proteins ##
#######################################
res2 = filter(res,method != 'perseus') %>% group_by(protein) %>% filter(n() == 4) %>%
  group_by(method) %>% mutate(qvalue = p.adjust(pvalue, 'BH')) %>% mutate(n_proteins = row_number(pvalue))
## filter(res2) %>% ggplot + geom_line(aes(qvalue, n_proteins, colour = method))
## res %>% ggplot + geom_line(aes(qvalue, n_proteins, colour = method)) + xlim(0, .2) + ylim(0,70)
pd = bind_rows(filter(res, method != 'perseus') %>% mutate(dataset = 'all proteins'), mutate(res2,dataset = 'common proteins'))
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
  ylab('Number of proteins returned')
p
ggsave('latosinka_comparisons_all_proteins.pdf',p,width = 5,height = 3.5)

dp = filter(pd,dataset == 'common proteins') %>% ggplot +
  geom_point(aes(qvalue, n_proteins, colour = method), size = .5)+
  geom_line(aes(qvalue, n_proteins, colour = method), size = .3) + xlim(0, .205) + ylim(0,100) +
  geom_vline(xintercept = c(.01,.05),colour = 'grey') +
  xlab('False Discovery Rate') +
  ylab('Number of proteins returned')
p
ggsave('latosinka_comparisons_common_proteins.pdf',p,width = 5,height = 3.5)

res2 = filter(res,method != 'perseus') %>% group_by(protein) %>% filter(n() == 6) %>%
  group_by(method) %>% mutate(qvalue = p.adjust(pvalue, 'BH')) %>% mutate(n_proteins = row_number(pvalue))
## filter(res2) %>% ggplot + geom_line(aes(qvalue, n_proteins, colour = method))
## res %>% ggplot + geom_line(aes(qvalue, n_proteins, colour = method)) + xlim(0, .2) + ylim(0,70)
pd = bind_rows(filter(res, method != 'perseus') %>% mutate(dataset = 'all proteins'), mutate(res2,dataset = 'common proteins'))
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


########################
## Z-tests comparison ##
########################
r = filter(res, method %in% c('msqrobsum','msqrob'))
r  = group_by(r, method) %>%
  mutate(pvalue = pnorm(-abs(t)) * 2) %>% arrange(pvalue) %>% 
  mutate(qvalue = p.adjust(pvalue, 'BH')
       , n_proteins = row_number(pvalue)) %>% 
  ungroup %>% mutate(method = str_glue('{method}_z')) %>%
  bind_rows(r)
r = filter(res,method %in% c('original')) %>% bind_rows(r)

p = r %>%
  ggplot + geom_line(aes(qvalue, n_proteins, colour = method)) +
  xlim(0, .1) + ylim(0,130) + xlab('q value') + ylab('Number of proteins')
p
ggsave('Latosinka_performance_ztest_pvalue.png', p, width = 5, height = 4);p

#####################################
## Check diferences msqrob and DEP ##
#####################################
p = res %>%
  ggplot + geom_line(aes(qvalue, n_proteins, colour = method)) +
  xlim(0, .1) + ylim(0,130) + xlab('q value') + ylab('Number of proteins')

sel = filter(res,method %in% c('msqrobsum', 'dep_miximp')) %>% arrange(qvalue) %>%
  ungroup
View(head(sel))
id = sel %>% select(protein) %>% slice(1:2)

out = read_rds('../analyses/output/msqrobsum')
pd1 = select(out,protein,data) %>% unnest %>% left_join(id,.)
pd2 = select(out,protein,data_summarized) %>% unnest %>% left_join(id,.)
ggplot(pd) + geom_point(aes(sample,expression,colour = feature)) + facet_grid(~protein) +
  geom_point(aes(sample, expression),data = pd2)

library(msqrobsum)
source('../analyses/functions.r')
out2 = read_maxlfq('../data/maxquant/') %>%
  preprocess(norm = 'vsn',logtrans = FALSE)
pd3 = msqrobsum::MSnSet2df(out2) %>% select(protein,sample,expression,condition) %>% unnest %>% left_join(id,.)


pd = bind_rows( mutate(pd1,summarisation = 'raw')
        , mutate(pd2,summarisation = 'robust')
        , mutate(pd3,summarisation = 'maxlfq')
          )
ggplot(pd) + geom_point(aes(sample,expression,colour = feature,group = condition)) + facet_grid(protein ~summarisation) + guides(colour = FALSE)

pd = mutate(pd, sample = str_glue('{str_extract(sample,"...$")}_{str_extract(sample,"^..")}'))
p = ggplot(pd) + geom_point(aes(sample,expression,colour = feature)) + facet_grid(protein ~summarisation) + guides(colour = FALSE)
p
ggsave('Latosinka_intensities_DEP_msqrobsum_topprots.png', p, width = 12, height = 8);p

count(res,method)

#####################################
## check msqrob z sample variances ##
#####################################
out = read_rds('../analyses/output/msqrob')
out = filter(out,str_detect(formula,'sample')) %>%
  mutate(sample_sd = map_dbl(model,~{a = attributes(summary(.)$varcor$sample)$stddev}))
count(out, formula)
table(out$sample_sd == 0)
## FALSE  TRUE 
## 1435   164

id = filter(out,sample_sd == 0) %>% select(protein)

r = filter(res, method %in% c('msqrobsum','msqrob'))
r  = group_by(r, method) %>%
  mutate(pvalue = pnorm(-abs(t)) * 2) %>% arrange(pvalue)

## r = bind_rows(
##   mutate(r, qvalue = p.adjust(pvalue, 'BH') , n_proteins = row_number(pvalue), subset = 'all')
## , left_join(id, r) %>% mutate(qvalue = p.adjust(pvalue, 'BH')
##                             , n_proteins = row_number(pvalue), subset = 'no sample variance')
## , anti_join(r, id) %>% mutate(qvalue = p.adjust(pvalue, 'BH')
##                             , n_proteins = row_number(pvalue), subset = 'sample variance')
## )

r = mutate(r, qvalue = p.adjust(pvalue, 'BH'))
r = bind_rows(
  mutate(r, n_proteins = row_number(pvalue), subset = 'all')
, left_join(id, r) %>% mutate( n_proteins = row_number(pvalue), subset = 'no sample variance')
, anti_join(r, id) %>% mutate(n_proteins = row_number(pvalue), subset = 'sample variance')
)

r = r %>% ungroup %>% mutate(method = str_glue('{method}_z'))
count(r, method, subset)

p = r %>%  ggplot + geom_line(aes(qvalue, n_proteins, colour = method)) +
  xlim(0, .5) + ylim(0,130) + xlab('q value') + ylab('Number of proteins') +
  facet_grid(~subset)
p

### same but look at theta

out = filter(out,str_detect(formula,'sample')) %>%
  mutate(sample_theta = map_dbl(model,~{a = lme4::getME(.,'theta')[2]}))
count(out, formula)
table(out$sample_sd == 0)
summary(out$sample_theta)


id = filter(out,sample_theta < .1) %>% select(protein)

r = filter(res, method %in% c('msqrobsum','msqrob'))
r  = group_by(r, method) %>%
  mutate(pvalue = pnorm(-abs(t)) * 2) %>% arrange(pvalue)

r = bind_rows(
  mutate(r, qvalue = p.adjust(pvalue, 'BH') , n_proteins = row_number(pvalue), subset = 'all')
, left_join(id, r) %>% mutate(qvalue = p.adjust(pvalue, 'BH')
                            , n_proteins = row_number(pvalue), subset = 'no sample variance')
, anti_join(r, id) %>% mutate(qvalue = p.adjust(pvalue, 'BH')
                            , n_proteins = row_number(pvalue), subset = 'sample variance')
)

r = r %>% ungroup %>% mutate(method = str_glue('{method}_z'))
count(r, method, subset)

p = r %>%  ggplot + geom_line(aes(qvalue, n_proteins, colour = method)) +
  xlim(0, .1) + ylim(0,130) + xlab('q value') + ylab('Number of proteins') +
  facet_grid(~subset)
p

##################################################################
## Check for differencess between msqrob and msqrobsum modellen ##
##################################################################
out1 = read_rds('../analyses/output/msqrob') %>% mutate(method = 'msqrob')
out2 = read_rds('../analyses/output/msqrobsum') %>% mutate(method = 'msqrobsum')

out1 = filter(out1,str_detect(formula,'sample'))
out2 = select(out1,protein) %>% left_join(out2)

r = bind_rows(out1,out2) %>% left_join(res, by = c("protein", "method"))
r = mutate(r, theta_condition = map_dbl(model,~{a = lme4::getME(.,'theta')['condition.(Intercept)']}))
r = group_by(r, method) %>% mutate(pvalue = pnorm(-abs(t)) * 2
                                   , qvalue = p.adjust(pvalue, 'BH'))

pd = select(r,protein,method,theta_condition) %>% spread(method, theta_condition)
pd %>% ggplot + geom_point(aes(msqrob,msqrobsum))
pd %>% ggplot + geom_point(aes(msqrob,msqrobsum))+ xlim(0,20) + geom_abline(slope = 1)
pd %>% ggplot + geom_point(aes(msqrob,msqrobsum))+ xlim(0,10) + geom_abline(slope = 1)

p = ggplot(r) + geom_point(aes(qvalue,theta_condition,colour = method)) + ylim(0,7) + xlim(0,.01)
ggsave('Latosinka_qvalue_vs_conditiontheta.png', p, width = 7, height = 6);p

## id =ungroup(r) %>% filter(method == 'msqrob',qvalue < .01) %>% select(protein)
id = group_by(r, protein) %>%
  summarise(msqrob = qvalue[method == 'msqrob'], msqrobsum = qvalue[method == 'msqrobsum']) %>%
  filter(msqrob < .01) %>% mutate(msqrobsum_signif = msqrobsum < .01)

pd = left_join(id,r) %>% select(protein,method,theta_condition,msqrobsum_signif) %>% spread(method, theta_condition)
## pd %>% ggplot + geom_point(aes(msqrob,msqrobsum, colour = msqrobsum_signif)) + geom_abline(slope = 1)
p = pd %>% ggplot + geom_point(aes(msqrob,msqrobsum, colour = msqrobsum_signif)) +
  xlim(0,7) + ylim(0,7) + geom_abline(slope = 1)
ggsave('Latosinka_msqrob_conditiontheta_vs_msqrobsum_conditiontheta.png', p, width = 7, height = 6);p
## pd %>% ggplot + geom_point(aes(log(msqrob),log(msqrobsum), colour = msqrobsum_signif)) + geom_abline(slope = 1)

pd = left_join(id,r) %>% select(protein,method,logFC,msqrobsum_signif) %>% spread(method, logFC)
p = pd %>% ggplot + geom_point(aes(abs(msqrob),abs(msqrobsum), colour = msqrobsum_signif)) +
  geom_abline(slope = 1)
p
ggsave('Latosinka_msqrob_absLogFC_vs_msqrobsum_absLogFC.png', p, width = 7, height = 6);p
## pd %>% ggplot + geom_point(aes(log(msqrob),log(msqrobsum), colour = msqrobsum_signif)) + geom_abline(slope = 1)

####
r = select(r, method,protein,data) %>% unnest %>% count(method,protein, sample) %>%
group_by(protein, method) %>% summarise(sample_size_sd = sd(n)) %>% left_join(r,.)

pd = left_join(id,r) %>% select(protein,method,sample_size_sd,msqrobsum_signif) %>% spread(method, sample_size_sd)
pd %>% ggplot + geom_point(aes(msqrob,msqrobsum, colour = msqrobsum_signif)) + geom_abline(slope = 1)

ungroup(r) %>% select(protein,sample_size_sd) %>% distinct %>% left_join(id,.) %>%
  ggplot + geom_boxplot(aes(msqrobsum_signif, sample_size_sd)) +
    geom_jitter(aes(msqrobsum_signif, sample_size_sd))

## ggplot(r) + geom_point(aes((sample_size_sd), qvalue,colour = method)) + ylim(0,.01)

pd = select(r, method,protein,data) %>% unnest %>% count(method,protein, sample) %>%
  group_by(protein, method) %>% summarise(sample_size_sd = sd(n),sample_size_mean = mean(n)) %>%
  ungroup %>% select(-method) %>% distinct %>% left_join(id,.)
p = pd %>% ggplot + geom_point(aes(sample_size_mean,sample_size_sd, colour = msqrobsum_signif)) +
  xlim(0,15) + ylim(0,5)
ggsave('Latosinka_samplesizeMean_vs_samplesizeSd.png', p, width = 7, height = 6);p
## pd %>% ggplot + geom_jitter(aes(sample_size_median,sample_size_mad, colour = msqrobsum_signif),width =.1,height =.1)
##############
## Take top 10 protein with the biggest difference in logFC between msqrob-msqrobsum
arrange(id, desc(msqrobsum)) %>% View
dat = read_rds('../analyses/output/msqrobsum')
intensities = dat %>% select(protein,data) %>% unnest
summaries = dat %>% select(protein,data_summarized) %>% unnest

topprots = arrange(id, desc(msqrobsum)) %>% head(10)
topprots = mutate(topprots,
                  title = str_replace(protein, 'OS=.*$',''),
                  title = str_glue('{title}\n qvalue:    msqrob: {round(msqrob,3)*100}%  -  msqrobsum: {round(msqrobsum,3)*100}%'))

pd = left_join(topprots,intensities) %>%
  mutate(sample2 = sample %>% str_split('_') %>% map_chr(~{paste(rev(.),collapse = '')}))
pd2 = left_join(topprots,summaries) %>%
  mutate(sample2 = sample %>% str_split('_') %>% map_chr(~{paste(rev(.),collapse = '')}))
p = ggplot(pd) + geom_point(aes(sample2,expression,colour = feature))  +
  geom_line(aes(sample2,expression,group = feature,colour = feature), linetype = 2,alpha = .5) +
  facet_wrap(~title,ncol = 2,scales = 'free_y') + guides(colour=FALSE) + xlab('') +
  scale_x_discrete(breaks = unique(pd$sample2),labels = unique(pd$sample) %>% str_split('_') %>% map_chr(last)) +
  geom_point(aes(sample2,expression),pd2)  +
  geom_line(aes(sample2,expression,group = protein),pd2, linetype = 2,alpha = .5) +
  geom_vline(xintercept = 4.5)
p
ggsave('Latosinka_top10highMsqrobsumQvalueProteins_Msqrob1PercSignificant.pdf',p,width = 10,height = 20)


topprots = filter(id, !msqrobsum_signif) %>% arrange(msqrob) %>% head(10)
topprots = mutate(topprots,
                  title = str_replace(protein, 'OS=.*$',''),
                  title = str_glue('{title}\n qvalue:    msqrob: {round(msqrob,3)*100}%  -  msqrobsum: {round(msqrobsum,3)*100}%'))

pd = left_join(topprots,intensities) %>%
  mutate(sample2 = sample %>% str_split('_') %>% map_chr(~{paste(rev(.),collapse = '')}))
pd2 = left_join(topprots,summaries) %>%
  mutate(sample2 = sample %>% str_split('_') %>% map_chr(~{paste(rev(.),collapse = '')}))
p = ggplot(pd) + geom_point(aes(sample2,expression,colour = feature))  +
  geom_line(aes(sample2,expression,group = feature,colour = feature), linetype = 2,alpha = .5) +
  facet_wrap(~title,ncol = 2,scales = 'free_y') + guides(colour=FALSE) + xlab('') +
  scale_x_discrete(breaks = unique(pd$sample2),labels = unique(pd$sample) %>% str_split('_') %>% map_chr(last)) +
  geom_point(aes(sample2,expression),pd2)  +
  geom_line(aes(sample2,expression,group = protein),pd2, linetype = 2,alpha = .5) +
  geom_vline(xintercept = 4.5)
p
ggsave('Latosinka_top10lowMsqrobQvalueProteins_MsqrobNotSignificant.pdf',p,width = 10,height = 20)


write_rds(id,'msqrob1percSignif_qvaluesMsqrobMsqrobsum.rds')
write_rds(intensities,'intensities_all_prot.rds')
write_rds(summaries,'robustsummaries_all_prot.rds')


ggplot(pd) + geom_point(aes(sample2,expression))  +
  geom_line(aes(sample2,expression,group = protein), linetype = 2,alpha = .5) +
  facet_wrap(~protein,ncol = 2) + guides(colour=FALSE) + xlab('') +
  scale_x_discrete(breaks = unique(pd$sample2),labels = unique(pd$sample) %>% str_split('_') %>% map_chr(last))


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


#####
id = group_by(r, protein) %>%
  summarise(msqrob = qvalue[method == 'msqrob'], msqrobsum = qvalue[method == 'msqrobsum']) %>%
  filter(msqrobsum < .01) %>% mutate(msqrob_signif = msqrobsum < .01)

pd = left_join(id,r) %>% select(protein,method,theta_condition,msqrob_signif) %>% spread(method, theta_condition)
pd %>% ggplot + geom_point(aes(msqrob,msqrobsum, colour = msqrob_signif)) + geom_abline(slope = 1)
pd %>% ggplot + geom_point(aes(msqrob,msqrobsum, colour = msqrob_signif)) +
  xlim(0,7) + ylim(0,7) + geom_abline(slope = 1)



summary(r$pvalue)


o = r$model[[2000]]
lme4::getME(o,'theta')['condition.(Intercept)']

glimpse(out)
glimpse(res)
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
  ## filter(msqrob < .05) %>%
  mutate(
   msqrob_signif = msqrob < .05
  , msqrobsum_signif = msqrobsum < .05
 , msqrobsum_signif_z = msqrobsum_z < .05
 , msqrobsum_signif_dfmsqrob = msqrobsum_dfmsqrob < .05
  )

pd = left_join(id,r) %>%
  transmute(protein,method,var = abs(logFC),msqrob_signif,msqrobsum_signif,msqrobsum_signif_z
           ,msqrobsum_signif_dfmsqrob
            ,statistic = 'absolute log fold change') %>%
  spread(method, var)
pd1 = pd
## p = pd %>% ggplot + geom_point(aes(abs(msqrob),abs(msqrobsum), colour = msqrobsum_signif)) + geom_abline(slope = 1)
## p
## p  + facet_grid(~msqrobsum_signif_z)

p = pd %>% ggplot + geom_point(aes(msqrob,msqrobsum, colour = msqrobsum_signif
                                 , shape = msqrobsum_signif_z)) +
  geom_abline(slope = 1) + scale_shape_manual(values = c(1,20))
p
ggsave('Latosinka_msqrob_absLogFC_vs_msqrobsum_absLogFC_ttest.png',p, width = 7, height = 6);p

###############
## SE
###########
pd = left_join(id,r) %>%
  transmute(protein,method,var = se,msqrob_signif,msqrobsum_signif,msqrobsum_signif_z
           ,msqrobsum_signif_dfmsqrob
           ,statistic = 'standard error') %>%
  spread(method, var)
pd2 = pd
p = pd %>% ggplot + geom_point(aes(msqrob,msqrobsum, colour = msqrobsum_signif
                                 , shape = msqrobsum_signif_z)) +
  geom_abline(slope = 1) + scale_shape_manual(values = c(1,20))
p
ggsave('Latosinka_msqrob_se_vs_msqrobsum_absLogFC_ttest.png',p, width = 7, height = 6);p

pd = left_join(id,r) %>%
  transmute(protein,method,var = df_post,msqrob_signif,msqrobsum_signif,msqrobsum_signif_z
           ,msqrobsum_signif_dfmsqrob
           ,statistic = 'degrees of freedom') %>%
  spread(method, var)
pd3 = pd

pd = left_join(id,r) %>%
  transmute(protein,method,var = abs(t),msqrob_signif,msqrobsum_signif,msqrobsum_signif_z
           ,msqrobsum_signif_dfmsqrob
           ,statistic = 'absolute t value') %>%
  spread(method, var)
pd4 = pd

pd = bind_rows(pd1, pd2, pd3, pd4) %>% filter(msqrob_signif)
p = pd %>% ggplot + geom_point(aes(msqrob,msqrobsum, colour = msqrobsum_signif
                                 , shape = msqrobsum_signif_z)) +
  geom_abline(slope = 1) + scale_shape_manual(values = c(1,16)) + facet_wrap(~statistic,nrow = 2,scales = 'free')
p
ggsave('Latosinka_msqrob_vs_msqrobsum_teststatistics_ttest_ztestInfo.png',p, width = 9, height = 7);p


pd = bind_rows(pd1, pd2, pd3, pd4) #%>% filter(msqrob_signif)
p = pd %>% ggplot + geom_point(aes(msqrob,msqrobsum, colour = msqrobsum_signif
                                 , shape = msqrobsum_signif_z)) +
  geom_abline(slope = 1) + scale_shape_manual(values = c(1,16)) + facet_grid(msqrob_signif~statistic,scales = 'free')
p

p = filter(pd,msqrob_signif) %>% ggplot + geom_point(aes(msqrob,msqrobsum, colour = msqrobsum_signif)) +
                               geom_abline(slope = 1) +
  facet_wrap(~statistic,nrow = 2,scales = 'free')
p
ggsave('Latosinka_msqrob_vs_msqrobsum_teststatistics_ttest.png',p, width = 9, height = 7);p


pd = bind_rows(pd1, pd2, pd3, pd4) %>%
  filter(!((statistic == 'standard error') & (msqrob > 5))) %>%
  filter(!((statistic == 'degrees of freedom') & (msqrob > 300))) 
p = filter(pd,!msqrob_signif) %>% ggplot + geom_point(aes(msqrob,msqrobsum), colour = 'grey') +
  geom_abline(slope = 1) +
  facet_wrap(~statistic,nrow = 2,scales = 'free')
p = p + geom_point(aes(msqrob,msqrobsum, colour = msqrobsum_signif
                                 , shape = msqrobsum_signif_z), filter(pd,msqrob_signif)) +
  scale_shape_manual(values = c(1,16))
p
ggsave('Latosinka_msqrob_vs_msqrobsum_teststatistics_ttest_ztestInfo_alldataGrey.png',p, width = 9, height = 7);p

pd = bind_rows(pd1, pd2, pd3, pd4) %>%
  filter(!((statistic == 'standard error') & (msqrob > 5))) %>%
  filter(!((statistic == 'degrees of freedom') & (msqrob > 300))) 
p = filter(pd,!msqrob_signif) %>% ggplot + geom_point(aes(msqrob,msqrobsum), colour = 'grey') +
  geom_abline(slope = 1) +
  facet_wrap(~statistic,nrow = 2,scales = 'free')
p = p + geom_point(aes(msqrob,msqrobsum, colour = msqrobsum_signif
                     , shape = msqrobsum_signif_dfmsqrob), filter(pd,msqrob_signif)) +
  scale_shape_manual(values = c(1,16))
p
ggsave('Latosinka_msqrob_vs_msqrobsum_teststatistics_ttest_dfmsqrobttestInfo_alldataGrey.png',p, width = 9, height = 7);p

######################333333
id = filter(out1, !is.na(formula)) %>%
  transmute(protein, theta_sample = map_dbl(model,~{a = lme4::getME(.,'theta')['sample.(Intercept)']}))

pd = left_join(pd2, id) %>%
  filter(theta_sample < 3)
##   filter(!((statistic == 'degrees of freedom') & (msqrob > 300))) 
p = filter(pd,!msqrob_signif) %>%
  ## ggplot + geom_point(aes(log(theta_sample),log(msqrobsum/msqrob)), colour = 'grey') +
  ggplot + geom_point(aes(theta_sample,log(msqrobsum/msqrob)), colour = 'grey') +
  geom_abline(slope = 0) +
  facet_wrap(~statistic,nrow = 2,scales = 'free')
p
p = p + geom_point(aes((theta_sample),log(msqrobsum/msqrob), colour = msqrobsum_signif
 ## p = p + geom_point(aes(log(theta_sample),log(msqrobsum/msqrob), colour = msqrobsum_signif
                     , shape = msqrobsum_signif_dfmsqrob), filter(pd,msqrob_signif)) +
  scale_shape_manual(values = c(1,16))
p
ggsave('Latosinka_msqrob_vs_msqrobsum_theta_logfcSE_ttest_dfmsqrobttestInfo_alldataGrey.png',p, width = 6, height = 4);p

pd = left_join(pd2, id)
  ## filter(theta_sample < 3000)
##   filter(!((statistic == 'degrees of freedom') & (msqrob > 300))) 
p = filter(pd,!msqrob_signif) %>%
  ggplot + geom_point(aes((theta_sample),(msqrobsum/msqrob)), colour = 'grey') +
  geom_abline(intercept = 1,slope = 0) +
  facet_wrap(~statistic,nrow = 2,scales = 'free')
p
p = p + geom_point(aes((theta_sample),(msqrobsum/msqrob), colour = msqrobsum_signif
                     , shape = msqrobsum_signif_dfmsqrob), filter(pd,msqrob_signif)) +
  scale_shape_manual(values = c(1,16)) + ylim(0,5) +xlim(0,5)
p
ggsave('Latosinka_msqrob_vs_msqrobsum_theta_fcSE_ttest_dfmsqrobttestInfo_alldataGrey.png',p, width = 8, height = 5);p

## p = filter(pd,msqrob_signif,abs(log(msqrobsum/msqrob)) < 10)  %>% ggplot +
##        ## geom_point(aes((theta_sample),log(msqrobsum/msqrob,base = 10), colour = msqrobsum_signif
##   geom_point(aes((theta_sample),msqrobsum/msqrob, colour = msqrobsum_signif
##                        ## p = p + geom_point(aes(log(theta_sample),log(msqrobsum/msqrob), colour = msqrobsum_signif
##                      , shape = msqrobsum_signif_dfmsqrob)) +
##        geom_abline(slope = 0) +
##   scale_shape_manual(values = c(1,16))
## p

## ggsave('Latosinka_msqrob_vs_msqrobsum_theta_logfcSE_ttest_dfmsqrobttestInfo_alldataGrey.png',p, width = 8, height = 5);p


## id = filter(out1, !is.na(formula)) %>%
##   transmute(protein
##           , theta_sample = map_dbl(model,~{a = lme4::getME(.,'theta')['sample.(Intercept)']})
##           , small_theta = theta_sample < .01
##           , small_theta = ifelse(is.na(small_theta),FALSE,small_theta))
## pd = bind_rows(pd1, pd2, pd3, pd4) %>% left_join(id)

## pd = pd  %>%
##   filter(!((statistic == 'standard error') & (msqrob > 1))) %>%
##   filter(!((statistic == 'degrees of freedom') & (msqrob > 300))) 
## p = filter(pd,!msqrob_signif) %>% ggplot + geom_point(aes(msqrob,msqrobsum), colour = 'grey') +
##   geom_abline(slope = 1) +
## facet_wrap(small_theta~statistic,scales = 'free',nrow = 2)
## p = p + geom_point(aes(msqrob,msqrobsum, colour = msqrobsum_signif
##                      , shape = msqrobsum_signif_dfmsqrob
##                        )
##                  , filter(pd,msqrob_signif)) +
##   scale_shape_manual(values = c(1,16))
## p = p + labs(colour = 'Theta sample < 0.01',shape = 'Significant at 5% FDR\nwith MSqRobSum')
## ggsave('Latosinka_msqrob_vs_msqrobsum_teststatistics_ttest_dfmsqrobttestInfo_smallTheta001_alldataGrey.png',p, width = 9, height = 7);p

pd = pd  %>%
  filter(!((statistic == 'standard error') & (msqrob > 1.5))) %>%
  filter(!((statistic == 'degrees of freedom') & (msqrob > 300))) %>%
  mutate(msqrobsum_significance = ifelse(msqrobsum_signif, 'Yes', ifelse(msqrobsum_signif_dfmsqrob,'Yes\nwhen using MSqRob df','No')))
p = filter(pd,!msqrob_signif) %>% ggplot + geom_point(aes(msqrob,msqrobsum), colour = 'grey') +
  geom_abline(slope = 1) +
  facet_wrap(~statistic,nrow = 2,scales = 'free')
  p = p + geom_point(aes(msqrob,msqrobsum
                       ## , colour = msqrobsum_signif
                       , colour = small_theta
                     , shape = msqrobsum_significance
                       )
                 , filter(pd,msqrob_signif)) +
  scale_shape_manual(values = c(1,8,16))
p
ggsave('Latosinka_msqrob_vs_msqrobsum_teststatistics_ttest_dfmsqrobttestInfo_smallTheta001_alldataGrey.png',p, width = 9, height = 7);p

;;

## aa = dplyr::filter(out1, grepl('O00154',out1$protein))$model[[1]]

## table(grepl('sp|O00154',out1$protein))
## aa = out1$model[[1]]
## aa
## aa = dplyr::filter(out1, grepl('P13798',out1$protein))$model[[1]]
## aa = dplyr::filter(out1, grepl('P54577',out1$protein))$model[[1]]

## aa = dplyr::filter(out1, grepl('O75131',out1$protein))$model[[1]]
## aa = dplyr::filter(out1, grepl('O75608',out1$protein))$model[[1]]
## aa = dplyr::filter(out1, grepl('O94905', out1$protein))$model[[1]]

## dplyr::filter(out1, grepl('P00505', out1$protein))$model[[1]]
## dplyr::filter(out1, grepl('P02656', out1$protein))$model[[1]]
## dplyr::filter(out1, grepl('P13489', out1$protein))$model[[1]]
## dplyr::filter(out1, grepl('P49720', out1$protein))$model[[1]]
## dplyr::filter(out1, grepl('P25205', out1$protein))$model[[1]]
## dplyr::filter(out1, grepl('Q5JTV8', out1$protein))$model[[1]]

## dplyr::filter(out2, grepl('Q7L5L3', out1$protein))$model[[1]]
## dplyr::filter(out2, grepl('Q7L5L3', out1$protein))$model[[1]]

## m = dplyr::filter(out1, grepl('Q9Y623', out1$protein))$model[[1]]
## ms = dplyr::filter(out2, grepl('Q9Y623', out1$protein))$model[[1]]
## msqrobsum::msqrobsum

## m
## ms
## summary(m)
## summary(ms)
## m@frame 
## ms@frame 
## m = dplyr::filter(out1, grepl('Q9Y623', out1$protein))
## glimpse(m)


## m = dplyr::filter(out1, grepl('Q9Y623', out1$protein))$model
## ms = dplyr::filter(out2, grepl('Q9Y623', out1$protein))$model

###################################
## make figures for all proteins ##
###################################
out1 = read_rds('../analyses/output/msqrob') %>% mutate(method = 'msqrob')
out2 = read_rds('../analyses/output/msqrobsum') %>% mutate(method = 'msqrobsum')

## r = bind_rows(out1,out2) %>% inner_join(res, by = c("protein", "method"))
## r = group_by(r, method) %>% mutate(pvalue_z = pnorm(-abs(t)) * 2
##                                  , qvalue_z = p.adjust(pvalue_z, 'BH'))
id = group_by(r, protein) %>%
  summarise(msqrob = qvalue[method == 'msqrob']
          , msqrobsum = qvalue[method == 'msqrobsum']
          , msqrobsum_z = qvalue_z[method == 'msqrobsum']
          , msqrobsum_dfmsqrob = qvalue_dfmsqrob[method == 'msqrobsum']
            ) %>%
  filter(msqrob < .05) %>%
  mutate(msqrobsum_signif = msqrobsum < .05
       , msqrobsum_signif_z = msqrobsum_z < .05
       , msqrobsum_signif_dfmsqrob = msqrobsum_dfmsqrob < .05
         )

## pd = left_join(id,r) %>% rowwise %>%
##   transmute(method,protein, data, data_summarized
##           , statistics = list(list(qvalue = qvalue, pvalue = pvalue,t = t, df = df_post, se = se, logFC = logFC))) %>%
##   group_by(protein) %>%
##   summarise(data = first(data), data_summarized = data_summarized[method == 'msqrobsum'])


## d = select(ranking, protein) %>% left_join(pd) %>% slice(1:2)

make_plot = function(d) {
  d1 = filter(d, method == 'msqrob')
  d2 = filter(d, method == 'msqrobsum')

  pd1 = d1 %>% pull(data) %>% first %>%
    mutate(sample2 = sample %>% str_split('_') %>% map_chr(~{paste(rev(.),collapse = '')}))

  pd2 = d2 %>% select(protein, data_summarized) %>% unnest %>%
    mutate(sample2 = sample %>% str_split('_') %>% map_chr(~{paste(rev(.),collapse = '')}))

  title = str_replace(d1$protein, 'OS=.*$','')

  ## add_info = function(d){
  ##   str_glue('{ifelse(d$method == "msqrob","MSqRob:  ","MSqRobSum:")} q-value: {round(d$qvalue,3)*100}%  log fold change: {round(d$logFC,2)}  t-value: {round(d$t,2)}  df: {round(d$df,1)}  se: {round(d$se,2)}')
  ## }
  ## subtitle = str_glue('{add_info(d1)}\n{add_info(d2)}')
  a = transmute(d, method, `q-value` = round(qvalue,3)*100
               ,`log fold change` = logFC
              , `t-value`=round(t,3)
              , df
              , se = round(se,3)
               ,`theta sample` = round(theta_sample,3))
  a[2,'theta sample'] = NA
  subtitle = capture.output(a)[-c(1,3)] %>% str_sub(3) %>%
    str_replace_all(' `','   ') %>%
    str_replace_all('`','') %>%
    paste0(collapse = '\n')

  p = ggplot(pd1) + geom_point(aes(sample2,expression,colour = feature))  +
    geom_line(aes(sample2,expression,group = feature,colour = feature), linetype = 2,alpha = .5) +
    ## facet_wrap(~title,ncol = 2,scales = 'free_y') +
    guides(colour=FALSE) + xlab('') +
    scale_x_discrete(breaks = unique(pd1$sample2)
                   , labels = unique(pd1$sample) %>% str_split('_') %>% map_chr(last)) +
    geom_point(aes(sample2,expression),pd2)  +
    geom_line(aes(sample2,expression, group = protein),pd2, linetype = 2,alpha = .5) +
    geom_vline(xintercept = count(pd2,condition) %>% arrange(condition) %>% {.$n[1] +.5}) +
    ggtitle(title,subtitle = subtitle) +
    theme(plot.subtitle = element_text(family = 'mono'))
  p
  ## ggsave('test.pdf',p,width = 8.5,height = 8.5)
}



pd = left_join(id,r) %>%
  transmute(method,protein, data, data_summarized, qvalue,pvalue, t, df = df_post, se, logFC)

pd = filter(out1, !is.na(formula)) %>% transmute(protein, theta_sample = map_dbl(model,~{a = lme4::getME(.,'theta')['sample.(Intercept)']})) %>% left_join(pd,.)


## I want the msqrobsum significant prots first
ranking = select(pd,protein,method,qvalue) %>% spread(method,qvalue) %>% arrange(msqrobsum > .05, msqrob)
pdf("Latosinka_Msqrob1PercSignifif.pdf", 8.5, 8.5)
## pdf("test.pdf", 8.5, 8.5)
mutate(ranking,id = row_number()) %>% select(id,protein) %>% left_join(pd) %>%
  group_by(id) %>% group_walk(~{print(make_plot(.))})
dev.off()

ranking = select(pd,protein,method,pvalue) %>% spread(method,pvalue) %>%
  mutate(diff = abs(msqrob-msqrobsum)) %>% arrange(diff)
pdf("Latosinka_Msqrob1PercSignifif_ranked_diff_pvalue.pdf", 8.5, 8.5)
mutate(ranking,id = row_number()) %>% select(id,protein) %>% left_join(pd) %>%
  group_by(id) %>% group_walk(~{print(make_plot(.))})
dev.off()

ranking = select(pd,protein,method,df) %>% spread(method,df) %>% mutate(diff = abs(msqrob-msqrobsum)) %>% arrange(diff)
pdf("Latosinka_Msqrob1PercSignifif_ranked_diff_df.pdf", 8.5, 8.5)
mutate(ranking,id = row_number()) %>% select(id,protein) %>% left_join(pd) %>%
  group_by(id) %>% group_walk(~{print(make_plot(.))})
dev.off()

ranking = select(pd,protein,method,qvalue) %>% spread(method,qvalue) %>%
  filter(msqrobsum < .05) %>% arrange(msqrobsum)
pdf("Latosinka_Msqrob5PercSignif_bothsignif.pdf", 8.5, 8.5)
mutate(ranking,id = row_number()) %>% select(id,protein) %>% left_join(pd) %>%
  group_by(id) %>% group_walk(~{print(make_plot(.))})
dev.off()

r = bind_rows(out1,out2) %>% inner_join(res, by = c("protein", "method"))
r = group_by(r, method) %>% mutate(pvalue_z = pnorm(-abs(t)) * 2
                                 , qvalue_z = p.adjust(pvalue_z, 'BH'))
id = group_by(r, protein) %>%
  summarise(msqrob = qvalue[method == 'msqrob'], msqrobsum = qvalue[method == 'msqrobsum']
           ,msqrobsum_z = qvalue_z[method == 'msqrobsum']) %>%
  filter(msqrob < .05) %>%
  mutate(msqrobsum_signif = msqrobsum < .05
       , msqrobsum_signif_z = msqrobsum_z < .05)

ranking = filter(id,!msqrobsum_signif, msqrobsum_signif_z)
pdf("Latosinka_Msqrob5PercSignif_msqrobsum_z_signif.pdf", 8.5, 8.5)
mutate(ranking,id = row_number()) %>% select(id,protein) %>% left_join(pd) %>%
  group_by(id) %>% group_walk(~{print(make_plot(.))})
dev.off()

ranking = filter(id,!msqrobsum_signif, !msqrobsum_signif_z)
pdf("Latosinka_Msqrob5PercSignif_msqrobsum_z_nonsignif.pdf", 8.5, 8.5)
mutate(ranking,id = row_number()) %>% select(id,protein) %>% left_join(pd) %>%
  group_by(id) %>% group_walk(~{print(make_plot(.))})
dev.off()
##################
ranking = filter(id,!msqrobsum_signif, msqrobsum_signif_dfmsqrob)
pdf("Latosinka_Msqrob5PercSignif_msqrobsum_dfmsqrob_signif.pdf", 8.5, 8.5)
mutate(ranking,id = row_number()) %>% select(id,protein) %>% left_join(pd) %>%
  group_by(id) %>% group_walk(~{print(make_plot(.))})
dev.off()

ranking = filter(id,!msqrobsum_signif, !msqrobsum_signif_dfmsqrob)
pdf("Latosinka_Msqrob5PercSignif_msqrobsum_dfmsqrob_nonsignif.pdf", 8.5, 8.5)
mutate(ranking,id = row_number()) %>% select(id,protein) %>% left_join(pd) %>%
  group_by(id) %>% group_walk(~{print(make_plot(.))})
dev.off()

################################################################################








ranking1 = select(pd,protein,method,qvalue) %>% spread(method,qvalue) %>%
  filter(msqrobsum > .05) %>% arrange(msqrobsum)
ranking = left_join(ranking1,pd) %>% select(protein,method,t) %>% spread(method,t) %>%
  mutate(diff = abs(msqrob-msqrobsum)) %>% arrange(diff)

pdf("Latosinka_Msqrob1PercSignifif_ranked_diff_t.pdf", 8.5, 8.5)
mutate(ranking,id = row_number()) %>% select(id,protein) %>% left_join(pd) %>%
  group_by(id) %>% group_walk(~{print(make_plot(.))})
dev.off()

ranking2 = left_join(ranking1,pd) %>% select(protein,method, t) %>% spread(method, t) %>%
  mutate(diff_t = abs(msqrob-msqrobsum))

ranking3 = left_join(ranking1,pd) %>% select(protein,method,df) %>% spread(method,df) %>%
  mutate(diff_df = abs(msqrob-msqrobsum))

ranking4 = left_join(ranking1,pd) %>% select(protein,method,se) %>% spread(method,se) %>%
  mutate(diff_se = abs(msqrob-msqrobsum))

p1 = inner_join(select(ranking2, protein,diff_t), select(ranking3, protein,diff_df)) %>% ggplot + geom_point(aes(diff_t,diff_df))

p2 = inner_join(select(ranking4, protein,diff_se), select(ranking3, protein,diff_df)) %>% ggplot + geom_point(aes(diff_se,diff_df))

cowplot::plot_grid(p1,p2)
pdf("Latosinka_Msqrob5PercSignif_mMsqrobSumNonSignif__diff_t_diff_se_diff_df.pdf")

##########################################
r = bind_rows(out1,out2) %>% inner_join(res, by = c("protein", "method"))
r = group_by(r, method) %>% mutate(pvalue_z = pnorm(-abs(t)) * 2
                                 , qvalue_z = p.adjust(pvalue_z, 'BH'))
id = group_by(r, protein) %>%
  summarise(msqrob = qvalue[method == 'msqrob'], msqrobsum = qvalue[method == 'msqrobsum']
           ,msqrobsum_z = qvalue_z[method == 'msqrobsum']) %>%
  filter(msqrob < .05) %>%
  mutate(msqrobsum_signif = msqrobsum < .05
       , msqrobsum_signif_z = msqrobsum_z < .05)

ranking = filter(id,!msqrobsum_signif, msqrobsum_signif_z)
pdf("Latosinka_Msqrob5PercSignif_msqrobsum_z_signif.pdf", 8.5, 8.5)
mutate(ranking,id = row_number()) %>% select(id,protein) %>% left_join(pd) %>%
  group_by(id) %>% group_walk(~{print(make_plot(.))})
dev.off()

ranking = filter(id,!msqrobsum_signif, !msqrobsum_signif_z)
pdf("Latosinka_Msqrob5PercSignif_msqrobsum_z_nonsignif.pdf", 8.5, 8.5)
mutate(ranking,id = row_number()) %>% select(id,protein) %>% left_join(pd) %>%
  group_by(id) %>% group_walk(~{print(make_plot(.))})
dev.off()
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
  facet_grid(~method) + guides(colour = FALSE) +
  xlab('Leading logFC dim 1') +
  ylab('Leading logFC dim 2')
pmds

dt = d %>% imap_dfr(~{MSnSet2df(.x) %>% mutate(method = .y)})
pdens = ggplot(dt) +  geom_density(aes(expression,group = sample,colour = condition)) +
  theme(legend.position = 'top') +
  facet_grid(~method)
pdens
p = cowplot::plot_grid(pdens,pmds,ncol = 1)
ggsave('normalization_density_mds_plots.png', p, width = 10, height = 8 );p
