library('tidyverse')
library('MSnbase')

### read mq dataset
#####################
read_mq = function(txt_path,fasta_path) {
  path_pep = paste0(txt_path,'/peptides.txt')
  path_prot = paste0(txt_path,'/proteinGroups.txt')

  print(str_glue('\n\nREADING MQ DATA'))

  exprs_col = grepEcols(path_pep, 'Intensity ',split = '\t')
  set = readMSnSet2(path_pep ,ecol = exprs_col,fnames = 'Sequence', sep = '\t',stringsAsFactors = FALSE)
  sampleNames(set) = str_replace(sampleNames(set),'Intensity.','')
  fd = fData(set) %>%
    transmute(protein = Proteins,
              sequence = Sequence,
              reverse = ifelse(is.na(Reverse),'',Reverse) == '+',
              contaminant = grepl('^CON_',protein))

  fd = add_id(fd,fasta_path)

  fd <-read_tsv(path_prot) %>%
    transmute(protein = `Protein IDs`, site_only = `Only identified by site`,
              ecoli = grepl('OS=Escherichia coli', `Fasta headers`)) %>%
    right_join(fd) %>%
    mutate(site_only = ifelse(is.na(site_only),'' ,site_only) == '+') %>%
    as.data.frame
  rownames(fd) = featureNames(set)

  pd = data.frame(condition = as.factor(str_extract(sampleNames(set),'^.'))
                , batch = factor(c(2,1,2,2,1,1,2,2,1,1,1,2,1,1,2,2,1,1,2,2)))
  rownames(pd) = sampleNames(set)

  set = MSnSet(exprs(set), fData =  AnnotatedDataFrame(fd), pData = AnnotatedDataFrame(pd))

  print(str_glue('dimension msnset: {paste(dim(set),collapse = " ")}'))
  print(str_glue('design'))
  print(count(pData(set),condition))
  print(str_glue('summary pData'))
  print(summary(pData(set)))
  print(str_glue('summary fData'))
  print(summary(fData(set)))

  set
}

### read maxlfq dataset
#####################
read_maxlfq = function(txt_path,fasta_path) {
  ## sometimes proteins with _REV in name are not labeled as Reverse because the majority protein is not Rev
  path_prot = paste0(txt_path,'proteinGroups.txt')

  print(str_glue('\n\nREADING MQ MAXLFQ DATA'))

    exprs_col = grepEcols(path_prot, 'LFQ.intensity ',split = '\t')
  set = readMSnSet2(path_prot ,ecol = exprs_col,fnames = 'Protein.IDs', sep = '\t',stringsAsFactors = FALSE)
  sampleNames(set) = str_replace(sampleNames(set),'LFQ.intensity.','')

  fd = fData(set) %>%
    transmute(protein = Protein.IDs,
              reverse = ifelse(is.na(Reverse),'',Reverse) == '+'
             ,contaminant = grepl('^CON_',protein)
            , site_only = Only.identified.by.site
            , site_only = ifelse(is.na(site_only),'' ,site_only) == '+'
              )

  fd = add_id(fd,fasta_path)

  rownames(fd) = featureNames(set)

  pd = data.frame(condition = as.factor(str_extract(sampleNames(set),'^.'))
                , batch = factor(c(2,1,2,2,1,1,2,2,1,1,1,2,1,1,2,2,1,1,2,2)))
  rownames(pd) = sampleNames(set)

  set = MSnSet(exprs(set), fData =  AnnotatedDataFrame(fd), pData = AnnotatedDataFrame(pd))

  print(str_glue('dimension msnset: {paste(dim(set),collapse = " ")}'))
  print(str_glue('design'))
  print(count(pData(set),condition))
  print(str_glue('summary pData'))
  print(summary(pData(set)))
  print(str_glue('summary fData'))
  print(summary(fData(set)))

  set
}

### preprocess data
###################
log_filter = function(msg,set,sel = FALSE){
  print(str_glue('\n\n- {msg}'))
  if(any(sel)) print(str_glue('{sum(!sel)} out of {length(sel)} sequences are filtered out'))
  print(str_glue('dimension msnset: {paste(dim(set),collapse = " ")}'))
  print(str_glue('{filter(fData(set),ecoli)$protein %>% unique %>% length} ecoli proteins found'))
}

preprocess = function(set, sep = ';', norm = 'quantiles', remove_condition_with_1_sample = TRUE,
                      remove_1_peptide_proteins = FALSE, logtrans = TRUE
                    , two_sample_condition = TRUE){
  log_filter('PREPROCESSING',set)

  exprs(set)[0 == (exprs(set))] <- NA
  nas = (is.na(exprs(set)))
  print(str_glue('{sum(nas)} out of {length(nas)} intensities are missing (= are NA)'))

  ## log transform data
#####################
  if (logtrans){
    print(str_glue('\n\n- log 2 transform (zero intensities assigned NA)'))
    set <- log(set, base=2)
  }
  ## quantile normalisation
  #########################
  if (!isFALSE(norm)){
    print(str_glue('\n\n- {norm} normalisation'))
    set <- normalize(set, norm)
  }

  ## keep smallest unique groups
  ##############################
  groups = tibble(protein = fData(set)$protein) %>% distinct %>%
      ## count the number of proteins in a protein group
      mutate(proteins = strsplit(protein, ';'), n = lengths(proteins)) %>% unnest %>%
      ## Check for every protein what the smallest protein group is, it belongs to
      ## remove larger protein groups
      group_by(proteins) %>% filter(n == min(n)) %>% ungroup %>%
      count(protein,n) %>% filter(n == nn) %>% pull(protein)

  sel <- fData(set)$protein %in% groups
  set <- set[sel,]
  log_filter('keep smallest unique groups',set,sel)


  ## Remove reverse
#################
  if('reverse' %in% fvarLabels(set)){
    sel <- !fData(set)$reverse
    summary(sel) %>% print
    set <- set[sel,]
    log_filter('remove reverse sequences',set,sel)
  }

  ## Remove contaminants
  ######################
  sel <- !fData(set)$contaminant
  set <- set[sel]
  log_filter('remove contaminants',set,sel)

  ## Remove proteins only identified with modification
  ######################################################
  if('site_only' %in% fvarLabels(set)){
    sel <- !fData(set)$site_only
    set <- set[sel]
    log_filter('remove proteins only identified with a modification',set,sel)
  }

  ## Remove peptides that assigned to both human and ecoli
########################################################
  sel <- fData(set) %>% transmute(flag = !(human & ecoli)) %>% pull(flag)
  set <- set[sel]
  log_filter('remove proteins assigned to both human and ecoli',set,sel)

  ## Remove peptides that assigned to unknown protein
###################################################
  sel <- fData(set)$protein != ''
  set <- set[sel]
  log_filter('remove `UNKNOWN` proteins assigned',set,sel)

  ## filter out conditions that have not at least 2 samples
  ### make them NA in msnset object
  ## and remove peptides with less then 2 observations
###########################################################
  while(two_sample_condition) {
    df = MSnSet2df(set)
    id <- df %>%  group_by(protein,condition) %>% summarise(n = length(unique(sample))) %>%
      filter(n < 2) %>% ungroup %>% select(protein, condition) %>%
      left_join(df, by = c('protein', 'condition')) %>% select(feature,sample)
    if(nrow(id) ==0) break
    exprs(set)[as.matrix(id)] = NA

    sel <- rowSums(!is.na(exprs(set))) >= 2
    set <- set[sel]
  }

  ## and remove peptides with less then 2 observations
  sel <- rowSums(!is.na(exprs(set))) >= 2
  set <- set[sel]
  log_filter('remove peptides with less then 2 observations and less then 2 samples per condition',set)

  ## remove proteins identified with only 1 peptide
####################
  if (remove_1_peptide_proteins){
    sel = fData(set) %>% group_by(protein) %>% mutate(flag = n()>1) %>% pull(flag)
    set <- set[sel]
    log_filter('remove proteins identified with only 1 peptide',set,sel)
  }
  fData(set) = fData(set)[,! colnames(fData(set)) %in%
                           c('contaminant','reverse', 'site_only','sequence')]
  set
}




########
rmm2res = function(obj){
  obj %>%
    transmute(protein,ecoli,human,contrasts) %>%
    filter(!map_lgl(contrasts,is.null)) %>%
    unnest %>%
    mutate(contrast = str_replace_all(contrast,'condition','')) %>%
    separate(contrast, c('sample1','sample2'),remove = FALSE)
}

limma2res = function(obj){
  obj$result %>% filter(!is.na(qvalue)) %>%
    mutate(contrast = str_replace_all(contrast,'condition','')) %>%
    separate(contrast, c('sample1','sample2'),remove = FALSE)
}


dep2res = function(res,fasta_path){
  res %>% 
    gather(statistic,value, contains('_vs_')) %>%
    mutate(statistic = str_replace(statistic, "_vs_",'vs')) %>%
    separate(statistic,c('contrast', 'statistic'),sep = '_') %>%
    spread(statistic,value) %>%
    transmute(protein = ID
            , feature = ID
              ## contrasts are differently specified
            , logFC = -ratio
            , pvalue = p.val
            , qvalue = p.adj
            , contrast) %>%
    separate(contrast, c('sample2', 'sample1'),sep = 'vs') %>%
    mutate(contrast = str_glue('{sample1}-{sample2}')) %>%
    add_id(fasta_path) %>% as_data_frame
}



perseus2res = function(perseus_file,fasta_path){
  perseus_file %>%
    select(protein, contains('T-test')) %>% #group_by(protein) %>%
    gather(statistic,value,contains('T-test')) %>%
    mutate(statistic = str_replace(statistic, "Student's T-test ",'')
         , statistic = str_replace(statistic, "Test statistic",'t')) %>%
    filter(statistic != 'significant') %>%
    separate(statistic,c('statistic','contrast'),sep = ' ') %>%
    mutate(contrast = str_replace(contrast,'_','-')) %>%
    spread(statistic,value) %>%
    transmute(protein
            , feature = protein
            , contrast = contrast
            , sample1 = str_replace(contrast,'-.','')
            , sample2 = str_replace(contrast,'.-','')
            , logFC = as.numeric(Difference)
            , t = as.numeric(t)
            , pvalue = as.numeric(`p-value`)
            , qvalue = as.numeric(`q-value`)
              ) %>%
    add_id(fasta_path) %>% as_data_frame
}



## add ecoli and human to protein id to fData dataframe with protein column
add_id = function(fd, fasta_path){
  id = list(ecoli = paste0(fasta_path, '/ecoli_up000000625_7_06_2018.fasta'),
            human = paste0(fasta_path, '/human_up000005640_sp_7_06_2018.fasta')) %>%
    map(~{read_lines(.x) %>%
            {.[str_detect(.,'^>')]} %>%
            str_extract(.,'(?<=\\|).*(?=\\|)')})
  fd2 = fd %>% transmute(protein = as.character(protein), proteins = strsplit(protein, ';')) %>% unnest %>%
    mutate(human = proteins %in% id$human, ecoli =  proteins %in% id$ecoli) %>% group_by(protein) %>%
    summarise(human = any(human), ecoli = any(ecoli)) %>%
    right_join(fd) %>% as.data.frame
  rownames(fd2) = rownames(fd)
  fd2
}
############################
DEP_read = function(txt_path){
    proteinGroups <- read.table(str_glue("{txt_path}proteinGroups.txt"), sep="\t", header=TRUE)
    annot = data_frame(condition = rep(letters[1:5], each = 4)
                     , replicate = rep(1:4,5)
                     , label = paste0(condition,replicate))

   import_MaxQuant(proteinGroups,annot)
}


  DEP_process = function(proteins,fasta_path){
  ## DEP preprocess on maxlfq
##################################
  ## txt_path = paste0("/home/st/Documents/quant_sum/ionstar_data/txt/")
  proteins <- filter_missval(proteins)
  norm <- normalize_vsn(proteins)

  proteins_MNAR <- get_df_long(norm) %>% as_data_frame %>%
    group_by(name, condition) %>%
    summarize(NAs = all(is.na(intensity))) %>%
    filter(NAs) %>%
    pull(name) %>%
    unique()

  ## Get a logical vector
  MNAR <- names(norm) %in% proteins_MNAR
  print('Number of proteins that are Missing Not At Random')
  print(table(MNAR))

  ## Perform a mixed imputation
  mixed_imputation <- DEP::impute(
                             norm, 
                             fun = "mixed",
                             randna = !MNAR, # we have to define MAR which is the opposite of MNAR
                             mar = "knn", # imputation function for MAR
                             mnar = "MinProb") # imputation function for MNAR

  protset_maxlfq_dep = as(mixed_imputation, 'MSnSet')
  featureNames(protset_maxlfq_dep) = fData(protset_maxlfq_dep)$ID

  fData(protset_maxlfq_dep) = fData(protset_maxlfq_dep) %>% select(protein = ID) %>%
    add_id(fasta_path)
  protset_maxlfq_dep
  }

DEP_preprocess = function(txt_path,fasta_path){
  DEP_read(txt_path) %>% DEP_process(fasta_path)
}
#####
msstats2res = function(obj, fasta_path){
  as_data_frame(obj$ComparisonResult) %>%
    transmute(protein = as.character(Protein)
             ,contrast = Label
            , sample1 = str_replace(contrast,'-.','')
            , sample2 = str_replace(contrast,'.-','')
            , logFC = log2FC
            , t = Tvalue
            , pvalue
            , qvalue = adj.pvalue
            , se = SE
            , df = DF
              ) %>%
    add_id(fasta_path) %>% as_data_frame
}
#############
DEP_impute = function(set_dep){
    fd = mutate(fData(set_dep)
              , name = featureNames(set_dep)
              , ID = name)
    pd = mutate(pData(set_dep)
              , replicate = rep(1:4,5)
              , label = paste0(condition,replicate))
    e = 2^exprs(set_dep) %>% as_data_frame %>% bind_cols(fd,.)
    colid = 6:25
    ## colid = 9:28
    print(colnames(e)[colid])
    set_dep = make_se(e, colid, pd)

    proteins_MNAR <- get_df_long(set_dep) %>% as_data_frame %>%
        group_by(name, condition) %>%
        summarize(NAs = all(is.na(intensity))) %>% 
        filter(NAs) %>% 
        pull(name) %>% 
        unique()

    ## Get a logical vector
    MNAR <- names(set_dep) %in% proteins_MNAR
    table(MNAR)

    ## Perform a mixed imputation
    mixed_imputation <- DEP::impute(
                                 set_dep, 
                                 fun = "mixed",
                                 randna = !MNAR, # we have to define MAR which is the opposite of MNAR
                                 mar = "knn", # imputation function for MAR
                                 mnar = "MinProb") # imputation function for MNAR

    set_dep = as(mixed_imputation, 'MSnSet')
    featureNames(set_dep) = fData(set_dep)$ID
   set_dep
}
