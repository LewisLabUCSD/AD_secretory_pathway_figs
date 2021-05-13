source('ppi_assist_func.R')
source('graph_assist.R')
# source('GSEA.R')
# library(BiocParallel)
library(igraph)
# library(plotly)
library(data.table)
library(magrittr)
library(ggplot2)
# 
# protein.data.source <- 'deep_proteome'
# secM.expr.col <- 'Proteins_iBAQ.sigmoid'
# secM.expr.col <- 'Proteins_iBAQ.rank.gene'
# secM.expr.col <- 'Proteins_iBAQ.sigmoid.byTissue'
# secM.expr.col <- 'Proteins_iBAQ.rank.byTissue'
# secM.expr.col <- 'Proteins_iBAQ.qtnorm.rank.byTissue'
# use.tissue.median <- T
# used.int.g <- 'Humannet'
load_components <- function(protein.data.source = NULL,
                            secM.expr.col = NULL,
                            use.tissue.median = NULL,
                            used.int.g = NULL,
                            return.omics.only = F, 
                            use.subgraph=T, n.cores=32,
                            debug=F) {
  
  
  if(is.null(protein.data.source)){
    message('protein.data.source not set, defaulting to HPA')
    protein.data.source <- 'HPA'
  }else{
    message(sprintf('using %s for protein.data.source', protein.data.source))
  }
  
  if(is.null(use.tissue.median)){
    message('use.tissue.median not set, defaulting to FALSE')
    use.tissue.median <- F
  }else{
    message(sprintf('using %s for use.tissue.median', use.tissue.median))
  }
  
  if(!return.omics.only){
    if(is.null(secM.expr.col)){
      if (protein.data.source=='HPA'){
        message('secM.expr.col not set, defaulting to HPA.protein.Level.numeric for HPA')
        secM.expr.col <- 'HPA.protein.Level.numeric'
      }else if (protein.data.source == 'deep_proteome'){
        message('secM.expr.col not set, defaulting to Proteins_iBAQ.rank for deep_proteome')
        secM.expr.col <- 'Proteins_iBAQ.rank.gene'
      }
    }else{
      message(sprintf('using %s for secM.expr.col', secM.expr.col))
    }
    if(is.null(used.int.g)){
      message('used.int.g not set, defaulting to PCNet')
      used.int.g = 'PCNet'
    }else{
      message(sprintf('using %s for used.int.g', used.int.g))
    }
  }
  
  secM.amir <<- fread(here::here('databases','Human_components_feizi.csv')) %>% setnames('gene_name', 'Symbol')
  
  ## fuzzy annotation for secretory subsystems
  secM.amir[grepl('trafficking', Subsystem, ignore.case = T), Subsystem.aggr:= 'trafficking']
  secM.amir[grepl('^COPI|^COPII', Subsystem, ignore.case = T), Subsystem.aggr:= 'COPI/COPII']
  secM.amir[grepl('^COPI|^COPII', Subsystem, ignore.case = T), Subsystem.aggr:= 'COPI/COPII']
  secM.amir[is.na(Subsystem.aggr), Subsystem.aggr:= Subsystem]
  
  
  ## could potentially use all membrane bound/ secretome in the HPA dataset
  
  # ```{r getting secP components}
  HPA.secretory_membrane <- fread(here::here('databases','HPA_secreted_proteins.tsv'))[, c('HPA.predicted.secretory', 
                                                                                           'HPA.predicted.membrane', 
                                                                                           'HPA.predicted.type'):=.(T, F, 'secreted')] %>% 
    rbind(fread(here::here('databases','HPA_membrane_proteins.tsv'))[, c('HPA.predicted.secretory', 
                                                                         'HPA.predicted.membrane', 
                                                                         'HPA.predicted.type'):=.(F, T, 'membrane')]) %>% 
    rbind(fread(here::here('databases','HPA_membrane_secreted.tsv'))[, c('HPA.predicted.secretory', 
                                                                         'HPA.predicted.membrane', 
                                                                         'HPA.predicted.type'):=.(T, T, 'membrane.secreted')]) %>% 
    .[, .(`Gene name` = Gene, Gene = Ensembl,
          HPA.predicted.secretory,
          HPA.predicted.membrane,
          HPA.predicted.type)] %>% 
    .[, .SD[1], by = c('Gene name', 'HPA.predicted.secretory', 'HPA.predicted.membrane', 'HPA.predicted.type')]## remove fuzzy id mapping 
  
  HPA.secretory_membrane <<- HPA.secretory_membrane[!Gene%in%c('ENSG00000283992', 'ENSG00000255292') ] # "LYNX1" "SDHD"  have ambiguous ensembl mappings. Removed.
  
  HPA.nx <<- fread(here::here('databases', 'rna_tissue_consensus.tsv'))
  
  ### getting secretory pathway resident genes
  subcellular_location <- fread(here::here('databases','HPA', 'subcellular_location.tsv'), na.strings = '') %>%
    .[, parse.localiation(.SD), by = c('Gene', 'Gene name')]
  HPA.secretory.loc <- subcellular_location[localization%in%c('Aggresome', 'Vesicles', 'Golgi apparatus', 'Endosome', 'Lysosomes', 'Endoplasmic reticulum')]
  uniprot.subcell <- fread(here::here('databases', 'uniprot_localization_keywork.tab')) %>% 
    .[, c('Subcellular.location', 
          'human_symbol'):= .(tolower(`Subcellular location [CC]`), 
                              gsub('_HUMAN', '', `Entry name`))]
  uniprot.subcell.keywords <- uniprot.subcell[, .(Keywords = tolower(parse.semicol.sep(Keywords))), by = c('Entry','human_symbol')]
  uniprot.secretory.loc <- uniprot.subcell.keywords[Keywords%like%c('aggresome|cytoplasmic vesicle|golgi apparatus|endosome|lysosome|endoplasmic reticulum')]
  all.secretory.resident.genes <<- c(HPA.secretory.loc$`Gene name`, uniprot.secretory.loc$human_symbol) %>% unique()
  
  # # getting proteomics datasets
  # 1. HPA; categorical
  # 2. http://www.nature.com/articles/nature13319#s1 coverage:6105 genes; protein expression values (log10 normalized iBAQ, +10 log10 units)
  # 3. Matthias Selbach more than 5,000 genes in mammalian cells
  # 
  # below data pre-processing from 171129_seqBasedFeatures.Rmd
  # 
  # Of all 723 unique genes, 522 have corresponding measurements in the HPA samples for which both RNA and proteins were quantified. 
  # of all 575 secMs, all of them have corresponding entries in the HPA RNA-seq dataset. However, only 468 of the secMs are in tissue(s) for which both RNA and protein levels were quantified.
  # ```{r HPA data load}
  
  if(protein.data.source == 'HPA'){
    HPA_normal_tissue <- fread("databases/HPA/normal_tissue.tsv") %>%  ## combat ambiguous tissue naming
      .[, Tissue:= gsub(pattern = '[0-9]*$', '', x = Tissue) %>% gsub(pattern = ' $', '', x=.)] %>% 
      .[,Level:=factor(Level, levels = c("Not detected","Low", "Medium", "High"), ordered = T)] %>% 
      .[, .(HPA.protein.Level = base::max(Level)), by = c('Gene name', 'Tissue')] ## for each sample in a tissue, chosse the expression from the cell type of the highest protein yield.
    
    ###
    
    HPA_rna_tissue <- fread('databases/HPA/rna_tissue.tsv')[, .(Gene, `Gene name`, Tissue = Sample, HPA.tpm = Value)] %>% # Sample column treated as tissues identifiers
      .[, .(HPA.tpm = max(HPA.tpm)), by = c('Gene name', 'Tissue')]
    ## discretize HPA_rna_tissue.2 based on HPA.tpm quantile
    
    HPA_rna_tissue[HPA.tpm==0, HPA.tpm.quantile := 0]
    HPA_rna_tissue[HPA.tpm>0, HPA.tpm.quantile:=cut(HPA.tpm, breaks=quantile(HPA.tpm, probs=seq(0,1, by=1/3), na.rm=TRUE), 
                                                    include.lowest=TRUE, labels = 1:3) %>% as.numeric]
    
    HPA.data <- merge(HPA_normal_tissue, HPA_rna_tissue, by = c('Gene name', 'Tissue'))
    # log transformation of tpms
    HPA.data[, HPA.tpm_log:=log(1+HPA.tpm)]
    
    HPA.data.all <- merge(HPA_normal_tissue, HPA_rna_tissue, by = c('Gene name', 'Tissue'), all = T) 
    stopifnot(nrow(HPA.data.all[is.na(HPA.tpm)&Tissue%in%HPA_rna_tissue$Tissue])==0)
    HPA.data.all<- HPA.data.all[Tissue%in%intersect(HPA_normal_tissue$Tissue, HPA_rna_tissue$Tissue)]
    # log transformation of tpms
    HPA.data.all[, c('HPA.tpm_log', 'HPA.protein.Level.numeric'):=.(log(1+HPA.tpm), as.numeric(HPA.protein.Level) - 1)]
    HPA.data.all[, imputed:=F]
    HPA.data.all[is.na(HPA.protein.Level), c('imputed', 'HPA.protein.Level.numeric'):=.(T, HPA.tpm.quantile)]
    
    ## transformations
    HPA.data.all[, HPA.tpm_sigmoid:= HPA.tpm_log %>% sigmoid %>% `-`(.5) %>% `*`(2)]
    HPA.data.all[, HPA.tpm_sigmoid_scaled:= scale(HPA.tpm_log)[,1] %>% sigmoid]
    HPA.data.all[, HPA.tpm_sigmoid_scaled_gene:= scale(HPA.tpm_log)[,1] %>% sigmoid, by = 'Gene name']
    # ```
    
    # ```{r getting secM components}
    secM.expr.HPA <- HPA.data.all %>% 
      .[secM.amir[, .(ensgid, Symbol, Subsystem.aggr, Subsystem)], on = c(`Gene name` = 'Symbol')] %>% 
      .[!is.na(HPA.tpm)]
    
    # secP.Johan.HPA <- merge(sample.info, HPA.data, by.x = 'protein', by.y = 'Gene name')
    
    HPA.intracellular <- HPA.data.all[, -c('imputed','HPA.tpm.quantile', 'HPA.protein.Level.numeric',
                                           'HPA.tpm_sigmoid',
                                           'HPA.tpm_sigmoid_scaled', 
                                           'HPA.tpm_sigmoid_scaled_gene'
    )] %>% 
      .[!HPA.secretory_membrane, on = c('Gene name')] %>% 
      .[, c('HPA.predicted.secretory', 
            'HPA.predicted.membrane', 
            'HPA.predicted.type',
            'Tissue',
            'HPA.protein.Level',
            'HPA.tpm',
            'HPA.tpm_log'):=.(
              F,
              F,
              'intracellular',
              NULL,
              NULL,
              NULL,
              NULL
            )] %>% 
      .[, .SD[1], by = c('Gene name')]
    secP.HPA <<- rbind(HPA.secretory_membrane[, -'Gene'],
                       HPA.intracellular)
    
    used.tissues <- HPA.data.all[,.N, by = .(Tissue)][N>10000, Tissue]
    if(use.tissue.median){
      predicted.secP.HPA <<- merge(secP.HPA, HPA.data.all[Tissue%in%used.tissues] %>%
                                     .[,lapply(.SD, median), by = 'Gene name', .SDcols = sapply(., is.numeric)] %>% 
                                     .[,Tissue:='tissue.median'], by = c('Gene name'))
      
      protein.expr.by.tissue <<- HPA.data.all[Tissue%in%used.tissues] %>%
        .[,lapply(.SD, median), by = 'Gene name', .SDcols = sapply(., is.numeric)] %>% 
        .[,Tissue:='tissue.median'] %>%
        split(by = 'Tissue')
    } else{
      predicted.secP.HPA <<- merge(secP.HPA, HPA.data.all, by = c('Gene name'))
      protein.expr.by.tissue <<- HPA.data.all[Tissue%in%used.tissues] %>% split(by = 'Tissue')
    }
    
    secP.components <<- secP.HPA$`Gene name` %>% unique()
    secM.components <<- secM.amir$Symbol %>% unique()
    
  }else if(protein.data.source=='deep_proteome'){
    library(readxl)
    
    proteins.wide <- read_excel("databases/deep_proteome_2019/Table_EV1.xlsx", 
                                sheet = "C. Genes") %>% as.data.table() %>% 
      .[, lapply(.SD, mean),.SDcols = `Adrenal gland`:`Urinary bladder`, by = .(`Gene ID`, `Gene name`)]
    
    
    proteins.genes <- proteins.wide %>%
      melt(id.vars = c('Gene ID','Gene name'), variable.name = 'Tissue', value.name = 'Proteins_iBAQ')
    
    proteins.genes.qtnorm <- 
      proteins.wide[, `Adrenal gland`:`Urinary bladder`] %>% as.matrix() %>% { 
        `colnames<-`( preprocessCore::normalize.quantiles(.), colnames(.))
      } %>% cbind(proteins.wide[, `Gene ID`:`Gene name`],.) %>% 
      melt(id.vars = c('Gene ID','Gene name'), variable.name = 'Tissue', value.name = 'Proteins_iBAQ.qtnorm')
    
    
    transcripts.genes <- read_excel("databases/deep_proteome_2019/Table_EV2.xlsx", 
                                    sheet = "B. Genes") %>% as.data.table() %>% 
      .[, lapply(.SD, mean),.SDcols = `Adrenal gland`:`Urinary bladder`, by = .(`Gene ID`, `Gene name`)] %>%
      melt(id.vars = c('Gene ID','Gene name'), variable.name = 'Tissue', value.name = 'Genes_FPKM')
    
    protein.tx.all <<- merge(transcripts.genes,
                             merge(proteins.genes, proteins.genes.qtnorm, by = c('Gene ID', 'Gene name', 'Tissue')),
                             by = c('Gene ID', 'Gene name', 'Tissue'))
    # rank transformations
    protein.tx.all[, Genes_FPKM.rank:= rank(Genes_FPKM), by = Tissue]
    protein.tx.all[, Proteins_iBAQ.rank.gene:= rank(Proteins_iBAQ), by = Tissue]
    protein.tx.all[, Proteins_iBAQ.rank.byTissue:= rank(Proteins_iBAQ), by = .(`Gene ID`, `Gene name`)]
    protein.tx.all[, Proteins_iBAQ.qtnorm.rank.byTissue:= rank(Proteins_iBAQ.qtnorm), by = .(`Gene ID`, `Gene name`)]
    
    protein.tx.all[, Proteins_iBAQ.sigmoid := sigmoid(scale(log(1+Proteins_iBAQ)))]
    protein.tx.all[ , Proteins_iBAQ.sigmoid.byTissue := sigmoid(scale(log(1+Proteins_iBAQ))), by = .(`Gene ID`, `Gene name`)]
    protein.tx.all[ , Proteins_iBAQ.qtnorm.sigmoid.byTissue := sigmoid(scale(log(1+Proteins_iBAQ.qtnorm))), by = .(`Gene ID`, `Gene name`)]
    
    protein.tx.all.sec <<- merge(protein.tx.all, HPA.secretory_membrane[, .('Gene ID' = Gene, HPA.predicted.type)],by = 'Gene ID', all.x = T)
    protein.tx.all.sec[is.na(HPA.predicted.type), HPA.predicted.type:='intracellular']
    predicted.secP.HPA <<- 
      protein.tx.all.sec %>% {
        merge(.[,lapply(.SD, mean), by = .(`Gene name`, Tissue), .SDcols = sapply(., is.numeric)],
              .[ , .SD[1], by = .(`Gene name`)] %>% .[, .(`Gene name`, HPA.predicted.type)], by = 'Gene name')
      }
    
    
    if(use.tissue.median){
      predicted.secP.HPA <<- predicted.secP.HPA %>% {
        merge(.[,lapply(.SD, median), by = 'Gene name', .SDcols = sapply(., is.numeric)] ,
              .[, .(`Gene name`, HPA.predicted.type)] %>% unique, by = 'Gene name')  ## retaining HPA.predicted.type
      } %>% 
        .[,Tissue:='tissue.median']
    }
    
    protein.expr.by.tissue <<- predicted.secP.HPA %>% split(by = 'Tissue')
    secP.components <<- predicted.secP.HPA$`Gene name` %>% unique()
    secM.components <<- protein.tx.all %>% .[secM.amir[, .(ensgid, Symbol, Subsystem.aggr)], on = c(`Gene ID` = 'ensgid')] %>% 
      .[!is.na(`Gene name`), unique(`Gene name`)]
  }
  
  
  
  if(return.omics.only) return()
  if(!use.subgraph) all.secretory.resident.genes <- NULL
  # ######
  #
  # start loading RWR3
  #
  # ######
  debug.RWR3 = F
  
  generate.new = F
  # n.cores = 32
  use.weights = F
  # scale.expr.RWR3 = F
  scale.expr.RWR3 = T
  # secM.expr.col = 'HPA.tpm_log'
  
  # secM.expr.col = 'HPA.tpm_sigmoid_scaled'
  # secM.expr.col = 'HPA.tpm_sigmoid_scaled'
  source.expr.override = NA ## use mean to override
  n.propagation = 20
  a.propagation = .1
  debug = debug
  if(debug){
    n.cores <- 1
  }
  
  # used.int.g = 'bioplex'
  ## setting up save directory
  version = sprintf('RWR4_190222_%s_%s%s', protein.data.source, used.int.g, ifelse(use.subgraph, yes = '', no='_fullGraph'))
  save.dir <- sprintf('output/180224_graphRWR/%s/%s', version, ifelse(use.tissue.median, 'tissue.median', 'allTissues'))
  if(!dir.exists(save.dir)) dir.create(save.dir, recursive = T)
  prop.int.score.dt.fn <- sprintf('%s/%s_prop.int.score.dt.csv.zip', save.dir, secM.expr.col)
  prop.intensity.all.save.fn <- sprintf('%s/%s_prop.intensity.all.fst', save.dir, secM.expr.col)
  
  if(!debug & file.exists(prop.int.score.dt.fn)){
    prop.int.score.dt <<- rio::import(prop.int.score.dt.fn) %>% as.data.table
  }else{
    if(!debug & file.exists(prop.intensity.all.save.fn)){
      RWR.dt <<- fst::read.fst(prop.intensity.all.save.fn) %>% as.data.table
      
    }else{
      
      int.g <- readRDS('output/int.db.g.RDS')[[used.int.g]] # generated by 171202_graphBasedClustering.Rmd
      write_graph(int.g, sprintf('auxiliary_file/graph_based/load_components_files/%s.txt', used.int.g))
      RWR.dt <<- parallel::mclapply(secP.components, get.RWR2.RWR3, secM.components = secM.components,
                                    source.expr.override = source.expr.override,
                                    all.secretory.resident.genes = all.secretory.resident.genes,
                                    protein.expr.by.tissue = protein.expr.by.tissue, secM.expr.col = secM.expr.col,
                                    scale.expr.RWR3 = scale.expr.RWR3,
                                    geneSymbol.col = 'Gene name',
                                    n.propagation = n.propagation,
                                    a.propagation = a.propagation,
                                    int.g = int.g, mc.cores = n.cores) %>% rbindlist()
      RWR.dt %>%  fst::write_fst(prop.intensity.all.save.fn)
    }
    
    prop.int.score.dt <<-
      rbind(RWR.dt[ , list(int.g = used.int.g,
                           use.weights = use.weights,
                           int.score.mean.RWR3 = mean(intensity.RWR3),
                           int.score.mean.RWR2 = mean(intensity.RWR),
                           int.score.diff.RWR3.RWR2 = mean(intensity.RWR3 - intensity.RWR),
                           int.score.diff.RWR3.RWR_half = mean(intensity.RWR3 - intensity.RWR_half),
                           int.score.diff.RWR3.stationary = mean(intensity.RWR3 - stationary.dist),
                           secM.expr.col = secM.expr.col,
                           collapse.by = 'secP'),
                    by = c('secP', 'Tissue')] %>%
              setnames('secP', 'Symbol'),
            RWR.dt[ , list(int.g = used.int.g,
                           use.weights = use.weights,
                           int.score.mean.RWR3 = mean(intensity.RWR3),
                           int.score.mean.RWR2 = mean(intensity.RWR),
                           int.score.diff.RWR3.RWR2 = mean(intensity.RWR3 - intensity.RWR),
                           int.score.diff.RWR3.RWR_half = mean(intensity.RWR3 - intensity.RWR_half),
                           int.score.diff.RWR3.stationary = mean(intensity.RWR3 - stationary.dist),
                           secM.expr.col = secM.expr.col,
                           collapse.by = 'secM'),
                    by = c('secM', 'Tissue')] %>%
              setnames('secM', 'Symbol'))
    
    rio::export(prop.int.score.dt, prop.int.score.dt.fn)
  }
}


