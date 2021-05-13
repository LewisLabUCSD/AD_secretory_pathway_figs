library(data.table)
data.table::setDTthreads(32)
library(magrittr)
library(ggplot2)
#### Assign multiple new variables 
# http://stackoverflow.com/questions/7519790/assign-multiple-new-variables-on-lhs-in-a-single-line-in-r
# Generic form
'%=%' = function(l, r, ...) UseMethod('%=%')
# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}

# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}

# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}


#' Title
#'
#' @param ll 
#'
#' @return
#' @export
#'
#' @examples 
#' 
# ## An example of an 'unbalanced' list
# z <- list(z1 = list(a = 1, b = 2), 
#           z2 = list(b = 4, a = 1, c = 0))
list.nested.flip <-  function(ll) {
  nms <- unique(unlist(lapply(ll, function(X) names(X))))
  ll <- lapply(ll, function(X) setNames(X[nms], nms))
  ll <- apply(do.call(rbind, ll), 2, as.list)
  lapply(ll, function(X) X[!sapply(X, is.null)])
}


#' fuzzy id conversion using the Genome wide annotation for Human 
#'
#' @param to.convert vector of ids to convert
#' @param to.form format to convert to
#' @param from.form can be inferred from number of matches with the annotation table
#'
#' @return
#' @export
#'
#' @examples
fuzzy.id.conversion <- function(to.convert, to.form, from.form = NULL, use.db = 'org.Hs.eg.db'){
  require(use.db, character.only = T)
  to.convert %<>% as.character() 
  
  ## fuzzy naming convention
  fuzzy.key.type <- function(key.type) {
    key.type <- toupper(key.type)
    if (grepl(pattern = 'SYMBOL', key.type)) return('SYMBOL')
    if (grepl(pattern = 'ID', key.type)&grepl(pattern = 'GENE|ENTREZ', key.type)) return('ENTREZID')
    if (grepl(pattern = 'ENSEMBL', key.type)) return('ENSEMBL')
    return(key.type)
  }
  
  to.form %<>% fuzzy.key.type()
  if(is.null(from.form)){
    message('inferring format from which to convert')
    from.form <- 
      sapply(setdiff(c('ENSEMBL', 'ENTREZID', 'SYMBOL', 'UNIPROT',  "ENSEMBLPROT",  "ENSEMBLTRANS",
                       'PROTEINID', 'GENEID', 'GENENAME'), to.form),
             function(try.from.form){
               tryCatch(select(get(use.db), to.convert, to.form, try.from.form) %>% nrow, error = function(x) 0)
             }
      ) %>% which.max %>% names
  }else{
    from.form %<>% fuzzy.key.type()
  }
  message(sprintf('converting %s to %s', from.form, to.form))
  select(get(use.db), to.convert, to.form, from.form) %>% as.data.table() %>% setnames(c('before', 'after')) %>% setkey('before')
}

id.conversion.table <- function(dt, id.col, to.form, from.form = NULL, use.db = 'org.Hs.eg.db', converted.col.name = paste0(id.col, '.', to.form)){
  to.convert <- unique(dt[[id.col]])
  conversion.dict <- fuzzy.id.conversion(to.convert = to.convert, to.form = to.form, from.form = from.form, use.db = use.db) %>% setnames('after', converted.col.name)
  merge(dt, conversion.dict, by.x = id.col, by.y = 'before', allow.cartesian=TRUE)
}

#' A vectorized string occurence counter
#'
#' @param x 
#' @param pattern 
#'
#' @return
#' @export
#'
#' @examples
count.string.match <- function(x, pattern, ...) {
  regmatches(x, gregexpr(pattern = pattern, text = x, ...)) %>% lengths
}

#' parse semi-colon separated strings.
#'
#' @param features 
#'
#' @return
#' @export
#'
#' @examples
parse.semicol.sep <- function(features, split.text = '; |;') {
  feature.split <- strsplit(features, split.text) %>% unlist %>% gsub(x = ., pattern = "^ |^\\?", replacement = '')
}



######


## implmentation of sensible complete case selection.
# imp_Vis <- function(dt, plot = F, completeness = .75) {
#   ## using numeric columns only
#   ## beware trade-off between number of complete cases and variables to retain
#   dt <- dt[,sapply(dt, is.numeric), with = F]
#   
#   imp_combo <- VIM::aggr(dt, cex.axis=.5,cex.lab = .5, numbers = T, col=c('navyblue','yellow'), gap = 0.5, sortVars = T)
#   imp_combo.dt <- data.table(count = imp_combo$count, combo = lapply( strsplit(imp_combo$combinations, ':'), as.numeric))
#   imp_combo.dt[,impsum:=sapply(combo, sum)]
#   selected_combo <- imp_combo.dt[impsum<ncol(dt)*imps_thresh,]
#   names(dt)[colSums(selected_combo[, transpose(combo)])==0]
# }
## chi-square calculation
chisq.test.pvalue <- function (x, y = NULL, correct = TRUE,
                               p = rep(1/length(x), length(x)), rescale.p = FALSE,
                               simulate.p.value = FALSE, B = 2000){
  if (NA %in% x | NaN %in% x| sum(x)==0){
    return(NA)
  }
  else{
    return (chisq.test(x, y, correct, p, rescale.p, simulate.p.value, B)$p.value)
  }
}
chisq.test.statistic <- function (x, y = NULL, correct = TRUE,
                                  p = rep(1/length(x), length(x)), rescale.p = FALSE,
                                  simulate.p.value = FALSE, B = 2000){
  if (NA %in% x | NaN %in% x | sum(x)==0){
    return(NA)
  }
  else{
    return (chisq.test(x, y, correct, p, rescale.p, simulate.p.value, B)$statistic)
  }
}




get_int_table <- function(dt_rows, dt_cols = NULL, int_database, int_db_cols, ID_col = 'uniprot_id', convert_to = NULL, use.Long = T){
  if(!'Weight'%in%names(int_database)) int_database[, Weight:= 1]
  int_db_colA = int_db_cols[1]
  int_db_colB = int_db_cols[2]
  row_IDs = if (is.vector(dt_rows)) unique(dt_rows) else unique(dt_rows[[ID_col]])
  col_IDs = if (is.vector(dt_cols)) unique(dt_cols) else if(is.null(dt_cols)) int_database[,unique(c(get(int_db_colA),get(int_db_colB)))] else unique(dt_cols[[ID_col]])
  int_database.a <- copy(int_database[(get(int_db_colA)%in%row_IDs & get(int_db_colB)%in%col_IDs)|
                                        (get(int_db_colB)%in%row_IDs & get(int_db_colA)%in%col_IDs)])
  int_database.b <- copy(int_database.a)
  int_database.b[[int_db_colA]] <- int_database.a[[int_db_colB]]
  int_database.b[[int_db_colB]] <- int_database.a[[int_db_colA]]
  int_combined <- rbind(int_database.a, int_database.b)
  int_combined <- int_combined[, c('Weight', int_db_colA, int_db_colB), with = F]
  int_combined <- int_combined[complete.cases(int_combined)]
  int_combined <- int_combined[!duplicated(int_combined)]
  
  if (use.Long){
    indexed.int <- as.data.table(expand.grid(row_IDs, col_IDs)) %>%
      ## all pairwise combinations of dt_cols and dt_rows with filled NAs
      setnames(1:2, c(int_db_colA,int_db_colB)) %>%
      ## semi-merge with int_combined
      int_combined[., on = c(int_db_colA, int_db_colB)]
  } else{
    indexed.int <- int_combined[get(int_db_colA)%in%row_IDs & get(int_db_colB)%in%col_IDs]
  }
  
  if(!is.null(convert_to)){
    if (anyDuplicated(dt_cols, by = convert_to) | anyDuplicated(dt_rows, by = convert_to)){
      print('duplication incurred during ID conversion, extracting first occurrence') ##todos
      dt_rows <- dt_rows[!duplicated(get(convert_to))]
      dt_cols <- dt_cols[!duplicated(get(convert_to))]
    }
    
    setkeyv(dt_cols, ID_col)
    setkeyv(dt_rows, ID_col)
    ## check_duplicates
    indexed.int[, eval(int_db_colA):= dt_rows[.(get(int_db_colA)), get(convert_to)]]
    indexed.int[, eval(int_db_colB):= dt_cols[.(get(int_db_colB)), get(convert_to)]]
  }
  if (nrow(indexed.int)==0) stop( sprintf('no interactions between geneList query: %s and secM\n', dt_cols))
  
  int_database.int <- dcast(indexed.int,
                            get(int_db_colA)~get(int_db_colB), 
                            value.var = 'Weight')
  
  
  ## factor to string for GSEA
  int.df <- indexed.int[!is.na(Weight)][,c(int_db_colA, int_db_colB):=list(as.character(get(int_db_colA)), as.character(get(int_db_colB)))] %>%
    .[,list(flattened_list=list(get(int_db_colB))) , by =.(get(int_db_colA))]
  flattened_PPI_map <- lapply(lapply(split(int.df[,2], int.df[,1], drop = T), unname), unlist)
  
  
  # return(int_combined)            ## convert to string
  return(list(indexed.int = indexed.int,
              int.table = int_database.int, 
              flattened_PPI_map = flattened_PPI_map))
}

merge.prefix <- function(df1, df2, prefixes = c('x__', 'y__'), postfix = F, ...){
  mrg <- data.table:::merge.data.table(x = df1, y = df2, ...)
  if (postfix){
    setnames(mrg, paste0(names(mrg), ifelse(names(mrg) %in% setdiff(names(df1),names(df2)),prefixes[1],"")))
    setnames(mrg, paste0(names(mrg), ifelse(names(mrg) %in% setdiff(names(df2),names(df1)),prefixes[2],"")))
  }else{
    setnames(mrg, paste0(ifelse(names(mrg) %in% setdiff(names(df1),names(df2)),prefixes[1],""),names(mrg)))
    setnames(mrg, paste0(ifelse(names(mrg) %in% setdiff(names(df2),names(df1)),prefixes[2],""),names(mrg)))
  }
  return(mrg)
}


ID_conversion <- function(ids,fromType = 'UNIPROT', OrgDb = org.Mm.eg.db, toType = 'SYMBOL', par = F){
  ## splitting ids
  if(length(ids)==1){
    # suppose ids are synonyms separated by wither ';', ' ' or ','
    ids = gsub(pattern = ';|,|\ ', x = ids, replacement = ';')
    converted = tryCatch(names(sort(table(bitr(unlist(strsplit(ids, ';')), fromType = fromType,
                                               OrgDb = OrgDb, toType = toType)[,2])))[1], error = function(e){NULL}) ## choose according to majority votes.
  } else{
    
    if(par){
      
      # Calculate the number of cores
      no_cores <- detectCores() - 1
      # Initiate cluster
      # cl <- makeCluster(no_cores, type="FORK")
      converted <- unlist(mclapply(X = ids, FUN = ID_conversion, fromType = fromType, OrgDb = OrgDb, toType = toType, par = F, par = F, mc.cores = no_cores))
      # stopCluster(cl)
      
    }else{
      converted <- sapply(ids, ID_conversion, fromType = fromType, OrgDb = OrgDb, toType = toType, par = F)
    }
  }
  if(toType == 'ENTREZID'){
    converted <- as.integer(converted)
  }
  return(converted)
}
# ids = 'Hist1h2af;Hist3h2a'

# orthologs_df <- fread('databases/human_mouse_hcop_fifteen_column.txt') ## for mouse <-> human
# setkeyv(orthologs_df, 'mouse_symbol')
orthologs_conversion <- function(ids, orthologs_df = orthologs_df, fromType = 'mouse_symbol', toType = 'human_entrez_gene'){
  ## splitting ids
  if(length(ids)==1){
    if (nchar(ids)==0){
      
      NA_to_use <- list('integer'= NA_integer_ , 
                        'character' = NA_character_,
                        'logical' =  NA)
      NA_ <- NA_to_use[[typeof(orthologs_df[[toType]])]]
      converted = setNames(data.table(NA_), eval(toType)) ## NA type must match that of the mapped. NA_character_ for chars.
      # converted = setNames(data.table(NA_character_), eval(toType)) ## NA type must match that of the mapped. NA_character_ for chars.
    }else{
      # suppose synonymous ids are separated by either ';', ' ' or ','
      ids = gsub(pattern = ';|,|\ ', x = ids, replacement = ';')
      converted = unique(orthologs_df[unique(tstrsplit(ids, ';')), eval(toType), with = F])
    }
    
  } else{
    converted <- sapply(ids, orthologs_conversion, orthologs_df = orthologs_df, fromType = fromType, toType = toType)
  }
  return(converted)
}

AAcomp_Cor <- function(queryComp, AAComp.dt, essential = F){
  
}



#' convert a data.table to a data.frame with unique rownames
#'
#' @param dt 
#' @param rn.col name of the column in the data table to be used as row nanes. If unspecified, the first column will be used
#'
#' @return
#' @export
#'
#' @examples
dt.to.df <- function(dt, rn.col = NULL, remove.rn = T) {
  if(is.null(rn.col)){
    rn.col <- names(dt)[1]
  }
  if (sum(names(dt) == rn.col)!=1){
    message('row name column matched to multiple columns in data.table')
    stop()
  }
  if (remove.rn){
    df <- as.data.frame(dt[, -which(names(dt)==rn.col),  with = F])
  }else{
    df <- as.data.frame(dt)
  }
  rownames(df) <- dt[, get(rn.col)]
  return(df)
}


dt.to.ExpSet = function(count.dt, count_col, gene_col = 'genes', sample_col = 'sample_ID', remove.0.Var= TRUE, log_ = FALSE,  ntop=Inf)
{
  
  # calculate the variance for each gene
  # select the ntop genes by variance
  var_cutoff <- ifelse(remove.0.Var, 0, -1)
  count.dt <- count.dt[!is.na(get(gene_col))] # remove entries with undefined gene_col
  select <- count.dt[count.dt[,var(get(count_col)), by = gene_col][V1>var_cutoff][order(V1, decreasing = T)[seq_len(min(ntop, .N))]], on = gene_col]
  if(any(select[,duplicated(get(gene_col)),by = sample_col]$V1)) message('duplication detected, using mean aggregation to bin counts/ fpkm from duplicate genes')
  select.df <- select %>% dcast(get(gene_col)~get(sample_col), value.var = count_col,  fun.aggregate = mean) %>% dt.to.df
  if (log_){
    select.df <- log(1+select.df)
  }
  return(select.df)
}

plotPCA.DT = function(ExpSet.df, sample_info, sample_col = 'sample_ID', 
                      group_col="Group.info", returnData=FALSE, biplot = F, scale = T, plot.labels = T)
{
  group.dt <- sample_info[names(ExpSet.df), c(sample_col, group_col), with = F,on = sample_col]
  # calculate the variance for each gene
  # select the ntop genes by variance
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t( ExpSet.df), scale = scale)
  if(biplot){
    p.biplot <- ggbiplot::ggbiplot(pcobj = pca, obs.scale = 1, var.scale = 1, ellipse = F, groups = group_col)
    # p.biplot <- ggbiplot::ggbiplot(pcobj = pca,
    #                                obs.scale = 1, var.scale = 1, ellipse = TRUE, circle = F,
    #                                groups = cut(group.dt[[group_col]], breaks = quantile(group.dt[[group_col]]),
    #                                             labels = c('0-25%', '25-50%', '50-75%', '75-100%'), include.lowest = T), 
    #                                alpha = 0.5)+
    #   scale_colour_discrete(name = 'protein yield') +
    #   theme(legend.direction = 'horizontal', legend.position = 'top') + theme_classic()
    return(p.biplot)
  }
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  # if (!all(group_col %in% names(colData(object)))) {
  #   stop("the argument 'group_col' should specify columns of colData(dds)")
  # }
  
  # add the group_col factors together to create a new grouping factor
  group <- group.dt[[group_col]]
  # assembly the data for the plot
  d <- data.table(PC1=pca$x[,1], PC2=pca$x[,2], group.dt) %>% merge(sample_info[, -group_col, with = F], by = sample_col)
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  p <- ggplot(data=d, aes_string(x="PC1", y="PC2", color=group_col, label = sample_col)) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed()
  if(plot.labels) return( p+ ggrepel::geom_text_repel() + geom_point() ) else return(p + geom_point(alpha = .5))
}


plotTSNE.DT = function(ExpSet.df, sample_info, tsne = NULL, subset.sample = NULL, sample_col = 'sample_ID', 
                       group_col="Group.info", returnData=FALSE, ...)
{
  
  # calculate the variance for each gene
  # select the ntop genes by variance
  # perform a PCA on the data in assay(x) for the selected genes
  if(is.null(tsne)) tsne <- Rtsne::Rtsne(t(ExpSet.df), ...)
  
  if(!is.null(subset.sample)){
    included.ids <- which(names(ExpSet.df)%in%subset.sample)
    ExpSet.df <- ExpSet.df[, included.ids]
  }else{
    included.ids <- seq_along(ExpSet.df)
  }
  group.dt <- sample_info[names(ExpSet.df), c(sample_col, group_col), with = F,on = sample_col]
  # add the group_col factors together to create a new grouping factor
  group <- group.dt[[group_col]]
  # assembly the data for the plot
  # assembly the data for the plot
  d <- data.table(TSNE1=tsne$Y[included.ids,1], TSNE2=tsne$Y[included.ids,2], group.dt) %>%
    merge(sample_info[, -group_col, with = F], by = sample_col)
  if (returnData) {
    # attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data=d, aes_string(x="TSNE1", y="TSNE2", color=group_col, label = sample_col)) + ggrepel::geom_text_repel() + geom_point()
}

plotUMAP.DT = function(ExpSet.df, sample_info, umap = NULL, subset.sample = NULL, sample_col = 'sample_ID', 
                       group_col="Group.info", returnData=FALSE, ...)
{
  
  # calculate the variance for each gene
  # select the ntop genes by variance
  # perform a PCA on the data in assay(x) for the selected genes
  if(is.null(umap)) umap <- umap::umap(t(ExpSet.df), ...)
  
  if(!is.null(subset.sample)){
    included.ids <- which(names(ExpSet.df)%in%subset.sample)
    ExpSet.df <- ExpSet.df[, included.ids]
  }else{
    included.ids <- seq_along(ExpSet.df)
  }
  group.dt <- sample_info[names(ExpSet.df), c(sample_col, group_col), with = F,on = sample_col]
  # add the group_col factors together to create a new grouping factor
  group <- group.dt[[group_col]]
  # assembly the data for the plot
  # assembly the data for the plot
  d <- data.table(UMAP1=umap$layout[included.ids,1], UMAP2=umap$layout[included.ids,2], group.dt) %>%
    merge(sample_info[, -group_col, with = F], by = sample_col)
  if (returnData) {
    # attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data=d, aes_string(x="UMAP1", y="UMAP2", color=group_col, label = sample_col)) + ggrepel::geom_text_repel() + geom_point()
}


parse.localiation <- function(localization.list) {
  localization = c()
  Reliability = c()
  Enhanced = localization.list[['Enhanced']]
  Supported = localization.list[['Supported']]
  Approved  = localization.list[['Approved']]
  Uncertain = localization.list[['Uncertain']]
  
  if(!is.na(Enhanced)){
    localization <- c(localization, unlist(strsplit(Enhanced, ';')))
    Reliability <- c(Reliability, 'Enhanced')
  }else if (!is.na(Supported)){
    localization <- c(localization, unlist(strsplit(Supported, ';')))
    Reliability <- c(Reliability, 'Supported')
  }else if (!is.na(Approved)){
    localization <- c(localization, unlist(strsplit(Approved, ';')))
    Reliability <- c(Reliability, 'Approved')
  }else if (!is.na(Uncertain)){
    localization <- c(localization, unlist(strsplit(gsub(pattern =  ' \\(GO:\\d\\d\\d\\d\\d\\d\\d\\)', x = Uncertain, replacement = ''), ';')))
    Reliability <- c(Reliability, 'Uncertain')
  }
  
  return( list(localization = unique(localization), 
               Reliability = factor(Reliability, levels = rev(c('Enhanced', 'Supported', 'Approved', 'Uncertain')))
  )
  )
  
}

#' due to gene synonyms that occurred inevitably during ortholog conversion, some proteins may have different localizations. This function summarizes and removed multiple localizations, leaving only the most confident localizations.
#'
#' @param localization.list.m 
#'
#' @return
#' @export
#'
#' @examples
parse.localiation.duplicates <- function(localization.list.m) {
  if(nrow(localization.list.m)!=length(unique(localization.list.m$localization))){
    return(localization.list.m[, unique(.SD[Reliability==max(Reliability)]), .SDcols = c('Reliability') ,by = localization])
  }else{
    return(localization.list.m[, .SD, .SDcols = c('localization', 'Reliability')])
  }
}


parse.anne.mathilde.feature <- function(features) {
  feature.split <- strsplit(features, '; |;') %>% unlist %>% gsub(x = ., pattern = "^ |^\\?", replacement = '')
}


#' generate AA composition table for a list of protein sequences in data.table
#' 
#'This function generates AA composition signatures from AA sequences using the seqinr package. Reference AA composition in CHO is 
#'compared to illustrate how far the AA comp of a given protein deviates from the regular CHO cell's.
#'
#' @param seq.col the name of the column containing the AA seq without signal peptide for the proteins
#' @param proteins.dt a data.tble of n rows with seq.col containing all AA sequences
#'
#' @return a data.table of n rows whose columns represent the AA features generated by seqinr package and the Mesia algorithm.
#' @export
#'
#' @examples
get.AA.composition <- function(proteins.dt, seq.col = 'seq', AAcomp = fread('databases/AA_composition.csv')) {
  # get.seq.occur.dt<- function(seq) {
  #   seq.occur.dt <- as.data.table(table(strsplit(toupper(seq), split = ''))) %>% setnames('V1', 'AA.code')
  #   seq.occur.dt[, freq:=N/sum(N)]
  #   return(seq.occur.dt)
  # }
  
  library(seqinr)
  AAcomp <- fread('databases/AA_composition.csv') ## empirical AA compostion of CHO
  AAcomp_cor_dt <- proteins.dt[,setNames(data.table(AAstat(getSequence(get(seq.col)), plot = F)[['Compo']]/ nchar(get(seq.col))),
                                         c('AA.Code', 'AA.Comp')), by =1:nrow(proteins.dt)] %>% 
    merge(AAcomp[, .(AA.Code, AA.comp.CHO, is.essential)],
          all.x = T) %>% 
    .[AA.Code!='*'] %>% 
    .[, AA.comp_CHO_ratio:= AA.Comp/ AA.comp.CHO]
  AAcomp_cor <- merge(AAcomp_cor_dt[, list(cor.all = cor(AA.Comp,AA.comp.CHO)),  by = 'nrow'],
                      AAcomp_cor_dt[is.essential==T, list(cor.essential = cor(AA.Comp,AA.comp.CHO)),
                                    by = 'nrow'], by = 'nrow')
  AAcomp_Mesia <- dcast(AAcomp_cor_dt, formula = nrow~AA.Code, value.var = c('AA.Comp', 'AA.comp_CHO_ratio'))
  AAcomp_final <- merge(AAcomp_cor, AAcomp_Mesia, by = 'nrow') %>% .[,nrow:=NULL]
  
  return(AAcomp_final)
}

get.AA.count <- function(seq, name_prefix = ''){
  if(is.na(seq) | nchar(seq)==0) seq <- 'X'
  AA.count.dt <- 
    seq %>% seqinr::getSequence() %>% seqinr::AAstat(plot = F) %>% .[['Compo']] %>% 
    data.table %>% melt(id.vars = 'V1') %>% dcast(variable~V1, value.var = 'value') %>%
    .[, -c(1:2)] #removing the variable column and *
  if(nchar(name_prefix) == 0) return(AA.count.dt)
  else{
    AA.count.dt %>% setnames(., paste0(name_prefix, '_', names(.)))
  }
}
get.AA.count.m <- function(seq, name_prefix = ''){
  if(is.na(seq) | nchar(seq)==0) return(NULL)
  AA.count.dt <- 
    seq %>% seqinr::getSequence() %>% seqinr::AAstat(plot = F) %>% .[['Compo']] %>% 
    data.table %>% melt(id.vars = 'V1')
  if(nchar(name_prefix) == 0) return(AA.count.dt)
  else{
    AA.count.dt %>% setnames(., paste0(name_prefix, '_', names(.)))
  }
}


#' Sigmoid transformation
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
sigmoid <- function(x) {
  1/(1 + exp(-x))
}

## if you want the so-called 'error function'
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
## (see Abramowitz and Stegun 29.2.29)
## and the so-called 'complementary error function'
erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
## and the inverses
erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
erfcinv <- function (x) qnorm(x/2, lower = FALSE)/sqrt(2)

txi.id.conversion <- function(df,
                              conversion.table = fread('databases/ConversionTable/180405CHO_HUMAN_SYMBOL_ID_MAPPING_all.tsv', na.strings = ''),
                              convert.to = 'humanSymbol', convert.from = 'choSymbol'){
  if(is.character(df)) return(df)
  dt <- df %>% as.data.table(keep.rownames = convert.from)
  ids <- colnames(df)
  dt.transgene <- dt[get(convert.from)%like%'Johan'][, .SD, .SDcols = c(convert.from, ids)]
  dt.endo <- dt[!get(convert.from)%like%'Johan']
  
  dt.endo.human <- conversion.table %>% 
    .[!is.na(get(convert.from)) & !is.na(get(convert.to))] %>% 
    .[dt.endo, on = convert.from] %>% 
    .[,lapply(.SD, mean) ,by = convert.to, .SDcols = sapply(., is.numeric)] %>% 
    .[!is.na(get(convert.to))] # remove unmapped genes
  dt.endo.human[, .SD, .SDcols = c(convert.to, ids)]
  df.combined <- rbind(dt.transgene %>% dt.to.df(), dt.endo.human[, .SD, .SDcols = c(convert.to, ids)] %>% dt.to.df())
  stopifnot(all(colnames(df)==colnames(df.combined)))
  return(as.matrix(df.combined))
}

txi.mean.function <- function(df, duplication.id.dt){
  if(is.character(df)) return(df)
  duplication.id.dt %>% setkeyv(names(.)[1])
  dt.m <- df %>% as.data.table(keep.rownames = T) %>% melt(id.vars = 'rn')
  dt.m.merged <- merge(dt.m, duplication.id.dt, by.x = c('variable'), by.y = names(duplication.id.dt)[1])
  stopifnot(nrow(dt.m) == nrow(dt.m.merged))
  dt.m.merged.mean <- dt.m.merged[, .(mean.value = mean(value)),by =c('rn', names(duplication.id.dt)[-1])]
  return(dt.m.merged.mean %>% dcast(
    as.formula(paste0('rn~',paste0(names(duplication.id.dt)[-1], collapse = '+'))),
    value.var = 'mean.value') %>% dt.to.df() %>% as.matrix)
}

txi.merge <- function(df1, df2, ...){
  if(is.character(df1)){
    if (df1==df2) return(df1) else stop('attributes do not match')}
  dt1 <- df1 %>% as.data.table(keep.rownames = T)
  dt2 <- df2 %>% as.data.table(keep.rownames = T)
  return(merge(dt1, dt2, by = 'rn', ...) %>% dt.to.df() %>% as.matrix())
}

txi.merge.techrep <- function(txi){
  lapply(txi, function(x){
    if(is.character(x)) return(x)
    x %>% melt() %>% as.data.table %>% .[, mean(value), by = .(Var1, Var2)] %>%
      dcast(Var1~Var2, value.var = 'V1') %>% dt.to.df() %>% as.matrix()
  }
  )
}


txi.filter.transgene <- function(df, sample.info = files.list, remove.NA = F) {
  dt <- df %>% as.data.table(keep.rownames = 'choSymbol')
  ids <- colnames(df)
  dt.endo <- dt[!choSymbol%like%'Johan']
  
  dt.transgene.m <- dt[choSymbol%like%'Johan'][, .SD, .SDcols = c('choSymbol', ids)] %>% 
    melt(measure.vars = ids, variable.name = 'sample_ID_unique')
  dt.transgene.renamed <- sample.info[dt.transgene.m, on = 'sample_ID_unique'] %>% 
    .[,.SD[choSymbol%like%sample_ID],by = c('sample_ID_unique', 'sample_ID')] %>%
    .[, choSymbol := 'transgene'] %>% 
    dcast(choSymbol ~sample_ID_unique, value.var = 'value')
  if(remove.NA){
    common.ids <- intersect(names(dt.transgene.renamed), names(dt.endo))
    return(rbind(dt.transgene.renamed[, common.ids, with = F] , dt.endo[, common.ids, with = F] , use.names = T) %>% dt.to.df %>% as.matrix())}
  
  return(rbind(dt.transgene.renamed , dt.endo , use.names = T, fill = T) %>% dt.to.df %>% as.matrix())
  
  
}

txi.filter <- function(txi, row.index) {
  lapply(txi, function(df, row.index){
    if(!is.data.frame(df)) return(df)
    return(df[row.index,])
  }, row.index =  row.index)
}

txi.filter.col <- function(txi, col.index) {
  lapply(txi, function(df, col.index){
    if(!is.data.frame(df) & !is.matrix(df)) return(df)
    return(df[,col.index])
  }, col.index = col.index)
}


pheatmap_0.center <- function(mat, ...) {
  paletteLength <- 50
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  # row.tree <- mat %>% dist %>% hclust %>% dendsort::dendsort()
  # col.tree <- mat %>% t %>% dist %>% hclust %>% dendsort::dendsort()
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks <- c(seq(min(mat), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(mat)/paletteLength, max(mat), length.out=floor(paletteLength/2)))
  pheatmap::pheatmap(mat, color=myColor, breaks=myBreaks, ...)
}

#' Count the number of unique entries in the vector
#'
#' @param variables 
#'
#' @return
#' @export
#'
#' @examples
n.unique <- function(variables) {
  variables %>% unique() %>% length
}



quantile_breaks <- function(xs, n = 20) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

pv.trans <- function(psuedo.pv=1e-8){
  return(scales::trans_new('pv.trans',
                           transform = 
                             function(x){ 
                               x[x<psuedo.pv] <- psuedo.pv
                               -log10(x) }, 
                           inverse = function(y){10**(-y)}))
}


#' Generate a sequence over the range of a vector
#'
#' @param x A numeric vector
#' @param n,by Specify the output sequence either by supplying the
#'   length of the sequence with `n`, or the spacing between value
#'   with `by`. Specifying both is an error.
#'
#'   I recommend that you name these arguments in order to make it clear to
#'   the reader.
#' @param pretty If `TRUE`, will generate a pretty sequence. If `n`
#'   is supplied, this will use [pretty()] instead of
#'   [seq()]. If `by` is supplied, it will round the first
#'   value to a multiple of `by`.
#' @param trim Optionally, trim values off the tails.
#'   `trim / 2 * length(x)` values are removed from each tail.
#' @param expand Optionally, expand the range by `expand * (1 + range(x)`
#'   (computed after trimming).
#' @export
#' @examples
#' x <- rcauchy(100)
#' seq_range(x, n = 10)
#' seq_range(x, n = 10, trim = 0.1)
#' seq_range(x, by = 1, trim = 0.1)
#'
#' # Make pretty sequences
#' y <- runif(100)
#' seq_range(y, n = 10)
#' seq_range(y, n = 10, pretty = TRUE)
#' seq_range(y, n = 10, expand = 0.5, pretty = TRUE)
#'
#' seq_range(y, by = 0.1)
#' seq_range(y, by = 0.1, pretty = TRUE)
seq_range <- function(x, n, by, trim = NULL, expand = NULL, pretty = FALSE) {
  if (!missing(n) && !missing(by)) {
    stop("May only specify one of `n` and `by`", call. = FALSE)
  }
  
  if (!is.null(trim)) {
    rng <- stats::quantile(x, c(trim / 2, 1 - trim / 2), na.rm = TRUE)
  } else {
    rng <- range(x, na.rm = TRUE)
  }
  
  if (!is.null(expand)) {
    rng <- rng + c(-expand / 2, expand / 2) * (rng[2] - rng[1])
  }
  
  if (missing(by)) {
    if (pretty) {
      pretty(rng, n)
    } else {
      seq(rng[1], rng[2], length.out = n)
    }
  } else {
    if (pretty) {
      rng[1] <- floor(rng[1] / by) * by
      rng[2] <- ceiling(rng[2] / by) * by
    }
    seq(rng[1], rng[2], by = by)
  }
}

fisher.method <- function(...){
  P.vec <- unlist(list(...)) %>% na.omit()
  if(length(P.vec)==0){ 
    stop()
  }
  if(length(P.vec)==1){
    return(P.vec)
  }
  if(length(P.vec)>1){
    return(pchisq(-2*sum(log(P.vec+10^-10)),df=2*length(P.vec),lower.tail = F))
  }
}


ggsave2 <- function(file.name, path.name = 'Figures_CHOHEK/', height, width, plot = last_plot()) {
  ggsave(plot = plot, paste0(path.name, file.name, '.png'),height = height, width = width, dpi = 150)
  ggsave(plot = plot, paste0(path.name, file.name, '.pdf'),height = height, width = width, device = cairo_pdf)
}


#' Geometric Mean
#'
#' @param x 
#' @param na.rm 
#' @param zero.propagate 
#'
#' @return
#' @export
#'
#' @examples
gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}