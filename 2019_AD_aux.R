library(enrichR)
##load genesets
ChEA.gmt <- fgsea::gmtPathways('databases/enrichrTFs/ChEA_2016.txt') %>% 
  .[grepl('_human$', names(.) , ignore.case = T)] %>% ## use human genesets only
  `names<-`(tstrsplit(names(.), '_')[[1]]) ## remove cell line information

ENCODE.gmt <- fgsea::gmtPathways('databases/enrichrTFs/ENCODE_TF_ChIP-seq_2015.txt') %>% 
  .[grepl('_hg19$', names(.) , ignore.case = T)] %>% ## use human genesets only
  `names<-`(tstrsplit(names(.), '_')[[1]]) ## remove cell line information
# litChip.gmt <- fgsea::gmtPathways('databases/enrichrTFs/Literature_ChIP-seq.gmt.txt')
# litChip.gmt <- litChip.gmt[grepl('_human$', names(litChip.gmt) , ignore.case = T)] ## use human genesets only
# ReMapChip.gmt <- fgsea::gmtPathways('databases/enrichrTFs/ReMap_ChIP-seq.gmt.txt')

eqtl.gmt <- fread('databases/enrichrTFs/eqtl.genesets.csv') %>% 
  .[, .SD[.N>5], by=risk_gene] %>% 
  split(by='risk_gene') %>% lapply(function(dt) dt$eqtl.gene)
enhancer.gwas.gmt <- fread('databases/2019_AD_GWAS/natGen_ALL_GWAS_enhancer.csv')[coverage>.1] %>% 
  split(by='SYMBOL') %>% lapply(function(dt) dt$attribute)
gmt.lists <- list(ChEA = ChEA.gmt,
                  ENCODE = ENCODE.gmt
                  # GeneHancer = enhancer.gwas.gmt
                  # eqtl = eqtl.gmt ## TODO
                  # LitChIP = litChip.gmt,
                  # ReMapCHIP = ReMapChip.gmt
)
combined.gmt.dt <- lapply(gmt.lists, function(gmt){lapply(gmt, function(geneList) data.table(target = geneList)) %>% rbindlist(idcol='regulator')}) %>% rbindlist(idcol='gmt.name')
combined.gmt <- combined.gmt.dt[, .(target=unique(target)),by=regulator] %>% split(by='regulator') %>% lapply(function(dt) dt$target)
# all.TFs <- unique(unlist(lapply(gmt.lists, names)))

## adapted from fisher.test; Override conf.int calculation when alternative =='g', conf.int is calculated using the default formulation for alternative=='two.sided'
fisher.test_confIntOverride <- 
  function (x, y = NULL, workspace = 2e+05, hybrid = FALSE, hybridPars = c(expect = 5, 
                                                                           percent = 80, Emin = 1), control = list(), or = 1, alternative = "two.sided", 
            conf.int = TRUE, conf.level = 0.95, simulate.p.value = FALSE, 
            B = 2000) 
  {
    DNAME <- deparse(substitute(x))
    METHOD <- "Fisher's Exact Test for Count Data"
    if (is.data.frame(x)) 
      x <- as.matrix(x)
    if (is.matrix(x)) {
      if (any(dim(x) < 2L)) 
        stop("'x' must have at least 2 rows and columns")
      if (!is.numeric(x) || any(x < 0) || anyNA(x)) 
        stop("all entries of 'x' must be nonnegative and finite")
      if (!is.integer(x)) {
        xo <- x
        x <- round(x)
        if (any(x > .Machine$integer.max)) 
          stop("'x' has entries too large to be integer")
        if (!identical(TRUE, (ax <- all.equal(xo, x)))) 
          warning(gettextf("'x' has been rounded to integer: %s", 
                           ax), domain = NA)
        storage.mode(x) <- "integer"
      }
    }
    else {
      if (is.null(y)) 
        stop("if 'x' is not a matrix, 'y' must be given")
      if (length(x) != length(y)) 
        stop("'x' and 'y' must have the same length")
      DNAME <- paste(DNAME, "and", deparse1(substitute(y)))
      OK <- complete.cases(x, y)
      x <- as.factor(x[OK])
      y <- as.factor(y[OK])
      if ((nlevels(x) < 2L) || (nlevels(y) < 2L)) 
        stop("'x' and 'y' must have at least 2 levels")
      x <- table(x, y)
    }
    con <- list(mult = 30)
    con[names(control)] <- control
    if ((mult <- as.integer(con$mult)) < 2) 
      stop("'mult' must be integer >= 2, typically = 30")
    nr <- nrow(x)
    nc <- ncol(x)
    have.2x2 <- (nr == 2) && (nc == 2)
    if (have.2x2) {
      alternative <- char.expand(alternative, c("two.sided", 
                                                "less", "greater"))
      if (length(alternative) > 1L || is.na(alternative)) 
        stop("alternative must be \"two.sided\", \"less\" or \"greater\"")
      if (!((length(conf.level) == 1L) && is.finite(conf.level) && 
            (conf.level > 0) && (conf.level < 1))) 
        stop("'conf.level' must be a single number between 0 and 1")
      if (!missing(or) && (length(or) > 1L || is.na(or) || 
                           or < 0)) 
        stop("'or' must be a single number between 0 and Inf")
    }
    PVAL <- NULL
    if (!have.2x2) {
      if (simulate.p.value) {
        sr <- rowSums(x)
        sc <- colSums(x)
        x <- x[sr > 0, sc > 0, drop = FALSE]
        nr <- as.integer(nrow(x))
        nc <- as.integer(ncol(x))
        if (is.na(nr) || is.na(nc) || is.na(nr * nc)) 
          stop("invalid nrow(x) or ncol(x)", domain = NA)
        if (nr <= 1L) 
          stop("need 2 or more non-zero row marginals")
        if (nc <= 1L) 
          stop("need 2 or more non-zero column marginals")
        METHOD <- paste(METHOD, "with simulated p-value\n\t (based on", 
                        B, "replicates)")
        STATISTIC <- -sum(lfactorial(x))
        tmp <- .Call(C_Fisher_sim, rowSums(x), colSums(x), 
                     B)
        almost.1 <- 1 + 64 * .Machine$double.eps
        PVAL <- (1 + sum(tmp <= STATISTIC/almost.1))/(B + 
                                                        1)
      }
      else if (hybrid) {
        if (!is.null(nhP <- names(hybridPars)) && !identical(nhP, 
                                                             c("expect", "percent", "Emin"))) 
          stop("names(hybridPars) should be NULL or be identical to the default's")
        stopifnot(is.double(hypp <- as.double(hybridPars)), 
                  length(hypp) == 3L, hypp[1] > 0, hypp[3] >= 
                    0, 0 <= hypp[2], hypp[2] <= 100)
        PVAL <- .Call(C_Fexact, x, hypp, workspace, mult)
        METHOD <- paste(METHOD, sprintf("hybrid using asym.chisq. iff (exp=%g, perc=%g, Emin=%g)", 
                                        hypp[1], hypp[2], hypp[3]))
      }
      else {
        PVAL <- .Call(C_Fexact, x, c(-1, 100, 0), workspace, 
                      mult)
      }
      RVAL <- list(p.value = max(0, min(1, PVAL)))
    }
    else {
      if (hybrid) 
        warning("'hybrid' is ignored for a 2 x 2 table")
      m <- sum(x[, 1L])
      n <- sum(x[, 2L])
      k <- sum(x[1L, ])
      x <- x[1L, 1L]
      lo <- max(0L, k - n)
      hi <- min(k, m)
      NVAL <- c(`odds ratio` = or)
      support <- lo:hi
      logdc <- dhyper(support, m, n, k, log = TRUE)
      dnhyper <- function(ncp) {
        d <- logdc + log(ncp) * support
        d <- exp(d - max(d))
        d/sum(d)
      }
      mnhyper <- function(ncp) {
        if (ncp == 0) 
          return(lo)
        if (ncp == Inf) 
          return(hi)
        sum(support * dnhyper(ncp))
      }
      pnhyper <- function(q, ncp = 1, upper.tail = FALSE) {
        if (ncp == 1) {
          return(if (upper.tail) phyper(x - 1, m, n, k, 
                                        lower.tail = FALSE) else phyper(x, m, n, k))
        }
        if (ncp == 0) {
          return(as.numeric(if (upper.tail) q <= lo else q >= 
                              lo))
        }
        if (ncp == Inf) {
          return(as.numeric(if (upper.tail) q <= hi else q >= 
                              hi))
        }
        sum(dnhyper(ncp)[if (upper.tail) support >= q else support <= 
                           q])
      }
      if (is.null(PVAL)) {
        PVAL <- switch(alternative, less = pnhyper(x, or), 
                       greater = pnhyper(x, or, upper.tail = TRUE), 
                       two.sided = {
                         if (or == 0) as.numeric(x == lo) else if (or == 
                                                                   Inf) as.numeric(x == hi) else {
                                                                     relErr <- 1 + 10^(-7)
                                                                     d <- dnhyper(or)
                                                                     sum(d[d <= d[x - lo + 1] * relErr])
                                                                   }
                       })
        RVAL <- list(p.value = PVAL)
      }
      mle <- function(x) {
        if (x == lo) 
          return(0)
        if (x == hi) 
          return(Inf)
        mu <- mnhyper(1)
        if (mu > x) 
          uniroot(function(t) mnhyper(t) - x, c(0, 1))$root
        else if (mu < x) 
          1/uniroot(function(t) mnhyper(1/t) - x, c(.Machine$double.eps, 
                                                    1))$root
        else 1
      }
      ESTIMATE <- c(`odds ratio` = mle(x))
      if (conf.int) {
        ncp.U <- function(x, alpha) {
          if (x == hi) 
            return(Inf)
          p <- pnhyper(x, 1)
          if (p < alpha) 
            uniroot(function(t) pnhyper(x, t) - alpha, 
                    c(0, 1))$root
          else if (p > alpha) 
            1/uniroot(function(t) pnhyper(x, 1/t) - alpha, 
                      c(.Machine$double.eps, 1))$root
          else 1
        }
        ncp.L <- function(x, alpha) {
          if (x == lo) 
            return(0)
          p <- pnhyper(x, 1, upper.tail = TRUE)
          if (p > alpha) 
            uniroot(function(t) pnhyper(x, t, upper.tail = TRUE) - 
                      alpha, c(0, 1))$root
          else if (p < alpha) 
            1/uniroot(function(t) pnhyper(x, 1/t, upper.tail = TRUE) - 
                        alpha, c(.Machine$double.eps, 1))$root
          else 1
        }
        CINT <- switch(alternative, less = c(0, ncp.U(x, 
                                                      1 - conf.level)), 
                       greater = {
                         alpha <- (1 - conf.level)/2
                         c(ncp.L(x, alpha), ncp.U(x, alpha))
                       }, 
                       two.sided = {
                         alpha <- (1 - conf.level)/2
                         c(ncp.L(x, alpha), ncp.U(x, alpha))
                       })
        attr(CINT, "conf.level") <- conf.level
      }
      RVAL <- c(RVAL, list(conf.int = if (conf.int) CINT, 
                           estimate = ESTIMATE, null.value = NVAL))
    }
    structure(c(RVAL, alternative = alternative, method = METHOD, 
                data.name = DNAME), class = "htest")
  }




#' Enrichment Query
#'
#' @param genes 
#' @param server use chea3 for online enrichment
#' @param encode 
#'
#' @return
#' @export
#'
#' @examples
enrichr.wrapper <- function(genes, server='chea3', fisher.gmt.lists= NULL, background=NULL, subset_background= T, mc.cores=1) {
  
  if(server=='chea3'){
    encode = "json"
    url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
    payload = list(query_name = "myQuery", gene_set = genes)
    
    #POST to ChEA3 server
    response = httr::POST(url = url, body = payload, encode = encode)
    json = httr::content(response, "text")
    #results as list of R dataframes
    results = jsonlite::fromJSON(json) 
    results.dt <- results %>% rbindlist(fill=T, idcol='library')					
    results.dt[, c('Rank', 'Score', 'Scaled Rank', 'Intersect', 'Set length', 'FET p-value', 'FDR', 'Odds Ratio'):=
                 .(as.numeric(Rank), as.numeric(Score), as.numeric(`Scaled Rank`), as.numeric(Intersect),
                   as.numeric(`Set length`), as.numeric(`FET p-value`), as.numeric(FDR), as.numeric(`Odds Ratio`)
                 )] 
  }else if (server=='enrichr'){
    capture.output(results.dt <- enrichR::enrichr(as.character(genes), c('ChEA_2016', 'ENCODE_TF_ChIP-seq_2015')) %>% rbindlist(idcol='library') %>% .[, Genes:=NULL])
    results.dt[, TF:=tstrsplit(Term, '_')[1]] %>% setnames(c('P.value', 'Adjusted.P.value'), c('FET p-value', 'FDR'))
  }else if (server== 'fisher.exact'){
    #load GMTs
    if(is.null(background)| is.null(fisher.gmt.lists)) 
      background <- c(unlist(ENCODE.gmt, use.names = F), unlist(ChEA.gmt, use.names = F))%>% unique
    
    
    enrich.genelist <- function(gl, gmt.background){
      results.dt.singleTF <- fisher.test( 
        table(factor(gmt.background%in%gl, levels = c(T,F)),
              factor(gmt.background%in%genes, levels = c(T,F))),
        alternative = 'g', conf.int = F)
      data.table(`FET p-value`=results.dt.singleTF$p.value, 
                 `Odds Ratio`=results.dt.singleTF$estimate)
    }
    get.gme.enrich <- function(gmt, gmt.background, subset_background) {
      if(!is.null(gmt.background)){
        if(subset_background){
          gmt.background=intersect(unlist(gmt, use.names = F),
                                   gmt.background)
        }
      }else gmt.background <- unique(unlist(gmt, use.names = F))
      results.dt.single.gmt <- parallel::mclapply( gmt, enrich.genelist,
                                                   gmt.background = gmt.background, 
                                                   mc.cores = mc.cores) %>% rbindlist(idcol = 'TF')
      return(results.dt.single.gmt)
    }
    
    results.dt <- lapply(fisher.gmt.lists, get.gme.enrich, gmt.background=background, subset_background = subset_background) %>% 
      rbindlist(idcol='gmt.source')
    # results.dt[, TF:=tstrsplit(TF, '_')[1]]
  }
  # if(!is.null(background)) results.dt <- results.dt[TF%in%background]
  return(results.dt)
}

#' Title
#'
#' @param enrich.results 
#' @param TF.col 
#' @param pv.col 
#' @param p.adj.method set to NA to disable
#'
#' @return
#' @export
#'
#' @examples
aggr.enrich.table <- function(enrich.results, TF.col='TF', pv.col='FET p-value', p.adj.method='fdr',
                              odds.ratio.col='Odds Ratio', filter.background=NULL){
  if(!is.null(filter.background)){
    enrich.results <- enrich.results[get(TF.col)%in%filter.background]
  }
  
  enrich.results.aggr <- ## pvalue aggregation with fisher's method
    enrich.results[!is.na(get(pv.col)), .(TF.pvalue = pchisq(-2*sum(log(get(pv.col) + 10^-10)),
                                                             df=2*.N,lower.tail = F),
                                          TF.OR = mean(get(odds.ratio.col))), by=TF.col]
  if(!is.na(p.adj.method)){
    enrich.results.aggr[, TF.pvalue:=p.adjust(TF.pvalue, method = p.adj.method)]
  }
  return(enrich.results.aggr)
}

quantify.enrichment <- function(enrich.results, AD_genes, 
                                # method='fisher.exact', #GSEA
                                TF.id.col = 'TF', TF.pv.col='TF.pvalue', TF.sig=.05,
                                background=NULL,
                                fisher.alternative='g'){
  
  ## diagnostic mosaic plot
  # data.table(risk_genes = ifelse(background%in%AD_genes,  'is_risk_genes', 'non_risk_genes'),
  #            TF.signif = ifelse(background%in%enrich.results[get(TF.pv.col)<TF.sig, get(TF.id.col)], 'is_TF.signif', 'non_TF.signif'),
  #            in.TF = ifelse(background%in%enrich.results$TF, 'is_in.TF', 'non_in.TF')) %>%
  #   ggplot() + ggmosaic::geom_mosaic(aes(x=ggmosaic::product(risk_genes), fill=TF.signif, conds=ggmosaic::product(in.TF)))
  
  
  if(is.null(background)) background <- enrich.results$TF
  sig.enrich <- enrich.results[get(TF.pv.col)<TF.sig, get(TF.id.col)]
  
  if(length(dim(AD_genes))==2){ ## hierarchical enrichment
    non.AD.genes <- setdiff(background, AD_genes$Gene)
    all.TF.locus <- AD_genes[Gene%in%background, unique(Locus)]
    sig.enrich <- c(setdiff(sig.enrich, AD_genes$Gene) , AD_genes[Gene%in%sig.enrich, unique(Locus)])
    background <- c(all.TF.locus, non.AD.genes)
    
    secondary.enrichment.table <- 
      table(TF.signif = factor(background%in%sig.enrich, levels = c(T,F)),
            risk_genes = factor(background%in%AD_genes$Locus, levels = c(T,F))) 
    hit.TF <- Reduce(intersect, list(background, sig.enrich, all.TF.locus))
    
  }else{
    secondary.enrichment.table <- table(TF.signif = factor(background%in%sig.enrich, levels = c(T,F)),
                                        risk_genes = factor(background%in%AD_genes, levels = c(T,F))
    )
    hit.TF <- Reduce(intersect, list(background, sig.enrich, AD_genes))
  }
  
  fisher.res <- tryCatch(
    {
      fisher.test_confIntOverride(secondary.enrichment.table, alternative = fisher.alternative, conf.level = .5)
      
    }, error=function(e){
      print(secondary.enrichment.table)
      return(list(estimate=NaN, 
                  p.value=NaN,
                  l.confInt=NaN,
                  r.confInt=NaN))
    }
  )
  
  return(list(OR= fisher.res$estimate,
              p.value=fisher.res$p.value,
              hit.TF = hit.TF,
              l.confInt=fisher.res$conf.int[1],
              r.confInt=fisher.res$conf.int[2])
  )
  # }else if(method=='GSEA'){ ## slow permutation
  #   pathway <- list('AD_risk_genes'=AD_risk_genes)
  #   stats=enrich.results[, .(get(TF.id.col), get(TF.pv.col))] %>% tibble::deframe()
  #   fgsea::fgsea(pathway, stats)
  # }
}


#' Title
#'
#' @param test.dt 
#' @param n.cores 
#' @param cumulative.t.statistic 
#' @param enrich 'targets' for direct target enrichment; 'TF' for enrichment of transcription factors as specified; 'TF_perm' for empirial hierarchical enrichment analysis
#' @param hyper.lookup 
#' @param fisher.alternative 
#' @param AD_genes 
#' @param TF.enrich.p.adj.method 
#' @param TF.enrich.sig 
#' @param enrich.server 
#' @param background 
#' @param secondary.background 
#' @param gwas.cutoff 
#' @param hyper.test.term 
#' @param match.term 
#' @param return.remaining.n 
#' @param x.col.name 
#' @param plot 
#' @param log.p 
#'
#' @return
#' @export
#'
#' @examples
enrich.step.1d <- function(test.dt, 
                           n.cores=1, cumulative.t.statistic = F, gmt.lists=gmt.lists,
                           ## begin TF enrichment parameters##############
                           enrich='targets', hyper.lookup=NULL, fisher.alternative='g', AD_genes=NULL, TF.enrich.p.adj.method = NA,
                           TF_nperm = 200,
                           TF.enrich.sig = .05, enrich.server= 'chea3', background=NULL, secondary.background=NULL,
                           ## end TF enrichment parameters################
                           # gwas.pv.colname = NULL, #'P.raw',
                           gwas.cutoff = 1e-5,
                           hyper.test.term, match.term, max_size= Inf,
                           return.remaining.n=F,
                           x.col.name = 'p_k_mean',
                           plot=T, log.p=T, gsea.dt=NULL) {
  
  # if(is.na(secondary.background)) secondary.background <- background
  ## target enrichment only
  # if(enrich=='targets'){
  #   stopifnot(!is.null(gwas.pv.colname))
  #   m.hyper <- test.dt[get(gwas.pv.colname) <= gwas.cutoff, .N]
  #   n.hyper <- test.dt[get(gwas.pv.colname) > gwas.cutoff, .N]
  # }
  
  test.dt[, x.rank:=rank(get(x.col.name))]
  if(cumulative.t.statistic){
    lfc.var <- test.dt[, var(t.statistic)]
  }
  
  if(enrich=='gmtTarget'){
    if(is.data.table(gmt.lists)){
      targets <- gmt.lists[padj<=TF.enrich.sig, unique(symbol)]
    }else{
      stop('must provide enrichment targets via gmt.list')
    }
    if(!is.null(AD_genes)) stop('AD genes are not used during target gene enrichment')
    
    # hyper.lookup <- data.table(k=numeric(), q=numeric(), hyper.pv=numeric()) %>% setkey(k,q)
  }
  
  if(is.null(hyper.lookup)){
    hyper.lookup <- pbapply::pblapply(#parallel::mclapply(
      1:nrow(test.dt),
      function(single.i){
        current.xrank = test.dt[single.i,x.rank]
        selected.geneSymbol <- test.dt[x.rank>=current.xrank, geneSymbol] %>% sort
        if(cumulative.t.statistic){
          ret.dt <- data.table(
            k=length(selected.geneSymbol),
            t.stat.mean = test.dt[x.rank>=current.xrank, mean(t.statistic)],
            t.stat.sd = test.dt[x.rank>=current.xrank, lfc.var/.N]**.5
          )
        }else if(enrich == 'gmtTarget'){
          if(is.numeric(background)){
            background.gmtTarget <- gmt.lists[padj<=background, unique(symbol)]
          }else if (is.character(background)){
            background.gmtTarget <- background
          }else{
            background.gmtTarget <- gmt.lists[, unique(symbol)]
          }
          
          # filtered.gmt.lists <- lapply(1, function(padj.cutoff){ gmt.lists[padj<=padj.cutoff]})
          # background <- c(unlist(ENCODE.gmt, use.names = F), unlist(ChEA.gmt, use.names = F))%>% unique
          # enrichr.wrapper(selected.geneSymbol, server = 'fisher.exact', fisher.gmt.lists = list(TF = targets))
          enrich.res <- enrichr.wrapper(selected.geneSymbol, server = 'fisher.exact',
                                        fisher.gmt.lists = list('differential_epigenome' = list('pv_cutoff_0.5' = targets)),
                                        background = union(background.gmtTarget, test.dt$geneSymbol),
                                        subset_background = F)
          ret.dt <- data.table(k=length(selected.geneSymbol),
                               sorted.TF=list(selected.geneSymbol),
                               OR = enrich.res$`Odds Ratio`,
                               hyper.pv = enrich.res$`FET p-value`,
                               TF.enrich.sig = TF.enrich.sig)
          
          
        }else if(enrich%in%c('GSEA','gsea')){
          pathway.dt <- data.table(sel.pathway = list(selected.geneSymbol), k = length(selected.geneSymbol))
          return(pathway.dt)
        }else if (enrich=='TF'){
          enrich.res <- enrichr.wrapper(selected.geneSymbol, server = enrich.server, fisher.gmt.lists = gmt.lists,
                                        background = background, mc.cores = 1) %>% 
            aggr.enrich.table(p.adj.method = TF.enrich.p.adj.method, filter.background = background)
          
          fisher.enrich <- quantify.enrichment(enrich.res, AD_genes= AD_genes, fisher.alternative = fisher.alternative, 
                                               TF.sig = TF.enrich.sig, background = secondary.background)
          
          ret.dt <- data.table(k=length(selected.geneSymbol),
                               sorted.TF=list(selected.geneSymbol),
                               OR = fisher.enrich$OR,
                               hyper.pv=fisher.enrich$p.value,
                               l.confInt = fisher.enrich$l.confInt,
                               r.confInt = fisher.enrich$r.confInt,
                               hit.TF = list(fisher.enrich$hit.TF))
          return(ret.dt)
        }else if (enrich=='TF_perm'){
          ### get AD focus genes (AD risk genes that are regulators)
          rAD <- if(length(dim(AD_genes))==2){AD_genes[, unique(Gene)] } else{AD_genes}## hierarchical enrichment
          rAD.TF <- intersect(unlist(lapply(gmt.lists, names), use.names = F), rAD)
          rAD.TF.gmts <- lapply(gmt.lists, function(g.list){
            g.list[which(names(g.list)%in%rAD)]
          })
          ### calculate enrichment odds ratio for all AD focus genes
          rAD.TF.enrich.dt <- 
            enrichr.wrapper(selected.geneSymbol, server=enrich.server, fisher.gmt.lists = rAD.TF.gmts, background = background) %>% 
            aggr.enrich.table(p.adj.method = TF.enrich.p.adj.method) %>% setkey(TF)
          
          perm.rAD.TF.enrich.dt <- lapply(1:TF_nperm, function(n.iter){
            selected.geneSymbol.perm <- test.dt[, sample(geneSymbol,length(selected.geneSymbol))]
            enrichr.wrapper(selected.geneSymbol.perm, server=enrich.server, fisher.gmt.lists = rAD.TF.gmts, background = background) %>% 
              aggr.enrich.table(p.adj.method = TF.enrich.p.adj.method)
          }
          ) %>% rbindlist(idcol='perm.id')
          ### summarize empirical distribution
          return(
            perm.rAD.TF.enrich.dt[, .(
              k=length(selected.geneSymbol),
              # sorted.TF = list(select)
              perm.quantile = ecdf(TF.OR)(rAD.TF.enrich.dt[TF, TF.OR]),
              null.05 = quantile(TF.OR, .05),
              null.25 = quantile(TF.OR, .25),
              null.50 = quantile(TF.OR, .50),
              null.75 = quantile(TF.OR, .75),
              null.95 = quantile(TF.OR, .95)), by=TF] %>% 
              merge(rAD.TF.enrich.dt, by='TF')
          )
        }else stop('valid enrich options: targets/ TF')
      }, cl=n.cores) %>% rbindlist()
    # hyper.lookup <- data.table(k=numeric(), sorted.TF=list(), hyper.pv=numeric())
  }
  if(enrich%in%c('GSEA', 'gsea')){
    gsea.ranks <- gsea.dt[padj<=TF.enrich.sig, .(symbol, lfc)] %>% tibble::deframe()
    # gsea.res <- hyper.lookup[, fgsea::fgsea(`names<-`(sel.pathway, 'sel.pathway'), stats= gsea.ranks, nperm = 5000),by = k]
    hyper.lookup %>% setkey(k)
    gsea.res <- hyper.lookup[k<max_size, {
      if(k%%100==0){message(sprintf('processing %s', k))}
      suppressWarnings(
      # fgsea::fgseaMultilevel(`names<-`(sel.pathway, 'sel.pathway'), stats= gsea.ranks))
      fgsea::fgsea(`names<-`(sel.pathway, 'sel.pathway'), stats= gsea.ranks, nperm = 100))
      },by = k]
    return(gsea.res)
  }else{
    return(hyper.lookup)
  }
  
  
}


#' Sensitivity analysis for 1-way (target) or 2-way enrichment (TF)
#'
#' @param test.dt target data.table for 1-way (target) or 2-way enrichment (TF) analysis
#' @param x.steps 
#' @param y.steps 
#' @param length.out 
#' @param enrich "TF" for for 2-way analysis (the regulatory targets for the reamining genes given the current target data.table cutoff are inferred by ChIP-X Enrichment Analysis. ); 
#' "target" for 1-way enrichment. 
#' @param hyper.lookup enrichment lookup table as returned by previous executions on coarser x/ysteps to speed things up. Optional
#' @param fisher.alternative 
#' @param AD_genes AD genes in list format
#' @param gwas.pv.colname 
#' @param gwas.cutoff 
#' @param hyper.test.term 
#' @param match.term 
#' @param return.remaining.n 
#' @param x.col.name 
#' @param y.col.name 
#' @param plot 
#' @param log.p 
#' @param TF.enrich.p.adj.method for TF enrichment analysis only
#' @param TF.enrich.sig TF FET enrichment cutoff. Only the remaining TFs will be tested for AD_genes enrichment. 
#'
#' @return
#' @export
#'
#' @examples
enrich.step <- function(test.dt, 
                        x.steps=NULL, y.steps=NULL, length.out=NULL, return.hyper.lookup.only=F,
                        ## begin TF enrichment parameters##############
                        enrich='targets',hyper.lookup=NULL, fisher.alternative='g', AD_genes=NULL, TF.enrich.p.adj.method = NA,
                        TF.enrich.sig = .05, enrich.server= 'chea3', background=NULL, secondary.background=NULL,
                        ## end TF enrichment parameters################
                        gwas.pv.colname = 'P.raw',
                        gwas.cutoff = 1e-5,
                        hyper.test.term, match.term, 
                        return.remaining.n=F,
                        x.col.name = 'p_k_mean',
                        y.col.name ='p_k_sd', 
                        plot=T, log.p=T) {
  
  
  if(is.null(x.steps))  x.steps <- seq(range(test.dt[[x.col.name]])[1], range(test.dt[[x.col.name]])[2], length.out=length.out)
  if(is.null(y.steps)) y.steps <- seq(range(test.dt[[y.col.name]])[1], range(test.dt[[y.col.name]])[2], length.out=length.out)
  
  # if(is.na(secondary.background)) secondary.background <- background
  
  ## target enrichment only
  if(enrich=='targets'){
    m.hyper <- test.dt[get(gwas.pv.colname) <= gwas.cutoff, .N]
    n.hyper <- test.dt[get(gwas.pv.colname) > gwas.cutoff, .N]
  }
  
  test.dt[, x.rank:=rank(get(x.col.name))]
  test.dt[, y.rank:=rank(get(y.col.name))]
  
  hyper.pv.m <- matrix(nrow = length(y.steps), ncol = length(x.steps),
                       dimnames = list(ifelse(seq_along(y.steps)%%10==0, formatC(y.steps, format = 'e', digits = 2),''),
                                       ifelse(seq_along(x.steps)%%10==0, formatC(x.steps, format = 'e', digits = 2),'')))
  remaining.m <- copy(hyper.pv.m)
  if(is.null(hyper.lookup)){
    if(enrich=='targets'){
      if(!is.null(AD_genes)) stop('AD genes are not used during target gene enrichment')
      hyper.lookup <- data.table(k=numeric(), q=numeric(), hyper.pv=numeric()) %>% setkey(k,q)
    }else if(enrich=='TF'){
      hyper.lookup <- parallel::mclapply(1:nrow(test.dt),
                                         function(single.i){
                                           current.xrank = test.dt[single.i,x.rank]
                                           current.yrank = test.dt[single.i,y.rank]
                                           selected.geneSymbol <- test.dt[x.rank>=current.xrank & y.rank>=current.yrank, geneSymbol] %>% sort
                                           enrich.res <- enrichr.wrapper(selected.geneSymbol, server = enrich.server, background = background, mc.cores = 1) %>% 
                                             aggr.enrich.table(p.adj.method = TF.enrich.p.adj.method)
                                           
                                           hyper.pv <- quantify.enrichment(enrich.res, AD_genes= AD_genes, fisher.alternative = fisher.alternative, 
                                                                           TF.sig = TF.enrich.sig, background = secondary.background)$p.value
                                           data.table(k=length(selected.geneSymbol),
                                                      sorted.TF=list(selected.geneSymbol),
                                                      hyper.pv=hyper.pv)
                                           
                                         }, mc.cores=16
      ) %>% rbindlist()
      # hyper.lookup <- data.table(k=numeric(), sorted.TF=list(), hyper.pv=numeric())
    }else stop('valid enrich options: taargets/ TF')
  }
  
  if(return.hyper.lookup.only) return(hyper.lookup)
  
  
  for(i in seq_along(y.steps)){
    if((100*i/length(y.steps))%%5==0){
      message(sprintf('running... %s percent', (100*i/length(y.steps))))
    }
    for(j in seq_along(x.steps)){
      select.test.dt <- test.dt[get(y.col.name)>=y.steps[i] & get(x.col.name)>= x.steps[j]]
      k.hyper <- select.test.dt[, .N]
      
      if(enrich=='targets'){
        remaining.m[i,j] <- k.hyper
        q.hyper <- select.test.dt[get(gwas.pv.colname) <= gwas.cutoff, .N]
        
        lookup.res <- hyper.lookup[k==k.hyper&q==q.hyper]
        if(nrow(lookup.res)==0){ ## no match, generate hyper.pv
          hyper.pv <- phyper(q.hyper, m.hyper, n.hyper, k.hyper, lower.tail = F)
          hyper.lookup %<>% rbind(data.table(k=k.hyper, q=q.hyper, hyper.pv=hyper.pv)) 
        }else{
          hyper.pv <- lookup.res$hyper.pv
        }
      }else if(enrich=='TF'){
        if(k.hyper<3){
          hyper.pv <- NaN
        }else{
          hyper.pv <- NA
          lookup.res <- hyper.lookup[k==k.hyper]
          if(nrow(lookup.res)>0){
            tmp.lookup.res <- lookup.res[sapply(sorted.TF, function(x) identical(x, select.test.dt[, sort(geneSymbol)])), hyper.pv]
            if(length(tmp.lookup.res)==1){
              hyper.pv <- tmp.lookup.res
            }
          }
          if(nrow(lookup.res)==0|is.na(hyper.pv)){
            enrich.res <- enrichr.wrapper(select.test.dt$geneSymbol, server = enrich.server, background = background, mc.cores = 32) %>% aggr.enrich.table(p.adj.method = TF.enrich.p.adj.method)
            hyper.pv <- quantify.enrichment(enrich.res, AD_genes= AD_genes, fisher.alternative = fisher.alternative, TF.sig = TF.enrich.sig, background = background)$p.value
            hyper.lookup %<>% rbind(data.table(k=k.hyper,
                                               sorted.TF= list(sort(select.test.dt$geneSymbol)),
                                               hyper.pv=hyper.pv))
          }
        }
      }
      hyper.pv.m[i,j] <- hyper.pv
    }
  }
  if(log.p) hyper.pv.m <- -log10(hyper.pv.m)
  if(plot){
    library(latticeExtra)
    
    test.dt <- copy(test.dt)
    test.dt[get(gwas.pv.colname) <= gwas.cutoff, signif:='AD_risk_genes']
    test.dt[get(gwas.pv.colname) > gwas.cutoff, signif:='non_AD_risk_genes']
    p <- hyper.pv.m %>% t %>%
      lattice::levelplot(col.regions = viridis::viridis_pal(direction=-1, option='A')(100), # gray(100:0/100),
                         pretty=F, cuts=50, scales=list(x=list(rot=45)),
                         xlab = sprintf('%s cutoff', x.col.name), ylab=sprintf('%s cutoff', y.col.name))
    x.range <- test.dt[[x.col.name]] %>% range
    y.range <- test.dt[[y.col.name]] %>% range
    x.scale <- length.out/(x.range[2] - x.range[1])
    y.scale <- length.out/(y.range[2] - y.range[1])
    
    scatter.plot.signif <- lattice::xyplot(as.formula(sprintf('y.scale* (%s - y.range[1]) ~ x.scale * (%s - x.range[1])', 
                                                              y.col.name, x.col.name)),
                                           test.dt,
                                           group= signif,
                                           par.settings = list(superpose.symbol = list(pch = c(19, 1),
                                                                                       col = c("orange", "blue"),
                                                                                       alpha= c(1, .2)
                                           )
                                           )
    )
    # plot(p)
    if(enrich=='TF'){
      return(list(hyper.pv.m = hyper.pv.m %>% t,
                  hyper.lookup=hyper.lookup,
                  plt= p + latticeExtra::as.layer(scatter.plot.signif) ))
    }else{
      return(p + latticeExtra::as.layer(scatter.plot.signif) )
    }
  }
  
  else if(return.remaining.n){
    return(remaining.m)
  }else{
    return(hyper.pv.m)
  }
}

pv.trans <- scales::trans_new('pv.trans',
                              transform = function(x){ 
                                x[x<1e-8] <- 1e-8
                                -log10(x) }, 
                              inverse = function(y){10**(-y)})


#' APOE allele to risk approximation
#'Apolipoprotein E (ApoE) is a major cholesterol carrier that supports lipid transport
#' and injury repair in the brain. APOE polymorphic alleles are the main genetic determinants
#'  of Alzheimer disease (AD) risk: individuals carrying the ε4 allele are at increased risk 
#'  of AD compared with those carrying the more common ε3 allele, whereas the ε2 allele 
#'  decreases risk. 
#' @param apoe.genotype 
#'
#' @return
#' @export
#'
#' @examples
apoe.genotype.risk.converter <- function(apoe.genotype) {
  stopifnot(is.integer(apoe.genotype))
  allele1 <- apoe.genotype%%10
  allele2 <- apoe.genotype%/%10
  
  return(sum(allele1, allele2))
}


draw.diff.network <- function(int.g.induced, highlightUID='sample.mean',focus.secP='CTHRC1',genExp.Symbol.dt, secM.components=NULL, layout= 'fruchtermanreingold',
                              plot.text = T, label.text.size=5, plt.toy.network=F,
                              RWR.grad, seed=123, neighbor.degree = 2, reset_APP_p_k=F,
                              plot.grad = 'p_k', plot.grad.scale = F, plot.expr = 'sigmoidExp', color.divergent=NULL,
                              scale.lim.prob = c(.1,.99)) {
  library(igraph)
  if(is.null(secM.components)){
    secM.components <- fread('databases/Human_components_feizi.csv')[, unique(gene_name)]
    
  }
  g <- igraph::induced_subgraph(int.g.induced,
                                igraph::neighborhood(int.g.induced, nodes = focus.secP, order = neighbor.degree)[[1]])
  RWR.grad <- RWR.grad[geneSymbol%in%names(V(g))]
  if(reset_APP_p_k){
    RWR.grad[ID==highlightUID & grad.type== plot.grad &secP==focus.secP & geneSymbol==focus.secP, value:=RWR.grad[ID==highlightUID & grad.type== plot.grad &secP==focus.secP & geneSymbol!=focus.secP, max(value)]]
  }
  
  value.var <- ifelse(plot.grad.scale, 'value.scaled', 'value') ## refer to RWR.grad colname
  
  
  V(g)$expression <- genExp.Symbol.dt[UID==highlightUID][names(V(g)), get(plot.expr), on='geneSymbol']
  V(g)$expression <- (V(g)$expression -min(V(g)$expression, na.rm = T))/ (max(V(g)$expression, na.rm = T) - min(V(g)$expression, na.rm = T))
  V(g)$color.var <- RWR.grad[ID==highlightUID & grad.type== plot.grad &secP==focus.secP][names(V(g)), get(value.var), on='geneSymbol']
  V(g)$dist_to_APP <- distances(g , focus.secP, V(g ), weights = NA)[focus.secP, V(g )]
  if(is.null(color.divergent)){
    color.divergent <- any(V(g)$color.var<0, na.rm = T)
  }
  if(color.divergent){
    color.limit <- RWR.grad[ geneSymbol!=secP& grad.type==plot.grad, quantile(get(value.var), na.rm = T,probs = scale.lim.prob)]
    ## c(-.3, 1) ==> c(-1,1)
    color.limit <- c(-abs(color.limit)[which.max(abs(color.limit))],
                     abs(color.limit)[which.max(abs(color.limit))])
    color.scale <- scale_colour_distiller(name='influence', limits=color.limit,
                                          breaks=c(-1e-5, 0, 1e-5), labels = c('negative', 'neutral', 'positive'),
                                          oob=scales::squish, palette='RdBu')
  }else{
    color.limit <- RWR.grad[ geneSymbol!=secP& grad.type==plot.grad, quantile(get(value.var), na.rm = T,probs = scale.lim.prob)]
    color.scale <- viridis::scale_color_viridis(name= 'component score',limits= color.limit, 
                                                breaks= seq(from=color.limit[1], to=color.limit[2], length.out = 5)[2:4],#RWR.grad[ geneSymbol!=secP& grad.type==plot.grad, quantile(get(value.var), probs = c(0,1/3, 2/3), na.rm = T)],
                                                labels = c('low', 'medium', 'high'),
                                                oob=scales::squish, option = 'magma')
  }
  
  
  ## edge direction
  for (e in E(g)){
    edge.dist.pair <- V(g)[.inc(E(g)[e])]$dist_to_APP
    edge.colorvar.pair <- V(g)[.inc(E(g)[e])]$color.var
    # E(g)[e]$e.alpha <- 1- (sum(edge.dist.pair)-1)/4
    E(g)[e]$e.alpha <- exp(sum(log(edge.colorvar.pair)/2))
    # E(g)[e]$e.alpha <- mean(edge.colorvar.pair)
    if(edge.dist.pair[1]<edge.dist.pair[2]){
      E(g)[e]$e.dir <- 1
    }else if(edge.dist.pair[1]>edge.dist.pair[2]){
      E(g)[e]$e.dir <- -1
    }else{
      E(g)[e]$e.dir <- 0
    }
  }
  
  library(ggnetwork)
  plot.ggnet.graph <- function(net.dt) {
    plt <- net.dt %>% {
      ggplot(., aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(data=net.dt[e.dir==0], aes(alpha=e.alpha), color='black', curvature = .2) +
        geom_edges(data=net.dt[e.dir==1], aes(alpha=e.alpha),
                   # arrow = arrow(length = unit(8, "pt"), type = "open"),
                   color='black') +
        geom_edges(data=net.dt[e.dir==-1], aes(x = xend, y = yend, xend = x, yend = y, alpha=e.alpha), ## reversed
                   # arrow = arrow(length = unit(8, "pt"), type = "open"), 
                   color='black') +
        geom_nodes(data = .[vertex.names!=focus.secP], 
                   aes(color = color.var, #size=expression, 
                       alpha= max(e.alpha, na.rm = T))
        ) +
        # geom_nodes(data=net.dt[vertex.names%in%secM.components], 
        #            aes(color = color.var, size = .5+ expression,  alpha= max(e.alpha, na.rm = T))
        #            # alpha = .5
        # ) +
        # scale_radius(name = 'expression', # range = c(0.5, 15), limits = c(0, 1.5), 
        #              breaks = c(0,.5, 1), 
        #              labels = c('low', 'medium', 'high')) +#, oob=scales::squish)+
        scale_alpha(range = c(0.02, .5)) +
        # geom_nodetext_repel(data = net.dt[!vertex.names%in%secM.components],
        #               aes(label = vertex.names), color='black', alpha=.5) +
        geom_nodelabel(data = net.dt[vertex.names==focus.secP],
                       aes(label = vertex.names),
                       color = 'red', fill='white', size=label.text.size) +
        color.scale +
        guides(alpha=F) +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_blank(),
              panel.background = element_rect(fill = "white"),
              panel.grid = element_blank())# +
      # coord_fixed()
      #legend.position = c(.98, .15))
      # legend.box = "vertical")
    }
  }
  if(plt.toy.network){
    g1 <- induced_subgraph(g, vids = names(V(g)[which(rank(-V(g)$color.var)<=21)]))
    # g2 <- induced_subgraph(g, vids = names(V(g)[which(rank(-V(g)$color.var)<vcount(g)*.8)]))
    set.seed(seed)
    net.dt <- ggnetwork(intergraph::asNetwork(g), layout=layout) %>% as.data.table %>% .[, graph:='n = 100']
    net.dt_small <- net.dt[vertex.names%in%names(V(g1))] %>% 
      .[, .SD[any(.$x%in%xend & .$y%in%yend)], by=.(vertex.names, xend, yend)] %>% .[, graph:='n = 20']
    toy.netowrk <- rbind(net.dt, net.dt_small) %>% 
      .[, graph:=factor(graph, levels = c('n = 20', 'n = 100'))] %>% 
      {
        plot.ggnet.graph(.) +  
          facet_grid(~graph)+
          # guides(alpha=F, color=F, size=F) +
          theme(strip.background = element_blank()
                # strip.text.x = element_blank()
          ) +
          geom_nodetext_repel(data = .[graph=='n = 20' & vertex.names!=focus.secP],
                              aes(label = vertex.names),
                              size=label.text.size)
      } 
    
    # small.graph <- net.dt[vertex.names%in%names(V(g1))] %>% 
    #   .[, .SD[any(.$x%in%xend & .$y%in%yend)], by=.(vertex.names, xend, yend)] %>% plot.ggnet.graph
    big.graph <- net.dt %>% plot.ggnet.graph
    big.graph.labels <-
      if(plot.text){
        big.graph+ geom_nodetext_repel(data = net.dt[vertex.names%in%secM.components],
                                       aes(label = vertex.names),size=label.text.size)
      }else big.graph
    
    
    return(list(big.graph.labels,
                toy.netowrk, rbind(net.dt, net.dt_small)))
  }else{
    set.seed(seed)
    net.dt <- ggnetwork(intergraph::asNetwork(g), layout=layout) %>% as.data.table
    plt <- net.dt %>% plot.ggnet.graph
    if(plot.text){
      plt <- plt+ geom_nodetext_repel(data = net.dt[vertex.names%in%secM.components],
                                      aes(label = vertex.names),size=label.text.size)
    }
    return(plt)
  }
}

CJ.dt <- function(...) {
  DTls <- list(...)
  rows <- do.call(CJ, lapply(DTls, function(x) x[, seq_len(.N)]))
  res <- DTls[[1L]][rows[[1L]]]
  for (n in seq_along(DTls)[-1L])
    res <- res[, c(.SD, DTls[[n]][rows[[n]]])]
  res
}


