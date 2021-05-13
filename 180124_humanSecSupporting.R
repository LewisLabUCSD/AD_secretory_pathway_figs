
#' Generate model significance statistics using reverse model selection
#'
#' @param m 
#' @param size 
#' @param formula_ 
#' @param scale.cov whether to normalize and reconstruct the model matrix
#'
#' @return
#' @export
#'
#' @examples
model.significance <- function(m, size, use='lm',
                               formula_ = 'HPA.protein.Level.numeric ~ HPA.tpm_log + int.score.mean.diff',
                               diagnostic.plot=F,
                               scale.cov = T, test = 'F', returnCI = F) {
  if(scale.cov){
    m <- cbind(m[, !sapply(m, is.numeric), with=F],m[, lapply(.SD, scale),.SDcols =sapply(m, is.numeric)])
  }
  # model.null <- glm(HPA.protein.Level.numeric ~ 1, data = m)
  # model.tpm <- glm(HPA.protein.Level.numeric ~ HPA.tpm_log, data = m)
  if(use=='lm'){
    model.tpm.intscore <-
      lm(formula = as.formula(formula_),
         data = m)
    model.significance.dt <-
      data.table(
        term = names(coef(model.tpm.intscore)),
        coefficient = coef(model.tpm.intscore)
        # as.data.table(drop1(model.tpm.intscore, test = 'F'
      ) %>%
      setkey('term') %>% 
      .[ as.data.table(drop1(model.tpm.intscore, test = test), keep.rownames = 'term')[1, term:='(Intercept)']]
    # model.significance.dt[, `-logPr(>F)` := -log(`Pr(>F)`)]
    model.significance.dt %>% setnames('F value', 'F.statistic')
    if(!'AIC'%in%names(model.significance.dt)){
      model.significance.dt[, AIC:=NA_real_]
    }
    if(returnCI){
      CI = confint(model.tpm.intscore) %>% as.data.table(keep.rownames = 'term')
      model.significance.dt <- merge(model.significance.dt, CI, by = 'term')
    }
    model.significance.dt[, rmse:= sum(model.tpm.intscore$residuals**2)]
  }else if(use=='lmer'){
    message('use mixed effects modelling')
    model.tpm.intscore <-
      lmerTest::lmer(formula = as.formula(formula_),data = m)
    model.significance.dt <- model.tpm.intscore %>% anova %>% as.data.table(keep.rownames = 'term') 
    model.significance.dt[, coefficient:=model.tpm.intscore@beta]
    if(returnCI){
      CI = model.tpm.intscore %>% confint %>% .[-c(1:4),] %>% as.data.table(keep.rownames = 'term')
      model.significance.dt <- merge(model.significance.dt, CI, by = 'term')
    }
    
    if(diagnostic.plot){
      gene.fitness <- coef(model.tpm.intscore)$Symbol %>% as.data.table(keep.rownames = 'Gene name') %>% 
        merge(predicted.secP.HPA[, lapply(.SD,mean), by=.(HPA.predicted.type,`Gene name`), .SDcols=Genes_FPKM:Proteins_iBAQ.sigmoid.byTissue], by='Gene name')
      gene.fitness %>%
        { nhighlight <- 20
        highlight.dt <- .[,.SD[order(-secretion.fitness)][1:nhighlight],by=HPA.predicted.type]
        ggplot(data = .,aes(secretion.fitness, color=HPA.predicted.type, fill=HPA.predicted.type)) + 
          geom_vline(data=highlight.dt, 
                     aes(xintercept=secretion.fitness, color=HPA.predicted.type), alpha=.5, linetype='longdash')+
          ggrepel::geom_text_repel(data=highlight.dt,
                                   aes(secretion.fitness, 0.1, label=`Gene name`))+
          geom_density(alpha=.3) +facet_grid(HPA.predicted.type~.)
        }
    }
  }
  return(model.significance.dt)
}


#' Plot coefficients, pvalues from MLR result tables generated from model.significance
#'
#' @param . 
#'
#' @return
#' @export
#'
#' @examples
plt.model.sig <- function(.) {
  p1 <- ggplot(data = .,
               aes(term, `Pr(>F)`, color = term, fill = term)) + geom_point() + scale_y_continuous(limits = c(NA,.05)) +
    coord_flip() +
    facet_grid(HPA.predicted.type~.) + theme(axis.title.y=element_blank(),
                                             axis.text.y=element_blank(),
                                             axis.ticks.y=element_blank())
  
  p2 <- ggplot(data = ., aes(term, coefficient, color = term)) + 
    geom_pointrange(aes(ymin=`2.5 %`, ymax=`97.5 %`)) +# ylim(c(-.5, .5)) + 
    geom_hline(yintercept = 0, color = 'green', linetype = "longdash") +
    facet_grid(HPA.predicted.type~.) + coord_flip()+ theme(axis.title.y=element_blank(),
                                                           axis.text.y=element_blank(),
                                                           axis.ticks.y=element_blank())
  
  ggpubr::ggarrange(p1, p2, common.legend = T)
  # geom_point() + geom_errorbar(aes(ymin=`2.5 %`, ymax=`97.5 %`), colour="black", width=.1, position=pd) +
  
  # p2 <- ggplot(data = ., aes(coefficient, color = term, fill = term)) +
  #   geom_density(alpha = .5) +
  #   facet_grid(HPA.predicted.type~.) + xlim(c(-3,3))
  # plot(p2)
  #
  # p3 <- ggplot(data = ., aes(term, coefficient, color = term, fill = term)) + geom_point() +
  #   facet_grid(HPA.predicted.type~my_q) + ylim(c(-1.1, 1.1)) + geom_hline(yintercept = 0, color = 'green', linetype = "longdash")
  # plot(p3)
}

plt.model.sig.multimodel <- function(.) {
  p1 <- ggplot(data = .,
               aes(term, `Pr(>F)`, color = term, fill = term)) + 
    geom_violin(alpha = .3) + geom_boxplot(alpha = .5, width = .1) +
    facet_grid(HPA.predicted.type~.) +
    # coord_flip() +
    # facet_grid(HPA.predicted.type~.) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  p2 <- ggplot(data = ., aes(term, coefficient, color = term, fill = term)) + 
    geom_violin(alpha = .3) + geom_boxplot(alpha = .5, width = .1) + 
    facet_grid(HPA.predicted.type~.) +
    ylim(c(-1.1, 1.1)) + geom_hline(yintercept = 0, color = 'green', linetype = "longdash") +
    # plot(p3)
    # geom_pointrange(aes(ymin=`2.5 %`, ymax=`97.5 %`)) +# ylim(c(-.5, .5)) + 
    # geom_hline(yintercept = 0, color = 'green', linetype = "longdash") +
    # facet_grid(HPA.predicted.type~.) + coord_flip()+ 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  ggpubr::ggarrange(p1, p2, common.legend = T)
  # geom_point() + geom_errorbar(aes(ymin=`2.5 %`, ymax=`97.5 %`), colour="black", width=.1, position=pd) +
  
  # p2 <- ggplot(data = ., aes(coefficient, color = term, fill = term)) +
  #   geom_density(alpha = .5) +
  #   facet_grid(HPA.predicted.type~.) + xlim(c(-3,3))
  # plot(p2)
  #
  # p3 <- ggplot(data = ., aes(term, coefficient, color = term, fill = term)) + geom_point() +
  #   facet_grid(HPA.predicted.type~my_q) + ylim(c(-1.1, 1.1)) + geom_hline(yintercept = 0, color = 'green', linetype = "longdash")
  # plot(p3)
}


#' scatter plot with coef
#'
#' @param quant.dt 
#' @param enrichment.score.col the col.name for the normalized enrichment score
#' @param expr.col the col.name for the expression to correlate with 'enrichment.score.col'
#' @param return.plot 
#' @param cor.use.methods a character string indicating which correlation coefficient (or covariance) is 
#' to be computed. One of "pearson" (default), "kendall", or "spearman": can be abbreviated. 
#' @param x_lim_cutoff_probs 
#' @param adj.bw defaults to NULL, where densities are estimated empirically through the MASS package. Additionally, a vector of length of two specifying
#' the adjustments factor to the x and y-axis bandwidth respectively.
#'
#' @return
#' @export
#'
#' @examples
plot.scatter.corr <- function(quant.dt, x_lim_cutoff_probs = c(.05, .95),
                              enrichment.score.col = 'enrichment.score',
                              expr.col = 'HPA.tpm',
                              return.plot = F, cor.use.methods = 'pearson', adj.bw = NULL) {
  quant.dt <- quant.dt[int.score>0]
  margin.max <- max(abs(quant.dt$enrichment.score %>% quantile(x_lim_cutoff_probs)))
  # scatter.value.dt <- quant.dt[, .(x = 0,
  #                                  y = 10000,
  #                                  rank.cor = paste0('r=', round(
  #                                    cor(
  #                                      get(enrichment.score.col),
  #                                      get(expr.col),
  #                                      method = 'spearman'),
  #                                    digits = 4))),
  #                              by = .(int.g, HPA.predicted.type)]
  
  plt.scatter <- quant.dt %>% 
    ggplot(aes_string(enrichment.score.col, expr.col))+
    geom_point(alpha = .05)+
    # geom_smooth(method = 'lm')+
    xlim(c(-margin.max, margin.max))+
    facet_grid(int.g ~ HPA.predicted.type)+
    stat_cor(method = cor.use.methods, color = 'red')+
    theme(legend.position="none")
  
  if(!is.null(adj.bw)){
    g(bw.adj.x, bw.adj.y) %=% adj.bw
    scatter.data <- plt.scatter %>% layer_data()
    bw.x <- scatter.data$x %>% na.omit %>% MASS::bandwidth.nrd() %>% `*`(bw.adj.x)
    bw.y <- scatter.data$y %>% na.omit %>% MASS::bandwidth.nrd() %>% `*`(bw.adj.y)
    plt.scatter <- plt.scatter + stat_density2d(aes_string(alpha='..level..', fill='..level..'), geom="polygon", h = c(bw.x, bw.y)) +
      scale_fill_gradient(low = "yellow", high = "red") # +ggtitle(sprintf('%s; %s', type, secM.expr.col))
  }else{
    plt.scatter <- plt.scatter + stat_density2d(aes_string(alpha='..level..', fill='..level..'), geom="polygon", h = NULL) +
      scale_fill_gradient(low = "yellow", high = "red") # +ggtitle(sprintf('%s; %s', type, secM.expr.col))
  }
  if(return.plot==F){
    plot(plt.scatter)
    return(NULL)
  }else return(plt.scatter)
}


## get plot margins
plot.trend <- function(int.score.protein.quant.dt, type, use.subgraph, secM.expr.col, use.weights, return.plot=F){
  plt.margin <-
    int.score.protein.quant.dt[, enrichment.score %>% quantile(probs = c(.1, .9))]
  subtitle <-
    sprintf(
      'propagation using %s with %s as subgraph; secM expression from %s; %s graph enrichment',
      type,
      use.subgraph,
      secM.expr.col,
      ifelse(use.weights, 'confidence level weighted', 'unweighted')
    )
  overall.trend.plot <- int.score.protein.quant.dt %>%
    # .[int.score != 0] %>% ## remove secPs that dont interact with any of the secMs
    ggplot(aes(enrichment.score, color = HPA.protein.Level)) +
    facet_grid(int.g ~ HPA.predicted.type, scales = 'free') +
    geom_vline(xintercept = 0, color = 'grey', linetype = "longdash") +
    xlim(plt.margin) +
    scale_color_brewer(palette = "RdYlBu", direction = -1)
  overall.trend.plot.density <-
    overall.trend.plot + geom_density(alpha = .025) +
    ggtitle('Density plot: int.score vs protein level', subtitle = subtitle)
  overall.trend.plot.ecdf <- overall.trend.plot + stat_ecdf() + geom_hline(yintercept = .5, color = 'grey', linetype = "longdash") +
    ggtitle('Cumulative density plot: int.score vs protein level',
            subtitle = subtitle)
  # plot(overall.trend.plot.density)
  # plot(overall.trend.plot.ecdf)
  # overall.box.plot <- int.score.protein.quant.dt %>%
  #   .[int.score != 0] %>% ## remove secPs that dont interact with any of the secMs
  #   ggplot(aes(HPA.protein.Level, int.score.mean.diff, color = HPA.protein.Level)) +
  #   facet_grid(int.g ~ HPA.predicted.type, scales = 'free') +
  #   geom_violin() + geom_boxplot(width = 0.1) + scale_color_brewer(palette = "RdYlBu", direction =
  #                                                                    -1) + ylim(c(-1, 1))
  p <- gridExtra::arrangeGrob(overall.trend.plot.density, overall.trend.plot.ecdf, nrow = 1)
  if(return.plot==F){
    plot(p)
    return(NULL)
  }else return(p)
  
}



plot.secP.model.sig <- function(secP.model.significance.by.secP, secM.expr.col = NULL, use.weights = NULL, type = NULL, scatter.only = F){
  # use.subgraph <- 'all.secretory.resident.genes'
  secP.model.significance.by.secP <- secP.model.significance.by.secP[, .SD[!any(is.na(coefficient), is.na(`Pr(>F)`))], by = .(HPA.predicted.type, Symbol)]
  plt.margin.max <- max(abs(secP.model.significance.by.secP[, coefficient %>% quantile(probs = c(.05, .95))]))
  if(!(is.null(secM.expr.col) | is.null(use.weights) | is.null(type))){
    subtitle <-
      sprintf(
        '%s propagation; secM expression from %s; %s graph enrichment',
        type,
        gsub('^HPA.', '', secM.expr.col),
        ifelse(use.weights, 'CI weighted', 'unweighted')
      )
  }else subtitle <- ''
  
  scatter.model.sig.secP <- secP.model.significance.by.secP %>% 
    ggplot(aes(coefficient, `Pr(>F)`, color = term, size = F.statistic)) + 
    # facet_grid(int.g ~ HPA.predicted.type, scales = 'free') +
    facet_grid( ~ HPA.predicted.type, scales = 'free') +
    geom_point(alpha =ifelse(scatter.only, 1, .2)) +
    geom_rug(alpha = .1) + 
    xlim(c(-plt.margin.max, plt.margin.max)) +
    ylim(c(0,1)) + # limit for pvalues
    geom_vline(xintercept = 0, color = 'green', linetype = "longdash") +
    labs(subtitle = 'coefficient vs significance')
  if(scatter.only) {
    gridExtra::grid.arrange(scatter.model.sig.secP, nrow = 1, top = subtitle)
    return(NULL)
  }
  box.model.sig.secP <- secP.model.significance.by.secP %>% 
    ggplot(aes(term, `Pr(>F)`, color = term)) +
    # facet_grid(int.g ~ HPA.predicted.type, scales = 'free') + 
    facet_grid( ~ HPA.predicted.type, scales = 'free') + 
    geom_hline(yintercept = 0.05, color = 'green', linetype = "longdash") + 
    geom_boxplot() +
    labs(subtitle = 'significance comparison')
  
  density.model.sig.secP <- secP.model.significance.by.secP %>% 
    ggplot(aes(coefficient, color = term)) +
    # facet_grid(int.g ~ HPA.predicted.type, scales = 'free') + 
    facet_grid( ~ HPA.predicted.type, scales = 'free') + 
    geom_vline(xintercept = 0, color = 'green', linetype = "longdash") + 
    geom_hline(yintercept = .5, color = 'green', linetype = "longdash") + 
    xlim(c(-plt.margin.max, plt.margin.max)) +
    stat_ecdf() +
    labs(subtitle = 'coefficient comparison')
  
  # plot(plot.model.sig.secP)
  gridExtra::grid.arrange(scatter.model.sig.secP, box.model.sig.secP, density.model.sig.secP, nrow = 2,  top = subtitle)
  return(NULL)
}


stat_cor <- function(mapping = NULL, data = NULL,
                     method = "pearson", label.sep = ", ",
                     label.x.npc = "left", label.y.npc = "top",
                     label.x = NULL, label.y = NULL, output.type = "expression",
                     geom = "text", position = "identity",  na.rm = FALSE, show.legend = NA,
                     inherit.aes = TRUE, ...) {
  parse <- ifelse(output.type == "expression", TRUE, FALSE)
  layer(
    stat = StatCor, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(label.x.npc  = label.x.npc , label.y.npc  = label.y.npc,
                  label.x = label.x, label.y = label.y, label.sep = label.sep,
                  method = method, output.type = output.type, parse = parse, na.rm = na.rm, ...)
  )
}


StatCor<- ggproto("StatCor", Stat,
                  required_aes = c("x", "y"),
                  default_aes = aes(hjust = ..hjust.., vjust = ..vjust..),
                  
                  compute_group = function(data, scales, method, label.x.npc, label.y.npc,
                                           label.x, label.y, label.sep, output.type)
                  {
                    if (length(unique(data$x)) < 2) {
                      # Not enough data to perform test
                      return(data.frame())
                    }
                    # Returns a data frame with estimate, p.value, label, method
                    .test <- .cor_test(
                      data$x, data$y, method = method, label.sep = label.sep,
                      output.type = output.type
                    )
                    # Returns a data frame with label: x, y, hjust, vjust
                    .label.pms <- ggpubr:::.label_params(data = data, scales = scales,
                                                         label.x.npc = label.x.npc, label.y.npc = label.y.npc,
                                                         label.x = label.x, label.y = label.y ) %>%
                      mutate(hjust = 0)
                    cbind(.test, .label.pms)
                  }
)





# Correlation test
#::::::::::::::::::::::::::::::::::::::::
# Returns a data frame: estimatel|p.value|method|label
.cor_test <- function(x, y, method = "pearson", label.sep = ", ", output.type = "expression"){
  
  .cor <- stats::cor.test(x, y, method = method, exact = FALSE, use = "complete.obs")
  estimate <- p.value <- p <- r <- rr <-  NULL
  z <- data.frame(estimate = .cor$estimate, p.value = .cor$p.value, method = method) %>%
    mutate(
      r = signif(estimate, 2),
      rr = signif(estimate^2, 2),
      p = signif(p.value, 2)
    )
  
  
  # Defining labels
  pval <- .cor$p.value
  
  if(output.type == "expression"){
    
    z <- z %>% dplyr::mutate(
      r.label = paste("italic(R)", r, sep = "~`=`~"),
      rr.label = paste("italic(R)^2", rr, sep = "~`=`~"),
      p.label = paste("italic(p)", p, sep = "~`=`~")
    )
    # Default label
    pvaltxt <- ifelse(pval < 2.2e-16, "italic(p)~`<`~2.2e-16",
                      paste("italic(p)~`=`~", signif(pval, 2)))
    cortxt <- paste0("italic(R)~`=`~", signif(.cor$estimate, 2),
                     "~`,`~",  pvaltxt)
    z$label <- cortxt
    
  }
  else if (output.type %in% c("latex", "tex", "text")){
    
    z <- z %>% dplyr::mutate(
      r.label = paste("R", r, sep = " = "),
      rr.label = paste("R2", rr, sep = " = "),
      p.label = paste("p", p, sep = "=")
    )
    
    # Default label
    pvaltxt <- ifelse(pval < 2.2e-16, "p < 2.2e-16",
                      paste("p =", signif(pval, 2)))
    cortxt <- paste0("R = ", signif(.cor$estimate, 2),
                     " , ",  pvaltxt)
    z$label <- cortxt
  }
  p.value.signif <- if(.cor$p.value<1e-4){
    '`****`' 
  }else if(.cor$p.value<1e-3){
    '`***`'
  }else if(.cor$p.value<1e-2){
    '`**`'
  }else if(.cor$p.value<.05){
    '`*`'
  }else{ 'ns'}
  
  z['p.signif'] <- p.value.signif
  return(z)
}



plot.scatter.corr.pub <- function(quant.dt, .xlim = NA,
                                  enrichment.score.col = 'enrichment.score', filter0 = T, facet = T, 
                                  ylab. = 'expression (tpm)',
                                  expr.col = 'HPA.tpm',
                                  return.plot = F, cor.use.methods = 'pearson', adj.bw = NULL, add.density = T, .ylim = NULL, ylog = T, norm.y = F) {
  if(filter0) quant.dt <- quant.dt[int.score>0] ## todo change col name
  if(norm.y) quant.dt[, expr:= (expr - min(expr))/ max(expr - min(expr))]
  # margin.max <- max(abs(quant.dt$enrichment.score %>% quantile(x_lim_cutoff_probs)))
  plt.scatter <- quant.dt %>% 
    ggplot(aes_string(enrichment.score.col, expr.col))+
    geom_point(alpha = .05)+
    # geom_smooth(method = 'lm')+
    xlim(.xlim)+
    # facet_grid(~ HPA.predicted.type)+
    stat_cor(aes(label = paste(..r.label.., ..p.signif.., sep = "~`,`~")),
             method = cor.use.methods, color = "black", label.x.npc = 'center', label.y.npc = 'top') + # t-test for beta https://stats.stackexchange.com/questions/270612/why-test-statistic-for-the-pearson-correlation-coefficient-is-frac-r-sqrtn-2
    xlab('machinery support score') + 
    ylab(ylab.) +
    theme_bw() %+replace%
    theme(legend.position="none",
          axis.title.x=element_blank())
  
  if(ylog) plt.scatter <- plt.scatter + scale_y_log10(limits = .ylim)
  
  if(add.density){
    if(!is.null(adj.bw)){
      g(bw.adj.x, bw.adj.y) %=% adj.bw
      scatter.data <- layer_data(plt.scatter)
      bw.x <- scatter.data$x %>% na.omit %>% MASS::bandwidth.nrd() %>% `*`(bw.adj.x)
      bw.y <- scatter.data$y %>% na.omit %>% MASS::bandwidth.nrd() %>% `*`(bw.adj.y)
      plt.scatter <- plt.scatter + stat_density2d(aes_string(alpha='..level..', fill='..level..'), geom="polygon", h = c(bw.x, bw.y)) +
        viridis::scale_fill_viridis(option = 'plasma') # +ggtitle(sprintf('%s; %s', type, secM.expr.col))
    }else{
      plt.scatter <- plt.scatter + stat_density2d(aes_string(alpha='..level..', fill='..level..'), geom="polygon") +
        viridis::scale_fill_viridis(option = 'plasma') # +ggtitle(sprintf('%s; %s', type, secM.expr.col))
    }}
  
  if(facet) plt.scatter <- plt.scatter + facet_grid( ~ HPA.predicted.type)
  
  if(return.plot==F){
    plot(plt.scatter)
    return(NULL)
  }else return(plt.scatter)
}



#' Compare models
#'
#' @param combo.row 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' ## total combinations
#' benchmark.combo <- 
#'   expand.grid(used.int.g = 'PCNet',
#'               use.tissue.median = F, # use all tissues
#'               protein.data.source = 'deep_proteome',
#'               secM.expr.col = 
#'                 c('Proteins_iBAQ.rank.gene', 'Proteins_iBAQ.rank.byTissue', 'Proteins_iBAQ.sigmoid'), #Proteins_iBAQ.sigmoid.byTissue
#'               Y.var = c('Proteins_iBAQ.rank.gene', 'Proteins_iBAQ.rank.byTissue', 'Proteins_iBAQ.sigmoid'),
#'               mRNA.var = c('log(1+Genes_FPKM)'),
#'               secretion.fitness.var = c('int.score.mean.RWR3', 'int.score.diff.RWR3.RWR2', 'int.score.diff.RWR3.RWR_half'),
#'               stringsAsFactors = F
#'   ) %>% as.data.table() %>% split(by = names(.),drop = T)
#' 
#' 
#' all.combo.stats <- 
#'  parallel::mclapply(benchmark.combo, load.and.model, mc.cores = 32) %>% rbindlist()

load.and.model <- function(combo.row, mixed_modelling=F, returnCI=F) {
  ## data load
  version = sprintf('RWR4_190222_%s_%s', combo.row$protein.data.source, combo.row$used.int.g)
  save.dir <- sprintf('output/180224_graphRWR/%s/%s', version, ifelse(combo.row$use.tissue.median, 'tissue.median', 'allTissues'))
  prop.int.score.dt.fn <- sprintf('%s/%s_prop.int.score.dt.csv.zip', save.dir, combo.row$secM.expr.col)
  prop.int.score.dt <- rio::import(prop.int.score.dt.fn) %>% as.data.table
  
  prop.int.score.dt.quant.tissues <- prop.int.score.dt %>%
    merge(predicted.secP.HPA, by.x= c('Symbol', 'Tissue'), by.y = c('Gene name', 'Tissue'))
  prop.int.score.dt.quant.tissues[HPA.predicted.type == 'membrane.secreted', HPA.predicted.type:= 'secreted']
  m.tissues <- prop.int.score.dt.quant.tissues[ collapse.by=='secP'] %>%
    # .[HPA.predicted.type=='intracellular'] %>% 
    # .[, .(Symbol, HPA.predicted.type, Y = HPA.protein.Level.numeric, X1 = log(1+HPA.tpm), X2 = int.score.diff.RWR3.RWR_half)] 
    .[, .(Symbol, HPA.predicted.type, 
          Y = eval(parse(text = combo.row$Y.var)),
          mRNA = eval(parse(text = combo.row$mRNA.var)),
          secretion.fitness = eval(parse(text = combo.row$secretion.fitness.var))
    )] 
  
  test.stat <- m.tissues[, .(var.Y = var(Y),
                             mean.Y = mean(Y),
                             median.Y = median(Y),
                             var.mRNA = var(mRNA),
                             mean.RNA = mean(mRNA),
                             median.mRNA = median(mRNA)),
                         by = .(HPA.predicted.type, Symbol)] 
  
  if(mixed_modelling){
    m0 <- m.tissues %>% 
      .[ ,.SD[var(Y)>0 & mean(Y) > 0 & sd(mRNA)>0 & var(secretion.fitness)>0], by = .(HPA.predicted.type, Symbol)]
    
    m0_rank <- m.tissues[, lapply(.SD, rank), .SDcols = Y:secretion.fitness, by = .(Symbol, HPA.predicted.type)] %>% 
      .[ ,.SD[var(Y)>0 & mean(Y) > 0 & sd(mRNA)>0 & var(secretion.fitness)>0], by = .(HPA.predicted.type, Symbol)]
    # 
    # d_1 <- m0 %>%
    #   .[,model.significance(.SD,formula_ = 'Y ~ 0 + mRNA + secretion.fitness + (1|Symbol)', test = 'F', use='lmer'), 
    #     by = .(HPA.predicted.type)] %>% 
    #   .[, lmer.type := 'rand.intercept'] %>% .[, var.rank := F]
    
    d_2 <- m0 %>%
      .[,model.significance(.SD,formula_ = 'Y ~ 0 + mRNA + secretion.fitness + (secretion.fitness|Symbol)', test = 'F', use='lmer')] %>% 
      .[, lmer.type := 'rand.slope_sec'] %>% .[, var.rank := F]
    
    
    return(cbind(combo.row, d_2))
    # 
    # d_3 <- m0 %>%
    #   .[,model.significance(.SD,formula_ = 'Y ~ 0 + mRNA + secretion.fitness + (mRNA|Symbol)', test = 'F', use='lmer'), 
    #     by = .(HPA.predicted.type)] %>% 
    #   .[, lmer.type := 'rand.slope_mRNA'] %>% .[, var.rank := F]
    # 
    # d_4 <- m0_rank %>%
    #   .[,model.significance(.SD,formula_ = 'Y ~ 0 + mRNA + secretion.fitness + (1|Symbol)', test = 'F', use='lmer'), 
    #     by = .(HPA.predicted.type)] %>% 
    #   .[, lmer.type := 'rand.intercept'] %>% .[, var.rank := T]
    # 
    # d_5 <- m0_rank %>%
    #   .[,model.significance(.SD,formula_ = 'Y ~ 0 + mRNA + secretion.fitness + (secretion.fitness|Symbol)', test = 'F', use='lmer'), 
    #     by = .(HPA.predicted.type)] %>% 
    #   .[, lmer.type := 'rand.slope_sec'] %>% .[, var.rank := T]
    # 
    # d_6 <- m0_rank %>%
    #   .[,model.significance(.SD,formula_ = 'Y ~ 0 + mRNA + secretion.fitness + (mRNA|Symbol)', test = 'F', use='lmer'), 
    #     by = .(HPA.predicted.type)] %>% 
    #   .[, lmer.type := 'rand.slope_mRNA'] %>% .[, var.rank := T]
    # 
    # return(cbind(combo.row, rbind(d_1, d_2, d_3, d_4, d_5, d_6)))
    
  }else{
    d2 <- m.tissues %>% 
      # .[Symbol %in% model.used.genes ] %>%
      .[ ,.SD[var(Y)>0 & mean(Y) > 0 & sd(mRNA)>0 & var(secretion.fitness)>0], by = .(HPA.predicted.type, Symbol)] %>%
      .[,model.significance(.SD,formula_ = 'Y ~ 0 + mRNA + secretion.fitness', test = 'F'), by = .(HPA.predicted.type, Symbol)] %>% 
      .[, var.rank := F]
    
    d2.rank <- m.tissues[, lapply(.SD, rank), .SDcols = Y:secretion.fitness, by = .(Symbol, HPA.predicted.type)] %>% 
      .[ ,.SD[var(Y)>0 & mean(Y) > 0 & sd(mRNA)>0 & var(secretion.fitness)>0], by = .(HPA.predicted.type, Symbol)] %>% 
      .[,model.significance(.SD,formula_ = 'Y ~ 0 + mRNA + secretion.fitness', test = 'F'), by = .(HPA.predicted.type, Symbol)] %>% 
      .[, var.rank := T]
    return(cbind(combo.row, rbind(d2, d2.rank)))
  }
}


plot_grad_targets_barplots <- 
  function(.SD, secP, n.top = 50){
    top.mask.dt <-
      .SD[targetGene!=secP][!targetGene%in%c(secM.amir$Symbol)][order(-abs(p_pk__p_mask))][1:n.top]
    top.mask.dt[, targetGene:=reorder(targetGene, abs(p_pk__p_mask))]
    top.mask.dt[, c('p_pk__p_mask.scale',
                    'effect'):= .(p_pk__p_mask/ max(abs(p_pk__p_mask)),
                                  ifelse(p_pk__p_mask>0, 'positive', 'negative')
                    )]
    top.mask.dt[, effect:=factor(effect, levels = c('negative', 'positive'))]
    p_mask <- ggplot(top.mask.dt,
                     aes(targetGene, p_pk__p_mask.scale, fill = HPA.predicted.type)) + geom_bar(stat = 'identity') +
      coord_flip() + ggtitle(secP, subtitle = 'other genes')+ scale_fill_discrete(drop=FALSE) + 
      theme(legend.position="bottom",
            axis.title.y = element_blank())
    
    top.sec.expr.dt <- 
      .SD[targetGene!=secP][targetGene%in%c(secM.amir$Symbol)][order(-abs(p_pk__p_expr))][1:n.top]#[targetGene%in%c(secM.amir$Symbol)]
    top.sec.expr.dt[, targetGene:=reorder(targetGene, abs(p_pk__p_expr))]
    top.sec.expr.dt[, c('p_pk__p_expr.scale',
                        'effect'):= .(p_pk__p_expr/ max(abs(p_pk__p_expr)),
                                      ifelse(p_pk__p_expr>0, 'positive', 'negative')
                        )]
    top.sec.expr.dt[, effect:=factor(effect, levels = c('negative', 'positive'))]
    p_expr <- ggplot(top.sec.expr.dt,aes(targetGene, p_pk__p_expr.scale, fill = Subsystem.aggr)) + geom_bar(stat = 'identity') +
      coord_flip() + ggtitle(secP, subtitle = 'secM genes') + scale_fill_discrete(drop=FALSE) + 
      theme(legend.position="bottom",
            axis.title.y = element_blank())
    
    gridExtra::grid.arrange(p_expr+ylim(c(-1,1)), p_mask+ylim(c(-1,1)), ncol = 2, nrow = 1)
    return(NULL)
    # .SD[targetGene!=secP][order(abs(p_pk__p_mask))]
  }
