suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(glmnet))

tpm = function(counts.mat, gene.length) {
    x <- counts.mat / gene.length
    tpm.mat <- t( t(x) * 1e6 / colSums(x) )
    return(x)
}

selectGenes <- function(counts, min.count=10, N=0.90){
 
  lib.size <- colSums(counts)
  MedianLibSize <- median(lib.size)
  CPM.Cutoff <- min.count / MedianLibSize*1e6
  CPM <- edgeR::cpm(counts,lib.size=lib.size)
 
  min.samples <- round(N * ncol(counts))
 
  f1 <- genefilter::kOverA(min.samples, CPM.Cutoff)
  flist <- genefilter::filterfun(f1)
  keep <- genefilter::genefilter(CPM, flist)
 
  return(keep)
}

run_ComBat = function(df_list, ann, batch_col, join_by="Geneid", cov_col=NULL, mean.only=F, ref.batch=NULL, int_counts=F) {
    col_names = df_list %>% lapply(FUN = function(x) {
        colnames(x)[-1]
    })
    full_df = df_list %>% reduce(inner_join, by = join_by) %>% column_to_rownames(var="Geneid")
    samples = colnames(full_df)
    ann = ann %>% slice(match(samples, sample_id))
    batch = ann %>% pull(batch_col)
    
    if (!is.null(cov_col)) {
        formula = paste0("~as.factor(", cov_col, ")")
        mod = model.matrix(as.formula(formula), data = ann)
    } else {
        mod=NULL
    }

    if (!int_counts) {
        combat_edata = ComBat(full_df, batch = batch, mean.only = mean.only, mod = mod, ref.batch=ref.batch)
    } else {
        full_df = as.matrix(full_df)
        combat_edata = ComBat_seq(full_df, batch = batch, covar_mod = mod)
    }
    out = lapply(col_names, FUN = function(x) {
        df = combat_edata %>% as.data.frame() %>% select(one_of(x)) %>% t() %>% as.data.frame()
        return(df)
    })
    return(out)
}

y_weights <- function(xy){
  cts <- xy %>% dplyr::count(y)
  tot <- sum(cts$n)
  cts$Weight <- 1 - cts$n / tot; cts$n <- NULL; names(cts) <- c("y", "Weight")

  wdf <- inner_join(data.frame(y=xy$y), cts, by="y")
  weights <- wdf$Weight

  return(weights)
}

build_model <- function(xy, alpha, nfolds = 5, family="gaussian", type.measure="mse", weighted=T, standardize=T, keep=T){
  xy <- xy %>% arrange(y)
  x <- as.matrix(xy %>% select(-y))
  y <- as.matrix(xy %>% select(y))

  if (weighted){weights <- y_weights(xy)}
  else{weights <- NULL}

  fit <- cv.glmnet(x, y, nfolds=nfolds, family=family, type.measure=type.measure, weights=weights, alpha=alpha, standardize=standardize, keep=keep)

  return(fit)
}

nzcs <- function(coefs){
  nzc <- coefs@Dimnames[[1]][coefs@i + 1]
  nzc <- nzc[2:length(nzc)]
  
  return(nzc)
}
                         
extract_nzc <- function(fit, lm, family){
  coefs <- coef(fit, s=lm)
  if (family == "multinomial"){nzc <- unique(unlist(lapply(coefs, nzcs)))}
  else(nzc <- nzcs(coefs))

  return(nzc)
}

make_boxplot = function(df, ann_df, x, y="QDS", byvar=c("sample_id"="sample_id"), color=NULL, title=NULL, 
                        levels=NULL, shape=NULL) {
    df = inner_join(df, ann_df, by=byvar) 
    if (!is.null(levels)) {
        df[[x]] = factor(df[[x]], ordered=T, levels=levels)
    }
    p = ggplot(df, aes_string(x = x, y = y)) + 
        geom_boxplot(outlier.shape = NA) + geom_jitter(
            height = 0, width = 0.2, aes_string(col=color, shape=shape)) + 
        ggtitle(title)
    return(p)
}