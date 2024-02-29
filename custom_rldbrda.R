#' Functions are adapted or copied (because they were not exported) from https://github.com/raeslab/RLdbRDA/tree/main
custom_rldbrda <- function(distmat, meta, p_cutoff=0.05) {
  
  r2 = get_r2(distmat, meta)
  
  sign_r2 = rownames(r2[which(r2$padj < p_cutoff),]) # selects only variables significant in step 1
  
  if (length(sign_r2) < 1) {
    stop("No significant features found, cannot continue!")
  }
  
  cumul <- get_cumul(distmat, meta[, sign_r2])
  
  out <- combine_data(r2, cumul)
  
  return(out)
}

get_r2_single <- function(distmat, meta, feature) {
  capsc <- capscale(distmat ~ meta[, feature], na.action=na.omit)
  an <- anova.cca(capsc)
  
  Fa <- an["F"][[1]][[1]]
  r2 <- RsquareAdj(capsc)[[1]]
  r2adj <- RsquareAdj(capsc)[[2]]
  N <- nrow(na.exclude(meta[,feature,drop=FALSE]))
  pval <- an["Pr(>F)"][[1]][[1]]
  
  output <- cbind(feature, Fa,r2,r2adj,N,pval)
  
  return(output)
}

get_r2 <- function(distmat, meta) {
  features = colnames(meta)
  
  all <- c()
  
  for (feature in features){
    es <- get_r2_single(distmat, meta, feature)
    all <- rbind(all, es)
  }
  
  all <- data.frame(all)
  all$padj <- p.adjust(all$pval,method="BH")
  
  rownames(all) <- all$feature
  
  return(all)
}

get_cumul <- function(distmat, meta) {
  mod0=capscale(distmat ~ 1) #H0: unconstrained ordination
  mod1=capscale(distmat ~ ., data=meta) #H1: full constrained ordination, all metadata
  
  
  attach(meta)
  
  step.res<-ordiR2step(mod0, scope=formula(mod1), data=meta ,direction="forward", Pin = 1, R2scope = TRUE, pstep = 100, perm.max = 999, permutations=9999, trace = F) #forward stepwise dbRDA
  res=step.res$anova
  
  
  row.names(res)=gsub(pattern="\\+ ", "",row.names(res))
  colnames(res)=gsub(pattern="Pr\\(>F\\)", "pval", colnames(res)) # replace column name
  colnames(res)=paste0("RDAcumul_",colnames(res))
  res[,"RDAcumul_N"]=nrow(meta)
  
  detach(meta)
  
  
  return(res)
}

combine_data <- function(r2, cumul){
  
  all=data.frame(merge(r2,cumul,by="row.names",all=T),row.names=1)
  all=all[order(all$r2,decreasing=TRUE),]
  all=all[order(all$RDAcumul_R2.adj),]
  
  return(all)
}

plot_dbrda <- function(plot_data) {
    if (length(unique(plot_data$significant)) > 1) {
      alpha_scale <- scale_alpha_continuous(range = c(0.5, 1))
    } else if (unique(plot_data$significant) == 1) {
      alpha_scale <- scale_alpha_continuous(range = 1)
    } else {
      alpha_scale <-scale_alpha_continuous(range = 0.5)
    }
    
    ggplot(data = plot_data, aes(x = value, y = rowname, fill = variable)) +
      geom_bar(aes(alpha = significant), stat = 'identity', position = position_dodge2()) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
      xlab('Effect size') +
      ylab('Feature') +
      labs(fill = "", alpha = "") + # Hide title of the legend and alpha legend
      guides(alpha = "none") + # Hide alpha legend
      scale_fill_manual(values = c("#60A68B", "#1F4068"),
                        labels = c(bquote(R^2), bquote('Cumulative' ~ R^2))) +
      alpha_scale
  }

