# Preliminary script to compute Figure S30 and S31 from the supplement
# - "classical" correlation matrix.

create_cor_mat <- function(corr_data) {
  cor_mat <- signif_mat <- array(NA,dim=c(max(floor(unique(c(corr_data$from,corr_data$to)))),max(floor(unique(c(corr_data$from,corr_data$to))))))
  colnames(cor_mat) <- rownames(cor_mat) <- colnames(signif_mat) <- rownames(signif_mat) <- 1:dim(cor_mat)[1]
  for (i in 1:length(corr_data$from)) {
    cor_mat[floor(corr_data$from[i]),floor(corr_data$to[i])] <- max(-1,min(1,corr_data$corr[i]))
    signif_mat[floor(corr_data$from[i]),floor(corr_data$to[i])] <- sign(corr_data$diff_sig[i])
  }
  keep_ind <- which(sapply(1:dim(cor_mat)[1],function(i) length(which(!is.na(c(cor_mat[i,],cor_mat[,i])))))!=0)
  cor_mat <- cor_mat[keep_ind,keep_ind]
  signif_mat <- signif_mat[keep_ind,keep_ind]
  diag(cor_mat) <- diag(signif_mat) <- 1
  return(list(cor_mat=cor_mat,signif_mat=signif_mat))
}

clrrmp <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                                 "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061")))

pdf(file.path(DIR_FIGURES,'corrplot_ACERAP.pdf'), width = 10.5, height = 9.5)
corr_obj <- create_cor_mat(ACERtrc$ACER_ap$corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_probit_arb_pollen_data)
corrplot::corrplot(corr_obj$cor_mat,type="upper",p.mat=1-corr_obj$signif_mat,pch.cex=0.75,col=clrrmp(200),
                   mar = c(0,0,1,0), number.cex = 0.5, title = 'ACER AP')
dev.off()

pdf(file.path(DIR_FIGURES,'corrplot_TRACETS.pdf'), width = 10.5, height = 9.5)
corr_obj <- create_cor_mat(ACERtrc$TRACE_ts$corr_set_summary_ratiotest_corr_list_gaussian_kernel_hres_identity_arb_pollen_data)
corrplot::corrplot(corr_obj$cor_mat,type="upper",p.mat=1-corr_obj$signif_mat,pch.cex=0.75,col=clrrmp(200),
                   mar = c(0,0,1,0), number.cex = 0.5, title = 'TraCE TS')
dev.off()
