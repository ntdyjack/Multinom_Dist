
#' sequence_cell
#' @description Generate count vector for a number of cells of a 
#' fixed subpopulation, given fixed number of reads and a probability mass.
#' @param ncells number of cells to generate count vectors for
#' @param reads_per_cell number of reads/UMIs per cell (assumed equal for all cells)
#' @param prob_dist probability mass (relative abundance of genes)
#' @param seed random seed
#'
#' @return The object returned is the MLE of probability mass 
#' (restrict to this for this study, as do not require actual 
#' count vectors for most analysis)
#' @export
#'
#' @examples
sequence_cell <- function(ncells, reads_per_cell, prob_dist,seed){
  #	draw samples from multinomial distribution 
#  set.seed(seed)
  obs_counts_cell <- rmultinom(n=ncells,prob=prob_dist,size=reads_per_cell)
# turn off library size scaling for now
#  estimate_prob_cell <- obs_counts_cell / reads_per_cell	
  return(obs_counts_cell)
#  return(t(estimate_prob_cell))
}

#' generate_probability_mass
#'
#' @param n_subtypes 
#' @param frac_hk_genes 
#' @param fold_changes 
#'
#' @return
#' @export
#'
#' @examples
generate_probability_mass <- function(n_subtypes, frac_hk_genes, fold_changes){
  #    Generate probability mass for a given cell/gene setup.
  
  #    Inputs:
  #    -- n_subtypes: number of cell subpopulations
  #    -- frac_hk_genes: (common) proportion of house-keeping genes
  #    -- fold_changes: list of fold-changes for marker-genes of length=n_subtypes
  
  #    Outputs:
  #    -- probability masses
  
  n_markers <- length(fold_changes)
  marker_density_sum <- 1 - frac_hk_genes #density left for marker genes to occupy
  
  pre_prob <- matrix(1,nrow=n_subtypes,ncol=n_subtypes)
  diag(pre_prob) <- fold_changes
  
  pre_densities <- solve(pre_prob)
  pre_densities <- (pre_densities %*% rep(1,n_markers)) * marker_density_sum
  
  almost_densities <- apply(pre_prob,1,function(x) x*pre_densities)
  prob_densities <- cbind(almost_densities,rep(frac_hk_genes,nrow(almost_densities)))
  
  return(prob_densities)
}


#' generate_probability_mass
#'
#' @param n_subtypes 
#' @param frac_hk_genes 
#' @param fold_changes 
#'
#' @return
#' @export
#'
#' @examples
generate_prob_mass <- function(n_subtypes, n_genes, gene_pcts, gene_lfcs){
  #    Generate probability mass for a given cell/gene setup.
  
  #    Inputs:
  #    -- n_subtypes: number of cell subpopulations
  #    -- n_genes: number of genes to simulate	
  #    -- gene_pcts: percent marker genes for each type, length=n_subtypes + 1
  #    -- gene_lfcs: fold changes of genes of each type, length=n_subtypes + 1 

  #    Outputs:
  #    -- probability masses
  if(sum(gene_pcts)!=1){stop("gene %s must sum to 1")}

  l <- length(gene_pcts) #last index for noise genes
#  marker_pct <- 1 -  gene_pcts[l] #pct left for marker genes to occupy
#  marker_sum <- n_genes * marker_pct# number of marker genes
  gene_sums <- ceiling(n_genes * gene_pcts)
#  noise_prob <-  gene_pcts[l] / gene_sums[l]
#  base_prob <- marker_pct /  (gene_lfcs[-l] + (l - 2))
#  base_prob <- base_prob / gene_sums[-l]
#  base_probs <- (1 / (gene_lfcs[-l] + 1)) / n_genes
#  marker_prob <- base_prob * gene_lfcs[i]
  idxs <- c(0,cumsum(gene_sums[-l])) 

  prob_densities <- t(sapply(1:(l-1), function(i) {
    base_prob <- 1 / ( gene_lfcs[i] * gene_sums[i] + sum(gene_sums[-i]) )
    tmp <- rep_len(base_prob,n_genes)
    tmp[(idxs[i]+1):(idxs[i+1])] <- rep_len(base_prob * gene_lfcs[i], gene_sums[i])
#    tmp[(tail(idxs,1)+1):(n_genes)] <- noise_prob
    return(tmp)
  }))  

  return(prob_densities)
}


#' log_likelihood_ratio_multinomial
#'
#' @param prob1 
#' @param prob2 
#' @param totcounts 
#'
#' @return
#' @export
#'
#' @examples
log_likelihood_ratio_multinomial <- function(counts1,counts2){
  #	Given count vectors for two cells, compute log-likelihood ratio comparing
  #    alternative vs null hypothesis, as described above.
  
  #    Inputs: 
  #    -- prob1: estimate of probability mass for cell 1
  #    -- prob2: estimate of probability mass for cell 2
  #    -- totcounts: total counts for each cell (assumed to be equal for both cells)
  
  #    Outputs:
  #    -- log-likelihood ratio 
  #  counts1 <- round(prob1*totcounts)
  #  counts2 <- round(prob2*totcounts) #rounding required for stability
  prob1 <- counts1/sum(counts1)
  prob2 <- counts2/sum(counts2)

  tmp_lhd_ha <- sum(log(prob1^counts1)) + sum(log(prob2^counts2))	
  prob_avg <- (prob1+prob2)/2
  tmp_lhd_h0 <- sum(log(prob_avg^counts1)) + sum(log(prob_avg^counts2))
  
  return(tmp_lhd_ha - tmp_lhd_h0)
}

#' compute_distance
#'
#' @param dat 
#' @param type 
#' @param totcounts 
#'
#' @return
#' @export
#'
#' @examples
compute_distance <- function(dat,type){
  #	Compute distance matrix from a probability data matrix
  
  #    Inputs:
  #    -- dat: probability data matrix (samples in rows)
  #    -- type: likelihood ratio test ("lrm") or euclidean ("euc")
  #    -- totcounts: total counts for each cell (assumed to be equal for both cells)
  
  #    Outputs:
  #    -- nxn distance matrix
  
  n <- ncol(dat)
  dist_mat <- matrix(0,nrow=n,ncol=n)
  dist_vec <- unlist(lapply(1:(n-1), function(i) {
    sapply((i+1):(n), function(j) {
      if(type=="lrm"){
        d <- log_likelihood_ratio_multinomial(dat[,i],dat[,j])
      }else if(type=="euc"){
        d <- sqrt(sum((dat[,i] - dat[,j])^2)) 
      }
      return(d)
    })
  }))
  
  dist_mat[lower.tri(dist_mat)] <- dist_vec
  dist_mat <- t(dist_mat)
  dist_mat[lower.tri(dist_mat)] <- dist_vec
  
  return(as.dist(dist_mat))
}

