context("all dissimilarities")

test_that("dissimilarities give expected results", {
  otus <- matrix(nrow = 4, ncol = 3) 
  rownames(otus) <- paste(rep(paste("subj", 1:2, sep = ""), each = 2), 
                          rep(c("a","b"), 2), sep = "")
  metadata <- data.frame(subjID = rep(paste("subj", 1:2, sep = ""), each = 2), 
                         sampID = paste(rep(paste("subj", 1:2, sep = ""), each = 2), 
                                        rep(c("a","b"), 2), sep = ""), 
                         time = rep(1:2, 2)) 
  otus[1, ] <- c(0, 0.2, 0.8)
  otus[2, ] <- c(0.1, 0.3, 0.6)
  otus[3, ] <- c(0.4, 0.4, 0.2) 
  otus[4, ] <- c(0.2, 0.8, 0) 
  otus[5, ] <- c(0.2, 0, 0.8) 
  otus[6, ] <- c(0, 0.4, 0.6) 
  colnames(otus) <- paste("otu", 1:3, sep = "")
  set.seed(1); sim.tree <- rtree(3, tip.label = paste("otu", 1:3, sep = ""))
  
  paired.dat <- pltransform(otus, metadata, paired = TRUE)
  longit.dat <- pltransform(otus, metadata, paired = FALSE) 
  
  ## Bray-Curtis 
  bray.pb <- sum(abs(c(0.5, 0.5, 0.5)))/3
  bray.pq <- sum(abs(c(0.5+1/6, -0.5-1/6, -1/14+1/2)))/3 
  bray.lb <- sum(abs(c(1,1,1)))/3
  bray.lq <- sum(abs(c(2/3, 2/3, 6/7)))/3

  expect_equal(pldist(otus, metadata, paired = TRUE, binary = TRUE, method = "bray")$D[1,2], bray.pb)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = FALSE, method = "Bray")$D[1,2], bray.pq)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = TRUE, method = "br")$D[1,2], bray.lb)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = FALSE, method = "bray")$D[1,2], bray.lq)
  
  ## Jaccard 
  jac.pb <- 1 
  jac.pq <- 1 - (0 + 0 + 1/14)/(0.5 + 0.5 + 0.5)
  jac.lb <- 1 
  jac.lq <- 1 - (1/3 + 1/3 + 1/7) / 3
  
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = TRUE, method = "jaccard")$D[1,2], jac.pb)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = FALSE, method = "jac")$D[1,2], jac.pq)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = TRUE, method = "j")$D[1,2], jac.lb)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = FALSE, method = "Jacc")$D[1,2], jac.lq)
  
  ## Kulczynski 
  kul.pb <- 1
  kul.pq <- 1 - 0.5 * (1/sum(c(0.5, 0.5, 1/14)) + 1/sum(c(1/6,1/6,1/2))) * sum(c(0, 0, 1/14))
  kul.lb <- 1 
  kul.lq <- 1 - 0.5 * (1/(2 + 1/7) + 1/(1 + 2/3)) * sum(c(1/3, 1/3, 1/7))
  
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = TRUE, method = "kul")$D[1,2], kul.pb)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = FALSE, method = "kul")$D[1,2], kul.pq)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = TRUE, method = "kul")$D[1,2], kul.lb)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = FALSE, method = "kul")$D[1,2], kul.lq)
  
  
  ## Gower 
  gow.pb <- 
  gow.pq <- 
  gow.lb <- 
  gow.lq <- 
  
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = TRUE, method = "gower")$D[1,2], gow.pb)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = FALSE, method = "Gower")$D[1,2], gow.pq)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = TRUE, method = "gow")$D[1,2], gow.lb)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = FALSE, method = "Gow")$D[1,2], gow.lq)
  
  
  ## UniFrac 
  uf.pb <- sum(abs(paired.dat$tsf.data$dat.binary[1,] - paired.dat$tsf.data$dat.binary[2,]))/3
  uf.pq <- sum(abs(paired.dat$tsf.data$dat.quant[1,] - paired.dat$tsf.data$dat.quant[2,]))/3
  uf.lb <- sum(abs(longit.dat$tsf.data$dat.binary[1,] - longit.dat$tsf.data$dat.binary[2,]))/3
  uf.lq <- sum(abs(longit.dat$tsf.data$dat.quant[1,] - longit.dat$tsf.data$dat.quant[2,]))/3
  
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = TRUE, method = "unifrac")$D[1,2,"d_UW"], uf.pb)
  expect_equal(pldist(otus, metadata, paired = TRUE, binary = FALSE, method = "unifrac")$D[1,2,"d_1"], uf.pq)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = TRUE, method = "unifrac")$D[1,2,"d_UW"], uf.lb)
  expect_equal(pldist(otus, metadata, paired = FALSE, binary = FALSE, method = "unifrac")$D[1,2,"d_1"], uf.lq)
  
})