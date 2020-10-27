# Andrea Sanchez-Tapia Santos Drumont CC Script
## This script contains the code developed by AST to run the niche equivalency and similarity tests on SDrumont computation cluster

# #!/usr/bin/env Rscript
# 
# install.packages("doParallel")
# 
# 
# load("./output/pca_grid.rdata")
# indexes <- expand.grid(1:length(pca.grid), 1:length(pca.grid))
# z1_list <- pca.grid[indexes$Var2[1:3]]
# z2_list <- pca.grid[indexes$Var1[1:3]]
# rm(pca.grid)
# message("cluster init...")
# library(parallel)
# library(doParallel)
# library(foreach)
# core <- detectCores()
# 
# #inicia el cluster
# cl <- parallel::makeCluster(4,
#                             outfile =
#                               paste("niche.log"), type = "SOCK")
# registerDoParallel(cl)
# message(paste("registercluster OK", getDoParWorkers()))
# 
# (ini = Sys.time())
# a <- foreach(lista1 = z1_list,
#              lista2 = z2_list,
#              .packages = c("ecospat")) %dopar%
#   ecospat::ecospat.niche.equivalency.test(z1 = lista1,
#                                           z2 = lista2,
#                                           rep = 100,
#                                           alternative = "lower")
# Sys.time() - ini
# stopCluster(cl)
# 