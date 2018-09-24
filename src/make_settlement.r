#
# create residential settlement types on the simulated grid data
# Code contributed by Claire Dooley
# 

make_settlement <- function(population_raster){
  # # select test pop  -- TESTING
  # population_raster <- countfields[[4]]
  #   plot(population_raster)
  
  # create blocks of cells
  block_mat <- matrix(0, nrow=regionsize, ncol=regionsize)
  block_size <- 2   # 2x2 blocks
  v1 <- seq(1, regionsize, by=block_size)
  block_id <- 0
  for(a in v1){
    for(b in v1){
      block_id <- block_id + 1
      block_mat[a:(a+(block_size-1)),b:(b+(block_size-1))] <- block_id
    }
  }
  block_raster <- raster(block_mat, xmn=0, xmx=50, ymn=0, ymx=50) 
    # plot(block_raster)
  
  # mean pop count per block 
  block_vec <- as.vector(t(block_mat))
  population_vec <- population_raster[]
  mean_blocks_vec <- rep(NA, length(block_vec))
  for(a in unique(block_vec)){
    mean_blocks_vec[which(block_vec==a)] <- mean(population_vec[which(block_vec==a)])
  }
  
  # urban, rural, unsettled raster for blocks 
  urban_threshold_pop <- 10  ### impoortant for determining how many 'urban' cells there are, therefore variation in resid type across whole grid
  
  urban_rural_raster <- raster(nrows=regionsize, ncols=regionsize, xmn=0, xmx=50, ymn=0, ymx=50)
  urban_rural_raster[which(mean_blocks_vec==0)] <- 0
  urban_rural_raster[which(mean_blocks_vec>0 & mean_blocks_vec<=urban_threshold_pop)] <- 1
  urban_rural_raster[which(mean_blocks_vec>urban_threshold_pop)] <- 2
    # plot(urban_rural_raster)
  
  #urban centres
  urban_only_raster <- urban_rural_raster
  urban_only_raster[!urban_only_raster==2] <- NA
  urban_clumps <- clump(urban_only_raster, directions = 4)
    # plot(urban_clumps)
  urban_only <- data.frame(rasterToPoints(urban_clumps))
  urban_centres <- sapply(split(urban_only[,c("x", "y")], urban_only$clumps), colMeans)
  
  # block centres
  block_df <- data.frame(rasterToPoints(block_raster))
  block_centres <- sapply(split(block_df[, c("x", "y")], block_df$layer), colMeans)
  
  ### CREATE RESIDENTIAL TYPE BY BLOCKS LAYER ###
  # define parameters
  small_sett_block_limit <- 30   # below & equal to this number of blocks (not cells!) in urban clump is small sett, above is large sett
  # thresholds for high/low pop
  pop_threshold_citycentre <- 40 # for city centre of large setts
  pop_threshold_edgeorsmall <- 20 # for small setts and edge of large setts 
  # threshold distance from centre of large settlements, at which resid types change
  dist_thrshold <- 2
  resid_types <- 2:5  # for urban only here, rural=1; unsettled=0
  
  resid_vec <- urban_rural_raster[]
  cells_per_clump <- table(urban_clumps[])
  sett_type <- NULL
  
  # loop through urban clumps as will only be a few
  for(i in colnames(urban_centres)){
      clump_cells <- cells_per_clump[names(cells_per_clump)==i]
      clump_blocks <- clump_cells/(block_size^2)  # number of blocks within clump
      blocks_in_clump <- unique(na.omit(block_vec[urban_clumps[] == i]))
      if(clump_blocks <= small_sett_block_limit){
          for(j in blocks_in_clump){
              if(mean_blocks_vec[which(block_vec == j)][1] <= pop_threshold_edgeorsmall){
                    resid_vec[which(block_vec == j)] <- resid_types[1]
              }else{
                    resid_vec[which(block_vec == j)] <- resid_types[2]
              }
      }
      }else{ 
          clump_centre <- urban_centres[,colnames(urban_centres)==i]
          for(j in blocks_in_clump){
              j_block_centre <- block_centres[,colnames(block_centres) == j]
              if(dist(rbind(j_block_centre, clump_centre)) <= dist_thrshold){
                  if(mean_blocks_vec[which(block_vec == j)][1] <= pop_threshold_edgeorsmall){
                    resid_vec[which(block_vec == j)] <- resid_types[1]
                  }else{
                    resid_vec[which(block_vec == j)] <- resid_types[2]
                  }
              }else{
                  if(mean_blocks_vec[which(block_vec == j)][1] <= pop_threshold_citycentre){
                    resid_vec[which(block_vec == j)] <- resid_types[3]
                  }else{
                    resid_vec[which(block_vec == j)] <- resid_types[4]
                  }
              }
          }
      }
  }     
  
  # create raster for resid types
  resid_mat <- matrix(resid_vec, nrow = regionsize, byrow=T)
  resid_raster <- raster(resid_mat, xmn=0, xmx=50, ymn=0, ymx=50)
    # plot(resid_raster)
  return(resid_raster)
}


