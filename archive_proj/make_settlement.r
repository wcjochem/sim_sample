#
# Create simulated "settlement strata" from population surface
# Similar to the Degree of Urbanisation model
#
#  Chris Jochem
#

require(raster)

make_settlement <- function(pop_raster){
  # create blocks of cells
  block_size <- rows <- cols <- 2 # 2x2 blocks
  nx <- dim(population_raster)[2]
  ny <- dim(population_raster)[1]
  
  blocks <- outer(1:ny, 1:nx, function(i,j) (i-1) %/% rows * ((ny+1) %/% cols) + (j-1) %/% cols + 1)
  block_raster <- raster(blocks, xmn=0, xmx=nx, ymn=0, ymx=ny)
  extent(block_raster) <- extent(population_raster)
  
  # population per block 
  block_pop <- as.data.frame(zonal(population_raster, block_raster, fun='sum'))
  block_pop_r <- reclassify(block_raster, block_pop)
  # smoothed 3x3 window
  block_pop_sm <- focal(block_pop_r, w=matrix(1/9,nrow=3,ncol=3))

  # classification thresholds - CHANGE HERE
  low_dens <- 50 # density to define clump
  sm_sett <- 250 # ppl per clump
  hi_dens <- 100
  lg_sett <- 800
  # find areas of density types
  # "low density"
  m <- matrix(c(0, low_dens, 0,
                low_dens, max(values(block_pop_sm),na.rm=T),1), 
              ncol=3, byrow=T)
  blocks_low <- reclassify(block_pop_sm, m)
  # make clumps
  low_clumps <- clump(blocks_low, directions=8, gaps=F)
    plot(low_clumps)
  # get population per clump
  clump_pop <- as.data.frame(zonal(population_raster, low_clumps, 'sum'))
  clump_pop$class <- ifelse(clump_pop$sum >= sm_sett, 1, NA)
  # reclassify to create small settlement areas
  small_sett <- reclassify(low_clumps, clump_pop[,c("zone","class")])
  small_sett[is.na(small_sett)] <- 0

  # "high density"
  m <- matrix(c(0, hi_dens, 0,
                hi_dens, max(values(block_pop_sm),na.rm=T),1), 
              ncol=3, byrow=T)
  blocks_hi <- reclassify(block_pop_sm, m)
  # make clumps
  hi_clumps <- clump(blocks_hi, directions=8, gaps=F)
  # get population per clump
  clump_pop <- as.data.frame(zonal(population_raster, hi_clumps, 'sum'))
  clump_pop$class <- ifelse(clump_pop$sum >= lg_sett, 2, NA)
  # reclassify to create small settlement areas
  lrg_sett <- reclassify(hi_clumps, clump_pop[,c("zone","class")])
  lrg_sett[is.na(lrg_sett)] <- 0
  
  return(max(small_sett, lrg_sett)) # raster
} 
# End
