#reproducible example for thinTrackTime issue

library(move)
library(dplyr)


load("sample_gannet_track.RData")

track_th <- x %>%
  thinTrackTime(interval = as.difftime(60, units='mins'),
                tolerance = as.difftime(30, units='mins')) 

#Error: C stack usage  7969860 is too close to the limit


# > sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Pop!_OS 20.10
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] dplyr_1.0.5      move_4.0.6       rgdal_1.5-23     raster_3.4-5     sp_1.4-5         geosphere_1.5-10
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.6       xml2_1.3.2       magrittr_2.0.1   tidyselect_1.1.0 lattice_0.20-41  R6_2.5.0         rlang_0.4.10     fastmap_1.1.0    fansi_0.4.2      httr_1.4.2      
# [11] tools_4.0.3      parallel_4.0.3   grid_4.0.3       packrat_0.6.0    utf8_1.2.1       DBI_1.1.1        ellipsis_0.3.1   digest_0.6.27    assertthat_0.2.1 tibble_3.1.0    
# [21] lifecycle_1.0.0  crayon_1.4.1     purrr_0.3.4      vctrs_0.3.7      codetools_0.2-18 memoise_2.0.0    glue_1.4.2       cachem_1.0.4     compiler_4.0.3   pillar_1.5.1    
# [31] generics_0.1.0   pkgconfig_2.0.3 