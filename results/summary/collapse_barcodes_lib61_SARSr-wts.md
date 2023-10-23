Collapse barcodes to final per-RBD/mutant phenotype scores for the
wildtype sarbecoviruses pool
================
Tyler Starr
9/27/2023

- <a href="#setup" id="toc-setup">Setup</a>
- <a href="#calculate-per-variant-mean-scores"
  id="toc-calculate-per-variant-mean-scores">Calculate per-variant mean
  scores</a>
- <a href="#heatmaps" id="toc-heatmaps">Heatmaps!</a>

This notebook reads in the per-barcode sera binding values and
previously measured expression for sarbecovirus homologs pool. It
synthesizes these two sets of results and calculates the final ‘mean’
phenotypes for each variant, and generates some coverage and QC
analyses.

``` r
#list of packages to install/load
packages = c("yaml","data.table","tidyverse","gridExtra")
#install any packages not already installed
installed_packages <- packages %in% rownames(installed.packages())
if(any(installed_packages == F)){
  install.packages(packages[!installed_packages],
                   lib=c(paste("/uufs/chpc.utah.edu/common/home/",Sys.getenv("USER"),"/RLibs/",Sys.getenv("R_VERSION"),sep="")),
                   repos=c("http://cran.us.r-project.org"))
}
#load packages
invisible(lapply(packages, library, character.only=T))

knitr::opts_chunk$set(echo = T)
knitr::opts_chunk$set(dev.args = list(png = list(type = "cairo")))

#read in config file
config <- read_yaml("config.yaml")

#make output directory
if(!file.exists(config$final_variant_scores_dir)){
  dir.create(file.path(config$final_variant_scores_dir))
}
```

Session info for reproducing environment:

``` r
sessionInfo()
```

    ## R version 4.1.3 (2022-03-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Rocky Linux 8.8 (Green Obsidian)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /uufs/chpc.utah.edu/sys/spack/linux-rocky8-nehalem/gcc-8.5.0/intel-oneapi-mkl-2021.4.0-h43nkmwzvaltaa6ii5l7n6e7ruvjbmnv/mkl/2021.4.0/lib/intel64/libmkl_rt.so.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] gridExtra_2.3     forcats_0.5.1     stringr_1.4.0     dplyr_1.0.8      
    ##  [5] purrr_0.3.4       readr_2.1.2       tidyr_1.2.0       tibble_3.1.6     
    ##  [9] ggplot2_3.4.1     tidyverse_1.3.1   data.table_1.14.2 yaml_2.3.5       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.1.2 xfun_0.30        haven_2.4.3      colorspace_2.0-3
    ##  [5] vctrs_0.5.2      generics_0.1.2   htmltools_0.5.2  utf8_1.2.2      
    ##  [9] rlang_1.0.6      pillar_1.7.0     glue_1.6.2       withr_2.5.0     
    ## [13] DBI_1.1.2        dbplyr_2.1.1     modelr_0.1.8     readxl_1.3.1    
    ## [17] lifecycle_1.0.3  munsell_0.5.0    gtable_0.3.0     cellranger_1.1.0
    ## [21] rvest_1.0.2      evaluate_0.15    knitr_1.37       tzdb_0.2.0      
    ## [25] fastmap_1.1.0    fansi_1.0.2      broom_0.7.12     Rcpp_1.0.11     
    ## [29] backports_1.4.1  scales_1.2.1     jsonlite_1.8.7   fs_1.5.2        
    ## [33] hms_1.1.1        digest_0.6.29    stringi_1.7.6    grid_4.1.3      
    ## [37] cli_3.6.0        tools_4.1.3      magrittr_2.0.2   crayon_1.5.0    
    ## [41] pkgconfig_2.0.3  ellipsis_0.3.2   xml2_1.3.3       reprex_2.0.1    
    ## [45] lubridate_1.8.0  rstudioapi_0.13  assertthat_0.2.1 rmarkdown_2.13  
    ## [49] httr_1.4.7       R6_2.5.1         compiler_4.1.3

## Setup

Read in tables of per-barcode AUC

``` r
dt <- data.table(read.csv(config$mAb_EC50_file),stringsAsFactors=F)[library=="lib61_SARSr-wts",]

#check all targets are in the targets_ordered config list
unique(dt$target) %in% config$targets_ordered
```

    ##   [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [31] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [46] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [61] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [76] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ##  [91] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

``` r
#assign target as a factor in my desired overall plotting order
dt[,target := factor(dt$target,levels=config$targets_ordered)]

#read in previously measured expression measurements for these variants? Note, don't have this metric for the sarbs that have been pooled in since the original library pool back in 2021
dt_expr <- data.table(read.csv(config$SARSr_lib61_mut_bind_expr),stringsAsFactors=F)
dt_expr[target=="SARS-CoV-2",target:="SARS-CoV-2_WH1"]
```

## Calculate per-variant mean scores

Unfiltered, look at distribution of AUC scores

``` r
p1 <- ggplot(dt[!is.na(EC50_S2V29_v37_2),],aes(x=target,y=log10(EC50_S2V29_v37_2)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S2V29_v37_2 mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p2 <- ggplot(dt[!is.na(EC50_S2V29),],aes(x=target,y=log10(EC50_S2V29)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S2V29 mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p3 <- ggplot(dt[!is.na(EC50_S2V29_UCA),],aes(x=target,y=log10(EC50_S2V29_UCA)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S2V29_UCA mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p4 <- ggplot(dt[!is.na(EC50_S2K146_v860),],aes(x=target,y=log10(EC50_S2K146_v860)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S2K146_v860 mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p5 <- ggplot(dt[!is.na(EC50_S2K146),],aes(x=target,y=log10(EC50_S2K146)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S2K146 mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p6 <- ggplot(dt[!is.na(EC50_S2K146_UCA),],aes(x=target,y=log10(EC50_S2K146_UCA)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S2K146_UCA mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p7 <- ggplot(dt[!is.na(EC50_S3L17),],aes(x=target,y=log10(EC50_S3L17)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S3L17 mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p8 <- ggplot(dt[!is.na(EC50_S3L17_UCA),],aes(x=target,y=log10(EC50_S3L17_UCA)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S3L17_UCA mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p9 <- ggplot(dt[!is.na(EC50_S309_887),],aes(x=target,y=log10(EC50_S309_887)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S309_887 mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=1)
```

    ## Warning: Groups with fewer than two data points have been dropped.
    ## Groups with fewer than two data points have been dropped.
    ## Groups with fewer than two data points have been dropped.

<img src="collapse_barcodes_lib61_SARSr-wts_files/figure-gfm/unfiltered_AUCs-1.png" style="display: block; margin: auto;" />

Let’s add a variable that flags the top and bottom 5% of expression
scores for each variant, and see how violin plots look when censoring
these top and bottom 10% of expressed barcodes. Let’s *only* do this for
targets that are in the many-barcode pool, not the spike-in targets that
had between 2-5 barcodes total (but with much higher per-barcode cell
counts for the titraitons!)

``` r
dt[,`S2V29_v37_2_censor_lower`:=quantile(EC50_S2V29_v37_2,0.05,na.rm=T,type=7),by=c("library","target")]
dt[,`S2V29_v37_2_censor_upper`:=quantile(EC50_S2V29_v37_2,0.95,na.rm=T,type=7),by=c("library","target")]

dt[,`S2V29_censor_lower`:=quantile(EC50_S2V29,0.05,na.rm=T,type=7),by=c("library","target")]
dt[,`S2V29_censor_upper`:=quantile(EC50_S2V29,0.95,na.rm=T,type=7),by=c("library","target")]

dt[,`S2V29_UCA_censor_lower`:=quantile(EC50_S2V29_UCA,0.05,na.rm=T,type=7),by=c("library","target")]
dt[,`S2V29_UCA_censor_upper`:=quantile(EC50_S2V29_UCA,0.95,na.rm=T,type=7),by=c("library","target")]

dt[,`S2K146_v860_censor_lower`:=quantile(EC50_S2K146_v860,0.05,na.rm=T,type=7),by=c("library","target")]
dt[,`S2K146_v860_censor_upper`:=quantile(EC50_S2K146_v860,0.95,na.rm=T,type=7),by=c("library","target")]

dt[,`S2K146_censor_lower`:=quantile(EC50_S2K146,0.05,na.rm=T,type=7),by=c("library","target")]
dt[,`S2K146_censor_upper`:=quantile(EC50_S2K146,0.95,na.rm=T,type=7),by=c("library","target")]

dt[,`S2K146_UCA_censor_lower`:=quantile(EC50_S2K146_UCA,0.05,na.rm=T,type=7),by=c("library","target")]
dt[,`S2K146_UCA_censor_upper`:=quantile(EC50_S2K146_UCA,0.95,na.rm=T,type=7),by=c("library","target")]

dt[,`S3L17_censor_lower`:=quantile(EC50_S3L17,0.05,na.rm=T,type=7),by=c("library","target")]
dt[,`S3L17_censor_upper`:=quantile(EC50_S3L17,0.95,na.rm=T,type=7),by=c("library","target")]

dt[,`S3L17_UCA_censor_lower`:=quantile(EC50_S3L17_UCA,0.05,na.rm=T,type=7),by=c("library","target")]
dt[,`S3L17_UCA_censor_upper`:=quantile(EC50_S3L17_UCA,0.95,na.rm=T,type=7),by=c("library","target")]

dt[,`S309_887_censor_lower`:=quantile(EC50_S309_887,0.05,na.rm=T,type=7),by=c("library","target")]
dt[,`S309_887_censor_upper`:=quantile(EC50_S309_887,0.95,na.rm=T,type=7),by=c("library","target")]

p1 <- ggplot(dt[!is.na(EC50_S2V29_v37_2) & ((EC50_S2V29_v37_2 > S2V29_v37_2_censor_lower & EC50_S2V29_v37_2 < S2V29_v37_2_censor_upper) | target %in% config$targets_low_bc) ,],
             aes(x=target,y=log10(EC50_S2V29_v37_2)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S2V29_v37_2 mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p2 <- ggplot(dt[!is.na(EC50_S2V29) & ((EC50_S2V29 > S2V29_censor_lower & EC50_S2V29 < S2V29_censor_upper) | target %in% config$targets_low_bc),],
             aes(x=target,y=log10(EC50_S2V29)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S2V29 mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p3 <- ggplot(dt[!is.na(EC50_S2V29_UCA) & ((EC50_S2V29_UCA > S2V29_UCA_censor_lower & EC50_S2V29_UCA < S2V29_UCA_censor_upper) | target %in% config$targets_low_bc),],
             aes(x=target,y=log10(EC50_S2V29_UCA)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S2V29_UCA mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p4 <- ggplot(dt[!is.na(EC50_S2K146_v860) & ((EC50_S2K146_v860 > S2K146_v860_censor_lower & EC50_S2K146_v860 < S2K146_v860_censor_upper) | target %in% config$targets_low_bc),],
             aes(x=target,y=log10(EC50_S2K146_v860)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S2K146_v860 mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p5 <- ggplot(dt[!is.na(EC50_S2K146) & ((EC50_S2K146 > S2K146_censor_lower & EC50_S2K146 < S2K146_censor_upper) | target %in% config$targets_low_bc),],
             aes(x=target,y=log10(EC50_S2K146)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S2K146 mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p6 <- ggplot(dt[!is.na(EC50_S2K146_UCA) & ((EC50_S2K146_UCA > S2K146_UCA_censor_lower & EC50_S2K146_UCA < S2K146_UCA_censor_upper) | target %in% config$targets_low_bc),],
             aes(x=target,y=log10(EC50_S2K146_UCA)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S2K146_UCA mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p7 <- ggplot(dt[!is.na(EC50_S3L17) & ((EC50_S3L17 > S3L17_censor_lower & EC50_S3L17 < S3L17_censor_upper) | target %in% config$targets_low_bc),],
             aes(x=target,y=log10(EC50_S3L17)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S3L17 mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p8 <- ggplot(dt[!is.na(EC50_S3L17_UCA) & ((EC50_S3L17_UCA > S3L17_UCA_censor_lower & EC50_S3L17_UCA < S3L17_UCA_censor_upper) | target %in% config$targets_low_bc),],
             aes(x=target,y=log10(EC50_S3L17_UCA)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S3L17_UCA mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))

p9 <- ggplot(dt[!is.na(EC50_S309_887) & ((EC50_S309_887 > S309_887_censor_lower & EC50_S309_887 < S309_887_censor_upper) | target %in% config$targets_low_bc),],
             aes(x=target,y=log10(EC50_S309_887)))+
  geom_violin(scale="width")+stat_summary(fun=median,geom="point",size=1)+
  ggtitle("S309_887 mAb binding")+xlab("homolog")+theme(axis.text.x=element_text(angle=-90,hjust=0))


grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=1)
```

    ## Warning: Groups with fewer than two data points have been dropped.
    ## Groups with fewer than two data points have been dropped.

<img src="collapse_barcodes_lib61_SARSr-wts_files/figure-gfm/censor_5_AUCs-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/lib61_vioplots_AUC-cens.pdf",sep="")))
```

Calculate the mean (geometric mean) per variant, the standard deviation,
and the number of (post-filter) barcodes on which a variant score was
determined

``` r
#apply the censors to NA out the phenotypes outside the range for the targets that aren't in the 2-5 higher-count barcode spike in list
dt[!(target %in% config$targets_low_bc) & (EC50_S2V29_v37_2 < `S2V29_v37_2_censor_lower` | EC50_S2V29_v37_2 > `S2V29_v37_2_censor_upper`), EC50_S2V29_v37_2:=NA]
dt[!(target %in% config$targets_low_bc) & (EC50_S2V29 < `S2V29_censor_lower` | EC50_S2V29 > `S2V29_censor_upper`), EC50_S2V29:=NA]
dt[!(target %in% config$targets_low_bc) & (EC50_S2V29_UCA < `S2V29_UCA_censor_lower` | EC50_S2V29_UCA > `S2V29_UCA_censor_upper`), EC50_S2V29_UCA:=NA]
dt[!(target %in% config$targets_low_bc) & (EC50_S2K146_v860 < `S2K146_v860_censor_lower` | EC50_S2K146_v860 > `S2K146_v860_censor_upper`), EC50_S2K146_v860:=NA]
dt[!(target %in% config$targets_low_bc) & (EC50_S2K146 < `S2K146_censor_lower` | EC50_S2K146 > `S2K146_censor_upper`), EC50_S2K146:=NA]
dt[!(target %in% config$targets_low_bc) & (EC50_S2K146_UCA < `S2K146_UCA_censor_lower` | EC50_S2K146_UCA > `S2K146_UCA_censor_upper`), EC50_S2K146_UCA:=NA]
dt[!(target %in% config$targets_low_bc) & (EC50_S3L17 < `S3L17_censor_lower` | EC50_S3L17 > `S3L17_censor_upper`), EC50_S3L17:=NA]
dt[!(target %in% config$targets_low_bc) & (EC50_S3L17_UCA < `S3L17_UCA_censor_lower` | EC50_S3L17_UCA > `S3L17_UCA_censor_upper`), EC50_S3L17_UCA:=NA]
dt[!(target %in% config$targets_low_bc) & (EC50_S309_887 < `S309_887_censor_lower` | EC50_S309_887 > `S309_887_censor_upper`), EC50_S309_887:=NA]

#calculate geomean, sd, and # bcs collapsed per genotype
dt[,mean_EC50_S2V29_v37_2:=10^mean(log10(EC50_S2V29_v37_2),na.rm=T),by=c("library","target")]
dt[,sd_EC50_S2V29_v37_2:=sd(log10(EC50_S2V29_v37_2),na.rm=T),by=c("library","target")]
dt[,n_bc_EC50_S2V29_v37_2:=sum(!is.na(EC50_S2V29_v37_2)),by=c("library","target")]

dt[,mean_EC50_S2V29:=10^mean(log10(EC50_S2V29),na.rm=T),by=c("library","target")]
dt[,sd_EC50_S2V29:=sd(log10(EC50_S2V29),na.rm=T),by=c("library","target")]
dt[,n_bc_EC50_S2V29:=sum(!is.na(EC50_S2V29)),by=c("library","target")]

dt[,mean_EC50_S2V29_UCA:=10^mean(log10(EC50_S2V29_UCA),na.rm=T),by=c("library","target")]
dt[,sd_EC50_S2V29_UCA:=sd(log10(EC50_S2V29_UCA),na.rm=T),by=c("library","target")]
dt[,n_bc_EC50_S2V29_UCA:=sum(!is.na(EC50_S2V29_UCA)),by=c("library","target")]

dt[,mean_EC50_S2K146_v860:=10^mean(log10(EC50_S2K146_v860),na.rm=T),by=c("library","target")]
dt[,sd_EC50_S2K146_v860:=sd(log10(EC50_S2K146_v860),na.rm=T),by=c("library","target")]
dt[,n_bc_EC50_S2K146_v860:=sum(!is.na(EC50_S2K146_v860)),by=c("library","target")]

dt[,mean_EC50_S2K146:=10^mean(log10(EC50_S2K146),na.rm=T),by=c("library","target")]
dt[,sd_EC50_S2K146:=sd(log10(EC50_S2K146),na.rm=T),by=c("library","target")]
dt[,n_bc_EC50_S2K146:=sum(!is.na(EC50_S2K146)),by=c("library","target")]

dt[,mean_EC50_S2K146_UCA:=10^mean(log10(EC50_S2K146_UCA),na.rm=T),by=c("library","target")]
dt[,sd_EC50_S2K146_UCA:=sd(log10(EC50_S2K146_UCA),na.rm=T),by=c("library","target")]
dt[,n_bc_EC50_S2K146_UCA:=sum(!is.na(EC50_S2K146_UCA)),by=c("library","target")]

dt[,mean_EC50_S3L17:=10^mean(log10(EC50_S3L17),na.rm=T),by=c("library","target")]
dt[,sd_EC50_S3L17:=sd(log10(EC50_S3L17),na.rm=T),by=c("library","target")]
dt[,n_bc_EC50_S3L17:=sum(!is.na(EC50_S3L17)),by=c("library","target")]

dt[,mean_EC50_S3L17_UCA:=10^mean(log10(EC50_S3L17_UCA),na.rm=T),by=c("library","target")]
dt[,sd_EC50_S3L17_UCA:=sd(log10(EC50_S3L17_UCA),na.rm=T),by=c("library","target")]
dt[,n_bc_EC50_S3L17_UCA:=sum(!is.na(EC50_S3L17_UCA)),by=c("library","target")]

dt[,mean_EC50_S309_887:=10^mean(log10(EC50_S309_887),na.rm=T),by=c("library","target")]
dt[,sd_EC50_S309_887:=sd(log10(EC50_S309_887),na.rm=T),by=c("library","target")]
dt[,n_bc_EC50_S309_887:=sum(!is.na(EC50_S309_887)),by=c("library","target")]
```

Collapse down to tables reporting just the summary statistics for each
genotype.

``` r
dt_final <- dt[,.(library,target,variant_class,
                  mean_EC50_S2V29_v37_2, sd_EC50_S2V29_v37_2, n_bc_EC50_S2V29_v37_2,
                  mean_EC50_S2V29, sd_EC50_S2V29, n_bc_EC50_S2V29,
                  mean_EC50_S2V29_UCA, sd_EC50_S2V29_UCA, n_bc_EC50_S2V29_UCA,
                  mean_EC50_S2K146_v860, sd_EC50_S2K146_v860, n_bc_EC50_S2K146_v860,
                  mean_EC50_S2K146, sd_EC50_S2K146, n_bc_EC50_S2K146,
                  mean_EC50_S2K146_UCA, sd_EC50_S2K146_UCA, n_bc_EC50_S2K146_UCA,
                  mean_EC50_S3L17, sd_EC50_S3L17, n_bc_EC50_S3L17,
                  mean_EC50_S3L17_UCA, sd_EC50_S3L17_UCA, n_bc_EC50_S3L17_UCA,
                  mean_EC50_S309_887, sd_EC50_S309_887, n_bc_EC50_S309_887)]

dt_final <- unique(dt_final); setkey(dt_final, target)
```

Let’s look how SEM is distributed. Can see that SEM is generally very,
very low. Also that it doesn’t really have a relationship with the EC50
metric, which is good.

``` r
par(mfrow=c(3,3))
#S2V29_v37_2
x <- log10(dt_final[,mean_EC50_S2V29_v37_2]); y <- dt_final[,sd_EC50_S2V29_v37_2/sqrt(n_bc_EC50_S2V29_v37_2)]; plot(x,y,pch=16,col="#00000090",xlab="EC50",ylab="SEM",main="S2V29_v37_2")

#S2V29
x <- log10(dt_final[,mean_EC50_S2V29]); y <- dt_final[,sd_EC50_S2V29/sqrt(n_bc_EC50_S2V29)]; plot(x,y,pch=16,col="#00000090",xlab="EC50",ylab="SEM",main="S2V29")

#S2V29_UCA
x <- log10(dt_final[,mean_EC50_S2V29_UCA]); y <- dt_final[,sd_EC50_S2V29_UCA/sqrt(n_bc_EC50_S2V29_UCA)]; plot(x,y,pch=16,col="#00000090",xlab="EC50",ylab="SEM",main="S2V29_UCA")

#S2K146_v860
x <- log10(dt_final[,mean_EC50_S2K146_v860]); y <- dt_final[,sd_EC50_S2K146_v860/sqrt(n_bc_EC50_S2K146_v860)]; plot(x,y,pch=16,col="#00000090",xlab="EC50",ylab="SEM",main="S2K146_v860")

#S2K146
x <- log10(dt_final[,mean_EC50_S2K146]); y <- dt_final[,sd_EC50_S2K146/sqrt(n_bc_EC50_S2K146)]; plot(x,y,pch=16,col="#00000090",xlab="EC50",ylab="SEM",main="S2K146")

#S2K146_UCA
x <- log10(dt_final[,mean_EC50_S2K146_UCA]); y <- dt_final[,sd_EC50_S2K146_UCA/sqrt(n_bc_EC50_S2K146_UCA)]; plot(x,y,pch=16,col="#00000090",xlab="EC50",ylab="SEM",main="S2K146_UCA")

#S3L17
x <- log10(dt_final[,mean_EC50_S3L17]); y <- dt_final[,sd_EC50_S3L17/sqrt(n_bc_EC50_S3L17)]; plot(x,y,pch=16,col="#00000090",xlab="EC50",ylab="SEM",main="S3L17")

#S3L17_UCA
x <- log10(dt_final[,mean_EC50_S3L17_UCA]); y <- dt_final[,sd_EC50_S3L17_UCA/sqrt(n_bc_EC50_S3L17_UCA)]; plot(x,y,pch=16,col="#00000090",xlab="EC50",ylab="SEM",main="S3L17_UCA")

#S309_887
x <- log10(dt_final[,mean_EC50_S309_887]); y <- dt_final[,sd_EC50_S309_887/sqrt(n_bc_EC50_S309_887)]; plot(x,y,pch=16,col="#00000090",xlab="EC50",ylab="SEM",main="S309_887")
```

<img src="collapse_barcodes_lib61_SARSr-wts_files/figure-gfm/plot_SEMs-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/SEM-v-EC50.pdf",sep=""),useDingbats=F))
```

Add in the previously-measured expression values.

``` r
dt_final[,expr:=as.numeric(NA)]
for(i in 1:nrow(dt_final)){
  bg <- as.character(dt_final[i,target])
  if(bg %in% dt_expr$target){
    dt_final[i,expr := dt_expr[target==bg,expression]]
  }
}

#compute a delta_expr relative to the median
dt_final[,delta_expr := expr - median(dt_final$expr,na.rm=T)]
```

Filter out the two backgrounds that were completely non-expressing. Most
barcodes were purged before the affinity measurements for these
backgrounds, so the affinities are determined from few barcodes and are
just generally unreliable because these are poorly folded/expressing
variants. (E.g. could see very high standard deviations)

``` r
dt_final[target %in% c("HKU3-8","AncSARS1a_alt"),c("mean_EC50_S2V29_v37_2","sd_EC50_S2V29_v37_2","n_bc_EC50_S2V29_v37_2",
                                                   "mean_EC50_S2V29","sd_EC50_S2V29","n_bc_EC50_S2V29",
                                                   "mean_EC50_S2V29_UCA","sd_EC50_S2V29_UCA","n_bc_EC50_S2V29_UCA",
                                                   "mean_EC50_S2K146_v860","sd_EC50_S2K146_v860","n_bc_EC50_S2K146_v860",
                                                   "mean_EC50_S2K146","sd_EC50_S2K146","n_bc_EC50_S2K146",
                                                   "mean_EC50_S2K146_UCA","sd_EC50_S2K146_UCA","n_bc_EC50_S2K146_UCA",
                                                   "mean_EC50_S3L17","sd_EC50_S3L17","n_bc_EC50_S3L17",
                                                   "mean_EC50_S3L17_UCA","sd_EC50_S3L17_UCA","n_bc_EC50_S3L17_UCA",
                                                   "mean_EC50_S309_887","sd_EC50_S309_887","n_bc_EC50_S309_887"):=NA]
```

Coverage stats on n_barcodes for different measurements in the final
pooled measurements.

``` r
par(mfrow=c(3,3))

hist(dt_final$n_bc_EC50_S2V29_v37_2,col="gray50",main=paste("S2V29_v37_2,\nmedian ",median(dt_final$n_bc_EC50_S2V29_v37_2,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)

hist(dt_final$n_bc_EC50_S2V29,col="gray50",main=paste("S2V29,\nmedian ",median(dt_final$n_bc_EC50_S2V29,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)

hist(dt_final$n_bc_EC50_S2V29_UCA,col="gray50",main=paste("S2V29_UCA,\nmedian ",median(dt_final$n_bc_EC50_S2V29_UCA,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)

hist(dt_final$n_bc_EC50_S2K146_v860,col="gray50",main=paste("S2K146_v860,\nmedian ",median(dt_final$n_bc_EC50_S2K146_v860,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)

hist(dt_final$n_bc_EC50_S2K146,col="gray50",main=paste("S2K146,\nmedian ",median(dt_final$n_bc_EC50_S2K146,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)

hist(dt_final$n_bc_EC50_S2K146_UCA,col="gray50",main=paste("S2K146_UCA,\nmedian ",median(dt_final$n_bc_EC50_S2K146_UCA,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)

hist(dt_final$n_bc_EC50_S3L17,col="gray50",main=paste("S3L17,\nmedian ",median(dt_final$n_bc_EC50_S3L17,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)

hist(dt_final$n_bc_EC50_S3L17_UCA,col="gray50",main=paste("S3L17_UCA,\nmedian ",median(dt_final$n_bc_EC50_S3L17_UCA,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)

hist(dt_final$n_bc_EC50_S309_887,col="gray50",main=paste("S309_887,\nmedian ",median(dt_final$n_bc_EC50_S309_887,na.rm=T),sep=""),xlab="number barcodes",ylab="number genotypes",breaks=20)
```

<img src="collapse_barcodes_lib61_SARSr-wts_files/figure-gfm/n_barcode_plots-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/lib61_hist_n_barcodes.pdf",sep="")))
```

Order factor variables for plotting

``` r
#order target by order given in config
dt_final$target <- factor(dt_final$target,levels=config$targets_ordered)

#rename some columns for convenience
setnames(dt_final,"mean_EC50_S2V29_v37_2","EC50_S2V29_v37_2")
setnames(dt_final,"mean_EC50_S2V29","EC50_S2V29")
setnames(dt_final,"mean_EC50_S2V29_UCA","EC50_S2V29_UCA")
setnames(dt_final,"mean_EC50_S2K146_v860","EC50_S2K146_v860")
setnames(dt_final,"mean_EC50_S2K146","EC50_S2K146")
setnames(dt_final,"mean_EC50_S2K146_UCA","EC50_S2K146_UCA")
setnames(dt_final,"mean_EC50_S3L17","EC50_S3L17")
setnames(dt_final,"mean_EC50_S3L17_UCA","EC50_S3L17_UCA")
setnames(dt_final,"mean_EC50_S309_887","EC50_S309_887")
```

## Heatmaps!

Output heatmaps illustrating all wildtype variants with separate columns
for each mAb.

``` r
#make temp long-form data frame
temp1 <- data.table::melt(dt_final[,.(target, EC50_S309_887, EC50_S3L17_UCA, EC50_S3L17, EC50_S2K146_UCA,EC50_S2K146,EC50_S2K146_v860,EC50_S2V29_UCA,EC50_S2V29 ,EC50_S2V29_v37_2)],
                          id.vars=c("target"),
                          measure.vars=c("EC50_S309_887", "EC50_S3L17_UCA", "EC50_S3L17", "EC50_S2K146_UCA","EC50_S2K146","EC50_S2K146_v860","EC50_S2V29_UCA","EC50_S2V29", "EC50_S2V29_v37_2"),
                          variable.name="mAb",value.name="EC50")

temp1$log10_EC50 <- log10(temp1$EC50)

p1 <- ggplot(temp1,aes(target,mAb))+geom_tile(aes(fill=log10_EC50),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#003366","#FFFFFF","#FFFFFF"),limits=c(-2,5),values=c(0,7/7,7/7),na.value="gray40")+
  #scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(5,12),values=c(0,1),na.value="yellow")+
  #scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,430,by=5)))+
  labs(x="RBD homolog",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold"))

grid.arrange(p1,nrow=1)
```

    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm):
    ## collapsing to unique 'x' values

    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm):
    ## collapsing to unique 'x' values

    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm):
    ## collapsing to unique 'x' values

<img src="collapse_barcodes_lib61_SARSr-wts_files/figure-gfm/heatmap_wildtypes_all_EC50-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/lib61_heatmap_EC50_all_wildtypes.pdf",sep="")))
```

Showing just extant sarbs.

``` r
#make temp long-form data frame
extant <- c(config$EurAf_extant,config$RsYN04_extant,config$SARS2_extant,config$SARS1_extant,config$Clade2_extant)

temp2 <- temp1[target %in% extant,];temp2$target <- factor(temp2$target,levels=extant)

p1 <- ggplot(temp2,aes(target,mAb))+geom_tile(aes(fill=log10_EC50),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#003366","#FFFFFF","#FFFFFF"),limits=c(-2,5),values=c(0,7/7,7/7),na.value="gray40")+
  #scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(5,12),values=c(0,1),na.value="yellow")+
  #scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,430,by=5)))+
  labs(x="RBD homolog",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold"))

p1
```

    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm):
    ## collapsing to unique 'x' values

    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm):
    ## collapsing to unique 'x' values

    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm):
    ## collapsing to unique 'x' values

<img src="collapse_barcodes_lib61_SARSr-wts_files/figure-gfm/heatmap_wildtypes_extants_EC50-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/lib61_heatmap_EC50_extant_wildtypes.pdf",sep="")))
```

And for ancestors.

``` r
#make temp long-form data frame
ancestors <- c(config$ancestors_MAP)

temp3 <- temp1[target %in% ancestors,];temp3$target <- factor(temp3$target,levels=ancestors)

p1 <- ggplot(temp3,aes(target,mAb))+geom_tile(aes(fill=log10_EC50),color="black",lwd=0.1)+
  scale_fill_gradientn(colours=c("#003366","#FFFFFF","#FFFFFF"),limits=c(-2,5),values=c(0,7/7,7/7),na.value="gray40")+
  #scale_fill_gradientn(colours=c("#FFFFFF","#003366"),limits=c(5,12),values=c(0,1),na.value="yellow")+
  #scale_x_continuous(expand=c(0,0),breaks=c(331,seq(335,430,by=5)))+
  labs(x="RBD homolog",y="")+theme_classic(base_size=9)+
  coord_equal()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.6,face="bold"))

p1
```

    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm):
    ## collapsing to unique 'x' values

    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm):
    ## collapsing to unique 'x' values

    ## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm):
    ## collapsing to unique 'x' values

<img src="collapse_barcodes_lib61_SARSr-wts_files/figure-gfm/heatmap_wildtypes_MAP-ancestors_EC50-1.png" style="display: block; margin: auto;" />

``` r
invisible(dev.print(pdf, paste(config$final_variant_scores_dir,"/lib61_heatmap_EC50_MAP-ancestors.pdf",sep="")))
```

Save output file.

``` r
dt_final %>%
  mutate_if(is.numeric, round, digits=5) %>%
  write.csv(file=config$final_variant_scores_lib61_file, row.names=F,quote=F)
```
