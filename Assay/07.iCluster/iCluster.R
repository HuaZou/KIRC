#!/usr/bin/R
library(argparser)
library(dplyr)
library(tibble)
library(ggplot2)
library(data.table)
library(iCluster)

rm(list = ls())
options(stringsAsFactors = F)
options(future.globals.maxSize = 1000 * 1024^2)

grp <- c("S1", "S2")
grp.col <- c("#6C326C", "#77A2D1")

# parameter input
parser <- arg_parser("iCluster function") %>%
    add_argument("-p", "--phen", 
        help = "phenotype") %>%
    add_argument("-c", "--copyNumber", 
        help = "copyNumber matrix") %>%
    add_argument("-g", "--geneExp", 
        help = "gene Expression matrix") %>%
    add_argument("-m", "--methylation", 
        help = "DNA methylation matrix")  %>%
    add_argument("-r", "--RPPA", 
        help = "protein_RPPA")  %>%
    add_argument("-o", "--out", 
        help = "result with director", default = "./")
    
args <- parse_args(parser)

# prepare for function 
phen <- fread(args$p)                                  
copyNumber <- fread(args$c)
geneExp <- fread(args$g)
methylation <- fread(args$m)
protein_RPPA <- fread(args$r)
out <- args$o	

get_profile <- function(dataset = copyNumber,
                        metadata = phen, 
                        tag = "copyNumber"){
  
  # dataset = protein_RPPA
  # metadata = phen
  # tag = "protein_RPPA"
  
  sid <- intersect(phen$Barcode, colnames(dataset)) 
  res <- dataset %>% dplyr::select(c("V1", sid)) %>% 
    mutate(Type=tag) %>%
    mutate(Name=paste(V1, Type, sep = "_")) %>%
    dplyr::select(Name, V1, Type, everything()) %>%
    dplyr::select(-c("V1", "Type")) %>%
    column_to_rownames("Name") %>%
    t()
  return(res)
}

datasets <- list(copyNumber=get_profile(dataset = copyNumber, tag = "copyNumber"),
                 geneExp=get_profile(dataset = geneExp, tag = "geneExp"),
                 methylation=get_profile(dataset = methylation, tag = "methylation"),
                 protein_RPPA=get_profile(dataset = protein_RPPA, tag = "protein_RPPA"))

print(names(datasets))


# icluster 
fit <- iCluster(datasets = datasets, k=2, lambda=rep(0.2, 4), max.iter = 50, epsilon = 1e-3)
#plotiCluster(fit=fit, label=rownames(datasets[[2]]))
#compute.pod(fit)

phen_new_icluster <- inner_join(phen,
              data.frame(SampleID=rownames(datasets[[2]]), Cluster=fit$clusters), 
             by = c("Barcode"="SampleID")) %>% 
  dplyr::select(Barcode, Cluster, everything()) %>%
  mutate(Cluster=paste0("S", Cluster))

name1 <- paste(out, "phenotype_cluster_iCluster.csv")
write.csv(phen_new_icluster, file = name1, row.names = F)


fit2 <- iCluster2(datasets = datasets, k=2, lambda=list(0.2, 0.2, 0.2, 0.2),
                 max.iter = 50, verbose = TRUE)

if(0){
    phen_new_icluster2 <- inner_join(phen,
                data.frame(SampleID=rownames(datasets[[2]]), Cluster=fit2$clusters), 
                by = c("Barcode"="SampleID")) %>% 
    dplyr::select(Barcode, Cluster, everything()) %>%
    mutate(Cluster=paste0("S", Cluster))

    name2 <- paste(out, "phenotype_cluster_iCluster2.csv")
    write.csv(phen_new_icluster2, file = name2, row.names = F)
}

save(fit, fit2, file = "iCluster_fit.RData")
