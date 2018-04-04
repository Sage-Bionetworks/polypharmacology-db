library(cluster)
library(synapser)
library(tidyverse)
library(pbmcapply)
synLogin()
sim.mat <- read.table(synGet("syn11948129")$path, header = T)

db <- readRDS(synGet("syn11712148")$path) %>% select(internal_id, common_name) %>% distinct()

names <- readRDS(synGet("syn11712154")$path) %>% 
  filter(database %in% c("drugbank")) %>% 
  select(internal_id) %>% 
  distinct() %>% 
  inner_join(db)

colnames(sim.mat) <- rownames(sim.mat)

sim.mat.reduced <- sim.mat[rownames(sim.mat) %in% names$internal_id,
                           colnames(sim.mat) %in% names$internal_id]

sim.mat.tidy <- sim.mat.reduced %>% rownames_to_column("a") %>% gather("b", "sim", -a)

dissim.mat <- 1-sim.mat.reduced



# agnes <- agnes(dissim.mat, diss = TRUE)
# names.sort <- names[order(names$internal_id, agnes$order.lab),]
# agnes$order.lab <- names$common_name
# 
# agdg <- as.dendrogram(agnes)
# plot(cut(agdg, h=0.9)$upper, 
#      main="Upper tree of cut at h=0.9")
# plot(cut(agdg, h=0.9)$lower[[74]], 
#      main="Upper tree of cut at h=0.9")
# 
# diana <- diana(dissim.mat, diss = TRUE)
# plot(diana)


hierarchical <- hclust(as.dist(dissim.mat), method = "ward.D2")
names.sort <- names[order(names$internal_id, hierarchical$labels),]
hierarchical$labels <- sapply(names.sort$common_name, function(x) strtrim(x, 30))
hcd <- as.dendrogram(hierarchical)
plot(hcd)
plot(cut(hcd, h=3.5)$upper)
plot(cut(hcd, h=3.5)$lower[[76]], 
     main="Second branch of lower tree with cut at h=75", xlab = "compound", asp = 2)

##Looking at histogram of similarities, the vast majority of similarities are very low
##however, many values still above 0.2. perhaps try this as a cutoff (everything >0.8 in dissim mat)
# hist(sim.mat.tidy$sim)
# hist(sim.mat.tidy$sim[sim.mat.tidy$sim>0.2])
# 
# dist.mat.cutoff <- 1-sim.mat.reduced
# dist.mat.cutoff[sim.mat.cutoff>0.8] <- 1
# 
# hierarchical <- hclust(as.dist(dist.mat.cutoff), method = "ward.D")
# names.sort <- names[order(names$internal_id, hierarchical$labels),]
# hierarchical$labels <- sapply(names.sort$common_name, function(x) strtrim(x, 30))
# hcd <- as.dendrogram(hierarchical)
# plot(hcd)
# plot(cut(hcd, h=2.5)$upper)
# plot(cut(hcd, h=2.5)$lower[[5]], 
#      main="Second branch of lower tree with cut at h=75", xlab = "compound")


db <- readRDS(synGet("syn11712148")$path) %>% distinct() %>% 
  filter(internal_id %in% names$internal_id) %>% 
  filter(mean_pchembl > 7 | n_qualitative > 0)
  
library(pbapply)
library(pbmcapply)
library(parallel)
library(plyr)
targs <- pblapply(unique(db$internal_id), function(x){
  foo <- filter(db, internal_id == x) %>% select(hugo_gene)
  foo$hugo_gene
})

names(targs) <- unique(db$internal_id)

simDrugGeneLists <- function(list.a, list.b){ 
  sims <- pbmclapply(list.a, function(i) {
    sim <- mclapply(list.b, function(j) {

      jaq <- length(intersect(i,j))/length(union(i,j))
      jaq
    }, mc.cores = detectCores())
    bar <- ldply(sim) 
    colnames(bar) <- c("match", "similarity")
    bar
  })
  sims <- ldply(sims) %>%  spread(.id, similarity) %>% remove_rownames() %>%  column_to_rownames("match")
  sims
}

target.jaccard <- simDrugGeneLists(targs,targs)

name.subset <- filter(names.sort, internal_id %in% colnames(target.jaccard))

target.jaccard.dist <- 1-target.jaccard
clust <- hclust(as.dist(target.jaccard.dist), method = "ward.D2")
names.sort.subset <- name.subset[order(name.subset$internal_id, clust$labels),]
clust$labels <- sapply(names.sort.subset$common_name, function(x) strtrim(x, 30))
plot(clust)
hcd.targ <- as.dendrogram(clust)
plot(cut(hcd.targ, h=0.5)$upper)
plot(cut(hcd.targ, h=2)$lower[[10]], 
     main="Second branch of lower tree with cut at h=75", xlab = "compound", asp = 2)


dissim.mat.min <- dissim.mat[rownames(dissim.mat) %in% rownames(target.jaccard.dist),
                                    colnames(dissim.mat) %in% colnames(target.jaccard.dist)]

hierarchical <- hclust(as.dist(dissim.mat.min), method = "ward.D2")

names.min <- names %>% filter(internal_id %in% rownames(dissim.mat.min))
names.sort <- names.min[order(names.min$internal_id, hierarchical$labels),]
hierarchical$labels <- sapply(names.sort$common_name, function(x) strtrim(x, 30))
hcd <- as.dendrogram(hierarchical)

library(dendextend)

cophenetic <- cor_cophenetic(hcd, hcd.targ)
fm <- Bk_plot(hcd, hcd.targ)
