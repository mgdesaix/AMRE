### Isolation by distance

library(geosphere)
library(ade4)
library(ggplot2)

pairwise.fst <- read.table("./data/pairwise.FST.txt")
colnames(pairwise.fst) <- c("Pairs", "UnweightedFST", "WeightedFST")
sites <- read.csv("./data/sites.csv")

### Unweighted FST matrix
# 14 total sites
fst.unweighted.matrix <- matrix(nrow = 14, ncol = 14)
fst.unweighted.matrix[lower.tri(fst.unweighted.matrix)] <- pairwise.fst$UnweightedFST
diag(fst.unweighted.matrix) <- 0
unweighted.gen.dist <- as.dist(fst.unweighted.matrix)

### Weighted FST matrix
# 14 total sites
fst.weighted.matrix <- matrix(nrow = 14, ncol = 14)
fst.weighted.matrix[lower.tri(fst.weighted.matrix)] <- pairwise.fst$WeightedFST
diag(fst.weighted.matrix) <- 0
weighted.gen.dist <- as.dist(fst.weighted.matrix)

### Distance matrix
xy <- sites[,c("Lon", "Lat")]
# Returns matrix of distances in kilometers
geo <- matrix(nrow = nrow(xy), ncol = nrow(xy))
for(i in 1:nrow(xy)) {
  for(j in 1:nrow(xy)) {
    geo[i,j] <- distGeo(xy[i,], xy[j,])
  }
}
# change to kilometers from meters
geo <- geo/1000
rownames(geo) <- sites$Sites
geo.dist <- as.dist(geo)


### Isolation by distance (Unweighted FST)
# log10 distance
ibd.unweighted <- mantel.randtest(unweighted.gen.dist, log10(geo.dist), nrepet = 99999)
ibd.unweighted

### Weighted FST
# 
ibd.weighted <- mantel.randtest(weighted.gen.dist, log10(geo.dist), nrepet = 99999)
ibd.weighted



df <- data.frame("Geodist" = as.vector(geo.dist), 
                 "WeigthedGendist" = as.vector(weighted.gen.dist),
                 "UnweightedGendist" = as.vector(unweighted.gen.dist))
ggplot(df) +
  geom_point(color = "red", alpha = 0.7, aes(x = log10(Geodist), y = WeigthedGendist)) +
  geom_point(color = "blue", alpha = 0.7, aes(x = log10(Geodist), y = UnweightedGendist)) +
  ylab("Genetic Distance") +
  xlab("Geographic Distance (log10)") +
  theme_bw()

