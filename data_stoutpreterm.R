## Stout (preterm birth)
# load OTUS 
otus <- read.csv("~/Documents/Research/Data/StoutPreterm/otu_table.csv")
rownames(otus) <- otus[,1]
otus$OTU_ID <- NULL 
otus <- t(otus) 

# change run IDs (SRR) to sample IDs (SRS) 
idkey <- read.csv("~/Documents/Research/Data/StoutPreterm/id_key.csv")
sampid <- c() 
for (i in 1:nrow(otus)) {
  sampid <- c(sampid, as.character(idkey$SampleName)[which(idkey$Run == rownames(otus)[i])])
}
rownames(otus) <- sampid

# match up metadata 
meta <- read.csv("~/Documents/Research/Data/StoutPreterm/deidentified_metadata.csv")
otus <- otus[which(rownames(otus) %in% meta$psaf), ]
meta <- meta[which(meta$psaf %in% rownames(otus)), ]

meta <- meta[order(meta$psaf), ]
otus <- otus[order(rownames(otus)), ]
all(meta$psaf == rownames(otus))

# prep metadata and otus (two samples per subject)
coll.id <- aggregate(meta$collection_number, by = list(meta$newid), FUN = function(x) x)
with.two <- c() 
for (i in 1:length(unique(meta$newid))) {
  this.id <- unique(meta$newid)[i] 
  if (all(c(1,2) %in% coll.id$x[[i]])) { with.two <- c(with.two, this.id) }
}
with.two  # 52 subj with times 1 and 2 

keep.meta <- meta[which(meta$newid %in% with.two & meta$collection_number %in% c(1,2)), ]
keep.meta$time <- NA 
for (i in 1:52) {
  this.idx <- which(keep.meta$newid == with.two[i])
  keep.meta$time[this.idx[which.min(keep.meta$ga_collection[this.idx])]] <- 1 
  keep.meta$time[this.idx[which.max(keep.meta$ga_collection[this.idx])]] <- 2 
}

pl.meta <- data.frame(subjID = keep.meta$newid, 
                      sampID = keep.meta$psaf, 
                      time = keep.meta$time)
otus <- otus[which(rownames(otus) %in% pl.meta$sampID), ]

# read tree 
tree <- read.tree("~/Documents/Research/Data/StoutPreterm/rooted_tree.nwk")

# outcomes, distances 
y.stout <- unique(data.frame(newid = keep.meta$newid, ga = keep.meta$ga_delivery))
y.stout$preterm <- as.numeric(y.stout$ga <= 37)
Ds <- pldist_all(otus, pl.meta, paired = TRUE, tree = tree)
Ks <- lapply(Ds, FUN = function(d) D2K(d))
all(y.stout$newid == rownames(Ds[[1]]))

MiRKAT(y = y.stout$preterm, X = NULL, Ks = Ks, 
       out_type = "D", nperm = 9999)
