## Romero (preterm birth - no differences detected in original paper)
library(ape)
#devtools::install_github("aplantin/pldist")
library(pldist)
library(MiRKAT)

##### All time points ###### 

meta <- read.csv("~/Documents/Research/Data/RomeroPreterm/metadata_processed.csv")
otus <- read.csv("~/Documents/Research/Data/RomeroPreterm/otu_tab.csv")
rownames(otus) <- otus[,1] 
otus$sample.ID <- NULL

meta <- meta[-which(meta$SubjectID %in% names(which(table(meta$SubjectID) == 1))), ]
otus <- otus[which(rownames(otus) %in% meta$SampleID), ]

pl.meta <- data.frame(subjID = meta$SubjectID, 
                      sampID = as.character(meta$SampleID), 
                      time = meta$GesAge)
pl.meta <- pl.meta[order(pl.meta$sampID), ]
otus <- otus[order(rownames(otus)), ]
Ds <- pldist_all(otus, pl.meta, paired = FALSE, method = c("b", "j", "g", "k"))
#Ds <- pldist_all(otus, pl.meta, paired = FALSE, method = c("b", "j"))
Ks <- lapply(Ds, FUN = function(d) D2K(d))

y.romero <- data.frame(subjID = meta$SubjectID, 
                       ptb = meta$PregnancyOutcome, 
                       bw = meta$BirthweightGms, 
                       race = meta$Race, 
                       bmi = meta$PrePregnancyBMI, 
                       age = meta$Age)
y.romero <- unique(y.romero)
y.romero$ptb <- as.numeric(y.romero$ptb == "PTD")
all(rownames(Ds[[1]]) == y.romero$subjID)

res1 <- MiRKAT(y = y.romero$ptb, X = cbind(y.romero$race, y.romero$bmi, y.romero$age), Ks = Ks, out_type = "D"); res1
res2 <- MiRKAT(y = y.romero$bw, X = cbind(y.romero$race, y.romero$bmi, y.romero$age), Ks = Ks, out_type = "C"); res2



#####  Two time points (13-27 and closest to 40)  #####
meta <- read.csv("~/Documents/Research/Data/RomeroPreterm/metadata_processed.csv")
otus <- read.csv("~/Documents/Research/Data/RomeroPreterm/otu_tab.csv")
rownames(otus) <- otus[,1] 
otus$sample.ID <- NULL
otus <- otus[order(rownames(otus)), ]
meta <- meta[order(meta$SampleID), ]
all(meta$SampleID == rownames(otus))


ga.sid <- aggregate(meta$GesAge, by = list(meta$SubjectID), FUN = function(x) x)
keep.idx <- c() 
time.bin <- c() 
for (i in 1:length(unique(meta$SubjectID))) {
  this.id <- ga.sid$Group.1[i]
  this.ga <- ga.sid$x[[i]]
  if (length(this.ga) > 1) {
    #tri2 <- this.ga[this.ga < 28 & this.ga > 12]
    #tri3 <- this.ga[this.ga >= 28] 
    #if (length(tri2) > 0 & length(tri3) > 0) {
      #tri2.idx <- which(meta$SubjectID == this.id & meta$GesAge == min(tri2)) 
      tri2.idx <- which(meta$SubjectID == this.id & meta$GesAge == min(this.ga)) 
      #tri3.idx <- which(meta$SubjectID == this.id & meta$GesAge == max(tri3))
      tri3.idx <- which(meta$SubjectID == this.id & meta$GesAge == max(this.ga))
      keep.idx <- c(keep.idx, tri2.idx, tri3.idx)
      time.bin <- c(time.bin, 0, 1)
    #}
  }
}

meta2 <- meta[keep.idx, ]
otus2 <- otus[keep.idx, ]
pl.meta2 <- data.frame(subjID = meta2$SubjectID, 
                       sampID = meta2$SampleID,
                       time = time.bin[order(keep.idx)])



#Ds <- pldist_all(otus2, pl.meta2, paired = TRUE, method = c("b", "j", "g", "k"))
Ds <- pldist_all(otus2, pl.meta2, paired = TRUE, method = c("j"))
Ks <- lapply(Ds, FUN = function(d) D2K(d))

y.romero2 <- data.frame(subjID = meta2$SubjectID, 
                       ptb = meta2$PregnancyOutcome, 
                       bw = meta2$BirthweightGms, 
                       race = meta2$Race, 
                       bmi = meta2$PrePregnancyBMI, 
                       age = meta2$Age)
y.romero2 <- unique(y.romero2)
y.romero2$ptb <- as.numeric(y.romero2$ptb == "PTD")
all(rownames(Ds[[1]]) == y.romero2$subjID)

res1 <- MiRKAT(y = y.romero2$ptb, X = cbind(y.romero2$race, y.romero2$bmi, y.romero2$age), Ks = Ks, out_type = "D"); res1
res2 <- MiRKAT(y = y.romero2$bw, X = cbind(y.romero2$race, y.romero2$bmi, y.romero2$age), Ks = Ks, out_type = "C"); res2



##### Single time point ##### 