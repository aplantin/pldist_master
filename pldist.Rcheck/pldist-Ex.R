pkgname <- "pldist"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('pldist')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("LUniFrac")
### * LUniFrac

flush(stderr()); flush(stdout())

### Name: LUniFrac
### Title: LUniFrac
### Aliases: LUniFrac

### ** Examples

data("bal.long.otus")
data("bal.long.meta")
data("sim.tree")
D2.unifrac <- LUniFrac(otu.tab = bal.long.otus, metadata = bal.long.meta, 
    tree = sim.tree, gam = c(0, 0.5, 1), paired = FALSE, check.input = TRUE)
D2.unifrac[, , "d_1"]   # gamma = 1 (quantitative longitudinal transformation)
D2.unifrac[, , "d_UW"]  # unweighted LUniFrac (qualitative/binary longitudinal transf.)




cleanEx()
nameEx("pldist")
### * pldist

flush(stderr()); flush(stdout())

### Name: pldist
### Title: pldist
### Aliases: pldist

### ** Examples

# Gower distance, paired & quantitative transformation 
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, method = "gower")$D

# Gower distance, paired & qualitative/binary transformation 
pldist(paired.otus, paired.meta, paired = TRUE, binary = TRUE, method = "gower")$D

# Gower distance, longitudinal & quantitative transformation 
pldist(bal.long.otus, bal.long.meta, paired = FALSE, binary = FALSE, method = "gower")$D

# Gower distance, longitudinal & qualitative/binary transformation 
pldist(bal.long.otus, bal.long.meta, paired = FALSE, binary = TRUE, method = "gower")$D

# Other distances 
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, method = "bray")$D
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, method = "kulczynski")$D
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, method = "jaccard")$D

# UniFrac additionally requires a phylogenetic tree and gamma values 
# (Gamma controls weight placed on abundant lineages) 
pldist(paired.otus, paired.meta, paired = TRUE, binary = FALSE, 
    method = "unifrac", tree = sim.tree, gam = c(0, 0.5, 1))$D 
    



cleanEx()
nameEx("pltransform")
### * pltransform

flush(stderr()); flush(stdout())

### Name: pltransform
### Title: pltransform
### Aliases: pltransform

### ** Examples

data("paired.otus")
data("paired.meta")
 # paired transformation
res1 <- pltransform(paired.otus, paired.meta, paired = TRUE, check.input = TRUE) 
 # longitudinal transformation 
res2 <- pltransform(paired.otus, paired.meta, paired = FALSE, check.input = TRUE) 
    



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
