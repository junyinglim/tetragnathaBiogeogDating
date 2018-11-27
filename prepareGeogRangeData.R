library(ape)

main.dir <- "/Users/junyinglim/Dropbox/spiders/phylo/ultros/"
output.dir <- "/Users/junyinglim/Dropbox/spiders/phylo/tetragnathaBiogeogDating/data"

biogeog<- read.csv(file.path(main.dir, "ultros_islands_SK.csv"))
biogeog_subset <- biogeog[c("Co","Ka","Oa","Ma", "Ha")]
rownames(biogeog_subset) <- biogeog$specimen_UID

write.nexus.data(as.matrix(biogeog_subset), file.path(output.dir, "spiders_biogeog.nex"), format = "Standard")
