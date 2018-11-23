
library(ape)

main.dir <- "~/Dropbox/spiders/phylo/jointBiogeogDating/"
output.dir <- file.path(main.dir, "output")
data.dir <- file.path(main.dir, "data")
state_labels <- read.csv(file.path(output.dir, "state_labels.txt"))
ref_data <- read.csv(file.path(data.dir, "vargas_islands_jyl.csv"))
state_labels$rangeLab <- c("Null", "C", "K", "O", "M", "H",
                           "C+K", "C+O", "C+M", "C+H",
                           "K+M", "O+M", "C+H", "K+H", "O+H", "M+H",
                           "C+K+O", "C+K+M", "C+O+M", "K+O+M", "C+K+H",
                           "C+O+H", "K+O+H", "C+M+H", "K+M+H", "O+M+H",
                           "C+K+O+M","C+K+O+H", "C+K+M+H", "C+O+M+H",
                           "K+O+M+H", "C+K+O+M+H")

trees <- read.tree(file.path(output.dir, "ucln_noconstraints_run_1.trees"))
anc_states <- read.delim(file.path(output.dir, "states_run_1.log"))

library(ggtree)
library(reshape2)
y <- melt(anc_states[100,], id.vars = "Iteration")
z <- y[grep(y$variable, pattern = "end_"),]
z$node <- gsub(z$variable, pattern = "end_", replacement = "")
rownames(z) <- z$node
z$value2 <- state_labels$rangeLab[match(z$value, state_labels$state)]


trees_prelim <- trees[[300]]
trees_prelim$tip.label <- as.vector(ref_data$Species[match(trees_prelim$tip.label, ref_data$specimen_UID)])

asd <- ggtree(tr = trees_prelim) %<+% z
asd2 <- asd + geom_nodepoint(aes(colour = factor(value2))) + geom_tippoint(aes(colour = factor(value2))) + theme(legend.position = "bottom") + geom_tiplab(size = 1)

ggsave(asd2, filename = file.path(main.dir, "anc_state.pdf"), device = cairo_pdf, height = 10, width = 10)
