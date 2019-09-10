## Plot ancestral states onto Tetragnatha phylogeny
# Author: Jun Ying Lim (with help from Michael Landis)

## Packages ============
rm(list = ls())
library(ape)
library(ggtree)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(treeio)

source("~/Dropbox/Projects/Resources/revbayes/RevGadgets/R/plot_ancestral_states.R") # RevGadgets local pull

# Notes:
# Stochastic character mapping performed on each tree sampled in the posterior (i.e., each stochastic map is an independent posterior sample of phylogeny, biogeographic histopry and the uncertainty in ages of islands)

## Input data ============
main.dir <- "~/Dropbox/spiders/phylo/tetragnathaBiogeogDating/"
output.dir <- file.path(main.dir, "output")
fig.dir <- file.path(main.dir, "figures")
# data.dir <- file.path(main.dir, "data")

# Import specimen metadata
ref_data <- read.csv("/Users/junyinglim/Dropbox/spiders/phylo/ultros/ultros_islands_SK.csv", stringsAsFactors = F)
ref_data$Species <- gsub(ref_data$Species, pattern = "\"", replacement = "")
ref_data$specimen_Island[ref_data$specimen_Island == "continental North America"] <- "N.Amer"
ref_data$Species_Island <- paste0(ref_data$Species, "_", ref_data$specimen_Island)

# import tree

stem <- "distscale2_root6"
tetragnathaTree <- read.beast(file.path(output.dir, paste0(stem, "_ancstate.tre")))

# Rename tips
tetragnathaTree@phylo$tip.label <- ref_data$Species_Island[match(tetragnathaTree@phylo$tip.label, ref_data$specimen_UID)]

# Create labels
state_labels <- read.csv(file.path(output.dir, "state_labels.txt"),colClasses="character")

area_names = c("R","K","O","M","H")

# map presence-absence ranges to area names
range_labels = sapply(state_labels$range[2:nrow(state_labels)],
                      function(x) {
                        present = as.vector(gregexpr(pattern="1", x)[[1]])
                        paste( area_names[present], collapse="")
                      })

st_lbl = list()
for (j in 1:length(range_labels)) {
  st_lbl[[ as.character(j) ]] = range_labels[j]
}

# Generate colours
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
st_colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:length(range_labels)]

# Use "..." to represent remaining low-prob states
st_colors[ length(st_colors)+1 ] = "lightgray"
st_lbl[["..."]] = "..."

pp = plot_ancestral_states(tree = tetragnathaTree,
                           include_start_states = T,
                           shoulder_label_size = 0,
                           summary_statistic = "PieRange",
                           state_labels = st_lbl,
                           state_colors = st_colors,
                           node_label_size = 0,
                           node_size_range = c(0.5,2),
                           node_label_nudge_x = 0.5,
                           tip_node_size = 3,
                           tip_label_size = 1.8,
                           tip_label_offset = 0.1,
                           xlim_visible = c(-10,50),
                           shoulder_label_nudge_x = 0.0,
                           show_posterior_legend = T,
                           pie_diameter = 3,
                           alpha = 1,
                           show_tree_scale = T)


ggsave(plot = pp + theme_tree2() + coord_cartesian(xlim = c(0,8)) + theme(legend.position = "bottom"),
       file=file.path(fig.dir, paste0(stem, "_ancstateplot.pdf")),
       device="pdf",height=8, width=8, useDingbats=F)



## CODE FROM PLOT_ANCESTRAL_RANGES ========================
t = tetragnathaTree
include_start_states = T
shoulder_label_size = 0
summary_statistic = "PieRange"
state_labels = st_lbl
state_colors = st_colors
node_label_size = 0
node_size_range = c(0.5,2)
node_label_nudge_x = 0.5
xlim_visible = c(-10,50)
ylim_visible = NULL
shoulder_label_nudge_x = 0.0
show_posterior_legend = T
node_pie_diameter = 4
tip_pie_diameter = 4
alpha = 1
fig_width = 10
fig_height = 5
show_tree_scale = T
show_state_legend = T
use_state_colors = T
tree_layout = "rectangular"


# Tip para
tip_node_size = 2
tip_label_size = 1.8
tip_label_offset = 0.1
tip_label_italics = F

# Pie para
pie_diameter=0.1
pie_nudge_x=0.0
pie_nudge_y=0.0
  
t = assign_state_labels(t, state_labels, include_start_states)

# add range for pp factors
t = set_pp_factor_range(t, include_start_states)

# add state colors
use_state_colors = !is.null(state_colors)
if (!is.null(state_colors) && !is.null(state_labels))
{
  names(state_colors) = state_labels
}

tree = attributes(t)$phylo
n_node = ggtree:::getNodeNum(tree)

# remove underscores from tip labels
attributes(t)$phylo$tip.label = gsub("_", " ", attributes(t)$phylo$tip.label)

if (tip_label_italics) {
  attributes(t)$phylo$tip.label = paste("italic('", attributes(t)$phylo$tip.label, "')", sep="")
}

# add tip labels
p = ggtree(t, layout=tree_layout, ladderize=TRUE)
p = p + geom_tiplab(size=tip_label_size, offset=tip_label_offset, parse=tip_label_italics)

# Start of the PieRange specific code
p = p + geom_tippoint(aes(colour=factor(end_state_1)), size=tip_node_size, alpha=alpha) 

# plot invisible node states (for legend)
p = p + geom_nodepoint(aes(colour=factor(start_state_1), size=0),na.rm=TRUE, alpha=0.0)
p = p + geom_nodepoint(aes(colour=factor(start_state_2), size=0),na.rm=TRUE, alpha=0.0)
p = p + geom_nodepoint(aes(colour=factor(start_state_3), size=0),na.rm=TRUE, alpha=0.0)

# set up the legend
if (show_state_legend) {
  p = p + guides(colour=guide_legend("State"), order=1)        
} else {
  p = p + guides(colour=FALSE, order=2)
}
p = p + guides(size=FALSE)
if (use_state_colors) {
  p = p + scale_color_manual(values=state_colors, breaks=state_labels)
}

#p = p + scale_radius(range = node_size_range)
p = p + theme(legend.position="left")

# get anc state matrices (for pie/bar charts)
dat_state_end = build_state_probs(t, state_labels, include_start_states)$end
dat_state_start = build_state_probs(t, state_labels, include_start_states)$start

# make pie objects
n_tips = length(tree$tip.label)
n_nodes = 2 * n_tips - 1
node_idx = (n_tips+1):n_nodes
node_idx_no_root = (n_tips+1):(n_nodes-1)
pies_end = nodepie(dat_state_end,cols=1:(ncol(dat_state_end)-1),color=state_colors,alpha=alpha)
pies_start = nodepie(dat_state_start,cols=1:(ncol(dat_state_start)-1),color=state_colors,alpha=alpha)

# print pies
hjust=-tree$edge.length/2

pie_diameter = 5

p_node = inset.revgadgets(tree_view=p,
                        insets=pies_end[node_idx],
                        x="node",
                        height=pie_diameter,
                        width=pie_diameter,
                        hjust=0,
                        vjust=0)


p_branch = inset.revgadgets(tree_view=p_node,
                            insets=pies_start,
                            x="parent_shoulder",
                            height=pie_diameter*0.85,
                            width=pie_diameter*0.85,
                            hjust=pie_nudge_x,
                            vjust=pie_nudge_y)

p <- p_branch
if (use_state_colors) {
  #print(state_colors)
  #print(state_labels)
  p = p + scale_color_manual(values=state_colors, breaks=as.vector(state_labels))
}

p = p + scale_radius(range = node_size_range)
p = p + theme(legend.position="left")

if (show_tree_scale)
{
  #p = p + theme_tree2()
}
p = p + ggtitle(title)
# set visible area
p = p + coord_cartesian(xlim = xlim_visible, ylim=ylim_visible, expand=TRUE)

#    if (!(summary_statistic %in% c("PieState", "PieRange"))) {
#        print(p)
#    }

ggsave("~/Desktop/pie.pdf",p)
