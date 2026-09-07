################################################################################
## 04_trees.R
## Load the VertLife master tree file and sample n_trees_use trees, once.
## Creates the global object `trees` (a list of 100 trees) used by every PGLS in
## 05_models.R. Requires 00_setup.R and 03_load_dataset.R (for coverage checks).
##
## The file `data/AllBirdsHackett1.tre` is 1000 Newick trees (9993 tips each);
## it is large and NOT stored in the repository -- download it from birdtree.org
## (Hackett backbone, "All species"); see README.
################################################################################

if (!file.exists(path_master))
  stop("Tree file not found: ", path_master,
       "\n  Download the Hackett 'All species' trees from birdtree.org (see README).")

## --- load all trees, then sample regularly-spaced ones ----------------------
trees_all <- ape::read.tree(path_master)
if (inherits(trees_all, "phylo")) trees_all <- list(trees_all)  # single-tree file
cat(sprintf("Tree file: %d trees | %d tips each\n",
            length(trees_all), length(trees_all[[1]]$tip.label)))

idx   <- unique(round(seq(1, length(trees_all), length.out = n_trees_use)))
trees <- trees_all[idx]
cat(sprintf("Sampled %d trees.\n", length(trees)))

rm(trees_all); invisible(gc())   # free the full 1000-tree object

## --- tip coverage of each backbone in the trees -----------------------------
tips <- trees[[1]]$tip.label
cat(sprintf("Tip coverage  | avo_base: %d/%d | merged: %d/%d\n",
            sum(avo_base$tip_label %in% tips), nrow(avo_base),
            sum(merged$tip_label   %in% tips), nrow(merged)))

cat("04_trees.R done.\n")
