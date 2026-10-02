################################################################################
## 06_figures.R
## Paper-ready figures, written to output/figures/ (PDF + PNG preview), each
## with a ready-to-paste LaTeX figure environment (*.tex, caption included).
##
##   Figure 1  : Block A main effects of 3D                        [main text]
##   Figure S1 : Block A within Passeriformes                      [appendix]
##   Figure S2 : Block B2, dim_bin x continuous parental care      [appendix]
##
## Needs only the model cache (output/model_results.rds, from 05_models.R) and
## 00_setup.R: no refitting, so figures can be restyled in seconds.
##
## Forest-plot conventions (all figures):
##   point      = Rubin-pooled estimate across the 100 trees
##   thick bar  = 90% CI  (excludes 0 <=> one-sided p < 0.05)
##   thin bar   = 95% CI  (excludes 0 <=> two-sided p < 0.05)
##   colour     = sign relative to the a priori prediction
##   right text = N species | % of trees with the sign of the pooled estimate
##   continuous responses in SD units of the response (= PGLS on z-scored
##   trait); binary responses in log-odds.
################################################################################

if (!exists("model_results")) model_results <- readRDS(path_model_cache)
mr <- model_results

## --- palette: Okabe-Ito (CVD-safe), doubled by shape --------------------------
agree_levels <- c("predicted", "opposite", "none")
agree_labels <- c(predicted = "as predicted",
                  opposite  = "opposite to prediction",
                  none      = "no directional prediction")
agree_cols   <- c(predicted = "#0072B2", opposite = "#D55E00", none = "#6E6E6E")
agree_shapes <- c(predicted = 16, opposite = 18, none = 15)

## --- summary rows for a set of results -----------------------------------------
## spec: rows of specA / specB2; res: named list of run_pgls() results.
## Response labels carry their source when needed to tell rows apart.
fig_rows <- function(spec, res, std = TRUE) {
  dup <- spec$label[duplicated(spec$label)]
  out <- lapply(seq_len(nrow(spec)), function(i) {
    r <- summarise_res(res[[spec$id[i]]], std = std)
    if (is.null(r)) return(NULL)
    short <- source_info$short[source_info$source == spec$source[i]]
    r$row_label <- if (spec$label[i] %in% dup || spec$block[i] == "ss")
      paste0(spec$label[i], " (", short, ")") else spec$label[i]
    r$block <- spec$block[i]
    r
  })
  df <- bind_rows(out)
  df$row_label <- gsub("$-$", "-", df$row_label, fixed = TRUE)
  df
}

## --- one forest panel ----------------------------------------------------------
## Returns two aligned plots: the forest itself, and a text strip on its right
## with N and % trees (a separate plot, so its position does not depend on the
## x-range or on the width of the response labels).
## xlim: common x-range (so panels sharing a scale are directly comparable).
forest_panel <- function(df, title, xlab = NULL, xlim, header = TRUE) {
  df$row_label <- factor(df$row_label, levels = rev(unique(df$row_label)))
  df$agrees    <- factor(df$agrees, levels = agree_levels)
  ny <- nlevels(df$row_label)
  forest <- ggplot(df, aes(y = row_label, colour = agrees, shape = agrees)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
    geom_linerange(aes(xmin = lo95, xmax = hi95), linewidth = 0.45,
                   show.legend = FALSE) +
    geom_linerange(aes(xmin = lo90, xmax = hi90), linewidth = 1.5,
                   show.legend = FALSE) +
    geom_point(aes(x = est), size = 2.3, show.legend = TRUE) +
    scale_colour_manual(values = agree_cols, labels = agree_labels,
                        drop = FALSE, name = NULL) +
    scale_shape_manual(values = agree_shapes, labels = agree_labels,
                       drop = FALSE, name = NULL) +
    scale_x_continuous(breaks = scales::breaks_pretty(5)) +
    coord_cartesian(xlim = xlim, ylim = c(0.5, ny + 0.5), expand = FALSE) +
    labs(title = title, x = xlab, y = NULL) +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 9, hjust = 0,
                                    margin = margin(b = 6)),
          plot.title.position = "plot",
          axis.text.y  = element_text(colour = "grey10"),
          axis.ticks.y = element_blank(),
          axis.line.y  = element_blank(),
          panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3),
          plot.margin = margin(4, 2, 4, 4))
  strip <- ggplot(df, aes(y = row_label)) +
    geom_text(aes(x = 1, label = format(N, big.mark = ",")), size = 2.6,
              hjust = 1, colour = "grey20") +
    geom_text(aes(x = 2.5, label = paste0(pct_same, "%")), size = 2.6,
              hjust = 1, colour = "grey20") +
    coord_cartesian(xlim = c(0, 2.6), ylim = c(0.5, ny + 0.5), expand = FALSE,
                    clip = "off") +
    labs(title = if (header) " " else NULL) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(size = 9, margin = margin(b = 6)),
          plot.margin = margin(4, 4, 4, 0))
  if (header)
    strip <- strip + annotate("text", x = c(1, 2.5), y = ny + 0.5, vjust = -0.6,
                              label = c("N", "% trees"), size = 2.5, hjust = 1,
                              fontface = "bold", colour = "grey20")
  list(forest = forest, strip = strip)
}

## common x-range over a set of panels (95% CIs), padded
x_range <- function(...) {
  d <- bind_rows(...)
  r <- range(c(0, d$lo95, d$hi95), na.rm = TRUE)
  r + c(-1, 1) * 0.06 * diff(r)
}

## stack panels (patchwork), heights proportional to the number of rows;
## one collected legend at the bottom
stack_panels <- function(panels, nrows) {
  cells <- unlist(lapply(panels, function(p) list(p$forest, p$strip)),
                  recursive = FALSE)
  patchwork::wrap_plots(cells, ncol = 2, widths = c(1, 0.24),
                        heights = nrows + 1.4, guides = "collect") &
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8, face = "bold"),
          legend.margin = margin(0, 0, 0, 0),
          legend.key.spacing.x = unit(2, "pt"))
}

save_fig <- function(g, name, width, height) {
  ggsave(file.path(dir_fig, paste0(name, ".pdf")), g, width = width, height = height)
  ggsave(file.path(dir_fig, paste0(name, ".png")), g, width = width, height = height,
         dpi = 200, bg = "white")
  cat("  saved", name, "(.pdf, .png)\n")
}

## write a ready-to-paste LaTeX figure environment
write_fig_tex <- function(name, caption, label, width = "0.95\\textwidth") {
  tex <- c("\\begin{figure}[htbp]", "\\centering",
           sprintf("\\includegraphics[width=%s]{%s.pdf}", width, name),
           sprintf("\\caption{%s}", caption),
           sprintf("\\label{%s}", label), "\\end{figure}")
  writeLines(tex, file.path(dir_fig, paste0(name, ".tex")))
}

bar_legend <- paste0(
  "Points are Rubin-pooled estimates across 100 phylogenies (combining within- and ",
  "between-tree variance); thick bars are 90\\% and thin bars 95\\% confidence ",
  "intervals. Because predictions are directional, a thick bar excluding zero ",
  "corresponds to a one-sided $p < 0.05$ (for intensity, which has no directional ",
  "prediction, only the 95\\% interval is relevant). Colour and symbol indicate whether ",
  "the estimate lies in the predicted direction. Right-hand columns: number of ",
  "species and percentage of trees in which the estimate has the same sign as the ",
  "pooled estimate.")

################################################################################
## FIGURE 1: Block A main effects ----------------------------------------------
################################################################################
cat("\n=== Figure 1 (Block A main effects) ===\n")

## one forest figure for a Block A-like result list (used for Fig 1 and Fig S1)
block_a_figure <- function(res) {
  sp  <- mr$specA %>% filter(in_fig)
  bin <- fig_rows(sp %>% filter(block == "ss", type == "binary"), res)
  int <- fig_rows(sp %>% filter(block == "ss", type == "gaussian"), res)
  con <- fig_rows(sp %>% filter(block == "contest"), res)
  cho <- fig_rows(sp %>% filter(block == "choice", type == "gaussian"), res)
  chb <- fig_rows(sp %>% filter(block == "choice", type == "binary"), res)
  xl_z <- x_range(int, con, cho)
  panels <- list(
    forest_panel(bin, "a   Sexual selection: mating system",
                 "Effect of 3D lifestyle (log-odds)", x_range(bin)),
    forest_panel(int, "b   Sexual selection: intensity", NULL, xl_z, header = FALSE),
    forest_panel(con, "c   Contest: sexual size dimorphism", NULL, xl_z, header = FALSE),
    forest_panel(cho, "d   Choice and display",
                 "Effect of 3D lifestyle (SD units of the response)", xl_z,
                 header = FALSE))
  nrows <- c(nrow(bin), nrow(int), nrow(con), nrow(cho))
  ## binary choice responses (allopreening), on the log-odds scale; the panel
  ## is dropped when none is estimable (e.g. within Passeriformes)
  if (nrow(chb) > 0) {
    panels <- c(panels, list(
      forest_panel(chb, "e   Choice: pair bond", "Effect of 3D lifestyle (log-odds)",
                   x_range(chb), header = FALSE)))
    nrows <- c(nrows, nrow(chb))
  }
  list(plot = stack_panels(panels, nrows), n = sum(nrows), n_panels = length(panels))
}

f1 <- block_a_figure(mr$resA)
save_fig(f1$plot, "fig1_main_effects", width = 6.8,
         height = 0.3 + 0.35 * f1$n_panels + 0.27 * f1$n)
write_fig_tex(
  "fig1_main_effects",
  paste0(
    "\\textbf{Effect of a three-dimensional lifestyle on the form and intensity of ",
    "sexual selection.} Each row is a separate phylogenetic model of the response ",
    "on the 2D/3D contrast (3D = Insessorial + Aerial vs 2D = Terrestrial). ",
    "\\textbf{(a)} Mating-system contrasts, each against monogamy, and ",
    "\\textbf{(e)} allopreening between pair members (penalised phylogenetic ",
    "logistic regression; log-odds). \\textbf{(b--d)} Continuous responses (PGLS, ",
    "Pagel's $\\lambda$), in standard-deviation units of the response so that ",
    "panels b--d share a common scale. ", bar_legend,
    " Bony spurs, which show near-complete separation between 2D and 3D, are ",
    "reported in ",
    "Table~\\ref{tab:descriptive} and Table~\\ref{tab:mainfull}."),
  "fig:effects")

################################################################################
## FIGURE S1: Block A within Passeriformes -------------------------------------
################################################################################
cat("\n=== Figure S1 (Passeriformes) ===\n")

resP <- lapply(mr$resC[grepl("_passer$", names(mr$resC))], `[[`, "ctl")
names(resP) <- sub("_passer$", "", names(resP))
fS1 <- block_a_figure(resP)
save_fig(fS1$plot, "figS1_passeriformes", width = 6.8,
         height = 0.3 + 0.35 * fS1$n_panels + 0.27 * fS1$n)
write_fig_tex(
  "figS1_passeriformes",
  paste0(
    "\\textbf{Main effects of a three-dimensional lifestyle within Passeriformes.} ",
    "Same models and conventions as Fig.~\\ref{fig:effects}, restricted to the one ",
    "order large enough for a within-clade test. Mating-system contrasts that are ",
    "not estimable within Passeriformes ((quasi-)complete separation) are omitted. ",
    bar_legend),
  "fig:passeriformes")

################################################################################
## FIGURE S2: Block B2, dim_bin x continuous parental care ---------------------
################################################################################
cat("\n=== Figure S2 (care interaction) ===\n")

spB <- mr$specB2 %>% filter(type == "gaussian")
b2_con <- fig_rows(spB %>% filter(block == "contest"), mr$resB2)
b2_cho <- fig_rows(spB %>% filter(block == "choice"),  mr$resB2)
xl_b2  <- x_range(b2_con, b2_cho)
fS2 <- stack_panels(list(
  forest_panel(b2_con, "a   Contest: sexual size dimorphism (predicted < 0)", NULL, xl_b2),
  forest_panel(b2_cho, "b   Choice and display (predicted > 0)",
               "Interaction 3D × parental care (SD units of the response per SD of care)",
               xl_b2, header = FALSE)),
  c(nrow(b2_con), nrow(b2_cho)))
save_fig(fS2, "figS2_care_interaction", width = 6.8,
         height = 1.5 + 0.27 * (nrow(b2_con) + nrow(b2_cho)))
write_fig_tex(
  "figS2_care_interaction",
  paste0(
    "\\textbf{Does parental-care asymmetry modulate the effect of a three-dimensional ",
    "lifestyle?} Interaction coefficient 3D $\\times$ care from ",
    "$\\mathrm{trait} \\sim \\mathrm{3D} \\times \\mathrm{care}$ ",
    "(care: relative care investment of the sexes, $z$-scored; positive = ",
    "female-biased). Under the channel-switch hypothesis, female-biased care should ",
    "strengthen the 3D reduction of contest traits (negative interaction, a) and the ",
    "3D increase of choice traits (positive interaction, b). ", bar_legend,
    " Allopreening (binary) is reported in Table~\\ref{tab:care}."),
  "fig:care")

cat("\n06_figures.R done. Figures in", dir_fig, "\n")
