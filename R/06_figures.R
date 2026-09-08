################################################################################
## 06_figures.R
## Figure 1 (form of sexual selection), Figure 2 (resource-sharing interaction)
## and Figure 3 (within-Passeriformes replication). Written to output/.
##
## Panel A (mating system) reuses the log-odds res_* from 05_models.R.
## Panel B (morphology + behaviour) is refitted here on z-scored responses so the
## standardised coefficients are comparable within the panel.
## Requires 00_setup.R, 03_load_dataset.R, 04_trees.R, 05_models.R.
################################################################################

## continuous traits shown, z-scored, in panel B (both Fig 1 and Fig 3)
cont_traits <- c("ssd_mass", "ssd_tarsus", "ssd_wing", "ssd_bill", "ssd_tail",
                 "display_num")

## --- shared figure helpers ---------------------------------------------------
## vivid colour if the across-tree range excludes zero (regardless of sign)
.fig_colour <- function(med, lo, hi) {
  sig <- (lo > 0) | (hi < 0)
  paste0(if_else(med < 0, "neg_", "pos_"), if_else(sig, "sig", "ns"))
}
cols <- c(neg_sig = "#D55E00", pos_sig = "#009E73",
          neg_ns  = "#E8A87C", pos_ns  = "#80CBC4")

## one-row summary of a log-odds res_* object (panel A)
summ_logodds <- function(res, label) {
  if (is.null(res)) return(NULL)
  data.frame(label = label, median = median(res$coefs),
             lo = unname(quantile(res$coefs, 0.025)),
             hi = unname(quantile(res$coefs, 0.975)))
}

## one-row summary of a z-scored PGLS across the 100 trees (panel B)
summ_z <- function(response, label, dat) {
  fits <- lapply(trees, .fit_gaussian, response = response, dat = dat,
                 drop_na_response = TRUE)
  ok   <- Filter(Negate(is.null), fits)
  if (!length(ok)) return(NULL)
  coefs <- sapply(ok, function(f) coef(f)[COEF_3D])
  data.frame(label = label, median = median(coefs),
             lo = unname(quantile(coefs, 0.025)),
             hi = unname(quantile(coefs, 0.975)))
}

## shared forest-plot panel
make_forest <- function(df, subtitle, xlab, caption = NULL) {
  ggplot(df, aes(median, label, xmin = lo, xmax = hi, colour = couleur)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey30") +
    geom_errorbarh(height = 0.22, linewidth = 0.9) +
    geom_point(size = 3.5) +
    scale_colour_manual(values = cols, guide = "none") +
    labs(subtitle = subtitle, x = xlab, y = NULL, caption = caption) +
    theme_bw(base_size = 11) +
    theme(plot.subtitle = element_text(face = "bold"),
          plot.caption = element_text(colour = "grey50", size = 8),
          panel.grid.minor = element_blank())
}

## stack two panels with a title (patchwork -> cowplot -> gridExtra fallback)
save_two_panels <- function(top, bottom, path, title) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    g <- patchwork::wrap_plots(top, bottom, ncol = 1, heights = c(1.5, 3.2)) +
      patchwork::plot_annotation(
        title = title,
        theme = ggplot2::theme(plot.title =
                                 ggplot2::element_text(face = "bold", size = 13)))
    ggsave(path, g, width = 7, height = 6.8)
  } else if (requireNamespace("cowplot", quietly = TRUE)) {
    body <- cowplot::plot_grid(top, bottom, ncol = 1,
                               rel_heights = c(1.5, 3.2), align = "v")
    ttl  <- cowplot::ggdraw() +
      cowplot::draw_label(title, fontface = "bold", size = 13, x = 0, hjust = 0)
    ggsave(path, cowplot::plot_grid(ttl, body, ncol = 1, rel_heights = c(0.08, 1)),
           width = 7, height = 6.8)
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    g <- gridExtra::arrangeGrob(
      top, bottom, ncol = 1, heights = c(1.5, 3.2),
      top = grid::textGrob(title, gp = grid::gpar(fontface = "bold", fontsize = 13)))
    ggsave(path, g, width = 7, height = 6.8)
  } else {
    warning("No layout package (patchwork/cowplot/gridExtra): saving panels separately.")
    ggsave(sub("\\.pdf$", "_A.pdf", path), top,    width = 7, height = 2.2)
    ggsave(sub("\\.pdf$", "_B.pdf", path), bottom, width = 7, height = 4.6)
  }
}

################################################################################
## FIGURE 1 -------------------------------------------------------------------
################################################################################
cat("\n=== Figure 1 ===\n")

## panel A: mating system (log-odds)
figA_df <- bind_rows(
  summ_logodds(res_barber_poly,      "Strong polygamy vs monogamy (Barber)"),
  summ_logodds(res_marc_rdp_vs_mono, "Resource-defense polygamy vs monogamy (Marcondes)"),
  summ_logodds(res_barber_lek,       "Lek vs monogamy (Barber)"),
  summ_logodds(res_marc_lek_vs_mono, "Lekking vs monogamy (Marcondes)")
) %>%
  mutate(label = factor(label, levels = rev(c(
           "Strong polygamy vs monogamy (Barber)",
           "Resource-defense polygamy vs monogamy (Marcondes)",
           "Lek vs monogamy (Barber)",
           "Lekking vs monogamy (Marcondes)"))),
         couleur = .fig_colour(median, lo, hi))
figA <- make_forest(figA_df,
                    "A. Mating system (log-odds, 3D vs 2D; monogamy-referenced)",
                    "Log-odds (3D vs 2D)")

## panel B: morphology + behaviour (z-scored)
dat_z <- merged %>%
  filter(!is.na(dim_bin), !is.na(tip_label)) %>%
  distinct(tip_label, .keep_all = TRUE) %>%
  mutate(across(all_of(cont_traits), ~ as.numeric(scale(.)))) %>%
  as.data.frame()
rownames(dat_z) <- dat_z$tip_label

dat_z_avo <- avo_base %>%
  filter(!is.na(dim_bin), !is.na(dichromatism), !is.na(tip_label)) %>%
  distinct(tip_label, .keep_all = TRUE) %>%
  mutate(across(c(dichromatism, male_plumage, female_plumage),
                ~ as.numeric(scale(.)))) %>%
  as.data.frame()
rownames(dat_z_avo) <- dat_z_avo$tip_label

## intensity (Barber) z-scored, for the panel-B row (an OUTCOME, shown for
## comparison alongside the form traits -- provisional placement, revisit w/ Jan)
dat_z_barber <- avo_base %>%
  filter(!is.na(dim_bin), !is.na(barber_ss), !is.na(tip_label)) %>%
  distinct(tip_label, .keep_all = TRUE) %>%
  mutate(barber_ss = as.numeric(scale(barber_ss))) %>%
  as.data.frame()
rownames(dat_z_barber) <- dat_z_barber$tip_label

figB_df <- bind_rows(
  summ_z("ssd_mass",       "Mass dimorphism",           dat_z),
  summ_z("ssd_tarsus",     "Tarsus dimorphism",         dat_z),
  summ_z("ssd_wing",       "Wing dimorphism",           dat_z),
  summ_z("ssd_bill",       "Bill dimorphism",           dat_z),
  summ_z("ssd_tail",       "Tail dimorphism",           dat_z),
  summ_z("display_num",    "Display agility",            dat_z),
  summ_z("dichromatism",   "Plumage dichromatism (M-F)", dat_z_avo),
  summ_z("male_plumage",   "Male plumage elaboration",   dat_z_avo),
  summ_z("female_plumage", "Female plumage elaboration", dat_z_avo),
  summ_z("barber_ss",      "Sexual-selection intensity (Barber)", dat_z_barber)
) %>%
  mutate(label = factor(label, levels = rev(c(
           "Mass dimorphism", "Tarsus dimorphism", "Tail dimorphism",
           "Wing dimorphism", "Bill dimorphism",
           "Display agility",
           "Plumage dichromatism (M-F)",
           "Male plumage elaboration", "Female plumage elaboration",
           "Sexual-selection intensity (Barber)"))),
         couleur = .fig_colour(median, lo, hi))
figB <- make_forest(
  figB_df, "B. Morphology and behaviour (PGLS, z-scored responses)",
  "Standardised coefficient (3D vs 2D)",
  caption = paste0(
    "Negative = trait reduced in 3D; positive = increased in 3D. ",
    "Vivid hue = across-tree range excludes 0.\n",
    "Panels use different x-scales (log-odds vs standardised) and are not ",
    "directly comparable in magnitude."))

save_two_panels(figA, figB, file.path(dir_out, "fig1_form_of_selection.pdf"),
                "Effect of 3D lifestyle on the form of sexual selection")
cat("  saved fig1_form_of_selection.pdf\n")

################################################################################
## FIGURE 2: resource-sharing interaction (marginal predictions, tree 1) ------
################################################################################
cat("\n=== Figure 2 ===\n")

## NB: this plot uses SCALED resource_c for both traits (so the two panels are
## drawn on a common resource axis). The manuscript's baked-in beta values were
## on inconsistent scales (see 05_models.R note); the caption here is qualitative
## to avoid re-embedding those numbers -- to be finalised in the revision.
marginal_grid <- function(response, trait_label) {
  d <- merged %>%
    filter(!is.na(dim_bin), !is.na(.data[[response]]),
           !is.na(resource), !is.na(tip_label)) %>%
    mutate(resource_c = as.numeric(scale(resource))) %>%
    distinct(tip_label, .keep_all = TRUE) %>% as.data.frame()
  rownames(d) <- d$tip_label
  tr <- ape::drop.tip(trees[[1]], setdiff(trees[[1]]$tip.label, rownames(d)))
  d1 <- d[tr$tip.label, ]
  fit <- phylolm::phylolm(as.formula(paste(response, "~ dim_bin * resource_c")),
                          data = d1, phy = tr, model = "lambda")
  co  <- coef(fit)
  res_levels <- quantile(d1$resource_c, c(0.15, 0.5, 0.85))
  g <- expand.grid(dim_bin = factor(DIM_LEVELS, levels = DIM_LEVELS),
                   resource_c = res_levels)
  g$dim3d <- as.integer(g$dim_bin == "3D")
  g$pred  <- co["(Intercept)"] + co[COEF_3D] * g$dim3d +
             co["resource_c"] * g$resource_c +
             co[COEF_INT_C] * g$dim3d * g$resource_c
  g$res_lab <- factor(rep(c("Low resource-sharing\n(no shared territory)",
                            "Medium",
                            "High resource-sharing\n(year-round territory)"),
                          each = 2),
                      levels = c("Low resource-sharing\n(no shared territory)",
                                 "Medium",
                                 "High resource-sharing\n(year-round territory)"))
  g$trait <- trait_label
  g
}

fig2_df <- bind_rows(
  marginal_grid("dichromatism", "Plumage dichromatism (male - female)"),
  marginal_grid("display_num",  "Display agility (1-5)")
)
fig2_df$trait <- factor(fig2_df$trait,
                        levels = c("Plumage dichromatism (male - female)",
                                   "Display agility (1-5)"))

fig2 <- ggplot(fig2_df, aes(dim_bin, pred, group = res_lab, colour = res_lab)) +
  geom_line(linewidth = 1.1) + geom_point(size = 3) +
  facet_wrap(~ trait, scales = "free_y") +
  scale_colour_manual(values = c("#D55E00", "grey55", "#0072B2"), name = NULL) +
  labs(
    title = "Both choice-related traits depend on between-sex resource sharing",
    subtitle = "Predicted trait value from the PGLS interaction model (tree 1)",
    x = "Mating-arena dimensionality (lifestyle proxy)",
    y = "Predicted trait value",
    caption = paste0(
      "In both panels, a 3D lifestyle raises the trait where the sexes share ",
      "little on a common territory, but not where they share resources ",
      "year-round.\nY-axes are trait-specific and not comparable in magnitude.")) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.caption = element_text(colour = "grey50", size = 8),
        strip.text = element_text(face = "bold"),
        legend.position = "right", panel.grid.minor = element_blank())

ggsave(file.path(dir_out, "fig2_interaction.pdf"), fig2, width = 9, height = 4.4)
cat("  saved fig2_interaction.pdf\n")

################################################################################
## FIGURE 3: within-Passeriformes replication of Figure 1 ---------------------
################################################################################
cat("\n=== Figure 3 ===\n")

## panel A: harem + RDP (lek not estimable within Passeriformes)
figA_pass_df <- bind_rows(
  summ_logodds(res_marc_rdp_vs_mono_p, "Resource-defense polygamy vs monogamy (Marcondes)"),
  summ_logodds(res_barber_poly_p,      "Strong polygamy vs monogamy (Barber)")
)
if (nrow(figA_pass_df) == 0)
  figA_pass_df <- data.frame(label = "(no estimable mating-system model)",
                             median = 0, lo = 0, hi = 0)
figA_pass_df <- figA_pass_df %>%
  mutate(label = factor(label, levels = rev(c(
           "Resource-defense polygamy vs monogamy (Marcondes)",
           "Strong polygamy vs monogamy (Barber)"))),
         couleur = .fig_colour(median, lo, hi))
figA_pass <- make_forest(figA_pass_df,
                         "A. Mating system (Passeriformes only; monogamy-referenced; lek not estimable)",
                         "Log-odds (3D vs 2D)")

## panel B: z-scored morphology + behaviour, Passeriformes only
dat_z_pass <- merged %>%
  filter(order == "Passeriformes", !is.na(dim_bin), !is.na(tip_label)) %>%
  distinct(tip_label, .keep_all = TRUE) %>%
  mutate(across(all_of(cont_traits), ~ as.numeric(scale(.)))) %>%
  as.data.frame()
rownames(dat_z_pass) <- dat_z_pass$tip_label

dat_z_avo_pass <- avo_base %>%
  filter(order == "Passeriformes", !is.na(dim_bin), !is.na(dichromatism),
         !is.na(tip_label)) %>%
  distinct(tip_label, .keep_all = TRUE) %>%
  mutate(across(c(dichromatism, male_plumage, female_plumage),
                ~ as.numeric(scale(.)))) %>%
  as.data.frame()
rownames(dat_z_avo_pass) <- dat_z_avo_pass$tip_label

dat_z_barber_pass <- avo_base %>%
  filter(order == "Passeriformes", !is.na(dim_bin), !is.na(barber_ss),
         !is.na(tip_label)) %>%
  distinct(tip_label, .keep_all = TRUE) %>%
  mutate(barber_ss = as.numeric(scale(barber_ss))) %>%
  as.data.frame()
rownames(dat_z_barber_pass) <- dat_z_barber_pass$tip_label

figB_pass_df <- bind_rows(
  summ_z("ssd_mass",       "Mass dimorphism",           dat_z_pass),
  summ_z("ssd_tarsus",     "Tarsus dimorphism",         dat_z_pass),
  summ_z("ssd_wing",       "Wing dimorphism",           dat_z_pass),
  summ_z("ssd_bill",       "Bill dimorphism",           dat_z_pass),
  summ_z("ssd_tail",       "Tail dimorphism",           dat_z_pass),
  summ_z("display_num",    "Display agility",            dat_z_pass),
  summ_z("dichromatism",   "Plumage dichromatism (M-F)", dat_z_avo_pass),
  summ_z("male_plumage",   "Male plumage elaboration",   dat_z_avo_pass),
  summ_z("female_plumage", "Female plumage elaboration", dat_z_avo_pass),
  summ_z("barber_ss",      "Sexual-selection intensity (Barber)", dat_z_barber_pass)
) %>%
  mutate(label = factor(label, levels = rev(c(
           "Mass dimorphism", "Tarsus dimorphism", "Tail dimorphism",
           "Wing dimorphism", "Bill dimorphism",
           "Display agility",
           "Plumage dichromatism (M-F)",
           "Male plumage elaboration", "Female plumage elaboration",
           "Sexual-selection intensity (Barber)"))),
         couleur = .fig_colour(median, lo, hi))
figB_pass <- make_forest(
  figB_pass_df, "B. Morphology and behaviour (Passeriformes only, z-scored)",
  "Standardised coefficient (3D vs 2D)",
  caption = paste0("Within-Passeriformes replication of Fig. 1. ",
                   "Vivid hue = across-tree range excludes 0."))

save_two_panels(figA_pass, figB_pass, file.path(dir_out, "fig3_passeriformes.pdf"),
                "Within-Passeriformes replication of Fig. 1")
cat("  saved fig3_passeriformes.pdf\n")

cat("\n06_figures.R done.\n")
