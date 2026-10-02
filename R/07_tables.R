################################################################################
## 07_tables.R
## Paper-ready tables, written to output/tables/ as ready-to-paste LaTeX table
## environments (booktabs; caption + label included) and as CSV:
##
##   Table 1  : descriptive values of every Block A response by foraging
##              lifestyle (variables in rows, lifestyles in columns)  [main text]
##   Table S1 : Block A main effects, full statistics                 [appendix]
##   Table S2 : Block A with controls / on subsets                    [appendix]
##   Table S3 : Block B2, dim_bin x continuous parental care          [appendix]
##
## Needs the model cache (output/model_results.rds) and, for Table 1, the
## dataset (03_load_dataset.R). LaTeX requirements: booktabs only.
################################################################################

if (!exists("model_results")) model_results <- readRDS(path_model_cache)
mr <- model_results
specA <- mr$specA

## ---------------------------------------------------------------- helpers ---
## decimals adapted to the precision of the estimate (from its SE)
dec_from_se <- function(se) pmax(1, pmin(5, 1 - floor(log10(se))))
fmt <- function(x, d) {
  x <- round(x, d); x[!is.na(x) & x == 0] <- 0      # no "-0.000"
  ifelse(is.na(x), "---", gsub("^-", "$-$", formatC(x, format = "f", digits = d)))
}
fmt_ci <- function(lo, hi, d) sprintf("[%s, %s]", fmt(lo, d), fmt(hi, d))
fmt_p <- function(p) ifelse(is.na(p), "---",
                     ifelse(p < 0.001, "$<$0.001", formatC(p, format = "f", digits = 3)))
fmt_n <- function(n) formatC(n, format = "d", big.mark = ",")

## superscript letter per data source, in order of first appearance
src_used    <- unique(specA$source)
src_letters <- setNames(letters[seq_along(src_used)], src_used)
src_sup     <- function(src) sprintf("$^{\\mathrm{%s}}$", src_letters[src])
src_notes   <- paste0(
  "Sources: ",
  paste(sprintf("$^{\\mathrm{%s}}$%s", src_letters,
                source_info$cite[match(names(src_letters), source_info$source)]),
        collapse = "; "), ".")

## write a LaTeX table environment (booktabs) + a CSV twin
write_tex_table <- function(name, header, body, caption, label, colspec,
                            notes = NULL, size = "\\footnotesize", csv = NULL,
                            tabcolsep = "3pt") {
  tex <- c("\\begin{table}[htbp]", "\\centering",
           sprintf("\\caption{%s}", caption), sprintf("\\label{%s}", label),
           size, sprintf("\\setlength{\\tabcolsep}{%s}", tabcolsep),
           sprintf("\\begin{tabular}{%s}", colspec), "\\toprule",
           header, "\\midrule", body, "\\bottomrule", "\\end{tabular}")
  if (!is.null(notes))
    tex <- c(tex, "", sprintf("\\par\\smallskip\\parbox{\\linewidth}{\\scriptsize %s}", notes))
  tex <- c(tex, "\\end{table}")
  writeLines(tex, file.path(dir_tab, paste0(name, ".tex")))
  if (!is.null(csv)) write_csv(csv, file.path(dir_tab, paste0(name, ".csv")))
  cat("  saved", name, "(.tex, .csv)\n")
}

## one table row; arguments may be scalars or vectors (cells, in order)
## wrap a (long) label cell into a ragged-right parbox (no extra package needed)
wrap <- function(x, width) sprintf("\\parbox[t]{%s}{\\raggedright %s}", width, x)

row_tex   <- function(...) paste0(paste(unlist(list(...)), collapse = " & "), " \\\\")
group_row <- function(title, ncol)
  sprintf("\\addlinespace\\multicolumn{%d}{l}{\\textit{%s}} \\\\", ncol, title)

block_title <- setNames(block_info$block_label, block_info$block)

################################################################################
## TABLE 1: descriptive values by foraging lifestyle ---------------------------
## Binary responses: % of species coded 1 (with the coding used in the models:
## mating-system contrasts are against monogamy). Continuous: mean (SD).
## Under each value, the number of species with data.
################################################################################
cat("\n=== Table 1 (descriptives by lifestyle) ===\n")

groups <- list(
  "Terrestrial"  = quote(primary_lifestyle == "Terrestrial"),
  "Insessorial"  = quote(primary_lifestyle == "Insessorial"),
  "Aerial"       = quote(primary_lifestyle == "Aerial"),
  "All 3D"       = quote(dim_bin == "3D"),
  "Generalist"   = quote(primary_lifestyle == "Generalist"),
  "Aquatic"      = quote(primary_lifestyle == "Aquatic"))

## decimals of the descriptive values, per response
desc_digits <- c(barber_ss = 2, ssd_mass = 3, ssd_tarsus = 3, ssd_wing = 3,
                 ssd_bill = 3, ssd_tail = 3, dichromatism = 2, male_plumage = 1,
                 female_plumage = 1, dichro_sr = 2, uv_cd = 2, display_num = 2)

desc_cell <- function(s, g) {
  d <- get(s$base) %>% filter(!!groups[[g]]) %>% distinct(key, .keep_all = TRUE)
  x <- d[[s$response]]; x <- x[!is.na(x)]
  n <- length(x)
  if (n == 0) return(c(val = "---", n = "0", num = NA, sd = NA, nn = 0))
  if (s$type == "binary") {
    c(val = formatC(100 * mean(x), format = "f", digits = 1), n = fmt_n(n),
      num = 100 * mean(x), sd = NA, nn = n)
  } else {
    dg <- desc_digits[[s$response]]
    c(val = sprintf("%s (%s)", fmt(mean(x), dg), fmt(sd(x), dg)), n = fmt_n(n),
      num = mean(x), sd = sd(x), nn = n)
  }
}

t1_body <- character(0); t1_csv <- list()
for (b in block_info$block) {
  t1_body <- c(t1_body, group_row(block_title[[b]], 1 + length(groups)))
  sp <- specA %>% filter(block == b)
  for (i in seq_len(nrow(sp))) {
    s     <- sp[i, ]
    cells <- lapply(names(groups), function(g) desc_cell(s, g))
    unit  <- if (s$type == "binary") "\\%" else "mean (SD)"
    lab   <- wrap(sprintf("%s%s%s%s", s$label,
                          if (s$type == "binary") " (\\%)" else "",
                          src_sup(s$source),
                          if (s$separated) "$^{\\dagger}$" else ""), "4.3cm")
    t1_body <- c(t1_body,
                 row_tex(lab, sapply(cells, `[[`, "val")),
                 row_tex("\\quad{\\scriptsize $n$}",
                         sprintf("{\\scriptsize %s}", sapply(cells, `[[`, "n"))))
    t1_csv[[s$id]] <- data.frame(
      block = b, response = s$label, source = s$source, stat = unit,
      group = names(groups),
      value = as.numeric(sapply(cells, `[[`, "num")),
      sd    = as.numeric(sapply(cells, `[[`, "sd")),
      n     = as.integer(sapply(cells, `[[`, "nn")))
  }
}
## N species per lifestyle (AVONET)
n_tot <- sapply(names(groups), function(g)
  avo_base %>% filter(!!groups[[g]]) %>% distinct(key) %>% nrow())

t1_header <- c(
  "& \\multicolumn{1}{c}{2D} & \\multicolumn{3}{c}{3D} & \\multicolumn{2}{c}{Excluded from the contrast} \\\\",
  "\\cmidrule(lr){2-2}\\cmidrule(lr){3-5}\\cmidrule(lr){6-7}",
  row_tex("", names(groups)),
  row_tex("\\textit{Species (AVONET)}", sprintf("\\textit{%s}", fmt_n(n_tot))))

fisher_note <- paste0(
  "$^{\\dagger}$Rare trait with (quasi-)complete separation between 2D and 3D; ",
  "exact Fisher test (2D vs 3D, non-phylogenetic): ",
  paste(sapply(names(mr$fisherA), function(id) {
    f <- mr$fisherA[[id]]
    sprintf("%s, 2D %d/%d vs 3D %d/%d, $p = %s$",
            tolower(specA$label[specA$id == id]),
            f$tab["2D", "1"], sum(f$tab["2D", ]), f$tab["3D", "1"], sum(f$tab["3D", ]),
            gsub("e(-?)0*(\\d+)", "\\\\times10^{\\1\\2}", formatC(f$p, format = "e", digits = 1)))
  }), collapse = "; "), ".")

write_tex_table(
  "table1_descriptives", t1_header, t1_body,
  caption = paste0(
    "\\textbf{Response variables by AVONET primary foraging lifestyle.} ",
    "Binary responses (\\%): percentage of species coded 1 (mating-system contrasts use ",
    "the same coding as the models, i.e.\\ each polygamous system against monogamy); ",
    "continuous responses: mean (SD). Below each value, the number of species with ",
    "data ($n$). The focal contrast is 2D (Terrestrial) versus 3D (Insessorial + ",
    "Aerial); Generalist and Aquatic species are shown for completeness but excluded ",
    "from the contrast. Size dimorphism is $\\log$(male/female); plumage dichromatism ",
    "is the male minus female plumage score."),
  label = "tab:descriptive",
  colspec = "lrrrrrr",
  notes = paste(src_notes, fisher_note),
  size = "\\scriptsize",
  csv = bind_rows(t1_csv))

################################################################################
## TABLE S1: Block A main effects, full statistics -----------------------------
################################################################################
cat("\n=== Table S1 (Block A, full statistics) ===\n")

s1_body <- character(0); s1_csv <- list()
for (b in block_info$block) {
  s1_body <- c(s1_body, group_row(block_title[[b]], 9))
  sp <- specA %>% filter(block == b)
  for (i in seq_len(nrow(sp))) {
    s <- sp[i, ]
    r <- summarise_res(mr$resA[[s$id]])
    if (is.null(r)) next
    z <- summarise_res(mr$resA[[s$id]], std = TRUE)
    d <- dec_from_se(r$se)
    pred <- switch(s$direction, negative = "$-$", positive = "$+$", none = "none")
    flag <- if (s$separated) "$^{\\dagger}$" else ""
    if (r$n_trees < 80 || r$n_btol > 0) flag <- paste0(flag, "$^{\\ddagger}$")
    s1_body <- c(s1_body, row_tex(
      wrap(sprintf("%s%s%s", s$label, src_sup(s$source), flag), "3.4cm"),
      pred,
      sprintf("%s {\\tiny(%s/%s)}", fmt_n(r$N), fmt_n(r$N2D), fmt_n(r$N3D)),
      fmt(r$est, d), fmt_ci(r$lo95, r$hi95, d),
      if (s$type == "gaussian") fmt(z$est, 2) else "---",
      sprintf("%d", as.integer(r$pct_same)),
      fmt_p(r$p_one), fmt_p(r$p_two)))
    s1_csv[[s$id]] <- cbind(block = b, id = s$id, source = s$source,
                            r, est_std = if (s$type == "gaussian") z$est else NA)
  }
}
s1_header <- row_tex("Response", "Pred.", "$N$ (2D/3D)", "Estimate", "95\\% CI",
                     "Std.", "\\% trees", "$p_{1}$", "$p_{2}$")
write_tex_table(
  "tableS1_main_effects", s1_header, s1_body,
  caption = paste0(
    "\\textbf{Main (bivariate) effects of a three-dimensional lifestyle.} One ",
    "phylogenetic model per response (response $\\sim$ 3D), fitted on each of 100 ",
    "trees and pooled with Rubin's rules \\citep{rubin1987multiple,nakagawa2019general}. ",
    "Binary responses: penalised phylogenetic logistic regression ",
    "\\citep{ives2010phylogenetic}, estimate in log-odds; continuous responses: PGLS ",
    "with Pagel's $\\lambda$, estimate in response units, and in SD units of the ",
    "response (Std., as in Fig.~\\ref{fig:effects}). Pred.: predicted sign of the 3D ",
    "effect. \\% trees: percentage of trees in which the estimate has the sign of the ",
    "pooled estimate. $p_{1}$: one-sided pooled $p$ for the predicted direction ",
    "(values above 0.5 indicate an estimate opposite to the prediction); $p_{2}$: ",
    "two-sided pooled $p$."),
  label = "tab:mainfull",
  colspec = "lcrrrrrrr",
  notes = paste(src_notes,
                "$^{\\dagger}$(Quasi-)complete separation between 2D and 3D: the",
                "penalised estimate is finite but reflects the concentration of the",
                "trait in a few clades (see Fisher test in Table~\\ref{tab:descriptive}).",
                "$^{\\ddagger}$Fewer than 80 trees converged, or some fits reached the",
                "bound of the linear predictor: interpret with caution."),
  size = "\\scriptsize",
  csv = bind_rows(s1_csv))

################################################################################
## TABLE S2: Block A with controls / on subsets --------------------------------
################################################################################
cat("\n=== Table S2 (controls) ===\n")

ctrl_groups <- list(
  "Body-size control (Rensch's rule)"    = paste0(c("ssd_mass", "ssd_tarsus", "ssd_wing",
                                                    "ssd_bill", "ssd_tail"), "_rensch"),
  "Display agility"                      = c("display_noaerial", "display_hwi"),
  "Passeriformes only"                   = grep("_passer$", names(mr$resC), value = TRUE))

s2_body <- character(0); s2_csv <- list()
for (gname in names(ctrl_groups)) {
  s2_body <- c(s2_body, group_row(gname, 7))
  for (k in ctrl_groups[[gname]]) {
    e   <- mr$resC[[k]]
    s   <- specA[specA$id == e$id, ]
    ref <- summarise_res(e$ref); ctl <- summarise_res(e$ctl)
    lab <- wrap(if (s$block == "ss") sprintf("%s%s", s$label, src_sup(s$source)) else s$label,
                "2.5cm")
    ## the group title already names the subset for the Passeriformes rows
    ctl_lab <- if (gname == "Passeriformes only") "" else wrap(e$control, "2.3cm")
    if (is.null(ctl)) {
      s2_body <- c(s2_body, row_tex(lab, ctl_lab, "---", "---",
                                    "\\multicolumn{3}{c}{not estimable (separation)}"))
      next
    }
    d <- dec_from_se(ref$se)
    s2_body <- c(s2_body, row_tex(
      lab, ctl_lab, fmt_n(ctl$N),
      sprintf("%s %s", fmt(ref$est, d), fmt_ci(ref$lo95, ref$hi95, d)),
      sprintf("%s %s", fmt(ctl$est, d), fmt_ci(ctl$lo95, ctl$hi95, d)),
      sprintf("%d", as.integer(ctl$pct_same)),
      if (s$direction == "none") paste0(fmt_p(ctl$p_two), "$^{*}$") else fmt_p(ctl$p_one)))
    s2_csv[[k]] <- data.frame(group = gname, id = e$id, control = e$control,
                              N_ctl = ctl$N, est_ref = ref$est, lo95_ref = ref$lo95,
                              hi95_ref = ref$hi95, N_ref = ref$N,
                              est_ctl = ctl$est, lo95_ctl = ctl$lo95,
                              hi95_ctl = ctl$hi95, pct_same_ctl = ctl$pct_same,
                              p_one_ctl = ctl$p_one, p_two_ctl = ctl$p_two)
  }
}
s2_header <- row_tex("Response", "Control / subset", "$N$",
                     "Bivariate [95\\% CI]", "Controlled [95\\% CI]",
                     "\\% trees", "$p_{1}$")
write_tex_table(
  "tableS2_controls", s2_header, s2_body,
  caption = paste0(
    "\\textbf{Main effects of a three-dimensional lifestyle with controls and on ",
    "subsets.} Pooled 3D estimates (response units; log-odds for binary responses) ",
    "and 95\\% CIs across 100 trees. For covariate controls (body size, ",
    "hand-wing index), the bivariate estimate is refitted on the same ",
    "species as the controlled model; for subsets (excluding aerial species, ",
    "Passeriformes only), the bivariate estimate is the main effect ",
    "on the full sample (Table~\\ref{tab:mainfull}). $N$, \\% trees and $p_{1}$ ",
    "(one-sided, predicted direction) refer to the controlled / subset model."),
  label = "tab:controls",
  colspec = "llrrrrr",
  notes = paste(src_notes,
                "$^{*}$No directional prediction: two-sided $p$.",
                "Plumage scores from \\citet{dale2015effects} cover passerines only, so",
                "their Passeriformes-only estimates equal the main effects."),
  size = "\\scriptsize",
  csv = bind_rows(s2_csv))

################################################################################
## TABLE S3: Block B2, dim_bin x continuous parental care -----------------------
################################################################################
cat("\n=== Table S3 (care interaction) ===\n")

spB <- mr$specB2
s3_body <- character(0); s3_csv <- list()
for (b in c("contest", "choice")) {
  s3_body <- c(s3_body, group_row(block_title[[b]], 8))
  sp <- spB %>% filter(block == b)
  for (i in seq_len(nrow(sp))) {
    s <- sp[i, ]
    r <- summarise_res(mr$resB2[[s$id]])
    if (is.null(r)) next
    z <- summarise_res(mr$resB2[[s$id]], std = TRUE)
    d <- dec_from_se(r$se)
    s3_body <- c(s3_body, row_tex(
      wrap(sprintf("%s%s", s$label, src_sup(s$source)), "4.2cm"),
      if (s$direction == "negative") "$-$" else "$+$",
      fmt_n(r$N), fmt(r$est, d), fmt_ci(r$lo95, r$hi95, d),
      if (s$type == "gaussian") fmt(z$est, 2) else "---",
      sprintf("%d", as.integer(r$pct_same)), fmt_p(r$p_one)))
    s3_csv[[s$id]] <- cbind(block = b, id = s$id, r,
                            est_std = if (s$type == "gaussian") z$est else NA)
  }
}
s3_header <- row_tex("Response", "Pred.", "$N$", "3D $\\times$ care", "95\\% CI",
                     "Std.", "\\% trees", "$p_{1}$")
write_tex_table(
  "tableS3_care_interaction", s3_header, s3_body,
  caption = paste0(
    "\\textbf{Modulation of the 3D effect by parental-care asymmetry.} Interaction ",
    "coefficient from response $\\sim$ 3D $\\times$ care, where care is the relative ",
    "care investment of the sexes ($z$-scored; positive = female-biased). Under the ",
    "channel-switch hypothesis, female-biased care frees males to compete, so the ",
    "interaction should be negative for contest traits and positive for choice traits ",
    "(Pred.). Estimates in response units (log-odds for allopreening) and, for ",
    "continuous responses, in SD units of the response (Std., as in ",
    "Fig.~\\ref{fig:care}). Pooling, \\% trees and $p_{1}$ as in ",
    "Table~\\ref{tab:mainfull}."),
  label = "tab:care",
  colspec = "lcrrrrrr",
  notes = src_notes,
  size = "\\scriptsize",
  csv = bind_rows(s3_csv))

cat("\n07_tables.R done. Tables in", dir_tab, "\n")
