################################################################################
## 00_setup.R
## Dimensionnalite de l'arene de reproduction & forme de la selection sexuelle.
## Chemins, librairies, fonctions utilitaires et constantes PARTAGEES par tous
## les modules. Ce fichier ne fait AUCUN calcul : il est source en premier par
## run_all.R, puis 01..06 s'appuient sur les objets definis ici.
################################################################################

## --------------------------------------------------------------- librairies ---
library(readr)
library(readxl)
library(data.table)
library(dplyr)
library(stringr)
library(janitor)
library(tidyr)
library(ape)
library(phylolm)
library(ggplot2)
## (fuzzyjoin : plus necessaire -- le diagnostic de matching a ete retire)

## ------------------------------------------------------------------ chemins ---
## Arborescence du depot :
##   <racine>/data    : toutes les donnees d'entree (cf. README pour le manifeste)
##   <racine>/output  : toutes les sorties (figures + tables)
## La racine est detectee via {here} si disponible, sinon le repertoire courant.
if (requireNamespace("here", quietly = TRUE)) {
  dir_root <- here::here()
} else {
  dir_root <- getwd()
  message("Package {here} absent : racine = getwd() = ", dir_root,
          "\n  -> lancer R depuis la racine du projet, ou installer {here}.")
}
dir_data <- file.path(dir_root, "data")
dir_out  <- file.path(dir_root, "output")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)

## --- chemins des fichiers d'entree (noms nettoyes ; cf. manifeste README) -----
path_lislevand <- file.path(dir_data, "lislevand_2007_avian_ssd.txt")
path_avonet    <- file.path(dir_data, "avonet3_birdtree.xlsx")
avonet_sheet   <- "AVONET3_BirdTree"
path_blio      <- file.path(dir_data, "birdtree_taxonomy.csv")   # BLIOCPhyloMasterTax
path_dale      <- file.path(dir_data, "dale_2015_plumage_scores.csv")
path_marcondes <- file.path(dir_data, "marcondes_douvas_2024_mating_systems.xlsx")
marcondes_sheet<- "Species_data"
path_spur      <- file.path(dir_data, "menezes_palaoro_2022_spurs.csv")
## Arbre VertLife (1000 arbres Hackett, 9993 tips). Fichier volumineux :
## NON versionne sur GitHub (cf. .gitignore) -- voir README pour le telechargement.
path_master    <- file.path(dir_data, "AllBirdsHackett1.tre")

## --- parametres globaux -------------------------------------------------------
n_trees_use   <- 100   # nombre d'arbres echantillonnes parmi les 1000
SEP_THRESHOLD <- 8     # seuil de (quasi-)separation pour les binaires rares

################################################################################
## UTILITAIRES DE DONNEES ------------------------------------------------------
################################################################################

## normalisation d'un binome (cle de jointure) : minuscules, espaces propres
clean_binom <- function(x) {
  x %>% str_replace_all("_", " ") %>% str_squish() %>% str_to_lower()
}

## 1re colonne dont le nom matche un motif (lecture robuste d'AVONET)
get_col <- function(df, pattern) {
  hit <- names(df)[str_detect(names(df), regex(pattern, ignore_case = TRUE))]
  if (length(hit)) hit[1] else NA_character_
}

################################################################################
## MOTEUR PGLS -----------------------------------------------------------------
## Toutes les fonctions supposent :
##   - un data.frame indexe par tip_label (rownames = feuilles de l'arbre) ;
##   - l'objet global `trees` (liste de 100 arbres), cree par 03_trees.R.
## Deux familles :
##   .fit_binary   -> phyloglm logistique  (reponses binaires : harem, lek, ...)
##   .fit_gaussian -> phylolm lambda       (reponses continues : SSD, dichrom, ...)
################################################################################

## prepare un data.frame indexe par tip_label, filtre sur la reponse non-NA.
## `data` par defaut = `merged` (cree par 02_build_data.R).
.prep_dat <- function(data = merged, response, extra_filter = NULL) {
  d <- data %>%
    filter(!is.na(dim_bin), !is.na(.data[[response]]), !is.na(tip_label))
  if (!is.null(extra_filter)) d <- d %>% filter(!!extra_filter)
  d <- d %>% distinct(tip_label, .keep_all = TRUE) %>% as.data.frame()
  rownames(d) <- d$tip_label
  d
}

## --- PGLS binaire (phyloglm) sur 1 arbre -------------------------------------
## Renvoie NULL si phyloglm echoue OU ne converge pas (warning de convergence),
## pour eviter de rapporter des coefficients issus de modeles instables.
.fit_binary <- function(tree, response, dat, formula = NULL) {
  keep <- intersect(tree$tip.label, rownames(dat))
  if (length(keep) < 10) return(NULL)
  f <- if (is.null(formula)) as.formula(paste(response, "~ dim_bin")) else formula
  fit <- tryCatch({
    tr <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
    d  <- dat[tr$tip.label, ]
    if (sum(d[[response]], na.rm = TRUE) < 5) return(NULL)
    withCallingHandlers(
      phylolm::phyloglm(f, data = d, phy = tr,
                        method = "logistic_MPLE", btol = 30),
      warning = function(w) {
        if (grepl("converge", conditionMessage(w), ignore.case = TRUE))
          invokeRestart("muffleWarning")   # avale le warning, gere via $convergence
      })
  }, error = function(e) NULL)
  if (!is.null(fit) && !is.null(fit$convergence) && fit$convergence != 0)
    return(NULL)
  fit
}

## --- PGLS continue (phylolm lambda) sur 1 arbre ------------------------------
## `drop_na_response` : si TRUE, retire d'abord les especes NA sur la reponse
## (utile quand `dat` porte plusieurs reponses z-scorees simultanement).
.fit_gaussian <- function(tree, response, dat, formula = NULL,
                          drop_na_response = FALSE) {
  if (drop_na_response) dat <- dat[!is.na(dat[[response]]), , drop = FALSE]
  keep <- intersect(tree$tip.label, rownames(dat))
  if (length(keep) < 20) return(NULL)
  f <- if (is.null(formula)) as.formula(paste(response, "~ dim_bin")) else formula
  tryCatch({
    tr <- ape::drop.tip(tree, setdiff(tree$tip.label, keep))
    d  <- dat[tr$tip.label, ]
    phylolm::phylolm(f, data = d, phy = tr, model = "lambda")
  }, error = function(e) NULL)
}

## --- boucle sur les 100 arbres + resume d'un coefficient ---------------------
## Sortie : mediane + IC across-trees, % de signe predit, % p<0.05 / p<0.10
## (Wald par arbre), et p combinee a la Rubin (intra + inter-arbres).
## Alimente directement la Table A2 (robustesse) du manuscrit.
run_pgls <- function(response, type = "gaussian",
                     coef_name = "dim_bintrois_D",
                     direction = c("negatif", "positif"),
                     formula = NULL, dat = NULL, label = response) {
  direction <- match.arg(direction)
  if (is.null(dat)) dat <- .prep_dat(merged, response)
  fitter <- if (type == "binary") .fit_binary else .fit_gaussian
  
  fits <- lapply(trees, fitter, response = response, dat = dat, formula = formula)
  ok   <- Filter(Negate(is.null), fits)
  
  cat(sprintf("\n== %-22s | N=%d | arbres OK=%d/%d | dir.predite=%s ==\n",
              label, nrow(dat), length(ok), length(trees), direction))
  if (length(ok) == 0) { cat("  Aucun arbre converge.\n"); return(invisible(NULL)) }
  if (length(ok) < 80)
    cat(sprintf("  /!\\ seulement %d/%d arbres converges : resultat NON FIABLE",
                length(ok), length(trees)),
        "(instabilite / quasi-separation), a ne pas rapporter.\n")
  
  coefs <- sapply(ok, function(f) coef(f)[coef_name])
  ses   <- sapply(ok, function(f) sqrt(diag(vcov(f)))[coef_name])
  sgn   <- if (direction == "negatif") -1 else 1
  
  ## p de Wald (deux queues) par arbre
  pvals     <- 2 * pnorm(-abs(coefs / ses))
  pct_sig   <- round(100 * mean(pvals < 0.05, na.rm = TRUE))
  pct_sig10 <- round(100 * mean(pvals < 0.10, na.rm = TRUE))
  
  ## combinaison a la Rubin (1987) : W intra-arbre, B inter-arbres,
  ## T = W + (1 + 1/M) B -> un seul z-test / p combinee.
  M        <- length(ok)
  W        <- mean(ses^2, na.rm = TRUE)
  B        <- stats::var(coefs)
  T_var    <- W + (1 + 1/M) * B
  se_rubin <- sqrt(T_var)
  z_rubin  <- mean(coefs) / se_rubin
  p_rubin  <- 2 * pnorm(-abs(z_rubin))
  
  cat(sprintf("  mediane = %+.4f | IC95%% = [%+.4f, %+.4f] | %d%% dans direction | %d%% p<0.05 | %d%% p<0.10 (Wald)\n",
              median(coefs), quantile(coefs, 0.025), quantile(coefs, 0.975),
              round(100 * mean(coefs * sgn > 0)), pct_sig, pct_sig10))
  cat(sprintf("  Rubin combine : SE=%.4f | z=%.3f | p=%.4g\n",
              se_rubin, z_rubin, p_rubin))
  
  invisible(list(coefs = coefs, ses = ses, pvals = pvals,
                 pct_sig = pct_sig, pct_sig10 = pct_sig10,
                 se_rubin = se_rubin, z_rubin = z_rubin, p_rubin = p_rubin,
                 n = nrow(dat), direction = direction, label = label))
}

## --- controle generique par UNE covariable continue z-scoree -----------------
## Generalise le controle de taille (Rensch) a n'importe quelle covariable
## (taille, HWI, ...). Retourne l'effet 3D brut, l'effet 3D controle, et imprime
## le coefficient de la covariable elle-meme.
## Alimente la Table A3 (Rensch, covar = body_size_log) et la Table A4 (HWI).
run_pgls_ctrl <- function(response, covar, direction, label,
                          type = "gaussian", data = merged) {
  base <- data %>%
    filter(!is.na(dim_bin), !is.na(.data[[response]]),
           !is.na(tip_label), !is.na(.data[[covar]])) %>%
    distinct(tip_label, .keep_all = TRUE)
  covz <- paste0(covar, "_z")
  base[[covz]] <- as.numeric(scale(base[[covar]]))
  dat <- as.data.frame(base); rownames(dat) <- dat$tip_label
  
  res_raw <- run_pgls(response, type = type, direction = direction, dat = dat,
                      label = paste0(label, " [brut]"))
  res_ctl <- run_pgls(
    response, type = type, direction = direction, coef_name = "dim_bintrois_D",
    formula = as.formula(paste(response, "~ dim_bin +", covz)),
    dat = dat, label = paste0(label, " [+ ", covar, "]"))
  
  ## coefficient de la covariable elle-meme (sens attendu : Rensch>0 ; HWI<0)
  if (!is.null(res_ctl)) {
    fitter <- if (type == "binary") .fit_binary else .fit_gaussian
    fits <- lapply(trees, fitter, response = response, dat = dat,
                   formula = as.formula(paste(response, "~ dim_bin +", covz)))
    okf  <- Filter(Negate(is.null), fits)
    if (length(okf)) {
      cc <- sapply(okf, function(f) coef(f)[covz])
      cat(sprintf("   -> coef %s : mediane = %+.4f [%+.4f, %+.4f]\n",
                  covar, median(cc), quantile(cc, .025), quantile(cc, .975)))
    }
  }
  invisible(list(raw = res_raw, ctl = res_ctl))
}

################################################################################
## HELPERS BINAIRES RARES (spurs, systemes Marcondes) --------------------------
## Table de contingence 2D/3D, garde-fou de separation, test exact de Fisher.
################################################################################

## table 2x2 (dim_bin x {0,1}) pour une variable binaire `var` d'un jeu `data`.
xtab_dim <- function(data, var) {
  d <- data %>%
    filter(dim_bin %in% c("deux_D", "trois_D"), !is.na(.data[[var]])) %>%
    distinct(key, .keep_all = TRUE) %>%
    as.data.frame()
  d$dim_bin <- droplevels(d$dim_bin)
  table(dim_bin = d$dim_bin, val = d[[var]])
}

## effectif de la categorie positive la plus rare (2D vs 3D) dans une table 2x2.
min_pos_cell <- function(tab) {
  if (!("1" %in% colnames(tab))) return(0)
  min(tab[, "1"])
}

## test exact de Fisher sur la table 2x2, avec impression formatee.
fisher_dim <- function(tab, label = "") {
  ft <- fisher.test(tab)
  cat(sprintf("  [%s] Fisher exact : OR=%.3f, p=%.3g\n",
              label, unname(ft$estimate), ft$p.value))
  ft
}

cat("00_setup.R charge : chemins, helpers et constantes prets.\n")