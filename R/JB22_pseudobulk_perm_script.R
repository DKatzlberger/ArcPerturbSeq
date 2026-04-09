## max. code width ============================================================

## Libraries =================================================================

# Scource R scripts
source("R/utils.R")
source("R/JB22_pseudobulk_utils.R")

# Activate renv + load libs.
.script_startup()

## Directories ================================================================

# Name of experiment
exp_name <- "JB22_pseudobulk_perm"

# Define paths
res_path <- file.path("results", exp_name)
fig_path <- file.path(res_path, "fig")
dge_path <- file.path(res_path, "dge")
tab_path <- file.path(res_path, "tab")

# Create output dirs
dir.create(res_path, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)
dir.create(dge_path, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_path, recursive = TRUE, showWarnings = FALSE)

# Shared directory (location: vscratch)
JB22_pb_dir  <- "/vscratch/wes/arc/share"


## Dataset ====================================================================

# Pseudobulk specs
JBB_pb_specs <- list(
  Peritoneum = list(
    base_dir = file.path(JB22_pb_dir, "JB22_PM"),
    files = c(
      Macrophages = "pseudobulk_Macrophages.h5ad"
    )
  ),
  Liver = list(
    base_dir = file.path(JB22_pb_dir, "JB22_Liver"),
    files = c(
      Bcells = "pseudobulk_B_cells.h5ad",
      Macrophages = "pseudobulk_Monocytes.h5ad",
      Neutrophils = "pseudobulk_Neutrophils.h5ad"
    )
  ),
  Spleen = list(
    base_dir = file.path(JB22_pb_dir, "JB22_Spleen"),
    files = c(
      Bcells = "pseudobulk_B_cells.h5ad"
    )
  )
)

# Experiment JB22 | SingleCellExperiment obj.
JB22_pb_sce <- .combine_pseudobulk(
  JBB_pb_specs, 
  exp_name = exp_name,
  output = "sce"
)

## JB22 bulk RNA-seq (perm.) ==================================================

dge_methods <- list(
  edgeR_QLF = .dge_fit_edgeR_QLF,
  DESeq2_Wald = .dge_fit_DESeq2_Wald,
  limma_voom = .dge_fit_limma_voom
)

dge_formulas <- list(
  genotype = ~ genotype,
  genotype_sample = ~ genotype + sample
)

.build_tasks <- function(
  sce_obj,
  formula,
  method,
  strata,
  n_perm = 10
) {

  meta <- as.data.table(colData(sce_obj))

  # Check strata
  miss <- setdiff(strata, names(meta))
  if (length(miss)) {
    stop("Missing strata columns: ", paste(miss, collapse = ", "))
  }

  # Unique strata + index
  combos <- unique(meta[, ..strata])[, row_id := .I]

  # Full expansion
  tasks <- CJ(
    row_id  = combos$row_id,
    formula = names(formula),
    method  = names(method),
    perm_id = 0:n_perm,
    sorted = FALSE
  )[combos, on = "row_id"][, row_id := NULL]

  # Add shuffled + seed
  tasks[, `:=`(
    shuffle = perm_id > 0,
    seed = {
      f_id <- match(formula, names(formula))
      m_id <- match(method,  names(method))
      s <- sprintf("%08d", perm_id * 1e6 + f_id * 1e4 + m_id * 1e2)
      s[perm_id == 0] <- sprintf("%0*d", nchar(s[1]), 0)
      s
    }
  )][, perm_id := NULL]

  # Return
  order <- c("tissue", "cell_type", "formula", "method", "shuffle", "seed")
  return(tasks[, ..order])
}



generate_permutations <- function(
  sce_obj,
  formula,
  strata,
  ref_level = NULL,
  n_perm = 10,
  verbose = TRUE
) {

  meta <- as.data.table(colData(sce_obj))

  # Check strata
  miss <- setdiff(strata, colnames(meta))
  if (length(miss)) {
    stop("Missing strata columns: ", paste(miss, collapse = ", "))
  }

  # Split once by strata
  split_keys <- unique(meta[, ..strata])
  n_strata <- nrow(split_keys)

  subset_list <- vector("list", n_strata)

  for (k in seq_len(n_strata)) {
    key <- split_keys[k]

    idx <- rep(TRUE, ncol(sce_obj))
    for (s in strata) {
      idx <- idx & (meta[[s]] == key[[s]])
    }

    subset_list[[k]] <- sce_obj[, idx]
  }

  names(subset_list) <- apply(split_keys, 1, paste, collapse = "__")

  # Outputs
  results <- list()
  perm_info <- list()
  counter <- 1

  # Main loop
  for (k in seq_len(n_strata)) {

    key <- split_keys[k]
    sce_sub <- subset_list[[k]]

    if (verbose) {
      message(sprintf(
        "\n[Stratum %d/%d] %s",
        k, n_strata,
        paste(key, collapse = " | ")
      ))
    }

    for (f_name in names(formula)) {

      form <- formula[[f_name]]
      f_id <- match(f_name, names(formula))

      if (verbose) {
        message(sprintf("  Formula: %s", f_name))
      }

      for (p in 0:n_perm) {

        shuffle <- (p > 0)

        # Deterministic seed
        if (shuffle) {
          seed_int <- p + f_id * 1e4 + k * 1e6
        } else {
          seed_int <- 0
        }

        set.seed(seed_int)
        seed_chr <- sprintf("%08d", seed_int)

        if (verbose) {
          message(sprintf("   Perm %d (seed=%s)", p, seed_chr))
        }

        # Permute
        sce_perm <- .prepare_metadata(
          sce_obj = sce_sub,
          formula = form,
          ref_level = ref_level,
          shuffle = shuffle,
          verbose = FALSE
        )

        # ID
        id <- paste(
          c(as.character(key),
            f_name,
            paste0("perm_", p)),
          collapse = "__"
        )

        # Store
        results[[id]] <- sce_perm

        perm_info[[counter]] <- cbind(
          key,
          data.table(
            formula = f_name,
            shuffle = shuffle,
            seed = seed_chr,
            id = id
          )
        )

        counter <- counter + 1
      }
    }
  }

  perm_info <- rbindlist(perm_info)

  return(list(
    perm_list = results,
    perm_info = perm_info
  ))
}

perm_list <- generate_permutations(
  sce_obj = JB22_pb_sce,
  formula = dge_formulas,
  strata = c("tissue", "cell_type"),
  n_perm = 10
)

perm_list$perm_info


factor_encoding <- function(
  formula,
  data,
  ref_levels = list(),
  verbose = TRUE
){
  dt <- data.table::as.data.table(data.table::copy(data))
  
  vars <- all.vars(formula)
  summary_list <- list()
  
  for (v in vars) {
    if (is.character(dt[[v]]) || is.factor(dt[[v]])) {
      
      dt[, (v) := factor(get(v))]
      
      original_levels <- levels(dt[[v]])
      ref_used <- original_levels[1]
      ref_type <- "default"
      
      if (v %in% names(ref_levels)) {
        if (!(ref_levels[[v]] %in% original_levels)) {
          stop(sprintf("Reference level '%s' not found in %s", ref_levels[[v]], v))
        }
        dt[, (v) := stats::relevel(get(v), ref = ref_levels[[v]])]
        ref_used <- ref_levels[[v]]
        ref_type <- "user"
      }
      
      summary_list[[v]] <- list(
        levels = levels(dt[[v]]),
        ref = ref_used,
        type = ref_type
      )
    }
  }
  
  # Verbose summary
  if (isTRUE(verbose) && length(summary_list) > 0) {
    
    cat("Factor encoding:\n")
    
    for (v in names(summary_list)) {
      info <- summary_list[[v]]
      levs <- info$levels
      n_lev <- length(levs)
      
      # Show max 3 levels
      if (n_lev > 3) {
        lev_preview <- paste(c(levs[1:2], "..."), collapse = ", ")
      } else {
        lev_preview <- paste(levs, collapse = ", ")
      }
      
      cat(sprintf(
        "  - '%s': ref = '%s' (%s); levels (n=%d) = [%s]\n",
        v,
        info$ref,
        info$type,
        n_lev,
        lev_preview
      ))
    }
    
    cat("\n")
  }
  
  return(dt[])
}

conditional_permutation <- function(
  formula,
  data,
  permute,
  within,
  verbose = TRUE
){
  dt <- data.table::as.data.table(data.table::copy(data))
  
  
  # Permute must exist
  if (!permute %in% names(dt)) {
    stop(sprintf("Variable '%s' not found in data", permute))
  }
  
  # Within must be provided
  if (is.null(within) || length(within) == 0) {
    stop("'within' must be specified explicitly")
  }
  
  # Within variables must exist
  missing_within <- setdiff(within, names(dt))
  if (length(missing_within) > 0) {
    stop(sprintf(
      "Variables not found in data: %s",
      paste(missing_within, collapse = ", ")
    ))
  }
  
  # Ensure factor
  if (!is.factor(dt[[permute]])) {
    stop(sprintf("Variable '%s' must be a factor", permute))
  }
  
  # Store levels (CRITICAL for model consistency)
  lvls <- levels(dt[[permute]])
  
  # Verbose message (H0 explanation)
  if (isTRUE(verbose)) {
    cat("Permutation scheme:\n")
    cat(sprintf("  - Variable permuted: '%s'\n", permute))
    cat(sprintf("  - Within strata: '%s'\n", paste(within, collapse = " × ")))
    cat(sprintf("  - Implied H0: '%s' cannot show true associations with response\n",
      permute
    ))
    cat("\n")
  }
  
  # Perform permutation
  dt[, (permute) := sample(get(permute)), by = within]
  
  # Restore factor levels
  dt[, (permute) := factor(get(permute), levels = lvls)]
  
  # Check permutation effectiveness
  chk <- dt[, data.table::uniqueN(get(permute)), by = within]
  if (any(chk$V1 < 2)) {
    warning("Some strata contain only one level — permutation ineffective there")
  }
  
  # Return
  return(dt[])
}

build_model <- function(
  formula,
  data,
  verbose = TRUE
){
  # Build model matrix
  mm <- model.matrix(formula, data = data)
  
  # Dimensions
  p <- ncol(mm)

  # Verbose message
  if (isTRUE(verbose)) {
    
    cat("Model design:\n")
    cat(sprintf("  - Total coefficients: %d\n", p))
    cat("\n")
  }
  
  return(mm)
}

explain_coefs <- function(
  formula,
  data,
  permute = NULL,
  verbose = TRUE
){

  # Build model
  mm <- build_model(
    formula = formula, 
    data = data,
    verbose = FALSE
  )

  mm_names <- colnames(mm)
  
  vars <- all.vars(formula)
  factor_vars <- vars[sapply(data[, ..vars], is.factor)]
  

  # Helper: parse coef
  parse_coef <- function(coef) {
    for (v in factor_vars) {
      if (startsWith(coef, v)) {
        lvl <- sub(paste0("^", v), "", coef)
        return(list(var = v, level = lvl))
      }
    }
    return(NULL)
  }
  
  # Build result table
  dt_res <- data.table::data.table(
    coef = mm_names
  )
  
  dt_res[, type := data.table::fifelse(
    coef == "(Intercept)", "intercept",
    data.table::fifelse(grepl(":", coef), "interaction", "main")
  )]
  
  dt_res[, permuted := if (!is.null(permute)) {
    grepl(permute, coef)
  } else FALSE]
  
  # Comparison
  dt_res[, comparison := sapply(coef, function(cf) {
    
    if (cf == "(Intercept)") {
      return("baseline")
    }
    
    # Main effect
    if (!grepl(":", cf)) {
      parsed <- parse_coef(cf)
      
      if (!is.null(parsed)) {
        ref <- levels(data[[parsed$var]])[1]
        return(sprintf("%s - %s", parsed$level, ref))
      } else {
        return(NA_character_)
      }
    }
    
    # Interaction
    parts <- strsplit(cf, ":")[[1]]
    parsed_parts <- lapply(parts, parse_coef)
    
    if (all(!sapply(parsed_parts, is.null))) {
      
      diffs <- sapply(parsed_parts, function(x) {
        ref <- levels(data[[x$var]])[1]
        sprintf("(%s - %s)", x$level, ref)
      })
      
      if (length(diffs) == 2) {
        return(sprintf(
          "%s in %s vs %s",
          diffs[1],
          parsed_parts[[2]]$level,
          levels(data[[parsed_parts[[2]]$var]])[1]
        ))
      }
      
      return(sprintf("interaction: %s", paste(diffs, collapse = " × ")))
    }
    
    return("interaction effect")
  })]
  
  # Verbose
  if (isTRUE(verbose)) {
    
    cat("Model:\n")
    
    # Intercept (only printed, not returned)
    if ("(Intercept)" %in% mm_names) {
      cat("  - Intercept:\n")
      for (v in factor_vars) {
        ref_level <- levels(data[[v]])[1]
        cat(sprintf("      %s: %s\n", v, ref_level))
      }
    }
    
    terms_dt <- dt_res[coef != "(Intercept)"]
    
    # Formatting helpers
    label_flag <- function(permuted) {
      if (permuted) "[assoc. broken under H0]" else "[structure preserved]"
    }
    
    # dynamic widths
    coef_width <- max(nchar(terms_dt$coef), na.rm = TRUE)
    comp_width <- max(nchar(terms_dt$comparison), na.rm = TRUE)
    
    fmt <- paste0("      %-", coef_width, "s → %-", comp_width, "s %s\n")
    
    # Main effects
    main_dt <- terms_dt[type == "main"]
    if (nrow(main_dt) > 0) {
      cat("  - Main effects:\n")
      for (i in seq_len(nrow(main_dt))) {
        cat(sprintf(
          fmt,
          main_dt$coef[i],
          main_dt$comparison[i],
          label_flag(main_dt$permuted[i])
        ))
      }
    }
    
    # Interaction effects
    int_dt <- terms_dt[type == "interaction"]
    if (nrow(int_dt) > 0) {
      cat("  - Interaction effects:\n")
      for (i in seq_len(nrow(int_dt))) {
        cat(sprintf(
          fmt,
          int_dt$coef[i],
          int_dt$comparison[i],
          label_flag(int_dt$permuted[i])
        ))
      }
    }
    
    cat("\n")
  }
  
  # Return
  return(dt_res[coef != "(Intercept)"][])
}


formula <- ~ genotype * cell_type * tissue + sample

meta <- factor_encoding(
  formula = formula,
  data = as.data.table(colData(JB22_pb_sce)),
  ref_levels = c(
    genotype = "NTC",
    cell_type = "Bcells"
  ),
  verbose = TRUE
)

meta_perm <- conditional_permutation(
  formula = formula,
  data = meta,
  permute = "genotype",
  within = c("cell_type", "tissue", "sample"),
  verbose = TRUE
)

mm <- build_model(
  formula = formula,
  data = meta_perm,
  verbose = TRUE
)

exp_coef <- explain_coefs(
  formula = formula,
  data = meta_perm,
  permute = "genotype",
  verbose = TRUE
)


