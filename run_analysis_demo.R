# ==============================================================================

# --- 1. Environment Setup & Libraries ---
message("[Step 1] Loading required libraries...")

# Install packages if missing (Reproducibility safety check)
required_packages <- c("survival", "survminer", "timeROC", "dplyr", "ggplot2", "data.table")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(survival)
library(survminer)
library(timeROC)
library(dplyr)
library(ggplot2)
library(data.table)

# Set seed for reproducibility (e.g., for K-means)
set.seed(123)

# ==============================================================================
# --- 2. Configuration Parameters (User Inputs) ---
# ==============================================================================
# Change these parameters to reproduce specific figures in the manuscript
CONFIG <- list(
  # Target Selection
  gene_symbol = "ABCD4",         # Target Gene
  cancer_code = "ACC",           # TCGA Cancer Code
  data_type   = "mRNA",          # Options: mRNA, miRNA, Methylation, etc.
  
  # Analysis Settings
  survival_type = "OS",          # Options: OS, DSS, DFI, PFI
  group_method  = "median",      # Options: median, mean, optimal, kmeans, custom
  custom_cutoff = 1.5,           # Only used if group_method = "custom"
  
  # ROC Settings
  roc_years = c(1, 3, 5),        # Time points for ROC
  
  # Multivariate Covariates (Column names in clinical data)
  covariates = c("age", "gender", "stage"), 
  
  # Data Path Settings
  use_demo_data = TRUE,          # Set to FALSE to use real local TCGA files
  work_dir      = getwd()        # Working directory
)

# ==============================================================================
# --- 3. Data Retrieval & Loading Functions ---
# ==============================================================================

#' Function to load TCGA data
#' @description Reproduces the data loading logic from the Shiny App
load_tcga_data <- function(config) {
  
  if (config$use_demo_data) {
    message("[Info] Generating SYNTHETIC DEMO DATA for reproducibility check...")
    # Simulation of merged Clinical + Expression data
    n_samples <- 200
    df <- data.frame(
      sample = paste0("TCGA-SAMPLE-", 1:n_samples),
      OS.time = runif(n_samples, 10, 3000),    # Survival time in days
      OS = rbinom(n_samples, 1, 0.4),          # Event status (0/1)
      age = sample(c("Young", "Old"), n_samples, replace = TRUE),
      gender = sample(c("Male", "Female"), n_samples, replace = TRUE),
      stage = sample(c("Stage I", "Stage II", "Stage III"), n_samples, replace = TRUE)
    )
    # Simulate Gene Expression (Normal distribution)
    df[[config$gene_symbol]] <- rnorm(n_samples, mean = 5, sd = 2)
    return(df)
  }
  
  # --- REAL DATA LOADING LOGIC (From app6.R) ---
  message("[Info] Attempting to load local TCGA files...")
  setwd(config$work_dir)
  
  # 1. Load Clinical Data
  if(!file.exists("clin")) stop("Error: 'clin' file not found.")
  clin <- read.table("clin", header = TRUE, sep = "\t", check.names = FALSE)
  
  # 2. Load Molecular Data Helper Function
  read_molecular_file <- function(file_name, target_id) {
    if(!file.exists(file_name)) return(NULL)
    # Fast reading using scan/grep logic to avoid loading GBs of data
    # (Simplified for script: using fread with grep command is often better in CLI, 
    # but sticking to R logic for cross-platform compatibility)
    con <- file(file_name, open = "r")
    header <- readLines(con, n = 1)
    target_data <- NULL
    pattern <- paste0("^", target_id, "\t")
    
    while(TRUE) {
      line <- readLines(con, n = 1)
      if(length(line) == 0) break
      if(grepl(pattern, line)) {
        target_data <- unlist(strsplit(line, split = "\t"))
        break
      }
    }
    close(con)
    
    if(is.null(target_data)) return(NULL)
    
    df_mol <- as.data.frame(target_data)
    colnames(df_mol) <- config$gene_symbol # Assign Gene Name
    
    # Process header
    samples <- unlist(strsplit(header, split = "\t"))
    df_mol <- cbind(sample = samples, df_mol)
    return(df_mol)
  }
  
  # 3. Handle Specific Data Types
  mol_df <- NULL
  target_lookup <- config$gene_symbol
  
  if (config$data_type == "Methylation") {
    # Methylation Logic: Gene -> Probe -> Data
    if(!file.exists("probe_met")) stop("probe_met file missing")
    probes <- read.table("probe_met", header = F, sep = "\t", stringsAsFactors = F)
    # Find probe for gene
    matched_probe <- probes[probes$V2 == target_lookup, 1][1]
    if(is.na(matched_probe)) stop(paste("No probe found for", target_lookup))
    message(paste("Map:", target_lookup, "->", matched_probe))
    mol_df <- read_molecular_file("met.gz", matched_probe)
    
  } else if (config$data_type == "miRNA") {
    # miRNA Logic: replace - with .
    target_lookup <- gsub("-", ".", target_lookup)
    mol_df <- read_molecular_file("mi_exp.gz", target_lookup)
    
  } else {
    # Default: mRNA, lncRNA, etc.
    mol_df <- read_molecular_file("exp.gz", target_lookup)
  }
  
  if(is.null(mol_df)) stop("Could not retrieve molecular data from files.")
  
  # 4. Merge
  # Filter clinical for specific cancer
  clin_sub <- clin[clin$`cancer type abbreviation` == config$cancer_code, ]
  
  # Merge logic
  final_df <- merge(clin_sub, mol_df, by.x = "sampleID", by.y = "sample") 
  
  return(final_df)
}

# ==============================================================================
# --- 4. Preprocessing Functions ---
# ==============================================================================

preprocess_data <- function(data, config) {
  message("[Step 2] Preprocessing data...")
  
  time_col <- paste0(config$survival_type, ".time")
  status_col <- config$survival_type
  gene <- config$gene_symbol
  
  # Ensure numeric
  data[[time_col]] <- as.numeric(data[[time_col]]) / 30 # Convert days to months
  data[[status_col]] <- as.numeric(data[[status_col]])
  data[[gene]] <- as.numeric(data[[gene]])
  
  # Remove NA
  data <- data[complete.cases(data[, c(time_col, status_col, gene)]), ]
  
  # Grouping (High/Low)
  cutoff <- NA
  if (config$group_method == "median") {
    cutoff <- median(data[[gene]])
  } else if (config$group_method == "mean") {
    cutoff <- mean(data[[gene]])
  } else if (config$group_method == "optimal") {
    # Surv_cutpoint
    cut.obj <- surv_cutpoint(data, time = time_col, event = status_col, variables = gene)
    cutoff <- cut.obj$cutpoint$cutpoint
  } else if (config$group_method == "kmeans") {
    km <- kmeans(data[[gene]], 2)
    centers <- tapply(data[[gene]], km$cluster, mean)
    cutoff <- mean(centers) # Approximate cutoff between clusters
  }
  
  data$group <- ifelse(data[[gene]] > cutoff, "High", "Low")
  data$group <- factor(data$group, levels = c("Low", "High"))
  
  message(paste0("   Cutoff value used: ", round(cutoff, 3), " (Method: ", config$group_method, ")"))
  message(paste0("   High group: n=", sum(data$group=="High")))
  message(paste0("   Low group:  n=", sum(data$group=="Low")))
  
  return(list(data = data, cutoff = cutoff))
}

# ==============================================================================
# --- 5. Analysis Execution ---
# ==============================================================================

# A. Run Workflow
data_raw <- load_tcga_data(CONFIG)
prep_res <- preprocess_data(data_raw, CONFIG)
df_clean <- prep_res$data

# Define Survival Object
surv_obj <- Surv(time = df_clean[[paste0(CONFIG$survival_type, ".time")]], 
                 event = df_clean[[CONFIG$survival_type]])

# --- Analysis 1: Kaplan-Meier Plot ---
message("[Step 3] Running Kaplan-Meier Analysis...")
fit_km <- survfit(surv_obj ~ group, data = df_clean)

# Plotting with App3 Style (Yellow/Blue)
p_km <- ggsurvplot(
  fit_km, 
  data = df_clean,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  palette = c("#2E9FDF", "#E7B800"), # Low=Blue, High=Yellow
  title = paste0("Survival Analysis: ", CONFIG$gene_symbol, " (", CONFIG$cancer_code, ")"),
  ggtheme = theme_bw(),
  surv.median.line = "hv"
)

# --- Analysis 2: Multivariate Cox Regression ---
message("[Step 4] Running Multivariate Cox Regression...")
# Construct formula: Surv ~ group + age + gender ...
if(length(CONFIG$covariates) > 0 & all(CONFIG$covariates %in% colnames(df_clean))) {
  f_multi <- as.formula(paste("surv_obj ~ group +", paste(CONFIG$covariates, collapse = "+")))
  fit_cox <- coxph(f_multi, data = df_clean)
  
  # Forest Plot
  p_forest <- ggforest(fit_cox, data = df_clean, 
                       main = "Multivariate Cox Hazard Ratio")
  
  # Save Cox Results Table
  cox_res <- summary(fit_cox)$coefficients
  write.csv(cox_res, file = "output_multivariate_results.csv")
} else {
  message("   [Warning] Covariates not found in data or demo mode. Skipping Multivariate.")
  p_forest <- NULL
}

# --- Analysis 3: Time-dependent ROC ---
message("[Step 5] Running Time-dependent ROC Analysis...")
# Convert years to months (since we converted time to months in preprocess)
times_months <- CONFIG$roc_years * 12 

roc_res <- timeROC(
  T = df_clean[[paste0(CONFIG$survival_type, ".time")]],
  delta = df_clean[[CONFIG$survival_type]],
  marker = df_clean[[CONFIG$gene_symbol]],
  cause = 1,
  times = times_months,
  iid = TRUE
)

# Plot ROC
roc_df <- data.frame()
for(i in 1:length(times_months)) {
  roc_df <- rbind(roc_df, data.frame(
    FP = roc_res$FP[, i], 
    TP = roc_res$TP[, i], 
    Time = paste0(CONFIG$roc_years[i], "-Year (AUC=", round(roc_res$AUC[i], 3), ")")
  ))
}

p_roc <- ggplot(roc_df, aes(x = FP, y = TP, color = Time)) +
  geom_line(size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
  theme_bw() +
  labs(title = "Time-dependent ROC", x = "1 - Specificity", y = "Sensitivity")

# ==============================================================================
# --- 6. Save Outputs ---
# ==============================================================================
message("[Step 6] Saving outputs to PDF...")
pdf("Analysis_Report_Reproducible.pdf", width = 14, height = 8)

# Page 1: KM Plot
print(p_km, newpage = FALSE)

# Page 2: Forest Plot (if exists)
if(!is.null(p_forest)) print(p_forest)

# Page 3: ROC Plot
print(p_roc)

dev.off()

message("Analysis Pipeline Completed Successfully.")
message("Results saved to: Analysis_Report_Reproducible.pdf")

# ==============================================================================
# --- 7. Session Info (For Versioning) ---
# ==============================================================================
message("\n--- Reproducibility Information ---")
sessionInfo()