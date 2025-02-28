# ================================================================
# R Script: Simulation of 300 Individuals with QTLs and Traits
#           ensuring each trait gets a unique QTL set
#           and final trait values scaled to [0..2].
# ================================================================

# --- 1. Setup and Directory ---
set.seed(123)  # Reproducible random seed

if (!dir.exists("data")) {
  dir.create("data")
  cat("Created folder 'data' for output files.\n")
} else {
  cat("Folder 'data' already exists. Files will be saved there.\n")
}

# --- 2. Genome and Genotype Simulation ---
chromosomes <- paste0("chr", 1:5)
n_qtls_per_chr <- 60
total_qtls <- length(chromosomes) * n_qtls_per_chr  # = 300

qtl_positions <- data.frame(
  chrom  = rep(chromosomes, each = n_qtls_per_chr),
  pos    = rep(seq_len(n_qtls_per_chr), times = length(chromosomes)),
  stringsAsFactors = FALSE
)
qtl_positions$qtl_id <- paste0(qtl_positions$chrom, "_", qtl_positions$pos)

cat("Created QTL positions for", total_qtls, "QTLs.\n")

n_ind <- 300  # number of individuals
genotype_options <- c("AA", "AC", "CC")
geno_probs <- c(0.25, 0.50, 0.25)

# Generate genotype matrix: (individuals x QTLs)
genotype_matrix <- matrix(
  sample(genotype_options, size = n_ind * total_qtls, replace = TRUE, prob = geno_probs),
  nrow = n_ind, ncol = total_qtls
)
colnames(genotype_matrix) <- qtl_positions$qtl_id
rownames(genotype_matrix) <- paste0("Ind", 1:n_ind)

cat("Simulated genotype matrix for", n_ind, "individuals and", total_qtls, "QTLs.\n")

# Save genotype matrix
genotype_df <- as.data.frame(genotype_matrix)
write.csv(genotype_df, "data/genotypes.csv", row.names = TRUE)
cat("Genotype matrix saved to 'data/genotypes.csv'.\n\n")

# --- 3. Define Traits and Assign QTLs with Random Effects ---
trait_names_raw <- c(
  "SL_n", "BRC1_2_n", "CK_n", "Aux_n", "ABA_n",
  "PIN1_3_4_7_n", "Aux_Bud_n", "SUC_n", "GA_n", "Decap_signal_n",
  "CXE_n", "SMXL6_7_8_n", "FHY3_FAR1_n", "SPL9_15_n", "MAX2_n",
  "D14_n", "D27_n", "MAX3_n", "MAX4_n", "MAX1_n",
  "LBO_n", "CLAMT_n", "PhyB_n", "PIFs_n", "HB21_40_53_n",
  "NCED3_n"
)
trait_names <- gsub("_n$", "", trait_names_raw)  # strip off "_n" at the end
n_traits <- length(trait_names)  # 26

qtl_per_trait <- 10
total_assigned_qtls <- n_traits * qtl_per_trait  # 26 * 10 = 260

if (total_assigned_qtls > total_qtls) {
  stop("Not enough QTLs available to assign to all traits.")
}

# Choose QTLs in one block (260 unique QTLs), then slice them 10 per trait
set.seed(456)  # separate seed for QTL assignment
selected_qtls <- sample(qtl_positions$qtl_id, size = total_assigned_qtls, replace = FALSE)

trait_qtl_mapping <- data.frame()
start_idx <- 1
for (i in seq_along(trait_names)) {
  # For this trait, pick a unique block of 10 QTLs
  qtls_for_trait <- selected_qtls[start_idx:(start_idx + qtl_per_trait - 1)]
  start_idx <- start_idx + qtl_per_trait
  
  # Assign random effect sizes for these 10 QTLs
  effects <- rnorm(qtl_per_trait, mean = 0, sd = 1)
  
  df <- data.frame(
    trait   = rep(trait_names[i], qtl_per_trait),
    qtl_id  = qtls_for_trait,
    effect  = effects,
    stringsAsFactors = FALSE
  )
  trait_qtl_mapping <- rbind(trait_qtl_mapping, df)
}

cat("Assigned", qtl_per_trait, "QTLs (with random effects) to each of the", n_traits, "traits.\n")
write.csv(trait_qtl_mapping, "data/trait_qtl_mapping.csv", row.names = FALSE)
cat("Trait-QTL mapping saved to 'data/trait_qtl_mapping.csv'.\n\n")

# --- 4. Convert Genotypes to Allele Counts for Additive Model ---
geno_to_count <- function(genotypes) {
  ifelse(genotypes == "AA", 2,
         ifelse(genotypes %in% c("AC", "CA"), 1,
                ifelse(genotypes == "CC", 0, NA)))
}

# --- 5. Compute Raw Additive Trait Values ---
trait_values_raw <- matrix(0, nrow = n_ind, ncol = n_traits,
                           dimnames = list(rownames(genotype_matrix), trait_names))

for (trait in trait_names) {
  # Grab the 10 QTLs/effects for this trait
  mapping <- subset(trait_qtl_mapping, trait == trait)
  
  # Start each individual's trait sum at 0
  trait_contrib <- rep(0, n_ind)
  
  # For each QTL, add effect * genotype_count
  for (j in seq_len(nrow(mapping))) {
    qtl_id  <- mapping$qtl_id[j]
    effect  <- mapping$effect[j]
    geno    <- genotype_matrix[, qtl_id]
    dosage  <- geno_to_count(geno)  # 0/1/2
    trait_contrib <- trait_contrib + (dosage * effect)
  }
  
  trait_values_raw[, trait] <- trait_contrib
  cat("Computed raw trait values for trait:", trait, "\n")
}

# --- 6. Scale Each Trait to [0, 2] ---
trait_values_scaled <- trait_values_raw

for (trait_idx in seq_along(trait_names)) {
  raw_vec <- trait_values_raw[, trait_idx]
  min_val <- min(raw_vec)
  max_val <- max(raw_vec)
  
  # Handle no-variation case
  if (max_val == min_val) {
    # If there's no variation, everyone gets 1 (the midpoint)
    trait_values_scaled[, trait_idx] <- 1
  } else {
    # Linear scaling to [0,2]
    scaled_vec <- 2 * (raw_vec - min_val) / (max_val - min_val)
    trait_values_scaled[, trait_idx] <- scaled_vec
  }
}

cat("\nAll traits scaled to [0, 2].\n")

# --- 7. Save Final Data ---
trait_values_scaled_df <- as.data.frame(trait_values_scaled)
write.csv(trait_values_scaled_df, "data/individual_trait_values_scaled.csv", row.names = TRUE)
cat("Scaled trait values saved to 'data/individual_trait_values_scaled.csv'.\n")
