library(AlphaSimR)
library(dplyr)
library(tidyr)
library(stringr)

# Set seed for reproducibility
set.seed(123)

# Define simulation parameters
nInd <- 200               # Number of individuals
nChr <- 5                 # Number of chromosomes
chrLen <- 1               # Length of each chromosome in Morgans
nQTL <- 50                # Number of QTLs per chromosome
nTraits <- 26             # Number of traits
qtlPerTrait <- 10         # Number of QTLs controlling each trait
traitNames <- c("SL","BRC1_2","CK","Aux","ABA","PIN1_3_4_7","Aux_Bud","SUC","GA",
                "Decap_signal","CXE","SMXL6_7_8","FHY3_FAR1","SPL9_15","MAX2",
                "D14","D27","MAX3","MAX4","MAX1","LBO","CLAMT","PhyB","PIFs",
                "HB21_40_53","NCED3")

# Initialize simulation parameters
founderPop <- runMacs(nInd = nInd, nChr = nChr, segSites = nQTL, inbred = TRUE, species = "GENERIC")

# Define the simulation parameters
SP <- SimParam$new(founderPop)
SP$setTrackPed(TRUE)

# Assign additive QTL effects
qtlEffects <- lapply(1:nTraits, function(i) {
  qtls <- sample(1:(nChr * nQTL), qtlPerTrait, replace = FALSE)  # Select 10 random QTLs per trait
  effects <- runif(qtlPerTrait, min = 0, max = 2)  # Ensure additive effects sum within [0,2]
  return(data.frame(QTL = qtls, Effect = effects))
})

# Convert list to a single dataframe
qtlEffectDF <- bind_rows(qtlEffects, .id = "Trait")
qtlEffectDF$Trait <- factor(qtlEffectDF$Trait, labels = traitNames)



# Assign QTL effects in AlphaSimR
SP$addTraitA(nQtlPerChr = nQTL, mean = rep(1, nTraits), var = rep(1, nTraits), corA = diag(nTraits))


SP$setVarE(h2 = rep(0.5, nTraits))  # Assuming heritability of 0.5 for all traits

# Simulate the population
pop <- newPop(founderPop, simParam = SP)

# Get genotype data
genoMat_numeric <- pullSegSiteGeno(pop, simParam = SP)  # Keep numeric genotype matrix

for (i in 1:nTraits) {
  selectedQTLs <- as.numeric(qtlEffects[[i]]$QTL)  # Ensure indices are numeric
  selectedEffects <- as.numeric(qtlEffects[[i]]$Effect)  # Ensure effects are numeric
  
  traitValues[, i] <- rowSums(genoMat_numeric[, selectedQTLs] * selectedEffects)
}


# Convert genotypes to "AA", "AC", or "CC"
genoMat[genoMat == 0] <- "AA"
genoMat[genoMat == 1] <- "AC"
genoMat[genoMat == 2] <- "CC"

# Convert genotype matrix to dataframe
genoDF <- as.data.frame(genoMat)
colnames(genoDF) <- paste0("QTL", 1:ncol(genoDF))
genoDF <- cbind(Individual = paste0("Ind", 1:nInd), genoDF)

# Compute trait values
traitValues <- matrix(0, nrow = nInd, ncol = nTraits)
colnames(traitValues) <- traitNames

for (i in 1:nTraits) {
  selectedQTLs <- qtlEffects[[i]]$QTL
  selectedEffects <- qtlEffects[[i]]$Effect
  traitValues[, i] <- rowSums(genoMat[, selectedQTLs] * selectedEffects)
}

# Ensure trait values are within range [0,2]
traitValues <- pmin(pmax(traitValues, 0), 2)

# Convert to dataframe
traitDF <- as.data.frame(traitValues)
traitDF <- cbind(Individual = paste0("Ind", 1:nInd), traitDF)

# Save to CSV
write.csv(genoDF, "data/genotypes.csv", row.names = FALSE)
write.csv(traitDF, "data/phenotypes.csv", row.names = FALSE)

cat("Genotype and phenotype data saved in the 'data/' folder.\n")

