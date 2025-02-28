# Load the AlphaSimR package
library(AlphaSimR)

# Set a seed for reproducibility
set.seed(123)

#### 1. Create the Founder Population ####

# We define a founder population of, say, 100 individuals.
nFounders <- 100

# We want 5 chromosomes. In the original description each chromosome had 1 Morgan
# and 20 causal loci (total = 100), but because we want 25 traits × 10 non‐pleiotropic QTL
# (i.e. 250 distinct QTL) we simulate 50 segregating sites per chromosome.
nChr <- 5
segSites <- rep(50, nChr)  # 50 sites per chromosome -> 5 x 50 = 250 sites in total

# Create the founder population. (newMapPop simulates haplotypes based on a recombination map;
# by default each chromosome is 1 Morgan long.)
founderPop <- newMapPop(nInd = nFounders, nChr = nChr, segSites = segSites)

#### 2. Set Up Simulation Parameters and Add Traits ####

# Initialize the simulation parameter object using the founder population.
SP <- SimParam$new(founderPop)

# We want to simulate 25 traits. For each trait we want 10 QTL in total.
# Here we assign 2 QTL per chromosome so that 2 x 5 = 10 QTL per trait.
nTraits <- 25
qtlPerChr <- 2  # per trait

# When adding multiple traits, setting pleiotropy = FALSE ensures that each trait’s QTL
# are drawn from non-overlapping sets of loci.
for(trait in 1:nTraits) {
  SP$addTraitA(nQtlPerChr = rep(qtlPerChr, nChr),
               mean = 0,     # initial intercept (will be adjusted later)
               var = 1,      # variance before scaling
               name = paste0("Trait", trait),
               pleiotropy = FALSE)
}

#### 3. Compute and Scale Breeding Values ####

# For each trait the raw breeding value is the sum over (allele copy × QTL effect).
# We now compute the raw breeding values for each founder and then “scale”
# (shift and multiply) so that for each trait the minimum is 0 and the maximum is 2.
# (This is our way to “constrain” the additive genetic values to the 0–2 range.)

# Create a matrix to store the raw breeding values.
rawBVs <- matrix(NA, nrow = nFounders, ncol = nTraits)
colnames(rawBVs) <- paste0("Trait", 1:nTraits)

for(trait in 1:nTraits) {
  rawBVs[, trait] <- getBV(founderPop, simParam = SP, trait = trait)
}

# Now scale each trait’s breeding values to [0, 2].
scaledBVs <- rawBVs  # to store the scaled values
scalingParams <- list()  # to record the scaling parameters for each trait

cat("Scaling breeding values for each trait:\n")
for(trait in 1:nTraits) {
  minVal <- min(rawBVs[, trait])
  maxVal <- max(rawBVs[, trait])
  rangeVal <- maxVal - minVal
  if(rangeVal == 0) rangeVal <- 1  # safeguard against division by zero
  scalingFactor <- 2 / rangeVal
  # New (scaled) breeding value: subtract minimum and multiply so that the range becomes 2.
  scaledBVs[, trait] <- (rawBVs[, trait] - minVal) * scalingFactor
  scalingParams[[trait]] <- list(min = minVal, max = maxVal, factor = scalingFactor)
  cat(sprintf("  Trait %d: raw range = [%.3f, %.3f] --> scaling factor = %.3f\n",
              trait, minVal, maxVal, scalingFactor))
}

#### 4. Prepare the Parents for Divergent Selection ####

# At this point the founder population (stored in 'founderPop') along with the scaled
# breeding values (in 'scaledBVs') represent our set of parents.
# For example, you might later select individuals from the top and bottom tails of a trait’s
# distribution for divergent selection.

# Here we package the information in a list.
parents <- list(population = founderPop,
                scaledBVs = scaledBVs,
                scalingParams = scalingParams)

# (Optional) Display a summary of the scaled breeding values for the first few individuals.
cat("\nScaled Breeding Values (first 6 individuals):\n")
print(round(scaledBVs[1:6, ], 3))

# The 'parents' object now contains the simulated founders with their trait values
# in the 0–2 range. You can now use these as the starting material for your divergent
# selection simulation.

#### End of Script ####
