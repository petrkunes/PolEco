#' Calculates the mean ecological index for the whole pollen assemblage
#'
#' The function requires a matrix of pollen counts and a matrix of pollen taxa linked with plant species associated with ecological index
#' and frequency of the species in the present-day vegetation. The function will assign plant species to each pollen grain (with permutations that can be set) and based on the species'
#' frequency in the vegetation (FreqVeg). Then it calculates the mean for each pollen assemblage including standard deviation.
#' See Kuneš et al. 2019 for details.
#'
#' @param pollen_counts data frame containing following columns: Entity, Sample, Varname, Count
#' \describe{
#' \item{Entity}{integer referring to site}
#' \item{Sample}{number referring  to depth or age}
#' \item{Varname}{character referring to pollen taxon name}
#' \item{Count}{number referring to pollen count}
#' }
#'
#' @param pollen_plants_eco data frame containing the following columns: Species, Varname, EcoIndex, FreqVeg
#' \describe{
#' \item{Species}{character referring to plant species name}
#' \item{Varname}{character referring to pollen taxon name}
#' \item{EcoIndex}{number referring to ecological index}
#' \item{FreqVeg}{number referring to frequency of the species in vegetation (use 1 if you do not wish to apply the weighting)}
#' }
#'
#' @param perm integer number of permutations (default perm = 100)
#'
#' @return Returns a data frame with mean EcoIndex and Standard Deviation for each sample.
#'
#' @usage
#' Eco.Index(pollen_counts, pollen_plants_eco, perm = 100)
#'
#' @examples
#'
#' #CALCULATE DF FOR THE FIRST 3 SAMPLES OF EXAMPLE DATA
#'
#' Eco.Index(PC_Prasilske[1:102,], pollen_plant_dist_2019)
#'
#' @references Kuneš, P., Abraham, V., & Herben, T. (2019). Changing disturbance-diversity relationships in temperate ecosystems over the past 12000 years. Journal of Ecology, 107(4), 1678–1688. doi: 10.1111/1365-2745.13136
#'
#' @note Citation of this programme: Kuneš, P., Abraham, V., & Herben, T. (2019). Changing disturbance-diversity relationships in temperate ecosystems over the past 12000 years. Journal of Ecology, 107(4), 1678–1688. doi: 10.1111/1365-2745.13136
#' @export



Eco.Index <- function(pollen_counts, pollen_plants_eco, perm = 100)

{
  DF <- data.frame(Entity = numeric(), Sample = numeric(), EcoIndex = numeric(), SD = numeric())
  pollen_counts <- aggregate(Count ~ Entity + Sample + Varname, data = pollen_counts, FUN = sum)
  pb <-
    txtProgressBar(min = 0,
                   max = sum(pollen_counts[, 4]),
                   style = 3)
  pbcount <- 0
  samples <- unique(pollen_counts[, 1:2])
  for (i in 1:length(rownames(samples))) {
    species <-  pollen_counts$Varname[which(pollen_counts[, 1] == samples[i, 1] & pollen_counts[, 2] == samples[i, 2])]
    freq_boot <- c()
    for (s in 1:length(species)) {
      species_sel <-  pollen_plants_eco[which(as.character(pollen_plants_eco[, 2]) == species[s]), ]
      for (j in 1:pollen_counts[which(pollen_counts$Entity == samples[i, 1] & pollen_counts$Sample == samples[i, 2] & pollen_counts$Varname == species[s]), 4]) {
        pbcount <- pbcount + 1
        if (length(rownames(species_sel)) > 0) {
          for (m in 1:perm) {
            random <- sample(1:length(rownames(species_sel)), 1, prob = species_sel[, 4])
            freq_boot <- c(freq_boot, species_sel[random, 3])
          }
          setTxtProgressBar(pb, pbcount)
        }
      }
    }
    DF <- rbind(DF, data.frame(Entity = samples[i, 1], Sample = samples[i, 2], EcoIndex = mean(freq_boot), SD = sd(freq_boot), SEM = sd(freq_boot)/sqrt(length(freq_boot))))
  }
  cat(" Done!\n")
  close(pb)

  return(DF)

}
