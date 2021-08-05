#' Data frame containing list of plant species with appropriate pollen type, disturbance frequency index and frequency of the species in the vegetation.
#'
#' @source Kuneš, P., Abraham, V., & Herben, T. (2019). Changing disturbance-diversity relationships in temperate ecosystems over the past 12000 years. Journal of Ecology, 107(4), 1678–1688. https://doi.org/10.1111/1365-2745.13136
#'
#' @format A data frame with four columns:
#' \describe{
#' \item{Species}{character referring to plant species name}
#' \item{Varname}{character referring to pollen taxon name according to Beug (2004)}
#' \item{DF}{number referring to whole-community disturbance frequency index according to Herben et al. (2016)}
#' \item{FreqVeg}{number referring to frequency of the species in vegetation}
#' }
#' @references Beug, H. J. (2004). Leitfaden der Pollenbestimmung für Mitteleuropa und angrenzende Gebiete. Verlag Dr. Friedrich Pfeil.
#' @references Herben, T., Chytrý, M., & Klimešová, J. (2016). A quest for species-level indicator values for disturbance. Journal of Vegetation Science, 27(3), 628–636. https://doi.org/10.1111/jvs.12384
#' @examples
#' data(pollen_plant_dist_2019)
#' \donttest{
#'  pollen_plant_dist_2019
#' }
"pollen_plant_dist_2019"

#' Pollen counts from Prášilské jezero in Šumava Mts., Czech Republic
#'
#' @source Carter, V.A., Chiverrell, R.C., Clear, J.L., Kuosmanen, N., Moravcová, A., Svoboda, M., Svobodová-Svitavská, H., Leeuwen, V., Van Leeuwen, J., van der Knaap, W.O., & Kuneš, P. 2018. Quantitative palynology informing conservation ecology in the Bohemian/Bavarian Forests of Central Europe. Frontiers in Plant Science 8: 1–14.
#' @source Carter, V.A., Moravcová, A., Chiverrell, R.C., Clear, J.L., Finsinger, W., Dreslerová, D., Halsall, K., & Kuneš, P. 2018. Holocene-scale fire dynamics of central European temperate spruce-beech forests. Quaternary Science Reviews 191: 15–30.
#' @format A data frame containing following columns: Entity, Sample, Varname, Count
#' \describe{
#' \item{Entity}{integer referring to site}
#' \item{Sample}{number referring  to depth or age}
#' \item{Varname}{character referring to pollen taxon name}
#' \item{Count}{number referring to pollen count}
#' }
#'
#' @examples
#' data(PC_Prasilske)
#' \donttest{
#'  PC_Prasilske
#' }
"PC_Prasilske"
