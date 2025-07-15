#' Functions to convert Tilia format into list of data and metadata and further formats for different analyses
#'
#'
#' @param tilia_file data frame in Tilia-format spreadsheet containing samples as columns and taxa as rows
#' @param tilia_t output object from tilia_trans() containing list of two data frames with data and metadata
#' @param age_min minimum age for the selection
#' @param age_max maximum age for the selection
#' @param mean.age list of mean ages of time slices for
#'
#' @return Returns list of two data frames containing data and metadata.
#'
#' @examples
#' #Transform data from Tilia format
#' tilia_trans(Tilia_Prasilske)
#' @import remotes
#' @import tidyr
#' @import janitor
#' @import vegan
#' @import paleotree
#' @import changepoint
#' @import stringr
#' @import dplyr
#' @import tibble
#' @name tilia_functions
NULL

#' @describeIn tilia_functions Convert Tilia format into list data frames with data and metadata
#' @export


tilia_trans <- function(tilia_file) {
  colnames(tilia_file)[1:7] <- tilia_file[1,1:7]
  tilia_file <- tilia_file[-1,]

  # Extract metadata
  metadata <- tilia_file[str_detect(tilia_file$Code, "#") | str_detect(tilia_file$Group, "LABO"),]
  metadata[which(metadata$Code == "Lyc.spik"),"Code"] <- metadata %>% dplyr::filter(Code == "Lyc.spik") %>%
    unite(Code, Code, Element, Units, sep = "_", remove = FALSE) %>% select(Code)
  metadata[which(metadata$Code == "Lyc.tab"),"Code"] <- metadata %>% dplyr::filter(Code == "Lyc.tab") %>%
    unite(Code, Code, Element, Units, sep = "_", remove = FALSE) %>% select(Code)
  metadata <- t(metadata[,-c(2:7)]) %>% row_to_names(1)
  metadata <- rownames_to_column(as.data.frame(metadata), var = "depth")
  metadata <- metadata %>%
    mutate_if(str_detect(colnames(.), "Samp.Keywords|#Samp.Name|#Samp.Analyst|#Anal.Notes", negate = TRUE), as.numeric)

  # pivot to long
  pollen_data <- tilia_file[-which(str_detect(tilia_file$Code, "#")),] %>%
    mutate(across(8:last_col(), ~as.numeric(.))) %>%
    pivot_longer(cols = -1:-7, names_to = 'depth', values_to = 'value', names_transform = list(value = as.numeric, depth = as.numeric))

  # output variables
  data_metadata <- list(data = pollen_data, metadata = metadata)
  return(data_metadata)
}

#' @describeIn tilia_functions Select data based on age
#' @export


tilia_select <- function(tilia_t, agemin, agemax) {
  samples <- which(tilia_t$metadata$`#Chron1`>agemin & tilia_t$metadata$`#Chron1`<agemax)
  tilia_t$metadata <- tilia_t$metadata[samples,]
  tilia_t$data <- tilia_t$data %>% filter(depth %in% tilia_t$metadata$depth[samples])
  return(tilia_t)
}

#' @describeIn tilia_functions Select and aggregate samples based on time slices for REVEALS analysis in LRA package
#' @export


tilia_REVEALS <- function(tilia_t, mean.age = list(), SD.age = list(), taxalist) {
  tilia_R <- data.frame(taxa = sort(unlist(taxalist)))
  colnames(tilia_R) <- "Name"
  for (i in 1:length(mean.age)){
    samples <- which(tilia_t$metadata$`#Chron1`>mean.age[[i]]-SD.age[[i]] & tilia_t$metadata$`#Chron1`<mean.age[[i]]+SD.age[[i]])
    tilia_R_sel <- tilia_t$data %>% dplyr::filter(depth %in% tilia_t$metadata$depth[samples] & Name %in% taxalist & Element == "pollen") %>%
      dplyr::group_by(Name) %>%
      dplyr::mutate(value = sum(value, na.rm = TRUE)) %>%
      dplyr::distinct(Name, value)
    colnames(tilia_R_sel)[2] <- as.character(mean.age)[i]
    tilia_R <- tilia_R %>% left_join(tilia_R_sel)
  }
  rownames(tilia_R) <- tilia_R$Name
  tilia_R <- tilia_R[,-1]
  return(tilia_R)
}

#' @describeIn tilia_functions Function to calculate percentages for selected groups
#' @export

tilia_percent <- function(tilia_t, groups_sum, groups, elements) {
  # filter data to plot
  pollen_data_f <- tilia_t$data %>%
    # filter(Element %in% elements & Group %in% groups) %>%
    mutate(value = as.numeric(value)) %>%
    mutate(value = replace_na(value, 0))

  # calculate sums
  pollen_data_sum <- pollen_data_f %>%
    filter(Element %in% elements & Group %in% groups_sum) %>%
    dplyr::group_by(depth) %>%
    dplyr::mutate(pollensum = sum(value, na.rm = TRUE)) %>%
    dplyr::distinct(depth, pollensum) %>%
    arrange(depth)

  # convert to proportions
  pollen_data_perc <- pollen_data_f %>%
    filter(Element %in% elements & Group %in% groups) %>%
    pivot_wider(id_cols = c(depth), names_from = Name, values_from = value, values_fill = 0) %>%
    arrange(depth)
  pollen_data_perc[,-1] <- pollen_data_perc[,-1]*100/pollen_data_sum$pollensum

  return(pollen_data_perc)
}

#' @describeIn tilia_functions Function to calculate pollen acummulation rates
#' @export

tilia_PAR <- function(tilia_t, groups = list("TRSH", "UPHE"), elements = list("pollen", "spores"), piwot.wider = TRUE) {
  # filter data to plot
  pollen_data_f <- tilia_t$data %>%
    mutate(value = as.numeric(value)) %>%
    mutate(value = replace_na(value, 0))

  # add concentration
  pollen_data_f <- pollen_data_f %>%
    dplyr::group_by(depth) %>%
    dplyr::mutate(conccoef = if_else(Element == "quantity added", value, NA_real_)) %>%
    tidyr::fill(conccoef, .direction = "updown") %>%
    dplyr::mutate(conccoef = if_else(Element == "counted", conccoef/value, NA_real_)) %>%
    tidyr::fill(conccoef, .direction = "updown") %>%
    dplyr::mutate(conccoef = if_else(Element == "volume", conccoef/value, NA_real_)) %>%
    tidyr::fill(conccoef, .direction = "updown") %>%
    dplyr::group_by(Name) %>%
    dplyr::mutate(conc = value * conccoef) %>%
    dplyr::arrange(depth)

  sample_deposition <- tilia_t$metadata %>%
    mutate(SDT = (lead(`#Chron1`) - `#Chron1`)/(lead(depth)-depth) * `samp.quant`^(1/3)) %>%
    fill(SDT, .direction = "down")

  pollen_data_f <- pollen_data_f %>%
    left_join(sample_deposition, by = join_by(depth)) %>%
    filter(SDT > 0) %>%
    mutate(PAR = conc / SDT) %>%
    filter(Element %in% elements & Group %in% groups) %>%
    select(Name, Element, Group, depth, PAR) %>%
    rename(value = PAR)

  if(piwot.wider == TRUE) {
    pollen_data_PAR <- pollen_data_f %>%
      pivot_wider(id_cols = c(depth), names_from = Name, values_from = value, values_fill = 0) %>%
      arrange(depth)
    # pollen_data_PAR <- pollen_data_PAR[is.na(rowSums(pollen_data_PAR[,-1])),]
    return(pollen_data_PAR)
  }
  else {
    return(list(data = pollen_data_f, metadata = tilia_t$metadata))
  }
}

#' @describeIn tilia_functions Function to calculate charcoal accumulation rate
#' @export


tilia_CHAR <- function(tilia_t, groups = list("CHAR"), elements = list("area"), unit = list("mm^2"), piwot.wider = TRUE) {
  # filter data to plot
  charcoal_data_f <- tilia_t$data %>%
    mutate(value = as.numeric(value)) %>%
    mutate(value = replace_na(value, 0))

  sample_deposition <- tilia_t$metadata %>%
    mutate(SDT = (lead(`#Chron1`) - `#Chron1`)/(lead(depth)-depth)) %>%
    fill(SDT, .direction = "down") %>%
    mutate(depthmin = depth-`#Anal.Thick`/2) %>%
    mutate(depthmax = depth+`#Anal.Thick`/2)

  charcoal_data_f <- charcoal_data_f %>%
    left_join(sample_deposition, by = join_by(depth)) %>%
    filter(SDT > 0) %>%
    mutate(CHAR = value / SDT) %>%
    filter(Element %in% elements & Group %in% groups & Units %in% unit) %>%
    select(Name, Element, Group, Units, depth, depthmin, depthmax, CHAR) %>%
    rename(value = CHAR)

  if(piwot.wider == TRUE) {
    charcoal_data_CHAR <- charcoal_data_f %>%
      pivot_wider(id_cols = c(depth, depthmin, depthmax), names_from = Units, values_from = value, values_fill = 0) %>%
      arrange(depth)
    return(charcoal_data_CHAR)
  }
  else {
    return(list(data = charcoal_data_f, metadata = tilia_t$metadata))
  }
}

#' @describeIn tilia_functions Function to extract ecological groups
#' @export

tilia_ecog <- function(tilia_t, groups, elements) {
  tilia_t$data %>% ungroup() %>% dplyr::distinct(Name, Group, Element) %>%
    filter(Group %in% groups & Element %in% elements) %>%
    select(Name, Group)
}

#' @describeIn tilia_functions Function to define taxa to plot - order taxa based on sum of values
#' @export

tilia_taxa <- function(tilia_wide) {
  tilia_wide %>% pivot_longer(cols = -1, names_to = "Name", values_to = "Values") %>%
    group_by(Name) %>%
    summarise(values = sum(Values)) %>%
    dplyr::arrange(desc(values))
}

#' @describeIn tilia_functions Harmonize pollen taxonomy with EPD accepted taxa
#' @export


tilia_EPD <- function(tilia_t, epd_p_vars) {
  # Join EPD taxonomy to counts
  p_counts_merge <- left_join(tilia_t$data %>% dplyr::filter(Element == "pollen" | Element == "spore"), epd_p_vars, by = "Name")
  p_counts_missing <- p_counts_merge[which(is.na(p_counts_merge$VarCode)),]
  p_counts_missing <- data.frame(unique(p_counts_missing$Name))
  if(length(rownames(p_counts_missing))>0) {p_counts_merge <- p_counts_merge[-which(is.na(p_counts_merge$Var.)),]}

  ### GET COUNTS OF ACCEPTED EPD NAMES
  p_counts_accept <- p_counts_merge %>%
    mutate(value = as.numeric(value)) %>%
    mutate(value = replace_na(value, 0)) %>%
    group_by(depth, AccVar.) %>%
    summarise(value = sum(value)) %>%
    rename(Var. = AccVar.)

  p_counts_accept <- left_join(p_counts_accept, epd_p_vars, by = "Var.") %>%
    rename(Group = group) %>%
    mutate(Element = if_else(Group %in% c("VACR", "BRYO"), "spore", "pollen")) %>%
    select("Name", "Group", "Element", "depth", "value") %>% ungroup()

  tilia_EPD <- list(data = p_counts_accept, metadata = tilia_t$metadata, missing_taxa = p_counts_missing)
  return(tilia_EPD)
}

#' @describeIn tilia_functions Calculate diversity indexes
#' @export


tilia_divers <- function(tilia_t, PAR = NULL, index = "rarefy", lowsum = 300, groups = list("TRSH", "UPHE"), elements = list("pollen", "spore")) {
  table_wide <- tilia_t$data %>%
    filter(Group %in% groups & Element %in% elements) %>%
    mutate(value = as.numeric(value)) %>%
    mutate(value = replace_na(value, 0)) %>%
    pivot_wider(id_cols = depth, names_from = Name, values_from = value, values_fill = 0) %>%
    arrange(depth)

  table_wide <- table_wide[rowSums(table_wide[,-1])>=lowsum,]


  if(!is.null(PAR)) {
    table_wide_PAR <- PAR$data %>%
      filter(Group %in% groups & Element %in% elements) %>%
      filter(depth %in% table_wide$depth) %>%
      mutate(value = as.numeric(value)) %>%
      mutate(value = replace_na(value, 0)) %>%
      pivot_wider(id_cols = depth, names_from = Name, values_from = value) %>%
      arrange(depth)
    delsamp <- rowSums(table_wide_PAR[,-1])!=Inf
    table_wide_PAR <- table_wide_PAR[delsamp,]
    table_wide <- table_wide %>% filter(depth %in% table_wide_PAR$depth)
    surface <- rowSums(table_wide[,-1])/rowSums(table_wide_PAR[,-1])
    sum_surf <- round(rowSums(table_wide_PAR[,-1])*min(surface))

    table_out <- table_wide %>%
      left_join(tilia_t$metadata) %>%
      select(depth, `#Chron1`)

    table_out <- table_out %>%
      ungroup() %>%
      mutate(AccRichness = 0)
    for (n in 1:length(rownames(table_wide_PAR))) {
      table_out$AccRichness[n] <- rarefy(round(table_wide[n,-1]), sample = sum_surf[n])
    }
    table_out$Surface <- surface

  } else {
    table_out <- table_wide %>%
      left_join(tilia_t$metadata) %>%
      select(depth, `#Chron1`)
  }

  if("rarefy" %in% index) {
    table_out <- table_out %>%
      ungroup() %>%
      mutate(Richness = rarefy(table_wide[,-1] %>%
                                 mutate(across(1:last_col(), ceiling)) %>%
                                 mutate(across(1:last_col(), ~as.integer(.)))
                               , sample = lowsum))
  }
  if("simpson" %in% index) {
    table_out <- table_out %>%
      ungroup() %>%
      mutate(Simpson = diversity(table_wide[,-1], index = "simpson"))
  }
  if("invsimpson" %in% index) {
    table_out <- table_out %>%
      ungroup() %>%
      mutate(invSimpson = diversity(table_wide[,-1], index = "invsimpson"))
  }
  if("PIE" %in% index) {
    table_out <- table_out %>%
      ungroup() %>%
      mutate(PIE = HurlbertPIE(table_wide[,-1]))
  }

  return(table_out)
}

#' @describeIn tilia_functions Check for the minimum sums for rarefaction
#' @export


tilia_check_rarefy <- function(tilia_t, PAR = null, groups = list("TRSH", "UPHE"), elements = list("pollen", "spore")) {
  tilia_sum <- tilia_t$data %>%
    filter(Group %in% groups & Element %in% elements) %>%
    mutate(value = as.numeric(value)) %>%
    mutate(value = replace_na(value, 0)) %>%
    select(depth, value) %>%
    group_by(depth) %>%
    summarise(total = sum(value)) %>%
    arrange(depth)
  return(tilia_sum)
}

#' @describeIn tilia_functions Calculate changepoint analysis
#' @export


tilia_parchp <- function(tilia_PAR_wide, taxon, method = "PELT", penalty = "BIC", minseglen,  all.samples = F) {
  stab.out <- data.frame(Taxon = character(), CPTS = numeric(), stability = numeric())
  for(tx in taxon) {
    CPA <- cpt.meanvar(rev(tilia_PAR_wide %>% select(all_of(tx)) %>% unlist()), method = method, penalty = penalty, minseglen = minseglen)
    stability <- param.est(CPA)$mean/sqrt(param.est(CPA)$variance)
    stability[is.na(stability)] <- 0
    stab.out <- rbind(stab.out, data.frame(Taxon = tx, CPTS = CPA@cpts, stability = stability, mean = param.est(CPA)$mean, var = param.est(CPA)$variance))
  }
  stab.out <- stab.out %>%
    left_join(tilia_PAR_wide %>%
                select(depth) %>%
                arrange(desc(depth)) %>%
                rownames_to_column("CPTS") %>%
                mutate(CPTS = as.numeric(CPTS)), by = "CPTS")
  if(all.samples == T) {
    stab.allsamples <- tilia_PAR_wide %>%
      select(depth) %>% arrange(depth) %>%
      left_join(stab.out %>% pivot_wider(id_cols = depth, names_from = Taxon, values_from = stability)) %>% fill(unlist(taxon))
    stab.allsamples.stand <- stab.allsamples
    stab.allsamples.stand[,-1] <- decostand(stab.allsamples[,-1], method = "max", margin = 2)
    return(list(Zones = stab.out, Samples = stab.allsamples, SamplesStand = stab.allsamples.stand))
  } else {
    return(stab.out)
  }
}

#' @describeIn tilia_functions Plot changepoint analysis
#' @export


tilia_plot_CPA <- function(tilia_t, taxon, method, penalty, minseglen) {
  rp1 <- riojaPlot(tilia_percent(tilia_t, list("TRSH", "UPHE"), list("TRSH", "UPHE"), list("pollen"))[,-1],tilia_t$metadata, groups = tilia_ecog(tilia_t, list("TRSH", "UPHE"), list("pollen")),
                   plot.sec.axis = F, yvar.name = "depth",
                   # sec.yvar.name = "agebp",
                   selVars = taxon, plot.cumul = TRUE,
                   # ytks1 = seq(1700, 2000, by = 10),
                   scale.percent = TRUE, xRight = 0.25)

  section <- 1
  shift <- 0.25
  for (tx in taxon) {
    assign(paste0("rp",section+1), riojaPlot(tilia_PAR(tilia_t, list("TRSH", "UPHE"), list("pollen"), piwot.wider = T) %>% select(all_of(tx)),
                                             tilia_PAR(tilia_t, list("TRSH", "UPHE"), list("pollen"), piwot.wider = T) %>% select(depth),
                                             yvar.name = "depth", scale.percent = FALSE, xRight = shift + 0.1, riojaPlot = get(paste0("rp",section))
    ))

    assign(paste0("rp",section+2), riojaPlot(tilia_parchp(tilia_PAR(tilia_t, list("TRSH", "UPHE"), list("pollen"), piwot.wider = T), taxon = tx, method = method, penalty = penalty, minseglen = minseglen) %>% filter(Taxon == tx) %>% select(stability),
                                             tilia_parchp(tilia_PAR(tilia_t, list("TRSH", "UPHE"), list("pollen"), piwot.wider = T), taxon = tx, method = method, penalty = penalty, minseglen = minseglen) %>% filter(Taxon == tx) %>% select(depth) + 2, yvar.name = "depth",
                                             scale.percent = FALSE,
                                             plot.poly = F, plot.bar = T, plot.line = F, riojaPlot = get(paste0("rp",section+1)),
                                             lwd.bar = 10,
                                             xRight = shift + 0.15)
    )

    addRPZone(get(paste("rp",section+1, sep = "")), tilia_parchp(tilia_PAR(tilia_t, list("TRSH", "UPHE"), list("pollen"), piwot.wider = T), taxon = tx, method = method, penalty = penalty, minseglen = minseglen) %>% filter(Taxon == tx) %>% select(depth) %>% unlist(), col = "red", xLeft = shift, xRight = shift + 0.15)
    section <-  section + 2
    shift <- shift + 0.15
  }

}

#' @describeIn tilia_functions Merge different proxies for ordination
#' @export

tilia_merge <- function(tilia_CHAR, tilia_percent, tilia_divers, tilia_CPA, taxa = list("Pinus", "Picea", "Fagus", "Abies"), CHARunits = "mm^2") {
  tilia_char_pollen <- tilia_percent %>%
    left_join(tilia_CHAR, by = join_by(depth >= depthmin, depth <depthmax))
  tilia_char_pollen <- tilia_char_pollen %>% select(-depth.y, -depthmax, -depthmin) %>%
    rename(depth = depth.x) %>% group_by(depth) %>% summarise(across(everything(), mean))
  # tilia.CPAw <- pivot_wider(tilia_CPA %>% select(-CPTS), names_from = Taxon, values_from = stability) %>%
  #   arrange(depth)
  tilia_envi <- tilia_char_pollen %>% select(depth, CHARunits) %>%
    left_join(tilia_CPA) %>% fill(unlist(taxa)) %>%
    rowwise() %>% mutate(stability = mean(c_across(all_of(unlist(taxa)))[c_across(all_of(unlist(taxa))) > 0], na.rm = TRUE)) %>%
    left_join(tilia_divers) %>%
    select(-c(unlist(taxa))) %>% na.omit()

  tilia_char_pollen <- tilia_char_pollen %>% select(-one_of(CHARunits)) %>% filter(depth %in% tilia_envi$depth)

  return(list(pollen = tilia_char_pollen, env = tilia_envi))
}


#' @describeIn tilia_functions Converts Neotoma sample() data into tilia_t
#' @export

neotoma_tilia <- function(n_samples) {
  if(length(unique(n_samples$collunitid))!=1) {return ("Samples contain more than one collection units!")
  }
  div_pollen <- list()
  div_pollen$data <- n_samples %>% select(taxonid, variablename, element, units, ecologicalgroup, depth, value) %>%
    rename(c(Name = variablename, Element = element, Units = units, Group = ecologicalgroup))
  div_pollen$metadata <- n_samples %>% select(sitename, depth, age, ageolder, ageyounger) %>% unique() %>%
    rename(c('#Chron1' = age, '#Chron1.Young' = ageyounger, '#Chron1.Old' = ageolder))

  div_temp <- n_samples %>% dplyr::filter(ecologicalgroup %in% "LABO") %>%
    unite(variablename, variablename, element, units, sep = "_", remove = FALSE) %>%
    select(sitename, depth, variablename, value) %>%
    pivot_wider(names_from = variablename, values_from = value)

  div_pollen$metadata <- div_pollen$metadata %>% left_join(div_temp)
  div_pollen$data <- div_pollen$data %>% dplyr::filter(Group != "LABO")
  return(div_pollen)
}
