#' KEGG orthologues involved in hydrocarbon degradation and related processes
#'
#' A curated set of KEGG orthologues (KOs) grouped into the processes that
#' matter when studying oil-polluted sediments. It is meant as a starting point
#' for [plot_pathway_presence()], not as an exhaustive list: pass your own
#' `data.frame` with the same two columns to use a different set.
#'
#' @format A data frame with two columns:
#' \describe{
#'   \item{pathway}{Process the orthologue belongs to}
#'   \item{ko}{KEGG orthologue identifier}
#' }
#'
#' @details
#' Marker genes are grouped as follows.
#'
#' * *Alkanes (aerobic)*: alkane 1-monooxygenase (`alkB`), long-chain alkane
#'   monooxygenase (`ladA`) and cytochrome P450 CYP153.
#' * *Aromatics (aerobic)*: ring-hydroxylating dioxygenases and the central
#'   benzoate and catechol routes.
#' * *PAHs (aerobic)*: naphthalene and related polycyclic routes.
#' * *Hydrocarbons (anaerobic)*: fumarate addition (`bssA`, `assA`) and the
#'   central benzoyl-CoA route (`bamA`, `bcr`).
#' * *Oxidative stress*: catalase, superoxide dismutase and peroxiredoxins.
#' * *Sulfur cycling* and *Nitrogen cycling*: environmental context, since both
#'   drive anaerobic degradation in sediments.
#'
#' @source KEGG orthology (https://www.genome.jp/kegg/ko.html)
#' @export
#'
#' @examples
#' head(hydrocarbon_kos())
#' table(hydrocarbon_kos()$pathway)
hydrocarbon_kos <- function() {
  data.frame(
    pathway = c(
      rep("Alkanes (aerobic)", 5),
      rep("Aromatics (aerobic)", 7),
      rep("PAHs (aerobic)", 5),
      rep("Hydrocarbons (anaerobic)", 5),
      rep("Oxidative stress", 6),
      rep("Sulfur cycling", 5),
      rep("Nitrogen cycling", 5)
    ),
    ko = c(
      # alkB, almA/ladA, CYP153, alcohol and aldehyde dehydrogenases
      "K00496", "K10944", "K18089", "K00114", "K00128",
      # benzoate dioxygenase, catechol 1,2- and 2,3-dioxygenase, muconate route
      "K05549", "K05550", "K03381", "K00446", "K01055", "K01821", "K01856",
      # naphthalene dioxygenase and downstream steps
      "K14579", "K14580", "K14581", "K14582", "K14583",
      # bssA, assA, bamA, benzoyl-CoA reductase
      "K07540", "K19746", "K19745", "K04112", "K04113",
      # katE/katG, sodA/sodB, ahpC, prx
      "K03781", "K03782", "K04564", "K04565", "K03386", "K24119",
      # dsrA/dsrB, aprA, sat, soxB
      "K11180", "K11181", "K00394", "K00958", "K17224",
      # nifH, nirS, nirK, nosZ, narG
      "K02588", "K15864", "K00368", "K00376", "K00370"
    ),
    stringsAsFactors = FALSE
  )
}
