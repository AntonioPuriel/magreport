mh_path <- function(file) {
  system.file("extdata", "marine_honour", file, package = "magreport")
}

skip_tree <- function() {
  skip_if_not_installed("ggplot2")
  for (pkg in c("ggtree", "ggnewscale", "ape")) skip_if_not_installed(pkg)
}

mh_inputs <- function() {
  abundance <- read_mag_abundance(mh_path("mag_abundance.tsv.gz"))
  list(
    mags = read_mimag_quality(mh_path("mimag_quality.tsv")),
    genes = read_hydrocarbon_genes(mh_path("hydrocarbon_per_mag.tsv")),
    abundance = abundance,
    samples = parse_sample_names(abundance$sample,
                                 into = c("time", "site", "condition", "replicate"))
  )
}

test_that("plot_mag_tree draws every ring with the full data", {
  skip_tree()
  x <- mh_inputs()

  expect_message(
    p <- plot_mag_tree(mh_path("bac120.tree"), x$mags,
                       genes = x$genes, abundance = x$abundance,
                       samples = x$samples, condition = "condition",
                       contrast = c(O = "H", C = "L")),
    "1 genome"  # the archaeon is not in the bacterial tree
  )
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))

  # Lineage, novelty, MIMAG, 4 hydrocarbon processes, the response ring and the bars
  rects <- vapply(p$layers, function(l) inherits(l$geom, "GeomRect"), logical(1))
  expect_equal(sum(rects), 3 + 4 + 1 + 1)

  # The tree has the 280 bacteria, rooted
  tips <- p$data[p$data$isTip, ]
  expect_equal(nrow(tips), 280)
})

test_that("plot_mag_tree works with only the tree and the MIMAG table", {
  skip_tree()
  x <- mh_inputs()

  p <- suppressMessages(plot_mag_tree(mh_path("bac120.tree"), x$mags, bar_scale = "linear"))
  expect_no_error(ggplot2::ggplot_build(p))
  rects <- vapply(p$layers, function(l) inherits(l$geom, "GeomRect"), logical(1))
  expect_equal(sum(rects), 3)
})

test_that("plot_mag_tree computes the response and checks its inputs", {
  skip_tree()
  x <- mh_inputs()
  tree <- ape::read.tree(mh_path("bac120.tree"))

  # The response ring is log2 of mean H over mean L, capped at +-6
  p <- suppressMessages(plot_mag_tree(tree, x$mags, abundance = x$abundance,
                                      samples = x$samples, condition = "condition",
                                      contrast = c(O = "H", C = "L")))
  rects <- Filter(function(l) inherits(l$geom, "GeomRect"), p$layers)
  response <- rects[[4]]$data # lineage, novelty, MIMAG, response
  tip_y <- p$data$y[p$data$isTip & p$data$label == "T0.metabat2.156"]
  resp <- response$value[response$ymin == tip_y - 0.5]
  ab <- x$abundance[x$abundance$bin == "T0.metabat2.156", ]
  cond <- substr(sub("^[^-]+-[^-]+-", "", ab$sample), 1, 1)
  expected <- log2((mean(ab$rel_abundance[cond == "O"]) + 0.01) /
                   (mean(ab$rel_abundance[cond == "C"]) + 0.01))
  expect_equal(resp, max(min(expected, 6), -6))

  expect_error(
    suppressMessages(plot_mag_tree(tree, x$mags, abundance = x$abundance, samples = x$samples,
                                   condition = "condition", contrast = c("O", "C"))),
    "named vector"
  )
  expect_error(
    suppressMessages(plot_mag_tree(tree, x$mags, abundance = x$abundance, samples = x$samples,
                                   condition = "condition", contrast = c(O = "H", X = "L"))),
    "No sample of condition X"
  )

  renamed <- tree
  renamed$tip.label <- paste0("x", seq_along(renamed$tip.label))
  expect_error(suppressMessages(plot_mag_tree(renamed, x$mags)), "tip names")
})

test_that("plot_mag_tree(interactive = TRUE) returns a widget with tooltips and ids", {
  skip_tree()
  skip_if_not_installed("ggiraph")
  x <- mh_inputs()

  w <- suppressMessages(plot_mag_tree(mh_path("bac120.tree"), x$mags,
                                      genes = x$genes, abundance = x$abundance,
                                      samples = x$samples, condition = "condition",
                                      contrast = c(O = "H", C = "L"), interactive = TRUE))
  expect_s3_class(w, "girafe")
  svg <- w$x$html
  # Each MAG has its own id, and the cells of a process share one
  expect_match(svg, "T0.metabat2.156", fixed = TRUE)
  expect_match(svg, "process:Alkanes (aerobic)", fixed = TRUE)
  expect_match(svg, "process:Anaerobic", fixed = TRUE)

  # Tooltip text of one MAG
  tips <- magreport:::mag_tooltip(
    data.frame(bin = "a<b", species = NA, genus = "Cycloclasticus", family = "x",
               lineage = "Gammaproteobacteria", novelty = "New species", quality = "high",
               completeness = 95.123, contamination = 1.2, markers = "NdoB, NdoC",
               top_abundance = 2.1234, log2_ratio = 3.456),
    contrast = c(O = "H", C = "L")
  )
  expect_match(tips, "a&lt;b", fixed = TRUE)
  expect_match(tips, "Cycloclasticus sp.", fixed = TRUE)
  expect_match(tips, "95.1 % complete", fixed = TRUE)
  expect_match(tips, "CANT-HYD: NdoB, NdoC", fixed = TRUE)
  expect_match(tips, "log2 H / L: 3.46", fixed = TRUE)
})
