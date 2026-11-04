#' Plot a circular tree of the MAGs with annotation rings
#'
#' Draws the MAG phylogeny as a circular tree, in the style of iTOL figures,
#' with rings around it for the lineage, the novelty against GTDB, the MIMAG
#' quality, the hydrocarbon degradation markers and the response to a
#' disturbance, and bars with the abundance of each MAG.
#'
#' Every ring after the lineage is optional and is drawn only when its data
#' are given, so the same function works with just a tree and the MIMAG table.
#'
#' @param tree A tree of the MAGs: an `ape::phylo` object or the path to a
#'   Newick file, for example the one written by `gtdbtk infer`. Tips must be
#'   named as the genomes in `mags`.
#' @param mags A data frame returned by [read_mimag_quality()].
#' @param genes Optional data frame returned by [read_hydrocarbon_genes()]:
#'   adds one ring per process of [hydrocarbon_groups()].
#' @param abundance Optional data frame returned by [read_mag_abundance()]:
#'   adds the outer bars with the highest relative abundance of each MAG in
#'   any sample.
#' @param samples,condition,contrast Optional, used with `abundance` to add a
#'   ring with the response to a disturbance: `samples` is a data frame with
#'   the columns `sample` and `condition` (for example from
#'   [parse_sample_names()]), and `contrast` a named character vector of length
#'   two with the condition codes and the names shown in the figure, numerator
#'   first, for example `c(O = "H", C = "L")`. The ring shows
#'   \eqn{log_2} of the mean relative abundance in the first condition over the
#'   second, with a pseudocount of 0.01 %.
#' @param rank Taxonomic rank used for the lineage colours.
#' @param split Lineages at `rank` that are split by class instead, such as
#'   `"Pseudomonadota"`. Use `NULL` to split none.
#' @param n_lineages Number of lineages with their own colour; the rest are
#'   grouped into "Other". At most 8.
#' @param outgroup Lineage at `rank` used to root the tree, or `NULL` to keep
#'   the root of the file. The default roots on the Patescibacteriota, as in
#'   GTDB.
#' @param bar_scale `"sqrt"` (the default) spreads the many small values so
#'   their bars are visible; `"linear"` keeps the bar length proportional.
#' @param ring_labels If `TRUE`, writes the name of each ring in the opening of
#'   the tree.
#' @param interactive If `TRUE`, returns an interactive figure (a ggiraph
#'   htmlwidget) instead of a `ggplot`. Pointing at a MAG shows its name,
#'   taxonomy, quality, markers and abundance, and highlights it in every ring;
#'   pointing at a coloured cell of a hydrocarbon ring highlights every MAG of
#'   that process. A search box above the tree highlights and lists the MAGs
#'   whose name, taxonomy, lineage or CANT-HYD markers contain the text typed.
#'   Save it with `htmlwidgets::saveWidget()`. Needs ggiraph.
#'
#' @return A `ggplot` object, or a `girafe` htmlwidget when `interactive = TRUE`.
#' @export
#'
#' @details
#' Genomes of `mags` that are not in the tree (for example archaea, when the
#' tree is bacterial) are left out with a message.
#'
#' Needs the package ggtree from Bioconductor, and ape and ggnewscale from
#' CRAN:
#'
#' ```r
#' install.packages(c("BiocManager", "ape", "ggnewscale"))
#' BiocManager::install("ggtree")
#' ```
#'
#' @examples
#' \donttest{
#' if (all(vapply(c("ggtree", "ggnewscale", "ape"),
#'                requireNamespace, logical(1), quietly = TRUE))) {
#'   dir <- system.file("extdata", "marine_honour", package = "magreport")
#'   mags <- read_mimag_quality(file.path(dir, "mimag_quality.tsv"))
#'   genes <- read_hydrocarbon_genes(file.path(dir, "hydrocarbon_per_mag.tsv"))
#'   abundance <- read_mag_abundance(file.path(dir, "mag_abundance.tsv.gz"))
#'   samples <- parse_sample_names(
#'     abundance$sample,
#'     into = c("time", "site", "condition", "replicate")
#'   )
#'
#'   plot_mag_tree(
#'     file.path(dir, "bac120.tree"), mags,
#'     genes = genes, abundance = abundance,
#'     samples = samples, condition = "condition", contrast = c(O = "H", C = "L")
#'   )
#'
#'   # Interactive version, to explore in a web browser
#'   if (requireNamespace("ggiraph", quietly = TRUE)) {
#'     widget <- plot_mag_tree(
#'       file.path(dir, "bac120.tree"), mags, genes = genes,
#'       abundance = abundance, interactive = TRUE
#'     )
#'     # htmlwidgets::saveWidget(widget, "mag_tree.html")
#'   }
#' }
#' }
plot_mag_tree <- function(tree,
                          mags,
                          genes = NULL,
                          abundance = NULL,
                          samples = NULL,
                          condition = NULL,
                          contrast = NULL,
                          rank = "phylum",
                          split = "Pseudomonadota",
                          n_lineages = 8,
                          outgroup = "Patescibacteriota",
                          bar_scale = c("sqrt", "linear"),
                          ring_labels = TRUE,
                          interactive = FALSE) {
  check_ggplot2()
  for (pkg in c("ggtree", "ggnewscale", "ape")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is needed for plot_mag_tree(): see ?plot_mag_tree", call. = FALSE)
    }
  }
  if (isTRUE(interactive)) {
    for (pkg in c("ggiraph", "htmlwidgets", "htmltools")) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package '", pkg, "' is needed for interactive = TRUE", call. = FALSE)
      }
    }
  }
  bar_scale <- match.arg(bar_scale)
  stopifnot(is.data.frame(mags), all(c("bin", "quality") %in% names(mags)))
  if (!rank %in% names(mags)) {
    stop("Column '", rank, "' not found in mags: use read_mimag_quality()", call. = FALSE)
  }
  n_lineages <- min(as.integer(n_lineages), length(tree_lineage_colours))

  if (is.character(tree)) {
    tree <- ape::read.tree(tree)
  }
  stopifnot(inherits(tree, "phylo"))

  # Tips and annotation must cover the same genomes
  absent <- setdiff(mags$bin, tree$tip.label)
  if (length(absent) > 0) {
    message(length(absent), " genome(s) of mags are not in the tree and are left out: ",
            paste(utils::head(absent, 3), collapse = ", "), if (length(absent) > 3) ", ..." else "")
  }
  if (sum(tree$tip.label %in% mags$bin) < 3) {
    stop("Fewer than 3 genomes of mags are in the tree: are the tip names the genome names?",
         call. = FALSE)
  }
  extra <- setdiff(tree$tip.label, mags$bin)
  if (length(extra) > 0) {
    tree <- ape::drop.tip(tree, extra)
  }
  ann <- mags[match(tree$tip.label, mags$bin), , drop = FALSE]

  # Lineage: the rank, or the class for the lineages in `split`
  lineage <- as.character(ann[[rank]])
  if (!is.null(split) && "class" %in% names(ann)) {
    lineage <- ifelse(lineage %in% split & !is.na(ann$class), ann$class, lineage)
  }
  lineage[is.na(lineage) | lineage == ""] <- "Other"
  counts <- sort(table(lineage[lineage != "Other"]), decreasing = TRUE)
  kept <- names(counts)[seq_len(min(n_lineages, length(counts)))]
  lineage[!lineage %in% kept] <- "Other"
  lineage_levels <- c(kept, if (any(lineage == "Other")) "Other")
  ann$lineage <- factor(lineage, levels = lineage_levels)
  lineage_cols <- stats::setNames(
    c(tree_lineage_colours[seq_along(kept)], "#B9B4A8")[seq_along(lineage_levels)],
    lineage_levels
  )

  if (!is.null(outgroup)) {
    out_tips <- ann$bin[ann[[rank]] %in% outgroup]
    if (length(out_tips) >= 2) {
      tree <- ape::root(tree, node = ape::getMRCA(tree, out_tips), resolve.root = TRUE)
    } else if (length(out_tips) == 1) {
      tree <- ape::root(tree, outgroup = out_tips, resolve.root = TRUE)
    } else {
      message("No genome of ", outgroup, " in the tree: kept the root of the file")
    }
  }

  # Novelty against the reference taxonomy
  rank_value <- function(column) {
    if (column %in% names(ann)) as.character(ann[[column]]) else rep(NA_character_, nrow(ann))
  }
  has <- function(x) !is.na(x) & x != ""
  ann$novelty <- factor(
    ifelse(has(rank_value("species")), "Known species",
      ifelse(has(rank_value("genus")), "New species",
        ifelse(has(rank_value("family")), "New genus", "New family or above")
      )
    ),
    levels = names(novelty_colours)
  )

  ann$quality <- factor(ann$quality, levels = quality_levels())

  # Hydrocarbon processes and markers of each MAG
  groups <- NULL
  if (!is.null(genes)) {
    groups <- hydrocarbon_groups(genes, any = FALSE)
    groups <- groups[groups$bin %in% ann$bin, , drop = FALSE]
    canthyd <- genes[genes$database == "CANT-HYD" & genes$bin %in% ann$bin, , drop = FALSE]
    markers <- tapply(canthyd$marker, canthyd$bin, function(m) paste(sort(unique(m)), collapse = ", "))
    ann$markers <- as.character(markers[ann$bin])
  }

  # Response to the disturbance
  if (!is.null(abundance) && !is.null(samples) && !is.null(condition) && !is.null(contrast)) {
    if (length(contrast) != 2 || is.null(names(contrast)) || any(names(contrast) == "")) {
      stop("contrast must be a named vector of length 2, such as c(O = \"H\", C = \"L\")", call. = FALSE)
    }
    if (!all(c("sample", condition) %in% names(samples))) {
      stop("samples needs the columns sample and ", condition, call. = FALSE)
    }
    cond <- stats::setNames(as.character(samples[[condition]]), samples$sample)[abundance$sample]
    missing_codes <- setdiff(names(contrast), cond)
    if (length(missing_codes) > 0) {
      stop("No sample of condition ", paste(missing_codes, collapse = ", "), call. = FALSE)
    }
    mean_in <- function(code) {
      x <- abundance[!is.na(cond) & cond == code, , drop = FALSE]
      as.numeric(tapply(x$rel_abundance, x$bin, mean)[ann$bin])
    }
    ann$log2_ratio <- log2((mean_in(names(contrast)[1]) + 0.01) / (mean_in(names(contrast)[2]) + 0.01))
    ann$response <- pmax(pmin(ann$log2_ratio, 6), -6)
  }

  # Highest relative abundance in any sample
  if (!is.null(abundance)) {
    top <- as.numeric(tapply(abundance$rel_abundance, abundance$bin, max)[ann$bin])
    top[is.na(top)] <- 0
    ann$top_abundance <- top
  }

  ann$tooltip <- mag_tooltip(ann, contrast)

  tree <- ggtree::groupOTU(tree, split(ann$bin, ann$lineage), group_name = "branch_lineage")
  p <- ggtree::ggtree(tree, layout = "fan", open.angle = 30, linewidth = 0.3,
                      ggplot2::aes(colour = branch_lineage)) +
    ggplot2::scale_colour_manual(values = c(lineage_cols, "0" = "#8C887D"), guide = "none")

  # Rings are drawn as rectangles in the tree's own coordinates: x is the
  # distance from the root, y the position of each tip. In the fan layout each
  # rectangle becomes a sector of a ring.
  tips <- p$data[p$data$isTip, c("label", "y")]
  ann$y <- tips$y[match(ann$bin, tips$label)]
  depth <- max(p$data$x, na.rm = TRUE)
  ring_width <- 0.048 * depth
  gap <- 0.006 * depth
  edge <- depth + 0.06 * depth # inner edge of the next ring
  rings <- data.frame(label = character(), centre = numeric())
  key <- function(order) ggplot2::guide_legend(order = order, keyheight = 0.7, keywidth = 0.7)

  # One rectangle per MAG; with interactive = TRUE each one carries a tooltip
  # and an id: MAGs sharing an id are highlighted together
  rects <- function(df) {
    if (isTRUE(interactive)) {
      ggiraph::geom_rect_interactive(
        data = df, inherit.aes = FALSE,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = value,
                     tooltip = tooltip, data_id = data_id)
      )
    } else {
      ggplot2::geom_rect(
        data = df, inherit.aes = FALSE,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = value)
      )
    }
  }

  add_ring <- function(p, values, scale, name, width = ring_width, extra_gap = 0,
                       ids = ann$bin, tooltips = ann$tooltip) {
    edge <<- edge + extra_gap
    df <- data.frame(
      xmin = edge, xmax = edge + width,
      ymin = ann$y - 0.5, ymax = ann$y + 0.5,
      value = values, tooltip = tooltips, data_id = ids
    )
    rings <<- rbind(rings, data.frame(label = name, centre = edge + width / 2))
    edge <<- edge + width + gap
    p + rects(df) + scale + ggnewscale::new_scale_fill()
  }

  p <- add_ring(p, ann$lineage, name = "Lineage",
                ggplot2::scale_fill_manual(values = lineage_cols, name = "Lineage", guide = key(1)))
  p <- add_ring(p, ann$novelty, name = "Novelty",
                ggplot2::scale_fill_manual(values = novelty_colours, name = "Novelty (GTDB)",
                                           drop = TRUE, guide = key(2)))
  p <- add_ring(p, ann$quality, name = "MIMAG",
                ggplot2::scale_fill_manual(
                  values = quality_colours_tree, drop = TRUE,
                  labels = c(
                    high = "High: >90 % complete, <5 % contamination,\n5S, 16S and 23S rRNA, >= 18 tRNA types",
                    medium = "Medium: >= 50 % complete, <10 % contamination",
                    low = "Low: below medium"
                  ),
                  name = "MIMAG quality (Bowers et al. 2017)", guide = key(3)
                ))

  # One ring per hydrocarbon process
  if (!is.null(groups)) {
    processes <- levels(droplevels(groups$group))
    short <- c(
      "Alkanes (aerobic)" = "Alkanes", "Short-chain alkanes (aerobic)" = "Short alkanes",
      "Aromatics and PAHs (aerobic)" = "Aromatics", "Anaerobic" = "Anaerobic", "Other" = "Other markers"
    )
    for (i in seq_along(processes)) {
      present <- ann$bin %in% groups$bin[groups$group == processes[i]]
      header <- paste0("<b>", processes[i], "</b>: ", sum(present), " MAGs<br/>")
      p <- add_ring(
        p, factor(ifelse(present, "yes", "no"), levels = c("yes", "no")),
        name = short[[processes[i]]], extra_gap = if (i == 1) 2 * gap else 0,
        ids = ifelse(present, paste0("process:", processes[i]), ann$bin),
        tooltips = ifelse(present, paste0(header, ann$tooltip), ann$tooltip),
        ggplot2::scale_fill_manual(values = c(yes = "#8A2E0F", no = "#F1EFEA"), guide = "none")
      )
    }
  }

  # Response to the disturbance
  if ("response" %in% names(ann)) {
    p <- add_ring(
      p, ann$response, name = paste(contrast[[1]], "/", contrast[[2]]),
      width = 1.3 * ring_width, extra_gap = 2 * gap,
      ggplot2::scale_fill_gradient2(
        low = "#01665e", mid = "#F4F1EA", high = "#8c510a", midpoint = 0,
        limits = c(-6, 6), breaks = c(-6, 0, 6),
        labels = c("<= -6", "0", ">= 6"), na.value = "#F1EFEA",
        name = paste0("log2 ", contrast[[1]], " / ", contrast[[2]], "\nmean abundance"),
        guide = ggplot2::guide_colourbar(order = 5, barheight = 4, barwidth = 0.7)
      )
    )
  }

  n_tips <- nrow(ann)

  # Outer bars: highest relative abundance in any sample
  if (!is.null(abundance)) {
    top <- ann$top_abundance
    transform <- if (bar_scale == "sqrt") sqrt else identity
    base <- edge + 2 * gap
    per_unit <- 0.45 * depth / transform(max(top))

    ticks <- if (bar_scale == "sqrt") c(0, 1, 5, 10, 20, 50, 100) else pretty(c(0, max(top)), 4)
    ticks <- ticks[ticks <= max(top) * 1.05]
    at <- base + per_unit * transform(ticks)

    bars <- data.frame(
      xmin = base, xmax = base + per_unit * transform(top),
      ymin = ann$y - 0.4, ymax = ann$y + 0.4,
      value = paste0("Highest relative abundance\nin any sample (%",
                     if (bar_scale == "sqrt") ", square-root scale" else "", ")"),
      tooltip = ann$tooltip, data_id = ann$bin
    )
    p <- p +
      ggplot2::annotate("segment", x = at, xend = at, y = 0.5, yend = n_tips + 0.5,
                        colour = "#E6E0D3", linewidth = 0.2) +
      rects(bars) +
      ggplot2::scale_fill_manual(values = "#7A7463", name = "Outer bars", guide = key(6)) +
      # Tick labels written across the opening, like the ring names, so close
      # ticks of the square-root scale do not overlap
      ggplot2::annotate("text", x = at, y = 0, label = paste0(ticks, " %"),
                        angle = min(p$data$angle[p$data$isTip], na.rm = TRUE) + 90,
                        size = 1.9, colour = "#4A473F", hjust = 1, vjust = 0.5)
  }

  # Ring names in the opening of the fan, centred on each ring
  if (isTRUE(ring_labels)) {
    last_angle <- max(p$data$angle[p$data$isTip], na.rm = TRUE)
    p <- p + ggplot2::annotate(
      "text", x = rings$centre, y = n_tips + 2, label = rings$label,
      angle = last_angle + 90, hjust = 0, vjust = 0.5,
      size = 2.3, colour = "#4A473F"
    )
  }

  p <- p + ggplot2::theme(
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 8, face = "bold", colour = "#2B3119"),
    legend.text = ggplot2::element_text(size = 7, colour = "#2B3119"),
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )

  if (!isTRUE(interactive)) {
    return(p)
  }

  widget <- ggiraph::girafe(
    ggobj = p, width_svg = 11, height_svg = 9,
    options = list(
      ggiraph::opts_hover(css = "stroke:#1F1D1A;stroke-width:1.2px;"),
      ggiraph::opts_hover_inv(css = "opacity:0.35;"),
      ggiraph::opts_tooltip(
        css = paste(
          "background:#FDFBF7;color:#2B3119;border:1px solid #D8C8A8;",
          "border-radius:4px;padding:6px 8px;font-family:sans-serif;font-size:12px;"
        ),
        use_fill = FALSE
      ),
      ggiraph::opts_selection(type = "none"),
      ggiraph::opts_toolbar(saveaspng = TRUE)
    )
  )
  add_tree_search(widget, ann)
}

# A search box above the interactive tree: MAGs whose name, taxonomy, lineage
# or markers contain the text are highlighted in every ring and listed
add_tree_search <- function(widget, ann) {
  # Empty text for missing values (GTDB names such as "JANAQB01" contain "NA")
  taxon <- function(column) {
    x <- if (column %in% names(ann)) as.character(ann[[column]]) else rep("", nrow(ann))
    ifelse(is.na(x), "", x)
  }
  text <- tolower(paste(
    ann$bin, taxon("phylum"), taxon("class"), taxon("order"), taxon("family"),
    taxon("genus"), taxon("species"), as.character(ann$lineage), taxon("markers")
  ))
  label <- ifelse(
    taxon("species") != "", taxon("species"),
    ifelse(taxon("genus") != "", paste(taxon("genus"), "sp."),
      ifelse(taxon("family") != "", paste(taxon("family"), "(family)"), as.character(ann$lineage))
    )
  )

  box <- paste0("magreport-search-", widget$x$uid)
  search <- htmltools::div(
    style = "font-family:sans-serif;font-size:13px;margin:4px 0 8px 0;",
    htmltools::tags$input(
      id = box, type = "search",
      placeholder = "Search a MAG, taxon or marker (e.g. Alcanivorax, AlkB)",
      style = paste("width:360px;max-width:95%;padding:5px 8px;border:1px solid #D8C8A8;",
                    "border-radius:4px;")
    ),
    htmltools::tags$div(id = paste0(box, "-out"), style = "color:#4A473F;margin-top:4px;")
  )

  widget <- htmlwidgets::prependContent(widget, search)
  htmlwidgets::onRender(widget, "
    function(el, x, data) {
      var input = document.getElementById(data.box);
      var out = document.getElementById(data.box + '-out');
      if (!input) return;
      function cells() { return el.querySelectorAll('svg [data-id]'); }
      input.addEventListener('input', function() {
        var q = input.value.trim().toLowerCase();
        var hits = {};
        var listed = [];
        if (q.length >= 2) {
          for (var i = 0; i < data.bin.length; i++) {
            if (data.text[i].indexOf(q) !== -1) {
              hits[data.bin[i]] = true;
              listed.push(data.bin[i] + ' (' + data.label[i] + ')');
            }
          }
        }
        var any = listed.length > 0;
        cells().forEach(function(c) {
          var id = c.getAttribute('data-id') || '';
          var hit = hits[id] === true;
          if (!hit && id.indexOf('process:') === 0) {
            var t = c.getAttribute('title') || '';
            for (var b in hits) { if (t.indexOf(b + '&lt;') !== -1 || t.indexOf(b + '<') !== -1) { hit = true; break; } }
          }
          if (!any) { c.style.opacity = ''; c.style.stroke = ''; c.style.strokeWidth = ''; }
          else if (hit) { c.style.opacity = '1'; c.style.stroke = '#1F1D1A'; c.style.strokeWidth = '1.2px'; }
          else { c.style.opacity = '0.15'; c.style.stroke = ''; c.style.strokeWidth = ''; }
        });
        if (q.length < 2) { out.textContent = ''; }
        else if (!any) { out.textContent = 'No MAG matches \u201c' + input.value + '\u201d'; }
        else {
          out.textContent = listed.length + ' MAG' + (listed.length > 1 ? 's' : '') + ': ' +
            listed.slice(0, 25).join(', ') + (listed.length > 25 ? ', ...' : '');
        }
      });
    }",
    data = list(box = box, bin = ann$bin, text = text, label = label)
  )
}

# Text shown when pointing at a MAG in the interactive tree
mag_tooltip <- function(ann, contrast = NULL) {
  esc <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    gsub(">", "&gt;", x, fixed = TRUE)
  }
  named <- function(column) {
    if (column %in% names(ann)) as.character(ann[[column]]) else rep(NA_character_, nrow(ann))
  }
  species <- named("species")
  genus <- named("genus")
  family <- named("family")
  taxon <- ifelse(!is.na(species) & species != "", species,
    ifelse(!is.na(genus) & genus != "", paste(genus, "sp."),
      ifelse(!is.na(family) & family != "", paste(family, "(family)"), "unclassified")
    )
  )

  lines <- paste0(
    "<b>", esc(ann$bin), "</b><br/>",
    "<i>", esc(taxon), "</i> (", esc(as.character(ann$lineage)), ")<br/>",
    esc(as.character(ann$novelty)), "<br/>",
    "MIMAG ", esc(as.character(ann$quality))
  )
  if (all(c("completeness", "contamination") %in% names(ann))) {
    lines <- paste0(lines, ": ", round(as.numeric(ann$completeness), 1), " % complete, ",
                    round(as.numeric(ann$contamination), 1), " % contamination")
  }
  if ("markers" %in% names(ann)) {
    m <- ifelse(is.na(ann$markers), "none", esc(ann$markers))
    lines <- paste0(lines, "<br/>CANT-HYD: ", m)
  }
  if ("top_abundance" %in% names(ann)) {
    lines <- paste0(lines, "<br/>Highest abundance: ", signif(ann$top_abundance, 3), " %")
  }
  if ("log2_ratio" %in% names(ann) && !is.null(contrast)) {
    lines <- paste0(lines, "<br/>log2 ", esc(contrast[[1]]), " / ", esc(contrast[[2]]), ": ",
                    sprintf("%.2f", ann$log2_ratio))
  }
  lines
}

# Lineage colours: a categorical order checked for colour-blind separation
tree_lineage_colours <- c(
  "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
  "#e87ba4", "#008300", "#6250d6", "#e34948"
)

novelty_colours <- c(
  "Known species" = "#ECEAE4", "New species" = "#B8B3A6",
  "New genus" = "#7D786C", "New family or above" = "#3A3731"
)

quality_colours_tree <- c(high = "#C9A227", medium = "#E4E1D9", low = "#B05A45")

utils::globalVariables(c(
  "branch_lineage", "xmin", "xmax", "ymin", "ymax", "value", "tooltip", "data_id"
))
