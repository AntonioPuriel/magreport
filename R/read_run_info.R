#' Read the run information of a mag-pipeline run
#'
#' Reads `run_info.txt`, which every mag-pipeline run writes in its results
#' directory: pipeline and Nextflow versions, dates, command line, reference
#' databases with their modification date, and every parameter. Runs of
#' mag-pipeline 2.1.0 or later also have `manifest.json`, with the same
#' information ([read_manifest()]); this function is for older runs.
#'
#' @param path Path to `run_info.txt`.
#'
#' @return A list with the same elements as [read_manifest()]: `pipeline`,
#'   `run`, `databases` and `params`; `schema_version` is `NA` and `tables` is
#'   empty. Parameter values are converted to logical or numeric when they
#'   look like one.
#' @export
#'
#' @examples
#' path <- tempfile(fileext = ".txt")
#' writeLines(c(
#'   "mag-pipeline run information",
#'   "",
#'   "pipeline         : AntonioPuriel/mag-pipeline 2.0.0",
#'   "nextflow         : 26.04.6",
#'   "success          : true",
#'   "",
#'   "command line:",
#'   "  nextflow run AntonioPuriel/mag-pipeline -profile test",
#'   "",
#'   "reference databases:",
#'   "  CheckM2 : /db/uniref100.KO.1.dmnd (last modified 2021-03-23)",
#'   "  CAT     : not used",
#'   "",
#'   "parameters:",
#'   "  binners                 = concoct,metabat2",
#'   "  min_contig_length       = 1500"
#' ), path)
#' info <- read_run_info(path)
#' info$pipeline$version
#' info$databases
read_run_info <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)

  # Sections start with a line ending in ":" and no leading spaces
  section <- rep("header", length(lines))
  current <- "header"
  for (i in seq_along(lines)) {
    if (grepl("^[a-z][a-z ]*:$", lines[i])) {
      current <- sub(":$", "", lines[i])
    }
    section[i] <- current
  }
  key_value <- function(x, sep) {
    x <- x[grepl(sep, x, fixed = TRUE)]
    keys <- trimws(sub(paste0(sep, ".*$"), "", x))
    values <- trimws(substring(x, regexpr(sep, x, fixed = TRUE) + 1))
    stats::setNames(as.list(values), keys)
  }

  header <- key_value(lines[section == "header"], ":")
  pipeline <- header$pipeline
  version <- if (!is.null(pipeline) && grepl(" ", pipeline)) sub(".* ", "", pipeline) else NULL

  command <- trimws(lines[section == "command line" & grepl("^\\s+\\S", lines)])

  databases <- list()
  db_lines <- lines[section == "reference databases" & grepl("^\\s+\\S", lines)]
  db_entries <- key_value(db_lines, ":")
  for (entry in names(db_entries)) {
    value <- db_entries[[entry]]
    if (identical(value, "not used")) next
    key <- tolower(gsub("[^A-Za-z0-9]", "", entry))
    modified <- regmatches(value, regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}(?=\\)$)", value, perl = TRUE))
    databases[[key]] <- list(
      path = sub("\\s*\\(last modified [^)]*\\)$", "", value),
      modified = if (length(modified) == 1) modified else NULL
    )
  }

  params <- lapply(
    key_value(lines[section == "parameters" & grepl("^\\s+\\S", lines)], "="),
    parse_param
  )

  list(
    schema_version = NA_character_,
    pipeline = list(
      name = if (is.null(version)) pipeline else sub(" [^ ]*$", "", pipeline),
      version = version,
      revision = NULL,
      commit = NULL,
      nextflow = header$nextflow
    ),
    run = list(
      started = header$started,
      completed = header$completed,
      duration = header$duration,
      success = parse_param(header$success),
      profile = header$profile,
      container_engine = header[["container engine"]],
      command_line = if (length(command) > 0) paste(command, collapse = " ") else NULL
    ),
    databases = databases,
    params = params,
    tables = list()
  )
}

parse_param <- function(x) {
  if (is.null(x) || identical(x, "null") || identical(x, "")) {
    return(NULL)
  }
  if (x %in% c("true", "false")) {
    return(x == "true")
  }
  number <- suppressWarnings(as.numeric(x))
  if (!is.na(number) && grepl("^-?[0-9.]+(e-?[0-9]+)?$", x)) {
    return(number)
  }
  x
}

#' Describe a mag-pipeline run for the methods of a paper
#'
#' Turns the run information ([read_manifest()] or [read_run_info()]) into
#' four tables: the pipeline and the run, the reference databases with their
#' dates, the versions of the tools (manifests of mag-pipeline 2.1.0 or later)
#' and the parameters that were set.
#'
#' @param info A list returned by [read_manifest()] or [read_run_info()], or a
#'   run returned by [read_pipeline()] (its attribute `"run_info"` is used).
#'
#' @return A list of four data frames: `run` (columns `item` and `value`),
#'   `databases` (`database`, `path`, `modified`), `software` (`tool`,
#'   `version`; empty for runs without versions) and `params` (`parameter`,
#'   `value`; parameters left at `null` are not shown).
#' @export
#'
#' @examples
#' path <- tempfile(fileext = ".txt")
#' writeLines(c(
#'   "pipeline         : AntonioPuriel/mag-pipeline 2.0.0",
#'   "nextflow         : 26.04.6",
#'   "",
#'   "reference databases:",
#'   "  GTDB-Tk : /db/gtdbtk (last modified 2026-05-31)",
#'   "",
#'   "parameters:",
#'   "  min_contig_length       = 1500"
#' ), path)
#' describe_run(read_run_info(path))
describe_run <- function(info) {
  if (!is.null(attr(info, "run_info"))) {
    info <- attr(info, "run_info")
  }
  if (is.null(info) || is.null(info$pipeline)) {
    stop("No run information: the run has neither manifest.json nor run_info.txt", call. = FALSE)
  }

  text <- function(x) {
    if (is.null(x) || length(x) == 0) {
      return(NA_character_)
    }
    if (is.numeric(x)) {
      x <- format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
    }
    paste(as.character(x), collapse = ", ")
  }
  pipeline <- info$pipeline
  version <- text(pipeline$version)
  if (!is.null(pipeline$revision) || !is.null(pipeline$commit)) {
    version <- sprintf("%s (%s)", version,
                       paste(c(pipeline$revision, substr(pipeline$commit %||% "", 1, 7)),
                             collapse = ", "))
  }
  run <- data.frame(
    item = c("Pipeline", "Version", "Nextflow", "Started", "Completed",
             "Duration", "Success", "Profile", "Command line"),
    value = c(text(pipeline$name), version, text(pipeline$nextflow),
              text(info$run$started), text(info$run$completed),
              text(info$run$duration), text(info$run$success),
              text(info$run$profile), text(info$run$command_line)),
    stringsAsFactors = FALSE
  )
  run <- run[!is.na(run$value), , drop = FALSE]

  databases <- data.frame(
    database = names(info$databases) %||% character(),
    path = vapply(info$databases, function(d) text(d$path), character(1)),
    modified = vapply(info$databases, function(d) text(d$modified), character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  set <- Filter(Negate(is.null), info$params)
  params <- data.frame(
    parameter = names(set) %||% character(),
    value = vapply(set, text, character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  software <- Filter(Negate(is.null), info$software)
  software <- data.frame(
    tool = names(software) %||% character(),
    version = vapply(software, text, character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  list(run = run, databases = databases, software = software, params = params)
}
