#!/usr/bin/env Rscript
#===============================================================================
# Run one marked block of a script file (e.g. a single figure) without running
# the rest of the script. Tooling helper; not part of the numbered scripts.
#
# Usage (from the handbook_R folder, or `make block name=fig8` from the root):
#   Rscript _run_block.R <block>
#
# Blocks are delimited in the script files by comment lines:
#   # ==== BLOCK: name DEPS: dep1,dep2 ====
#   ... code ...
#   # ==== END BLOCK ====
# The runner finds the script that contains <block>, runs that script's
# "setup" block, then the declared DEPS (transitively, each once), then the
# block itself. Because markers are comments, running a whole script with
# source() or Rscript is unaffected.
#===============================================================================

# The runner keeps its own objects in a local environment because the setup
# blocks start with rm(list=ls()). Blocks are evaluated in the global
# environment, as when a script is run from top to bottom.
local({

  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("usage: Rscript _run_block.R <block>", call. = FALSE)
  target <- args[1]

  # Parse the blocks of a script file
  parse_blocks <- function(file) {
    src    <- readLines(file, warn = FALSE)
    starts <- grep("^#\\s*=+\\s*BLOCK:", src)
    ends   <- grep("^#\\s*=+\\s*END BLOCK", src)
    if (length(starts) != length(ends)) stop("unbalanced BLOCK markers in ", file, call. = FALSE)
    blocks <- list()
    for (i in seq_along(starts)) {
      hdr  <- src[starts[i]]
      name <- sub("^#\\s*=+\\s*BLOCK:\\s*([A-Za-z0-9_.-]+).*$", "\\1", hdr)
      deps <- character(0)
      if (grepl("DEPS:", hdr)) {
        deps <- trimws(strsplit(sub("^.*DEPS:\\s*(.*?)\\s*=+\\s*$", "\\1", hdr), ",")[[1]])
        deps <- deps[nzchar(deps)]
      }
      body <- if (ends[i] > starts[i] + 1L) src[(starts[i] + 1L):(ends[i] - 1L)] else character(0)
      blocks[[name]] <- list(deps = deps, body = body)
    }
    blocks
  }

  # Find the script that contains the block
  files <- list.files(".", pattern = "^[0-9].*\\.R$")
  found <- Filter(function(f) target %in% names(parse_blocks(f)), files)
  if (length(found) == 0) {
    all <- unlist(lapply(files, function(f) setdiff(names(parse_blocks(f)), "setup")))
    stop(sprintf("block '%s' not found. Available: %s", target, paste(all, collapse = ", ")), call. = FALSE)
  }
  script <- found[1]
  blocks <- parse_blocks(script)

  # Run setup, then dependencies, then the block
  before <- file.info(list.files("../figures", full.names = TRUE))[, "mtime", drop = FALSE]
  done <- character(0)
  run <- function(name) {
    if (name %in% done) return(invisible())
    b <- blocks[[name]]
    if (is.null(b)) stop(sprintf("block '%s' not found in %s", name, script), call. = FALSE)
    for (d in b$deps) run(d)
    message(sprintf("==> BLOCK %s (%s)", name, script))
    eval(parse(text = b$body), envir = globalenv())
    done <<- c(done, name)
  }
  if ("setup" %in% names(blocks)) run("setup")
  run(target)

  # Report figures written by the block
  after <- file.info(list.files("../figures", full.names = TRUE))[, "mtime", drop = FALSE]
  new <- rownames(after)[is.na(before[rownames(after), "mtime"]) | after$mtime > before[rownames(after), "mtime"]]
  for (f in new) message("Saved: ", normalizePath(f))

})
