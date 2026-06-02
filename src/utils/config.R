#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The R package 'yaml' is required to read configs/config.yaml.")
  }
})

find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, mustWork = TRUE)
  while (TRUE) {
    if (file.exists(file.path(current, "configs", "config.yaml"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find repository root containing configs/config.yaml.")
    }
    current <- parent
  }
}

load_config <- function(config_path = NULL) {
  if (is.null(config_path)) {
    config_path <- file.path(find_repo_root(), "configs", "config.yaml")
  }
  yaml::read_yaml(config_path)
}

repo_path <- function(..., config = NULL) {
  root <- find_repo_root()
  file.path(root, ...)
}

config_path <- function(config, key, default = NULL) {
  parts <- strsplit(key, "\\.")[[1]]
  value <- config
  for (part in parts) {
    if (!is.list(value) || is.null(value[[part]])) {
      return(default)
    }
    value <- value[[part]]
  }
  value
}

resolve_config_path <- function(config, key, default = NULL) {
  value <- config_path(config, key, default)
  if (is.null(value)) {
    return(NULL)
  }
  if (grepl("^/", value)) {
    return(value)
  }
  file.path(find_repo_root(), value)
}

ensure_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  path
}

