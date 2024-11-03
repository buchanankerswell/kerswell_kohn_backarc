#######################################################
## General helper functions                          ##
#######################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
with_error_handling <- function(expr, default = stop, quiet = FALSE) {
  expr_sub <- substitute(expr)
  caller <- sys.call(-1)[[1]] |> as.character()

  tryCatch(
    suppressMessages(
      suppressWarnings(
        eval(expr_sub, parent.frame())
      )
    ),
    error = function(e) {
      if (!quiet) cat(" !! ERROR in ", caller, ": ", conditionMessage(e), "\n", sep = "")
      if (is.function(default)) default(e) else default
    }
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
create_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
create_path_dir <- function(path) {
  if (!dir.exists(dirname(path))) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load_data <- function(path) {
  with_error_handling({
    cat("    Loading data: ", path, "\n", sep = "")
    load(path, envir = .GlobalEnv)
  })
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
check_data_path <- function(path) {
  if (file.exists(path)) {
    cat(" -- Found data: ", path, "\n", sep = "")
    TRUE
  } else {
    create_path_dir(path)
    FALSE
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
check_plot_path <- function(path) {
  if (file.exists(path)) {
    cat(" -- Found plot: ", path, "\n", sep = "")
    TRUE
  } else {
    create_path_dir(path)
    cat(" -> Drawing: ", path, "\n", sep = "")
    FALSE
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
compile_jsons_to_df <- function(files) {
  read_json <- function(file) {
    with_error_handling(
      {
        fromJSON(file)
      },
      default = NULL
    )
  }
  bind_rows(map(files, read_json))
}
