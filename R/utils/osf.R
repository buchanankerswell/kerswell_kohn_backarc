#######################################################
## Sync data with OSF repo ca6zu                     ##
#######################################################
upload_simulation_results_to_osf <- function(data_dir = "out", conflicts = "overwrite", verbose = TRUE) {
  if (Sys.getenv("OSF_PAT") == "") {
    stop(
      "\n ======================================================================\n",
      " [Access Denied] Unauthorized Operation.\n",
      " Only the OSF repository owner can upload data to this project.\n",
      " ======================================================================\n",
      call. = FALSE
    )
  }
  if (!dir.exists(data_dir)) {
    stop(" -- Error: Local directory '", data_dir, "' not found. Nothing to upload.\n", sep = "")
  }

  node <- osf_retrieve_node("ca6zu")

  all_files <- list.files(data_dir, recursive = TRUE, full.names = TRUE)
  files_to_upload <- all_files[!grepl("\\.(tif|tiff|geoid)$", all_files, ignore.case = TRUE)]
  print(files_to_upload)

  osf_upload(node, path = files_to_upload, conflicts = conflicts, verbose = verbose)
}

download_simulation_results_from_osf <- function(data_dir = "out", conflicts = "overwrite", verbose = TRUE) {
  if (!dir.exists(data_dir)) {
    parent_dir <- dirname(data_dir)
    if (!dir.exists(parent_dir)) dir.create(parent_dir, recursive = TRUE)

    osf_retrieve_node("ca6zu") |>
      osf_ls_files() |>
      subset(name == basename(data_dir)) |>
      osf_download(path = parent_dir, recurse = TRUE, conflicts = conflicts, verbose = verbose)
  }
}
