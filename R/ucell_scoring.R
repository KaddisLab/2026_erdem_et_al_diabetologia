# UCell-based gene signature scoring
# Rank-based method with proper parallelization for large datasets

#' Calculate UCell scores for gene signatures
#'
#' @param seurat_obj Seurat object
#' @param signatures Named list of gene signatures
#' @param signature_names Character vector of signature names to calculate.
#'   Can be a single name or multiple names for batch calculation.
#'   Default: "ku_ductal"
#' @return Seurat object with UCell scores added
#' @details
#' UCell is a rank-based method that calculates per-cell "signatures scores" for
#' each defined signature. Works with raw or normalized data. Parallelized using
#' BiocParallel for efficient computation on large datasets.
#'
#' When multiple signature names are provided, all signatures are calculated
#' in a single UCell call, which is more efficient than calling separately
#' for each signature (BiocParallel overhead is paid only once).
#'
#' @examples
#' \dontrun{
#' # Single signature
#' obj <- calculate_ucell_scores(seurat_obj, signatures, "ku_ductal")
#'
#' # Multiple signatures in one call (more efficient)
#' obj <- calculate_ucell_scores(seurat_obj, signatures,
#'   c("ku_ductal", "ifng_chemokine_effector", "pan_interferon_inflammatory"))
#' }
#' @export
calculate_ucell_scores <- function(seurat_obj, signatures, signature_names = "ku_ductal") {
  # Detect cores from SLURM allocation
  slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
  if (slurm_cpus == "") slurm_cpus <- Sys.getenv("SLURM_CPUS_ON_NODE", unset = "")
  ncores <- if (slurm_cpus != "") as.integer(slurm_cpus) else parallel::detectCores()
  ncores <- min(ncores, 8)  # Cap at 8 cores

  ncells <- ncol(seurat_obj)

  # Optimize chunk.size: larger chunks = fewer parallel operations = more efficient
  # With high memory available, we can use larger chunks
  # Target: ~8-16 chunks for good parallelization without overhead
  chunk_size <- max(1000, ceiling(ncells / 16))

  # Report signatures being calculated
  n_sigs <- length(signature_names)
  if (n_sigs == 1) {
    message("Calculating UCell scores for signature: ", signature_names)
  } else {
    message("Calculating UCell scores for ", n_sigs, " signatures: ",
            paste(signature_names, collapse = ", "))
  }
  message("Total cells: ", ncells, ", chunk size: ", chunk_size)

  # Use MulticoreParam (fork-based) for better memory tracking via autometric
  # SnowParam spawns separate processes that autometric cannot track
  # MulticoreParam keeps all memory in same process tree
  hostname <- Sys.info()["nodename"]
  bp_param <- BiocParallel::MulticoreParam(workers = ncores)

  message("Using MulticoreParam (fork) with ", ncores, " workers on ", hostname)

  # Use BPPARAM for proper parallelization with optimized chunk size
  # UCell accepts multiple signatures in a single call
  seurat_obj <- UCell::AddModuleScore_UCell(
    seurat_obj,
    features = signatures[signature_names],
    BPPARAM = bp_param,
    chunk.size = chunk_size
  )

  message("UCell scores calculated")
  return(seurat_obj)
}
