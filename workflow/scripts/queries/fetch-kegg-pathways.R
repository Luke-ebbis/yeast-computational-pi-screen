get_kegg_data <- function(kegg_gene_id) {
  keggset <-tibble::tibble()
  for (i in 1:length(kegg_gene_id)) {
    id <- kegg_gene_id[i]
    message(stringr::str_glue("retrieving {id}"))
    kegg_gene <- stringr::str_glue("sce:{id}")
    
    tryCatch(expr = {
      obj <- KEGGREST::keggGet(kegg_gene)
      current_obj <- obj
      paths <- current_obj[[1]]$PATHWAY
      paths_df <- paths |> data.frame() |>
        tibble::rownames_to_column() |> 
        dplyr::mutate(KEGG = id)
      keggset <- dplyr::bind_rows(keggset, paths_df) },
      error = function(e) {
        message(stringr::str_glue("fetch failed for {id} because of {e}"))
      })
    
    
  }
  return(keggset)
}


meta <- readr::read_csv(snakemake@input[[1]],
                        na = c("na"))

kegg_indentifiers <- meta |> dplyr::select(dplyr::matches("KEGG\\d")) |>
  tidyr::pivot_longer(matches("KEGG\\d")) |>
  dplyr::select(value) |> unique() |> dplyr::pull(value)

paths <- kegg_indentifiers |> na.omit() |> get_kegg_data()

paths |> readr::write_csv(snakemake@output[[1]])