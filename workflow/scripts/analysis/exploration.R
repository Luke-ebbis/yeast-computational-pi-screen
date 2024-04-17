library(igraph)
library(tidyverse)
library(ggnetwork)
#' Turn a protein interaction table into a igraph object.
#'
#' @param df A 
#' @param directed 
#'
#' @return an igraph::graph()
#' @export
#'
#' @examples
make_edge_list <- function(path_records, directed=FALSE) {
  edges <- path_records |> dplyr::rename(gene_1 = `gene name1`,
                               gene_2 = `gene name2`,
                               score = `PPI Score`) |> 
    dplyr::select(dplyr::matches("gene_"))#,score)
  if(any(is.na(edges))) {
    stop("the data must not contain Na values")
  } else {
    g <- igraph::graph_from_edgelist(as.matrix(edges), directed = directed)
    igraph::E(g)$weight <- path_records$`PPI Score`  
    return(g)
  }
}

make_graph_of_filtered_data <- function(path_records, colname, path_pattern){
  colname1 <- stringr::str_glue("{colname}1")
  colname2 <- stringr::str_glue("{colname}2")
  filted_paths <- path_records |> 
    dplyr::filter(!!rlang::sym(colname1) == path_pattern |
                  !!rlang::sym(colname2) == path_pattern)
  message(path_pattern)
  print(filted_paths |> dplyr::select(dplyr::matches("gene name")))
  make_edge_list(filted_paths)
}


plot_igraph <- function(graph) {
  plot(g,
       edge.width = round(igraph::E(g)$weight * 2))
}


merge_with_kegg_table <- function(meta, paths){
  path_records <- meta |> 
    dplyr::left_join(paths |> dplyr::rename(KEGG1 = KEGG,
                                            kegg_path1 = paths)) |> 
    dplyr::left_join(paths |> dplyr::rename(KEGG2 = KEGG,
                                            kegg_path2 = paths)) |> 
    dplyr::filter(!is.na(`gene name1`)) |> dplyr::filter(!is.na(`gene name2`))
  path_records
}


merge_with_go_table <- function(meta, paths){
  path_records <- meta |> 
    dplyr::left_join(paths |> dplyr::rename(Uniprot1 = id,
                                            uniprot_go1 = go_name)) |> 
    dplyr::left_join(paths |> dplyr::rename(Uniprot2 = id,
                                            uniprot_go2 = go_name)) |> 
    dplyr::filter(!is.na(`gene name1`)) |> dplyr::filter(!is.na(`gene name2`))
  path_records
}

if(!exists("snakemake")) {
  meta <- readr::read_csv("resources/data/metadata.csv",
                          na = c("na"))
  go_paths <- readr::read_csv("results/data/queries/uniprot/primary-go-annotation.csv")
  kegg_paths <- readr::read_csv("results/data/queries/kegg/kegg-pathways.csv")
  
} else {
  meta <- readr::read_csv(snakemake@input[1],
                          na = c("na"))
  kegg_paths <- readr::read_csv(snakemake@input[2])
  go_paths <- readr::read_csv(snakemake@input[3])
}

metabolism_interactions <- meta


go_paths |> dplyr::count(go_name) |> dplyr::arrange(go_name) |> print(n = Inf)

kegg_paths |> dplyr::count(paths) |> dplyr::arrange(paths) |> print(n = Inf)

path_records <- merge_with_go_table(meta, go_paths) |> merge_with_kegg_table(kegg_paths)


kegg_paths_to_plot <- c("Pyruvate metabolism",
                        "Glycolysis / Gluconeogenesis",
                        "Citrate cycle (TCA cycle)",
                        "Carbon metabolism")
go_paths_to_plot <- c("carbohydrate metabolic process",
                      "tricarboxylic acid cycle"
  
)



for (path_name in kegg_paths_to_plot) {
  path <-  make_graph_of_filtered_data(path_records, 
                                            colname = "kegg_path",
                                            path_pattern = path_name) 
  file_path_name <- stringr::str_replace_all(path_name, "\\s/", "_") |> 
    stringr::str_replace_all("\\s", "_")
  png(stringr::str_glue("results/analysis/figures/prelim/kegg/{file_path_name}.png"),
      600, 600)
  igraph::V(path)$label.cex = 1.5
  plot(path,
       cex = 5,
       edge.width = round(igraph::E(g)$weight * 2))
  title(path_name)
  dev.off()
}

for (path_name in go_paths_to_plot) {
  path <-  make_graph_of_filtered_data(path_records, 
                                       colname = "uniprot_go",
                                       path_pattern = path_name) 
  file_path_name <- stringr::str_replace_all(path_name, "\\s/", "_") |> 
    stringr::str_replace_all("\\s", "_")
  png(stringr::str_glue("results/analysis/figures/prelim/go/{file_path_name}.png"),
      600, 600)
  igraph::V(path)$label.cex = 1.5
  plot(path,
       cex = 5,
       edge.width = round(igraph::E(g)$weight * 2))
  title(path_name)
  dev.off()
}


