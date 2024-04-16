library(igraph)
library(tidyverse)

#' Turn a protein interaction table into a igraph object.
#'
#' @param df A 
#' @param directed 
#'
#' @return an igraph::graph()
#' @export
#'
#' @examples
make_edge_list <- function(df=path_records, directed=FALSE) {
  edges <- df |> dplyr::rename(gene_1 = `gene name1`,
                               gene_2 = `gene name2`,
                               score = `PPI Score`) |> 
    dplyr::select(dplyr::matches("gene_"))#,score)
  if(any(is.na(edges))) {
    stop("the data must not contain Na values")
  }
  g <- igraph::graph_from_edgelist(as.matrix(edges), directed = directed)
  igraph::E(g)$weight <- df$`PPI Score`
  return(g)
}



meta <- readr::read_csv("resources/data/metadata.csv",
                        na = c("na"))
metabolism_interactions <- meta |> 
  dplyr::filter(stringr::str_detect(KEGGpath1, "metabolism") |
           stringr::str_detect(KEGGpath2, "metabolism"))


paths <- readr::read_csv("pathways.csv")

path_records <- metabolism_interactions |> 
  dplyr::left_join(paths |> dplyr::rename(KEGG1 = KEGG,
                            kegg_path1 = paths)) |> 
  dplyr::left_join(paths |> dplyr::rename(KEGG2 = KEGG,
                            kegg_path2 = paths)) |> 
  dplyr::filter(!is.na(`gene name1`)) |> dplyr::filter(!is.na(`gene name2`))


path_records |> dplyr::filter(stringr::str_detect(kegg_path1, "Glycolysis") | 
                         stringr::str_detect(kegg_path2, "Glycolysis")) |> 
  dplyr::select(dplyr::matches("gene name"),
         dplyr::matches("Uniprot"),
         dplyr::matches("Loc"),
         dplyr::matches("kegg_path")) 
paths |> dplyr::count(paths) |> print(n = Inf)
paths |> dplyr::filter(stringr::str_detect(paths, "Glycolysis"))






paths |> dplyr::count(paths) |> 
  dplyr::arrange(dplyr::desc(n)) |>
  dplyr::filter(n >= 2) |> 
  print(n = Inf)






  
path_records |> dplyr::filter(KEGG2 == "YOR065W")
metabolism_interactions |> dplyr::filter(KEGG2 == "YOR065W")


g <- path_records |> dplyr::filter(stringr::str_detect(kegg_path1, "Pyruvate") | 
                         stringr::str_detect(kegg_path2, "Pyruvate")) |> 
  make_edge_list()
plot(g,
     edge.width = round(igraph::E(g)$weight * 2))



g <- path_records |> dplyr::filter(stringr::str_detect(kegg_path1, "Purine") | 
                              stringr::str_detect(kegg_path2, "Purine")) |> 
  make_edge_list()
plot(g,
     edge.width = round(igraph::E(g)$weight * 2))



g <- path_records |> dplyr::filter(stringr::str_detect(kegg_path1, "Citrate") | 
                              stringr::str_detect(kegg_path2, "Citrate")) |> 
  make_edge_list()
plot(g,
     edge.width = round(igraph::E(g)$weight * 2))

g <- path_records |> dplyr::filter(stringr::str_detect(kegg_path1, "Citrate") | 
                                     stringr::str_detect(kegg_path2, "Citrate")) |> 
  make_edge_list()
plot(g,
     edge.width = round(igraph::E(g)$weight * 2))