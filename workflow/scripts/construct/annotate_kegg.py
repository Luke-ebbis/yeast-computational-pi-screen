#! /usr/bin/env python
from collections.abc import Iterable
import glob
import os
from pathlib import Path
from typing import List
from polars import dataframe
import requests
import rdflib
from rdflib import graph
import yaml
from rdflib.namespace import CSVW, DC, DCAT, DCTERMS, DOAP, FOAF, ODRL2, ORG, OWL, \
                           PROF, PROV, RDF, RDFS, SDO, SH, SKOS, SOSA, SSN, TIME, \
                           VOID, XMLNS, XSD
from rdflib import Graph, Namespace, URIRef, Literal, BNode
import polars as pl
from Bio.KEGG import Gene
from Bio.KEGG import REST
OVERWRITE = False


def read_metadata(filepath: str):
    return pl.read_csv(filepath,
                      null_values="na")


def collect_kegg_identifiers(metadata: pl.DataFrame):
   id_list = [metadata.get_column("KEGG1").to_list(),
              metadata.get_column("KEGG1").to_list()]
   return set([i for i in id_list][0])


def download_kegg_information(kegg_genes: List[str], 
                             organism_code: str = "sce",
                             folder = 'results/data/queries'):
    """Function to download all the KEGG identifiers

    Parameters
    ----------
    kegg_genes : List[str]
        The kegg genes to look up
    organism_code : str = "sce"
        The organism code to look up. Yeast by default.
    """
    for kegg_gene in kegg_genes:
        rq = f"http://rest.kegg.jp/get/{organism_code}:{kegg_gene}"
        result = requests.get(rq)
        outpath = Path(f"{folder}/{organism_code}/{kegg_gene}.txt")
        outpath.parent.mkdir(parents=True, exist_ok=True)
        write_text(result.text, output=outpath)


def write_text(text: str, output: Path, overwrite=False):
    print(f"writing to {output}")
    output.write_text(text)



def main():
   input_file = snakemake.input[0] 
   metadata = read_metadata(input_file)
   kegg_identifiers = collect_kegg_identifiers(metadata=metadata)
   print(len(kegg_identifiers))
   download_kegg_information(["YLL041C"])



if __name__ == "__main__":

    main()
