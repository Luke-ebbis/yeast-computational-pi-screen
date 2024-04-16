#! /usr/bin/env python
from collections.abc import Iterable
import glob
import os
from typing import List
from polars import dataframe
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



def main():
   input_file = snakemake.input[0] 
   metadata = read_metadata(input_file)
   kegg_identifiers = collect_kegg_identifiers(metadata=metadata)
   print(len(kegg_identifiers))



if __name__ == "__main__":

    main()
