#! /usr/bin/env python
from collections.abc import Iterable
import glob
import os
from pathlib import Path
from typing import List
from polars import DataFrame, dataframe
from rdflib import namespace
from rdflib.query import Result
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
import time
import re
from rdflib import Namespace
OVERWRITE = False


def read_metadata(filepath: str):
    return pl.read_csv(filepath,
                      null_values="na")


def collect_uniprot_identifiers(metadata: pl.DataFrame) -> set[str]:
   id_list = [metadata.get_column("Uniprot1").to_list(),
              metadata.get_column("Uniprot2").to_list()]
   required_identifiers = set([i for i in id_list][0])
   return required_identifiers

def metabolic_query():
    query = """
        PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
        PREFIX wikibase: <http://wikiba.se/ontology#>
        PREFIX wdt: <http://www.wikidata.org/prop/direct/>
        PREFIX wd: <http://www.wikidata.org/entity/>
        PREFIX vg: <http://biohackathon.org/resource/vg#>
        PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>
        PREFIX uberon: <http://purl.obolibrary.org/obo/uo#>
        PREFIX sp: <http://spinrdf.org/sp#>
        PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
        PREFIX sio: <http://semanticscience.org/resource/>
        PREFIX sh: <http://www.w3.org/ns/shacl#>
        PREFIX sd: <http://www.w3.org/ns/sparql-service-description#>
        PREFIX schema: <http://schema.org/>
        PREFIX sachem: <http://bioinfo.uochb.cas.cz/rdf/v1.0/sachem#>
        PREFIX rh: <http://rdf.rhea-db.org/>
        PREFIX pubmed: <http://rdf.ncbi.nlm.nih.gov/pubmed/>
        PREFIX ps: <http://www.wikidata.org/prop/statement/>
        PREFIX pq: <http://www.wikidata.org/prop/qualifier/>
        PREFIX patent: <http://data.epo.org/linked-data/def/patent/>
        PREFIX p: <http://www.wikidata.org/prop/>
        PREFIX owl: <http://www.w3.org/2002/07/owl#>
        PREFIX orthodbGroup: <http://purl.orthodb.org/odbgroup/>
        PREFIX orthodb: <http://purl.orthodb.org/>
        PREFIX orth: <http://purl.org/net/orth#>
        PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
        PREFIX obo: <http://purl.obolibrary.org/obo/>
        PREFIX np: <http://nextprot.org/rdf#>
        PREFIX nextprot_cv: <http://nextprot.org/rdf/terminology/>
        PREFIX nextprot: <http://nextprot.org/rdf/entry/>
        PREFIX mnx: <https://rdf.metanetx.org/schema/>
        PREFIX mnet: <https://rdf.metanetx.org/mnet/>
        PREFIX mesh: <http://id.nlm.nih.gov/mesh/>
        PREFIX lscr: <http://purl.org/lscr#>
        PREFIX lipidmaps: <https://www.lipidmaps.org/rdf/>
        PREFIX keywords: <http://purl.uniprot.org/keywords/>
        PREFIX insdcschema: <http://ddbj.nig.ac.jp/ontologies/nucleotide/>
        PREFIX insdc: <http://identifiers.org/insdc/>
        PREFIX identifiers: <http://identifiers.org/>
        PREFIX glyconnect: <https://purl.org/glyconnect/>
        PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
        PREFIX genex: <http://purl.org/genex#>
        PREFIX foaf: <http://xmlns.com/foaf/0.1/>
        PREFIX faldo: <http://biohackathon.org/resource/faldo#>
        PREFIX eunisSpecies: <http://eunis.eea.europa.eu/rdf/species-schema.rdf#>
        PREFIX ensembltranscript: <http://rdf.ebi.ac.uk/resource/ensembl.transcript/>
        PREFIX ensemblterms: <http://rdf.ebi.ac.uk/terms/ensembl/>
        PREFIX ensemblprotein: <http://rdf.ebi.ac.uk/resource/ensembl.protein/>
        PREFIX ensemblexon: <http://rdf.ebi.ac.uk/resource/ensembl.exon/>
        PREFIX ensembl: <http://rdf.ebi.ac.uk/resource/ensembl/>
        PREFIX ec: <http://purl.uniprot.org/enzyme/>
        PREFIX dcterms: <http://purl.org/dc/terms/>
        PREFIX dc: <http://purl.org/dc/terms/>
        PREFIX chebislash: <http://purl.obolibrary.org/obo/chebi/>
        PREFIX chebihash: <http://purl.obolibrary.org/obo/chebi#>
        PREFIX cco: <http://rdf.ebi.ac.uk/terms/chembl#>
        PREFIX busco: <http://busco.ezlab.org/schema#>
        PREFIX bibo: <http://purl.org/ontology/bibo/>
        PREFIX allie: <http://allie.dbcls.jp/>
        PREFIX SWISSLIPID: <https://swisslipids.org/rdf/SLM_>
        PREFIX GO: <http://purl.obolibrary.org/obo/GO_>
        PREFIX ECO: <http://purl.obolibrary.org/obo/ECO_>
        PREFIX CHEBI: <http://purl.obolibrary.org/obo/CHEBI_>
        PREFIX up: <http://purl.uniprot.org/core/>
        PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        select distinct ?protein ?name ?biotype ?protein_classification
        WHERE {
      		?protein a up:Protein .
            service <https://sparql.uniprot.org/> {
            ?protein up:recommendedName ?rec_name .
            ?protein up:classifiedWith ?protein_classification .
            ?rec_name up:fullName ?name .
            ?protein_classification up:database <http://purl.uniprot.org/database/go> .
            # a primary metabolism
            ?protein_classification rdfs:subClassOf <http://purl.obolibrary.org/obo/GO_0044238>.
            ?protein_classification rdfs:label ?biotype .
            }
        } limit 10
    """
    return query


def make_graph(identifiers: set[str]) -> rdflib.Graph:
    g = Graph()
    protein = URIRef("http://purl.uniprot.org/core/Protein")
    for identifier in identifiers:
        if identifier:
            prot = URIRef(f"http://purl.uniprot.org/uniprot/{identifier}")
            g.add((prot, RDF.type, protein))
    return g

def result_to_df(query_result: Result) -> DataFrame:
    identifiers = []
    names = []
    golink = []
    goname = []
    for row in query_result:
        identifiers.append(str(row.protein))
        names.append(str(row.name))
        golink.append(str(row.protein_classification))
        goname.append(str(row.biotype))
    df = DataFrame({"id" :
                    [i.replace("http://purl.uniprot.org/uniprot/", "")
                     for i in identifiers],
                    "uniprot" : identifiers,
                    "name": names,
                    "go" : golink,
                    "go_name" : goname })
    return df


def main():
   input_file = snakemake.input[0] 
   output_folder = snakemake.output[0]
   metadata = read_metadata(input_file)
   identifiers = collect_uniprot_identifiers(metadata)
   proteins = make_graph(identifiers)
   qres = proteins.query(metabolic_query())
   df = result_to_df(qres)
   df.write_csv(output_folder)


if __name__ == "__main__":

    main()
