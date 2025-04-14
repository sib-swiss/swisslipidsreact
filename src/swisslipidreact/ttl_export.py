# Full script that processes all rows and exports RDF to a file using rdflib
# If alter namespace inform Marco!

from rdflib import Graph, Namespace, RDF, RDFS, URIRef, Literal
import pandas as pd
import re

# Define namespaces
RH = Namespace("http://rdf.rhea-db.org/")
#EX = Namespace("http://rdf.example.org/")
SLM = Namespace("https://www.swisslipids.org/#/entity/SLM:") #fixme check with Jerven and if change announce to Marco and Sebastien
CHEBI = Namespace("http://purl.obolibrary.org/obo/CHEBI_")

ANASVES = Namespace("http://rdf.example.org/anasves/")
ANASVES_str = "http://rdf.example.org/anasves/"

GO = Namespace("http://purl.obolibrary.org/obo/GO_")
RECONX = Namespace("https://reconx.vital-it.ch/kg/")
                   
# Load full dataset from user TSV
# tsv_data = """
# MASTER_ID\treaction_smiles\tchebi_equation\tswisslipids_equation\tRInChI\tWeb-RInChIKey
# 10272\t...\t1 15378 + 1 72859 + 1 456216 = 1 82929 + 1 30616\t1 NA + 1 SLM:000000808 + 1 NA = 1 SLM:000118523 + 1 NA\t...\tWeb-RInChIKey=FNWAMYNOSKDNAKXHE-EYGIABUZYIHKWSA
# 10520\t...\t1 63774 + 1 66914 = 1 15378 + 1 NA + 1 58223\t1 SLM:000500471 + 1 NA = 1 NA + 1 SLM:000500571 + 1 NA\t...\tWeb-RInChIKey=ZHPIQPDSWAVVBBYRS-IKGLJMJDLTMWWSA
# """
# df = pd.read_csv(StringIO(tsv_data), sep="\t")
df = pd.read_csv('../../results/reaction_per_palmitic_acid.tsv', sep="\t")

df.dropna(subset=['Web-RInChIKey'], inplace=True)
print(f'{len(df)} unique RInChIKey reactions')

# Parse equation sides
def parse_equation_sides(equation):
    if pd.isna(equation):
        return [], []
    left, right = equation.split(" = ")
    parse_side = lambda s: [(int(x.split()[0]), x.split()[1]) for x in s.split(" + ")]
    return parse_side(left), parse_side(right)

# Merge compound sources, preferring SwissLipids
def merge_compounds(primary, fallback):
    result = []
    for p, f in zip(primary, fallback):
        stoich_coef = p[0] if p else f[0]
        comp_id = p[1] if p and p[1] != "NA" else f[1]
        result.append((stoich_coef, comp_id))
    return result

# Create RDF graph
g = Graph()
g.bind("rdf", RDF)
g.bind("rdfs", RDFS)
g.bind("rh", RH)
g.bind("chebi", CHEBI)
g.bind("slm", SLM)
g.bind("anasves", ANASVES)

MNetIRI = ANASVES["lipid_C14"]
# g.add((MNetIRI, RDF.type, ANASVES.MNet)) # Marco does not need it but I can keep it
g.add((MNetIRI, RDFS.label, Literal("SwissLipids ontology * Rhea for C14")))

# Process each row
for _, row in df.iterrows():
    accession = row['Web-RInChIKey'].replace("Web-RInChIKey=", "")
    reaction_uri = f"{ANASVES_str}{accession.replace('-', '_')}"
    left_uri = URIRef(f"{reaction_uri}_L")
    right_uri = URIRef(f"{reaction_uri}_R")

    # Reaction triple
    reaction_ref = URIRef(reaction_uri)
    g.add((MNetIRI, ANASVES.reac, reaction_ref))
    g.add((reaction_ref, RDFS.subClassOf, RH[str(row['MASTER_ID'])]))
    g.add((reaction_ref, ANASVES["generatedFrom"], RH[str(row['MASTER_ID'])]))
    # add g.add((reaction_ref, RDFS.subClassOf, RH[str(row['MASTER_ID'])])) for all MASTER_ID
    g.add((reaction_ref, RH.accession, Literal(accession)))
    g.add((reaction_ref, RH.side, left_uri))
    g.add((reaction_ref, RH.side, right_uri))
    g.add((left_uri, RH.curatedOrder, Literal(1)))
    g.add((right_uri, RH.curatedOrder, Literal(2)))

    # Parse and merge equations
    chebi_l, chebi_r = parse_equation_sides(row['chebi_equation'])
    slm_l, slm_r = parse_equation_sides(row['swisslipids_equation'])
    left_compounds = merge_compounds(slm_l, chebi_l)
    right_compounds = merge_compounds(slm_r, chebi_r)

    # Add compounds for each side
    def add_compounds(side_uri, compounds):
        for stoich_coef, comp_id in compounds:
            if comp_id == "NA":
                continue
            part_uri = URIRef(f"{reaction_uri}_compound_{comp_id}")

            g.add((part_uri, RH.location, GO["0005575"]))

            comp_uri = URIRef(f"{ANASVES_str}Compound_{comp_id}")

            g.add((side_uri, RH[f"contains{stoich_coef}"], part_uri))
            g.add(( RH[f"contains{stoich_coef}"], RH.coefficient, Literal(stoich_coef)))
            g.add(( RH[f"contains{stoich_coef}"], RDFS.subPropertyOf, RH.contains))
            g.add((part_uri, RH.compound, comp_uri))
            #g.add((part_uri, RDFS.subClassOf, ANASVES.ReactionParticipant))

            full_accession = comp_id if "SLM:" in comp_id else f"CHEBI:{comp_id}"
            
            if full_accession.startswith( "SLM:"):
                chem_iri = SLM[ full_accession.replace("SLM:", "" )]
            elif full_accession.startswith( "CHEBI:"):
                chem_iri = CHEBI[ full_accession.replace("CHEBI:", "" )]
            else:
                raise Exception( "Unexpected chem name: " + full_accession )
            g.add((comp_uri, RH.accession, Literal(full_accession) ))
            g.add((comp_uri, RH.chebi, chem_iri)) # correct semantically for SLM entries since SwissLipids are extesion of ChEBI

    add_compounds(left_uri, left_compounds)
    add_compounds(right_uri, right_compounds)


# Save RDF to file
output_path = "../../results/reactions_v2.ttl"
g.serialize(destination=output_path, format="ttl")
output_path
