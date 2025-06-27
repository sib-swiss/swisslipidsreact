def export_ttl(full_scope=True, input_path=None, output_dir=None):
    from rdflib import Graph, Namespace, RDF, RDFS, URIRef, Literal
    import pandas as pd
    import os
    import re
    from datetime import datetime
    import glob
    
    # determine default
    if full_scope:
        MNetIRI_literal = "SwissLipids ontology * Rhea for curated FA list"
        MNetIRI_string = "lipid_curated_FA_list"
    else:
        MNetIRI_literal = "SwissLipids ontology * Rhea for C16"
        MNetIRI_string = "lipid_C16"
    
    # fallback input if none given
    if input_path is None:
        search_dir = output_dir if output_dir is not None else os.getcwd()
        pattern = os.path.join(search_dir, "*_enumerated_reactions.tsv")
        candidates = glob.glob(pattern)
        if not candidates:
            raise FileNotFoundError(
                f"No enumerated_reactions TSV files found in {search_dir}. "
                f"Please provide --input explicitly."
            )
        # sort candidates by modification time
        candidates.sort(key=os.path.getmtime, reverse=True)
        input_path = candidates[0]
        print(f"Auto-detected latest enumerated reactions file: {input_path}")

    # fallback output if none given
    if output_dir is None:
        output_dir = os.getcwd()
    else:
        os.makedirs(output_dir, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = os.path.join(
        output_dir,
        f"reactions_{'curated_FA' if full_scope else 'PAL_C16'}_{timestamp}.ttl"
    )

    # load dataframe
    df = pd.read_csv(input_path, sep="\t")
    df.dropna(subset=['Web-RInChIKey'], inplace=True)

    # namespaces
    RH = Namespace("http://rdf.rhea-db.org/")
    SLM = Namespace("https://www.swisslipids.org/#/entity/SLM:")
    CHEBI = Namespace("http://purl.obolibrary.org/obo/CHEBI_")
    ANASVES = Namespace("http://rdf.example.org/anasves/")
    ANASVES_str = "http://rdf.example.org/anasves/"
    GO = Namespace("http://purl.obolibrary.org/obo/GO_")

    # RDF graph
    g = Graph()
    g.bind("rdf", RDF)
    g.bind("rdfs", RDFS)
    g.bind("rh", RH)
    g.bind("chebi", CHEBI)
    g.bind("slm", SLM)
    g.bind("anasves", ANASVES)

    MNetIRI = ANASVES[MNetIRI_string]
    g.add((MNetIRI, RDFS.label, Literal(MNetIRI_literal)))
    accessions_ready = []

    def parse_equation_sides(equation):
        if pd.isna(equation):
            return [], []
        left, right = equation.split(" = ")
        parse_side = lambda s: [(int(x.split()[0]), x.split()[1]) for x in s.split(" + ")]
        return parse_side(left), parse_side(right)

    def merge_compounds(primary, fallback):
        result = []
        for p, f in zip(primary, fallback):
            stoich_coef = p[0] if p else f[0]
            comp_id = p[1] if p and p[1] != "NA" else f[1]
            result.append((stoich_coef, comp_id))
        return result

    def add_compounds(side_uri, compounds, reaction_uri):
        for stoich_coef, comp_id in compounds:
            if comp_id == "NA":
                continue
            part_uri = URIRef(f"{reaction_uri}_compound_{comp_id}")
            g.add((part_uri, RH.location, GO["0005575"]))
            comp_uri = URIRef(f"{ANASVES_str}Compound_{comp_id}")
            g.add((side_uri, RH[f"contains{stoich_coef}"], part_uri))
            g.add((RH[f"contains{stoich_coef}"], RH.coefficient, Literal(stoich_coef)))
            g.add((RH[f"contains{stoich_coef}"], RDFS.subPropertyOf, RH.contains))
            g.add((part_uri, RH.compound, comp_uri))
            full_accession = comp_id if "SLM:" in comp_id else f"CHEBI:{comp_id}"
            if full_accession.startswith("SLM:"):
                chem_iri = SLM[full_accession.replace("SLM:", "")]
            elif full_accession.startswith("CHEBI:"):
                chem_iri = CHEBI[full_accession.replace("CHEBI:", "")]
            else:
                raise Exception("Unexpected chem name: " + full_accession)
            g.add((comp_uri, RH.accession, Literal(full_accession)))
            g.add((comp_uri, RH.chebi, chem_iri))

    for _, row in df.iterrows():
        accession = row['Web-RInChIKey'].replace("Web-RInChIKey=", "")
        reaction_uri = f"{ANASVES_str}{accession.replace('-', '_')}"
        left_uri = URIRef(f"{reaction_uri}_L")
        right_uri = URIRef(f"{reaction_uri}_R")
        reaction_ref = URIRef(reaction_uri)
        g.add((reaction_ref, RDFS.subClassOf, RH[str(row['MASTER_ID'])]))
        g.add((reaction_ref, ANASVES["generatedFrom"], RH[str(row['MASTER_ID'])]))
        if accession not in accessions_ready:
            g.add((MNetIRI, ANASVES.reac, reaction_ref))
            g.add((reaction_ref, RH.accession, Literal(accession)))
            g.add((reaction_ref, RH.side, left_uri))
            g.add((reaction_ref, RH.side, right_uri))
            g.add((left_uri, RH.curatedOrder, Literal(1)))
            g.add((right_uri, RH.curatedOrder, Literal(2)))

            chebi_l, chebi_r = parse_equation_sides(row['chebi_equation'])
            slm_l, slm_r = parse_equation_sides(row['swisslipids_equation'])
            left_compounds = merge_compounds(slm_l, chebi_l)
            right_compounds = merge_compounds(slm_r, chebi_r)

            add_compounds(left_uri, left_compounds, reaction_uri)
            add_compounds(right_uri, right_compounds, reaction_uri)

            accessions_ready.append(accession)

    print(f"Saving RDF with {len(g)} triples to: {output_file}")
    g.serialize(destination=output_file, format="ttl")