def build_rdf(input_file=None, output_dir=None, output_format='nt'):
    
    import os
    import gzip
    import shutil
    from pathlib import Path
    from rdflib import Graph, Namespace, RDF, RDFS, URIRef, Literal
    import pandas as pd
    
    MNetIRI_literal = "SwissLipids Reactions C16"
    MNetIRI_string = "SlrC16"
    
    # Set output directory (default is current working directory).
    if output_dir is None:
        output_dir = os.getcwd()
    else:
        os.makedirs(output_dir, exist_ok=True)

    # Set input file (default is "enumerated_reactions.tsv" in output directory).
    if input_file is None:
        input_file = os.path.join(output_dir, "enumerated_reactions.tsv")
        if not os.path.isfile(input_file) or not os.access(input_file, os.R_OK):
            raise FileNotFoundError(
                f"Default file {input_file} does not exist or is not readable. Please provide file with --input option."
            )
    # Set output file name to be the same as input file name, but with file extension for output format.
    input_file_base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_file = os.path.join(output_dir, f"{input_file_base_name}.{output_format}")

    # Load dataframe.
    df = pd.read_csv(input_file, sep="\t")
    df.dropna(subset=['Web-RInChIKey'], inplace=True)

    # Namespaces
    RH = Namespace("http://rdf.rhea-db.org/")
    SLM = Namespace("https://swisslipids.org/rdf/SLM_")
    CHEBI = Namespace("http://purl.obolibrary.org/obo/CHEBI_")
    SLR_purl = "https://purl.expasy.org/lipid/reaction#"
    SLR = Namespace(SLR_purl)
    GO = Namespace("http://purl.obolibrary.org/obo/GO_")

    # RDF graph
    g = Graph()
    g.bind("rdf", RDF)
    g.bind("rdfs", RDFS)
    g.bind("rh", RH)
    g.bind("chebi", CHEBI)
    g.bind("slm", SLM)
    g.bind("slr", SLR)

    MNetIRI = SLR[MNetIRI_string]
    g.add((MNetIRI, RDFS.label, Literal(MNetIRI_literal)))
    accessions_ready = []

    def parse_equation_sides(equation):
        if pd.isna(equation):
            return [], []
        left, right = equation.split(" = ")
        parse_side = lambda s: [(int(x.split()[0]), x.split()[1]) for x in s.split(" + ")]
        return parse_side(left), parse_side(right)

    # Maps participants from 'slm_id_equation' to participants from 'chebi_equation'.
    # This is required to add a subClassOf <CHEBI> to <SLM> participants.
    def map_participants(slm, chebi):
        result = []
        for s, c in zip(slm, chebi):
            stoich_coef = s[0]
            slm_id      = s[1]
            chebi_id    = c[1]
            result.append((stoich_coef, slm_id, chebi_id))
        return result

    def add_participants(side_uri, compounds, reaction_uri):
        
        # NOTE: 'slm_id' is an ID from 'slm_id_equation', which may contain both
        # SLM IDs and CHEBI IDs!
        for stoich_coef, slm_id, chebi_id in compounds:
            # 1. Add participant.
            part_uri = URIRef(f"{reaction_uri}_compound_{slm_id}")
            g.add((part_uri, RH.location, GO["0005575"]))
            g.add((side_uri, RH[f"contains{stoich_coef}"], part_uri))
            g.add((RH[f"contains{stoich_coef}"], RH.coefficient, Literal(stoich_coef)))
            g.add((RH[f"contains{stoich_coef}"], RDFS.subPropertyOf, RH.contains))

            # 2. Add compound.
            comp_uri = URIRef(f"{SLR_purl}Compound_{slm_id}")
            g.add((part_uri, RH.compound, comp_uri))
            full_accession = slm_id if slm_id.startswith("SLM:") else f"CHEBI:{chebi_id}"
            if full_accession.startswith("SLM:"):
                chem_iri = SLM[full_accession.replace("SLM:", "")]
            elif full_accession.startswith("CHEBI:"):
                chem_iri = CHEBI[chebi_id]
            else:
                raise Exception("Unexpected chem name: " + full_accession)
            g.add((comp_uri, RH.accession, Literal(full_accession)))
            g.add((comp_uri, RH.chebi, chem_iri))
            g.add((comp_uri, RDFS.subClassOf, CHEBI[chebi_id]))

    for _, row in df.iterrows():
        
        accession = row['Web-RInChIKey'].replace("Web-RInChIKey=", "")
        reaction_uri = f"{SLR_purl}{accession.replace('-', '_')}"
        uri_l = URIRef(f"{reaction_uri}_L")
        uri_r = URIRef(f"{reaction_uri}_R")
        reaction_ref = URIRef(reaction_uri)
        
        g.add((reaction_ref, RDFS.subClassOf, RH[str(row['MASTER_ID'])]))
        g.add((reaction_ref, SLR["generatedFrom"], RH[str(row['MASTER_ID'])]))
        
        if accession not in accessions_ready:
            g.add((MNetIRI, SLR.reac, reaction_ref))
            g.add((reaction_ref, RH.accession, Literal(accession)))
            g.add((reaction_ref, RH.side, uri_l))
            g.add((reaction_ref, RH.side, uri_r))
            g.add((uri_l, RH.curatedOrder, Literal(1)))
            g.add((uri_r, RH.curatedOrder, Literal(2)))

            slm_l, slm_r     = parse_equation_sides(row['slm_id_equation'])
            chebi_l, chebi_r = parse_equation_sides(row['chebi_equation'])

            participants_l = map_participants(slm_l, chebi_l)
            participants_r = map_participants(slm_r, chebi_r)

            add_participants(uri_l, participants_l, reaction_uri)
            add_participants(uri_r, participants_r, reaction_uri)

            accessions_ready.append(accession)

    print(f"Saving RDF with {len(g)} triples to {output_file}.gz")
    g.serialize(destination=output_file, format=output_format, encoding='utf-8')

    # Compress the output file.
    with open(output_file, "rb") as f_in:
        with gzip.open(output_file + ".gz", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(output_file)
