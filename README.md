# SwissLipidsReact

Expands Rhea reaction patterns into complete lipid reactions, resolving structures and assigning Web-RInChIs.

[![License](https://img.shields.io/github/license/sib-swiss/pyrheadb)](LICENSE)
![OS Linux](https://img.shields.io/badge/OS-Linux-green)
![OS Windows](https://img.shields.io/badge/OS-Windows-blue)
![OS macOS](https://img.shields.io/badge/OS-macOS-lightgrey)

## Description

This code combines the [Rhea](https://www.rhea-db.org) database of biochemical reactions and the [SwissLipids](https://www.swisslipids.org/#/) database of lipid structures to enumerate the hypothetically possible space of biochemical reactions with fully defined lipid structures.

The subset of Rhea reactions that define the lipid reaction mechanisms are represented using the [ChEBI](https://www.ebi.ac.uk/chebi/) identifiers of the reacting lipid classes in the Rhea database.

SwissLipids provides connections between a lipid class - a hypothetical entity aiming to represent many lipids present in nature that share a particular substructure - and all of the hypothetically possible lipid structures with isomeric subspecies level of compound structure definition, i.e. 2.5D structure definition, allowing to recognise precisely atom composition and bond order, as well as stereochemical tags of the atoms of every molecule.

This code transforms each Rhea reaction that is defined in terms of lipid classes into a set of reactions where each reactant and product has a defined 2.5D structure, and checks the correspondance between reactants and products to ensure that the resulting reactions are atomically balanced and biochemically feasible.

## Data

It is necessary to download lipids.tsv (~700MB) from [SwissLipids](https://www.swisslipids.org/#/downloads) and copy it to `src/swisslipidsreact/package_data` before starting the execution.

## Installation

```bash
pip install .
```

## pyrheadb dependency

This package is dependent on [pyrheadb](https://github.com/sib-swiss/pyrheadb/tree/main).

To avoid downloading and preprocessing the full Rhea reaction data for every potential new execution, follow these [instructions](https://github.com/sib-swiss/pyrheadb/wiki) on how to set up the RHEADB_LOC environment variable.

## Run

```bash
# Enumerate reactions.
swisslipidsreact run

# Build RDF from enumeration results for integration into the RDF knowledge graph.
swisslipidsreact build-rdf
```

## Options

Explanation of fatty acid (FA) options:
  
Options | Meaning | Runtime | Usage |
--- | --- | --- | --- |
-filter-fa c16 --test | Use only SwissLipids compounds whose FAs are all palmitate | minutes | Testing with reduced dataset |
-filter-fa c16 | Use only SwissLipids compounds with maximum one FA that is not palmitate | hours | Integration in RDF knowledge graph |
-filter-fa none | Use all SwissLipids compounds | ∞ | Not recommended (too slow), but can be used in combination with the --rhea-id option |

*Reaction enumeration*

  ```bash
"--output-dir",
help = "Output directory (default: current working directory)"

"--filter-fa",
 help = "Filter the fatty acids: c16 (default), curated, none (use only in combination with --rhea-id option)"

"--filter-rhea",
help = "Filter Rhea by having a direct SLM parent class of an isomeric subspecies on at least one or both sides of the reaction: two-sides (default), one-side"

"--rhea-id",
help = "Enumerate reactions only for the given Rhea ID"

"--rhea-version",
help = "Use the given Rhea release version (default: latest release)"

"--test",
help = "Use only SwissLipids compounds whose FAs are all palmitate (default: False)"
```

*RDF build*

```bash
"--input",
 help = "Input TSV file (default: <output-dir>/enumerated_reactions.tsv)"
 
"--output-dir",
help = "Output directory (default: current working directory)"

"--output-format",
help = "RDF serialization format (default: nt)"
```

## Usage examples

To learn more about the options, check `swisslipidsreact --help`.

* Enumerate with SwissLipids compounds whose FAs are all palmitate (test set):
  ```bash
  swisslipidsreact run --filter-fa c16 --output-dir results-test-c16 --test
  ```

* Enumerate with SwissLipids compounds with maximum one FA that is not palmitate (production set):
  ```bash
  swisslipidsreact run --filter-fa c16 --output-dir results-prod-c16
  ```

* Enumerate with all SwissLipids compounds for one rhea ID:
  ```bash
  swisslipidsreact run --filter-fa none --rhea-id 78071 --output-dir results-rhea-78071
  ```

* Build RDF for test set:
  ```bash
  swisslipidsreact build-rdf --output-dir results-test-c16
  ```
* Build RDF for production set:
  ```bash
  swisslipidsreact build-rdf --output-dir results-prod-c16
  ```

## Debugging

Use the environment variable SLR_DEBUG to get more detailed debug information, e.g.:

```bash
SLR_DEBUG=1 swisslipidsreact run --filter-fa curated --output-dir results-test-curated --test
```

* SLR_DEBUG=1 prints debug messages.
* SLR_DEBUG=2 serializes various dataframes into DEBUG_...tsv files (this will take disk space, use only in test mode).


## Profiling

```bash
pip install pyinstrument
pyinstrument --from-path swisslipidsreact build-rdf -input ... --output-dir ...
```
