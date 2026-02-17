# SwissLipidsReact

Expands Rhea reaction patterns into complete lipid reactions, resolving structures and assigning RInChIs.

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
# Run enumeration
swisslipidsreact run

# Export .ttl (turtle) format for integration into the RDF knowledge graph.
swisslipidsreact export-ttl

# Analyse rhea reaction template usage.
swisslipidsreact master-id-analysis
```

## Options

Explanation of fatty acid options:
  
Option | Meaning | Runtime | Usage |
--- | --- | --- | --- |
none (default) | Only palmitate allowed as a fatty acid in any position | minutes | Testing with reduced dataset |
--curated-fa | Filter SwissLipids based on allowed FA per position | hours | Filtered for integration in RDF knowledge grap |
-all-fa | all SwissLipids considered | ∞ | not recommended (too slow), but can be used for an individual Rhea ID |

*Reaction enumeration*

  ```bash
"--output-dir",
type=str,
default=None,
help="Output directory (default: current working directory)"

"--curated-fa",
action="store_true",
help="Use curated fatty acid list (default: False for C16)"

"--all-fa",
action="store_true",
default=False,
help="No restrictions of FA per position"

"--rheaid",
type=int,
default=None,
help="run pipeline for only one rhea id"
```

*RDF export*

```bash
"--curated-fa",
action="store_true",
help="Use curated fatty acid list for TTL export (default: False for C16)"

"--input",
type=str,
default=None,
help="Input TSV file (default: inferred from mode)"
```

## Usage

By default, the pipeline will generate results for palmitate only.
To generate results for the whole list of fatty acids in human and enumerated classes, use the `--curated-fa` option.

To learn more about the options, check `swisslipidsreact --help`.

* Enumerate with C16 fatty acids test set:
  ```bash
  swisslipidsreact run --output-dir results_C16/
  ```

* Enumerate with curated list of fatty acids (execution time: several hours):
  ```bash
  swisslipidsreact run --curated-fa --output-dir results_curated_fatty_acids/
  ```

* Enumerate with all fatty acids (WARNING: execution time: ∞):
  ```bash
  swisslipidsreact run --all-fa --output-dir results_all_fatty_acids/
  ```

* Enumerate with all fatty acids for one rhea id:
  ```bash
  swisslipidsreact run --all-fa --output-dir results_78071/ --rheaid 78071
  ```

* Export RDF for C16 test set:
  ```bash
  swisslipidsreact export-ttl --output-dir results_C16/
  ```
* Export RDF for curated list of fatty acids (execution time: several hours):
  ```bash
  swisslipidsreact export-ttl --curated-fa --output-dir results_curated_fatty_acids/
  ```

* Analyse the Rhea reaction master id usage:
  ```bash
  swisslipidsreact master-id-analysis --input "results_merged/merged_enumerated_reactions.tsv" --all-fa
  ```

## Debugging

Use the environment variable SLR_DEBUG to get more detailed debug information, e.g.:

```bash
SLR_DEBUG=1 swisslipidsreact run --output-dir results_C16
```

* SLR_DEBUG=1 prints debug messages.
* SLR_DEBUG=2 serializes various dataframes into DEBUG_...tsv files (this will take disk space, use only in test mode).


## Profiling

```bash
pip install pyinstrument
pyinstrument --from-path swisslipidsreact export-ttl -input ... --output-dir ...
```
