It is necessary to download lipids.tsv from [SwissLipids](https://www.swisslipids.org/#/downloads).

lipids.tsv (~700MB) should be located in src/swisslipidsreact/package_data before starting the execution.

## Usage
By default, the pipeline will generate results for palmitate only.
To generate results for the whole list of fatty acids in human and enumerated classes, set --curated-fa option.

To learn more about the options, check swisslipidsreact --help.

* run pipeline for C16:
```bash
swisslipidsreact run --output-dir results_C16/
```

* run pipeline for curated list of fatty acids (execution time: several hours):
```bash
swisslipidsreact run --curated-fa --output-dir results_curated_fatty_acids/
```

* run pipeline for all fatty acids (execution time: ∞):
```bash
swisslipidsreact run --all-fa --output-dir results_all_fatty_acids/
```

* run pipeline for all fatty acids per rhea id (execution time: ∞):
```bash
swisslipidsreact run --all-fa --output-dir results_78071/ --rheaid 78071
```

* turtle (RDF type) export for palmitate:
```bash
swisslipidsreact export-ttl --output-dir results_C16/
```
* turtle (RDF type) export for curated list of fatty acids (execution time: several hours):
```bash
swisslipidsreact export-ttl --curated-fa --output-dir results_curated_fatty_acids/
```
