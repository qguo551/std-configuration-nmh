# Structure–function configuration in adolescent subthreshold depression

Analysis code and aggregate results accompanying the manuscript “Individualized Atypicality of Structure-Function Configuration Marks Subthreshold Depression in Adolescents”.

## Contents

- `discovery_symptom_associations.m`: the discovery-sample symptom analysis, including the six-covariate model.
- `validation_symptom_associations.m`: the validation StD symptom analysis adjusted for mean framewise displacement and total intracranial volume.
- `discovery_symptom_summary.csv` and `validation_symptom_summary.csv`: aggregate results from the verified analysis outputs, retaining their original numerical precision.
- `regional_group_results.csv`: the saved 210-region discovery ANCOVA result table, including the nine regions with FDR-adjusted P < 0.05.
- `demographics_discovery.csv` and `demographics_validation.csv`: the aggregate demographic tables as reported in the manuscript.
- `DATA_DICTIONARY.md`: variable definitions and input requirements.
- `verify_aggregate_data.py`: checks the structure and consistency of the released aggregate tables using Python's standard library.

This release covers symptom-association code and the aggregate tables listed above. Image preprocessing, network construction, molecular analyses and spin-test implementations are not included in this release. The regional statistics are saved analysis outputs; their model-fitting code is not included here.

## Data access

This repository contains aggregate results only. Participant-level MRI, clinical measurements, network matrices, residuals and identifiers are not included. No participant-level access arrangement is offered by this release.

The Allen Human Brain Atlas and receptor-map resources used in the manuscript are maintained by their original providers: https://human.brain-map.org/ and https://github.com/netneurolab/hansen_receptors. These external datasets are not redistributed here.

## Running the symptom analyses

The code uses MATLAB and the Statistics and Machine Learning Toolbox; validation was performed with MATLAB R2023b. Running these analyses requires the private participant-level workbooks described in `DATA_DICTIONARY.md`.

```matlab
addpath('/path/to/this/repository');
discovery = discovery_symptom_associations('/path/to/discovery_workbook.xlsx');
validation = validation_symptom_associations('/path/to/validation_workbook.xlsx');
```

Both functions return aggregate statistics and print aggregate results. They do not export participant rows. The discovery calculation retains PHQ as the third test in the established FDR correction family; the released symptom tables report SCL-90 and SDS, as in the manuscript.

The discovery reference uses the 32 healthy controls, leaving out the target participant for each healthy control. The validation reference uses its own 16 healthy controls. The symptom correlations are two-sided Spearman correlations of separate ordinary-least-squares residuals. The discovery main model adjusts for age, sex, only-child status, single-parent status, mean framewise displacement and total intracranial volume. Validation adjusts for mean framewise displacement and total intracranial volume within the 29 participants with StD.

## Checking aggregate tables

```sh
python verify_aggregate_data.py
```

This check needs no participant data and does not re-estimate the models.
