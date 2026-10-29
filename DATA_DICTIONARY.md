# Data dictionary

## Symptom summary tables

`model`: `six_covariates` is the primary discovery model; `FD_TIV` is the validation model. `scale` identifies SCL-90 or SDS. `n` is the number of participants. `rho` is Spearman's correlation of model residuals. `p` is the two-sided MATLAB Spearman P value. `p_BH` is the Benjamini–Hochberg adjusted P value across SCL-90, SDS and PHQ within the corresponding discovery model. Validation P values are unadjusted.

## Regional group results

Rows correspond to the 210 cortical regions in the original BNA-order result table. `ROI` is the original region-column identifier. `n_used` is the analysis sample size. `mean_Dep` and `mean_HC` are the saved group means in the original configuration scale; `mean_diff` is their StD-minus-HC difference. `t_group`, `F_group`, `p_group`, `partial_eta2`, `cohen_d`, `beta_group` and `p_FDR` retain the names and numerical values in the saved ANCOVA output. The signed coefficient follows the original group coding (StD = 1; HC = 2); it should not be treated as a StD-minus-HC mean difference. The accompanying mean difference supplies that direction directly.

## Demographic tables

Values are transcribed directly from the manuscript tables. Continuous measures use mean ± standard deviation, and categorical measures use count (percentage). The discovery sample contains 70 StD and 32 HC participants. The validation demographic table describes the 29 StD participants.

## Private inputs required by the MATLAB functions

The functions accept workbooks with the original column names. Participant rows are not distributed in this repository.

- `sub`: a unique participant identifier used internally for consistency checks.
- `组别（抑郁组=1，健康组=2）`: StD = 1; HC = 2.
- `sex`: the original binary coding (1 or 2).
- `age`: age in years.
- `独生子女（是=1，否=0）`: only-child status (yes = 1; no = 0).
- `单亲家庭（是=1.否=0）`: single-parent status (yes = 1; no = 0).
- `mean_FD_Power_mm`: mean framewise displacement in millimetres.
- `TIV_mm3`: total intracranial volume in cubic millimetres.
- `SCL90`, `SDS`, `PHQ`: original symptom scores. PHQ is required only by the discovery function for the established three-test FDR family.
- `MINDFIND_node_001` through `MINDFIND_node_210`: the regional structure–function configuration values, in the original BNA order.

The functions preserve the complete-case selection and fixed healthy-reference definitions of the source analyses. No new exclusion rule is introduced.
