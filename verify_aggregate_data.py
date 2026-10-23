from pathlib import Path
import csv, math

base = Path(__file__).resolve().parent
def read(name):
    with (base / name).open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))

discovery = read("discovery_symptom_summary.csv")
validation = read("validation_symptom_summary.csv")
regional = read("regional_group_results.csv")
assert len(discovery) == 4 and len(validation) == 2 and len(regional) == 210
assert {r["model"] for r in discovery} == {"six_covariates", "six_plus_group"}
assert all(int(r["n"]) == 102 for r in discovery)
assert all(int(r["n"]) == 29 for r in validation)
for r in discovery + validation:
    assert -1 <= float(r["rho"]) <= 1 and 0 <= float(r["p"]) <= 1
for r in discovery:
    assert float(r["p"]) <= float(r["p_BH"]) <= 1
assert len({r["ROI"] for r in regional}) == 210
assert all(int(r["n_used"]) == 102 for r in regional)
assert sum(float(r["p_FDR"]) < .05 for r in regional) == 9
for r in regional:
    assert math.isclose(float(r["mean_Dep"]) - float(r["mean_HC"]), float(r["mean_diff"]), abs_tol=1e-12)
print("Aggregate tables passed: 4 discovery associations, 2 validation associations, 210 regions; 9 FDR-significant regions.")
