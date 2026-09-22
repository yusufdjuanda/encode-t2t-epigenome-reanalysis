# Authors' ENCODE reanalysis results

Unmodified copies of three files from
[msauria/T2T_Encode_Analysis](https://github.com/msauria/T2T_Encode_Analysis)
at commit `3837e8a3d0` (GPL-3.0), the pipeline behind the ENCODE results in
Gershman et al. (2022). `workflow/scripts/06_compare_with_paper.py` compares
this project's results against them. See `docs/paper-methods.md` for what the
pipeline did.

| File | Source path | sha256 |
|---|---|---|
| `encode_samples.txt` | `data/encode_samples.txt` | `501355ab49d58a1fe82d0c0545097124405a88a07cb1c47a6189c3306bda7598` |
| `mapping_stats.txt` | `results/mapping_stats.txt` | `51f43a3e72bf730dd13c2773883dcd850954d3270b8eca77c7b52bca444a2c4a` |
| `macs2_peak_overlap.txt` | `results/macs2_peak_overlap.txt` | `52ba9a6cafcdb3963881a4eb76ba91881e05df7578a20c7e354a076374c626fe` |

`mapping_stats.txt` counts read pairs. For experiments whose libraries were
sequenced in several runs, its counts are inflated by the number of runs per
library (deviation D4 in `docs/replication-plan.md`); the comparison script
corrects for that.
