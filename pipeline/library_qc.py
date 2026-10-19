"""Collect per-library salmon QC for every quantified study into one table.
Run on Hazel:  python3 library_qc.py /rsstu/users/r/rrellan/sara/MaizeExpressN/results > library_qc.tsv
Columns: study_id, sample_id, reads_processed, reads_mapped, percent_mapped, library_type, compatible_ratio
(libraries quantified with STAR + salmon, i.e. Huang2025, report salmon's view of the transcriptome BAM)."""
import glob
import json
import os
import sys

R = sys.argv[1]
print("\t".join(["study_id", "sample_id", "reads_processed", "reads_mapped", "percent_mapped", "library_type", "compatible_ratio"]))
for study in sorted(os.listdir(R)):
    if study.startswith("_"):
        continue
    for mi in sorted(glob.glob(os.path.join(R, study, "*", "*", "aux_info", "meta_info.json"))):
        d = os.path.dirname(os.path.dirname(mi))
        if os.path.basename(os.path.dirname(d)) not in ("salmon", "star_salmon"):
            continue
        m = json.load(open(mi))
        lf = os.path.join(d, "lib_format_counts.json")
        l = json.load(open(lf)) if os.path.exists(lf) else {}
        print("\t".join(map(str, [study, os.path.basename(d), m.get("num_processed", ""), m.get("num_mapped", ""),
                                  round(m.get("percent_mapped", 0), 2), l.get("expected_format", ""),
                                  round(l.get("compatible_fragment_ratio", 0), 3)])))
