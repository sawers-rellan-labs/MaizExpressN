#!/usr/bin/env python3
"""Check metadata/contrasts.tsv against metadata/samples.tsv.

Checks: columns; class vocabulary; sample ids exist, are included and belong
to the contrast's study; numerator and denominator are disjoint and n counts
match; sign convention (lower N in the numerator; nitrate in the denominator
for form contrasts); held-constant factors really are equal on both sides;
no two pooled contrasts of a study use the same samples. Also lists included
samples that no contrast uses, for review. Usage:
python3 metadata/build/validate_contrasts.py
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_contrasts import COLS, ROOT  # noqa: E402

CLASSES = {"deficiency", "starvation", "resupply", "excess", "field_rate", "local", "form"}


def num(x):
    try:
        return float(x)
    except ValueError:
        return None


def main():
    samples = {s["sample_id"]: s for s in csv.DictReader(
        open(ROOT / "metadata" / "samples.tsv"), delimiter="\t")}
    rows = list(csv.reader(open(ROOT / "metadata" / "contrasts.tsv"), delimiter="\t"))
    problems = []
    if rows[0] != COLS:
        problems.append(f"header differs: {rows[0]}")
    contrasts = [dict(zip(rows[0], r)) for r in rows[1:]]
    used = set()
    seen_sets = defaultdict(dict)
    for c in contrasts:
        cid = c["contrast_id"]
        tags = set(filter(None, c["tags"].split(";")))
        if c["contrast_class"] not in CLASSES:
            problems.append(f"{cid}: class {c['contrast_class']}")
        ns = c["numerator_samples"].split(";")
        ds = c["denominator_samples"].split(";")
        if set(ns) & set(ds):
            problems.append(f"{cid}: numerator and denominator overlap")
        if len(ns) != int(c["n_numerator"]) or len(ds) != int(c["n_denominator"]):
            problems.append(f"{cid}: n counts do not match sample lists")
        side = {}
        for name, ids in (("num", ns), ("den", ds)):
            for i in ids:
                s = samples.get(i)
                if s is None:
                    problems.append(f"{cid}: unknown sample {i}")
                    continue
                if s["include"] != "TRUE":
                    problems.append(f"{cid}: excluded sample {i}")
                if s["study_id"] != c["study_id"]:
                    problems.append(f"{cid}: sample {i} belongs to {s['study_id']}")
                used.add(i)
            side[name] = [samples[i] for i in ids if i in samples]
        if not side.get("num") or not side.get("den"):
            continue
        a, b = side["num"][0], side["den"][0]
        # every sample of one side shares the design fields
        for name in ("num", "den"):
            fields = ("genotype", "tissue_detail", "N_label_paper", "N_amount", "phase",
                      "time_h", "other_factor_level")
            if len({tuple(s[f] for f in fields) for s in side[name]}) > 1:
                problems.append(f"{cid}: {name} samples are not one condition group")
        # held constant
        relaxed = "split_vs_homogeneous" in tags  # split-root half vs homogeneous roots
        for f in ("genotype",) if relaxed else ("genotype", "tissue_detail", "other_factor_level"):
            if a[f] != b[f]:
                problems.append(f"{cid}: {f} differs ({a[f]!r} vs {b[f]!r})")
        if "not_time_matched" not in tags and a["time_h"] != b["time_h"]:
            problems.append(f"{cid}: time differs ({a['time_h']} vs {b['time_h']}) without tag")
        # sign convention
        if c["contrast_class"] == "form":
            if b["N_form"] != "nitrate" or a["N_form"] == "nitrate":
                problems.append(f"{cid}: form contrast must be alternative / nitrate")
            if a["N_amount"] != b["N_amount"]:
                problems.append(f"{cid}: form contrast with different total N")
        elif "systemic_signal" not in tags:
            if int(a["N_rank_in_study"]) >= int(b["N_rank_in_study"]):
                problems.append(f"{cid}: numerator rank {a['N_rank_in_study']} not below "
                                f"denominator {b['N_rank_in_study']}")
            x, y = num(a["N_amount"]), num(b["N_amount"])
            if x is not None and y is not None and a["N_unit"] == b["N_unit"] and x >= y:
                problems.append(f"{cid}: numerator N {x} not below denominator N {y}")
        if c["pool"] == "TRUE":
            key = (frozenset(ns), frozenset(ds))
            if key in seen_sets[c["study_id"]]:
                problems.append(f"{cid}: same samples as {seen_sets[c['study_id']][key]}")
            seen_sets[c["study_id"]][key] = cid

    unused = defaultdict(list)
    for i, s in samples.items():
        if s["include"] == "TRUE" and i not in used:
            unused[s["study_id"]].append(i)
    print(f"{len(contrasts)} contrasts; {len(problems)} problems")
    for p in problems:
        print("  -", p)
    print(f"included samples used by no contrast: {sum(map(len, unused.values()))}")
    for st, ids in sorted(unused.items()):
        print(f"  {st} ({len(ids)}): {', '.join(sorted(ids))[:200]}")
    return problems


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
