#!/usr/bin/env python3
"""Merge per-study sample sheets into metadata/samples.tsv and needs_review.tsv.

Reads metadata/build/samples/<study_id>.tsv and <study_id>.review.tsv (format in
curation_spec.md), validates them against metadata/runs.tsv and the controlled
vocabularies, and prints every problem found. Writes the merged tables only
when called with --write.

Usage: python3 metadata/build/merge_samples.py [--write]
"""
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SAMPLES = ROOT / "metadata" / "build" / "samples"
COLS = ("sample_id study_id run_accessions genotype genotype_class background tissue_group "
        "tissue_detail po_term stage stage_class age_days growth_system N_label_paper N_form "
        "N_amount N_unit N_rank_in_study phase time_h other_factor other_factor_level replicate "
        "pooled_plants include curation_note").split()
REVIEW_COLS = ["study_id", "item", "question", "proposed_answer", "evidence"]
VOCAB = {
    "genotype_class": {"B73", "other inbred", "hybrid", "mutant", "transgenic", "NR"},
    "tissue_group": {"Root", "Leaf", "Shoot", "Ear/kernel", "Whole seedling"},
    "stage_class": {"seedling", "vegetative", "reproductive", "NR"},
    "N_form": {"nitrate", "ammonium", "ammonium nitrate", "urea", "mixed", "none", "NR"},
    "phase": {"steady", "starvation", "resupply"},
    "include": {"TRUE", "FALSE"},
}
NUMERIC = {"N_amount": True, "time_h": True, "age_days": True, "pooled_plants": True}


def is_num(x):
    return bool(re.fullmatch(r"-?\d+(\.\d+)?", x))


def main(write=False):
    runs = list(csv.DictReader(open(ROOT / "metadata" / "runs.tsv"), delimiter="\t"))
    want = defaultdict(set)
    run_study = {}
    for r in runs:
        run_study[r["run_accession"]] = (r["study_id"], r["include"])
        if r["include"] == "TRUE":
            want[r["study_id"]].add(r["run_accession"])

    problems, merged, reviews = [], [], []
    seen_ids = Counter()
    files = sorted(p for p in SAMPLES.glob("*.tsv") if not p.name.endswith(".review.tsv"))
    for p in files:
        study = p.stem
        rows = list(csv.reader(open(p), delimiter="\t"))
        if not rows or rows[0] != COLS:
            problems.append(f"{study}: header differs from spec: {rows[0] if rows else 'empty'}")
            continue
        got = Counter()
        for i, vals in enumerate(rows[1:], 2):
            if len(vals) != len(COLS):
                problems.append(f"{study} line {i}: {len(vals)} fields, expected {len(COLS)}")
                continue
            r = dict(zip(COLS, vals))
            seen_ids[r["sample_id"]] += 1
            if r["study_id"] != study:
                problems.append(f"{study} line {i}: study_id {r['study_id']}")
            if not r["sample_id"].startswith(study + "_") or " " in r["sample_id"]:
                problems.append(f"{study} line {i}: sample_id {r['sample_id']!r}")
            for col, vocab in VOCAB.items():
                if r[col] not in vocab:
                    problems.append(f"{study} {r['sample_id']}: {col} = {r[col]!r}")
            for col in NUMERIC:
                if r[col] and r[col] != "NR" and not is_num(r[col]):
                    problems.append(f"{study} {r['sample_id']}: {col} = {r[col]!r} not numeric")
            if r["N_rank_in_study"] != "NR" and not r["N_rank_in_study"].isdigit():
                problems.append(f"{study} {r['sample_id']}: N_rank_in_study = {r['N_rank_in_study']!r}")
            if not r["replicate"].isdigit():
                problems.append(f"{study} {r['sample_id']}: replicate = {r['replicate']!r}")
            if r["N_unit"] not in ("mM N", "kg N/ha") and r["N_amount"] not in ("", "NR"):
                problems.append(f"{study} {r['sample_id']}: N_unit = {r['N_unit']!r} (expect 'mM N' or 'kg N/ha')")
            if r["phase"] in ("starvation", "resupply") and r["time_h"] == "":
                problems.append(f"{study} {r['sample_id']}: {r['phase']} without time_h")
            for run in filter(None, r["run_accessions"].split(";")):
                got[run] += 1
                owner = run_study.get(run)
                if owner is None:
                    problems.append(f"{study} {r['sample_id']}: unknown run {run}")
                elif owner[0] != study:
                    problems.append(f"{study} {r['sample_id']}: run {run} belongs to {owner[0]}")
                elif owner[1] != "TRUE" and r["include"] == "TRUE":
                    problems.append(f"{study} {r['sample_id']}: run {run} is excluded in runs.tsv")
            merged.append(r)
        for run in want[study] - set(got):
            problems.append(f"{study}: included run {run} not assigned to a sample")
        for run, n in got.items():
            if n > 1:
                problems.append(f"{study}: run {run} assigned {n} times")
        rp = SAMPLES / f"{study}.review.tsv"
        if rp.exists():
            rr = list(csv.DictReader(open(rp), delimiter="\t"))
            if rr and list(rr[0]) != REVIEW_COLS:
                problems.append(f"{study}: review header {list(rr[0])}")
            reviews += [{k: x.get(k, "") or "" for k in REVIEW_COLS} for x in rr]
        else:
            problems.append(f"{study}: no review file")
    for sid, n in seen_ids.items():
        if n > 1:
            problems.append(f"duplicate sample_id {sid}")
    missing = sorted(set(want) - {p.stem for p in files})
    for s in missing:
        problems.append(f"{s}: no sample sheet yet ({len(want[s])} included runs)")

    print(f"{len(files)} sample sheets, {len(merged)} samples, {len(reviews)} review items")
    print(f"{len(problems)} problems")
    for p in problems:
        print("  -", p)
    if write:
        with open(ROOT / "metadata" / "samples.tsv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t")
            w.writeheader()
            w.writerows(merged)
        extra = list(csv.DictReader(open(ROOT / "metadata" / "build" / "review_global.tsv"),
                                    delimiter="\t"))
        with open(ROOT / "metadata" / "needs_review.tsv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=REVIEW_COLS, delimiter="\t")
            w.writeheader()
            w.writerows(extra + reviews)
        print("wrote metadata/samples.tsv and metadata/needs_review.tsv")
    return problems


if __name__ == "__main__":
    main("--write" in sys.argv)
