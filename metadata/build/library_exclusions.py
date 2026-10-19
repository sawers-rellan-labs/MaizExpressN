#!/usr/bin/env python3
"""Flag low-quality libraries for the DESeq2 step (M4). Strict by design: quality over quantity
(Rubén, 19 Oct 2026).

Reads metadata/library_qc.tsv (pipeline/library_qc.py), metadata/samples.tsv and the gene
count tables in results/<study>/, and writes metadata/library_exclusions.tsv
(study_id, sample_id, reason). A library is excluded when any of these holds:
  1. depth: fewer than 5 M reads mapped by salmon (2 M for 3' tag studies);
  2. mapping: below 60%, or more than 10 points below the median of the same study, tissue group and genotype
     (per genotype, because non-B73 lines map worse to the B73 reference for genetic, not technical, reasons);
  3. library type: differs from the study's majority type (a different or failed prep);
  4. replicates (same group key as build_contrasts.py; log2-CPM correlation over genes with mean CPM > 5):
     - in a group of >= 3, its mean r with the other replicates is > 0.03 below the study's
       median within-group r AND > 0.03 below the other replicates' r among themselves;
     - in a group of 2, the pair's r is > 0.05 below the study's median within-group r
       (both go: there is no telling which one is wrong).
Excluded libraries are left out of every DESeq2 fit; the samples stay in samples.tsv.
Contrasts that lose replication become display-only (de_per_study.R treats them as unreplicated).
Usage: python3 metadata/build/library_exclusions.py
"""
import collections
import csv
import glob
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MIN_MAPPED = 5_000_000
MIN_MAPPED_TAG = 2_000_000
TAG_STUDIES = {"Huang2025"}           # 3' tag-seq: fewer reads needed; replicate rule not applied
MIN_PCT = 60.0
MAX_BELOW_MEDIAN = 10.0
R_MARGIN_3 = 0.03
R_MARGIN_2 = 0.05
GROUP_KEY = ["study_id", "genotype", "tissue_group", "tissue_detail", "N_label_paper", "N_form", "N_amount",
             "N_unit", "N_rank_in_study", "phase", "time_h", "other_factor", "other_factor_level"]

qc_rows = list(csv.DictReader(open(ROOT / "metadata" / "library_qc.tsv"), delimiter="\t"))
smp = {s["sample_id"]: s for s in csv.DictReader(open(ROOT / "metadata" / "samples.tsv"), delimiter="\t")
       if s["include"] == "TRUE"}
qc_rows = [r for r in qc_rows if r["sample_id"] in smp]
group = {k: "|".join(v[c] for c in GROUP_KEY) for k, v in smp.items()}
reasons = collections.defaultdict(list)

# rules 1-3
by_study = collections.defaultdict(list)
for r in qc_rows:
    by_study[r["study_id"]].append(r)
for study, rs in by_study.items():
    major = collections.Counter(r["library_type"] for r in rs).most_common(1)[0][0]
    by_tissue = collections.defaultdict(list)
    tg = lambda sid: (smp[sid]["tissue_group"], smp[sid]["genotype"])
    for r in rs:
        by_tissue[tg(r["sample_id"])].append(float(r["percent_mapped"]))
    for r in rs:
        sid, pct, mapped = r["sample_id"], float(r["percent_mapped"]), int(r["reads_mapped"])
        floor = MIN_MAPPED_TAG if study in TAG_STUDIES else MIN_MAPPED
        if mapped < floor:
            reasons[sid].append(f"{mapped / 1e6:.1f} M mapped reads (< {floor / 1e6:.0f} M)")
        med = statistics.median(by_tissue[tg(sid)])
        if pct < MIN_PCT:
            reasons[sid].append(f"{pct:.1f}% mapped (< {MIN_PCT:.0f}%)")
        elif pct < med - MAX_BELOW_MEDIAN:
            reasons[sid].append(f"{pct:.1f}% mapped vs {med:.1f}% median for its study, tissue and genotype")
        if r["library_type"] != major:
            reasons[sid].append(f"library type {r['library_type']} vs study {major}")

# rule 4
study_of = {r["sample_id"]: r["study_id"] for r in qc_rows}
for path in sorted(glob.glob(str(ROOT / "results" / "*" / "salmon.merged.gene_counts.tsv"))):
    study = Path(path).parent.name
    if study in TAG_STUDIES:
        continue
    tab = list(csv.reader(open(path), delimiter="\t"))
    head = tab[0]
    cols = [j for j in range(2, len(head)) if head[j] in smp and head[j] not in reasons]
    if len(cols) < 2:
        continue
    tot = {j: sum(float(r[j]) for r in tab[1:]) for j in cols}
    lcpm = []
    for r in tab[1:]:
        c = [float(r[j]) / tot[j] * 1e6 for j in cols]
        if sum(c) / len(c) > 5:
            lcpm.append([math.log2(x + 1) for x in c])
    X = list(zip(*lcpm)); names = [head[j] for j in cols]
    mu = [sum(x) / len(x) for x in X]
    sd = [math.sqrt(sum((a - m) ** 2 for a in x)) for x, m in zip(X, mu)]
    memo = {}
    def cor(i, k):
        key = (min(i, k), max(i, k))
        if key not in memo:
            memo[key] = sum((a - mu[i]) * (b - mu[k]) for a, b in zip(X[i], X[k])) / (sd[i] * sd[k])
        return memo[key]
    groups = collections.defaultdict(list)
    for i, nm in enumerate(names):
        groups[group[nm]].append(i)
    within = [cor(a, b) for g in groups.values() for x, a in enumerate(g) for b in g[x + 1:]]
    if not within:
        continue
    wmed = statistics.median(within)
    for g in groups.values():
        if len(g) == 2:
            r = cor(g[0], g[1])
            if r < wmed - R_MARGIN_2:
                for i in g:
                    reasons[names[i]].append(f"its 2 replicates disagree: r = {r:.3f} vs study within-group median {wmed:.3f}")
        elif len(g) >= 3:
            for i in g:
                others = [k for k in g if k != i]
                r_i = statistics.mean(cor(i, k) for k in others)
                r_oth = statistics.mean(cor(a, b) for x, a in enumerate(others) for b in others[x + 1:])
                if r_i < wmed - R_MARGIN_3 and r_oth - r_i > R_MARGIN_3:
                    reasons[names[i]].append(f"disagrees with its replicates: r = {r_i:.3f} with them vs "
                                             f"{r_oth:.3f} among them (study within-group median {wmed:.3f})")

out = sorted([study_of.get(s, smp[s]["study_id"]), s, "; ".join(v)] for s, v in reasons.items())
with open(ROOT / "metadata" / "library_exclusions.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t", lineterminator="\n")
    w.writerow(["study_id", "sample_id", "reason"])
    w.writerows(out)

# impact on contrasts
excluded = {o[1] for o in out}
n_lib = collections.Counter(r["study_id"] for r in qc_rows)
n_out = collections.Counter(o[0] for o in out)
print(f"{len(out)} of {len(qc_rows)} libraries excluded -> metadata/library_exclusions.tsv")
lost = []
for c in csv.DictReader(open(ROOT / "metadata" / "contrasts.tsv"), delimiter="\t"):
    num = [s for s in c["numerator_samples"].split(";") if s and s not in excluded]
    den = [s for s in c["denominator_samples"].split(";") if s and s not in excluded]
    if not num or not den:
        lost.append((c["study_id"], c["contrast_id"], "no samples left"))
    elif c["replicated"] == "TRUE" and (len(num) < 2 or len(den) < 2):
        lost.append((c["study_id"], c["contrast_id"], "now unreplicated (display only)"))
print("per study (excluded / libraries):", ", ".join(f"{s} {n_out[s]}/{n_lib[s]}" for s in sorted(n_out)))
print(f"{len(lost)} contrasts affected:")
for l in lost:
    print("  " + "\t".join(l))
