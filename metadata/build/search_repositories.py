#!/usr/bin/env python3
"""List maize nitrogen records in NCBI BioProject and GEO, for title matching.

Writes metadata/build/repo_bioproject_maize_N.tsv and repo_geo_maize_N.tsv.
Studies without an accession in the paper are matched against these by hand
(title words, genotypes, tissue, submitter, dates).

Usage: python3 metadata/build/search_repositories.py
"""
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harvest_accessions as h  # noqa: E402

BUILD = h.BUILD
N_TERMS = ('(nitrogen OR nitrate OR ammonium OR "low N" OR "N deficiency" OR '
           '"N starvation" OR "N stress" OR "N supply" OR urea OR NUE)')


def esearch_all(db, term):
    ids, start = [], 0
    while True:
        txt = h.eutil("esearch", db=db, term=term, retmax=500, retstart=start, retmode="json")
        res = json.loads(txt)["esearchresult"]
        ids += res.get("idlist", [])
        start += 500
        if start >= int(res.get("count", 0)):
            return ids


def esummary(db, ids):
    out = []
    for i in range(0, len(ids), 200):
        txt = h.eutil("esummary", db=db, id=",".join(ids[i:i + 200]), retmode="json")
        res = json.loads(txt)["result"]
        out += [res[u] for u in res.get("uids", [])]
    return out


def clean(s):
    return " ".join(str(s).split())


def main():
    ids = esearch_all("bioproject", f'{N_TERMS} AND "Zea mays"[Organism]')
    rows = []
    for d in esummary("bioproject", ids):
        rows.append({
            "accession": d.get("project_acc", ""),
            "registration_date": d.get("registration_date", ""),
            "title": clean(d.get("project_title", "")),
            "submitter": clean(d.get("submitter_organization", "")),
            "data_type": clean(d.get("project_data_type", "")),
            "description": clean(d.get("project_description", ""))[:600],
        })
    path = BUILD / "repo_bioproject_maize_N.tsv"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: r["registration_date"]))
    print(f"{len(rows)} BioProjects -> {path}", file=sys.stderr)

    ids = esearch_all("gds", f'{N_TERMS} AND "Zea mays"[Organism] AND gse[Entry Type]')
    rows = []
    for d in esummary("gds", ids):
        rel = [x.get("targetobject", "") for x in d.get("extrelations", [])]
        rows.append({
            "accession": d.get("accession", ""),
            "pdat": d.get("pdat", ""),
            "n_samples": d.get("n_samples", ""),
            "gdstype": clean(d.get("gdstype", "")),
            "title": clean(d.get("title", "")),
            "bioproject": d.get("bioproject", ""),
            "sra": ";".join(rel),
            "pubmed": ";".join(str(x) for x in d.get("pubmedids", [])),
            "summary": clean(d.get("summary", ""))[:600],
        })
    path = BUILD / "repo_geo_maize_N.tsv"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: r["pdat"]))
    print(f"{len(rows)} GEO series -> {path}", file=sys.stderr)

    # SRA experiments, grouped by BioProject: catches projects whose titles are
    # generic ("Zea mays Raw sequence reads") but whose experiments mention N.
    ids = esearch_all("sra", f'{N_TERMS} AND "Zea mays"[Organism] AND "transcriptomic"[Source]')
    by_proj = {}
    for d in esummary("sra", ids):
        x = d.get("expxml", "")
        proj = re.search(r"<Bioproject>([^<]+)</Bioproject>", x)
        study = re.search(r'<Study acc="([^"]+)" name="([^"]*)"', x)
        title = re.search(r"<Title>(.*?)</Title>", x)
        center = re.search(r'center_name="([^"]*)"', x)
        key = proj.group(1) if proj else (study.group(1) if study else "")
        p = by_proj.setdefault(key, {"bioproject": key, "srp": study.group(1) if study else "",
                                     "study_name": clean(study.group(2)) if study else "",
                                     "center": center.group(1) if center else "",
                                     "n_experiments": 0, "create_date": d.get("createdate", ""),
                                     "example_titles": []})
        p["n_experiments"] += 1
        p["create_date"] = min(p["create_date"], d.get("createdate", "")) or p["create_date"]
        if title and len(p["example_titles"]) < 4:
            p["example_titles"].append(clean(title.group(1)))
    rows = []
    for p in by_proj.values():
        p["example_titles"] = " | ".join(p["example_titles"])
        rows.append(p)
    path = BUILD / "repo_sra_maize_N.tsv"
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: r["create_date"]))
    print(f"{len(rows)} SRA projects ({len(ids)} experiments) -> {path}", file=sys.stderr)


if __name__ == "__main__":
    main()
