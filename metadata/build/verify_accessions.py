#!/usr/bin/env python3
"""Fetch repository records for candidate accessions and summarise them.

Input:  metadata/build/accession_candidates.tsv (study_id, accession, source)
Output: metadata/build/accession_records.tsv (one row per accession)
        metadata/build/ena_runs/<accession>.tsv (ENA run report, metadata only)

Record types handled: NCBI BioProject (PRJNA/PRJEB/PRJDB), SRA study or
experiment (SRP/ERP/DRP/SRX), GEO series (GSE), CNCB GSA (CRA, PRJCA; HTML
page scrape). No sequence data is downloaded.

Usage: python3 metadata/build/verify_accessions.py
"""
import csv
import html
import json
import re
import sys
import urllib.parse
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import harvest_accessions as h  # noqa: E402

BUILD = h.BUILD
RUNS = BUILD / "ena_runs"
ENA_FIELDS = ("run_accession,study_accession,secondary_study_accession,sample_accession,"
              "secondary_sample_accession,experiment_accession,scientific_name,tax_id,"
              "library_strategy,library_source,library_selection,library_layout,"
              "instrument_model,read_count,base_count,sample_title,experiment_title,"
              "study_title,first_public,fastq_ftp,fastq_md5,fastq_bytes")


def clean(s):
    return " ".join(html.unescape(str(s)).split())


def ena_runs(acc):
    url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?"
           + urllib.parse.urlencode({"accession": acc, "result": "read_run",
                                     "fields": ENA_FIELDS, "format": "tsv", "limit": 0}))
    txt = h.fetch(url, "ebi")
    if not txt.strip():
        return []
    rows = list(csv.DictReader(txt.splitlines(), delimiter="\t"))
    if rows:
        RUNS.mkdir(exist_ok=True)
        (RUNS / f"{acc}.tsv").write_text(txt)
    return rows


def summarise_runs(rows):
    if not rows:
        return {"n_runs": 0}

    def cnt(key):
        return "; ".join(f"{k or 'NA'}:{v}" for k, v in Counter(r[key] for r in rows).most_common())
    pub = sorted(r["first_public"] for r in rows if r["first_public"])
    gb = sum(int(r["base_count"] or 0) for r in rows) / 1e9
    return {
        "n_runs": len(rows),
        "organisms": cnt("scientific_name"),
        "strategies": cnt("library_strategy"),
        "sources": cnt("library_source"),
        "layouts": cnt("library_layout"),
        "instruments": cnt("instrument_model"),
        "studies_ena": "; ".join(sorted({f"{r['study_accession']}/{r['secondary_study_accession']}"
                                     for r in rows})),
        "first_public": f"{pub[0]}..{pub[-1]}" if pub else "",
        "gbases": f"{gb:.1f}",
        "example_samples": " | ".join(dict.fromkeys(
            clean(r["sample_title"] or r["experiment_title"]) for r in rows))[:400],
    }


def bioproject(acc):
    txt = h.eutil("esearch", db="bioproject", term=f"{acc}[PRJA]", retmode="json")
    ids = json.loads(txt)["esearchresult"].get("idlist", []) if txt else []
    if not ids:
        return {"found": "no"}
    d = json.loads(h.eutil("esummary", db="bioproject", id=ids[0], retmode="json"))["result"][ids[0]]
    return {"found": "yes", "title": clean(d.get("project_title", "")),
            "date": d.get("registration_date", "")[:10],
            "submitter": clean(d.get("submitter_organization", "")),
            "organism": clean(d.get("organism_name", "")),
            "description": clean(d.get("project_description", ""))[:400]}


def geo(acc):
    txt = h.eutil("esearch", db="gds", term=f"{acc}[ACCN] AND gse[ETYP]", retmode="json")
    ids = json.loads(txt)["esearchresult"].get("idlist", []) if txt else []
    if not ids:
        return {"found": "no"}
    d = json.loads(h.eutil("esummary", db="gds", id=ids[0], retmode="json"))["result"][ids[0]]
    rel = ";".join(x.get("targetobject", "") for x in d.get("extrelations", []))
    return {"found": "yes", "title": clean(d.get("title", "")), "date": d.get("pdat", ""),
            "organism": clean(d.get("taxon", "")), "n_samples": d.get("n_samples", ""),
            "gdstype": clean(d.get("gdstype", "")), "bioproject": d.get("bioproject", ""),
            "sra_link": rel, "pubmed": ";".join(map(str, d.get("pubmedids", []))),
            "description": clean(d.get("summary", ""))[:400]}


def gsa(acc):
    kind = "gsa" if acc.startswith("CRA") else "bioproject"
    url = f"https://ngdc.cncb.ac.cn/{kind}/browse/{acc}"
    txt = h.fetch(url, "ebi")
    t = re.sub(r"<script.*?</script>|<style.*?</style>", " ", txt, flags=re.S)
    t = clean(re.sub(r"<[^>]+>", " ", t))
    if not t or "not found" in t.lower()[:400]:
        return {"found": "no", "url": url}
    m = re.search(r"Title\s*:?\s*(.{10,200}?)\s+(?:BioProject|Description|Release|Submitter|Organism)", t)
    runs = len(set(re.findall(r"\bCRR\d{6,}\b", t)))
    exps = len(set(re.findall(r"\bCRX\d{6,}\b", t)))
    org = re.search(r"Organism\s*:?\s*([A-Z][a-z]+ [a-z]+)", t)
    links = sorted(set(re.findall(r"\b(PRJCA\d+|CRA\d+|SAMC\d+)\b", t)) - {acc})
    return {"found": "yes", "url": url, "title": m.group(1) if m else "",
            "organism": org.group(1) if org else "",
            "gsa_runs_on_page": runs, "gsa_experiments_on_page": exps,
            "linked": "; ".join(links[:12]), "description": t[:400]}


def main():
    cands = list(csv.DictReader(open(BUILD / "accession_candidates.tsv"), delimiter="\t"))
    out = []
    for acc in dict.fromkeys(c["accession"] for c in cands):
        studies = "; ".join(dict.fromkeys(c["study_id"] for c in cands if c["accession"] == acc))
        print(acc, studies, file=sys.stderr)
        rec = {"accession": acc, "studies": studies}
        if acc.startswith(("PRJNA", "PRJEB", "PRJDB")):
            rec["type"] = "BioProject"
            rec.update(bioproject(acc))
            rec.update(summarise_runs(ena_runs(acc)))
        elif acc.startswith(("SRP", "ERP", "DRP", "SRX", "ERX", "DRX")):
            rec["type"] = "SRA study" if acc[2] == "P" else "SRA experiment"
            runs = ena_runs(acc)
            rec.update(summarise_runs(runs))
            rec["found"] = "yes" if runs else "no"
            if runs:
                rec["title"] = clean(runs[0]["study_title"])
        elif acc.startswith("GSE"):
            rec["type"] = "GEO series"
            rec.update(geo(acc))
            bp = rec.get("bioproject")
            if bp:
                rec.update(summarise_runs(ena_runs(bp)))
        elif acc.startswith(("CRA", "PRJCA")):
            rec["type"] = "GSA run set" if acc.startswith("CRA") else "CNCB BioProject"
            rec.update(gsa(acc))
        out.append(rec)
    fields = ["accession", "studies", "type", "found", "title", "date", "organism", "submitter",
              "n_samples", "gdstype", "bioproject", "sra_link", "pubmed", "n_runs", "organisms",
              "strategies", "sources", "layouts", "instruments", "studies_ena", "first_public",
              "gbases", "gsa_runs_on_page", "gsa_experiments_on_page", "linked", "url",
              "example_samples", "description"]
    with open(BUILD / "accession_records.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out)
    print(f"wrote {len(out)} records", file=sys.stderr)


if __name__ == "__main__":
    main()
