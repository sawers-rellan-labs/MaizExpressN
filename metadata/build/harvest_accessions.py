#!/usr/bin/env python3
"""Harvest candidate data accessions for each catalogue study.

Sources, per study:
  1. NCBI elink from the PMID to BioProject, SRA and GEO DataSets.
  2. Europe PMC text-mined accession annotations (PMC and MED ids).
  3. Europe PMC database cross-references for the article.
  4. Accession-pattern scan of full text (NCBI PMC efetch and Europe PMC XML).

Writes metadata/build/accession_evidence.tsv (one row per study x accession x
source). Responses are cached under metadata/build/cache/ so reruns are cheap.
Standard library only.

Usage: python3 metadata/build/harvest_accessions.py
"""
import csv
import hashlib
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CATALOGUE = ROOT / "catalogue" / "maize_N_rnaseq_studies.tsv"
BUILD = ROOT / "metadata" / "build"
CACHE = BUILD / "cache"
OUT = BUILD / "accession_evidence.tsv"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
EPMC = "https://www.ebi.ac.uk/europepmc/"

ACC_RE = re.compile(
    r"\b(PRJ(?:NA|EB|DB|CA)\d{3,}|[SED]R[PXRS]\d{5,}|GSE\d{3,}|GSM\d{3,}|"
    r"CRA\d{5,}|CRR\d{5,}|CRX\d{5,}|SAM[NED][A-Z]?\d{5,}|SAMC\d{5,}|"
    r"E-MTAB-\d{3,}|E-GEOD-\d{3,})\b"
)

_last_call = {"ncbi": 0.0, "ebi": 0.0}


def fetch(url, kind):
    """GET with on-disk cache and polite rate limiting."""
    CACHE.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha1(url.encode()).hexdigest()
    path = CACHE / key
    if path.exists():
        return path.read_text(errors="replace")
    wait = 0.4 if kind == "ncbi" else 0.2
    dt = time.time() - _last_call[kind]
    if dt < wait:
        time.sleep(wait - dt)
    req = urllib.request.Request(url, headers={"User-Agent": "MaizExpressN-accessions/0.1"})
    text = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                text = r.read().decode("utf-8", errors="replace")
            break
        except urllib.error.HTTPError as e:
            if e.code in (404, 400):
                text = ""
                break
            time.sleep(2 ** attempt)
        except Exception:
            time.sleep(2 ** attempt)
    _last_call[kind] = time.time()
    if text is None:
        print(f"  failed: {url}", file=sys.stderr)
        return ""
    path.write_text(text)
    (CACHE / (key + ".url")).write_text(url)
    return text


def eutil(name, **params):
    params.setdefault("tool", "MaizExpressN")
    return fetch(EUTILS + name + ".fcgi?" + urllib.parse.urlencode(params), "ncbi")


def elink_pubmed(pmid, db):
    """Return linked UIDs in db for a PMID."""
    txt = eutil("elink", dbfrom="pubmed", db=db, id=pmid, retmode="json")
    if not txt:
        return []
    uids = []
    for ls in json.loads(txt).get("linksets", []):
        for lsdb in ls.get("linksetdbs", []):
            uids += lsdb.get("links", [])
    return sorted(set(uids))


def uid_to_acc(db, uids):
    """Convert BioProject / SRA / GDS UIDs to accessions via esummary."""
    if not uids:
        return []
    txt = eutil("esummary", db=db, id=",".join(uids), retmode="json")
    if not txt:
        return []
    res = json.loads(txt).get("result", {})
    out = []
    for uid in res.get("uids", []):
        d = res[uid]
        if db == "bioproject":
            out.append((d.get("project_acc", uid), d.get("project_title", "")))
        elif db == "gds":
            out.append((d.get("accession", uid), d.get("title", "")))
        elif db == "sra":
            expxml = d.get("expxml", "")
            m = re.search(r'Study acc="([^"]+)"', expxml)
            t = re.search(r"<Title>(.*?)</Title>", expxml)
            out.append((m.group(1) if m else uid, t.group(1) if t else ""))
    return out


def epmc_annotations(article_id):
    url = (EPMC + "annotations_api/annotationsByArticleIds?"
           + urllib.parse.urlencode({"articleIds": article_id,
                                     "type": "Accession Numbers", "format": "JSON"}))
    txt = fetch(url, "ebi")
    if not txt:
        return []
    out = []
    for art in json.loads(txt):
        for a in art.get("annotations", []):
            exact = a.get("exact", "")
            prefix = (a.get("prefix") or "")[-60:]
            postfix = (a.get("postfix") or "")[:60]
            section = a.get("section", "")
            out.append((exact, f"[{section}] ...{prefix}{exact}{postfix}..."))
    return out


def epmc_dblinks(src, ext_id):
    url = EPMC + f"webservices/rest/{src}/{ext_id}/databaseLinks?format=json&pageSize=100"
    txt = fetch(url, "ebi")
    if not txt:
        return []
    out = []
    try:
        data = json.loads(txt)
    except json.JSONDecodeError:
        return []
    for xr in data.get("dbCrossReferenceList", {}).get("dbCrossReference", []):
        for info in xr.get("dbCrossReferenceInfo", []):
            out.append((info.get("info1", ""), xr.get("dbName", "")))
    return out


def text_from_xml(xml):
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        return re.sub(r"<[^>]+>", " ", xml)
    parts = []
    for el in root.iter():
        if el.text:
            parts.append(el.text)
        # keep hyperlink targets: accession numbers are often only in xlink:href
        for k, v in el.attrib.items():
            if k.endswith("href"):
                parts.append(v)
        if el.tail:
            parts.append(el.tail)
    return " ".join(parts)


def scan_text(text):
    out = {}
    for m in ACC_RE.finditer(text):
        acc = m.group(1)
        if acc not in out:
            s = max(0, m.start() - 80)
            out[acc] = re.sub(r"\s+", " ", text[s:m.end() + 60])
    return out


def main():
    studies = list(csv.DictReader(open(CATALOGUE), delimiter="\t"))
    rows = []

    def add(study, source, acc, context, url=""):
        rows.append({"study_id": study, "source": source, "accession": acc,
                     "context": context.replace("\t", " ").replace("\n", " ")[:400],
                     "url": url})

    for s in studies:
        sid, pmid, pmcid = s["study_id"], s["pmid"], s["pmcid"]
        pmid = pmid if pmid.isdigit() else ""
        pmcid = pmcid if pmcid.startswith("PMC") else ""
        print(sid, pmid, pmcid, file=sys.stderr)

        if pmid:
            for db in ("bioproject", "sra", "gds"):
                for acc, title in uid_to_acc(db, elink_pubmed(pmid, db)):
                    add(sid, f"elink_pubmed_{db}", acc, title,
                        f"https://www.ncbi.nlm.nih.gov/{db}/?term={acc}")
            for acc, ctx in epmc_annotations(f"MED:{pmid}"):
                add(sid, "epmc_annotation_MED", acc, ctx)
            for acc, db in epmc_dblinks("MED", pmid):
                add(sid, "epmc_dblink", acc, db)
        if pmcid:
            for acc, ctx in epmc_annotations(f"PMC:{pmcid}"):
                add(sid, "epmc_annotation_PMC", acc, ctx)
            xml = eutil("efetch", db="pmc", id=pmcid.replace("PMC", ""), retmode="xml")
            for acc, ctx in scan_text(text_from_xml(xml)).items():
                add(sid, "pmc_fulltext_scan", acc, ctx,
                    f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/")
            xml = fetch(EPMC + f"webservices/rest/{pmcid}/fullTextXML", "ebi")
            for acc, ctx in scan_text(text_from_xml(xml)).items():
                add(sid, "epmc_fulltext_scan", acc, ctx,
                    f"https://europepmc.org/article/PMC/{pmcid}")

    # collapse duplicates (same study, source, accession)
    seen, uniq = set(), []
    for r in rows:
        k = (r["study_id"], r["source"], r["accession"])
        if k not in seen:
            seen.add(k)
            uniq.append(r)
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["study_id", "source", "accession", "context", "url"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(uniq)
    print(f"wrote {len(uniq)} rows to {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
