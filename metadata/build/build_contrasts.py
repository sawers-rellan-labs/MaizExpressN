#!/usr/bin/env python3
"""Build metadata/contrasts.tsv from metadata/samples.tsv.

Each contrast compares two condition groups of one study (samples sharing
every design field except the replicate). Sign convention (CLAUDE.md):
log2FC = log2(numerator / denominator) with the lower-N group as numerator,
so positive = higher under N limitation; for N-form contrasts the
alternative form is the numerator and nitrate the denominator.

Most studies use the generic rule: inside each stratum of held-constant
factors, every lower-N group is contrasted with every higher-N group. The
STUDIES table says, per study, the contrast class, the factors held
constant and any special rule. Usage: python3 metadata/build/build_contrasts.py
"""
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "metadata" / "contrasts.tsv"
COLS = ["contrast_id", "study_id", "contrast_class", "numerator_group", "denominator_group",
        "held_constant", "N_numerator", "N_denominator", "time_h", "replicated",
        "n_numerator", "n_denominator", "sign_note", "notes",
        # additions beyond the HANDOFF schema
        "numerator_samples", "denominator_samples", "genotype", "tissue_group",
        "tissue_detail", "other_factor_level", "tags", "pool", "batch_term"]
GROUP_KEY = ["study_id", "genotype", "tissue_group", "tissue_detail", "N_label_paper",
             "N_form", "N_amount", "N_unit", "N_rank_in_study", "phase", "time_h",
             "other_factor", "other_factor_level"]
BY_DEFAULT = ["genotype", "tissue_detail", "other_factor_level", "time_h"]


def load_groups():
    groups = defaultdict(list)
    for s in csv.DictReader(open(ROOT / "metadata" / "samples.tsv"), delimiter="\t"):
        if s["include"] == "TRUE":
            groups[tuple(s[k] for k in GROUP_KEY)].append(s)
    out = defaultdict(list)
    for key, samples in groups.items():
        g = dict(zip(GROUP_KEY, key))
        g["samples"] = sorted(samples, key=lambda s: int(s["replicate"]))
        out[g["study_id"]].append(g)
    return out


def n_text(g):
    amt = g["N_amount"]
    unit = "" if amt in ("NR", "") else " " + g["N_unit"]
    return f"{amt}{unit} ({g['N_form']})"


def label(g):
    parts = [g["genotype"], g["tissue_group"], g["N_label_paper"]]
    if g["other_factor_level"]:
        parts.append(g["other_factor_level"])
    return " | ".join(parts)


def slug(*xs):
    s = "_".join(re.sub(r"[^A-Za-z0-9_.+-]+", "", str(x)) for x in xs if x not in ("", None))
    return s.replace("+", "plus")


def fmt_h(x):
    return ("%g" % float(x)) if x not in ("", None) else ""


def id_time(t):
    """Time label for ids: minutes below 1 h, else hours."""
    if t in ("", None):
        return ""
    t = float(t)
    return f"{round(t * 60)}min" if t < 1 else f"{t:g}h"


def id_level(x):
    """First two words of a factor level, for ids."""
    words = re.findall(r"[A-Za-z0-9.]+", x or "")
    return "".join(w[:1].upper() + w[1:] for w in words[:2])


def make(study, cls, num, den, *, time_h="", tags=(), notes="", pool=True, batch="",
         held=None, sid=None, drop=(), sign=None):
    ns = [s["sample_id"] for s in num["samples"] if s["sample_id"] not in drop]
    ds = [s["sample_id"] for s in den["samples"] if s["sample_id"] not in drop]
    held = held or {}
    tags = list(tags)
    if len(ns) < 2 or len(ds) < 2:
        tags.append("unreplicated" if (len(ns) < 2 and len(ds) < 2) else "single_library_arm")
    if sign:
        pass
    elif cls == "form":
        sign = (f"log2FC = log2({num['N_form']} / nitrate); positive = higher with "
                f"{num['N_form']} than with nitrate at the same total N")
    else:
        sign = (f"log2FC = log2({num['N_label_paper']} / {den['N_label_paper']}); "
                "positive = higher under the lower N supply")
    return {
        "contrast_id": sid,
        "study_id": study, "contrast_class": cls,
        "numerator_group": label(num), "denominator_group": label(den),
        "held_constant": "; ".join(f"{k}={v}" for k, v in held.items() if v != ""),
        "N_numerator": n_text(num), "N_denominator": n_text(den),
        "time_h": time_h, "replicated": "TRUE" if len(ns) >= 2 and len(ds) >= 2 else "FALSE",
        "n_numerator": len(ns), "n_denominator": len(ds),
        "sign_note": sign, "notes": notes,
        "numerator_samples": ";".join(ns), "denominator_samples": ";".join(ds),
        "genotype": num["genotype"] if num["genotype"] == den["genotype"] else
        f"{num['genotype']} / {den['genotype']}",
        "tissue_group": num["tissue_group"],
        "tissue_detail": num["tissue_detail"] if num["tissue_detail"] == den["tissue_detail"]
        else f"{num['tissue_detail']} / {den['tissue_detail']}",
        "other_factor_level": num["other_factor_level"] if
        num["other_factor_level"] == den["other_factor_level"] else
        f"{num['other_factor_level']} / {den['other_factor_level']}",
        "tags": ";".join(dict.fromkeys(tags)), "pool": "TRUE" if pool else "FALSE",
        "batch_term": batch,
    }


def tissue_tag(groups):
    """Short tissue label for ids when one tissue_group holds several tissue_details."""
    details = defaultdict(set)
    for g in groups:
        details[g["tissue_group"]].add(g["tissue_detail"])
    tags = {}
    for tg, ds in details.items():
        for d in ds:
            words = re.findall(r"[A-Za-z]+", d)
            tags[d] = tg if len(ds) == 1 else tg + "-" + "".join(w.capitalize() for w in words[:2])
    return tags


def generic(study, groups, cls, by=BY_DEFAULT, keep=lambda g: True, **kw):
    """Every lower-rank group vs every higher-rank group within strata of `by`."""
    ttag = tissue_tag(groups)
    strata = defaultdict(list)
    for g in groups:
        if keep(g):
            strata[tuple(g[b] for b in by)].append(g)
    out = []
    for key, gs in strata.items():
        for a in gs:
            for b in gs:
                if a["N_rank_in_study"] == "NR" or b["N_rank_in_study"] == "NR":
                    continue
                if int(a["N_rank_in_study"]) < int(b["N_rank_in_study"]):
                    held = dict(zip(by, key))
                    held = {("tissue" if k == "tissue_detail" else
                             "other" if k == "other_factor_level" else k): v
                            for k, v in held.items()}
                    t = fmt_h(a["time_h"]) if a["time_h"] == b["time_h"] else ""
                    sid = slug(study, cls, a["genotype"] if "genotype" in by else "",
                               ttag[a["tissue_detail"]], id_level(a["other_factor_level"]) if
                               "other_factor_level" in by else "",
                               id_time(t), id_level(b["N_label_paper"]) if
                               len([x for x in gs if x["N_rank_in_study"] ==
                                    b["N_rank_in_study"]]) > 1 else "")
                    out.append(make(study, cls, a, b, held=held, time_h=t, sid=sid, **kw))
    return out


def pick(groups, **kv):
    res = [g for g in groups if all(re.fullmatch(v, g[k]) for k, v in kv.items())]
    return res


# ---------------------------------------------------------------- per-study rules
def rules():
    R = {}

    def simple(study, cls, by=BY_DEFAULT, **kw):
        R[study] = lambda gs: generic(study, gs, cls, by=by, **kw)

    # steady-state level contrasts
    for st in ["Du2021", "Guo2022", "Mu2017", "Li2022", "Ma2021a"]:
        simple(st, "deficiency")
    simple("Lv2016", "deficiency", notes="rRNA-depleted total RNA libraries.")
    simple("Singh2022", "deficiency",
           notes="N+ concentration not reported (NR); ranks from the paper's N+/N- labels.")
    simple("Bi2014", "deficiency", tags=["batch_confounded"],
           notes="All low-N libraries in one GA-II run (36-nt) and all high-N in another (40-nt): "
                 "N is confounded with sequencing run. Trim reads to 36 nt before quantification.")
    simple("Mager2018", "deficiency",
           notes="-P samples (other_factor low P) have no N partner and are not used.")
    simple("Buzzicotti2026", "deficiency", keep=lambda g: g["other_factor_level"] == "control",
           notes="Water-stress samples got nutrient solution only at transplanting and are not "
                 "used in N contrasts.")
    simple("TorresRodriguez2021", "deficiency",
           notes="Fertigated sand. Within P level only: leaf libraries were sequenced in lane "
                 "batches that coincide with P level.")
    simple("AlonsoNieves2022", "deficiency",
           notes="pho1;2a' mutant N contrasts (wild type is TorresRodriguez2021's, counted there). "
                 "Within P level only (leaf lane batch = P level).")
    simple("Rosman2024", "deficiency", batch="replicate (sequencing run: FAS = rep 1, ADSF = rep 2)",
           tags=["batch_blocked"],
           notes="N supplied as 5:1 nitrate:ammonium. Model ~ batch + group; replicate equals "
                 "sequencing run, balanced across groups.")
    simple("Guo2023", "deficiency",
           notes="Both groups moved from 4 mM N to 0.04 (ND) or 8 mM N (NL) 3 d before sampling; "
                 "the sufficient 4 mM level was not sequenced.")
    simple("Wang2025b", "deficiency",
           notes="Greenhouse, area-based fertilizer for the whole season (100 vs 200 kg N/ha "
                 "equivalent); N form NR.")
    simple("Li2025", "excess",
           notes="HN = 5 mM KNO3, most likely added to the 6.5 mM N control (11.5 mM N); amount "
                 "kept NR until confirmed.")
    # field N rates
    simple("Chen2015", "field_rate", tags=["merged_replicates"],
           notes="Each run merges the paper's 2 biological replicates (8 plants); NN and LN plots "
                 "were in different fields.")
    simple("Ning2023", "field_rate", tags=["low_depth"],
           notes="Runs hold 0.95-10.6 M spots, far below the paper's ~47 M; check library size "
                 "after quantification.")
    simple("Li2023a", "field_rate", tags=["unexplained_sample_sets"],
           notes="SRA holds two sets of 12 libraries (numeric and letter aliases), probably two "
                 "sampling stages; kept as separate contrasts, never pooled as 6 replicates.")
    simple("Pan2015", "field_rate", tags=["deferred_processing", "label_check", "batch_blocked"],
           batch="sequencing run (replicate 3 = 35-nt run, replicates 1-2 = 49-nt run)",
           notes="21-nt DGE tags (processing deferred). GEO treatment field contradicts titles; "
                 "labels follow titles, verify the LN direction after quantification.")

    def cheng(gs):
        drop = {s["sample_id"] for g in gs for s in g["samples"]
                if "Rep3-2" in s["curation_note"] or s["sample_id"].endswith("R3b")}
        return generic("Cheng2021", gs, "field_rate", drop=drop,
                       notes="Field, Urbana IL 2014; 0 vs 200 kg N/ha ammonium sulfate. PH207 lowN "
                             "'Rep3-2' (technical duplicate of Rep3) left out. FR1064 has no "
                             "highN library and no contrast.")
    R["Cheng2021"] = cheng

    # short-term changes with time-matched controls
    simple("He2016", "starvation",
           notes="Transfer at the two-leaf stage from normal N to low N (0.04 mM, not zero); "
                 "normal-N kept as the time-matched control. rtcs background NR.")
    for st, note in [
        ("Cao2024", "2 h nitrate induction vs KCl mock after 4 d N-free."),
        ("Ge2020", "30 min 15 mM KNO3 vs 15 mM KCl after 14 d without N; whole seedlings."),
        ("Huang2025", "QuantSeq 3' tag libraries (UMI-aware processing). +N = 20 mM KNO3 + 20 mM "
                      "NH4NO3 after 2 d N starvation; -N = KCl."),
        ("Buoso2021", "First N supply (24 h) to 5-d-old seedlings germinated without N; no "
                      "nitrate treatment, so no form contrast against nitrate."),
        ("Ravazzolo2021a", "Inhibitor samples (TIS108 at -N only, PCIB at +NO3 only) have no N "
                           "partner and are not used."),
        ("Wang2022", "First N supply after N-free germination (CK 2 mM N vs LN 0), per stage."),
    ]:
        simple(st, "resupply", notes=note)
    simple("Du2020", "resupply",
           notes="15 mM KNO3 vs 15 mM KCl after pre-growth on 0.15 mM nitrate; one pooled sample "
                 "per time point (consider a time-course model across the 8 points).")
    simple("Wu2024", "resupply",
           notes="Mo17 series only; the B73 series is Du2020's (same libraries, counted once). "
                 "One pooled sample per time point.")
    simple("Ravazzolo2020", "resupply",
           notes="24 h after 24 h without N; 1 mM KNO3 or 1 mM NH4Cl. Its -N and +NO3 groups also "
                 "serve as Ravazzolo2021b's controls.")
    R["Ravazzolo2020"] = (lambda base: (lambda gs: base(gs) + ravazzolo_form(gs)))(R["Ravazzolo2020"])

    # time courses without time-matched controls
    def meng(gs):
        return generic("Meng2020", gs, "resupply", by=["genotype", "tissue_detail"],
                       tags=["not_time_matched"],
                       notes="2 h after 0.5 mM KNO3 vs the 0 h baseline (no 2 h mock).")

    def wang2020a(gs):
        return generic("Wang2020a", gs, "resupply", by=["genotype", "tissue_detail"],
                       tags=["not_time_matched"],
                       notes="30 min after nitrate supply vs the untreated 0.15 mM baseline; "
                             "resupply concentration NR.")
    R["Meng2020"], R["Wang2020a"] = meng, wang2020a

    def jia(gs):
        out = []
        for tissue in ("Leaf", "Root"):
            base = pick(gs, tissue_group=tissue, N_label_paper=r"DN0.*")[0]
            dn96 = pick(gs, tissue_group=tissue, N_label_paper="DN96")[0]
            for g in pick(gs, tissue_group=tissue, phase="starvation"):
                t = fmt_h(g["time_h"])
                out.append(make("Jia2023", "starvation", g, base, time_h=t,
                                tags=["not_time_matched"], sid=slug("Jia2023_starvation", tissue, id_time(t)),
                                held={"genotype": "Xu178", "tissue": g["tissue_detail"]},
                                notes="N deprivation (0 N) vs the 4 mM time-0 baseline; only one "
                                      "HN control was sequenced, so time and age are confounded."))
            for g in pick(gs, tissue_group=tissue, phase="resupply"):
                t = fmt_h(float(g["time_h"]) - 96)
                out.append(make("Jia2023", "resupply", dn96, g, time_h=t,
                                tags=["not_time_matched"], sid=slug("Jia2023_resupply", tissue, id_time(t)),
                                held={"genotype": "Xu178", "tissue": g["tissue_detail"]},
                                notes="4 mM nitrate resupplied after 96 h deprivation vs the 96 h "
                                      "deprived sample (time 0 of resupply); time_h = hours since "
                                      "resupply."))
        return out
    R["Jia2023"] = jia

    def braynen(gs):
        out = []
        for tissue in ("Leaf", "Root"):
            nn0 = pick(gs, tissue_group=tissue, phase="steady")[0]
            ln24 = pick(gs, tissue_group=tissue, phase="starvation", time_h="24")[0]
            for g in pick(gs, tissue_group=tissue, phase="starvation"):
                t = fmt_h(g["time_h"])
                out.append(make("Braynen2025", "starvation", g, nn0, time_h=t,
                                tags=["not_time_matched"], sid=slug("Braynen2025_starvation", tissue, id_time(t)),
                                held={"genotype": "B73", "tissue": g["tissue_detail"]},
                                notes="Nitrate removed vs the normal-N time-0 sample (the only "
                                      "normal-N sample); normal-N concentration NR."))
            for g in pick(gs, tissue_group=tissue, phase="resupply"):
                t = fmt_h(float(g["time_h"]) - 24)
                out.append(make("Braynen2025", "resupply", ln24, g, time_h=t,
                                tags=["not_time_matched"], sid=slug("Braynen2025_resupply", tissue, id_time(t)),
                                held={"genotype": "B73", "tissue": g["tissue_detail"]},
                                notes="Nitrate re-supplied after 24 h without N vs the 24 h "
                                      "N-deprived sample; time_h = hours since re-supply."))
        return out
    R["Braynen2025"] = braynen

    def wang2019(gs):
        out = []
        for tissue in ("Root", "Shoot"):
            ck = pick(gs, tissue_group=tissue, N_form="nitrate")[0]
            for g in pick(gs, tissue_group=tissue, N_form="mixed|ammonium"):
                out.append(make("Wang2019", "form", g, ck, sid=slug("Wang2019_form", tissue, g["N_form"]),
                                held={"genotype": "Zhengdan958", "tissue": tissue, "total N": "1 mM N"},
                                notes=f"{g['N_label_paper']} vs sole nitrate (CK), 1 mM total N, "
                                      "hydroponic, 12 d."))
        return out
    R["Wang2019"] = wang2019

    def li2026a(gs):
        out = []
        root = [g for g in gs if g["tissue_group"] == "Root"]
        shoot = [g for g in gs if g["tissue_group"] == "Shoot"]
        held = {"genotype": "B73", "time_h": "8"}
        cn, cc = pick(root, N_label_paper=r"C\.NO3")[0], pick(root, N_label_paper=r"C\.Cl")[0]
        sn, sc = pick(root, N_label_paper=r"Sp\.NO3")[0], pick(root, N_label_paper=r"Sp\.Cl")[0]
        note0 = "2 d N starvation, then 8 h; nitrate = 2 mM Ca(NO3)2 (4 mM N), control CaCl2."
        out.append(make("Li2026a", "resupply", cc, cn, time_h="8", held={**held, "tissue": "Root"},
                        sid="Li2026a_resupply_Root_8h", notes="Homogeneous: both root halves nitrate "
                        "vs both chloride. " + note0))
        out.append(make("Li2026a", "resupply", pick(shoot, N_label_paper=r"C\.Cl")[0],
                        pick(shoot, N_label_paper=r"C\.NO3")[0], time_h="8",
                        held={**held, "tissue": "Shoot"}, sid="Li2026a_resupply_Shoot_8h",
                        notes="Shoots of homogeneous nitrate vs chloride plants. " + note0))
        out.append(make("Li2026a", "local", cc, sn, time_h="8", held={**held, "tissue": "Root"},
                        tags=["local_supply", "split_vs_homogeneous"], sid="Li2026a_local_Root_8h",
                        notes="Root half in a local nitrate patch (split plant) vs homogeneous "
                              "chloride roots. " + note0))
        out.append(make("Li2026a", "local", sn, cn, time_h="8", held={**held, "tissue": "Root"},
                        tags=["local_supply", "systemic_signal", "split_vs_homogeneous"], pool=False,
                        sid="Li2026a_systemic_nitrateSide_Root_8h",
                        sign="log2FC = log2(Sp.NO3 / C.NO3); positive = higher in the split plant "
                             "(same local nitrate, less whole-plant N)",
                        notes="Systemic signal: nitrate-fed half of split plants vs homogeneous "
                              "nitrate (same local N, less whole-plant N). Display only, not pooled."))
        out.append(make("Li2026a", "local", cc, sc, time_h="8", held={**held, "tissue": "Root"},
                        tags=["local_supply", "systemic_signal", "split_vs_homogeneous"], pool=False,
                        sid="Li2026a_systemic_chlorideSide_Root_8h",
                        sign="log2FC = log2(C.Cl / Sp.Cl); positive = higher in the homogeneous "
                             "chloride plant (same local supply, less whole-plant N)",
                        notes="Systemic signal: chloride half of split plants vs homogeneous "
                              "chloride (same local N, more whole-plant N). Display only, not pooled."))
        return out
    R["Li2026a"] = li2026a

    for st, note in [("Yu2015", "Local 4 mM nitrate on one shoot-borne root (24 h) vs homogeneous "
                                "0.5 mM; root stele."),
                     ("Yu2016", "Local 4 mM nitrate (24 h) vs homogeneous 0.5 mM; LCM phloem-pole "
                                "pericycle; compare within root type only.")]:
        simple(st, "local", tags=["local_supply"], notes=note)
    return R


def ravazzolo_form(gs):
    nh4 = pick(gs, N_form="ammonium")[0]
    no3 = pick(gs, N_form="nitrate")[0]
    return [make("Ravazzolo2020", "form", nh4, no3, time_h="24", sid="Ravazzolo2020_form_Root_24h_ammonium",
                 held={"genotype": "B73", "tissue": nh4["tissue_detail"], "total N": "1 mM N"},
                 notes="1 mM NH4Cl vs 1 mM KNO3, 24 h after 24 h without N.")]


NO_CONTRASTS = {
    "Ravazzolo2021b": "Only NO3- + cPTIO samples are new; its -N and +NO3 controls are "
                      "Ravazzolo2020's, so it adds no N contrast (cPTIO effects are not N contrasts).",
}


def main():
    groups = load_groups()
    R = rules()
    contrasts = []
    missing = []
    for study in sorted(groups):
        if study in NO_CONTRASTS:
            continue
        if study not in R:
            missing.append(study)
            continue
        contrasts += R[study](groups[study])
    if missing:
        sys.exit(f"no contrast rule for: {missing}")
    ids = [c["contrast_id"] for c in contrasts]
    dup = {i for i in ids if ids.count(i) > 1}
    if dup:
        sys.exit(f"duplicate contrast ids: {sorted(dup)}")
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        w.writerows(contrasts)
    print(f"{len(contrasts)} contrasts from {len({c['study_id'] for c in contrasts})} studies "
          f"-> {OUT}", file=sys.stderr)
    for s, why in NO_CONTRASTS.items():
        print(f"  {s}: no contrasts ({why})", file=sys.stderr)


if __name__ == "__main__":
    main()
