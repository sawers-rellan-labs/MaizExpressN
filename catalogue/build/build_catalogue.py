#!/usr/bin/env python3
"""Merge agent extraction tables into the maize N RNA-seq catalogue.

Inputs : agent_out/pubmed_batch_*.tsv, agent_out/web_literature.tsv
Outputs: out/maize_N_rnaseq_studies.tsv, out/maize_N_related.tsv,
         out/maize_N_to_check.tsv, out/catalogue.json
"""
import csv, glob, json, re, sys, os

BASE = sys.argv[1] if len(sys.argv) > 1 else '.'
OUT = os.path.join(BASE, 'out')
os.makedirs(OUT, exist_ok=True)

rows = []
for f in sorted(glob.glob(os.path.join(BASE, 'agent_out/pubmed_batch_*.tsv'))) + [os.path.join(BASE, 'agent_out/web_literature.tsv')]:
    with open(f, encoding='utf-8') as fh:
        rows += list(csv.DictReader(fh, delimiter='\t', quoting=csv.QUOTE_NONE))

def key(d):
    return d['pmid'] if d['pmid'] not in ('', 'none', 'NR') else d['doi_url']

# ---- manual curation (decisions made while reviewing agent reports) ----
OV = {
 '41803685': dict(study_id='Li2026a'),
 '42210095': dict(study_id='Li2026b', fit='borderline', fit_note='All sequenced samples got 180 kg N/ha; the contrast is urea vs controlled-release urea.'),
 '37891910': dict(study_id='Li2023a'),
 '38005765': dict(study_id='Li2023b', acc_status='on request', acc_note='Data available on request (no deposit).'),
 '32393171': dict(study_id='Wang2020a'),
 '32155833': dict(study_id='Wang2020b'),
 '33579187': dict(study_id='Ma2021a', system_group='Hydroponic', acc_note='Two BioProjects (lncRNA-seq and small RNA-seq); paper does not say which is which.'),
 '33932084': dict(study_id='Ma2021b', related_to='Same group and design as Ma2021a (P178, HN vs LN); may reuse those libraries.'),
 '40941853': dict(study_id='Wang2025a', system_group='Solid substrate'),
 '39757743': dict(study_id='Wang2025b', tissue_group=['Ear/kernel']),
 '33508105': dict(study_id='Ravazzolo2021a', maizegdb=True, details_level='full text',acc_note='Cited with a GEO reviewer token in the accepted manuscript; confirm it is public.'),
 '34502437': dict(study_id='Ravazzolo2021b', leads='GSE173102', acc_status='lead', related_to='The -N and +NO3- libraries are probably the Ravazzolo2020 ones (GSE141860); the +cPTIO libraries are new.'),
 '31968691': dict(acc_status='repository match', maizegdb=True, acc_note='GEO series title matches the paper.'),
 '26990640': dict(acc_status='repository match', maizegdb=True, acc_note='From MaizeGDB run metadata (B73 roots); check that the rtcs runs are in the same BioProject.'),
 '29206208': dict(maizegdb=True),
 '37786281': dict(accession='GSE227472', acc_status='repository match', acc_note='GEO series title matches the paper; paper states GEO deposit.'),
 '34090337': dict(accession='SRP287300; PRJNA669055', related_to='Same BioProject as AlonsoNieves2022 (wild-type libraries shared).'),
 '35844781': dict(accession='SRP287300; PRJNA669055', related_to='Same BioProject as TorresRodriguez2021.'),
 '26136266': dict(leads='GSE59613; PRJNA255755', acc_status='lead', acc_note='Paper states GEO deposit; these IDs came up in searches but were not confirmed.'),
 '26198256': dict(tissue_group=['Root'], leads='SRX1020537; PRJNA283053', acc_status='lead', acc_note='SRA experiment titled "Homogeneous low nitrate, control, biol. Repl.2" matches the design; not confirmed.'),
 '34561450': dict(tissue_group=['Leaf'], acc_note='Paper says GEO; the number was not captured.', related_to='Coruzzi lab (NYU) with Moose lab field trials (UIUC); see Huang2025.'),
 '40365911': dict(tissue_group=['Root', 'Shoot'], acc_note='Paper lists accessions in its supplementary tables.', related_to='Same lab as Cheng2021; NUE models trained on Cheng2021 field data.'),
 '33793942': dict(fit='borderline', fit_note='N was removed together with P, K, Mg and S (multi-nutrient deprivation).', acc_note='BioProject stated; the number was not captured.'),
 '35860530': dict(fit='borderline', fit_note='Observational multi-site field survey; N was not experimentally varied.', acc_note='NCBI deposit stated (251 root mRNA-seq samples); the number was not captured.'),
 '40082800': dict(fit='borderline', fit_note='Fertilizer regimes; N is confounded with P and organic inputs.'),
 '40891006': dict(fit='borderline', fit_note='Dual RNA-seq of U. maydis galls; maize genes were not the focus.'),
 '39395686': dict(fit='tentative', fit_note='Abstract only; RNA-seq method and N design not confirmed.'),
 '30478358': dict(fit='tentative', fit_note='Abstract only; transcriptome method and N design not confirmed.'),
 '40305906': dict(fit='tentative', fit_note='Abstract only; which N x inoculant treatments were sequenced is unclear.'),
 '35736701': dict(related_to='Same HN/LN protocol and lab as Mu2017 (China Agricultural University).'),
 '37629121': dict(system_group='Hydroponic'),
 '42424137': dict(tissue_group=['Root']),
 'https://doi.org/10.1016/S2095-3119(20)63220-8': dict(tissue_group=['Root', 'Leaf']),
 '38545384': dict(acc_note='Data statement points to the article/supplement; raw reads not located.'),
 'https://doi.org/10.3390/agronomy14071541': dict(acc_status='not deposited', acc_note='Data statement: data are in the article and supplementary files.'),
}

TOK = re.compile(r'\b(PRJNA\d+|PRJEB\d+|PRJDB\d+|PRJCA\d+|SRP\d+|SRX\d+|SRR\d+|ERP\d+|GSE\d+|CRA\d+|E-MTAB-\d+)\b')

def acc_url(a):
    if a.startswith(('PRJNA', 'PRJEB', 'PRJDB')):
        return f'https://www.ncbi.nlm.nih.gov/bioproject/{a}'
    if a.startswith(('SRP', 'SRX', 'SRR', 'ERP')):
        return f'https://www.ncbi.nlm.nih.gov/sra/{a}'
    if a.startswith('GSE'):
        return f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={a}'
    if a.startswith('CRA'):
        return f'https://ngdc.cncb.ac.cn/gsa/browse/{a}'
    if a.startswith('PRJCA'):
        return f'https://ngdc.cncb.ac.cn/bioproject/browse/{a}'
    if a.startswith('E-MTAB'):
        return f'https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{a}'
    return ''

def tissue_group(t):
    t = t.lower(); g = []
    if re.search(r'root|pericycle|stele', t): g.append('Root')
    if re.search(r'leaf|leaves', t): g.append('Leaf')
    if re.search(r'shoot', t): g.append('Shoot')
    if re.search(r'whole seedling', t): g.append('Whole seedling')
    if re.search(r'ear|cob|floret|kernel|endosperm|gall', t) and 'ear leaf' not in t and 'leaf above the primary ear' not in t: g.append('Ear/kernel')
    if 'gall' in t: g = ['Leaf gall']
    return g or ['NR']

def system_group(s):
    s = s.lower()
    if 'field' in s: return 'Field'
    if re.search(r'sand|coir|coconut|vermiculite|perlite', s): return 'Solid substrate'
    if re.search(r'hydropon|nutrient solution|aerated|ms hydroponic', s): return 'Hydroponic'
    if re.search(r'pot|soil|greenhouse', s): return 'Soil/pots'
    return 'NR'

out = []
seen = set()
for d in rows:
    k = key(d)
    if k in seen:
        continue
    seen.add(k)
    o = dict(d)
    ov = OV.get(k, {})
    for f in ('study_id', 'accession'):
        if f in ov: o[f] = ov[f]
    toks = TOK.findall(o['accession'])
    lead_toks = TOK.findall(ov.get('leads', ''))
    if 'acc_status' in ov:
        st = ov['acc_status']
    elif toks:
        st = 'confirmed'
    elif 'request' in o['accession'].lower():
        st = 'on request'
    else:
        st = 'to resolve'
    note = ov.get('acc_note', '')
    if st == 'to resolve' and not note:
        nl = o['notes'].lower()
        if re.search(r'stated|deposit|states', nl) and re.search(r'stripp|missing|not captured|not recoverable|not legible', nl):
            note = 'Deposit stated in the paper; the number was not captured.'
        else:
            note = 'No data statement captured; look up by title/authors.'
    src = o['source_of_details'].lower()
    details = ov.get('details_level') or ('abstract only' if (src.startswith('abstract') and 'pmc' not in src) or 'abstract only' in src else 'full text')
    o.update(
        study_id=o['study_id'],
        fit=ov.get('fit', 'core'),
        fit_note=ov.get('fit_note', ''),
        details_level=details,
        design_tokens=[t.strip() for t in o['design_category'].split('+') if t.strip()],
        tissue_group=ov.get('tissue_group') or tissue_group(o['tissue']),
        system_group=ov.get('system_group') or system_group(o['growth_system']),
        accessions=[{'id': a, 'url': acc_url(a)} for a in dict.fromkeys(toks)],
        accession_leads=[{'id': a, 'url': acc_url(a)} for a in lead_toks],
        accession_status=st,
        accession_note=note,
        maizegdb=bool(ov.get('maizegdb')),
        related_to=ov.get('related_to', ''),
    )
    out.append(o)

# keep only core studies (decision 16 Oct 2026); borderline/tentative rows are dropped
dropped = [o for o in out if o['fit'] != 'core']
out = [o for o in out if o['fit'] == 'core']
print('dropped (non-core):', ', '.join(f"{o['study_id']} [{o['fit']}]" for o in dropped))

ids = [o['study_id'] for o in out]
dups = {i for i in ids if ids.count(i) > 1}
if dups:
    sys.exit(f'duplicate study_id: {dups}')

out.sort(key=lambda o: (-int(o['year']) if o['year'].isdigit() else 0, o['study_id']))

cols = ['study_id', 'details_level', 'year', 'first_author', 'journal', 'title', 'doi_url', 'pmid', 'pmcid',
        'design_category', 'genotypes', 'n_genotypes', 'tissue', 'tissue_group', 'stage_age', 'growth_system', 'system_group',
        'N_treatments', 'other_factors', 'timepoints_duration', 'bio_reps', 'n_libraries', 'platform_layout',
        'accession', 'accession_status', 'accession_leads', 'accession_note', 'repository', 'maizegdb', 'related_to',
        'other_omics', 'key_finding', 'notes', 'source_of_details']
with open(os.path.join(OUT, 'maize_N_rnaseq_studies.tsv'), 'w', encoding='utf-8', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(cols)
    for o in out:
        r = []
        for c in cols:
            v = o.get(c, '')
            if c == 'accession':
                v = '; '.join(a['id'] for a in o['accessions']) or o['accession']
            elif c == 'accession_leads':
                v = '; '.join(a['id'] for a in o['accession_leads'])
            elif c == 'tissue_group':
                v = '; '.join(v)
            elif c == 'maizegdb':
                v = 'yes' if v else ''
            r.append(str(v).replace('\t', ' ').replace('\n', ' '))
        w.writerow(r)

# ---- related datasets and items still to check (from agent reports + abstract checks) ----
RELATED = [
 # type, study, pmid, doi, design, accession, note
 ('Microarray', 'Schlüter2012', '22972706', '10.1104/pp.112.204420', 'A188 and B73 source leaves; 15 vs 0.15 mM nitrate', 'GSE40678', 'Series title matches design'),
 ('Microarray', 'Schlüter2013', '23822863', '10.1186/1471-2164-14-442', 'B73 source leaves; low N, low P and cold', 'GEO (number not captured)', 'N arrays shared with Schlüter2012'),
 ('Microarray', 'Yang2011', '21980173', '10.1104/pp.111.187898', 'N-status biomarker genes; limiting vs sufficient N, recovery, field', 'GSE32361', ''),
 ('Microarray', 'Coll2010', '20349115', '10.1007/s11103-010-9624-5', 'MON810 vs non-GM varieties; conventional vs low-N field', '', 'Affymetrix'),
 ('Microarray', 'Humbert2013', '23324127', '10.1186/1471-2164-14-3', 'Water x N stress combinations; SRG-200; several tissues', '', 'Custom Affymetrix'),
 ('Microarray', 'Zamboni2014', '24805158', '10.1111/jipb.12214', 'Nitrate induction time course in roots; Lo5 (high NUE) vs T250', '', ''),
 ('Microarray', 'Mascia2019', '31574923', '10.3390/ijms20194856', 'Lo5 vs T250 roots; 1 and 4 d N starvation', 'GEO (number not captured)', 'Agilent 4x44K'),
 ('Microarray', 'Plett2016a', '26038196', '10.1111/pbi.12388', 'Gaspé Flint leaf and root; 0.5 vs 2.5 mM nitrate across 7 life-cycle stages', 'GSE61820', ''),
 ('Microarray', 'Plett2016b', '27511191', '10.1007/s11103-016-0512-5', 'Same plants as Plett2016a; N-assimilation enzymes vs transcripts', 'GSE61820', 'Reuses Plett2016a arrays'),
 ('Microarray', 'Amiour2012', '22936829', '10.1093/jxb/ers186', 'Leaves under long-term N deficiency at two stages; + proteome, metabolome', '', ''),
 ('Microarray', 'Zanin2014', '25524070', '10.1093/pcp/pcu202', 'Roots: urea, nitrate, urea + nitrate vs no N', '', ''),
 ('Microarray', 'Zanin2016', '27446099', '10.3389/fpls.2016.00845', 'Urea-fed seedlings +/- urease inhibitor NBPT', '', 'No N-level contrast'),
 ('Microarray', 'Chowdhury2022', '34554248', '10.1093/jxb/erab435', 'B73 root metabolic model; 5 vs 0.1 mM nitrate', '', 'cDNA arrays; data in Table S2'),
 ('cDNA-AFLP', 'Trevisan2011', '21762167', '10.1111/j.1469-8137.2011.03822.x', 'Seedling roots, nitrate supply; hemoglobin-NR co-regulation', '', ''),
 ('Small RNA', 'Zhao2012', '22235323', '10.1371/journal.pone.0029669', 'B73 shoot and root miRNAs after 2 d N deprivation; + degradome', '', ''),
 ('Small RNA', 'Xu2011', '22132192', '10.1371/journal.pone.0028009', 'Ye478 leaf and root miRNAs; chronic and 2-h low nitrate (miRNA arrays)', 'GEO (number not captured)', ''),
 ('Small RNA', 'Zhao2013', '23788746', '10.1093/aob/mct133', 'Seedling miRNAs and targets under low nitrate; + degradome', '', ''),
 ('Small RNA', 'Trevisan2012', '22211437', '10.1111/j.1365-3040.2011.02478.x', 'Roots +/- nitrate; miRNA arrays + in situ', '', ''),
 ('Small RNA', 'Yang2019', '31313105', '10.1007/s13258-019-00848-0', 'miRNAs under N deficiency vs resupply vs control', '', 'sRNA-seq'),
 ('Reprocessed / meta-analysis', 'Hayford2024', '38816789', '10.1186/s12864-024-10443-7', 'MaizeGDB 2024: 24 B73 stress datasets remapped to B73 v5, incl. He2016, Mu2017, Ravazzolo2020, Ravazzolo2021a', '', 'Quick lookups in MaizeGDB/qTeller'),
 ('Reprocessed / meta-analysis', 'Pardo2025', '41092016', '10.1093/g3journal/jkaf223', 'Core abiotic-stress genes from 39 BioProjects (5 low-N)', 'PRJNA904734 (low-N project named in text)', 'Full list in supplement; preprint PMID 40027706'),
 ('Reprocessed / meta-analysis', 'Jiang2018', '29518046', '10.3390/genes9030151', 'Gene networks from public N microarrays (90 Affymetrix arrays; Schlüter set)', '', ''),
 ('Reprocessed / meta-analysis', 'Xue2024', '39529933', '10.3389/fpls.2024.1477383', 'ZmCASPL family; low-N time course data', 'PRJNA304223 (He2016)', ''),
 ('Reprocessed / meta-analysis', 'Ma2017', '28764653', '10.1186/s12870-017-1077-4', 'Co-expression network from 787 SRA runs (36 studies); nitrate-assimilation module', '', 'Study list in supplement'),
 ('Reprocessed / meta-analysis', 'Luo2022', '36388550', '10.3389/fpls.2022.1013598', 'WGCNA + meta-QTL across low-N, low-P and disease datasets', '', 'Dataset list in supplement'),
 ('Reprocessed / meta-analysis', 'Yang2026', '42452185', '10.3390/plants15131977', 'Stress co-expression modules from 15 datasets (286 samples), incl. low N/P', 'GSE142477 (validation set)', ''),
 ('Reprocessed / meta-analysis', 'Mia2026', '41746965', '10.1371/journal.pone.0342990', 'Sugar transporters; reused public HN/LN hybrid + parents RNA-seq', '', 'Probably Bi2014 data'),
 ('Reprocessed / meta-analysis', 'Liu2026', '41977721', '10.3390/plants15071062', 'AHL family; reused group\'s V33/V37 drought and low-N RNA-seq', '', ''),
 ('Reprocessed / meta-analysis', 'Yamamoto2022', '35958195', '10.3389/fpls.2022.905894', 'PEPC family; reanalysed maize/sorghum N-stress RNA-seq', '', ''),
 ('Reprocessed / meta-analysis', 'Xing2023', '37967146', '10.1111/tpj.16550', 'NUE candidate genes mined from published omics', '', ''),
 ('Reprocessed / meta-analysis', 'Yan2025', '41514965', '10.3390/plants15010020', 'Low-N root aerenchyma GWAS; reused published 4F1 vs Zheng32 transcriptome', '', ''),
 ('Reprocessed / meta-analysis', 'Lin2014', '25209110', '10.1007/s11103-014-0246-1', 'Reference genes stable across genotypes and N treatments', '', 'N tested by qPCR only'),
 ('Reprocessed / meta-analysis', 'Simons2014', '25248718', '10.1104/pp.114.245787', 'Leaf genome-scale model; N-limited vs N-complete; gln1 mutants', '', ''),
 ('RNA-seq, no N-supply contrast', 'Gong2020', '32392584', '10.1093/jxb/eraa229', 'N remobilization after silking; PH4CV vs PH6WC leaves; 15N', '', 'Useful for remobilization, not N levels'),
 ('RNA-seq, no N-supply contrast', 'Yang2023', '37922646', '10.1016/j.plaphy.2023.108147', 'Leaf N fluxes in stay-green hybrids XY335 vs XY696; 15N', '', ''),
 ('RNA-seq, no N-supply contrast', 'Sun2018', '29597009', '10.1016/j.molp.2018.03.013', 'ZmmiR528 transgenics under N luxury; lignin and lodging', '', 'N contrast in RNA-seq not confirmed'),
]
TOCHECK = [
 # kind, label, pmid, link, why
 ('Paper', 'Li2015 (Plant Cell)', '25944100', 'https://doi.org/10.1105/tpc.15.00158', 'atg12 autophagy mutants under N starvation; MeSH lists RNA sequencing'),
 ('Paper', 'Zheng2025 (J Exp Bot)', '40365707', 'https://doi.org/10.1093/jxb/eraf202', 'ZmNRL1 mutant/overexpressor; comparative transcriptome under N stress'),
 ('Paper', 'Yu2021 (Nat Plants)', '33833418', 'https://doi.org/10.1038/s41477-021-00897-y', 'Root-axis transcriptome; N-deprivation experiments with lrt1 and flavones'),
 ('Paper', 'Wang2024 (Plant Cell)', '38917216', 'https://doi.org/10.1093/plcell/koae184', 'ZmNLP3.2 knockout/overexpressor under low N'),
 ('Paper', 'Gu2025 (Plant J)', '40474668', 'https://doi.org/10.1111/tpj.70263', 'ZmNLP8 mutants under different nitrate conditions'),
 ('Paper', 'Wang2026 (J Integr Plant Biol)', '41566850', 'https://doi.org/10.1111/jipb.70147', 'JA / ZmbHLH99 / ZmNLP3.2 under low nitrate'),
 ('Paper', 'Li2026c (Plant Cell Environ)', '41555148', 'https://doi.org/10.1111/pce.70371', 'ZmPT7 mutants sensitive to low nitrate'),
 ('Paper', 'Hu2025 (New Phytol)', '39888029', 'https://doi.org/10.1111/nph.20422', 'ZmNPF7.10; K and N distribution'),
 ('Paper', 'Wang2020c (Plant Cell Environ)', '33215716', 'https://doi.org/10.1111/pce.13956', 'ZmAER overexpression; low-N tolerance'),
 ('Paper', 'Huang2022 (Nature)', '36385527', 'https://doi.org/10.1038/s41586-022-05441-2', 'THP9 teosinte allele; NUE under low N'),
 ('Paper', 'Hui2022 (Plant Cell)', '35880836', 'https://doi.org/10.1093/plcell/koac225', 'ZmAMT3;1 RNAi; mycorrhiza-dependent N uptake'),
 ('Paper', 'Pii2019 (J Exp Bot)', '30715422', 'https://doi.org/10.1093/jxb/ery433', 'Azospirillum inoculation x nitrate induction'),
 ('Paper', 'Hsu2025 (Plant J)', '40616845', 'https://doi.org/10.1111/tpj.70319', 'Rhizosphere N dynamics in Andropogoneae grasses; RNA-seq of roots and leaves deposited on USDA Ag Data Commons'),
 ('Paper', 'Liang2024 (Physiol Plant)', '38472657', 'https://doi.org/10.1111/ppl.14249', 'Fulvic acid under low N; DEGs (N-level contrast unclear)'),
 ('Paper', 'Li2021 (J Sci Food Agric)', '34888887', 'https://doi.org/10.1002/jsfa.11711', 'Controlled-release urea + fulvic acid; leaf transcriptome'),
 ('Paper', 'Plant and Soil (N/P/K deficiency, roots)', '', 'https://doi.org/10.1007/s11104-019-04385-3', 'Transcriptome x N in title; not read'),
 ('Paper', 'J Soil Sci Plant Nutr (ear leaves, long-term N)', '', 'https://doi.org/10.1007/s42729-021-00637-3', 'Transcriptome x N in title; not read'),
 ('Paper', 'Crop J 2023 (water x N)', '', 'https://www.sciencedirect.com/science/article/pii/S2214514123001162', 'Transcriptome x N in title; not read'),
 ('Paper', 'Crop J', '', 'https://www.sciencedirect.com/science/article/pii/S2214514124000990', 'RNA-seq status unknown'),
 ('Paper', 'Crop J', '', 'https://www.sciencedirect.com/science/article/pii/S2214514124001855', 'RNA-seq status unknown'),
 ('Paper', 'Crop J', '', 'https://www.sciencedirect.com/science/article/pii/S2214514125002028', 'RNA-seq status unknown'),
 ('Paper', 'Plant Stress', '', 'https://www.sciencedirect.com/science/article/pii/S2667064X25001617', 'RNA-seq status unknown'),
 ('Paper', 'Environ Exp Bot', '', 'https://www.sciencedirect.com/science/article/pii/S0098847221003622', 'RNA-seq status unknown'),
 ('Paper', 'Appl Soil Ecol', '', 'https://www.sciencedirect.com/science/article/pii/S0929139325000502', 'RNA-seq status unknown'),
 ('Paper', 'Plant and Soil 2024', '', 'https://doi.org/10.1007/s11104-024-06561-6', 'RNA-seq status unknown'),
 ('Paper', 'J Plant Biol 2024', '', 'https://doi.org/10.1007/s12374-024-09428-6', 'RNA-seq status unknown'),
 ('Preprint', 'Xu 2022 (bioRxiv)', '', 'https://doi.org/10.1101/2022.02.09.479784', 'Unpublished preprint; RNA-seq status unknown'),
 ('Paper', 'Scientia Agricultura Sinica 2023 (Chinese)', '', 'https://doi.org/10.3864/j.issn.0578-1752.2023.06.002', 'Title not retrieved'),
 ('Paper', 'Scientia Agricultura Sinica 2022 (Chinese)', '', 'https://doi.org/10.3864/j.issn.0578-1752.2022.03.002', 'Title not retrieved'),
 ('Paper', 'Acta Agronomica Sinica 2023 (Chinese)', '', 'https://doi.org/10.3724/SP.J.1006.2023.23045', 'Title not retrieved'),
 ('Dataset', 'pepck1::Ds vs W22 leaf gradient x N (USDA Ag Data Commons)', '', 'https://agdatacommons.nal.usda.gov/articles/dataset/Transcriptome_of_maize_leaf_gradient_and_nitrogen_response_in_pepck1_Ds_mutant_and_W22_control_plants/25089188/1', 'No linked paper found; BioProject not retrieved'),
 ('Dataset', 'PRJNA904734', '', 'https://www.ncbi.nlm.nih.gov/bioproject/PRJNA904734', 'Low-N BioProject named in Pardo2025; source study unknown'),
]
for a in ['E-MTAB-8633', 'GSE111425', 'GSE49681', 'GSE135613', 'PRJNA707066', 'PRJNA299645', 'PRJCA044958', 'PRJCA010624',
          'GSE17932', 'GSE125417', 'PRJNA1116564', 'PRJNA208608', 'E-MTAB-4246', 'E-MTAB-4208']:
    TOCHECK.append(('Dataset', a, '', acc_url(a), 'Came up in N-related repository searches; title not seen, may not be N'))

with open(os.path.join(OUT, 'maize_N_related.tsv'), 'w', encoding='utf-8', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['data_type', 'study_id', 'pmid', 'doi_url', 'design', 'accession', 'note'])
    for t in RELATED:
        w.writerow([t[0], t[1], t[2], 'https://doi.org/' + t[3], t[4], t[5], t[6]])
with open(os.path.join(OUT, 'maize_N_to_check.tsv'), 'w', encoding='utf-8', newline='') as fh:
    w = csv.writer(fh, delimiter='\t', lineterminator='\n')
    w.writerow(['kind', 'label', 'pmid', 'link', 'why'])
    for t in TOCHECK:
        w.writerow(t)

pub = []
for o in out:
    pub.append({k: o.get(k, '') for k in ['study_id', 'details_level', 'year', 'first_author', 'journal', 'title',
        'doi_url', 'pmid', 'design_tokens', 'genotypes', 'n_genotypes', 'tissue', 'tissue_group', 'stage_age', 'growth_system',
        'system_group', 'N_treatments', 'other_factors', 'timepoints_duration', 'bio_reps', 'n_libraries', 'platform_layout',
        'accessions', 'accession_leads', 'accession_status', 'accession_note', 'repository', 'maizegdb', 'related_to',
        'other_omics', 'key_finding', 'notes']})
json.dump({'studies': pub,
           'related': [dict(zip(['type', 'study', 'pmid', 'doi', 'design', 'accession', 'note'], t)) for t in RELATED],
           'tocheck': [dict(zip(['kind', 'label', 'pmid', 'link', 'why'], t)) for t in TOCHECK]},
          open(os.path.join(OUT, 'catalogue.json'), 'w', encoding='utf-8'), ensure_ascii=False)

from collections import Counter
print('studies:', len(out))
print('fit:', Counter(o['fit'] for o in out))
print('acc status:', Counter(o['accession_status'] for o in out))
print('system:', Counter(o['system_group'] for o in out))
print('tissue:', Counter(t for o in out for t in o['tissue_group']))
print('design:', Counter(t for o in out for t in o['design_tokens']))
print('details:', Counter(o['details_level'] for o in out))
print('related:', len(RELATED), 'tocheck:', len(TOCHECK))
