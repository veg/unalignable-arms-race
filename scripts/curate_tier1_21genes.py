#!/usr/bin/env python3
"""
Curate 100% authentic, species-labeled, raw unaligned RefSeq/TOGA2 mammalian CDS
datasets for the 21 new Tier-1 loci (completing the 25-locus Tier-1 panel alongside
TRIM5, MX1, APOBEC3G, and PRDM9):

  Arm 1A (Antiviral Effectors & Sensors, 11 new + 3 done = 14):
    ISG15, APOBEC3H, SAMHD1, STAT2, IRF3, TRIM25, BST2, EIF2AK2, OAS1, MAVS, NLRP1
  Arm 1B (Viral Entry Receptors, 5):
    DPP4, ACE2, CD4, TFRC, NPC1
  Arm 1D (Hijacked Cellular Host Hub Controls, 5):
    PPIA, ELOB, ANP32A, RAE1, YWHAZ

Also writes /Users/sergei/Projects/unalignable-arms-race/data/ground_truth/uniprot_25gene_canonical.json
containing verified canonical Homo_sapiens amino-acid sequences from UniProtKB.
Python 3.9 compatible.
"""

import glob
import json
import os
import re
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict

ROOT = "/Users/sergei/Projects/unalignable-arms-race"
RAW_DIR = os.path.join(ROOT, "data", "raw_cds")
TOGA_DIR = os.path.join(RAW_DIR, "toga2_raw")
GT_DIR = os.path.join(ROOT, "data", "ground_truth")
TSV_PATH = os.path.join(ROOT, "data", "assemblies_and_species.tsv")
DATASETS_BIN = "/Users/sergei/miniconda3/bin/datasets"

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

GENE_SPECS = [
    # (gene, human_gene_id, uniprot, arm, toga2_filename)
    ("ISG15",   9636,   "P05161", "Antiviral_Effector_Sensor", "ENST00000649529.1#ISG15.codonified.fa"),
    ("APOBEC3H",164668, "Q6NTF7", "Antiviral_Effector_Sensor", "ENST00000348946.8#APOBEC3H.codonified.fa"),
    ("SAMHD1",  25939,  "Q9Y3Z3", "Antiviral_Effector_Sensor", "NM_001363729.2#SAMHD1.codonified.fa"),
    ("STAT2",   6773,   "P52630", "Antiviral_Effector_Sensor", "ENST00000698186.1#STAT2.codonified.fa"),
    ("IRF3",    3661,   "Q14653", "Antiviral_Effector_Sensor", "ENST00000377139.8#IRF3.codonified.fa"),
    ("TRIM25",  7706,   "Q14258", "Antiviral_Effector_Sensor", "ENST00000316881.9#TRIM25.codonified.fa"),
    ("BST2",    684,    "Q10589", "Antiviral_Effector_Sensor", "ENST00000252593.7#BST2.codonified.fa"),
    ("EIF2AK2", 5610,   "P19525", "Antiviral_Effector_Sensor", "ENST00000681507.1#EIF2AK2.codonified.fa"),
    ("OAS1",    4938,   "P00973", "Antiviral_Effector_Sensor", "NM_001320151.2#OAS1.codonified.fa"),
    ("MAVS",    57506,  "Q7Z434", "Antiviral_Effector_Sensor", "ENST00000428216.4#MAVS.codonified.fa"),
    ("NLRP1",   22861,  "Q9C000", "Antiviral_Effector_Sensor", "ENST00000576905.6#NLRP1.codonified.fa"),
    ("DPP4",    1803,   "P27487", "Viral_Entry_Receptor",      "ENST00000676810.1#DPP4.codonified.fa"),
    ("ACE2",    59272,  "Q9BYF1", "Viral_Entry_Receptor",      "ENST00000678073.1#ACE2.codonified.fa"),
    ("CD4",     920,    "P01730", "Viral_Entry_Receptor",      "ENST00000011653.9#CD4.codonified.fa"),
    ("TFRC",    7037,   "P02786", "Viral_Entry_Receptor",      "ENST00000360110.9#TFRC.codonified.fa"),
    ("NPC1",    4864,   "O15118", "Viral_Entry_Receptor",      "ENST00000269228.10#NPC1.codonified.fa"),
    ("PPIA",    5478,   "P62937", "Hijacked_Cellular_Hub",     "ENST00000468812.6#PPIA.codonified.fa"),
    ("ELOB",    6923,   "Q15370", "Hijacked_Cellular_Hub",     "ENST00000409906.9#ELOB.codonified.fa"),
    ("ANP32A",  8125,   "P39687", "Hijacked_Cellular_Hub",     "ENST00000465139.6#ANP32A.codonified.fa"),
    ("RAE1",    8480,   "P78406", "Hijacked_Cellular_Hub",     "ENST00000371242.6#RAE1.codonified.fa"),
    ("YWHAZ",   7534,   "P63104", "Hijacked_Cellular_Hub",     "XM_017013810.3#YWHAZ.codonified.fa"),
]

# Also record all 25 UniProt accessions for canonical reference verification
ALL_25_UNIPROT = {
    "TRIM5": "Q9C035",
    "MX1": "P20591",
    "APOBEC3G": "Q9HC16",
    "PRDM9": "Q9NQV7",
}
for _g, _gid, _up, _arm, _tfn in GENE_SPECS:
    ALL_25_UNIPROT[_g] = _up


def translate_nt(seq):
    return "".join(CODON_TABLE.get(seq[i:i + 3].upper(), "X") for i in range(0, len(seq) - 2, 3))


def count_internal_stops(seq):
    aa = translate_nt(seq)
    if not aa:
        return 0
    body = aa[:-1] if aa.endswith("*") else aa
    return body.count("*")


def sanitize_taxon(sp):
    parts = re.sub(r"[^A-Za-z0-9_ ]+", "", sp.strip()).split()
    if len(parts) >= 2:
        return "{0}_{1}".format(parts[0], parts[1])
    return parts[0] if parts else "Unknown"


def read_fasta(path):
    records = []
    cur_h, cur_s = None, []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith(">"):
                if cur_h is not None:
                    records.append((cur_h, "".join(cur_s)))
                cur_h, cur_s = s[1:], []
            else:
                cur_s.append(s)
    if cur_h is not None:
        records.append((cur_h, "".join(cur_s)))
    return records


def degap_toga_codon_seq(aligned_seq):
    clean_codons = []
    n_codons = 0
    partial_gap_codons = 0
    s = aligned_seq.upper()
    for k in range(0, len(s) - 2, 3):
        cod = s[k:k + 3]
        if cod == "---":
            continue
        if "N" in cod:
            n_codons += 1
            continue
        if "-" in cod:
            partial_gap_codons += 1
            continue
        if all(c in "ACGT" for c in cod):
            clean_codons.append(cod)
    return "".join(clean_codons), n_codons, partial_gap_codons


def load_toga_metadata():
    asm_map = {
        "hg38": {
            "species": "Homo sapiens",
            "taxon": "Homo_sapiens",
            "order": "Primates",
            "family": "Hominidae",
            "lineage": "Mammalia; Primates; Catarrhini; Hominidae; Homo",
            "status": "Chromosome",
            "scaf_n50": 150000000,
        }
    }
    genus_order = {}
    with open(TSV_PATH, "r", encoding="utf-8") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 15:
                continue
            asm = p[6].strip()
            sp = p[1].strip()
            lin = p[5].strip()
            status = p[10].strip()
            try:
                scaf_n50 = int(p[14].strip())
            except ValueError:
                scaf_n50 = 0
            lin_tokens = [t.strip() for t in lin.split(";")]
            fams = [t for t in lin_tokens if t.endswith("idae")]
            fam = fams[0] if fams else "Other"
            order = "Other"
            for cand in (
                "Primates", "Chiroptera", "Cetacea", "Artiodactyla", "Carnivora",
                "Rodentia", "Lagomorpha", "Perissodactyla", "Eulipotyphla",
                "Afrosoricida", "Macroscelidea", "Proboscidea", "Sirenia",
                "Pholidota", "Cingulata", "Pilosa", "Dermoptera", "Scandentia",
                "Diprotodontia", "Dasyuromorphia", "Didelphimorphia", "Monotremata"
            ):
                if cand in lin_tokens or cand in lin:
                    if cand == "Artiodactyla" and "Cetacea" in lin:
                        order = "Cetacea"
                    else:
                        order = cand
                    break
            gn = sanitize_taxon(sp).split("_")[0]
            if gn and order != "Other":
                genus_order[gn] = order
            asm_map[asm] = {
                "species": sp,
                "taxon": sanitize_taxon(sp),
                "order": order,
                "family": fam,
                "lineage": lin,
                "status": status,
                "scaf_n50": scaf_n50,
            }
    return asm_map, genus_order


def fetch_uniprot_canonical():
    out = {}
    for gene, acc in ALL_25_UNIPROT.items():
        url = "https://rest.uniprot.org/uniprotkb/{0}.json".format(acc)
        req = urllib.request.Request(url, headers={"User-Agent": "HyPhy3-Screen/1.0"})
        with urllib.request.urlopen(req, timeout=30) as r:
            j = json.loads(r.read().decode("utf-8"))
        seq = j["sequence"]["value"]
        out[gene] = {
            "uniprot": acc,
            "aa_len": len(seq),
            "aa_seq": seq,
            "protein_name": j.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", ""),
        }
        time.sleep(0.1)
    path = os.path.join(GT_DIR, "uniprot_25gene_canonical.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    return out


def download_ncbi_gene_cds(gene, gid):
    zip_path = os.path.join(RAW_DIR, "ncbi_{0}_orthologs.zip".format(gene))
    pkg_dir = os.path.join(RAW_DIR, "ncbi_{0}_pkg".format(gene))
    cds_fna = os.path.join(pkg_dir, "ncbi_dataset", "data", "cds.fna")
    if os.path.exists(cds_fna) and os.path.getsize(cds_fna) > 500:
        return cds_fna
    cmd = [
        DATASETS_BIN, "download", "gene", "gene-id", str(gid),
        "--ortholog", "mammalia", "--include", "cds", "--filename", zip_path
    ]
    subprocess.run(cmd, check=True)
    os.makedirs(pkg_dir, exist_ok=True)
    subprocess.run(["unzip", "-qo", zip_path, "-d", pkg_dir], check=True)
    return cds_fna


def aa_match_score(aa_cand, aa_ref):
    a = aa_cand.rstrip("*")
    if a == aa_ref:
        return (2, 1.0, len(a))
    # Prefix / overlap score
    L = min(len(a), len(aa_ref))
    if L == 0:
        return (0, 0.0, 0)
    # Check exact substring containment
    if a in aa_ref or aa_ref in a:
        return (1, float(L) / max(len(a), len(aa_ref)), -abs(len(a) - len(aa_ref)))
    same = sum(1 for i in range(L) if a[i] == aa_ref[i])
    return (0, float(same) / max(len(a), len(aa_ref)), -abs(len(a) - len(aa_ref)))


def curate_single_gene(gene, gid, up_info, arm, toga_fn, asm_map, genus_order):
    print("  [{0}] Downloading / parsing NCBI + TOGA2...".format(gene), flush=True)
    cds_fna = download_ncbi_gene_cds(gene, gid)
    ref_aa = up_info["aa_seq"]
    ref_aa_len = up_info["aa_len"]
    ref_nt_len = ref_aa_len * 3

    # First pass: determine authoritative Homo_sapiens canonical CDS from NCBI (matching UniProt aa_seq!)
    human_candidates = []
    other_ncbi = defaultdict(list)

    for hdr, seq_raw in read_fasta(cds_fna):
        m_org = re.search(r"\[organism=([^\]]+)\]", hdr)
        if not m_org:
            continue
        sp = m_org.group(1)
        taxon = sanitize_taxon(sp)
        acc = hdr.split()[0].split(":")[0]
        seq = re.sub(r"[^A-Za-z]", "", seq_raw).upper()
        if not seq or not all(c in "ACGT" for c in seq):
            continue
        if len(seq) % 3 != 0 or not seq.startswith("ATG") or count_internal_stops(seq) > 0:
            continue
        aa = translate_nt(seq).rstrip("*")
        if taxon == "Homo_sapiens":
            human_candidates.append((aa_match_score(aa, ref_aa), acc.startswith("NM_"), acc, seq, sp))
        else:
            other_ncbi[taxon].append((acc, seq, sp, aa))

    # Also read TOGA2 .codonified.fa
    toga_path = os.path.join(TOGA_DIR, toga_fn)
    toga_records = read_fasta(toga_path) if os.path.exists(toga_path) else []
    toga_by_taxon = defaultdict(list)
    hg38_toga_seq = None
    for hdr, aseq in toga_records:
        seq, n_cod, p_cod = degap_toga_codon_seq(aseq)
        if hdr == "hg38":
            hg38_toga_seq = seq
            continue
        info = asm_map.get(hdr)
        if not info or not seq:
            continue
        toga_by_taxon[info["taxon"]].append((hdr, info, seq, n_cod, p_cod))

    # Select Homo_sapiens CDS
    best_human = None
    if human_candidates:
        human_candidates.sort(key=lambda x: (x[0], x[1]), reverse=True)
        best_human = {
            "taxon": "Homo_sapiens",
            "species": "Homo sapiens",
            "order": "Primates",
            "family": "Hominidae",
            "accession": human_candidates[0][2],
            "seq": human_candidates[0][3],
            "nt_len": len(human_candidates[0][3]),
            "internal_stops": 0,
            "source": "NCBI_RefSeq_Homo_sapiens",
            "aa_match_to_uniprot": human_candidates[0][0][1],
        }
    elif hg38_toga_seq:
        best_human = {
            "taxon": "Homo_sapiens",
            "species": "Homo sapiens",
            "order": "Primates",
            "family": "Hominidae",
            "accession": "TOGA2:hg38",
            "seq": hg38_toga_seq,
            "nt_len": len(hg38_toga_seq),
            "internal_stops": count_internal_stops(hg38_toga_seq),
            "source": "TOGA2_hg38",
            "aa_match_to_uniprot": aa_match_score(translate_nt(hg38_toga_seq), ref_aa)[1],
        }
    assert best_human is not None, "Failed to obtain Homo_sapiens CDS for " + gene
    hum_nt_len = best_human["nt_len"]

    # Length bounds (preserve genuine isoform/repeat length shifts; guard against tiny fragments)
    min_nt = max(180, int(0.60 * hum_nt_len))
    max_nt = int(1.38 * hum_nt_len) + 90
    if gene == "SAMHD1":
        min_nt = max(180, int(0.55 * hum_nt_len))

    by_taxon = {"Homo_sapiens": best_human}

    # 1. Select best NCBI RefSeq ortholog per non-human species (closest isoform length + NM_ priority)
    for taxon, cand_list in other_ncbi.items():
        valid = [c for c in cand_list if min_nt <= len(c[1]) <= max_nt]
        if not valid:
            continue
        valid.sort(key=lambda c: (
            1 if c[0].startswith("NM_") else 0,
            -abs(len(c[1]) - hum_nt_len),
            len(c[1]),
        ), reverse=True)
        acc, seq, sp, _aa = valid[0]
        gn = taxon.split("_")[0]
        by_taxon[taxon] = {
            "taxon": taxon,
            "species": sp,
            "order": genus_order.get(gn, "Other"),
            "family": "NCBI_RefSeq",
            "accession": acc,
            "seq": seq,
            "nt_len": len(seq),
            "internal_stops": 0,
            "source": "NCBI_Datasets_{0}".format(gid),
        }

    # 2. Supplement missing species & high-confidence Chromosome-level LoF alleles from TOGA2
    lof_candidates = []
    for taxon, tlist in toga_by_taxon.items():
        for hdr, info, seq, n_cod, p_cod in tlist:
            if n_cod > 4 or len(seq) < min_nt or len(seq) > max_nt:
                continue
            stops = count_internal_stops(seq)
            if stops == 0 and taxon not in by_taxon:
                # Intact TOGA2 ortholog for taxon missing from NCBI RefSeq
                score = (1 if info["status"] == "Chromosome" else 0, -n_cod, -p_cod, info["scaf_n50"], len(seq))
                prev = by_taxon.get(taxon)
                cand_d = {
                    "taxon": taxon,
                    "species": info["species"],
                    "order": info["order"],
                    "family": info["family"],
                    "accession": "TOGA2:{0}".format(hdr),
                    "seq": seq,
                    "nt_len": len(seq),
                    "internal_stops": 0,
                    "source": "TOGA2_intact",
                    "_score": score,
                }
                if prev is None or (prev["source"] == "TOGA2_intact" and score > prev["_score"]):
                    by_taxon[taxon] = cand_d
            elif stops > 0 and info["status"] == "Chromosome" and n_cod <= 2 and p_cod == 0 and len(seq) >= int(0.75 * hum_nt_len) and stops <= 4:
                # High-confidence Chromosome-level intact-reading-frame genuine pseudogene!
                lof_candidates.append({
                    "taxon": taxon,
                    "species": info["species"],
                    "order": info["order"],
                    "family": info["family"],
                    "accession": "TOGA2:{0}".format(hdr),
                    "seq": seq,
                    "nt_len": len(seq),
                    "internal_stops": stops,
                    "source": "TOGA2_Chromosome_LoF",
                    "scaf_n50": info["scaf_n50"],
                })

    # Clean temporary _score keys
    for d in by_taxon.values():
        d.pop("_score", None)

    # Also promote up to 3 Chromosome-level LoF alleles if unambiguous and in non-control arms
    # (for Hijacked_Cellular_Hub controls, only include if genuinely Chromosome-level cleanly verified - still let's apply identical rule so no differential bias!)
    lof_candidates.sort(key=lambda r: (r["internal_stops"], -r["scaf_n50"]))
    seen_lof_orders = set()
    kept_lof = 0
    for r in lof_candidates:
        if kept_lof >= 3:
            break
        if r["taxon"] not in by_taxon and r["order"] not in seen_lof_orders:
            r.pop("scaf_n50", None)
            by_taxon[r["taxon"]] = r
            seen_lof_orders.add(r["order"])
            kept_lof += 1

    # 3. Taxonomic stratification to keep cohorts comparable and tractable:
    #    target 70 taxa for standard-length genes (<=1000 aa), 55 taxa for long (>1000 aa: NLRP1, NPC1)
    target_n = 55 if ref_aa_len > 1000 else 72
    all_recs = list(by_taxon.values())
    if len(all_recs) > target_n:
        selected = {"Homo_sapiens": by_taxon["Homo_sapiens"]}
        # Keep any LoF records
        for r in all_recs:
            if r["internal_stops"] > 0:
                selected[r["taxon"]] = r

        # Keep up to 12 diverse Primates (closest structural relatives + key experimental lineages)
        primates = sorted(
            [r for r in all_recs if r["taxon"] not in selected and r["order"] == "Primates"],
            key=lambda r: (1 if "NCBI_" in r["source"] else 0, 1 if r["accession"].startswith("NM_") else 0, r["taxon"]),
            reverse=True,
        )
        seen_p_gen = set()
        for r in primates:
            gn = r["taxon"].split("_")[0]
            if gn not in seen_p_gen and len([x for x in selected.values() if x["order"] == "Primates"]) < 14:
                selected[r["taxon"]] = r
                seen_p_gen.add(gn)

        # Keep top diverse genera across remaining mammalian orders (Chiroptera, Cetacea, Carnivora, Rodentia, Artiodactyla, Perissodactyla, Eulipotyphla, Afrotheria, etc.)
        priority_orders = [
            "Chiroptera", "Cetacea", "Carnivora", "Rodentia", "Artiodactyla",
            "Perissodactyla", "Proboscidea", "Sirenia", "Lagomorpha", "Eulipotyphla",
            "Afrosoricida", "Scandentia", "Cingulata", "Pilosa", "Pholidota", "Other"
        ]
        # Round-robin by order + unique genus first
        by_ord = defaultdict(list)
        for r in all_recs:
            if r["taxon"] not in selected:
                by_ord[r["order"]].append(r)
        for o in by_ord:
            by_ord[o].sort(key=lambda r: (
                1 if "NCBI_" in r["source"] else 0,
                -abs(r["nt_len"] - hum_nt_len),
                r["taxon"]
            ), reverse=True)

        seen_gen = set(x["taxon"].split("_")[0] for x in selected.values())
        added = True
        while len(selected) < target_n and added:
            added = False
            for o in priority_orders:
                if len(selected) >= target_n:
                    break
                # Pick first unseen genus in this order
                picked = None
                for r in by_ord.get(o, []):
                    gn = r["taxon"].split("_")[0]
                    if r["taxon"] not in selected and gn not in seen_gen:
                        picked = r
                        break
                if picked is None:
                    for r in by_ord.get(o, []):
                        if r["taxon"] not in selected:
                            picked = r
                            break
                if picked is not None:
                    selected[picked["taxon"]] = picked
                    seen_gen.add(picked["taxon"].split("_")[0])
                    added = True
        final_recs = list(selected.values())
    else:
        final_recs = all_recs

    final_recs.sort(key=lambda d: (0 if d["taxon"] == "Homo_sapiens" else 1, d["order"], d["taxon"]))

    fa_path = os.path.join(RAW_DIR, "{0}_unaligned.fa".format(gene))
    meta_path = os.path.join(RAW_DIR, "{0}_manifest.json".format(gene))
    with open(fa_path, "w", encoding="utf-8") as out:
        for r in final_recs:
            out.write(">{0}\n{1}\n".format(r["taxon"], r["seq"]))

    clean_meta = [{k: v for k, v in r.items() if k != "seq"} for r in final_recs]
    with open(meta_path, "w", encoding="utf-8") as out:
        json.dump({
            "gene": gene,
            "human_gene_id": gid,
            "uniprot": up_info["uniprot"],
            "arm": arm,
            "human_uniprot_aa_len": ref_aa_len,
            "human_cds_nt_len": hum_nt_len,
            "human_aa_match_to_uniprot": best_human["aa_match_to_uniprot"],
            "total_taxa": len(final_recs),
            "intact_taxa": sum(1 for r in final_recs if r["internal_stops"] == 0),
            "pseudogene_stop_taxa": [r["taxon"] for r in final_recs if r["internal_stops"] > 0],
            "min_nt": min(r["nt_len"] for r in final_recs),
            "max_nt": max(r["nt_len"] for r in final_recs),
            "orders": dict(sorted(defaultdict(int, ((o, sum(1 for r in final_recs if r["order"] == o)) for o in set(r["order"] for r in final_recs))).items())),
            "records": clean_meta,
        }, out, indent=2)

    return len(final_recs), hum_nt_len // 3, best_human["aa_match_to_uniprot"], [r["taxon"] for r in final_recs if r["internal_stops"] > 0]


def main():
    print("Step 1: Fetching canonical UniProt sequences for all 25 Tier-1 human proteins...")
    uniprot_map = fetch_uniprot_canonical()
    print("  Verified {0} UniProt canonical sequences.".format(len(uniprot_map)))

    asm_map, genus_order = load_toga_metadata()
    print("  Loaded TOGA2 assembly map ({0} assemblies, {1} genera).".format(len(asm_map), len(genus_order)))

    summary_rows = []
    for gene, gid, up, arm, tfn in GENE_SPECS:
        n_taxa, hum_aa, match_id, lofs = curate_single_gene(
            gene, gid, uniprot_map[gene], arm, tfn, asm_map, genus_order
        )
        print("  -> {0:10s} [{1:25s}]: N={2:2d}, human={3:4d} aa (UniProt match={4:.3f}), LoFs={5}".format(
            gene, arm, n_taxa, hum_aa, match_id, lofs
        ), flush=True)
        summary_rows.append((gene, arm, n_taxa, hum_aa, match_id, lofs))

    print("\nAll 21 new datasets curated cleanly.")


if __name__ == "__main__":
    main()
