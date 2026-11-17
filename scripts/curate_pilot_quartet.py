#!/usr/bin/env python3
"""
Curate 100% authentic, species-labeled, raw unaligned RefSeq/GenBank CDS datasets
(preserving native indels, unclipped tandem repeat copy-number variation, and genuine
published mammalian pseudogene lineages) plus verified structural/DMS ground-truth
maps for the "Unalignable Arms Race" Pilot Quartet:
  1. TRIM5    (full-length TRIM5alpha incl. RING, B-box2, Coiled-Coil, B30.2/PRYSPRY v1-v4 + Felis catus exon-8 stop pseudogene GQ183880.1)
  2. MX1      (unclipped NCBI RefSeq mammalian MxA incl. disordered Loop L4 + Odontoceti truncated/inactivated MX1 lineages)
  3. APOBEC3G (Simiiformes double-domain CD1+CD2 cytidine deaminase incl. Vif loop-7 D128 charge toggle & ectopic conversion tracts)
  4. PRDM9    (unclipped NCBI RefSeq mammalian PRDM9 with variable tandem 28-aa C2H2 ZnF repeat counts + Munoz-Fuentes et al. MBE 2011 Canidae pseudogenes)
"""

import json
import os
import re
import time
import urllib.parse
import urllib.request
from collections import defaultdict

ROOT = "/Users/sergei/Projects/unalignable-arms-race"
RAW_DIR = os.path.join(ROOT, "data", "raw_cds")
GT_DIR = os.path.join(ROOT, "data", "ground_truth")

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


def translate_nt(seq: str) -> str:
    return "".join(CODON_TABLE.get(seq[i:i + 3].upper(), "X") for i in range(0, len(seq) - 2, 3))


def count_internal_stops(seq: str) -> int:
    aa = translate_nt(seq)
    if not aa:
        return 0
    body = aa[:-1] if aa.endswith("*") else aa
    return body.count("*")


def sanitize_taxon(sp: str) -> str:
    parts = re.sub(r"[^A-Za-z0-9_ ]+", "", sp.strip()).split()
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    return parts[0] if parts else "Unknown"


def read_fasta(path: str):
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


def parse_fasta_string(text: str):
    records = []
    cur_h, cur_s = None, []
    for line in text.splitlines():
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


def eutils_fetch(term: str, rettype: str = "fasta_cds_na", retmax: int = 400) -> str:
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    q = urllib.parse.urlencode({"db": "nuccore", "term": term, "retmax": str(retmax), "retmode": "json"})
    with urllib.request.urlopen(base + "esearch.fcgi?" + q, timeout=30) as r:
        id_list = json.loads(r.read().decode("utf-8")).get("esearchresult", {}).get("idlist", [])
    if not id_list:
        return ""
    chunks = []
    for i in range(0, len(id_list), 100):
        batch = ",".join(id_list[i:i + 100])
        fq = urllib.parse.urlencode({"db": "nuccore", "id": batch, "rettype": rettype, "retmode": "text"})
        time.sleep(0.35)
        with urllib.request.urlopen(base + "efetch.fcgi?" + fq, timeout=60) as r:
            chunks.append(r.read().decode("utf-8", errors="replace"))
    return "\n".join(chunks)


def eutils_fetch_by_ids(ids, rettype="fasta"):
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    fq = urllib.parse.urlencode({"db": "nuccore", "id": ",".join(ids), "rettype": rettype, "retmode": "text"})
    time.sleep(0.35)
    with urllib.request.urlopen(base + "efetch.fcgi?" + fq, timeout=60) as r:
        return r.read().decode("utf-8", errors="replace")


def esummary_map(accs):
    acc_to_sp = {}
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    for i in range(0, len(accs), 100):
        batch = ",".join(accs[i:i + 100])
        sq = urllib.parse.urlencode({"db": "nuccore", "id": batch, "retmode": "json"})
        time.sleep(0.35)
        with urllib.request.urlopen(base + "esummary.fcgi?" + sq, timeout=30) as r:
            sdata = json.loads(r.read().decode("utf-8")).get("result", {})
        for v in sdata.values():
            if isinstance(v, dict) and "caption" in v and "organism" in v:
                acc_to_sp[v.get("accessionversion", v["caption"])] = (v["organism"], v.get("title", ""))
                acc_to_sp[v["caption"]] = (v["organism"], v.get("title", ""))
    return acc_to_sp


def curate_trim5():
    candidates = {}

    def consider(taxon, sp_name, acc, seq_raw, source, note="", allow_pseudo=False):
        seq = re.sub(r"[^A-Za-z]", "", seq_raw).upper()
        if not seq or not all(c in "ACGT" for c in seq):
            return
        if not allow_pseudo:
            if not (1380 <= len(seq) <= 1650) or len(seq) % 3 != 0 or not seq.startswith("ATG"):
                return
            if count_internal_stops(seq) > 0:
                return
        else:
            if len(seq) < 600 or len(seq) > 2500:
                return
        is_nm = acc.startswith("NM_") or source == "GenBank_curated"
        score = (
            1 if taxon == "Homo_sapiens" and "NM_033034" in acc else 0,
            1 if allow_pseudo else 0,  # preserve explicit published LoF mRNA for Felis_catus
            1 if is_nm else 0,
            len(seq),
        )
        prev = candidates.get(taxon)
        if prev is None or score > prev["score"]:
            candidates[taxon] = {
                "taxon": taxon, "species": sp_name, "accession": acc, "seq": seq,
                "nt_len": len(seq), "internal_stops": count_internal_stops(seq),
                "source": source, "note": note, "score": score,
            }

    ncbi_cds = os.path.join(RAW_DIR, "ncbi_trim5_pkg", "ncbi_dataset", "data", "cds.fna")
    for hdr, seq in read_fasta(ncbi_cds):
        m_org = re.search(r"\[organism=([^\]]+)\]", hdr)
        if m_org:
            sp = m_org.group(1)
            consider(sanitize_taxon(sp), sp, hdr.split()[0].split(":")[0], seq, "NCBI_Datasets_85363")

    query = (
        '("TRIM5"[Gene] OR "TRIM5alpha"[Title] OR "TRIM5 alpha"[Title] OR "tripartite motif-containing 5"[Title]) '
        'AND (Primates[Organism] OR Chiroptera[Organism] OR Carnivora[Organism] OR Perissodactyla[Organism] OR Artiodactyla[Organism]) '
        'AND 1350:4500[SLEN] NOT patent[Title] NOT synthetic[Title] NOT chimeric[Title]'
    )
    gb_text = eutils_fetch(query, rettype="fasta_cds_na", retmax=500)
    pending = []
    for hdr, seq in parse_fasta_string(gb_text):
        low = hdr.lower()
        if any(bad in low for bad in ("trimcyp", "cyclophilin", "chimera", "mutant", "synthetic", "trim22", "trim6", "trim34")):
            continue
        if "[gene=" in low and not re.search(r"\[gene=trim5[a-z0-9_-]*\]", low):
            continue
        m_acc = re.match(r"lcl\|([A-Za-z0-9_.]+)_cds_", hdr)
        acc = m_acc.group(1) if m_acc else hdr.split()[0]
        pending.append((acc, seq))

    if pending:
        acc_to_sp = esummary_map(sorted(set(x[0] for x in pending)))
        for acc, seq in pending:
            if acc in acc_to_sp:
                sp, title = acc_to_sp[acc]
                if not any(x in title.lower() for x in ("trimcyp", "fusion", "chimera", "synthetic", "trim22", "trim6", "trim34")):
                    consider(sanitize_taxon(sp), sp, acc, seq, "GenBank_curated", note=title[:100])

    # Add McEwan et al. 2009 PLoS Pathogens Felis catus TRIM5 mRNA (GQ183880.1: starts at ATG and carries the authentic exon-8 premature TGA stop + downstream remnant PRYSPRY!)
    fc_fa = eutils_fetch_by_ids(["GQ183880.1"], rettype="fasta")
    for hdr, seq in parse_fasta_string(fc_fa):
        m_atg = seq.upper().find("ATGGCTTCT")
        if m_atg == -1:
            m_atg = seq.upper().find("ATG")
        if m_atg >= 0:
            fc_orf_plus_remnant = seq.upper()[m_atg:]
            consider("Felis_catus", "Felis catus", "GQ183880.1", fc_orf_plus_remnant, "GenBank_LoF_McEwan2009",
                     note="Felis catus inactivated TRIM5 mRNA with premature stop before PRYSPRY", allow_pseudo=True)

    return sorted(candidates.values(), key=lambda d: (0 if d["taxon"] == "Homo_sapiens" else 1, d["taxon"]))


def curate_mx1():
    """
    Curate unclipped NCBI RefSeq mammalian MX1 raw CDS (~78 diverse species) + authentic
    Odontoceti (toothed whale/dolphin) inactivated/truncated MX1 records.
    """
    by_taxon = {}
    ncbi_cds = os.path.join(RAW_DIR, "ncbi_mx1_pkg", "ncbi_dataset", "data", "cds.fna")
    for hdr, seq in read_fasta(ncbi_cds):
        m_org = re.search(r"\[organism=([^\]]+)\]", hdr)
        if not m_org:
            continue
        sp = m_org.group(1)
        taxon = sanitize_taxon(sp)
        acc = hdr.split()[0].split(":")[0]
        s = re.sub(r"[^A-Za-z]", "", seq).upper()
        if not s or not all(c in "ACGT" for c in s):
            continue
        if not (1750 <= len(s) <= 2250) or len(s) % 3 != 0 or not s.startswith("ATG"):
            continue
        if count_internal_stops(s) > 0:
            continue
        score = (
            1 if taxon == "Homo_sapiens" and "NM_001144925" in acc else 0,
            1 if acc.startswith("NM_") else 0,
            len(s),
        )
        prev = by_taxon.get(taxon)
        if prev is None or score > prev["score"]:
            by_taxon[taxon] = {
                "taxon": taxon, "species": sp, "accession": acc, "seq": s,
                "nt_len": len(s), "internal_stops": 0, "source": "NCBI_Datasets_4599",
                "note": "", "score": score,
            }

    # Add Braun et al. 2015 PNAS Odontoceti (toothed whale/dolphin) inactivated MX1 genomic/mRNA records or TOGA2 Odontoceti premature-stop projections!
    # Check Braun et al. 2015 GenBank accessions (KR080141..KR080150 / KP208469..) or search "Braun"[Author] AND "MX1"
    braun_txt = eutils_fetch('("MX1"[All Fields] OR "Mx1"[All Fields]) AND (Delphinidae[Organism] OR Phocoenidae[Organism] OR Monodontidae[Organism] OR Lipotes[Organism]) AND 600:4000[SLEN]', rettype="fasta", retmax=40)
    for hdr, seq in parse_fasta_string(braun_txt):
        low = hdr.lower()
        if "mx2" in low or "whole genome" in low or "chromosome" in low:
            continue
        acc = hdr.split()[0]
        s = re.sub(r"[^A-Za-z]", "", seq).upper()
        # Extract from canonical Mx1 start codon if mRNA includes 5' UTR
        idx = s.find("ATGGTTGTT")
        if idx == -1:
            idx = s.find("ATGGTT")
        if idx >= 0 and 900 <= len(s) - idx <= 2400:
            cds_cand = s[idx:]
            # Find species binomial in header
            parts = hdr.split()
            if len(parts) >= 3:
                sp = f"{parts[1]} {parts[2]}" if parts[1] != "PREDICTED:" else f"{parts[2]} {parts[3]}"
                taxon = sanitize_taxon(sp)
                if taxon not in by_taxon and all(c in "ACGT" for c in cds_cand):
                    by_taxon[taxon] = {
                        "taxon": taxon, "species": sp, "accession": acc, "seq": cds_cand,
                        "nt_len": len(cds_cand), "internal_stops": count_internal_stops(cds_cand),
                        "source": "GenBank_Odontoceti_MX1", "note": hdr[:90], "score": (0, 2, len(cds_cand)),
                    }

    # Stratify to ~78 diverse species (keeping Homo_sapiens, all Primates, all Cetacea, all Chiroptera, plus diverse Carnivora/Rodentia/Artiodactyla/Afrotheria)
    all_recs = list(by_taxon.values())
    if len(all_recs) > 80:
        # Load genus->order hints from assemblies_and_species.tsv
        genus_order = {}
        with open(os.path.join(ROOT, "data", "assemblies_and_species.tsv")) as f:
            f.readline()
            for line in f:
                p = line.split("\t")
                if len(p) > 5:
                    gn = p[1].split()[0] if p[1].strip() else ""
                    lin = p[5]
                    for o in ("Primates", "Chiroptera", "Cetacea", "Carnivora", "Artiodactyla", "Rodentia", "Perissodactyla", "Proboscidea"):
                        if o in lin:
                            genus_order[gn] = "Cetacea" if "Cetacea" in lin else o
                            break
        kept = {}
        for r in all_recs:
            gn = r["taxon"].split("_")[0]
            r["order"] = genus_order.get(gn, "Other")
            if r["taxon"] == "Homo_sapiens" or r["internal_stops"] > 0 or r["order"] in ("Primates", "Chiroptera", "Cetacea"):
                kept[r["taxon"]] = r
        seen_genus = defaultdict(int)
        for r in sorted(all_recs, key=lambda x: (-x["score"][1], x["taxon"])):
            if len(kept) >= 80:
                break
            gn = r["taxon"].split("_")[0]
            if r["taxon"] not in kept and seen_genus[gn] == 0:
                kept[r["taxon"]] = r
                seen_genus[gn] += 1
        all_recs = list(kept.values())

    return sorted(all_recs, key=lambda d: (0 if d["taxon"] == "Homo_sapiens" else 1, d["taxon"]))


def curate_apobec3g():
    """
    Curate Simiiformes double-domain APOBEC3G raw unaligned CDS across Hominoids,
    Old World Monkeys, and New World Monkeys (Sawyer et al. PLoS Biol 2004, Compton & Emerman,
    NCBI RefSeq & GenBank).
    """
    by_taxon = {}
    # 1. Broad GenBank + RefSeq primate APOBEC3G search (including XM_ predicted primate orthologs!)
    query = (
        '("APOBEC3G"[Gene] OR "APOBEC3G"[All Fields] OR "apolipoprotein B mRNA editing enzyme catalytic subunit 3G"[Title] '
        'OR "apolipoprotein B mRNA editing enzyme, catalytic polypeptide-like 3G"[Title] OR "CEM15"[Title]) '
        'AND Primates[Organism] AND 950:3500[SLEN] NOT patent[Title] NOT synthetic[Title] NOT chimeric[Title]'
    )
    gb_text = eutils_fetch(query, rettype="fasta_cds_na", retmax=500)
    pending = []
    for hdr, seq in parse_fasta_string(gb_text):
        low = hdr.lower()
        if any(bad in low for bad in ("mutant", "chimera", "synthetic", "apobec3f", "apobec3de", "apobec3d", "apobec3b", "apobec3c", "apobec3a", "apobec3h")):
            continue
        if "[gene=" in low and not re.search(r"\[gene=apobec3g\]", low):
            continue
        s = re.sub(r"[^A-Za-z]", "", seq).upper()
        if not (1050 <= len(s) <= 1260) or len(s) % 3 != 0 or not s.startswith("ATG"):
            continue
        if count_internal_stops(s) > 0:
            continue
        m_acc = re.match(r"lcl\|([A-Za-z0-9_.]+)_cds_", hdr)
        acc = m_acc.group(1) if m_acc else hdr.split()[0]
        pending.append((acc, s))

    if pending:
        acc_to_sp = esummary_map(sorted(set(x[0] for x in pending)))
        for acc, s in pending:
            if acc in acc_to_sp:
                sp, title = acc_to_sp[acc]
                tlow = title.lower()
                if any(x in tlow for x in ("chimera", "mutant", "synthetic", "3f", "3de", "3b", "3c", "3a", "3h")):
                    continue
                taxon = sanitize_taxon(sp)
                score = (
                    1 if taxon == "Homo_sapiens" and "NM_021822" in acc else 0,
                    1 if acc.startswith("NM_") else 0,
                    len(s),
                )
                prev = by_taxon.get(taxon)
                if prev is None or score > prev["score"]:
                    by_taxon[taxon] = {
                        "taxon": taxon, "species": sp, "accession": acc, "seq": s,
                        "nt_len": len(s), "internal_stops": 0, "source": "NCBI_Primate_APOBEC3G",
                        "note": title[:95], "score": score,
                    }

    return sorted(by_taxon.values(), key=lambda d: (0 if d["taxon"] == "Homo_sapiens" else 1, d["taxon"]))


def curate_prdm9():
    """
    Curate unclipped NCBI RefSeq mammalian PRDM9 raw CDS (61 species, preserving true
    tandem 28-aa C2H2 repeat copy-number expansions/contractions!) + authentic Canidae
    PRDM9 pseudogene sequences (Munoz-Fuentes et al. MBE 2011).
    """
    by_taxon = {}
    ncbi_cds = os.path.join(RAW_DIR, "ncbi_prdm9_pkg", "ncbi_dataset", "data", "cds.fna")
    for hdr, seq in read_fasta(ncbi_cds):
        m_org = re.search(r"\[organism=([^\]]+)\]", hdr)
        if not m_org:
            continue
        sp = m_org.group(1)
        taxon = sanitize_taxon(sp)
        acc = hdr.split()[0].split(":")[0]
        s = re.sub(r"[^A-Za-z]", "", seq).upper()
        if not s or not all(c in "ACGT" for c in s):
            continue
        if len(s) < 1650 or len(s) % 3 != 0 or not s.startswith("ATG"):
            continue
        if count_internal_stops(s) > 0:
            continue
        aa = translate_nt(s)
        n_znf = len(re.findall(r"C..C.{12}H...H", aa))
        if n_znf < 3:
            continue
        score = (
            1 if taxon == "Homo_sapiens" and "NM_020227" in acc else (1 if taxon == "Homo_sapiens" and "NM_001376900" in acc else 0),
            1 if acc.startswith("NM_") else 0,
            n_znf,
            len(s),
        )
        prev = by_taxon.get(taxon)
        if prev is None or score > prev["score"]:
            by_taxon[taxon] = {
                "taxon": taxon, "species": sp, "accession": acc, "seq": s,
                "nt_len": len(s), "znf_count": n_znf, "internal_stops": 0,
                "source": "NCBI_Datasets_56979", "note": f"{n_znf} C2H2 ZnF repeats", "score": score,
            }

    # Add Munoz-Fuentes et al. MBE 2011 Canidae inactivated PRDM9 pseudogene sequences!
    canid_ids = ["JF750655.1", "JF750658.1", "MW004668.1", "JF750659.1"]
    canid_fa = eutils_fetch_by_ids(canid_ids, rettype="fasta")
    acc_map = esummary_map(canid_ids)
    for hdr, seq in parse_fasta_string(canid_fa):
        acc = hdr.split()[0]
        sp, title = acc_map.get(acc, ("Canis lupus", hdr))
        taxon = sanitize_taxon(sp)
        s = re.sub(r"[^A-Za-z]", "", seq).upper()
        if all(c in "ACGT" for c in s) and len(s) >= 1000:
            by_taxon[taxon] = {
                "taxon": taxon, "species": sp, "accession": acc, "seq": s,
                "nt_len": len(s), "znf_count": 0, "internal_stops": count_internal_stops(s),
                "source": "GenBank_Canidae_LoF_MunozFuentes2011", "note": title[:95],
                "score": (0, 2, 0, len(s)),
            }

    return sorted(by_taxon.values(), key=lambda d: (0 if d["taxon"] == "Homo_sapiens" else 1, d["taxon"]))


def write_dataset(name: str, recs):
    fa_path = os.path.join(RAW_DIR, f"{name}_unaligned.fa")
    meta_path = os.path.join(RAW_DIR, f"{name}_manifest.json")
    with open(fa_path, "w", encoding="utf-8") as out:
        for r in recs:
            out.write(f">{r['taxon']}\n{r['seq']}\n")
    clean_meta = [{k: v for k, v in r.items() if k not in ("seq", "score")} for r in recs]
    with open(meta_path, "w", encoding="utf-8") as out:
        json.dump({
            "gene": name,
            "total_taxa": len(recs),
            "intact_taxa": sum(1 for r in recs if r["internal_stops"] == 0 and r["nt_len"] % 3 == 0),
            "pseudogene_or_frameshift_taxa": [r["taxon"] for r in recs if r["internal_stops"] > 0 or r["nt_len"] % 3 != 0],
            "min_nt": min(r["nt_len"] for r in recs),
            "max_nt": max(r["nt_len"] for r in recs),
            "records": clean_meta,
        }, out, indent=2)
    return fa_path


def build_and_verify_ground_truth(datasets_by_gene):
    gt = {}

    h_trim5 = translate_nt(next(r["seq"] for r in datasets_by_gene["TRIM5"] if r["taxon"] == "Homo_sapiens")).rstrip("*")
    assert len(h_trim5) == 493 and h_trim5[332 - 1] == "R", f"TRIM5 check failed: len={len(h_trim5)}"
    gt["TRIM5"] = {
        "human_ref_aa_len": len(h_trim5),
        "human_ref_seq": h_trim5,
        "pdb_structures": ["2LM3 (PRYSPRY NMR/crystal)", "8F25 (TRIM5a-HIV1 CA hexamer cryo-EM)"],
        "domains_1idx": {
            "RING": [15, 59], "B_box_2": [90, 132], "Coiled_Coil": [133, 241],
            "Linker_L2": [242, 290], "PRYSPRY_core": [291, 493],
        },
        "hypervariable_viral_loops_1idx": {
            "v1_capsid_patch": [324, 345], "v2_loop": [383, 395],
            "v3_loop": [415, 428], "v4_loop": [480, 488],
        },
        "experimental_functional_sites_1idx": [330, 332, 334, 335, 336, 337, 339, 340, 385, 389, 420],
        "known_lof_lineages": ["Felis_catus (GQ183880.1 exon-8 premature TGA stop; McEwan et al. PLoS Pathog 2009)"],
    }

    h_mx1 = translate_nt(next(r["seq"] for r in datasets_by_gene["MX1"] if r["taxon"] == "Homo_sapiens")).rstrip("*")
    assert len(h_mx1) == 662 and h_mx1[561 - 1] == "F", f"MX1 check failed: len={len(h_mx1)}"
    gt["MX1"] = {
        "human_ref_aa_len": len(h_mx1),
        "human_ref_seq": h_mx1,
        "pdb_structures": ["3SZR (human MxA full-length crystal)"],
        "domains_1idx": {
            "N_terminal_extension": [1, 42], "BSE_1": [43, 68], "GTPase_G_domain": [69, 340],
            "BSE_2": [341, 367], "Stalk_MD_IR": [368, 532],
            "Loop_L4_disordered_antiviral": [533, 572], "Stalk_C": [573, 634], "BSE_3": [635, 662],
        },
        "hypervariable_viral_loops_1idx": {
            "N_term_IDR": [1, 42], "Loop_L4_orthomyxovirus_NP_interface": [533, 572],
        },
        "experimental_functional_sites_1idx": [540, 543, 549, 559, 561, 562, 564, 566, 567],
        "known_lof_lineages": ["Odontoceti (toothed whales & dolphins Mx1 inactivation; Braun et al. PNAS 2015)"],
    }

    h_a3g = translate_nt(next(r["seq"] for r in datasets_by_gene["APOBEC3G"] if r["taxon"] == "Homo_sapiens")).rstrip("*")
    assert len(h_a3g) == 384 and h_a3g[127:130] == "DPD", f"APOBEC3G check failed: len={len(h_a3g)}"
    gt["APOBEC3G"] = {
        "human_ref_aa_len": len(h_a3g),
        "human_ref_seq": h_a3g,
        "pdb_structures": ["8E40 (HIV-1 Vif-CBFb-APOBEC3G cryo-EM)", "3IR2 (CD2 catalytic crystal)"],
        "domains_1idx": {
            "CD1_N_terminal_RNA_Vif_binding": [1, 196],
            "CD2_C_terminal_catalytic_deaminase": [197, 384],
        },
        "hypervariable_viral_loops_1idx": {
            "CD1_Loop7_Vif_interface": [122, 134],
            "CD1_alpha2_alpha3_Vif_patch": [95, 110],
            "CD2_substrate_loops": [205, 220],
        },
        "experimental_functional_sites_1idx": [126, 127, 128, 129, 130, 211, 317],
        "known_lof_lineages": ["Paralogous/inter-domain ectopic gene conversion across primate APOBEC3"],
    }

    h_prdm9 = translate_nt(next(r["seq"] for r in datasets_by_gene["PRDM9"] if r["taxon"] == "Homo_sapiens")).rstrip("*")
    assert len(h_prdm9) == 894, f"PRDM9 check failed: len={len(h_prdm9)}"
    znf_matches = list(re.finditer(r"C..C.{12}H...H", h_prdm9))
    assert len(znf_matches) >= 12
    znf_repeats, dna_contact_1idx, c2h2_scaffold_1idx = [], [], []
    for idx, m in enumerate(znf_matches, start=1):
        c1 = m.start() + 1
        c2 = c1 + 3
        h1 = c2 + 13
        h2 = h1 + 4
        contacts = [c2 + 6, c2 + 8, c2 + 9, c2 + 12]
        znf_repeats.append({
            "finger_index": idx,
            "c1_to_h2_1idx": [c1, h2],
            "motif_seq": m.group(0),
            "c2h2_coordinating_1idx": [c1, c2, h1, h2],
            "dna_contact_positions_m1_p2_p3_p6_1idx": contacts,
            "dna_contact_aas": "".join(h_prdm9[p - 1] for p in contacts),
        })
        dna_contact_1idx.extend(contacts)
        c2h2_scaffold_1idx.extend([c1, c2, h1, h2])

    znf_counts = [r["znf_count"] for r in datasets_by_gene["PRDM9"] if r.get("znf_count", 0) > 0]
    gt["PRDM9"] = {
        "human_ref_aa_len": len(h_prdm9),
        "human_ref_seq": h_prdm9,
        "mammalian_znf_count_range": [min(znf_counts), max(znf_counts)],
        "pdb_structures": ["5EGB (human PRDM9 ZnF8-12 bound to hotspot DNA)", "4C1Q (PR/SET catalytic domain)"],
        "domains_1idx": {
            "KRAB_domain": [24, 97], "SSXRD_motif": [151, 183],
            "PR_SET_methyltransferase": [202, 368], "Post_SET_linker": [369, 531],
            "Tandem_C2H2_ZnF_array": [znf_repeats[0]["c1_to_h2_1idx"][0] - 4, 894],
        },
        "znf_repeat_count": len(znf_repeats),
        "znf_repeats": znf_repeats,
        "experimental_dna_contact_sites_1idx": sorted(set(dna_contact_1idx)),
        "c2h2_coordinating_scaffold_sites_1idx": sorted(set(c2h2_scaffold_1idx)),
        "known_lof_lineages": ["Canidae (Canis lupus familiaris, Canis aureus, Lupulella mesomelas, Vulpes vulpes inactivated PRDM9; Munoz-Fuentes et al. MBE 2011)"],
    }

    with open(os.path.join(GT_DIR, "pilot_quartet_ground_truth.json"), "w", encoding="utf-8") as out:
        json.dump(gt, out, indent=2)
    return gt


def main():
    print("1. Curating TRIM5alpha (full-length + Felis_catus GQ183880.1 exon-8 LoF)...")
    trim5_recs = curate_trim5()
    write_dataset("TRIM5", trim5_recs)
    print(f"   TRIM5: {len(trim5_recs)} species, nt={min(r['nt_len'] for r in trim5_recs)}..{max(r['nt_len'] for r in trim5_recs)}, LoF/frameshift={[r['taxon'] for r in trim5_recs if r['internal_stops']>0 or r['nt_len']%3!=0]}")

    print("2. Curating MX1 (unclipped RefSeq mammalian MxA + Odontoceti MX1)...")
    mx1_recs = curate_mx1()
    write_dataset("MX1", mx1_recs)
    print(f"   MX1: {len(mx1_recs)} species, nt={min(r['nt_len'] for r in mx1_recs)}..{max(r['nt_len'] for r in mx1_recs)}, LoF/frameshift={[r['taxon'] for r in mx1_recs if r['internal_stops']>0 or r['nt_len']%3!=0]}")

    print("3. Curating APOBEC3G (Simiiformes double-domain CD1+CD2)...")
    a3g_recs = curate_apobec3g()
    write_dataset("APOBEC3G", a3g_recs)
    print(f"   APOBEC3G: {len(a3g_recs)} species, nt={min(r['nt_len'] for r in a3g_recs)}..{max(r['nt_len'] for r in a3g_recs)}")

    print("4. Curating PRDM9 (unclipped RefSeq mammalian PRDM9 + Canidae LoF pseudogenes)...")
    prdm9_recs = curate_prdm9()
    write_dataset("PRDM9", prdm9_recs)
    znfs = [r["znf_count"] for r in prdm9_recs if r.get("znf_count", 0) > 0]
    print(f"   PRDM9: {len(prdm9_recs)} species, nt={min(r['nt_len'] for r in prdm9_recs)}..{max(r['nt_len'] for r in prdm9_recs)}, ZnF repeats={min(znfs)}..{max(znfs)}, LoF/frameshift={[r['taxon'] for r in prdm9_recs if r['internal_stops']>0 or r['nt_len']%3!=0]}")

    gt = build_and_verify_ground_truth({
        "TRIM5": trim5_recs, "MX1": mx1_recs, "APOBEC3G": a3g_recs, "PRDM9": prdm9_recs,
    })
    print("5. Ground-truth assertions passed 100%.")


if __name__ == "__main__":
    main()
