#!/usr/bin/env python3
"""
Build exact atomic-resolution (heavy-atom Euclidean distance d_min <= 4.5 A)
host-partner contact maps and literature-verified disordered/functional loop
ground-truth annotations for all 25 Tier-1 loci.

Every PDB coordinate file is downloaded from RCSB (.cif.gz), every chain is mapped
to Homo_sapiens canonical codon coordinates (1..L_hum) via banded Smith-Waterman,
partner chains are restricted strictly to authentic viral polypeptides (or dsDNA/dsRNA
for PRDM9 5EGB and OAS1 4IG8), and minimum heavy-atom distances are computed
directly from atomic Cartesian coordinates.
Python 3.9 compatible.
"""

import gzip
import json
import math
import os
import re
import time
import urllib.request
from collections import defaultdict

import numpy as np

ROOT = "/Users/sergei/Projects/unalignable-arms-race"
RAW_DIR = os.path.join(ROOT, "data", "raw_cds")
GT_DIR = os.path.join(ROOT, "data", "ground_truth")
PDB_CACHE = os.path.join(GT_DIR, "pdb_cache")
os.makedirs(PDB_CACHE, exist_ok=True)

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

# Pre-registered representative structures for the 25 Tier-1 loci.
# Some loci have primary PDB + optional backup PDB or assembly-1 fallback checked automatically.
TIER1_STRUCTURE_CONFIG = {
    # Arm 1A: 14 Antiviral Effectors & Sensors (8 PDB complexes + 6 mutagenesis/loop maps)
    "ISG15":    {"arm": "Antiviral_Effector_Sensor", "pdbs": ["3SDL", "5TL6"], "partner_mode": "viral", "loops": {}},
    "APOBEC3G": {"arm": "Antiviral_Effector_Sensor", "pdbs": ["8CX0"],         "partner_mode": "viral",
                 "loops": {"CD1_Loop7_Vif": [122, 134], "CD1_alpha2_3": [95, 110]}},
    "APOBEC3H": {"arm": "Antiviral_Effector_Sensor", "pdbs": ["8FVI"],         "partner_mode": "viral", "loops": {}},
    "SAMHD1":   {"arm": "Antiviral_Effector_Sensor", "pdbs": ["4CC9"],         "partner_mode": "viral", "loops": {}},
    "STAT2":    {"arm": "Antiviral_Effector_Sensor", "pdbs": ["8T12"],         "partner_mode": "viral", "loops": {}},
    "IRF3":     {"arm": "Antiviral_Effector_Sensor", "pdbs": ["5JER"],         "partner_mode": "viral", "loops": {}},
    "TRIM25":   {"arm": "Antiviral_Effector_Sensor", "pdbs": ["5NT1"],         "partner_mode": "viral", "loops": {}},
    "OAS1":     {"arm": "Antiviral_Effector_Sensor", "pdbs": ["4IG8"],         "partner_mode": "nucleic_acid",
                 "loops": {"C_term_prenylation_tail": [360, 400]}},
    "BST2":     {"arm": "Antiviral_Effector_Sensor", "pdbs": [],             "partner_mode": "functional_only",
                 "loops": {"N_tail_Nef": [1, 22], "Transmembrane_Vpu": [22, 46]},
                 "functional_sites": [14, 15, 16, 17, 18, 22, 26, 30, 33, 34, 37, 41, 45],
                 "dms_ref": "McNatt et al. PLoS Pathog 2009 (TM Vpu) + Jia et al. 2009/2014 (N-tail Nef)"},
    "TRIM5":    {"arm": "Antiviral_Effector_Sensor", "pdbs": [],             "partner_mode": "functional_only",
                 "loops": {"v1_capsid_patch": [324, 345], "v2_loop": [383, 395], "v3_loop": [415, 428], "v4_loop": [480, 488]},
                 "functional_sites": [330, 332, 334, 335, 336, 337, 339, 340, 385, 389, 420],
                 "dms_ref": "Tenthorey et al. eLife 2020 + PDB 8F25"},
    "MX1":      {"arm": "Antiviral_Effector_Sensor", "pdbs": [],             "partner_mode": "functional_only",
                 "loops": {"N_term_IDR": [1, 42], "Loop_L4_NP": [533, 572]},
                 "functional_sites": [540, 543, 549, 559, 561, 562, 564, 566, 567],
                 "dms_ref": "Mitchell et al. Cell Host Microbe 2012 + PDB 3SZR"},
    "EIF2AK2":  {"arm": "Antiviral_Effector_Sensor", "pdbs": [],             "partner_mode": "functional_only",
                 "loops": {"Kinase_insert_L": [255, 278], "AlphaG_helix_K3L": [482, 502], "AlphaD_helix": [365, 385]},
                 "functional_sites": [376, 485, 486, 488, 489, 490, 492, 496, 499],
                 "dms_ref": "Dar et al. Cell 2005 PDB 2A1A + Elde et al. Nature 2009"},
    "MAVS":     {"arm": "Antiviral_Effector_Sensor", "pdbs": [],             "partner_mode": "functional_only",
                 "loops": {"Picorna_3ABC_region": [420, 440], "HCV_NS34A_region": [500, 514]},
                 "functional_sites": [427, 428, 429, 506, 507, 508, 509],
                 "dms_ref": "Patel et al. PLoS Biol 2012 HCV NS3/4A cleavage escape"},
    "NLRP1":    {"arm": "Antiviral_Effector_Sensor", "pdbs": [],             "partner_mode": "functional_only",
                 "loops": {"Disordered_3Cpro_tripwire": [100, 165], "FIIND_autoprocess": [1190, 1225]},
                 "functional_sites": [128, 129, 130, 131, 132],
                 "dms_ref": "Tsu et al. eLife 2021 + Chavarria-Smith et al. PLoS Biol 2016"},

    # Arm 1B: 5 Viral Entry Receptors
    "DPP4":     {"arm": "Viral_Entry_Receptor", "pdbs": ["4L72"], "partner_mode": "viral", "loops": {}},
    "ACE2":     {"arm": "Viral_Entry_Receptor", "pdbs": ["6M0J"], "partner_mode": "viral", "loops": {}},
    "CD4":      {"arm": "Viral_Entry_Receptor", "pdbs": ["1GC1"], "partner_mode": "viral", "loops": {}},
    "TFRC":     {"arm": "Viral_Entry_Receptor", "pdbs": ["3KAS"], "partner_mode": "viral", "loops": {}},
    "NPC1":     {"arm": "Viral_Entry_Receptor", "pdbs": ["5F1B"], "partner_mode": "viral", "loops": {"DomainC_Loop2": [495, 512]}},

    # Arm 1C: 1 Tandem-Repeat DNA Co-Crystal Benchmark
    "PRDM9":    {"arm": "Tandem_Repeat_Control", "pdbs": ["5EGB"], "partner_mode": "nucleic_acid", "loops": {}},

    # Arm 1D: 5 Matched Hijacked Cellular Host Hub Controls
    "PPIA":     {"arm": "Hijacked_Cellular_Hub", "pdbs": ["1AK4"],         "partner_mode": "viral", "loops": {}},
    "ELOB":     {"arm": "Hijacked_Cellular_Hub", "pdbs": ["2MA9", "4N9F"], "partner_mode": "viral", "loops": {}},
    "ANP32A":   {"arm": "Hijacked_Cellular_Hub", "pdbs": ["8RN0"],         "partner_mode": "viral", "loops": {}},
    "RAE1":     {"arm": "Hijacked_Cellular_Hub", "pdbs": ["7VPG"],         "partner_mode": "viral", "loops": {}},
    "YWHAZ":    {"arm": "Hijacked_Cellular_Hub", "pdbs": ["7ZIT"],         "partner_mode": "viral", "loops": {}},
}


def read_fasta_dict(path):
    d = {}
    h, s = None, []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    d[h.split()[0]] = "".join(s).upper()
                h, s = line[1:], []
            else:
                s.append(line)
    if h is not None:
        d[h.split()[0]] = "".join(s).upper()
    return d


def translate_nt(seq):
    return "".join(CODON_TABLE.get(seq[i:i + 3].upper(), "X") for i in range(0, len(seq) - 2, 3)).rstrip("*")


def is_viral_organism(sci_name):
    if not sci_name:
        return False
    s = sci_name.lower()
    if any(b in s for b in ("homo sapiens", "mus musculus", "bos taurus", "synthetic", "escherichia", "saccharomyces")):
        return False
    tokens = (
        "virus", "viral", "viridae", "hiv", "siv", "sars", "mers", "influenza",
        "ebola", "sudan", "marburg", "machupo", "mammarenavirus", "dengue",
        "herpes", "papillomavir", "coronavirus", "adenovir"
    )
    return any(t in s for t in tokens)


def query_rcsb_entries(pdb_ids):
    gql_url = "https://data.rcsb.org/graphql"
    q = """
    query($ids: [String!]!) {
      entries(entry_ids: $ids) {
        rcsb_id
        exptl { method }
        rcsb_entry_info { resolution_combined }
        polymer_entities {
          rcsb_polymer_entity_container_identifiers {
            entity_id
            asym_ids
            auth_asym_ids
            uniprot_ids
          }
          entity_poly {
            pdbx_seq_one_letter_code_can
            type
          }
          rcsb_entity_source_organism {
            ncbi_scientific_name
          }
        }
      }
    }
    """
    payload = json.dumps({"query": q, "variables": {"ids": sorted(set(pdb_ids))}}).encode("utf-8")
    req = urllib.request.Request(
        gql_url, data=payload,
        headers={"Content-Type": "application/json", "User-Agent": "HyPhy3-Screen/1.0"}
    )
    with urllib.request.urlopen(req, timeout=45) as r:
        resp = json.loads(r.read().decode("utf-8"))
    return {e["rcsb_id"]: e for e in (resp.get("data", {}).get("entries") or []) if e}


def download_cif(pdb_id, assembly=False):
    suffix = "-assembly1.cif.gz" if assembly else ".cif.gz"
    local_p = os.path.join(PDB_CACHE, "{0}{1}".format(pdb_id, suffix))
    if os.path.exists(local_p) and os.path.getsize(local_p) > 200:
        return local_p
    url = "https://files.rcsb.org/download/{0}{1}".format(pdb_id, suffix)
    req = urllib.request.Request(url, headers={"User-Agent": "HyPhy3-Screen/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        data = r.read()
    with open(local_p, "wb") as out:
        out.write(data)
    return local_p


def tokenize_mmcif_line(line):
    tokens = []
    n = len(line)
    i = 0
    while i < n:
        while i < n and line[i] in (" ", "\t"):
            i += 1
        if i >= n:
            break
        if line[i] in ("'", '"'):
            quote = line[i]
            i += 1
            start = i
            while i < n:
                if line[i] == quote and (i + 1 >= n or line[i + 1] in (" ", "\t", "\n")):
                    break
                i += 1
            tokens.append(line[start:i])
            i += 1
        else:
            start = i
            while i < n and line[i] not in (" ", "\t", "\n"):
                i += 1
            tokens.append(line[start:i])
    return tokens


def parse_cif_atom_coords(cif_gz_path, target_asym_ids):
    """
    Extract heavy-atom coordinates from Model 1 (or single model) of a .cif.gz file,
    returning: asym_id -> seq_id(int) -> list of (x, y, z).
    """
    coords = defaultdict(lambda: defaultdict(list))
    cols = []
    in_atom_loop = False
    reading_atom_rows = False
    idx_group = idx_asym = idx_seq = idx_x = idx_y = idx_z = idx_elem = idx_model = -1

    with gzip.open(cif_gz_path, "rt", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\r\n")
            if not line or line.startswith("#"):
                if reading_atom_rows:
                    break
                continue
            if line.startswith("loop_"):
                if reading_atom_rows:
                    break
                cols = []
                in_atom_loop = True
                continue
            if in_atom_loop and line.startswith("_atom_site."):
                cols.append(line.strip())
                continue
            if in_atom_loop and line.startswith("_"):
                # Some other loop!
                if cols:
                    in_atom_loop = False
                    cols = []
                continue
            if in_atom_loop and cols and not reading_atom_rows:
                col_idx = {c: k for k, c in enumerate(cols)}
                idx_group = col_idx.get("_atom_site.group_PDB", -1)
                idx_asym = col_idx.get("_atom_site.label_asym_id", -1)
                idx_seq = col_idx.get("_atom_site.label_seq_id", -1)
                idx_x = col_idx.get("_atom_site.Cartn_x", -1)
                idx_y = col_idx.get("_atom_site.Cartn_y", -1)
                idx_z = col_idx.get("_atom_site.Cartn_z", -1)
                idx_elem = col_idx.get("_atom_site.type_symbol", -1)
                idx_model = col_idx.get("_atom_site.pdbx_PDB_model_num", -1)
                if min(idx_group, idx_asym, idx_seq, idx_x, idx_y, idx_z) < 0:
                    in_atom_loop = False
                    cols = []
                    continue
                reading_atom_rows = True

            if reading_atom_rows:
                if line.startswith("loop_") or line.startswith("_") or line.startswith("data_"):
                    break
                toks = tokenize_mmcif_line(line)
                if len(toks) < len(cols):
                    continue
                asym = toks[idx_asym]
                if asym not in target_asym_ids:
                    continue
                if idx_model >= 0 and toks[idx_model] not in ("1", ".", "?"):
                    continue
                if toks[idx_group] not in ("ATOM", "HETATM"):
                    continue
                if idx_elem >= 0 and toks[idx_elem].upper() in ("H", "D"):
                    continue
                seq_s = toks[idx_seq]
                if not seq_s.isdigit():
                    continue
                try:
                    x = float(toks[idx_x])
                    y = float(toks[idx_y])
                    z = float(toks[idx_z])
                except ValueError:
                    continue
                coords[asym][int(seq_s)].append((x, y, z))

    return coords


def smith_waterman_map(pdb_seq, hum_seq):
    """
    Local Smith-Waterman alignment mapping 1-indexed PDB label_seq_id (1..len(pdb_seq))
    to 1-indexed Homo_sapiens codon/residue position (1..len(hum_seq)).
    Returns (dict pdb_1idx -> hum_1idx, identity_frac, aligned_len).
    """
    m, n = len(pdb_seq), len(hum_seq)
    if m == 0 or n == 0:
        return {}, 0.0, 0

    MATCH = 3
    MISMATCH = -2
    GAP = -4

    H = np.zeros((m + 1, n + 1), dtype=np.int32)
    B = np.zeros((m + 1, n + 1), dtype=np.int8)  # 0=stop, 1=diag, 2=up(gap hum), 3=left(gap pdb)

    best_score = 0
    best_i, best_j = 0, 0
    for i in range(1, m + 1):
        ci = pdb_seq[i - 1]
        hi = H[i]
        hi_prev = H[i - 1]
        bi = B[i]
        for j in range(1, n + 1):
            sc_diag = hi_prev[j - 1] + (MATCH if ci == hum_seq[j - 1] else MISMATCH)
            sc_up = hi_prev[j] + GAP
            sc_left = hi[j - 1] + GAP
            best = sc_diag
            bdir = 1
            if sc_up > best:
                best = sc_up
                bdir = 2
            if sc_left > best:
                best = sc_left
                bdir = 3
            if best <= 0:
                hi[j] = 0
                bi[j] = 0
            else:
                hi[j] = best
                bi[j] = bdir
                if best > best_score:
                    best_score = best
                    best_i, best_j = i, j

    mapping = {}
    matches = 0
    aligned = 0
    i, j = best_i, best_j
    while i > 0 and j > 0 and B[i, j] != 0:
        d = B[i, j]
        if d == 1:
            mapping[i] = j
            aligned += 1
            if pdb_seq[i - 1] == hum_seq[j - 1]:
                matches += 1
            i -= 1
            j -= 1
        elif d == 2:
            i -= 1
        else:
            j -= 1

    ident = float(matches) / max(1, aligned)
    return mapping, ident, aligned


def min_heavy_atom_distances(host_res_coords, partner_atom_array):
    """
    Given host_res_coords: dict hum_1idx -> list of (x,y,z) heavy atoms,
    and partner_atom_array: Nx3 float numpy array of partner heavy atoms,
    return dict hum_1idx -> minimum heavy-atom Euclidean distance (Angstroms).
    """
    if partner_atom_array.shape[0] == 0:
        return {}
    out = {}
    for hpos, atom_list in host_res_coords.items():
        if not atom_list:
            continue
        A = np.asarray(atom_list, dtype=np.float64)  # Kx3
        # Compute K x N squared Euclidean distances
        diff = A[:, None, :] - partner_atom_array[None, :, :]
        d2 = np.sum(diff * diff, axis=2)
        dmin = float(math.sqrt(float(np.min(d2))))
        out[hpos] = dmin
    return out


def compute_structure_contacts(gene, up_acc, hum_aa, cfg, rcsb_entries):
    partner_mode = cfg["partner_mode"]
    combined_dmin = {}  # hum_1idx -> best (smallest) distance across selected PDB(s)
    pdb_summaries = []

    for pdb_id in cfg["pdbs"]:
        ent = rcsb_entries.get(pdb_id)
        if not ent:
            continue
        host_chains = []      # (asym_id, seq_can, mapping_dict, ident)
        partner_chains = []   # (asym_id, desc, organism)
        hum_self_seq_fragments = set()

        for pe in (ent.get("polymer_entities") or []):
            ids_block = pe.get("rcsb_polymer_entity_container_identifiers") or {}
            asym_ids = ids_block.get("asym_ids") or []
            uids = ids_block.get("uniprot_ids") or []
            ep = pe.get("entity_poly") or {}
            seq_can = ep.get("pdbx_seq_one_letter_code_can") or ""
            ptype = (ep.get("type") or "").lower()
            orgs = [o.get("ncbi_scientific_name") or "" for o in (pe.get("rcsb_entity_source_organism") or [])]
            org_str = ", ".join(orgs)

            # Check if this entity matches our target human protein
            matched_host = False
            if ("polypeptide" in ptype) and (up_acc in uids or (not uids and len(seq_can) >= 40)):
                sw_map, ident, aln_len = smith_waterman_map(seq_can, hum_aa)
                if aln_len >= 35 and ident >= 0.72:
                    matched_host = True
                    hum_self_seq_fragments.add(seq_can)
                    for aid in asym_ids:
                        host_chains.append((aid, seq_can, sw_map, ident, aln_len))

            if matched_host:
                # Could be a single-chain fusion? Only if non-aligned portion is >=30 aa viral, which is rare in our chosen crystals.
                continue

            if partner_mode == "nucleic_acid":
                if "ribonucleotide" in ptype or "deoxyribonucleotide" in ptype:
                    for aid in asym_ids:
                        partner_chains.append((aid, ptype, org_str or "Nucleic_acid"))
            elif partner_mode == "viral":
                if any(is_viral_organism(o) for o in orgs):
                    for aid in asym_ids:
                        partner_chains.append((aid, ptype, org_str))

        if not host_chains or not partner_chains:
            print("    WARNING: {0} {1} host_chains={2} partner_chains={3}".format(
                gene, pdb_id, len(host_chains), len(partner_chains)), flush=True)
            continue

        # Try asymmetric unit first; fall back to assembly1 if 0 direct contacts
        all_target_asyms = set(h[0] for h in host_chains) | set(p[0] for p in partner_chains)
        cif_path = download_cif(pdb_id, assembly=False)
        coords = parse_cif_atom_coords(cif_path, all_target_asyms)

        # Evaluate best asymmetric unit contacts across all partner chains
        partner_pts_list = []
        for aid, _pt, _org in partner_chains:
            for res_atoms in coords.get(aid, {}).values():
                partner_pts_list.extend(res_atoms)

        def eval_contacts_for_coords(c_dict, p_pts):
            if not p_pts:
                return {}, {}
            P = np.asarray(p_pts, dtype=np.float64)
            # For homomeric multiple host copies in AU, compute min distance per hum_1idx across copies
            h_res_atoms = defaultdict(list)
            best_map = {}
            resolved_h = set()
            for aid, _sc, sw_map, _ident, _al in host_chains:
                for pdb_seq_i, atom_list in c_dict.get(aid, {}).items():
                    if pdb_seq_i in sw_map:
                        hpos = sw_map[pdb_seq_i]
                        resolved_h.add(hpos)
                        h_res_atoms[hpos].extend(atom_list)
            dists = min_heavy_atom_distances(h_res_atoms, P)
            for hpos in resolved_h:
                if hpos in dists:
                    best_map[hpos] = dists[hpos]
            return best_map, resolved_h

        d_au, res_au = eval_contacts_for_coords(coords, partner_pts_list)
        n45 = sum(1 for v in d_au.values() if v <= 4.5)
        used_assembly = False
        if n45 < 4:
            # Check assembly1 (e.g. symmetric oligomer interface)
            try:
                cif_asm = download_cif(pdb_id, assembly=True)
                coords_asm = parse_cif_atom_coords(cif_asm, all_target_asyms)
                p_pts_asm = []
                for aid, _pt, _org in partner_chains:
                    for res_atoms in coords_asm.get(aid, {}).values():
                        p_pts_asm.extend(res_atoms)
                d_asm, res_asm = eval_contacts_for_coords(coords_asm, p_pts_asm)
                if sum(1 for v in d_asm.values() if v <= 4.5) > n45:
                    d_au, res_au = d_asm, res_asm
                    n45 = sum(1 for v in d_au.values() if v <= 4.5)
                    used_assembly = True
            except Exception as exc:
                print("    assembly1 fallback skipped for {0} {1}: {2}".format(gene, pdb_id, exc))

        for hpos, dval in d_au.items():
            if hpos not in combined_dmin or dval < combined_dmin[hpos]:
                combined_dmin[hpos] = dval

        res_list = (ent.get("rcsb_entry_info") or {}).get("resolution_combined") or [None]
        pdb_summaries.append({
            "pdb_id": pdb_id,
            "resolution_A": res_list[0],
            "used_assembly1": used_assembly,
            "host_asyms": sorted(set(h[0] for h in host_chains)),
            "host_ident": round(max(h[3] for h in host_chains), 3),
            "partner_asyms": sorted(set(p[0] for p in partner_chains)),
            "partner_organisms": sorted(set(p[2] for p in partner_chains)),
            "resolved_human_residues": len(res_au),
            "contacts_4p5A_this_pdb": n45,
            "contacts_8p0A_this_pdb": sum(1 for v in d_au.values() if v <= 8.0),
        })

    return combined_dmin, pdb_summaries


def main():
    canon_path = os.path.join(GT_DIR, "uniprot_25gene_canonical.json")
    uniprot_map = json.load(open(canon_path))

    # Collect all PDB IDs needed
    all_pdbs = []
    for g, cfg in TIER1_STRUCTURE_CONFIG.items():
        all_pdbs.extend(cfg["pdbs"])
    print("Querying RCSB GraphQL for {0} representative complexes...".format(len(all_pdbs)))
    rcsb_entries = query_rcsb_entries(all_pdbs)
    print("Retrieved {0} RCSB entries.".format(len(rcsb_entries)))

    gt_out = {}
    print("\n{0:10s} | {1:25s} | {2:12s} | {3:6s} | {4:8s} | {5:8s} | {6:8s} | {7}".format(
        "Gene", "Arm", "PDB(s)", "HumAA", "Resolved", "<=4.5A", "<=8.0A", "Partner"
    ))
    print("-" * 112)

    for gene, cfg in TIER1_STRUCTURE_CONFIG.items():
        fa_path = os.path.join(RAW_DIR, "{0}_unaligned.fa".format(gene))
        seqs = read_fasta_dict(fa_path)
        assert "Homo_sapiens" in seqs, "Missing Homo_sapiens in " + fa_path
        hum_aa = translate_nt(seqs["Homo_sapiens"])
        L_hum = len(hum_aa)
        up_acc = uniprot_map[gene]["uniprot"]

        dmin_map, pdb_sums = compute_structure_contacts(gene, up_acc, hum_aa, cfg, rcsb_entries)
        c45 = sorted(h for h, d in dmin_map.items() if d <= 4.5)
        c80 = sorted(h for h, d in dmin_map.items() if d <= 8.0)
        resolved = sorted(dmin_map.keys())
        scaffold_gt8 = sorted(h for h, d in dmin_map.items() if d > 8.0)
        unmodeled = sorted(set(range(1, L_hum + 1)) - set(resolved))

        # Special verified C2H2 annotation for PRDM9 DNA contacts & scaffold
        extra = {}
        if gene == "PRDM9":
            znf_matches = list(re.finditer(r"C..C.{12}H...H", hum_aa))
            dna_contacts, c2h2_scaff, reps = [], [], []
            for idx, m in enumerate(znf_matches, start=1):
                c1 = m.start() + 1
                c2 = c1 + 3
                h1 = c2 + 13
                h2 = h1 + 4
                contacts = [c2 + 6, c2 + 8, c2 + 9, c2 + 12]
                dna_contacts.extend(contacts)
                c2h2_scaff.extend([c1, c2, h1, h2])
                reps.append({"finger": idx, "c1_h2": [c1, h2], "dna_contacts": contacts, "c2h2": [c1, c2, h1, h2]})
            extra = {
                "znf_repeats": reps,
                "znf_dna_contact_1idx": sorted(set(dna_contacts)),
                "znf_c2h2_scaffold_1idx": sorted(set(c2h2_scaff)),
            }

        gt_out[gene] = {
            "gene": gene,
            "uniprot": up_acc,
            "arm": cfg["arm"],
            "partner_mode": cfg["partner_mode"],
            "human_ref_aa_len": L_hum,
            "human_ref_seq": hum_aa,
            "pdb_summaries": pdb_sums,
            "resolved_residues_1idx": resolved,
            "contact_4p5A_1idx": c45,
            "proximal_4p5_8p0A_1idx": sorted(h for h, d in dmin_map.items() if 4.5 < d <= 8.0),
            "contact_8p0A_1idx": c80,
            "resolved_scaffold_gt8A_1idx": scaffold_gt8,
            "unmodeled_disordered_1idx": unmodeled,
            "dmin_by_human_pos": {str(k): round(v, 3) for k, v in sorted(dmin_map.items())},
            "functional_loops_1idx": cfg.get("loops", {}),
            "literature_functional_sites_1idx": cfg.get("functional_sites", []),
            "dms_or_functional_ref": cfg.get("dms_ref", ""),
            "extra": extra,
        }

        pdb_lbl = ",".join(s["pdb_id"] for s in pdb_sums) if pdb_sums else "DMS/Loop"
        partner_lbl = ", ".join(sorted(set(o for s in pdb_sums for o in s["partner_organisms"])))[:34] if pdb_sums else cfg.get("dms_ref", "")[:34]
        print("{0:10s} | {1:25s} | {2:12s} | {3:6d} | {4:8d} | {5:8d} | {6:8d} | {7}".format(
            gene, cfg["arm"], pdb_lbl, L_hum, len(resolved), len(c45), len(c80), partner_lbl
        ), flush=True)

    out_path = os.path.join(GT_DIR, "tier1_25gene_ground_truth.json")
    with open(out_path, "w", encoding="utf-8") as out:
        json.dump(gt_out, out, indent=2)
    print("\nWrote 25-gene atomic ground-truth map to:", out_path)


if __name__ == "__main__":
    main()
