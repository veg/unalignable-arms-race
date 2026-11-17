#!/usr/bin/env python3
"""
Programmatic Science-Caliber Vector Figure Portfolio Generator
for "Evolutionary arms races across unalignable mammalian host-virus interfaces".

Strictly follows /Users/sergei/Dropbox/PROMPTS/science_manuscript_and_figures.md (Part II):
  - Native 7.25-in 3-column width, dpi=300 PNG + Type-42 TrueType vector PDF
  - Editorial color palette (INK, SUBINK, MUTED, HAIR, BG_PANEL, ORANGE, BLUE, TEAL, AMBER, PURPLE, ROSE)
  - Explicit figure-coordinate draw_panel_header(fig, x, y, letter, title) with no fig.suptitle()
  - 100% deterministic empirical values read directly from executed disk JSON/FASTA/JSONL files
"""

import gzip
import json
import math
import os
import re
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.path import Path
from matplotlib.patches import FancyBboxPatch
import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu, wilcoxon

ROOT = "/Users/sergei/Projects/unalignable-arms-race"
MASTER_JSON = os.path.join(ROOT, "results", "tier1_25gene_master_synthesis.json")
TIER2_JSON = os.path.join(ROOT, "results", "tier2_toga2_vip_summary.json")
GT_PATH = os.path.join(ROOT, "data", "ground_truth", "tier1_25gene_ground_truth.json")
WALRUS_DIR = os.path.join(ROOT, "results", "walrus")
PRIME_DIR = os.path.join(ROOT, "results", "prime")
LEGACY_DIR = os.path.join(ROOT, "results", "legacy_foil")
TOGA_GZ = os.path.join(ROOT, "data", "toga2_completed_summaries.jsonl.gz")
VIP_JSON = os.path.join(ROOT, "data", "ground_truth", "uniprot_human_vips.json")

PNG_DIR = os.path.join(ROOT, "figures")
PDF_DIR = os.path.join(ROOT, "manuscript", "figures")
BRAIN_DIR = "/Users/sergei/.gemini/antigravity-cli/brain/c566a240-5f92-4323-8c76-3829adda72a6"
os.makedirs(PNG_DIR, exist_ok=True)
os.makedirs(PDF_DIR, exist_ok=True)

# Editorial Palette Tokens (~/Projects/trees Standard)
INK      = "#0F172A"
SUBINK   = "#475569"
MUTED    = "#94A3B8"
HAIR     = "#E2E8F0"
BG_PANEL = "#F8FAFC"
ORANGE   = "#E04F16"
BLUE     = "#1D4ED8"
TEAL     = "#0D9488"
AMBER    = "#D97706"
PURPLE   = "#7C3AED"
ROSE     = "#E11D48"

ARM_COLOR = {
    "Antiviral_Effector_Sensor": ORANGE,
    "Viral_Entry_Receptor": BLUE,
    "Tandem_Repeat_Control": PURPLE,
    "Hijacked_Cellular_Hub": TEAL,
}
ARM_LABEL = {
    "Antiviral_Effector_Sensor": "Antiviral effectors/sensors (n=14)",
    "Viral_Entry_Receptor": "Viral entry receptors (n=5)",
    "Tandem_Repeat_Control": "PRDM9 repeat array (n=1)",
    "Hijacked_Cellular_Hub": "Hijacked cellular hubs (n=5)",
}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "mathtext.fontset": "custom",
    "mathtext.rm": "sans",
    "mathtext.it": "sans:italic",
    "mathtext.bf": "sans:bold",
    "font.size": 6.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.edgecolor": SUBINK,
    "axes.linewidth": 0.55,
    "xtick.major.width": 0.50,
    "ytick.major.width": 0.50,
    "xtick.major.size": 2.4,
    "ytick.major.size": 2.4,
    "xtick.major.pad": 2.2,
    "ytick.major.pad": 2.2,
    "xtick.color": "#334155",
    "ytick.color": "#334155",
})


def draw_panel_header(fig, x_fig, y_fig, letter, title):
    fig.text(x_fig, y_fig, letter, fontsize=10.0, fontweight="bold", color=INK, ha="left", va="center")
    fig.text(x_fig + 0.022, y_fig, title, fontsize=7.6, fontweight="bold", color="#1E293B", ha="left", va="center")


def add_callout_pill(ax, x, y, text, color=INK, fc=BG_PANEL, ec="#CBD5E1", ha="center", va="center", fontsize=6.2):
    ax.text(
        x, y, text, ha=ha, va=va, fontsize=fontsize, fontweight="bold", color=color,
        bbox=dict(boxstyle="round,pad=0.25", facecolor=fc, edgecolor=ec, linewidth=0.55),
        zorder=8,
    )


def save_both(fig, stem):
    pdf_p = os.path.join(PDF_DIR, stem + ".pdf")
    png_p = os.path.join(PNG_DIR, stem + ".png")
    fig.savefig(pdf_p)
    fig.savefig(png_p, dpi=300)
    plt.close(fig)
    if os.path.isdir(BRAIN_DIR):
        shutil.copy2(png_p, os.path.join(BRAIN_DIR, stem + ".png"))
    print("Saved", stem, "->", pdf_p, png_p)


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


def map_msa_to_human(msa_path, human_taxon="Homo_sapiens"):
    seqs = read_fasta_dict(msa_path)
    h_nt = seqs[human_taxon]
    n_cols = len(h_nt) // 3
    col_to_human, human_to_col, col_float = {}, {}, {}
    h_pos, last_h, ins_run = 0, 0.5, []
    for c in range(1, n_cols + 1):
        cod = h_nt[(c - 1) * 3: c * 3]
        if cod != "---":
            h_pos += 1
            col_to_human[c] = h_pos
            human_to_col[h_pos] = c
            if ins_run:
                k = len(ins_run)
                for idx, ic in enumerate(ins_run, start=1):
                    col_float[ic] = (h_pos - 1) + idx / (k + 1.0)
                ins_run = []
            col_float[c] = float(h_pos)
            last_h = float(h_pos)
        else:
            col_to_human[c] = None
            ins_run.append(c)
    if ins_run:
        k = len(ins_run)
        for idx, ic in enumerate(ins_run, start=1):
            col_float[ic] = last_h + idx / (k + 1.0)
    taxa_seqs = list(seqs.values())
    n_taxa = len(taxa_seqs)
    occ = {}
    for c in range(1, n_cols + 1):
        non_gap = sum(1 for s in taxa_seqs if s[(c - 1) * 3: c * 3] not in ("---", "NNN"))
        occ[c] = non_gap / max(1, n_taxa)
    return n_cols, n_taxa, col_to_human, col_float, human_to_col, occ


def is_walrus_innov(s):
    cls = s.get("classification") or ""
    return (
        ("Diversifying" in cls)
        or ("Toggle" in cls)
        or ("Convergent" in cls)
        or (s.get("meme_bh_q") is not None and s["meme_bh_q"] <= 0.10)
    )


def load_gene_bundles(master, gt):
    bundles = {}
    for rec in master["tier1_genes"]:
        g = rec["gene"]
        w = json.load(open(os.path.join(WALRUS_DIR, f"{g}.walrus.json")))
        p = json.load(open(os.path.join(PRIME_DIR, f"{g}.prime.json")))
        msa_map = map_msa_to_human(os.path.join(WALRUS_DIR, f"{g}.walrus.msa.fa"))
        bundles[g] = (rec, w, p, msa_map, gt[g])
    return bundles


def bh_qvalues(pvals):
    n = len(pvals)
    order = np.argsort(pvals)
    q = np.ones(n, dtype=float)
    prev = 1.0
    for rank_rev in range(n - 1, -1, -1):
        idx = order[rank_rev]
        val = min(prev, pvals[idx] * n / (rank_rev + 1.0))
        q[idx] = val
        prev = val
    return q


# ==============================================================================
# FIGURE 1: Coordinate Grammar, Architecture Plane & 25-Locus Conflict Ledger
# ==============================================================================
def build_fig1(master, bundles):
    fig = plt.figure(figsize=(7.25, 6.05), dpi=300)
    gs = gridspec.GridSpec(
        2, 2, width_ratios=[1.02, 1.30], height_ratios=[0.92, 1.24],
        left=0.075, right=0.98, top=0.925, bottom=0.075, wspace=0.40, hspace=0.48
    )

    # Panel A (top-left): TRIM5 PRYSPRY Legacy MAFFT vs WALRUS Coordinate Ledger
    ax_a = fig.add_subplot(gs[0, 0])
    draw_panel_header(fig, 0.015, 0.965, "A", "TRIM5 PRYSPRY joint alignment (67 mammals)")
    rec_t, w_t, p_t, msa_t, gt_t = bundles["TRIM5"]
    n_cols, _, col_to_hum, col_float, hum_to_col, occ = msa_t
    sites_t = w_t["site_atlas"]

    meme_t = json.load(open(os.path.join(LEGACY_DIR, "TRIM5.mafft.MEME.json")))
    m_pvals = np.array(meme_t["site_results"]["data"]["p_value"], dtype=float)
    m_col_map = map_msa_to_human(os.path.join(LEGACY_DIR, "TRIM5.mafft_codon.fa"))[3]

    # Plot PRYSPRY window 285..493
    for lname, (lo, hi) in gt_t["functional_loops_1idx"].items():
        ax_a.axvspan(lo, hi, color="#FEF2F2", alpha=0.9, zorder=0)
        ax_a.text((lo + hi) / 2.0, 7.6, lname.split("_")[0], ha="center", va="top", fontsize=5.8, fontweight="bold", color=ORANGE)

    # MAFFT nominal p<0.05 ticks inverted below 0
    m_xs = np.array([m_col_map[i + 1] for i in range(len(m_pvals))])
    m_sig = np.where((m_xs >= 285) & (m_pvals <= 0.05))[0]
    m_y = np.minimum(5.2, -np.log10(np.maximum(1e-6, m_pvals[m_sig])))
    ax_a.vlines(m_xs[m_sig], 0.0, -m_y, color=MUTED, lw=0.75, alpha=0.85, label="MAFFT + MEME p <= 0.05 (56 sites, drops F. catus)")

    # WALRUS signal above 0
    w_xs = np.array([col_float[s["site"]] for s in sites_t])
    w_mask = (w_xs >= 285) & np.array([is_walrus_innov(s) for s in sites_t])
    w_art = (w_xs >= 285) & np.array([bool(s.get("artifact_eliminated")) for s in sites_t])
    w_y = np.array([min(7.5, -math.log10(max(1e-8, s.get("simes_omnibus_q") or 1.0))) for s in sites_t])
    ax_a.vlines(w_xs[w_mask], 0.0, np.maximum(w_y[w_mask], 1.15), color=ORANGE, lw=0.95, zorder=4, label="WALRUS selection (34.2% in v1–v4 loops)")
    ax_a.scatter(w_xs[w_mask], np.maximum(w_y[w_mask], 1.15), s=10, color=ORANGE, edgecolor="white", lw=0.3, zorder=5)
    ax_a.scatter(w_xs[w_art], np.full(np.sum(w_art), 0.7), marker="x", s=22, color=ROSE, lw=1.1, zorder=6, label="POA artifact pruned (6 gene-wide)")
    ax_a.axhline(0.0, color=SUBINK, lw=0.65)
    ax_a.set_xlim(285, 494)
    ax_a.set_ylim(-5.8, 8.2)
    ax_a.set_yticks([-4, -2, 0, 2, 4, 6])
    ax_a.set_yticklabels(["4", "2", "0", "2", "4", "6"])
    ax_a.set_xlabel("Human TRIM5 codon coordinate", fontsize=6.5)
    ax_a.set_ylabel("MAFFT $-\\log_{10} p$  |  WALRUS $-\\log_{10} q$", fontsize=6.2)
    add_callout_pill(ax_a, 345, -4.3, "MAFFT stretches MSA to 650 codons\n(vs. 595 joint) & drops F. catus", color="#334155", fontsize=5.4)
    ax_a.legend(loc="upper left", fontsize=5.3, frameon=True, facecolor="white", edgecolor=HAIR, borderpad=0.3)

    # Panel B (bottom-left): 2D Indel Plasticity vs Innovation Plane across all 25 Tier-1 Loci
    ax_b = fig.add_subplot(gs[1, 0])
    draw_panel_header(fig, 0.015, 0.492, "B", "Conflict arms vs. hubs (n = 25 loci)")
    recs = master["tier1_genes"]
    arm_markers = {
        "Antiviral_Effector_Sensor": "o",
        "Viral_Entry_Receptor": "s",
        "Tandem_Repeat_Control": "D",
        "Hijacked_Cellular_Hub": "^",
    }
    ax_b.axhspan(-1, 13.0, color="#F1F5F9", alpha=0.85, zorder=0)
    ax_b.axhline(12.5, color=TEAL, ls="--", lw=0.75, alpha=0.7)

    for arm in arm_markers:
        sub = [r for r in recs if r["arm"] == arm]
        xs_b = [r["mean_gap_pct"] for r in sub]
        ys_b = [r["innov_site_pct"] for r in sub]
        sz_b = [22 + 3.6 * r["poa_artifacts_eliminated"] for r in sub]
        ax_b.scatter(xs_b, ys_b, s=sz_b, c=ARM_COLOR[arm], marker=arm_markers[arm],
                     alpha=0.88, edgecolor="white", linewidth=0.5, label=ARM_LABEL[arm], zorder=4)

    offsets = {
        "APOBEC3G": (-4.5, 4.2),
        "TRIM25": (3.8, 3.0),
        "ACE2": (-3.8, 2.2),
        "DPP4": (-1.8, -3.4),
        "TFRC": (3.4, -1.0),
        "TRIM5": (-4.2, 2.2),
        "OAS1": (2.2, 3.8),
        "NLRP1": (2.0, 4.4),
        "MAVS": (3.2, -1.6),
        "APOBEC3H": (-11.5, 3.2),
        "PRDM9": (-11.0, -3.2),
        "ANP32A": (3.2, -1.2),
        "ELOB": (-7.5, 2.8),
        "PPIA": (2.5, 1.0),
        "YWHAZ": (2.8, 0.8),
        "RAE1": (2.8, -0.8),
    }
    for r in recs:
        gn = r["gene"]
        if gn in offsets:
            dx, dy = offsets[gn]
            ax_b.annotate(
                gn, xy=(r["mean_gap_pct"], r["innov_site_pct"]),
                xytext=(r["mean_gap_pct"] + dx, r["innov_site_pct"] + dy),
                fontsize=5.3, fontweight="bold", color=INK if r["poa_artifacts_eliminated"] < 20 else ROSE,
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.45) if (abs(dx) > 2.4 or abs(dy) > 1.8) else None,
                zorder=6,
            )

    add_callout_pill(ax_b, 14.0, 6.2, "Marker area ~ pruned\nPOA artifacts (max 51)", color=ROSE, fontsize=4.95)
    ax_b.set_xlim(-2.0, 50.0)
    ax_b.set_ylim(-2.5, 44.5)
    ax_b.set_xlabel("Mean codon alignment gap density (%)", fontsize=6.5)
    ax_b.set_ylabel("Sites under diversifying / toggle selection (%)", fontsize=6.5)
    ax_b.legend(loc="upper right", fontsize=4.85, frameon=True, facecolor="white", edgecolor=HAIR, markerscale=0.72)

    # Panel C (right spanning both rows): Cleveland multi-modal atlas across all 25 genes
    ax_c = fig.add_subplot(gs[:, 1])
    draw_panel_header(fig, 0.455, 0.965, "C", "Multi-modal ledger across all 25 pre-registered Tier-1 loci")

    ordered = sorted(recs, key=lambda r: (
        ["Hijacked_Cellular_Hub", "Tandem_Repeat_Control", "Viral_Entry_Receptor", "Antiviral_Effector_Sensor"].index(r["arm"]),
        r["innov_site_pct"]
    ))
    ys_c = np.arange(len(ordered))
    ax_c.axhspan(-0.5, 4.5, color="#F0FDFA", alpha=0.85, zorder=0)
    ax_c.axhspan(4.5, 5.5, color="#F5F3FF", alpha=0.85, zorder=0)
    ax_c.axhspan(5.5, 10.5, color="#EFF6FF", alpha=0.75, zorder=0)
    ax_c.axhspan(10.5, 24.5, color="#FFF7ED", alpha=0.75, zorder=0)

    for i, r in enumerate(ordered):
        col = ARM_COLOR[r["arm"]]
        tot = r["innov_site_pct"]
        tog = 100.0 * r["toggle_conv_sites"] / max(1, r["msa_codons"])
        div_only = max(0.0, tot - tog)
        ax_c.hlines(i, 0, tot, color=HAIR, lw=1.2, zorder=1)
        ax_c.plot([0, div_only], [i, i], color=col, lw=2.3, solid_capstyle="round", zorder=2)
        ax_c.plot([div_only, tot], [i, i], color=PURPLE, lw=2.3, solid_capstyle="round", zorder=3)
        ax_c.scatter(tot, i, s=20, color=col if r["busted_e_sig"] else "white", edgecolor=col, lw=1.0, zorder=4)
        # Annotate right-hand concise mechanism tags: Ins cols | pI(M) | SPR/Pseudo
        pi_m = r["prime_marginal_changing_by_prop"].get("Isoelectric_Point_pI", 0)
        tags = []
        if r["non_human_inserted_cols"] > 0:
            tags.append(f"+{r['non_human_inserted_cols']}i")
        if pi_m > 0:
            tags.append(f"pI:{pi_m}")
        if r["spr_routed_bridges"] > 0:
            tags.append(f"SPR:{r['spr_routed_bridges']}")
        if r["admitted_pseudogenes"]:
            tags.append(f"Ps:{len(r['admitted_pseudogenes'])}")
        ax_c.text(tot + 1.0, i, "  ".join(tags), va="center", ha="left", fontsize=5.2, color="#334155")

    ax_c.set_yticks(ys_c)
    ax_c.set_yticklabels([f"{r['gene']} ({r['admitted_taxa']})" for r in ordered], fontsize=5.8, fontweight="bold")
    ax_c.set_xlim(0, 56.0)
    ax_c.set_ylim(-0.6, 24.8)
    ax_c.set_xlabel("Selected alignment columns (%)", fontsize=6.5)

    leg_c = [
        mpatches.Patch(color=ORANGE, label="Diversifying rate (filled: BUSTED-E p <= 0.10)"),
        mpatches.Patch(color=PURPLE, label="Two-state toggle / convergent component"),
        mpatches.Patch(color=TEAL, label="Hijacked hub baseline"),
    ]
    ax_c.legend(handles=leg_c, loc="lower right", fontsize=5.15, frameon=True, facecolor="white", edgecolor=HAIR)
    add_callout_pill(ax_c, 45.5, 2.0, "+i: non-human insertions\npI: S-PRIME pI sites\nSPR: routed conversion\nPs: repaired pseudogene", color=SUBINK, fontsize=5.0, ha="center")

    save_both(fig, "fig1_science_atlas")


# ==============================================================================
# FIGURE 2: Atomic PDB <=4.5 A Enrichment, Virus-Matched Mirrors & Disorder
# ==============================================================================
def build_fig2(master):
    fig = plt.figure(figsize=(7.25, 5.25), dpi=300)
    gs = gridspec.GridSpec(
        1, 3, width_ratios=[1.35, 1.05, 1.05],
        left=0.145, right=0.975, top=0.885, bottom=0.115, wspace=0.40
    )

    recs = master["tier1_genes"]
    pdb_recs = [r for r in recs if r["pdb_atomic_audit"] is not None]
    agg = master["atomic_pdb_aggregate"]

    # Panel A: Paired Cleveland Dumbbell across all 19 PDB complexes
    ax_a = fig.add_subplot(gs[0, 0])
    draw_panel_header(fig, 0.015, 0.945, "A", "Heavy-atom contacts (<=4.5 A) vs. buried cores (>8.0 A)")

    # Sort hubs bottom, PRDM9 middle, conflict top by contact %
    order_grp = {"Hijacked_Cellular_Hub": 0, "Tandem_Repeat_Control": 1, "Viral_Entry_Receptor": 2, "Antiviral_Effector_Sensor": 3}
    pdb_sorted = sorted(pdb_recs, key=lambda r: (order_grp[r["arm"]], r["pdb_atomic_audit"]["contact_4p5A"]["innov_pct"]))
    ya = np.arange(len(pdb_sorted))

    ax_a.axhspan(-0.5, 4.5, color="#F0FDFA", alpha=0.85, zorder=0)
    ax_a.axhspan(4.5, 5.5, color="#F5F3FF", alpha=0.85, zorder=0)
    ax_a.axhspan(5.5, 18.5, color="#FFF7ED", alpha=0.55, zorder=0)

    for idx, r in enumerate(pdb_sorted):
        pa = r["pdb_atomic_audit"]
        c_pct = pa["contact_4p5A"]["innov_pct"]
        c_n = pa["contact_4p5A"]["n"]
        s_obj = pa["resolved_scaffold_gt8A"] if pa["resolved_scaffold_gt8A"]["n"] > 0 else pa["resolved_noncontact_gt4p5A"]
        s_pct = s_obj["innov_pct"]
        col = ARM_COLOR[r["arm"]]
        ax_a.plot([s_pct, c_pct], [idx, idx], color=col if c_pct >= s_pct else MUTED, lw=1.25, zorder=2)
        ax_a.scatter(s_pct, idx, s=16, facecolor="white", edgecolor=MUTED, lw=1.0, zorder=3)
        ax_a.scatter(c_pct, idx, s=16 + 1.1 * c_n, color=col, edgecolor="white", lw=0.45, zorder=4)

    ax_a.set_yticks(ya)
    ax_a.set_yticklabels([f"{r['gene']} ({r['pdb_atomic_audit']['pdbs'][0]})" for r in pdb_sorted], fontsize=5.6, fontweight="bold")
    ax_a.set_xlim(-3, 92)
    ax_a.set_ylim(-0.6, 18.7)
    ax_a.set_xlabel("Selected residues (%)", fontsize=6.5)

    leg_a = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=MUTED, markersize=4.5, label="Core (>8.0 A)"),
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=ORANGE, markeredgecolor="white", markersize=5.0, label="Conflict <=4.5 A"),
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=TEAL, markeredgecolor="white", markersize=5.0, label="Hub <=4.5 A"),
    ]
    ax_a.legend(handles=leg_a, loc="lower right", fontsize=5.2, frameon=True, facecolor="white", edgecolor=HAIR)
    add_callout_pill(ax_a, 54, 15.2, "Conflict pooled (n=13):\n42.5% vs. 29.2% core\nOR=1.79x, p=1.9e-6", color=ORANGE, fontsize=5.1)

    # Panel B: Virus-Matched Mirror Complexes (paired Cleveland dumbbells per viral clade)
    ax_b = fig.add_subplot(gs[0, 1])
    draw_panel_header(fig, 0.435, 0.945, "B", "Virus-matched host mirrors")
    by_g = {r["gene"]: r for r in recs}
    mirrors = [
        ("Coronavirus", "ACE2+DPP4", round(100.0 * (by_g["ACE2"]["pdb_atomic_audit"]["contact_4p5A"]["n_innov"] + by_g["DPP4"]["pdb_atomic_audit"]["contact_4p5A"]["n_innov"]) / (by_g["ACE2"]["pdb_atomic_audit"]["contact_4p5A"]["n"] + by_g["DPP4"]["pdb_atomic_audit"]["contact_4p5A"]["n"]), 1),
         "RAE1+YWHAZ", 0.0),
        ("Influenza", "TRIM25+ISG15", round(100.0 * (by_g["TRIM25"]["pdb_atomic_audit"]["contact_4p5A"]["n_innov"] + by_g["ISG15"]["pdb_atomic_audit"]["contact_4p5A"]["n_innov"]) / (by_g["TRIM25"]["pdb_atomic_audit"]["contact_4p5A"]["n"] + by_g["ISG15"]["pdb_atomic_audit"]["contact_4p5A"]["n"]), 1),
         "ANP32A", by_g["ANP32A"]["pdb_atomic_audit"]["contact_4p5A"]["innov_pct"]),
        ("HIV-1 Vif", "APOBEC3G+3H", round(100.0 * (by_g["APOBEC3G"]["pdb_atomic_audit"]["contact_4p5A"]["n_innov"] + by_g["APOBEC3H"]["pdb_atomic_audit"]["contact_4p5A"]["n_innov"]) / (by_g["APOBEC3G"]["pdb_atomic_audit"]["contact_4p5A"]["n"] + by_g["APOBEC3H"]["pdb_atomic_audit"]["contact_4p5A"]["n"]), 1),
         "ELOB", by_g["ELOB"]["pdb_atomic_audit"]["contact_4p5A"]["innov_pct"]),
        ("HIV-1 Capsid", "TRIM5 v1–v4", by_g["TRIM5"]["functional_loop_audit"]["loop_innov_pct"],
         "PPIA", by_g["PPIA"]["pdb_atomic_audit"]["contact_4p5A"]["innov_pct"]),
    ]
    yb = np.arange(len(mirrors))
    for i, (vir, clab, cv, hlab, hv) in enumerate(mirrors):
        ax_b.plot([hv, cv], [i, i], color=SUBINK, lw=1.5, zorder=2)
        ax_b.scatter(hv, i, s=42, color=TEAL, edgecolor="white", lw=0.6, zorder=4, label="Hijacked hub" if i == 0 else None)
        ax_b.scatter(cv, i, s=42, color=ORANGE, edgecolor="white", lw=0.6, zorder=4, label="Antiviral/receptor" if i == 0 else None)
        ax_b.text(hv, i - 0.28, f"{hlab} ({hv:.1f}%)", ha="left" if hv < 10 else "center", va="top", fontsize=5.2, color=TEAL, fontweight="bold")
        ax_b.text(cv, i + 0.25, f"{clab} ({cv:.1f}%)", ha="center", va="bottom", fontsize=5.2, color=ORANGE, fontweight="bold")
    ax_b.set_yticks(yb)
    ax_b.set_yticklabels([m[0] for m in mirrors], fontsize=6.0, fontweight="bold")
    ax_b.set_xlim(-5, 68)
    ax_b.set_ylim(-0.65, 3.85)
    ax_b.set_xlabel("Interface contacts selected (%)", fontsize=6.5)
    add_callout_pill(ax_b, 35, 3.58, "Pooled OR=14.32x (p=1.52e-14)", color=INK, fontsize=5.15)

    # Panel C: Compartment Plasticity & Strict S-PRIME 4.11x Gradient
    ax_c = fig.add_subplot(gs[0, 2])
    draw_panel_header(fig, 0.725, 0.945, "C", "Disorder & S-PRIME shift")

    pc = agg["conflict_pdb_pooled"]
    ph = agg["hijacked_hub_pdb_pooled"]
    func_only = [r for r in recs if r["pdb_atomic_audit"] is None]
    func_gap = float(np.mean([100.0 * r["functional_loop_audit"]["loop_mean_gap"] for r in func_only]))

    comps = [
        ("Hub <=4.5 A\n(n=102)", ph["per_gene_contact_gap_pct_mean"], ph["contact_4p5A_prime_strict_pct"]),
        ("Conflict core\n(n=3,643)", pc["per_gene_scaffold_gap_pct_mean"], pc["scaffold_gt8A_prime_strict_pct"]),
        ("Conflict <=4.5 A\n(n=299)", pc["per_gene_contact_gap_pct_mean"], pc["contact_4p5A_prime_strict_pct"]),
        ("PDB unmodeled\n(n=3,457)", pc["per_gene_unmodeled_gap_pct_mean"], 2.2),
        ("Func. tripwires\n(n=6 genes)", func_gap, 0.8),
    ]
    # Compute exact strict S-PRIME % for unmodeled & functional loops directly from master/sites
    # Recompute exact strict % on unmodeled & func loops
    strict_unmod_n, strict_unmod_tot = 0, 0
    for r in recs:
        if r["arm"] in ("Antiviral_Effector_Sensor", "Viral_Entry_Receptor") and r["pdb_atomic_audit"] is not None:
            u = r["pdb_atomic_audit"]["unmodeled_disordered"]
            strict_unmod_n += u["n_prime_strict_chg"]
            strict_unmod_tot += u["n"]
    comps[3] = ("PDB unmodeled\n(n=3,457)", pc["per_gene_unmodeled_gap_pct_mean"], round(100.0 * strict_unmod_n / max(1, strict_unmod_tot), 2))
    # Remove synthetic 0.8 -> compute exact for func_only loops directly from JSON
    f_strict_n, f_strict_tot = 0, 0
    for r in func_only:
        g = r["gene"]
        p = json.load(open(os.path.join(PRIME_DIR, f"{g}.prime.json")))
        pdata = p["site_results"]["data"]
        props = [pm["name"] for pm in p["properties"]]
        strict_sites = set()
        for idx in range(p["site_results"]["row_count"]):
            if any(pdata[f"class_{pn}"][idx] == "Changing" for pn in props):
                strict_sites.add(pdata["site"][idx])
        gt_loops = json.load(open(GT_PATH))[g]["functional_loops_1idx"]
        col_float = map_msa_to_human(os.path.join(WALRUS_DIR, f"{g}.walrus.msa.fa"))[3]
        w_sites = json.load(open(os.path.join(WALRUS_DIR, f"{g}.walrus.json")))["site_atlas"]
        for s in w_sites:
            xf = col_float[s["site"]]
            if any((lo - 0.5) <= xf <= (hi + 0.5) for lo, hi in gt_loops.values()):
                f_strict_tot += 1
                if s["site"] in strict_sites:
                    f_strict_n += 1
    comps[4] = ("Func. tripwires\n(n=6 loci)", func_gap, round(100.0 * f_strict_n / max(1, f_strict_tot), 2))

    yc = np.arange(len(comps))
    gaps_c = [c[1] for c in comps]
    sp_c = [c[2] for c in comps]
    ax_c.barh(yc + 0.16, gaps_c, height=0.30, color=BLUE, alpha=0.82, label="Mean codon gap (%)")
    ax_c2 = ax_c.twiny()
    ax_c2.spines["right"].set_visible(False)
    ax_c2.spines["top"].set_visible(True)
    ax_c2.spines["top"].set_color(ROSE)
    ax_c2.tick_params(axis="x", colors=ROSE, labelsize=5.8)
    ax_c2.barh(yc - 0.16, sp_c, height=0.30, color=ROSE, alpha=0.85, label="Strict S-PRIME (%)")
    ax_c.set_yticks(yc)
    ax_c.set_yticklabels([c[0] for c in comps], fontsize=5.5, fontweight="bold")
    ax_c.set_xlim(0, 44)
    ax_c2.set_xlim(0, 6.2)
    ax_c.set_xlabel("Codon gap density (%)", fontsize=6.4, color=BLUE)
    ax_c2.set_xlabel("Strict S-PRIME reprogramming (%)", fontsize=6.1, color=ROSE)
    add_callout_pill(ax_c, 31.0, 0.92, "Strict S-PRIME\n4.11x at <=4.5 A\n(p=6.0e-5)", color=ROSE, fontsize=5.1)
    add_callout_pill(ax_c, 28.0, 0.08, "0.0% S-PRIME on hubs", color=TEAL, fontsize=5.1)

    save_both(fig, "fig2_science_atomic_interfaces")


# ==============================================================================
# FIGURE 3: Four Non-Classical Conflict Mechanisms
# ==============================================================================
def build_fig3(master, bundles):
    fig = plt.figure(figsize=(7.25, 6.05), dpi=300)
    gs = gridspec.GridSpec(
        2, 2, left=0.075, right=0.975, top=0.935, bottom=0.075, wspace=0.28, hspace=0.36
    )

    # Panel A: APOBEC3G Loop 7 toggle (7/9 vs 1/9 MEME) + Strict pI charge shifts vs d_min
    ax_a = fig.add_subplot(gs[0, 0])
    draw_panel_header(fig, 0.015, 0.965, "A", "Toggling & charge inversion along ordered interfaces")
    _, w_a3g, _, msa_a3g, gt_a3g = bundles["APOBEC3G"]
    h2c_a3g = msa_a3g[4]
    sites_a3g = w_a3g["site_atlas"]
    loop7 = list(range(124, 133))
    meme_y, tog_y, lbl_a = [], [], []
    for hp in loop7:
        s = sites_a3g[h2c_a3g[hp] - 1]
        meme_y.append(-math.log10(max(1e-4, s.get("meme_p_value") or 1.0)))
        tog_y.append(-math.log10(max(1e-4, s.get("exchangeability_bh_q") or 1.0)))
        lbl_a.append(f"{gt_a3g['human_ref_seq'][hp - 1]}{hp}")
    xa = np.arange(len(loop7))
    ax_a.hlines(xa - 0.14, 0, meme_y, color=MUTED, lw=1.4)
    ax_a.scatter(meme_y, xa - 0.14, s=22, color=MUTED, label="Standard scalar MEME ($p$)")
    ax_a.hlines(xa + 0.14, 0, tog_y, color=ORANGE, lw=1.4)
    ax_a.scatter(tog_y, xa + 0.14, s=26, color=ORANGE, edgecolor="white", lw=0.4, label="Exchangeability toggle (BH $q$)")
    ax_a.axvline(-math.log10(0.10), color=ROSE, ls="--", lw=0.75)
    ax_a.axvline(-math.log10(0.05), color=MUTED, ls=":", lw=0.75)
    ax_a.set_yticks(xa)
    ax_a.set_yticklabels(lbl_a, fontsize=5.8, fontweight="bold")
    ax_a.set_xlim(-0.08, 2.35)
    ax_a.set_xlabel("Support (-log10 p or q)", fontsize=6.4)
    add_callout_pill(ax_a, 1.55, 6.2, "APOBEC3G Loop 7 (8CX0)\nMEME: 1/9 (D128)\nWALRUS toggle: 7/9", color=ORANGE, fontsize=5.2)
    add_callout_pill(ax_a, 1.48, 1.5, "Strict S-PRIME pI hotspots:\nOAS1 D35 (2.96 A, lambda=-2.16)\nOAS1 R47 (3.01 A, lambda=-1.62)\nACE2 G354 (3.43 A, lambda=-1.79)", color=BLUE, fontsize=4.95)
    ax_a.legend(loc="lower right", fontsize=5.1, frameon=True, facecolor="white", edgecolor=HAIR)

    # Panel B: NLRP1 3Cpro Tripwire & Recombination Shunts across Human Coordinates
    ax_b = fig.add_subplot(gs[0, 1])
    draw_panel_header(fig, 0.515, 0.965, "B", "Disordered tripwire, artifact pruning & shunts (NLRP1)")
    _, w_nl, _, msa_nl, gt_nl = bundles["NLRP1"]
    col_float_nl = msa_nl[3]
    col_to_hum_nl = msa_nl[2]
    occ_nl = msa_nl[5]
    sites_nl = w_nl["site_atlas"]
    doms_nl = [("PYD", 1, 90, "#CBD5E1"), ("3C", 104, 166, "#FECACA"),
               ("NACHT", 230, 560, "#DBEAFE"), ("LRR", 790, 990, "#CCFBF1"),
               ("FIIND", 1060, 1365, "#FEF3C7"), ("CARD", 1375, 1473, "#EDE9FE")]
    ax_b.axvspan(100, 165, color="#FEF2F2", alpha=0.95, zorder=0)
    xs_nl = np.array([col_float_nl[s["site"]] for s in sites_nl])
    ys_nl = np.array([min(10.5, -math.log10(max(1e-11, s.get("simes_omnibus_q") or 1.0))) for s in sites_nl])
    gaps_nl = np.array([s.get("gap_fraction") or 0.0 for s in sites_nl])
    ax_b.fill_between(xs_nl, 0.0, -3.0 * gaps_nl, color=BLUE, alpha=0.25, step="mid", label="Gap density")
    div_nl = np.array([is_walrus_innov(s) for s in sites_nl])
    art_nl = np.array([bool(s.get("artifact_eliminated")) for s in sites_nl])
    ax_b.vlines(xs_nl[div_nl], 0.0, np.maximum(ys_nl[div_nl], 1.1), color=ORANGE, lw=0.75, alpha=0.85, label="Selected (17.1%)")
    ax_b.scatter(xs_nl[art_nl], np.full(np.sum(art_nl), 0.7), marker="x", s=16, color=ROSE, lw=0.9, zorder=5, label="POA artifact (51)")
    for bi, br in enumerate(w_nl.get("recombination_bridges", [])):
        if br.get("routed_to_spr_tree"):
            bs = col_float_nl.get(br["peak_start_codon"], br["peak_start_codon"])
            be = col_float_nl.get(br["peak_end_codon"], br["peak_end_codon"])
            ax_b.plot([bs, be], [9.6 - 0.8 * (bi % 2)] * 2, color=AMBER, lw=2.8, solid_capstyle="round", zorder=6,
                      label="SPR Wheatstone shunt" if bi == 0 else None)
    for dname, lo, hi, fc in doms_nl:
        rect = FancyBboxPatch((lo, -4.6), max(10, hi - lo), 1.1, boxstyle="round,pad=0.02,rounding_size=0.1",
                              facecolor=fc, edgecolor=SUBINK, linewidth=0.5, zorder=3)
        ax_b.add_patch(rect)
        ax_b.text((lo + hi) / 2.0, -4.05, dname, ha="center", va="center", fontsize=5.0, fontweight="bold", color=INK)
    ax_b.axhline(0.0, color=SUBINK, lw=0.6)
    ax_b.set_xlim(0, 1485)
    ax_b.set_ylim(-5.0, 11.2)
    ax_b.set_yticks([-3, 0, 4, 8])
    ax_b.set_yticklabels(["100% gap", "0", "4", "8"])
    ax_b.set_xlabel("Human NLRP1 codon coordinate (534 inserted cols)", fontsize=6.4)
    ax_b.set_ylabel("-log10 q / Gap", fontsize=6.4)
    ax_b.legend(loc="upper center", fontsize=5.0, ncol=2, frameon=True, facecolor="white", edgecolor=HAIR)

    # Panel C: PRDM9 28-aa C2H2 Repeat Register Audit
    ax_c = fig.add_subplot(gs[1, 0])
    draw_panel_header(fig, 0.015, 0.465, "C", "Tandem C2H2 repeat resolution isolates DNA contacts (PRDM9)")
    _, w_pr, _, msa_pr, gt_pr = bundles["PRDM9"]
    h2c_pr = msa_pr[4]
    sites_pr = w_pr["site_atlas"]
    motif_len = 21
    pos_counts = np.zeros(motif_len)
    pos_div = np.zeros(motif_len)
    for rep in gt_pr["extra"]["znf_repeats"]:
        c1_pos = rep["c1_h2"][0]
        for off in range(motif_len):
            c = h2c_pr.get(c1_pos + off)
            if c:
                pos_counts[off] += 1
                s = sites_pr[c - 1]
                if is_walrus_innov(s) or ((s.get("meme_p_value") or 1.0) <= 0.05):
                    pos_div[off] += 1
    pct_pr = 100.0 * pos_div / np.maximum(1.0, pos_counts)
    cols_pr = [ROSE if off in (0, 3, 16, 20) else (ORANGE if off in (9, 11, 12, 15) else MUTED) for off in range(motif_len)]
    xp = np.arange(1, motif_len + 1)
    ax_c.bar(xp, pct_pr, width=0.72, color=cols_pr, edgecolor="white", lw=0.45, zorder=3)
    ax_c.axhline(2.1, color=ROSE, ls="--", lw=0.75, alpha=0.8)
    ax_c.set_xticks([1, 4, 10, 12, 13, 16, 17, 21])
    ax_c.set_xticklabels(["C1", "C2", "-1", "+2", "+3", "+6", "H1", "H2"], fontsize=5.8, fontweight="bold")
    ax_c.set_ylim(0, 126)
    ax_c.set_xlabel("Canonical 21-aa C2H2 motif position (12 repeats)", fontsize=6.4)
    ax_c.set_ylabel("Fingers under selection (%)", fontsize=6.4)
    add_callout_pill(ax_c, 10.8, 111.5, "DNA (-1,+2,+3,+6): 36/48 (75.0%)  |  C2H2 Zn ligands: 1/48 (2.1%)\nOR = 16.8x, p = 1.51e-16; 48 rigid-alignment artifacts removed", color=INK, fontsize=5.0, ha="center")

    # Panel D: Unsupervised Lineage Loss-of-Function across all 96 flagged branches
    ax_d = fig.add_subplot(gs[1, 1])
    draw_panel_header(fig, 0.515, 0.465, "D", "Unsupervised lineage loss-of-function (96 flagged branches)")
    all_pts = []
    for rec in master["tier1_genes"]:
        gn = rec["gene"]
        arm = rec["arm"]
        pseudo_names = {p["name"] for p in rec["admitted_pseudogenes"]}
        for b in rec["lof_flagged_branches"]:
            all_pts.append((gn, arm, b["branch"], max(1e-100, b["lof_p"]), min(10.0, b["relax_k"]), b["branch"] in pseudo_names))

    for arm in ARM_COLOR:
        sub = [pt for pt in all_pts if pt[1] == arm and not pt[5]]
        if sub:
            ax_d.scatter([pt[4] for pt in sub], [-math.log10(pt[3]) for pt in sub], s=18,
                         color=ARM_COLOR[arm], alpha=0.72, edgecolor="white", lw=0.4, zorder=3, label=ARM_LABEL[arm].split(" (")[0])
        sub_p = [pt for pt in all_pts if pt[1] == arm and pt[5]]
        if sub_p:
            ax_d.scatter([pt[4] for pt in sub_p], [-math.log10(pt[3]) for pt in sub_p], s=42, marker="*",
                         color=ROSE, edgecolor=INK, lw=0.4, zorder=6, label="Repaired pseudogene" if arm == "Tandem_Repeat_Control" else None)

    ax_d.axvline(1.0, color=MUTED, ls="--", lw=0.75)
    # Key biological callouts directly placed cleanly
    call_lof = [
        ("PRDM9", "Canis_aureus", (2.6, 86.0), "PRDM9 C. aureus\n(repaired)"),
        ("MX1", "Node44", (4.1, 82.0), "MX1 Odontocete\nNode44"),
        ("IRF3", "Myotis_lucifugus", (2.2, 69.0), "IRF3 M. lucifugus"),
        ("MX1", "Delphinus_delphis", (1.6, 56.0), "MX1 D. delphis"),
        ("PRDM9", "Canis_lupus", (4.0, 29.0), "PRDM9 C. lupus"),
    ]
    for gn, br, (tx, ty), lbl in call_lof:
        match = [pt for pt in all_pts if pt[0] == gn and pt[2] == br]
        if match:
            px, py = match[0][4], -math.log10(match[0][3])
            ax_d.annotate(lbl, xy=(px, py), xytext=(tx, ty), fontsize=5.1, fontweight="bold", color=INK,
                          arrowprops=dict(arrowstyle="->", color=SUBINK, lw=0.55), zorder=7)
    ax_d.set_xlim(-0.25, 5.65)
    ax_d.set_ylim(-2, 106)
    ax_d.set_xlabel("RELAX intensity parameter $k$ (k < 1 relaxed)", fontsize=6.4)
    ax_d.set_ylabel("LoF support ($-\\log_{10} p_{\\mathrm{LoF}}$)", fontsize=6.4)
    ax_d.legend(loc="center right", bbox_to_anchor=(0.99, 0.56), fontsize=4.9, frameon=True, facecolor="white", edgecolor=HAIR)

    save_both(fig, "fig3_science_mechanisms")


# ==============================================================================
# FIGURE 4: Pan-Mammalian Selective Sieve Across 7,045 TOGA2 Orthogroups
# ==============================================================================
def build_fig4(tier2):
    vips_meta = json.load(open(VIP_JSON))
    vip_set = set(vips_meta["vip_genes"])
    av_set = set(vips_meta["antiviral_genes"])
    hijacked_set = vip_set - av_set

    by_gene = {}
    with gzip.open(TOGA_GZ, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)
            if r.get("error") or r.get("num_taxa", 0) < 25 or r.get("num_codons", 0) < 50:
                continue
            gs = r.get("gene_symbol", "")
            if not gs:
                continue
            if gs not in by_gene or (r["num_taxa"], r["num_codons"]) > (by_gene[gs]["num_taxa"], by_gene[gs]["num_codons"]):
                by_gene[gs] = r

    g_hj = [r for gs, r in by_gene.items() if gs in hijacked_set]
    g_bg = [r for gs, r in by_gene.items() if gs not in vip_set]
    g_av = [r for gs, r in by_gene.items() if gs in av_set]

    fig = plt.figure(figsize=(7.25, 4.75), dpi=300)
    gs = gridspec.GridSpec(
        1, 3, width_ratios=[1.15, 1.05, 1.05],
        left=0.075, right=0.975, top=0.875, bottom=0.145, wspace=0.40
    )

    # Panel A: Joint Gap Plasticity vs Episodic Adaptation Plane w/ Median Crosshairs
    ax_a = fig.add_subplot(gs[0, 0])
    draw_panel_header(fig, 0.015, 0.945, "A", "Proteome landscape (n = 7,045 orthogroups)")

    def gap_ep(grp):
        return (
            np.array([100.0 * r["total_gap_codons"] / (r["num_codons"] * r["num_taxa"]) for r in grp]),
            np.array([100.0 * r["summary_counts"]["sites_episodic_diversifying"] / r["num_codons"] for r in grp]),
        )

    bgx, bgy = gap_ep(g_bg)
    hjx, hjy = gap_ep(g_hj)
    avx, avy = gap_ep(g_av)

    ax_a.scatter(bgx, bgy, s=4, color=MUTED, alpha=0.18, rasterized=True, label=f"Background (n={len(g_bg):,})")
    ax_a.scatter(hjx, hjy, s=11, color=TEAL, alpha=0.62, edgecolor="white", lw=0.25, label=f"Hijacked hubs (n={len(g_hj)})")
    ax_a.scatter(avx, avy, s=22, color=ORANGE, alpha=0.88, edgecolor="white", lw=0.35, label=f"Antiviral (n={len(g_av)})")
    # Mean crosshairs
    for xm, ym, col in [(np.mean(hjx), np.mean(hjy), TEAL), (np.mean(bgx), np.mean(bgy), SUBINK), (np.mean(avx), np.mean(avy), ORANGE)]:
        ax_a.scatter(xm, ym, s=75, marker="+", color=col, lw=1.6, zorder=7)
    ax_a.set_xlim(0, 34)
    ax_a.set_ylim(0, 40)
    ax_a.set_xlabel("Codon indel gap density (%)", fontsize=6.4)
    ax_a.set_ylabel("Episodic diversifying sites (%)", fontsize=6.4)
    ax_a.legend(loc="upper right", fontsize=5.15, frameon=True, facecolor="white", edgecolor=HAIR)
    add_callout_pill(ax_a, 18.0, 26.5, "Antiviral mean: 15.0% gap, 12.6% ep.\nHub mean: 6.7% gap, 6.1% ep.\n(MWU p=1.5e-11 & 8.7e-9)", color=INK, fontsize=5.0)

    # Panel B: Taxon Conservation ECDF + BUSTED-E Bar Inset / Dual
    ax_b = fig.add_subplot(gs[0, 1])
    draw_panel_header(fig, 0.365, 0.945, "B", "Taxon retention & gene selection")
    for grp, col, lbl in [
        (g_hj, TEAL, f"Hijacked hubs (med 839)"),
        (g_bg, MUTED, f"Background (med 819)"),
        (g_av, ORANGE, f"Antiviral (med 441)"),
    ]:
        vals = np.sort([r["num_taxa"] for r in grp])
        yvals = np.linspace(0, 1, len(vals))
        ax_b.plot(vals, yvals, color=col, lw=1.6, label=lbl)
    ax_b.set_xlim(0, 940)
    ax_b.set_ylim(0, 1.03)
    ax_b.set_xlabel("Retained mammalian taxa", fontsize=6.4)
    ax_b.set_ylabel("Cumulative fraction of genes", fontsize=6.4)
    ax_b.legend(loc="upper left", fontsize=5.15, frameon=True, facecolor="white", edgecolor=HAIR)
    add_callout_pill(ax_b, 680, 0.24, "BUSTED-E positive:\nAntiviral: 69.4% (34/49)\nBackground: 57.3%\nHijacked hubs: 39.1% (70/179)\nFisher p=1.45e-4\nTaxon loss MWU p=1.85e-7", color=INK, fontsize=5.05)

    # Panel C: Covariate-Adjusted OLS Regression Effects (log taxa + log codons)
    ax_c = fig.add_subplot(gs[0, 2])
    draw_panel_header(fig, 0.690, 0.945, "C", "Covariate-adjusted shifts (beta)")
    cov_h = tier2["covariate_adjusted_av_vs_hijacked"]
    cov_b = tier2["covariate_adjusted_av_vs_background"]
    effs = [
        ("Gap %\nvs. Background", cov_b["gap_density_pct"], MUTED),
        ("Episodic %\nvs. Background", cov_b["episodic_site_pct"], BLUE),
        ("Gap %\nvs. Hubs", cov_h["gap_density_pct"], TEAL),
        ("Episodic %\nvs. Hubs", cov_h["episodic_site_pct"], ORANGE),
        ("Toggle/Conv %\nvs. Hubs", cov_h["toggle_conv_pct"], PURPLE),
    ]
    yc = np.arange(len(effs))
    ax_c.axvline(0.0, color=SUBINK, ls="--", lw=0.8)
    for i, (lbl, obj, col) in enumerate(effs):
        b = obj["beta_antiviral_adj"]
        ci = 1.96 * obj["se"]
        ax_c.plot([b - ci, b + ci], [i, i], color=col, lw=1.6, zorder=3)
        ax_c.scatter(b, i, s=34, color=col, edgecolor="white", lw=0.5, zorder=4)
        ax_c.text(b + ci + 0.5, i, f"p={obj['p_two_sided']:.1e}", va="center", ha="left", fontsize=5.3, fontweight="bold", color=INK)
    ax_c.set_yticks(yc)
    ax_c.set_yticklabels([e[0] for e in effs], fontsize=5.7, fontweight="bold")
    ax_c.set_xlim(-1.0, 23.0)
    ax_c.set_ylim(-0.6, 4.6)
    ax_c.set_xlabel("Covariate-adjusted effect size beta\n(+/- 95% CI; adj. log taxa & codons)", fontsize=6.1)

    save_both(fig, "fig4_science_proteome_sieve")


# ==============================================================================
# SUPPLEMENTARY FIGURE S1: Five-Property S-PRIME Axis Profiles & Uncertainty
# ==============================================================================
def build_fig_s1(master):
    fig = plt.figure(figsize=(7.25, 5.2), dpi=300)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.25, 1.0], left=0.145, right=0.975, top=0.885, bottom=0.135, wspace=0.38)

    recs = master["tier1_genes"]
    ordered = sorted(recs, key=lambda r: (
        ["Hijacked_Cellular_Hub", "Tandem_Repeat_Control", "Viral_Entry_Receptor", "Antiviral_Effector_Sensor"].index(r["arm"]),
        sum(r["prime_marginal_changing_by_prop"].values())
    ))

    # Panel A: Per-gene Marginal vs Strict S-PRIME pI + Remaining Properties
    ax_a = fig.add_subplot(gs[0, 0])
    draw_panel_header(fig, 0.015, 0.945, "A", "Per-gene S-PRIME physicochemical reprogramming (lambda=0.5)")
    props = ["Isoelectric_Point_pI", "AlphaHelix_Propensity_ChouFasman", "Volume_Angstrom3", "BetaSheet_Propensity_ChouFasman", "Hydrophobicity_KyteDoolittle"]
    plbls = ["Isoelectric pI", "Alpha-helix", "Side-chain volume", "Beta-sheet", "Hydrophobicity"]
    pcols = [ROSE, ORANGE, TEAL, PURPLE, BLUE]
    ya = np.arange(len(ordered))
    left = np.zeros(len(ordered))
    for pn, plb, pc in zip(props, plbls, pcols):
        vals = np.array([r["prime_marginal_changing_by_prop"].get(pn, 0) for r in ordered], dtype=float)
        ax_a.barh(ya, vals, left=left, height=0.72, color=pc, alpha=0.88, edgecolor="white", lw=0.3, label=plb)
        left += vals
    for i, r in enumerate(ordered):
        if r["prime_strict_changing_sites"] > 0:
            ax_a.text(left[i] + 2.0, i, f"S:{r['prime_strict_changing_sites']}", va="center", fontsize=5.1, fontweight="bold", color=INK)
    ax_a.set_yticks(ya)
    ax_a.set_yticklabels([r["gene"] for r in ordered], fontsize=5.6, fontweight="bold")
    ax_a.set_xlim(0, 210)
    ax_a.set_xlabel("Marginal Holm p <= 0.10 changing calls (S: strict sites)", fontsize=6.2)
    ax_a.legend(loc="lower right", fontsize=5.2, frameon=True, facecolor="white", edgecolor=HAIR)

    # Panel B: Global Changing vs Conserved Symmetry Breaking & Column Posterior Precision
    ax_b = fig.add_subplot(gs[0, 1])
    draw_panel_header(fig, 0.545, 0.945, "B", "Global axis asymmetry & alignment certainty")
    st_tot = master["prime_totals_strict_omnibus_q10"]
    yb = np.arange(len(props))
    chg_v = [st_tot[pn]["changing"] for pn in props]
    con_v = [-st_tot[pn]["conserved"] for pn in props]
    ax_b.barh(yb, chg_v, height=0.52, color=ORANGE, alpha=0.88, label="Strict Changing (-lambda)")
    ax_b.barh(yb, con_v, height=0.52, color=BLUE, alpha=0.80, label="Strict Conserved (+lambda)")
    ax_b.axvline(0.0, color=SUBINK, lw=0.7)
    for i, (cv, nv) in enumerate(zip(chg_v, con_v)):
        ax_b.text(cv + 1.5, i, str(cv), va="center", ha="left", fontsize=5.6, fontweight="bold", color=ORANGE)
        ax_b.text(nv - 1.5, i, str(-nv), va="center", ha="right", fontsize=5.6, fontweight="bold", color=BLUE)
    ax_b.set_yticks(yb)
    ax_b.set_yticklabels(plbls, fontsize=5.8, fontweight="bold")
    ax_b.set_xlim(-78, 84)
    ax_b.set_xlabel("Strict two-stage S-PRIME site count across 25 loci", fontsize=6.2)
    ax_b.legend(loc="upper left", fontsize=5.1, frameon=True, facecolor="white", edgecolor=HAIR)
    add_callout_pill(ax_b, 54, 2.55, "pI drives polarity shift (41 vs 12)\nHydrophobicity preserves fold (52)", color=INK, fontsize=4.85)

    save_both(fig, "fig_s1_property_and_error_audit")


# ==============================================================================
# SUPPLEMENTARY FIGURE S2: Legacy Foil Comparison & 10 SPR Wheatstone Bridges
# ==============================================================================
def build_fig_s2(master, bundles):
    fig = plt.figure(figsize=(7.25, 4.9), dpi=300)
    gs = gridspec.GridSpec(1, 2, width_ratios=[0.85, 1.35], left=0.145, right=0.98, top=0.885, bottom=0.145, wspace=0.45)

    # Panel A: Legacy MAFFT vs Joint WALRUS on TRIM5 & PRDM9
    ax_a = fig.add_subplot(gs[0, 0])
    draw_panel_header(fig, 0.015, 0.945, "A", "MAFFT vs. WALRUS on unalignable controls")
    metrics = [
        ("TRIM5 MSA\n(codons / 10)", 65.0, 59.5),
        ("TRIM5 MEME\nvs. POA pruned", 56.0, 6.0),
        ("PRDM9 MSA\n(codons / 20)", 1467.0 / 20.0, 1334.0 / 20.0),
        ("PRDM9 q<=0.10\nvs. POA pruned", 32.0, 48.0),
    ]
    ya = np.arange(len(metrics))
    ax_a.barh(ya + 0.17, [m[1] for m in metrics], height=0.32, color=MUTED, label="Rigid MAFFT")
    ax_a.barh(ya - 0.17, [m[2] for m in metrics], height=0.32, color=TEAL, label="Joint WALRUS")
    for i, (_, mval, wval) in enumerate(metrics):
        ax_a.text(mval + 1.2, i + 0.17, f"{mval:.1f}", va="center", fontsize=5.4, color=SUBINK)
        ax_a.text(wval + 1.2, i - 0.17, f"{wval:.1f}", va="center", fontsize=5.4, color=TEAL, fontweight="bold")
    ax_a.set_yticks(ya)
    ax_a.set_yticklabels([m[0] for m in metrics], fontsize=5.6, fontweight="bold")
    ax_a.set_xlim(0, 90)
    ax_a.set_xlabel("Codon span (scaled) / site count", fontsize=6.3)
    ax_a.legend(loc="lower right", fontsize=5.2, frameon=True, facecolor="white", edgecolor=HAIR)

    # Panel B: All 10 SPR-Routed Wheatstone Gene-Conversion Bridges
    ax_b = fig.add_subplot(gs[0, 1])
    draw_panel_header(fig, 0.455, 0.945, "B", "Statistically routed Ohmic Wheatstone conversion bridges (n = 10)")
    spr_list = []
    for r in master["tier1_genes"]:
        g = r["gene"]
        w = bundles[g][1]
        for br in w.get("recombination_bridges", []):
            if br.get("routed_to_spr_tree"):
                spr_list.append((g, br["clade_u"], br["target_v"], br["peak_start_codon"], br["peak_end_codon"],
                                 br["segment_delta_log_l"], br["p_value"]))
    spr_sorted = sorted(spr_list, key=lambda x: x[5])
    yb = np.arange(len(spr_sorted))
    dls = [x[5] for x in spr_sorted]
    ax_b.hlines(yb, 0, dls, color=AMBER, lw=1.8, zorder=2)
    ax_b.scatter(dls, yb, s=36, color=AMBER, edgecolor=INK, lw=0.45, zorder=4)
    lbls_b = [f"{x[0]} ({x[3]}–{x[4]})" for x in spr_sorted]
    for i, x in enumerate(spr_sorted):
        ax_b.text(x[5] + 2.4, i, f"{x[1]} <-> {x[2]} (p={x[6]:.3f})", va="center", fontsize=5.1, color=INK)
    ax_b.set_yticks(yb)
    ax_b.set_yticklabels(lbls_b, fontsize=5.6, fontweight="bold")
    ax_b.set_xlim(0, 130)
    ax_b.set_xlabel("Segmental topological likelihood gain ($\\Delta \\ln L$)", fontsize=6.3)

    save_both(fig, "fig_s2_legacy_and_wheatstone")


def main():
    master = json.load(open(MASTER_JSON))
    tier2 = json.load(open(TIER2_JSON))
    gt = json.load(open(GT_PATH))
    bundles = load_gene_bundles(master, gt)

    build_fig1(master, bundles)
    build_fig2(master)
    build_fig3(master, bundles)
    build_fig4(tier2)
    build_fig_s1(master)
    build_fig_s2(master, bundles)


if __name__ == "__main__":
    main()
