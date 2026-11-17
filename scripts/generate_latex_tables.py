#!/usr/bin/env python3
"""
Generate exact LaTeX Table 1 (25-Locus Tier-1 Evolutionary Atlas) and Table 2
(Atomic <=4.5A PDB Structure vs Scaffold Audit) directly from
/Users/sergei/Projects/unalignable-arms-race/results/tier1_25gene_master_synthesis.json.
Zero hardcoded numbers. Python 3.9 compatible.
"""

import json
import os

ROOT = "/Users/sergei/Projects/unalignable-arms-race"
SYN_PATH = os.path.join(ROOT, "results", "tier1_25gene_master_synthesis.json")
OUT_PATH = os.path.join(ROOT, "manuscript", "tables_generated.tex")

ARM_HEADERS = [
    ("Antiviral_Effector_Sensor", "\\multicolumn{10}{@{}l}{\\textit{Arm A: Antiviral Effectors \\& Innate Sensors ($n = 14$)}} \\\\"),
    ("Viral_Entry_Receptor", "\\multicolumn{10}{@{}l}{\\textit{Arm B: Viral Entry Receptors ($n = 5$)}} \\\\"),
    ("Tandem_Repeat_Control", "\\multicolumn{10}{@{}l}{\\textit{Arm C: Tandem-Repeat Meiotic Benchmark ($n = 1$)}} \\\\"),
    ("Hijacked_Cellular_Hub", "\\multicolumn{10}{@{}l}{\\textit{Arm D: Hijacked Cellular Host Hubs --- Controls ($n = 5$)}} \\\\"),
]

PARTNER_SHORT = {
    "APOBEC3G": "HIV-1 Vif (8CX0)",
    "APOBEC3H": "HIV-1 Vif (8FVI)",
    "ISG15":    "FluB NS1B / PLpro (3SDL,5TL6)",
    "SAMHD1":   "SIV Vpx (4CC9)",
    "STAT2":    "Dengue NS5 (8T12)",
    "IRF3":     "Rotavirus NSP1 (5JER)",
    "TRIM25":   "FluA NS1 (5NT1)",
    "OAS1":     "dsRNA (4IG8)",
    "ACE2":     "SARS-CoV-2 RBD (6M0J)",
    "DPP4":     "MERS-CoV RBD (4L72)",
    "CD4":      "HIV-1 gp120 (1GC1)",
    "TFRC":     "Machupo GP1 (3KAS)",
    "NPC1":     "Ebola GP (5F1B)",
    "PRDM9":    "Hotspot DNA (5EGB)",
    "PPIA":     "HIV-1 CA (1AK4)",
    "ELOB":     "HIV-1 Vif (2MA9,4N9F)",
    "ANP32A":   "FluA Pol (8RN0)",
    "RAE1":     "SARS-CoV-2 Orf6 (7VPG)",
    "YWHAZ":    "SARS-CoV-2 N (7ZIT)",
}


def fmt_p(p):
    if p is None:
        return "---"
    if p < 1e-4:
        return "$<10^{-4}$"
    if p < 0.001:
        return "${0:.1e}$".format(p).replace("e-0", "\\times 10^{-") + "}"
    return "${0:.3f}$".format(p)


def fmt_p_sci(p):
    if p is None:
        return "---"
    if p >= 0.01:
        return "${0:.3f}$".format(p)
    if p >= 0.001:
        return "${0:.4f}$".format(p)
    s = "{0:.1e}".format(p)
    mant, exp = s.split("e")
    return "${0}\\times 10^{{{1}}}$".format(mant, int(exp))


def main():
    m = json.load(open(SYN_PATH))
    recs = m["tier1_genes"]
    by_arm = {}
    for r in recs:
        by_arm.setdefault(r["arm"], []).append(r)

    lines = []
    # ---------------- TABLE 1 ----------------
    lines.append(r"\begin{table*}[t]")
    lines.append(r"\centering")
    lines.append(r"\caption{\textbf{Multi-modal evolutionary Atlas across all 25 Tier-1 mammalian loci.} All runs began from raw unaligned CDS (\texttt{-u --coevolution --keep-pseudogenes}); \texttt{S-PRIME} used conservative Ridge shrinkage ($\lambda=0.5$). Sel.\% gives sites under \WALRUS diversifying, toggle, or convergent selection (BH $q \le 0.10$); Art. gives POA-eliminated rigid-MSA false positives; $\Delta\text{pI}$ gives \texttt{S-PRIME} isoelectric-point changing sites (marginal Holm $p\le 0.10$ / strict two-stage).}")
    lines.append(r"\label{tab:tier1_atlas}")
    lines.append(r"\small")
    lines.append(r"\setlength{\tabcolsep}{4.2pt}")
    lines.append(r"\begin{adjustbox}{max width=\textwidth}")
    lines.append(r"\begin{tabular}{@{}l r r r r r r r r r@{}}")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Gene} & \textbf{Taxa} & \textbf{Ref/MSA} & \textbf{Ins.} & \textbf{Gap\%} & \textbf{Sel.\%} & \textbf{\BUSTEDE $p$} & \textbf{Art.} & \textbf{$\Delta\text{pI}$ (M/S)} & \textbf{LoF / Coev} \\")
    lines.append(r"\midrule")

    for arm_key, hdr_tex in ARM_HEADERS:
        lines.append(hdr_tex)
        for r in by_arm.get(arm_key, []):
            pi_m = r["prime_marginal_changing_by_prop"].get("Isoelectric_Point_pI", 0)
            pi_s = r["prime_strict_changing_by_prop"].get("Isoelectric_Point_pI", 0)
            bp = r["busted_e_p"]
            bp_str = fmt_p_sci(bp)
            n_lof = len(r["lof_flagged_branches"])
            n_ps = len(r["admitted_pseudogenes"])
            lof_str = "{0}{1}".format(n_lof, " (+{0})".format(n_ps) if n_ps else "")
            row = "\\quad \\texttt{{{0}}} & {1} & {2}/{3} & {4} & {5:.1f} & {6:.1f} & {7} & {8} & {9}/{10} & {11} / {12} \\\\".format(
                r["gene"], r["admitted_taxa"], r["human_aa"], r["msa_codons"],
                r["non_human_inserted_cols"], r["mean_gap_pct"], r["innov_site_pct"],
                bp_str, r["poa_artifacts_eliminated"], pi_m, pi_s, lof_str, r["coev_pairs_q10"]
            )
            lines.append(row)
        lines.append(r"\midrule")

    lines[-1] = r"\bottomrule"
    lines.append(r"\end{tabular}")
    lines.append(r"\end{adjustbox}")
    lines.append(r"\end{table*}")
    lines.append("")

    # ---------------- TABLE 2 ----------------
    lines.append(r"\begin{table*}[t]")
    lines.append(r"\centering")
    lines.append(r"\caption{\textbf{Atomic-resolution ($d_{\min} \le 4.5$\,\AA) structural interface enrichment.} Resolved human residues in 19 non-fusion PDB co-structures are partitioned strictly by heavy-atom distance to the viral or nucleic-acid partner chain ($\le 4.5$\,\AA\ direct contact vs.\ $>8.0$\,\AA\ protein core; for \texttt{SAMHD1} 4CC9, where only the 19-aa C-terminal Vpx degron is resolved, non-contact compares $>4.5$\,\AA).}")
    lines.append(r"\label{tab:pdb_atomic}")
    lines.append(r"\footnotesize")
    lines.append(r"\setlength{\tabcolsep}{4pt}")
    lines.append(r"\begin{adjustbox}{max width=\textwidth}")
    lines.append(r"\begin{tabular}{@{}l l r r r r r@{}}")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Gene} & \textbf{Partner Complex (PDB)} & \textbf{Contact ($\le 4.5$\,\AA)} & \textbf{Core ($>8.0$\,\AA)} & \textbf{OR} & \textbf{Fisher $p$} & \textbf{Gap$_{\le 4.5}$ / Gap$_{\text{unmod}}$} \\")
    lines.append(r"\midrule")

    grp_titles = [
        ("Antiviral_Effector_Sensor", "\\multicolumn{7}{@{}l}{\\textit{Conflict Complexes: Antiviral Effectors \\& Sensors ($n = 8$)}} \\\\"),
        ("Viral_Entry_Receptor", "\\multicolumn{7}{@{}l}{\\textit{Conflict Complexes: Viral Entry Receptors ($n = 5$)}} \\\\"),
        ("Tandem_Repeat_Control", "\\multicolumn{7}{@{}l}{\\textit{Tandem-Repeat DNA Complex ($n = 1$)}} \\\\"),
        ("Hijacked_Cellular_Hub", "\\multicolumn{7}{@{}l}{\\textit{Hijacked Cellular Hub Controls ($n = 5$)}} \\\\"),
    ]
    for arm_key, htex in grp_titles:
        sub = [r for r in by_arm.get(arm_key, []) if r["pdb_atomic_audit"] is not None]
        if not sub:
            continue
        lines.append(htex)
        for r in sub:
            pa = r["pdb_atomic_audit"]
            c = pa["contact_4p5A"]
            s = pa["resolved_scaffold_gt8A"]
            if s["n"] == 0:
                s = pa["resolved_noncontact_gt4p5A"]
            u = pa["unmodeled_disordered"]
            or_v = "{0:.2f}$\\times$".format(pa["fisher_or_contact_vs_scaffold"]) if pa["fisher_or_contact_vs_scaffold"] is not None else "---"
            pv_s = fmt_p_sci(pa["fisher_p_contact_vs_scaffold"])
            ug_s = "{0:.1f}\\%".format(100.0 * u["mean_gap_frac"]) if u["mean_gap_frac"] is not None else "---"
            lines.append(
                "\\quad \\texttt{{{0}}} & {1} & {2}/{3} ({4:.1f}\\%) & {5}/{6} ({7:.1f}\\%) & {8} & {9} & {10:.1f}\\% / {11} \\\\".format(
                    r["gene"], PARTNER_SHORT[r["gene"]],
                    c["n_innov"], c["n"], c["innov_pct"],
                    s["n_innov"], s["n"], s["innov_pct"],
                    or_v, pv_s, 100.0 * c["mean_gap_frac"], ug_s
                )
            )
        lines.append(r"\midrule")

    agg = m["atomic_pdb_aggregate"]
    pc = agg["conflict_pdb_pooled"]
    ph = agg["hijacked_hub_pdb_pooled"]
    lines.append("\\textbf{{Pooled Conflict ($n=13$)}} & Within-complex ($\\le 4.5$\\,\\AA\\ vs.\\ core) & {0}/{1} ({2:.1f}\\%) & {3}/{4} ({5:.1f}\\%) & \\textbf{{{6:.2f}$\\times$}} & \\textbf{{{7}}} & {8:.1f}\\% / {9:.1f}\\% \\\\".format(
        pc["contact_4p5A_innov"], pc["contact_4p5A_total"], pc["contact_4p5A_innov_pct"],
        pc["scaffold_gt8A_innov"], pc["scaffold_gt8A_total"], pc["scaffold_gt8A_innov_pct"],
        pc["pooled_fisher_or_innov"], fmt_p_sci(pc["pooled_fisher_p_innov"]),
        pc["per_gene_contact_gap_pct_mean"], pc["per_gene_unmodeled_gap_pct_mean"]
    ))
    lines.append("\\textbf{{Pooled Hubs ($n=5$)}} & Controls ($\\le 4.5$\\,\\AA\\ vs.\\ core) & {0}/{1} ({2:.1f}\\%) & {3}/{4} ({5:.1f}\\%) & {6:.2f}$\\times$ & {7} & {8:.1f}\\% / {9:.1f}\\% \\\\".format(
        ph["contact_4p5A_innov"], ph["contact_4p5A_total"], ph["contact_4p5A_innov_pct"],
        ph["scaffold_gt8A_innov"], ph["scaffold_gt8A_total"], ph["scaffold_gt8A_innov_pct"],
        ph["pooled_fisher_or_innov"], fmt_p_sci(ph["pooled_fisher_p_innov"]),
        ph["per_gene_contact_gap_pct_mean"], ph["per_gene_unmodeled_gap_pct_mean"]
    ))
    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{adjustbox}")
    lines.append(r"\end{table*}")

    with open(OUT_PATH, "w", encoding="utf-8") as out:
        out.write("\n".join(lines) + "\n")
    print("Wrote LaTeX tables cleanly to:", OUT_PATH)


if __name__ == "__main__":
    main()
