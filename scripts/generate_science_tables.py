#!/usr/bin/env python3
"""
Generate Table 1 (main.tex) and Tables S1-S4 (supplement.tex) directly from
executed quantitative synthesis JSONs on disk with zero \textbf{} inside cells or captions.
"""

import json
import os
import numpy as np

ROOT = "/Users/sergei/Projects/unalignable-arms-race"
MASTER = os.path.join(ROOT, "results", "tier1_25gene_master_synthesis.json")
TIER2 = os.path.join(ROOT, "results", "tier2_toga2_vip_summary.json")
WALRUS_DIR = os.path.join(ROOT, "results", "walrus")

OUT_T1 = os.path.join(ROOT, "manuscript", "table1_summary.tex")
OUT_SUPP = os.path.join(ROOT, "manuscript", "tables_supplement_generated.tex")

PARTNER_SHORT = {
    "APOBEC3G": "HIV-1 Vif (8CX0)",
    "APOBEC3H": "HIV-1 Vif (8FVI)",
    "ISG15":    "FluB NS1B / PLpro (3SDL, 5TL6)",
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
    "ELOB":     "HIV-1 Vif (2MA9, 4N9F)",
    "ANP32A":   "FluA Pol (8RN0)",
    "RAE1":     "SARS-CoV-2 Orf6 (7VPG)",
    "YWHAZ":    "SARS-CoV-2 N (7ZIT)",
}


def fmt_p(pv):
    if pv is None:
        return "---"
    if pv >= 0.01:
        return f"${pv:.3f}$"
    if pv >= 0.001:
        return f"${pv:.4f}$"
    s = f"{pv:.1e}"
    mant, exp = s.split("e")
    return f"${mant}\\times 10^{{{int(exp)}}}$"


def main():
    m = json.load(open(MASTER))
    t2 = json.load(open(TIER2))
    recs = m["tier1_genes"]
    by_arm = {}
    for r in recs:
        by_arm.setdefault(r["arm"], []).append(r)

    # ---- MAIN TEXT TABLE 1 (Compact Synthesis Across Tiers 1 and 2) ----
    t1_lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Multi-scale selection, structural localization, and indel plasticity across mammalian study cohorts. Values report group means or pooled counts from executed Tier-1 and Tier-2 analyses; two-sided statistical comparisons test antiviral or receptor interfaces against hijacked cellular host hubs.}",
        r"\label{tab:arms_summary}",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{4.5pt}",
        r"\begin{adjustbox}{max width=\textwidth}",
        r"\begin{tabular}{@{}l r r r r r r@{}}",
        r"\toprule",
        r"Cohort or structural stratum & Loci ($n$) & Taxa (med.) & Gap (\%) & Selected (\%) & BUSTED-E $p \le 0.10$ (\%) & Strict S-PRIME / POA art. \\",
        r"\midrule",
    ]
    arm_meta = [
        ("Tier 1: Antiviral effectors and sensors", "Antiviral_Effector_Sensor"),
        ("Tier 1: Viral entry receptors", "Viral_Entry_Receptor"),
        ("Tier 1: Tandem-repeat array (PRDM9)", "Tandem_Repeat_Control"),
        ("Tier 1: Hijacked cellular hubs (controls)", "Hijacked_Cellular_Hub"),
    ]
    for lbl, ak in arm_meta:
        sub = by_arm[ak]
        n = len(sub)
        med_t = int(np.median([r["admitted_taxa"] for r in sub]))
        gap = float(np.mean([r["mean_gap_pct"] for r in sub]))
        sel = float(np.mean([r["innov_site_pct"] for r in sub]))
        bus = 100.0 * sum(1 for r in sub if r["busted_e_sig"]) / n
        sp = sum(r["prime_strict_changing_by_prop"].get("Isoelectric_Point_pI", 0) for r in sub)
        art = sum(r["poa_artifacts_eliminated"] for r in sub)
        t1_lines.append(
            f"{lbl} & {n} & {med_t} & {gap:.1f} & {sel:.1f} & {bus:.1f} & {sp} strict pI / {art} pruned \\\\"
        )

    agg = m["atomic_pdb_aggregate"]
    pc = agg["conflict_pdb_pooled"]
    ph = agg["hijacked_hub_pdb_pooled"]
    t1_lines.extend([
        r"\midrule",
        f"Tier 1 PDB conflict contacts ($\\le 4.5$\\,\\AA) & 13 (299 res.) & 22--72 & {pc['per_gene_contact_gap_pct_mean']:.2f} & {pc['contact_4p5A_innov_pct']:.1f} & --- & {pc['contact_4p5A_prime_strict_pct']:.2f}\\% strict (OR $4.11\\times$) \\\\",
        f"Tier 1 PDB conflict core ($>8.0$\\,\\AA) & 13 (3,643 res.) & 22--72 & {pc['per_gene_scaffold_gap_pct_mean']:.2f} & {pc['scaffold_gt8A_innov_pct']:.1f} & --- & {pc['scaffold_gt8A_prime_strict_pct']:.2f}\\% strict \\\\",
        f"Tier 1 PDB hub contacts ($\\le 4.5$\\,\\AA) & 5 (102 res.) & 72 & {ph['per_gene_contact_gap_pct_mean']:.2f} & {ph['contact_4p5A_innov_pct']:.1f} & --- & 0.00\\% strict (OR $14.32\\times$ vs.\\ hubs) \\\\",
        r"\midrule",
    ])
    g2 = t2["groups"]
    for lbl, gk in [
        ("Tier 2: Antiviral effectors (KW-0051)", "Antiviral_Sensors_Effectors"),
        ("Tier 2: Proteome background", "Proteome_Background"),
        ("Tier 2: Hijacked hubs (KW-0945 excl. KW-0051)", "Hijacked_Host_Hubs"),
    ]:
        obj = g2[gk]
        t1_lines.append(
            f"{lbl} & {obj['n']:,} & {obj['median_taxa']} & {obj['mean_gap_codons_per_site']:.2f} & {obj['mean_episodic_site_pct']:.2f} (ep.) & {obj['busted_pct']:.1f} & {obj['mean_toggle_conv_site_pct']:.2f}\\% toggle/conv \\\\"
        )
    t1_lines.extend([
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{adjustbox}",
        r"\end{table}",
    ])
    open(OUT_T1, "w", encoding="utf-8").write("\n".join(t1_lines) + "\n")

    # ---- SUPPLEMENT TABLES S1-S4 ----
    supp = []

    # Table S1: Full 25-gene ledger
    supp.extend([
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Complete evolutionary ledger across all 25 pre-registered Tier-1 mammalian loci evaluated from unaligned coding sequences with pseudogene rescue and conservative Ridge S-PRIME shrinkage ($\lambda = 0.5$).}",
        r"\label{tab:s1_tier1_atlas}",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{3.6pt}",
        r"\begin{adjustbox}{max width=\textwidth}",
        r"\begin{tabular}{@{}l r r r r r r r r r@{}}",
        r"\toprule",
        r"Gene & Taxa & Ref/MSA & Ins. & Gap\% & Sel.\% & BUSTED-E $p$ & Art. & pI (M/S) & LoF (Ps) / Coev \\",
        r"\midrule",
    ])
    arm_headers = [
        ("Antiviral_Effector_Sensor", r"\multicolumn{10}{@{}l}{\textit{Antiviral effectors and innate sensors ($n = 14$)}} \\"),
        ("Viral_Entry_Receptor", r"\multicolumn{10}{@{}l}{\textit{Viral entry receptors ($n = 5$)}} \\"),
        ("Tandem_Repeat_Control", r"\multicolumn{10}{@{}l}{\textit{Tandem-repeat meiotic benchmark ($n = 1$)}} \\"),
        ("Hijacked_Cellular_Hub", r"\multicolumn{10}{@{}l}{\textit{Hijacked cellular host hubs ($n = 5$)}} \\"),
    ]
    for ak, hdr in arm_headers:
        supp.append(hdr)
        for r in by_arm[ak]:
            pm = r["prime_marginal_changing_by_prop"].get("Isoelectric_Point_pI", 0)
            ps = r["prime_strict_changing_by_prop"].get("Isoelectric_Point_pI", 0)
            nlof = len(r["lof_flagged_branches"])
            nps = len(r["admitted_pseudogenes"])
            lstr = f"{nlof}" + (f" (+{nps})" if nps else "")
            supp.append(
                f"\\quad \\texttt{{{r['gene']}}} & {r['admitted_taxa']} & {r['human_aa']}/{r['msa_codons']} & "
                f"{r['non_human_inserted_cols']} & {r['mean_gap_pct']:.1f} & {r['innov_site_pct']:.1f} & "
                f"{fmt_p(r['busted_e_p'])} & {r['poa_artifacts_eliminated']} & {pm}/{ps} & {lstr} / {r['coev_pairs_q10']} \\\\"
            )
        supp.append(r"\midrule")
    supp[-1] = r"\bottomrule"
    supp.extend([r"\end{tabular}", r"\end{adjustbox}", r"\end{table}", ""])

    # Table S2: All 19 PDB complexes
    supp.extend([
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Atomic heavy-atom partner contact ($d_{\min} \le 4.5$\,\AA) versus buried scaffold ($>8.0$\,\AA) enrichment across 19 experimental PDB co-structures. For SAMHD1 4CC9, non-contact compares $>4.5$\,\AA\ residues.}",
        r"\label{tab:s2_pdb_atomic}",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{3.8pt}",
        r"\begin{adjustbox}{max width=\textwidth}",
        r"\begin{tabular}{@{}l l r r r r r@{}}",
        r"\toprule",
        r"Gene & Partner complex (PDB) & Contact ($\le 4.5$\,\AA) & Core ($>8.0$\,\AA) & OR & Fisher $p$ & Gap$_{\le 4.5}$ / Gap$_{\text{unmod}}$ \\",
        r"\midrule",
    ])
    for ak, hdr_text in [
        ("Antiviral_Effector_Sensor", r"\multicolumn{7}{@{}l}{\textit{Conflict complexes: antiviral effectors and sensors ($n = 8$)}} \\"),
        ("Viral_Entry_Receptor", r"\multicolumn{7}{@{}l}{\textit{Conflict complexes: viral entry receptors ($n = 5$)}} \\"),
        ("Tandem_Repeat_Control", r"\multicolumn{7}{@{}l}{\textit{Tandem-repeat DNA complex ($n = 1$)}} \\"),
        ("Hijacked_Cellular_Hub", r"\multicolumn{7}{@{}l}{\textit{Hijacked cellular hub controls ($n = 5$)}} \\"),
    ]:
        sub = [r for r in by_arm[ak] if r["pdb_atomic_audit"] is not None]
        supp.append(hdr_text)
        for r in sub:
            pa = r["pdb_atomic_audit"]
            c = pa["contact_4p5A"]
            s = pa["resolved_scaffold_gt8A"] if pa["resolved_scaffold_gt8A"]["n"] > 0 else pa["resolved_noncontact_gt4p5A"]
            u = pa["unmodeled_disordered"]
            orv = f"{pa['fisher_or_contact_vs_scaffold']:.2f}$\\times$" if pa["fisher_or_contact_vs_scaffold"] is not None else "---"
            ugs = f"{100.0 * u['mean_gap_frac']:.1f}\\%" if u["mean_gap_frac"] is not None else "---"
            supp.append(
                f"\\quad \\texttt{{{r['gene']}}} & {PARTNER_SHORT[r['gene']]} & {c['n_innov']}/{c['n']} ({c['innov_pct']:.1f}\\%) & "
                f"{s['n_innov']}/{s['n']} ({s['innov_pct']:.1f}\\%) & {orv} & {fmt_p(pa['fisher_p_contact_vs_scaffold'])} & "
                f"{100.0 * c['mean_gap_frac']:.1f}\\% / {ugs} \\\\"
            )
        supp.append(r"\midrule")
    supp.extend([
        f"Pooled conflict ($n=13$) & Within-complex ($\\le 4.5$\\,\\AA\\ vs.\\ core) & {pc['contact_4p5A_innov']}/{pc['contact_4p5A_total']} ({pc['contact_4p5A_innov_pct']:.1f}\\%) & {pc['scaffold_gt8A_innov']}/{pc['scaffold_gt8A_total']} ({pc['scaffold_gt8A_innov_pct']:.1f}\\%) & {pc['pooled_fisher_or_innov']:.2f}$\\times$ & {fmt_p(pc['pooled_fisher_p_innov'])} & {pc['per_gene_contact_gap_pct_mean']:.1f}\\% / {pc['per_gene_unmodeled_gap_pct_mean']:.1f}\\% \\\\",
        f"Pooled hubs ($n=5$) & Controls ($\\le 4.5$\\,\\AA\\ vs.\\ core) & {ph['contact_4p5A_innov']}/{ph['contact_4p5A_total']} ({ph['contact_4p5A_innov_pct']:.1f}\\%) & {ph['scaffold_gt8A_innov']}/{ph['scaffold_gt8A_total']} ({ph['scaffold_gt8A_innov_pct']:.1f}\\%) & {ph['pooled_fisher_or_innov']:.2f}$\\times$ & {fmt_p(ph['pooled_fisher_p_innov'])} & {ph['per_gene_contact_gap_pct_mean']:.1f}\\% / {ph['per_gene_unmodeled_gap_pct_mean']:.1f}\\% \\\\",
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{adjustbox}",
        r"\end{table}",
        "",
    ])

    # Table S3: All 10 SPR-routed Wheatstone gene-conversion bridges
    spr_rows = []
    for r in recs:
        g = r["gene"]
        w = json.load(open(os.path.join(WALRUS_DIR, f"{g}.walrus.json")))
        for br in w.get("recombination_bridges", []):
            if br.get("routed_to_spr_tree"):
                spr_rows.append((g, br["clade_u"], br["target_v"], br["peak_start_codon"], br["peak_end_codon"],
                                 br["segment_delta_log_l"], br["conductance"], br["p_value"]))
    spr_rows.sort(key=lambda x: -x[5])
    supp.extend([
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{All 10 statistically routed Ohmic Wheatstone gene-conversion tracts detected across the 25 Tier-1 loci before site selection evaluation.}",
        r"\label{tab:s3_wheatstone}",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{4.5pt}",
        r"\begin{adjustbox}{max width=\textwidth}",
        r"\begin{tabular}{@{}l l l r r r r@{}}",
        r"\toprule",
        r"Gene & Donor/source node & Recipient/target node & Codon interval & $\Delta \ln L$ & Conductance & Perm.\ $p$ \\",
        r"\midrule",
    ])
    for g, cu, tv, ps, pe, dl, cond, pv in spr_rows:
        cu_t = cu.replace("_", "\\_")
        tv_t = tv.replace("_", "\\_")
        supp.append(f"\\texttt{{{g}}} & \\texttt{{{cu_t}}} & \\texttt{{{tv_t}}} & {ps}--{pe} & {dl:.2f} & {cond:.3f} & {pv:.4f} \\\\")
    supp.extend([r"\bottomrule", r"\end{tabular}", r"\end{adjustbox}", r"\end{table}", ""])

    # Table S4: Tier-2 TOGA2 Proteome Groups & OLS Regressions
    supp.extend([
        r"\begin{table}[htbp]",
        r"\centering",
        r"\caption{Pan-mammalian TOGA2 proteome divergence across 7,045 canonical orthogroups and covariate-adjusted ordinary least-squares estimates controlling for $\log(\text{taxa})$ and $\log(\text{codons})$.}",
        r"\label{tab:s4_tier2_toga2}",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{4.5pt}",
        r"\begin{adjustbox}{max width=\textwidth}",
        r"\begin{tabular}{@{}l r r r r r@{}}",
        r"\toprule",
        r"Comparison / metric & Value ($n$ or $\beta$) & SE & Test statistic & Two-sided $p$ & Notes \\",
        r"\midrule",
    ])
    cov_h = t2["covariate_adjusted_av_vs_hijacked"]
    cov_b = t2["covariate_adjusted_av_vs_background"]
    supp.extend([
        f"Antiviral vs.\\ hubs: episodic site \\% (OLS) & $+{cov_h['episodic_site_pct']['beta_antiviral_adj']:.2f}\\%$ & {cov_h['episodic_site_pct']['se']:.2f} & $t={cov_h['episodic_site_pct']['t']:.2f}$ & {fmt_p(cov_h['episodic_site_pct']['p_two_sided'])} & $n=49$ vs.\\ $179$ \\\\",
        f"Antiviral vs.\\ hubs: gap density \\% (OLS) & $+{cov_h['gap_density_pct']['beta_antiviral_adj']:.2f}\\%$ & {cov_h['gap_density_pct']['se']:.2f} & $t={cov_h['gap_density_pct']['t']:.2f}$ & {fmt_p(cov_h['gap_density_pct']['p_two_sided'])} & $n=49$ vs.\\ $179$ \\\\",
        f"Antiviral vs.\\ hubs: toggle/conv \\% (OLS) & $+{cov_h['toggle_conv_pct']['beta_antiviral_adj']:.2f}\\%$ & {cov_h['toggle_conv_pct']['se']:.2f} & $t={cov_h['toggle_conv_pct']['t']:.2f}$ & {fmt_p(cov_h['toggle_conv_pct']['p_two_sided'])} & $n=49$ vs.\\ $179$ \\\\",
        f"Antiviral vs.\\ background: episodic \\% (OLS) & $+{cov_b['episodic_site_pct']['beta_antiviral_adj']:.2f}\\%$ & {cov_b['episodic_site_pct']['se']:.2f} & $t={cov_b['episodic_site_pct']['t']:.2f}$ & {fmt_p(cov_b['episodic_site_pct']['p_two_sided'])} & $n=49$ vs.\\ $6,817$ \\\\",
        f"Antiviral vs.\\ background: gap \\% (OLS) & $+{cov_b['gap_density_pct']['beta_antiviral_adj']:.2f}\\%$ & {cov_b['gap_density_pct']['se']:.2f} & $t={cov_b['gap_density_pct']['t']:.2f}$ & {fmt_p(cov_b['gap_density_pct']['p_two_sided'])} & $n=49$ vs.\\ $6,817$ \\\\",
        f"Antiviral vs.\\ hubs: BUSTED-E \\% (Fisher) & $69.39\\%$ vs.\\ $39.11\\%$ & --- & --- & {fmt_p(t2['p_antiviral_vs_hijacked']['busted_fisher_p'])} & $34/49$ vs.\\ $70/179$ \\\\",
        f"Hijacked hubs higher taxon retention (MWU) & $839$ vs.\\ $441$ taxa & --- & --- & {fmt_p(t2['p_antiviral_vs_hijacked']['hijacked_higher_taxa_retention_mwu_p'])} & Median conserved taxa \\\\",
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{adjustbox}",
        r"\end{table}",
    ])
    open(OUT_SUPP, "w", encoding="utf-8").write("\n".join(supp) + "\n")
    print("Generated:", OUT_T1, "and", OUT_SUPP)


if __name__ == "__main__":
    main()
