#!/usr/bin/env python3
"""
Orchestrate end-to-end HyPhy-3 inference on a single Pilot Quartet gene on silverback:
  1. WALRUS end-to-end from raw unaligned CDS (--coevolution --keep-pseudogenes, emitting JSON, HTML, MAP MSA, Stockholm PP, POA GFA)
  2. S-PRIME 5-property biophysical selection + Ohmic network on the co-estimated MAP codon MSA (after replacing any NNN-masked pseudogene stops with --- for strict sense-codon PRIME if needed)
  3. Illustrative legacy foil (MAFFT protein-guided back-translated codon alignment + standard FEL & MEME without POA marginalization or reticulation routing)
"""

import os
import subprocess
import sys
import time

CODON_STOPS = {"TAA", "TAG", "TGA"}
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


def read_fa(path):
    recs = []
    h, s = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    recs.append((h, "".join(s)))
                h, s = line[1:], []
            else:
                s.append(line)
    if h is not None:
        recs.append((h, "".join(s)))
    return recs


def sanitize_sense_msa(in_msa, out_msa):
    recs = read_fa(in_msa)
    with open(out_msa, "w") as out:
        for h, s in recs:
            s = s.upper()
            cods = []
            for k in range(0, len(s) - 2, 3):
                c = s[k:k + 3]
                if c in CODON_STOPS:
                    cods.append("NNN")
                else:
                    cods.append(c)
            out.write(f">{h}\n{''.join(cods)}\n")


def run_mafft_codon_foil(raw_fa, out_codon_fa, threads):
    """Standard legacy workflow: drop premature-stop pseudogenes, translate to AA, align with MAFFT --auto, back-translate to codons."""
    recs = read_fa(raw_fa)
    intact_nt = {}
    aa_path = out_codon_fa + ".tmp.aa.fa"
    aa_aln_path = out_codon_fa + ".tmp.aa.aln.fa"
    with open(aa_path, "w") as out:
        for h, s in recs:
            s = s.upper()
            # Trim trailing stop
            if len(s) >= 3 and s[-3:] in CODON_STOPS:
                s = s[:-3]
            if len(s) % 3 != 0:
                s = s[: len(s) - (len(s) % 3)]
            aa = "".join(CODON_TABLE.get(s[i:i + 3], "X") for i in range(0, len(s), 3))
            if "*" in aa:
                # Standard codon pipelines reject internal stops and drop pseudogene taxa
                continue
            intact_nt[h] = s
            out.write(f">{h}\n{aa}\n")

    subprocess.run(
        ["mafft", "--auto", "--thread", str(threads), "--quiet", aa_path],
        stdout=open(aa_aln_path, "w"),
        check=True,
    )
    aa_aln = read_fa(aa_aln_path)
    with open(out_codon_fa, "w") as out:
        for h, aseq in aa_aln:
            nt = intact_nt[h]
            pos = 0
            cods = []
            for ch in aseq:
                if ch == "-":
                    cods.append("---")
                else:
                    cods.append(nt[pos:pos + 3])
                    pos += 3
            out.write(f">{h}\n{''.join(cods)}\n")
    for p in (aa_path, aa_aln_path):
        if os.path.exists(p):
            os.remove(p)


def run_cmd(cmd, label):
    t0 = time.time()
    print(f"[{label}] Running: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True)
    print(f"[{label}] Completed in {time.time() - t0:.1f}s", flush=True)


def main():
    gene = sys.argv[1]
    work_dir = sys.argv[2]
    hyphy3 = sys.argv[3]
    threads = int(sys.argv[4]) if len(sys.argv) > 4 else 16

    raw_fa = os.path.join(work_dir, "raw_cds", f"{gene}_unaligned.fa")
    walrus_dir = os.path.join(work_dir, "results", "walrus")
    prime_dir = os.path.join(work_dir, "results", "prime")
    foil_dir = os.path.join(work_dir, "results", "legacy_foil")
    for d in (walrus_dir, prime_dir, foil_dir):
        os.makedirs(d, exist_ok=True)

    w_json = os.path.join(walrus_dir, f"{gene}.walrus.json")
    w_html = os.path.join(walrus_dir, f"{gene}.walrus.html")
    w_msa = os.path.join(walrus_dir, f"{gene}.walrus.msa.fa")
    w_sto = os.path.join(walrus_dir, f"{gene}.walrus.sto")
    w_poa = os.path.join(walrus_dir, f"{gene}.walrus.gfa")

    # 1. WALRUS end-to-end with --coevolution and --keep-pseudogenes
    run_cmd([
        hyphy3, "walrus",
        "-u", raw_fa,
        "--coevolution",
        "--keep-pseudogenes",
        "-o", w_json,
        "--html", w_html,
        "--out-msa", w_msa,
        "--out-stockholm", w_sto,
        "--out-poa", w_poa,
        "--threads", str(threads),
    ], f"{gene}:WALRUS")

    # 2. S-PRIME (5PROP + Ohmic) on WALRUS co-estimated codon MSA
    sense_msa = os.path.join(prime_dir, f"{gene}.walrus.sense.fa")
    sanitize_sense_msa(w_msa, sense_msa)
    p_json = os.path.join(prime_dir, f"{gene}.prime.json")
    p_html = os.path.join(prime_dir, f"{gene}.prime.html")
    run_cmd([
        hyphy3, "prime",
        "--alignment", sense_msa,
        "--tree", "auto",
        "--ohmic",
        "--property-set", "5PROP",
        "--threads", str(threads),
        "--output", p_json,
        "--html", p_html,
        "--no-progress",
    ], f"{gene}:PRIME")

    # 3. Legacy foil (MAFFT codon alignment + standalone MEME & FEL without POA/reticulation/pseudogene rescue)
    if subprocess.call(["which", "mafft"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0:
        m_msa = os.path.join(foil_dir, f"{gene}.mafft_codon.fa")
        run_mafft_codon_foil(raw_fa, m_msa, threads)
        run_cmd([
            hyphy3, "meme",
            "--alignment", m_msa,
            "--tree", "auto",
            "--threads", str(threads),
            "--output", os.path.join(foil_dir, f"{gene}.mafft.MEME.json"),
        ], f"{gene}:MAFFT+MEME")
        run_cmd([
            hyphy3, "fel",
            "--alignment", m_msa,
            "--tree", "auto",
            "--threads", str(threads),
            "--output", os.path.join(foil_dir, f"{gene}.mafft.FEL.json"),
        ], f"{gene}:MAFFT+FEL")


if __name__ == "__main__":
    main()
