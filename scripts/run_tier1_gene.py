#!/usr/bin/env python3
"""
Orchestrate end-to-end HyPhy-3 inference on a single Tier-1 gene on silverback:
  1. WALRUS end-to-end from raw unaligned CDS (--coevolution --keep-pseudogenes,
     emitting JSON, standalone HTML, MAP MSA, Stockholm PP, POA GFA)
  2. S-PRIME 5-property biophysical selection + Ohmic network on the co-estimated
     MAP codon MSA (replacing any stop codons in admitted pseudogenes with NNN)
     using conservative Empirical Bayes Ridge shrinkage (--shrinkage 0.5).
"""

import os
import subprocess
import sys
import time

CODON_STOPS = {"TAA", "TAG", "TGA"}


def read_fa(path):
    recs = []
    h, s = None, []
    with open(path, "r", encoding="utf-8") as f:
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
    with open(out_msa, "w", encoding="utf-8") as out:
        for h, s in recs:
            s = s.upper()
            cods = []
            for k in range(0, len(s) - 2, 3):
                c = s[k:k + 3]
                if c in CODON_STOPS:
                    cods.append("NNN")
                else:
                    cods.append(c)
            out.write(">{0}\n{1}\n".format(h, "".join(cods)))


def run_cmd(cmd, label):
    t0 = time.time()
    print("[{0}] Running: {1}".format(label, " ".join(cmd)), flush=True)
    subprocess.run(cmd, check=True)
    print("[{0}] Completed in {1:.1f}s".format(label, time.time() - t0), flush=True)


def main():
    gene = sys.argv[1]
    work_dir = sys.argv[2]
    hyphy3 = sys.argv[3]
    threads = int(sys.argv[4]) if len(sys.argv) > 4 else 32

    raw_fa = os.path.join(work_dir, "raw_cds", "{0}_unaligned.fa".format(gene))
    walrus_dir = os.path.join(work_dir, "results", "walrus")
    prime_dir = os.path.join(work_dir, "results", "prime")
    for d in (walrus_dir, prime_dir):
        os.makedirs(d, exist_ok=True)

    w_json = os.path.join(walrus_dir, "{0}.walrus.json".format(gene))
    w_html = os.path.join(walrus_dir, "{0}.walrus.html".format(gene))
    w_msa = os.path.join(walrus_dir, "{0}.walrus.msa.fa".format(gene))
    w_sto = os.path.join(walrus_dir, "{0}.walrus.sto".format(gene))
    w_poa = os.path.join(walrus_dir, "{0}.walrus.gfa".format(gene))

    # 1. WALRUS end-to-end with --coevolution and --keep-pseudogenes
    if not (os.path.exists(w_json) and os.path.getsize(w_json) > 1000 and os.path.exists(w_msa)):
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
        ], "{0}:WALRUS".format(gene))
    else:
        print("[{0}:WALRUS] Reusing existing outputs.".format(gene), flush=True)

    # 2. S-PRIME (5PROP + Ohmic, conservative --shrinkage 0.5) on WALRUS co-estimated codon MSA
    sense_msa = os.path.join(prime_dir, "{0}.walrus.sense.fa".format(gene))
    sanitize_sense_msa(w_msa, sense_msa)
    p_json = os.path.join(prime_dir, "{0}.prime.json".format(gene))
    if not (os.path.exists(p_json) and os.path.getsize(p_json) > 1000):
        run_cmd([
            hyphy3, "prime",
            "--alignment", sense_msa,
            "--tree", "auto",
            "--ohmic",
            "--property-set", "5PROP",
            "--shrinkage", "0.5",
            "--threads", str(threads),
            "--output", p_json,
            "--no-progress",
        ], "{0}:PRIME".format(gene))
    else:
        print("[{0}:PRIME] Reusing existing output.".format(gene), flush=True)


if __name__ == "__main__":
    main()
