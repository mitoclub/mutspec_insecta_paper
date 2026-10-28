#!/usr/bin/env python3
"""Prepare the mammalian CYTB / generation-length PGLS dataset."""

from __future__ import annotations

import argparse
import gzip
import re
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import Phylo


STOP_CODONS_MT2 = {"TAA", "TAG", "AGA", "AGG"}
PRO_CODONS = {"CCT", "CCC", "CCA", "CCG"}
PHE_CODONS = {"TTT", "TTC"}
LEU_TT_CODONS = {"TTA", "TTG"}

DIRECT_GL_SOURCES = {
    "GMA",
    "Rspan-AFR(SM+Gest)",
    "Rspan-AFB",
    "Rspan-ASMmales",
}

TAXON_TOKEN_RE = re.compile(r"([^_]+)_(.+)_(\d+)$")


def fasta_records_gz(path: Path):
    """Yield (header, sequence) records from a gzipped FASTA."""
    with gzip.open(path, "rt", errors="replace") as handle:
        header = None
        seq = []
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq).upper()
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header is not None:
            yield header, "".join(seq).upper()


def parse_midori_header(header: str):
    """Return accession and MIDORI rank -> (name, taxid) mapping."""
    accession, taxonomy = (header.split("\t", 1) + [""])[:2]
    fields = {}
    for token in taxonomy.split(";"):
        match = TAXON_TOKEN_RE.match(token)
        if match:
            rank, name, taxid = match.groups()
            fields[rank] = (name.replace("_", " "), int(taxid))
    return accession, fields


def cytb_metric(sequence: str):
    """Calculate sequence QC metrics and Pro/(Phe+LeuTT).

    MIDORI2 mitochondrial gene records are treated in the supplied coding
    orientation/frame. Internal stops are checked under vertebrate mitochondrial
    genetic code (NCBI table 2).
    """
    codons = [sequence[i:i + 3] for i in range(0, len(sequence) - 2, 3)]

    internal_stops = sum(codon in STOP_CODONS_MT2 for codon in codons[:-1])
    pro = sum(codon in PRO_CODONS for codon in codons)
    phe = sum(codon in PHE_CODONS for codon in codons)
    leu_tt = sum(codon in LEU_TT_CODONS for codon in codons)
    denominator = phe + leu_tt

    return {
        "sequence_length_nt": len(sequence),
        "ambiguous_nt_count": sum(base not in "ACGT" for base in sequence),
        "internal_stop_count_mt2": internal_stops,
        "pro_count": pro,
        "phe_count": phe,
        "leu_tt_count": leu_tt,
        "denominator_count": denominator,
        "pro_over_phe_plus_leutt": pro / denominator if denominator else None,
    }


def read_generation_length(path: Path):
    if path.suffix.lower() == ".csv":
        frame = pd.read_csv(path)
    else:
        frame = pd.read_excel(
            path,
            sheet_name="Generation_Lenght_for_Mammals",
        )

    required = {
        "Scientific_name",
        "GenerationLength_d",
        "Sources_GL",
        "Order",
        "Family",
        "TaxID",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing generation-length columns: {sorted(missing)}")

    frame = frame.copy()
    frame["Scientific_name"] = frame["Scientific_name"].astype(str).str.strip()

    # Conservative exact-name merge: names duplicated in the ecological source
    # are excluded rather than resolved silently.
    counts = Counter(frame["Scientific_name"])
    duplicate_names = {name for name, n in counts.items() if n > 1}
    unique = frame[~frame["Scientific_name"].isin(duplicate_names)].copy()

    return frame, unique, duplicate_names


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--midori", type=Path, required=True)
    parser.add_argument("--generation-length", type=Path, required=True)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--qc-report", type=Path, required=True)
    parser.add_argument("--min-length", type=int, default=1000)
    args = parser.parse_args()

    gl_all, gl_unique, duplicate_names = read_generation_length(
        args.generation_length
    )
    gl_by_name = gl_unique.set_index("Scientific_name", drop=False)

    mammal_records = 0
    exact_gl_matches = 0
    qc_rows = []

    for header, sequence in fasta_records_gz(args.midori):
        # Fast rejection before parsing the complete taxonomy string.
        if "class_Mammalia_40674" not in header:
            continue

        mammal_records += 1
        accession, taxonomy = parse_midori_header(header)
        if "species" not in taxonomy:
            continue

        species, ncbi_taxid = taxonomy["species"]
        if species not in gl_by_name.index:
            continue

        exact_gl_matches += 1
        metrics = cytb_metric(sequence)

        if metrics["sequence_length_nt"] < args.min_length:
            continue
        if metrics["internal_stop_count_mt2"] != 0:
            continue
        if not metrics["denominator_count"]:
            continue

        gl = gl_by_name.loc[species]

        qc_rows.append({
            "species": species,
            "tree_name": species.replace(" ", "_"),
            "midori_accession": accession,
            "midori_ncbi_taxid": ncbi_taxid,
            **metrics,
            "generation_length_d": float(gl["GenerationLength_d"]),
            "generation_length_source": gl["Sources_GL"],
            "generation_length_direct": gl["Sources_GL"] in DIRECT_GL_SOURCES,
            "generation_table_taxid": gl["TaxID"],
            "order": gl["Order"],
            "family": gl["Family"],
        })

    prepared = pd.DataFrame(qc_rows)

    tree = Phylo.read(str(args.tree), "newick")
    tree_tips = {tip.name for tip in tree.get_terminals()}

    prepared["matched_timetree"] = prepared["tree_name"].isin(tree_tips)
    final = prepared[prepared["matched_timetree"]].copy().reset_index(drop=True)
    final["zero_ambiguous_nt"] = final["ambiguous_nt_count"].eq(0)
    final["log10_generation_length"] = final["generation_length_d"].map(
        __import__("math").log10
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(args.output, index=False)

    report = "\n".join([
        "MAMMAL CYTB PREPROCESSING QC",
        "============================",
        f"MIDORI mammal CYTB records: {mammal_records}",
        f"Generation-length rows: {len(gl_all)}",
        (
            "Duplicate Scientific_name values in generation-length table "
            f"(excluded): {len(duplicate_names)}"
        ),
        f"Exact MIDORI x generation-length name matches before sequence QC: {exact_gl_matches}",
        f"After sequence QC: {len(prepared)}",
        f"Exact matches to TimeTree tips: {len(final)}",
        f"Not matched to TimeTree tips: {len(prepared) - len(final)}",
        (
            "Final species with zero ambiguous nucleotides: "
            f"{int(final['zero_ambiguous_nt'].sum())}"
        ),
        (
            "Final species with direct/observed-calculated generation length: "
            f"{int(final['generation_length_direct'].sum())}"
        ),
        "",
        "Primary QC:",
        f"- sequence length >= {args.min_length} nt",
        "- no internal stops under vertebrate mitochondrial code",
        "- Phe + LeuTT > 0",
        "- exact species-name matching only",
    ])
    args.qc_report.parent.mkdir(parents=True, exist_ok=True)
    args.qc_report.write_text(report + "\n", encoding="utf-8")

    print(report)


if __name__ == "__main__":
    main()
