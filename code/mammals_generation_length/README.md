# Mammalian CYTB composition and generation length

This directory contains the analysis used to test whether mammalian generation
length is associated with the CYTB amino-acid composition metric:

Pro / (Phe + LeuTT)

where LeuTT denotes leucine encoded by the TTA and TTG codons.

## Analysis

The workflow:

1. extracts mammalian CYTB records from MIDORI2;
2. performs sequence quality control;
3. calculates Pro/(Phe+LeuTT) for each species;
4. matches species to the mammalian generation-length dataset;
5. matches the resulting species to a TimeTree mammalian phylogeny;
6. fits PGLS models with Pagel's lambda estimated by maximum likelihood;
7. generates the Figure 3A panel.

The primary model is:

    Pro/(Phe+LeuTT) ~ log10(generation length in days)

The primary dataset contains 3000 mammalian species.

Current primary result:

- beta = 0.029225
- 95% CI = 0.007342 to 0.051108
- p = 0.008874
- Pagel's lambda = 0.964162

Sensitivity analyses are also performed using:

- direct/observed generation-length estimates only;
- sequences without ambiguous nucleotides;
- both restrictions simultaneously;
- untransformed generation length.

## Files

- scripts/01_prepare_mammals.py - sequence processing, QC and dataset assembly
- scripts/02_run_pgls.py - phylogenetic generalized least-squares analysis
- scripts/03_plot_figure.py - Figure 3A generation
- tests/test_metric.py - basic test of the Pro/(Phe+LeuTT) calculation
- run_all.sh - complete workflow

Processed data and statistical results are stored in:

    data/mammals_generation_length/

The final figure is stored in:

    figures/Fig3A_mammals.pdf

## Input data

The workflow requires three source files:

- MIDORI2_LONGEST_NUC_GB270_Cytb_RAW.fasta.gz
- Generation+Lenght+for+Mammals.xlsx
- Mammalia_species.nwk

By default they are expected under:

    data/mammals_generation_length/raw/

Raw external datasets are not included in this repository.

Alternative input paths can be supplied through the environment variables
MIDORI, GENERATION_LENGTH, and TREE.

## Requirements

Install the Python dependencies:

    python -m pip install -r code/mammals_generation_length/requirements.txt

Run the complete analysis from the repository root:

    bash code/mammals_generation_length/run_all.sh
