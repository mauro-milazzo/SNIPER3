# 🧬 SNIPER3

**S**NP-**N**eutralizing **I**ntelligent **P**rimer3 **E**xon & **R**egion Template Generator

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io)
[![Genome Assembly](https://img.shields.io/badge/Genome-GRCh38%2Fhg38-blue.svg)](https://genome.ucsc.edu/)
[![Python](https://img.shields.io/badge/Python-3.9%2B-green.svg)](https://www.python.org/)

**SNIPER3** is a web application designed to simplify and standardize the preparation of DNA sequence templates for **Primer3Web** and **Primer3Plus**. By directly leveraging genomic APIs, it automates flanking sequence retrieval, exonic capitalization, target region segmentation, and population-scale SNP masking (`< >`) using either **gnomAD v4** or **UCSC dbSNP (snp151)**.

---

## ✨ Key Features

* **GRCh38/hg38 Native Integration:** Fetches genomic reference sequences directly via the UCSC Genome Browser REST API.
* **Flexible Input Methods:**
  * **Gene Symbol + Exon Number:** Retrieves RefSeq transcripts (e.g., `CFTR` Exon 11) with automatic intron/exon mapping.
  * **Genomic Coordinates:** Direct input of target coordinates (`chr:start-end`).
  * **dbSNP rsID:** Resolves single-nucleotide variants directly to exact genomic coordinates (e.g., `rs113993960`).
* **Automated SNP Neutralization (`< >` Masking):**
  * Choose between **gnomAD v4** (GraphQL API with >800k exomes/genomes) or **UCSC snp151** (dbSNP).
  * Filter variants dynamically by **Minor Allele Frequency (MAF %)** and **Minimum Allele Count (AN)** to avoid common PCR mismatch artifacts.
  * Multi-nucleotide variants (MNVs), indels, and deletions are fully spanned within SNP brackets.
* **Smart Target Segmentation:** Automatically splits large exons or target regions exceeding maximum amplicon length into overlapping sub-amplicons.
* **Interactive Sequence Visualizer:**
  * Interactive HTML preview highlighting target boundaries `[ ]`, exonic regions (UPPERCASE), and flagged SNPs (<span style="color:red; font-weight:bold;">red underlined</span>).
  * Hover tooltips displaying exact rsIDs, Allele Frequencies (AF %), and sample size metrics (AN).
  * Direct hyperlinks to variant detail pages on gnomAD or NCBI dbSNP.
* **One-Click Clipboard & Export:** Fast interactive copy button with visual feedback and individual fragment `.txt` downloads formatted for Primer3.

---
