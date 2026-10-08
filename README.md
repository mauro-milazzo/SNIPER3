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

## 📖 How to Use SNIPER3

### 1. Generating a Template Sequence

1. **Select an Input Method:**
   * **Gene Symbol + Exon Number:** Enter a gene symbol (e.g., `CFTR`) and select the desired exon number. SNIPER3 automatically fetches RefSeq transcripts, maps the exon boundaries, and capitalizes exonic sequence while keeping intronic flanking sequence in lowercase.
   * **Genomic Coordinates (GRCh38/hg38):** Enter exact coordinates (`Chromosome`, `Start`, and `End`).
   * **rs-Number (dbSNP ID):** Enter a dbSNP ID (e.g., `rs113993960`) to automatically resolve its position and center the template around the variant.
2. **Configure Parameters (Sidebar):**
   * **Flanking Region & Padding:** Adjust upstream/downstream intronic sequence length and target padding around exons.
   * **Variant Filtering:** Toggle SNP highlighting, select the database (**gnomAD v4** or **UCSC snp151**), and set MAF % thresholds and minimum sample size (AN) to filter out rare artifacts.
   * **Segmentation:** Enable auto-splitting for large target regions (>600 bp) into overlapping sub-amplicons.
3. **Generate:** Click **🚀 Generate Sequence Templates**.

---

### 2. Understanding the Interactive Output

* **Sequence Visualizer:**
  * **`[ ]` Blue Brackets:** Mark the target region + padding.
  * **UPPERCASE vs lowercase:** UPPERCASE bases represent the target/exonic region; lowercase bases represent flanking/intronic regions.
  * **<span style="color:red; font-weight:bold;">Red Underlined Bases</span>:** Indicate flagged population variants (SNPs, MNVs, or indels) passing your MAF and sample size thresholds.
  * **Interactive Tooltips & Links:** Hover over any underlined variant to inspect its exact Allele Frequency (AF %) and total sample size (AN). Click on the variant to open its dedicated page on gnomAD or NCBI dbSNP.

---

### 3. Using Output in Primer3 / Primer3Plus

1. Click **📋 Copy Sequence** under your desired fragment. The text area will flash green to confirm the sequence is copied to your clipboard.
2. Open **[Primer3Web](https://primer3.ut.ee/)** or **[Primer3Plus](https://www.primer3plus.com/)**.
3. Paste the copied sequence directly into the **`SEQUENCE_TEMPLATE`** input field.
4. **How Primer3 Interprets SNIPER3 Formatting:**
   * **Brackets `[ ]`:** Define the target region that primers must flank.
   * **Angle Brackets `< >`:** Tell Primer3 to **avoid placing primer binding sites** on these positions (neutralizing SNP mismatch risks).
   * **Capitalization:** Useful when configuring Primer3 parameters like `Included Region` or setting specific conditions for exonic boundaries.
5. Click **Pick Primers** in Primer3 to receive optimal primer pairs that completely avoid high-frequency population variants!
