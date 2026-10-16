import math
import re
import time
import requests
import pandas as pd
import streamlit as st

# ==============================================================================
# 1. APPLICATION CONFIGURATION & CONSTANTS
# ==============================================================================
# Defines the maximum coordinate span (in base pairs) permitted for direct coordinate input.
# Safeguards against oversized genomic queries that could overload REST APIs or cause execution timeouts.
MAX_COORDINATE_REGION_SIZE = 20000

# Enforces a mandatory cooldown duration (in seconds) between user sequence generation requests.
# Prevents rate-limit exhaustion (HTTP 429) on public genomic REST endpoints (gnomAD GraphQL & UCSC REST).
COOLDOWN_SECONDS = 10.0

# Configure page-level layout, browser tab title, and favicon icon for Streamlit
st.set_page_config(
    page_title="SNIPER3 - Primer3 Template Generator",
    page_icon="🧬",
    layout="wide"
)

# Render main application header and subtitle explanation
st.title("🧬 SNIPER3")
st.caption("**S**NP-**N**eutralizing **I**ntelligent **P**rimer3 **E**xon & **R**egion Template Generator")

st.markdown(
    "Generate standardized sequence templates formatted for Primer3Web and Primer3Plus using **GRCh38/hg38** coordinates. "
    "Features include large exon auto-splitting, exonic/target sequence capitalization, "
    "and automated SNP filtering (< >) by MAF and population sample size via **gnomAD v4** or **dbSNP (151)**."
)

# ==============================================================================
# 2. STREAMLIT SESSION STATE INITIALIZATION
# ==============================================================================
# Streamlit re-executes the python script top-to-bottom on every user interaction.
# Initializing default state variables ensures data persists across rerun cycles.

if "results" not in st.session_state:
    # Storage container holding sequence fragments, variant tables, and parameters after execution
    st.session_state["results"] = None

if "last_button_click_time" not in st.session_state:
    # Stores the Unix epoch timestamp of the most recent generation run for cooldown calculation
    st.session_state["last_button_click_time"] = 0.0

if "trigger_generate" not in st.session_state:
    # Operational flag used to signal that the sequence generation pipeline should run
    st.session_state["trigger_generate"] = False

def clear_results():
    """
    Callback function attached to UI input controls.
    Clears stale result data whenever the user alters search terms or filter parameters.
    """
    st.session_state["results"] = None


# ==============================================================================
# 3. CACHED EXTERNAL API DATA FETCHING FUNCTIONS
# ==============================================================================
# Annotating network calls with @st.cache_data caches response payloads for 1 hour (3600s),
# significantly speeding up repeated queries and reducing redundant network overhead.

@st.cache_data(ttl=3600, show_spinner=False)
def fetch_refseq_transcripts(gene_symbol: str):
    """
    Fetches RefSeq transcript annotations for a given gene symbol using the UCSC REST API (GRCh38/hg38).

    Args:
        gene_symbol (str): HGNC official gene symbol (e.g., 'CFTR', 'BRCA1').

    Returns:
        list: Sorted list of dictionary objects representing RefSeq transcripts,
              ordered in descending order by exon count so coding NM_ transcripts take precedence.
    """
    gene_clean = gene_symbol.strip().upper()
    
    # Query UCSC search endpoint to obtain chromosome and genomic boundary coordinates
    search_url = f"https://api.genome.ucsc.edu/search?genome=hg38;search={gene_clean}"
    r_search = requests.get(search_url, timeout=10)
    r_search.raise_for_status()
    
    search_data = r_search.json()
    pos_str = None
    
    # Inspect search response structure to isolate coordinate match string
    if 'positionMatches' in search_data:
        for cat in search_data['positionMatches']:
            for match in cat.get('matches', []):
                if 'position' in match:
                    pos_str = match['position']
                    break
            if pos_str: break
    elif 'results' in search_data and len(search_data['results']) > 0:
        res = search_data['results'][0]
        pos_str = res.get('position') or f"{res.get('chrom')}:{res.get('chromStart')}-{res.get('chromEnd')}"

    # Query the ncbiRefSeq track for detailed transcript annotations across the gene locus
    if pos_str:
        match = re.match(r"(chr[0-9XY]+):(\d+)-(\d+)", pos_str, re.IGNORECASE)
        if match:
            c_name, s_pos, e_pos = match.group(1), int(match.group(2)), int(match.group(3))
            track_url = f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=ncbiRefSeq;chrom={c_name};start={s_pos};end={e_pos}"
            r_track = requests.get(track_url, timeout=10)
            r_track.raise_for_status()
            
            refseq_items = r_track.json().get('ncbiRefSeq', [])
            
            # Prioritize primary coding transcripts beginning with 'NM_' matching the query symbol
            gene_transcripts = [
                t for t in refseq_items 
                if t.get('name2', '').upper() == gene_clean and t.get('name', '').startswith('NM_')
            ]
            if not gene_transcripts:
                gene_transcripts = [t for t in refseq_items if t.get('name2', '').upper() == gene_clean]
            
            # Sort transcripts by total exon count descending to feature full-length isoforms
            return sorted(gene_transcripts, key=lambda x: x.get('exonCount', 0), reverse=True)
    return []


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_ucsc_sequence(chrom: str, start_pos: int, end_pos: int):
    """
    Retrieves raw GRCh38/hg38 reference DNA sequence from the UCSC REST API.

    Note:
        Converts 1-based human coordinate notation to UCSC's 0-based start index (start_pos - 1).

    Args:
        chrom (str): Chromosome identifier (e.g., '7', 'X').
        start_pos (int): 1-based start genomic position.
        end_pos (int): 1-based end genomic position.

    Returns:
        str: Uppercase nucleotide DNA sequence string.
    """
    c_clean = str(chrom).lower().replace('chr', '')
    if c_clean == '23': c_clean = 'X'
    elif c_clean == '24': c_clean = 'Y'
    
    ucsc_seq_url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr{c_clean};start={start_pos-1};end={end_pos}"
    r_seq = requests.get(ucsc_seq_url, timeout=15)
    r_seq.raise_for_status()
    
    res_json = r_seq.json()
    if 'dna' in res_json:
        return res_json['dna'].upper()
    else:
        raise RuntimeError(f"Failed to fetch reference sequence from UCSC for chr{c_clean}:{start_pos}-{end_pos}.")


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_ucsc_snps(chrom: str, start_pos: int, end_pos: int):
    """
    Retrieves variant records from the dbSNP build 151 track (snp151) via UCSC REST API.

    Args:
        chrom (str): Chromosome identifier.
        start_pos (int): 1-based start coordinate.
        end_pos (int): 1-based end coordinate.

    Returns:
        list: Raw list of variant dictionaries returned by UCSC track API.
    """
    c_clean = str(chrom).lower().replace('chr', '')
    if c_clean == '23': c_clean = 'X'
    elif c_clean == '24': c_clean = 'Y'
    
    snp_url = f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=snp151;chrom=chr{c_clean};start={start_pos-1};end={end_pos}"
    r_snp = requests.get(snp_url, timeout=15)
    r_snp.raise_for_status()
    
    return r_snp.json().get('snp151', [])


@st.cache_data(ttl=3600, show_spinner=False)
def fetch_gnomad_v4_snps(chrom: str, start_pos: int, end_pos: int):
    """
    Queries gnomAD v4 GraphQL endpoint in a single call to obtain population allele counts.
    Computes exact weighted Joint Allele Frequency (exome AC + genome AC) / (exome AN + genome AN).

    Args:
        chrom (str): Chromosome identifier.
        start_pos (int): 1-based start position.
        end_pos (int): 1-based end position.

    Returns:
        list: Parsed dictionaries containing variant details, calculated AF, and sample size.
    """
    c_clean = str(chrom).lower().replace('chr', '')
    if c_clean == '23': c_clean = 'X'
    elif c_clean == '24': c_clean = 'Y'

    url = "https://gnomad.broadinstitute.org/api"
    
    # GraphQL query string selecting essential allele frequency fields in gnomAD v4
    query = """
    query GnomadRegionVariants($chrom: String!, $start: Int!, $stop: Int!, $dataset: DatasetId!) {
      region(chrom: $chrom, start: $start, stop: $stop, reference_genome: GRCh38) {
        variants(dataset: $dataset) {
          variant_id
          rsids
          pos
          ref
          alt
          genome {
            ac
            an
          }
          exome {
            ac
            an
          }
        }
      }
    }
    """
    variables = {
        "chrom": c_clean,
        "start": int(start_pos),
        "stop": int(end_pos),
        "dataset": "gnomad_r4"
    }

    r = requests.post(url, json={"query": query, "variables": variables}, timeout=20)
    r.raise_for_status()

    data = r.json()
    if "errors" in data:
        raise RuntimeError(f"gnomAD GraphQL API returned errors: {data['errors'][0].get('message', 'Unknown error')}")

    raw_variants = data.get("data", {}).get("region", {}).get("variants", [])
    
    parsed_variants = []
    for v in raw_variants:
        # Extract allele counts (AC) and total alleles called (AN) for exomes and genomes
        exome_ac = v.get("exome", {}).get("ac") or 0 if v.get("exome") else 0
        exome_an = v.get("exome", {}).get("an") or 0 if v.get("exome") else 0

        genome_ac = v.get("genome", {}).get("ac") or 0 if v.get("genome") else 0
        genome_an = v.get("genome", {}).get("an") or 0 if v.get("genome") else 0

        # Combine counts across exomes and genomes for true joint population statistics
        total_ac = exome_ac + genome_ac
        total_an = exome_an + genome_an

        joint_af = (total_ac / total_an) if total_an > 0 else 0.0

        rsids = v.get("rsids") or []
        rs_name = rsids[0] if rsids else v["variant_id"]

        parsed_variants.append({
            "name": rs_name,
            "variant_id": v["variant_id"],
            "pos": v["pos"],
            "ref": v["ref"],
            "alt": v["alt"],
            "maf": joint_af,
            "joint_af": joint_af,
            "total_n": total_an
        })
    return parsed_variants


@st.cache_data(ttl=3600, show_spinner=False)
def resolve_rsid(rs_id: str):
    """
    Resolves a dbSNP rs-number (e.g., 'rs113993960') to exact GRCh38/hg38 genomic coordinates
    using the UCSC Search REST service.

    Returns:
        tuple: (chrom, start_pos, end_pos, formatted_rsid) or (None, None, None, formatted_rsid) if not found.
    """
    rs_clean = rs_id.strip().lower()
    if not rs_clean.startswith('rs'):
        rs_clean = f"rs{rs_clean}"

    try:
        search_url = f"https://api.genome.ucsc.edu/search?genome=hg38;search={rs_clean}"
        r_search = requests.get(search_url, timeout=10)
        
        if r_search.ok:
            data = r_search.json()
            position_matches = data.get('positionMatches', [])
            
            for category in position_matches:
                if category.get('trackName') == 'snp151' or category.get('name') == 'snp151':
                    for match in category.get('matches', []):
                        highlight_str = match.get('highlight', '')
                        m = re.search(r"chr([0-9XY]+):(\d+)-(\d+)", highlight_str, re.IGNORECASE)
                        if m:
                            c = m.group(1).upper()
                            start_p = int(m.group(2))
                            end_p = int(m.group(3))
                            return c, start_p, end_p, rs_clean

            results = data.get('results', [])
            for res in results:
                if res.get('trackName') == 'snp151' or res.get('name') == 'snp151':
                    highlight_str = res.get('highlight', '')
                    m = re.search(r"chr([0-9XY]+):(\d+)-(\d+)", highlight_str, re.IGNORECASE)
                    if m:
                        c = m.group(1).upper()
                        start_p = int(m.group(2))
                        end_p = int(m.group(3))
                        return c, start_p, end_p, rs_clean

    except Exception as e:
        st.warning(f"Error querying UCSC search API: {e}")

    return None, None, None, rs_clean


# ==============================================================================
# 4. SIDEBAR CONFIGURATION - USER PARAMETERS & OPTIONS
# ==============================================================================
st.sidebar.header("⚙️ Flanking & Target Parameters")

# Length of intronic/genomic sequence retrieved upstream and downstream of target
flank_size = st.sidebar.number_input(
    "Flanking Region Length (bp)", 
    min_value=100, max_value=2000, value=400, step=50,
    on_change=clear_results,
    help="Length of flanking (intronic/genomic) sequence retrieved upstream and downstream of the target region."
)

# Number of bases adjacent to target exon included INSIDE target brackets [ ]
target_padding = st.sidebar.number_input(
    "Target Padding (bp)", 
    min_value=0, max_value=200, value=40, step=5,
    on_change=clear_results,
    help="Number of bases directly adjacent to the target included INSIDE target brackets [ ]."
)

st.sidebar.header("🔬 Variant Filtering Options")

enable_snp_filtering = st.sidebar.checkbox(
    "Filter SNPs in sequence", 
    value=True, 
    on_change=clear_results
)

if enable_snp_filtering:
    snp_source = st.sidebar.radio(
        "Variant Database Source:",
        ["gnomAD v4", "dbSNP (151)"],
        index=0,
        on_change=clear_results,
        help="Choose between gnomAD v4 (>800k exomes/genomes) or the classic dbSNP151 track."
    )

    maf_option = st.sidebar.selectbox(
        "Minimum MAF Threshold (%)",
        options=["0.01%", "0.10%", "1.00%", "Manual input"],
        index=0,
        on_change=clear_results,
        help="Select a predefined Minor Allele Frequency threshold or enter a custom value."
    )

    if maf_option == "0.01%":
        maf_pct_val = 0.01
    elif maf_option == "0.10%":
        maf_pct_val = 0.10
    elif maf_option == "1.00%":
        maf_pct_val = 1.00
    else:  # Custom manual numeric entry
        maf_pct_val = st.sidebar.number_input(
            "Custom MAF Threshold (%)",
            min_value=0.0,
            max_value=50.0,
            value=0.0100,
            step=0.1000,
            format="%.4f",
            on_change=clear_results,
            help="Enter a custom MAF threshold (up to 4 decimal places)."
        )

    # Convert percentage to decimal ratio (e.g., 0.01% -> 0.0001)
    snp_maf_threshold = maf_pct_val / 100.0

    # Minimum total observed alleles required in population sample size
    min_allele_n_threshold = st.sidebar.number_input(
        "Minimum Allele Count Threshold (alleleNs / AN)", 
        min_value=10, max_value=1000000, value=1000, step=1000,
        on_change=clear_results,
        help="Minimum total observed alleles in database to filter artifacts from small cohorts."
    )
else:
    snp_source = "None"
    snp_maf_threshold = 0.0
    min_allele_n_threshold = 0

st.sidebar.header("✂️ Segmentation Options")

enable_segmentation = st.sidebar.checkbox(
    "Segment Target into Sub-amplicons", 
    value=True, 
    on_change=clear_results
)

if enable_segmentation:
    split_threshold = st.sidebar.number_input(
        "Maximum Target Length Before Split (bp)",
        min_value=300, max_value=2000, value=600, step=50,
        on_change=clear_results,
        help="Target regions exceeding this length (including padding) will be automatically segmented."
    )

    min_overlap = st.sidebar.number_input(
        "Minimum Sub-Amplicon Overlap (bp)",
        min_value=10, max_value=200, value=50, step=10,
        on_change=clear_results,
        help="Required overlapping base pairs between consecutive sub-amplicon fragments."
    )
else:
    split_threshold = 999999
    min_overlap = 0

# External tool navigation links rendered as sidebar buttons (updated with width="stretch")
st.sidebar.header("🔗 External Tools")
st.sidebar.link_button("🌐 Open Primer3web", "https://primer3.ut.ee/", width="stretch")
st.sidebar.link_button("🌐 Open Primer3Plus", "https://www.primer3plus.com/", width="stretch")
st.sidebar.link_button("🌐 Open NCBI Primer-BLAST", "https://www.ncbi.nlm.nih.gov/tools/primer-blast/", width="stretch")
st.sidebar.link_button("🌐 Open UNAFold", "https://www.unafold.org/mfold/applications/dna-folding-form.php", width="stretch")
st.sidebar.link_button("🌐 Open UCSC In-Silico PCR", "https://genome.ucsc.edu/cgi-bin/hgPcr", width="stretch")


# ==============================================================================
# 5. INPUT METHOD SELECTION & PARSING VALIDATION
# ==============================================================================
input_type = st.radio(
    "Select Input Method:", 
    ["Genomic Coordinates (GRCh38/hg38)", "Gene Symbol + Exon Number", "rs-Number (dbSNP ID)"], 
    index=0,
    horizontal=True,
    on_change=clear_results
)

chrom = None
exon_start = None
exon_end = None
header_label = ""
input_validation_error = None

# OPTION 1: Direct genomic coordinate parsing
if input_type == "Genomic Coordinates (GRCh38/hg38)":
    coord_raw = st.text_input(
        "Enter Genomic Region (GRCh38/hg38):", 
        value="chr7:117559464-117559655", 
        placeholder="e.g., chr7:117559464-117559655 or chr7:117559464-117559464",
        on_change=clear_results,
        help="Accepts formats like chr4:158679508-158679699 or single nucleotides like chr7:117559464-117559464. Max 20,000 bp."
    ).strip()

    if coord_raw:
        # Regex capturing chromosome string and start/end integer strings
        m = re.match(r"^(?:chr)?([0-9XY]+)\s*:\s*([\d,\.]+)\s*-\s*([\d,\.]+)$", coord_raw, re.IGNORECASE)
        if m:
            chrom_candidate = m.group(1).upper()
            start_str = re.sub(r"[,\.]", "", m.group(2))
            end_str = re.sub(r"[,\.]", "", m.group(3))
            
            try:
                start_val = int(start_str)
                end_val = int(end_str)
                
                if end_val < start_val:
                    input_validation_error = f"Invalid coordinates: End position ({end_val:,}) cannot be smaller than Start position ({start_val:,})."
                else:
                    region_len = end_val - start_val + 1
                    if region_len > MAX_COORDINATE_REGION_SIZE:
                        input_validation_error = f"Region length ({region_len:,} bp) exceeds the maximum allowed limit of {MAX_COORDINATE_REGION_SIZE:,} bp."
                    else:
                        chrom = chrom_candidate
                        exon_start = start_val
                        exon_end = end_val
                        header_label = f"chr{chrom}:{exon_start}-{exon_end} (GRCh38/hg38)"
                        st.info(f"**Target Region:** `chr{chrom}:{exon_start:,}-{exon_end:,}` | **Length:** {region_len:,} bp")
            except ValueError:
                input_validation_error = "Could not parse start or end positions as valid integers."
        else:
            input_validation_error = "Invalid coordinate format. Please use the format 'chr4:158679508-158679699' or 'chr7:117559464-117559464'."
    else:
        st.caption("💡 *Example format: `chr7:117559464-117559655` or single base `chr7:117559464-117559464`*")

# OPTION 2: Gene symbol + exon lookup
elif input_type == "Gene Symbol + Exon Number":
    col1, col2 = st.columns(2)
    with col1:
        gene_input = st.text_input("Gene Symbol (e.g., CFTR, COMT, ANO5):", value="CFTR", placeholder="e.g., CFTR", on_change=clear_results).strip().upper()
    with col2:
        exon_num = st.number_input("Exon Number:", min_value=1, value=11, step=1, on_change=clear_results)
        
    if gene_input:
        try:
            gene_transcripts = fetch_refseq_transcripts(gene_input)
            
            if gene_transcripts:
                transcript_options = {f"{t['name']} ({t.get('exonCount', 0)} exons)": t for t in gene_transcripts}
                selected_label = st.selectbox("Select RefSeq Transcript (GRCh38/hg38):", list(transcript_options.keys()), on_change=clear_results)
                selected_t = transcript_options[selected_label]
                
                transcript_acc = selected_t['name']
                chrom = selected_t['chrom'].replace('chr', '')
                strand = selected_t['strand']
                
                # Convert comma-delimited UCSC exon start/end arrays to 1-based tuples
                exon_starts = [int(x) + 1 for x in selected_t['exonStarts'].strip(',').split(',') if x]
                exon_ends = [int(x) for x in selected_t['exonEnds'].strip(',').split(',') if x]
                exons = list(zip(exon_starts, exon_ends))
                
                # Sort exon order according to transcription strand direction (5' -> 3')
                if strand == '-':
                    exons = sorted(exons, key=lambda x: x[0], reverse=True)
                else:
                    exons = sorted(exons, key=lambda x: x[0], reverse=False)
                    
                if exon_num <= len(exons):
                    exon_start, exon_end = exons[exon_num - 1]
                    header_label = f"{gene_input} ({transcript_acc}) Exon {exon_num}"
                    st.info(
                        f"**Exon {exon_num}/{len(exons)} Coordinates (GRCh38/hg38):** `chr{chrom}:{exon_start}-{exon_end}` | "
                        f"**Exon Length:** {exon_end - exon_start + 1} bp"
                    )
                else:
                    input_validation_error = f"Exon {exon_num} does not exist in selected transcript (Total exons: {len(exons)})."
            else:
                input_validation_error = f"Could not find RefSeq transcripts for gene symbol **{gene_input}** on GRCh38/hg38."
        except Exception as e:
            input_validation_error = f"Error fetching gene transcripts: {e}"

# OPTION 3: dbSNP rs-number lookup
else:
    rs_input = st.text_input("Enter dbSNP rs-number (e.g., rs113993960, rs4680):", value="rs113993960", placeholder="e.g., rs113993960", on_change=clear_results).strip()
    if rs_input:
        c, s, e, formatted_rs = resolve_rsid(rs_input)
        if c and s and e:
            chrom = c
            exon_start = s
            exon_end = e
            header_label = formatted_rs
            st.info(
                f"**Resolved {formatted_rs} Position:** `chr{chrom}:{exon_start}-{exon_end}` (GRCh38/hg38)"
            )
        else:
            input_validation_error = f"Could not locate **{formatted_rs}** on GRCh38/hg38 genome assembly. Please check the rsID."

# Display error box if input coordinate parsing or lookup fails
if input_validation_error:
    st.error(f"❌ **Input Error:** {input_validation_error}")


# ==============================================================================
# 6. BUTTON & RATE-LIMIT COOLDOWN TIMER FRAGMENT
# ==============================================================================
# @st.fragment(run_every="1s") creates a localized UI block that ticks every 1 second,
# enabling a real-time countdown timer without triggering a full page reload.
@st.fragment(run_every="1s")
def render_button_and_cooldown():
    time_since_last_click = time.time() - st.session_state["last_button_click_time"]
    remaining_time = COOLDOWN_SECONDS - time_since_last_click
    
    col_btn, col_timer = st.columns([1, 3], vertical_alignment="center")
    
    with col_btn:
        btn_disabled = (remaining_time > 0) or (input_validation_error is not None)
        if st.button("🚀 Generate Sequence Templates", type="primary", disabled=btn_disabled):
            st.session_state["last_button_click_time"] = time.time()
            st.session_state["trigger_generate"] = True
            st.rerun()

    with col_timer:
        if remaining_time > 0:
            rem_sec = int(math.ceil(remaining_time))
            st.caption(f"⏳ Cooldown active: Please wait **{rem_sec} second(s)** before making another request.")

render_button_and_cooldown()


# ==============================================================================
# 7. MAIN SEQUENCE PROCESSING & TEMPLATE GENERATION PIPELINE
# ==============================================================================
if st.session_state.get("trigger_generate", False):
    st.session_state["trigger_generate"] = False
    
    if chrom and exon_start and exon_end and not input_validation_error:
        with st.spinner(f"Fetching genomic data via {snp_source if enable_snp_filtering else 'GRCh38'}..."):
            try:
                # Step 1: Calculate full target boundaries including target padding
                full_target_start = exon_start - target_padding
                full_target_end = exon_end + target_padding
                full_target_len = full_target_end - full_target_start + 1

                # Step 2: Determine unified coordinate boundaries covering full target AND flanking regions
                overall_fetch_start = full_target_start - flank_size
                overall_fetch_end = full_target_end + flank_size

                c_clean = str(chrom).lower().replace('chr', '')
                if c_clean == '23': c_clean = 'X'
                elif c_clean == '24': c_clean = 'Y'

                # Step 3: Execute EXACTLY ONE API query to fetch the entire reference sequence
                full_raw_seq = fetch_ucsc_sequence(c_clean, overall_fetch_start, overall_fetch_end)

                # Step 4: Execute EXACTLY ONE API query to fetch all SNPs across the unified locus
                global_snp_map = {}
                global_all_snps_dict = {}

                if enable_snp_filtering and full_raw_seq:
                    if "gnomAD" in snp_source:
                        raw_snps = fetch_gnomad_v4_snps(c_clean, overall_fetch_start, overall_fetch_end)
                        for snp in raw_snps:
                            total_n = snp['total_n']
                            maf = snp['maf']
                            rs_name = snp['name']
                            var_id = snp['variant_id']
                            var_pos = snp['pos']
                            ref_allele = snp.get('ref', '')
                            alt_allele = snp.get('alt', '')
                            ref_len = max(1, len(ref_allele))
                            ref_alt_str = f"{ref_allele} / {alt_allele}"

                            # Apply user-defined quality thresholds (AN sample size and MAF ratio)
                            if total_n >= min_allele_n_threshold and maf >= snp_maf_threshold:
                                gnomad_url = f"https://gnomad.broadinstitute.org/variant/{var_id}?dataset=gnomad_r4"
                                dbsnp_url = f"https://www.ncbi.nlm.nih.gov/snp/{rs_name}" if str(rs_name).startswith('rs') else "#"
                                
                                global_all_snps_dict[var_id] = {
                                    "RS / Variant ID": rs_name,
                                    "Position (GRCh38)": f"chr{chrom}:{var_pos}",
                                    "Ref / Alt": ref_alt_str,
                                    "Joint AF (%)": f"{snp['joint_af'] * 100:.4f}%",
                                    "Allele Count (AN)": f"{int(total_n):,}",
                                    "Source": "gnomAD v4",
                                    "gnomAD v4 Link": gnomad_url,
                                    "dbSNP Link": dbsnp_url,
                                    "raw_pos": var_pos,
                                    "var_id": var_id
                                }

                                # If input was rs-Number, dynamically synchronize uppercase exon coordinates to exact target SNP boundaries
                                if "rs-Number" in input_type and rs_name.lower() in header_label.lower():
                                    exon_start = var_pos
                                    exon_end = var_pos + ref_len - 1

                                # Map variant metadata to every base offset spanned by reference allele
                                for bp_offset in range(ref_len):
                                    curr_pos = var_pos + bp_offset
                                    rel_pos = curr_pos - overall_fetch_start
                                    if 0 <= rel_pos < len(full_raw_seq):
                                        global_snp_map[rel_pos] = {
                                            'name': rs_name,
                                            'variant_id': var_id,
                                            'maf': maf * 100,
                                            'total_n': int(total_n),
                                            'pos': var_pos,
                                            'ref_alt': ref_alt_str,
                                            'source': 'gnomAD v4'
                                        }

                    else:  # dbSNP (151) track mapping
                        raw_snps = fetch_ucsc_snps(c_clean, overall_fetch_start, overall_fetch_end)
                        for snp in raw_snps:
                            allele_ns_raw = snp.get('alleleNs', '')
                            allele_freqs_raw = snp.get('alleleFreqs', '')
                            rs_name = snp.get('name', 'rsID Unknown')
                            
                            ref_allele = snp.get('ref', '')
                            observed_str = snp.get('observed', '') or snp.get('alleles', '')
                            
                            if not ref_allele and observed_str and '/' in observed_str:
                                parts = observed_str.split('/')
                                ref_allele = parts[0]
                                alt_str = "/".join(parts[1:])
                            elif observed_str:
                                alt_str = observed_str
                            else:
                                alt_str = "-"
                            
                            ref_alt_display = f"{ref_allele if ref_allele else '-'} / {alt_str}"
                            
                            total_n = 0.0
                            maf = 0.0
                            
                            if allele_ns_raw:
                                try:
                                    ns = [float(n) for n in str(allele_ns_raw).strip(',').split(',') if n]
                                    total_n = sum(ns)
                                except ValueError:
                                    total_n = 0.0
                                    
                            if allele_freqs_raw:
                                try:
                                    freqs = [float(f) for f in str(allele_freqs_raw).strip(',').split(',') if f]
                                    if len(freqs) > 1:
                                        sorted_freqs = sorted(freqs)
                                        maf = sorted_freqs[-2]
                                except ValueError:
                                    maf = 0.0
                            
                            if total_n >= min_allele_n_threshold and maf >= snp_maf_threshold:
                                rel_s = max(0, snp['chromStart'] - (overall_fetch_start - 1))
                                rel_e = min(len(full_raw_seq), snp['chromEnd'] - (overall_fetch_start - 1))
                                var_pos = overall_fetch_start + rel_s
                                
                                dbsnp_url = f"https://www.ncbi.nlm.nih.gov/snp/{rs_name}" if str(rs_name).startswith('rs') else "#"
                                gnomad_url = f"https://gnomad.broadinstitute.org/variant/{rs_name}?dataset=gnomad_r4" if str(rs_name).startswith('rs') else "#"
                                    
                                global_all_snps_dict[rs_name] = {
                                    "RS / Variant ID": rs_name,
                                    "Position (GRCh38)": f"chr{chrom}:{var_pos}",
                                    "Ref / Alt": ref_alt_display,
                                    "MAF (%)": f"{maf * 100:.4f}%",
                                    "Allele Count (AN)": f"{int(total_n):,}",
                                    "Source": "dbSNP (151)",
                                    "dbSNP Link": dbsnp_url,
                                    "gnomAD v4 Link": gnomad_url,
                                    "raw_pos": var_pos,
                                    "var_id": rs_name
                                }

                                # If input was rs-Number, dynamically synchronize uppercase exon coordinates to exact target SNP boundaries
                                if "rs-Number" in input_type and rs_name.lower() in header_label.lower():
                                    exon_start = snp['chromStart'] + 1
                                    exon_end = snp['chromEnd']

                                snp_info = {
                                    'name': rs_name,
                                    'variant_id': rs_name,
                                    'maf': maf * 100,
                                    'total_n': int(total_n),
                                    'pos': var_pos,
                                    'ref_alt': ref_alt_display,
                                    'source': 'dbSNP (151)'
                                }
                                
                                for pos_idx in range(rel_s, max(rel_s + 1, rel_e)):
                                    if 0 <= pos_idx < len(full_raw_seq):
                                        global_snp_map[pos_idx] = snp_info

                # Step 5: Perform target segmentation on the pre-fetched genomic span
                sub_targets = []
                if enable_segmentation and full_target_len > split_threshold:
                    num_fragments = math.ceil((full_target_len - min_overlap) / (split_threshold - min_overlap))
                    effective_step = math.ceil((full_target_len - min_overlap) / num_fragments)
                    
                    curr_start = full_target_start
                    for idx in range(num_fragments):
                        curr_end = min(curr_start + effective_step + min_overlap - 1, full_target_end)
                        if idx == num_fragments - 1:
                            curr_end = full_target_end
                        sub_targets.append((curr_start, curr_end))
                        curr_start = curr_end - min_overlap + 1
                else:
                    sub_targets.append((full_target_start, full_target_end))

                # Step 6: Assign sequence slices and mapped SNPs to each sub-amplicon fragment
                processed_fragments = []

                for sub_idx, (t_start, t_end) in enumerate(sub_targets):
                    frag_fetch_start = t_start - flank_size
                    frag_fetch_end = t_end + flank_size

                    # Slice raw sequence substring from pre-fetched full sequence
                    rel_slice_start = frag_fetch_start - overall_fetch_start
                    rel_slice_end = frag_fetch_end - overall_fetch_start + 1
                    raw_seq = full_raw_seq[rel_slice_start:rel_slice_end]

                    if raw_seq:
                        snp_map = {}
                        detailed_snps_list = []
                        unique_snps_count = set()

                        # Collect global SNPs that fall inside this sub-amplicon's coordinate window
                        for g_var_id, g_snp in global_all_snps_dict.items():
                            v_pos = g_snp["raw_pos"]
                            if frag_fetch_start <= v_pos <= frag_fetch_end:
                                unique_snps_count.add(g_var_id)
                                clean_dict = {k: v for k, v in g_snp.items() if k not in ["raw_pos", "var_id"]}
                                detailed_snps_list.append(clean_dict)

                        # Translate global SNP coordinate indices to local sub-amplicon base positions
                        for rel_p, snp_info in global_snp_map.items():
                            abs_p = overall_fetch_start + rel_p
                            if frag_fetch_start <= abs_p <= frag_fetch_end:
                                local_rel_p = abs_p - frag_fetch_start
                                if 0 <= local_rel_p < len(raw_seq):
                                    snp_map[local_rel_p] = snp_info

                        # Construct Primer3 template string (< > [ ]) and HTML visualization tags
                        target_rel_s = exon_start - frag_fetch_start
                        target_rel_e = exon_end - frag_fetch_start
                        
                        sub_target_rel_s = t_start - frag_fetch_start
                        sub_target_rel_e = t_end - frag_fetch_start

                        p3_str_list = []
                        html_list = []
                        in_p3_snp = False

                        for i in range(len(raw_seq)):
                            base = raw_seq[i]
                            
                            # Capitalize target bases; lowercase flanking regions
                            is_target_region = (target_rel_s <= i <= target_rel_e)
                            display_base = base if is_target_region else base.lower()

                            # Insert target start bracket '['
                            if i == sub_target_rel_s:
                                p3_str_list.append('[')
                                html_list.append('<b style="color: #1E88E5; font-size: 1.2em;">[</b>')

                            is_snp_base = i in snp_map
                            
                            # Insert Primer3 SNP masking brackets '<' and '>'
                            if enable_snp_filtering:
                                if is_snp_base and not in_p3_snp:
                                    p3_str_list.append('<')
                                    in_p3_snp = True
                                elif not is_snp_base and in_p3_snp:
                                    p3_str_list.append('>')
                                    in_p3_snp = False

                            p3_str_list.append(display_base)

                            # Build HTML visualizer markup with red underlining and hover tooltip (including Ref/Alt)
                            if is_snp_base and enable_snp_filtering:
                                snp_meta = snp_map[i]
                                rs_id = snp_meta['name']
                                var_id = snp_meta['variant_id']
                                ref_alt_info = snp_meta.get('ref_alt', '')

                                if "gnomAD" in snp_meta['source']:
                                    link_url = f"https://gnomad.broadinstitute.org/variant/{var_id}?dataset=gnomad_r4"
                                else:
                                    link_url = f"https://gnomad.broadinstitute.org/variant/{rs_id}?dataset=gnomad_r4" if str(rs_id).startswith('rs') else "#"

                                maf_text = f" | AF: {snp_meta['maf']:.4f}% (AN={snp_meta['total_n']:,})" if snp_meta['maf'] > 0 else ""
                                ref_alt_text = f" ({ref_alt_info})" if ref_alt_info else ""
                                tooltip_text = f"{rs_id}{ref_alt_text}{maf_text} | Pos: chr{chrom}:{snp_meta['pos']} ({snp_meta['source']})"
                                
                                html_list.append(
                                    f'<a href="{link_url}" target="_blank" title="{tooltip_text}" '
                                    f'style="color: #E53935; font-weight: bold; text-decoration: underline; '
                                    f'padding: 0; margin: 0; display: inline;">{display_base}</a>'
                                )
                            else:
                                html_list.append(display_base)

                            # Insert target end bracket ']'
                            if i == sub_target_rel_e:
                                p3_str_list.append(']')
                                html_list.append('<b style="color: #1E88E5; font-size: 1.2em;">]</b>')

                        if enable_snp_filtering and in_p3_snp:
                            p3_str_list.append('>')

                        primer3_text = "".join(p3_str_list)

                        processed_fragments.append({
                            "sub_idx": sub_idx,
                            "t_start": t_start,
                            "t_end": t_end,
                            "primer3_text": primer3_text,
                            "html_list": html_list,
                            "unique_snps_count": unique_snps_count,
                            "detailed_snps_list": detailed_snps_list,
                            "total_bp": t_end - t_start + 1
                        })

                # Save processed pipeline results to session state and trigger UI re-render
                st.session_state["results"] = {
                    "chrom": chrom,
                    "full_target_len": full_target_len,
                    "fragments": processed_fragments,
                    "input_type": input_type,
                    "header_label": header_label,
                    "enable_snp_filtering": enable_snp_filtering,
                    "snp_source": snp_source if enable_snp_filtering else "None",
                    "min_allele_n_threshold": min_allele_n_threshold,
                    "snp_maf_threshold_pct": snp_maf_threshold * 100
                }
                st.rerun()

            except Exception as err:
                st.session_state["results"] = None
                err_str = str(err)
                if "Read timed out" in err_str or "HTTPSConnectionPool" in err_str:
                    user_error_msg = f"The requested service ({snp_source if enable_snp_filtering else 'UCSC'}) timed out or is currently unresponsive. Please try again or select a different database in the sidebar."
                else:
                    user_error_msg = f"An error occurred while querying genomic data: {err}"

                st.error(f"❌ **Execution Error:** {user_error_msg}")
    else:
        st.warning("Please specify valid target coordinates, gene parameters, or dbSNP rs-number.")


# ==============================================================================
# 8. RESULTS DISPLAY & RENDERING SECTION
# ==============================================================================
if st.session_state["results"]:
    res = st.session_state["results"]
    fragments = res["fragments"]
    lbl = res["header_label"]
    
    st.success(
        f"Target Region: {res['full_target_len']} bp (including padding). "
        f"Divided into **{len(fragments)} sub-amplicon(s)**. Variant Source: **{res['snp_source']}**."
    )

    # Render a tab view for each generated sub-amplicon fragment
    tabs = st.tabs([f"Fragment {f['sub_idx']+1} ({f['total_bp']} bp)" for f in fragments])

    for f in fragments:
        sub_idx = f["sub_idx"]
        with tabs[sub_idx]:
            st.markdown(f"**Fragment {sub_idx+1} Coordinates (GRCh38/hg38):** `chr{res['chrom']}:{f['t_start']}-{f['t_end']}` ({f['total_bp']} bp)")
            if res["enable_snp_filtering"]:
                st.metric(
                    f"Identified Critical SNPs ({res['snp_source']} | MAF ≥ {res['snp_maf_threshold_pct']:.4f}% | AN ≥ {res['min_allele_n_threshold']:,})", 
                    len(f["unique_snps_count"])
                )

            st.subheader(f"👁️ Sequence Visualizer ({lbl})")
            
            if "Gene Symbol" in res["input_type"]:
                case_legend = "UPPERCASE = Exonic Sequence | lowercase = Flanking Sequence"
            elif "rs-Number" in res["input_type"]:
                case_legend = "UPPERCASE = Target SNP Base | lowercase = Flanking Sequence"
            else:
                case_legend = "UPPERCASE = Target Region | lowercase = Flanking Sequence"

            st.caption(
                f"**Legend:** Blue brackets **`[`** **`]`** = Target Region + Padding | {case_legend} | "
                " <span style='color: #E53935; font-weight: bold; text-decoration: underline;'>Red Underlined</span> = Flagged SNP. "
                f" *Hover over a flagged SNP for variant/AF details, or click to open {res['snp_source']}.*", 
                unsafe_allow_html=True
            )

            html_content = "".join(f["html_list"])

            # Render HTML sequence visualizer box
            st.markdown(
                f"""
                <div style="
                    font-family: 'Courier New', Courier, monospace; 
                    background-color: var(--background-color);
                    color: var(--text-color);
                    padding: 18px; 
                    border-radius: 8px; 
                    border: 1px solid rgba(128, 128, 128, 0.3); 
                    line-height: 1.6; 
                    font-size: 14px;
                    word-break: break-all;
                ">
                    {html_content}
                </div>
                """, 
                unsafe_allow_html=True
            )

            # Render expandable data frame containing detailed SNP metrics and external database links
            if res["enable_snp_filtering"] and f.get("detailed_snps_list"):
                with st.expander(f"🔽 Show Detailed List of Identified SNPs ({len(f['detailed_snps_list'])} variants found)", expanded=False):
                    df_snps = pd.DataFrame(f["detailed_snps_list"])
                    
                    # Shorten column labels to "gnomAD" and "dbSNP" with clean "Link" button text
                    col_config = {}
                    
                    if "gnomAD v4 Link" in df_snps.columns:
                        col_config["gnomAD v4 Link"] = st.column_config.LinkColumn(
                            "gnomAD",
                            display_text="Link",
                            width="small"
                        )
                    
                    if "dbSNP Link" in df_snps.columns:
                        col_config["dbSNP Link"] = st.column_config.LinkColumn(
                            "dbSNP",
                            display_text="Link",
                            width="small"
                        )

                    st.dataframe(
                        df_snps,
                        column_config=col_config,
                        width="stretch",
                        hide_index=True
                    )

            st.subheader(f"📋 Primer3 Copy-Paste Output ({lbl})")
            st.text_area(
                f"Primer3 Template Sequence (Fragment {sub_idx+1}):",
                value=f["primer3_text"],
                height=150,
                key=f"p3_text_{sub_idx}"
            )

            # Embedded iframe copy button to allow direct browser clipboard interaction
            col_copy, col_spacer = st.columns([2.5, 5])
            with col_copy:
                st.iframe(
                    f"""
                    <!DOCTYPE html>
                    <html>
                    <head>
                    <style>
                        body {{
                            margin: 0;
                            padding: 0;
                            background: transparent;
                        }}
                        button {{
                            background-color: #f0f2f6;
                            color: #31333F;
                            border: 1px solid rgba(49, 51, 63, 0.2);
                            padding: 6px 16px;
                            border-radius: 8px;
                            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
                            font-size: 14px;
                            font-weight: 500;
                            cursor: pointer;
                            width: 100%;
                            height: 38px;
                            white-space: nowrap;
                            display: flex;
                            align-items: center;
                            justify-content: center;
                            gap: 6px;
                            box-sizing: border-box;
                        }}
                    </style>
                    </head>
                    <body>
                    <button id="copy-btn-{sub_idx}" onclick="copySequence()">📋 Copy Sequence</button>
                    <script>
                    function copySequence() {{
                        var seqText = {repr(f["primer3_text"])};
                        navigator.clipboard.writeText(seqText).then(function() {{
                            var parentDoc = window.parent.document;
                            var textareas = parentDoc.querySelectorAll('textarea');
                            var targetArea = textareas[{sub_idx}];
                            if (targetArea) {{
                                var originalBg = targetArea.style.backgroundColor;
                                var originalTransition = targetArea.style.transition;
                                targetArea.style.transition = 'background-color 0.2s ease';
                                targetArea.style.backgroundColor = '#d4edda';
                                setTimeout(function() {{
                                    targetArea.style.backgroundColor = originalBg;
                                    setTimeout(function() {{
                                        targetArea.style.transition = originalTransition;
                                    }}, 300);
                                }}, 800);
                            }}
                        }});
                    }}
                    </script>
                    </body>
                    </html>
                    """,
                    height=40
                )
