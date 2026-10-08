import math
import re
import time
import requests
import streamlit as st

# ==========================================
# PAGE CONFIGURATION
# ==========================================
st.set_page_config(
    page_title="SNIPER3 - Primer3 Template Generator",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 SNIPER3")
st.caption("**S**NP-**N**eutralizing **I**ntelligent **P**rimer3 **E**xon & **R**egion Template Generator")

st.markdown(
    "Generate standardized sequence templates formatted for Primer3Web and Primer3Plus using **GRCh38/hg38** coordinates. "
    "Features include large exon auto-splitting, exonic/target sequence capitalization, "
    "and automated SNP filtering (< >) by MAF and population sample size."
)

# Initialize Session State
if "results" not in st.session_state:
    st.session_state["results"] = None

if "last_button_click_time" not in st.session_state:
    st.session_state["last_button_click_time"] = 0.0

if "trigger_generate" not in st.session_state:
    st.session_state["trigger_generate"] = False

def clear_results():
    """Clears generated template outputs while keeping input form selections intact."""
    st.session_state["results"] = None

# ==========================================
# CACHED API FUNCTIONS
# ==========================================
@st.cache_data(ttl=3600, show_spinner=False)
def fetch_refseq_transcripts(gene_symbol: str):
    """Fetch and cache RefSeq transcript items for a gene symbol via UCSC REST API (GRCh38/hg38)."""
    gene_clean = gene_symbol.strip().upper()
    search_url = f"https://api.genome.ucsc.edu/search?genome=hg38;search={gene_clean}"
    r_search = requests.get(search_url, timeout=10)
    
    if r_search.ok:
        search_data = r_search.json()
        pos_str = None
        
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

        if pos_str:
            match = re.match(r"(chr[0-9XY]+):(\d+)-(\d+)", pos_str, re.IGNORECASE)
            if match:
                c_name, s_pos, e_pos = match.group(1), int(match.group(2)), int(match.group(3))
                track_url = f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=ncbiRefSeq;chrom={c_name};start={s_pos};end={e_pos}"
                r_track = requests.get(track_url, timeout=10)
                if r_track.ok:
                    refseq_items = r_track.json().get('ncbiRefSeq', [])
                    gene_transcripts = [t for t in refseq_items if t.get('name2', '').upper() == gene_clean and t.get('name', '').startswith('NM_')]
                    if not gene_transcripts:
                        gene_transcripts = [t for t in refseq_items if t.get('name2', '').upper() == gene_clean]
                    return sorted(gene_transcripts, key=lambda x: x.get('exonCount', 0), reverse=True)
    return []

@st.cache_data(ttl=3600, show_spinner=False)
def fetch_ucsc_sequence(chrom: str, start_pos: int, end_pos: int):
    """Fetch genomic DNA sequence from UCSC REST API (GRCh38/hg38)."""
    c_clean = str(chrom).lower().replace('chr', '')
    if c_clean == '23': c_clean = 'X'
    elif c_clean == '24': c_clean = 'Y'
    
    ucsc_seq_url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr{c_clean};start={start_pos-1};end={end_pos}"
    r_seq = requests.get(ucsc_seq_url, timeout=10)
    if r_seq.ok and 'dna' in r_seq.json():
        return r_seq.json()['dna'].upper()
    return None

@st.cache_data(ttl=3600, show_spinner=False)
def fetch_ucsc_snps(chrom: str, start_pos: int, end_pos: int):
    """Fetch snp151 track data from UCSC REST API (GRCh38/hg38)."""
    c_clean = str(chrom).lower().replace('chr', '')
    if c_clean == '23': c_clean = 'X'
    elif c_clean == '24': c_clean = 'Y'
    
    snp_url = f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=snp151;chrom=chr{c_clean};start={start_pos-1};end={end_pos}"
    r_snp = requests.get(snp_url, timeout=10)
    if r_snp.ok:
        return r_snp.json().get('snp151', [])
    return []

@st.cache_data(ttl=3600, show_spinner=False)
def resolve_rsid(rs_id: str):
    """
    Resolves an rsID to exact single-nucleotide GRCh38/hg38 coordinates using the UCSC search API
    by parsing the 'highlight' field inside the 'snp151' track match.
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

# ==========================================
# SIDEBAR CONFIGURATION
# ==========================================
st.sidebar.header("⚙️ Flanking & Target Parameters")

flank_size = st.sidebar.number_input(
    "Flanking Region Length (bp)", 
    min_value=100, max_value=2000, value=400, step=50,
    on_change=clear_results,
    help="Length of intronic/genomic flanking sequence retrieved upstream and downstream of the target region."
)

target_padding = st.sidebar.number_input(
    "Target Padding (bp)", 
    min_value=0, max_value=200, value=40, step=5,
    on_change=clear_results,
    help="Number of base pairs immediately adjacent to the target included inside the target brackets [ ]."
)

st.sidebar.header("🔬 Variant Filtering Options")

enable_snp_filtering = st.sidebar.checkbox(
    "Highlight SNPs in Sequence", 
    value=True, 
    on_change=clear_results
)

if enable_snp_filtering:
    snp_maf_threshold = st.sidebar.number_input(
        "Minimum MAF Threshold (%)", 
        min_value=0.0, max_value=50.0, value=0.1, step=0.1,
        on_change=clear_results,
        help="Variants with a Minor Allele Frequency (MAF) at or above this threshold will be flagged."
    ) / 100.0

    min_allele_n_threshold = st.sidebar.number_input(
        "Minimum Allele Count Threshold (alleleNs)", 
        min_value=10, max_value=100000, value=1000, step=100,
        on_change=clear_results,
        help="Minimum total observed alleles in dbSNP. Filters out small-cohort artifacts (e.g., SGDP)."
    )
else:
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
        help="Target regions exceeding this length (including padding) will be split into overlapping sub-amplicons."
    )

    min_overlap = st.sidebar.number_input(
        "Minimum Sub-Amplicon Overlap (bp)",
        min_value=10, max_value=200, value=50, step=10,
        on_change=clear_results,
        help="Required base pair overlap between consecutive sub-amplicon fragments."
    )
else:
    split_threshold = 999999
    min_overlap = 0

st.sidebar.header("🔗 External Tools")
st.sidebar.link_button("🌐 Open Primer3web", "https://primer3.ut.ee/", use_container_width=True)
st.sidebar.link_button("🌐 Open Primer3Plus", "https://www.primer3plus.com/", use_container_width=True)
st.sidebar.link_button("🌐 Open NCBI Primer-BLAST", "https://www.ncbi.nlm.nih.gov/tools/primer-blast/", use_container_width=True)
st.sidebar.link_button("🌐 Open UNAFold", "https://www.unafold.org/mfold/applications/dna-folding-form.php", use_container_width=True)
st.sidebar.link_button("🌐 Open UCSC In-Silico PCR", "https://genome.ucsc.edu/cgi-bin/hgPcr", use_container_width=True)

# ==========================================
# INPUT SELECTION
# ==========================================
input_type = st.radio(
    "Select Input Method:", 
    ["Gene Symbol + Exon Number", "Genomic Coordinates (GRCh38/hg38)", "rs-Number (dbSNP ID)"], 
    horizontal=True,
    on_change=clear_results
)

chrom = None
exon_start = None
exon_end = None
header_label = ""

if input_type == "Gene Symbol + Exon Number":
    col1, col2 = st.columns(2)
    with col1:
        gene_input = st.text_input("Gene Symbol (e.g., COMT, FKRP, ANO5):", value="COMT", on_change=clear_results).strip().upper()
    with col2:
        exon_num = st.number_input("Exon Number:", min_value=1, value=4, step=1, on_change=clear_results)
        
    if gene_input:
        gene_transcripts = fetch_refseq_transcripts(gene_input)
        
        if gene_transcripts:
            transcript_options = {f"{t['name']} ({t.get('exonCount', 0)} exons)": t for t in gene_transcripts}
            selected_label = st.selectbox("Select RefSeq Transcript (GRCh38/hg38):", list(transcript_options.keys()), on_change=clear_results)
            selected_t = transcript_options[selected_label]
            
            transcript_acc = selected_t['name']
            chrom = selected_t['chrom'].replace('chr', '')
            strand = selected_t['strand']
            
            exon_starts = [int(x) + 1 for x in selected_t['exonStarts'].strip(',').split(',') if x]
            exon_ends = [int(x) for x in selected_t['exonEnds'].strip(',').split(',') if x]
            exons = list(zip(exon_starts, exon_ends))
            
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
                st.error(f"Exon {exon_num} does not exist in selected transcript (Total exons: {len(exons)}).")
        else:
            st.error(f"Could not find RefSeq transcripts for gene symbol **{gene_input}** on GRCh38/hg38.")

elif input_type == "Genomic Coordinates (GRCh38/hg38)":
    col1, col2, col3 = st.columns(3)
    with col1:
        chrom = st.text_input("Chromosome (e.g., 22, X):", value="22", on_change=clear_results).strip()
    with col2:
        exon_start = st.number_input("Target Start Position (GRCh38/hg38):", min_value=1, value=19963748, on_change=clear_results)
    with col3:
        exon_end = st.number_input("Target End Position (GRCh38/hg38):", min_value=1, value=19963748, on_change=clear_results)
    
    if chrom and exon_start and exon_end:
        header_label = f"chr{chrom}:{exon_start}-{exon_end} (GRCh38/hg38)"

else:  # rs-Number input
    rs_input = st.text_input("Enter dbSNP rs-number (e.g., rs4680, rs113744932):", value="rs4680", on_change=clear_results).strip()
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
            st.error(f"Could not locate **{formatted_rs}** on GRCh38/hg38 genome assembly. Please check the rsID.")

# ==========================================
# BUTTON & LIVE TIMER FRAGMENT
# ==========================================
COOLDOWN_SECONDS = 10.0

@st.fragment(run_every="1s")
def render_button_and_cooldown():
    time_since_last_click = time.time() - st.session_state["last_button_click_time"]
    remaining_time = COOLDOWN_SECONDS - time_since_last_click
    
    col_btn, col_timer = st.columns([1, 3], vertical_alignment="center")
    
    with col_btn:
        btn_disabled = remaining_time > 0
        if st.button("🚀 Generate Sequence Templates", type="primary", disabled=btn_disabled):
            st.session_state["last_button_click_time"] = time.time()
            st.session_state["trigger_generate"] = True
            st.rerun()

    with col_timer:
        if remaining_time > 0:
            rem_sec = int(math.ceil(remaining_time))
            st.caption(f"⏳ Cooldown active: Please wait **{rem_sec} second(s)** before making another request.")

render_button_and_cooldown()

# ==========================================
# SEQUENCE PROCESSING & EXECUTION
# ==========================================
if st.session_state.get("trigger_generate", False):
    st.session_state["trigger_generate"] = False
    
    if chrom and exon_start and exon_end:
        with st.spinner("Processing genomic sequence and variant data (GRCh38/hg38)..."):
            
            full_target_start = exon_start - target_padding
            full_target_end = exon_end + target_padding
            full_target_len = full_target_end - full_target_start + 1
            
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

            processed_fragments = []

            for sub_idx, (t_start, t_end) in enumerate(sub_targets):
                fetch_start = t_start - flank_size
                fetch_end = t_end + flank_size
                
                c_clean = str(chrom).lower().replace('chr', '')
                if c_clean == '23': c_clean = 'X'
                elif c_clean == '24': c_clean = 'Y'
                
                raw_seq = fetch_ucsc_sequence(c_clean, fetch_start, fetch_end)
                
                if raw_seq:
                    snp_map = {}
                    unique_snps_count = set()
                    
                    if enable_snp_filtering:
                        raw_snps = fetch_ucsc_snps(c_clean, fetch_start, fetch_end)
                        
                        for snp in raw_snps:
                            allele_ns_raw = snp.get('alleleNs', '')
                            allele_freqs_raw = snp.get('alleleFreqs', '')
                            rs_name = snp.get('name', 'rsID Unknown')
                            
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
                                rel_s = max(0, snp['chromStart'] - (fetch_start - 1))
                                rel_e = min(len(raw_seq), snp['chromEnd'] - (fetch_start - 1))
                                
                                unique_snps_count.add(rs_name)
                                
                                snp_info = {
                                    'name': rs_name,
                                    'maf': maf * 100,
                                    'total_n': int(total_n),
                                    'pos': fetch_start + rel_s
                                }
                                
                                for pos_idx in range(rel_s, rel_e):
                                    snp_map[pos_idx] = snp_info

                    target_rel_s = exon_start - fetch_start
                    target_rel_e = exon_end - fetch_start
                    
                    sub_target_rel_s = t_start - fetch_start
                    sub_target_rel_e = t_end - fetch_start

                    p3_str_list = []
                    html_list = []
                    in_p3_snp = False

                    for i in range(len(raw_seq)):
                        base = raw_seq[i]
                        
                        is_target_region = (target_rel_s <= i <= target_rel_e)
                        display_base = base if is_target_region else base.lower()

                        if i == sub_target_rel_s:
                            p3_str_list.append('[')
                            html_list.append('<b style="color: #1E88E5; font-size: 1.2em;">[</b>')

                        is_snp_base = i in snp_map
                        
                        if enable_snp_filtering:
                            if is_snp_base and not in_p3_snp:
                                p3_str_list.append('<')
                                in_p3_snp = True
                            elif not is_snp_base and in_p3_snp:
                                p3_str_list.append('>')
                                in_p3_snp = False

                        p3_str_list.append(display_base)

                        if is_snp_base and enable_snp_filtering:
                            snp_meta = snp_map[i]
                            rs_id = snp_meta['name']
                            ncbi_url = f"https://www.ncbi.nlm.nih.gov/snp/{rs_id}" if str(rs_id).startswith('rs') else "#"
                            maf_text = f" | MAF: {snp_meta['maf']:.2f}% (N={snp_meta['total_n']})" if snp_meta['maf'] > 0 else ""
                            tooltip_text = f"{rs_id}{maf_text} | Position: chr{chrom}:{snp_meta['pos']} (GRCh38/hg38)"
                            
                            html_list.append(
                                f'<a href="{ncbi_url}" target="_blank" title="{tooltip_text}" '
                                f'style="color: #E53935; font-weight: bold; text-decoration: underline; '
                                f'padding: 0; margin: 0; display: inline;">{display_base}</a>'
                            )
                        else:
                            html_list.append(display_base)

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
                        "total_bp": t_end - t_start + 1
                    })

            st.session_state["results"] = {
                "chrom": chrom,
                "full_target_len": full_target_len,
                "fragments": processed_fragments,
                "input_type": input_type,
                "header_label": header_label,
                "enable_snp_filtering": enable_snp_filtering,
                "min_allele_n_threshold": min_allele_n_threshold
            }
            st.rerun()
    else:
        st.warning("Please specify valid target coordinates, gene parameters, or dbSNP rs-number.")

# ==========================================
# DISPLAY PERSISTED RESULTS
# ==========================================
if st.session_state["results"]:
    res = st.session_state["results"]
    fragments = res["fragments"]
    lbl = res["header_label"]
    
    st.success(
        f"Target Region: {res['full_target_len']} bp (including padding). "
        f"Divided into **{len(fragments)} sub-amplicon(s)**."
    )

    tabs = st.tabs([f"Fragment {f['sub_idx']+1} ({f['total_bp']} bp)" for f in fragments])

    for f in fragments:
        sub_idx = f["sub_idx"]
        with tabs[sub_idx]:
            st.markdown(f"**Fragment {sub_idx+1} Coordinates (GRCh38/hg38):** `chr{res['chrom']}:{f['t_start']}-{f['t_end']}` ({f['total_bp']} bp)")
            if res["enable_snp_filtering"]:
                st.metric("Identified Critical SNPs (N ≥ " + f"{res['min_allele_n_threshold']})", len(f["unique_snps_count"]))

            st.subheader(f"👁️ Sequence Visualizer ({lbl})")
            
            if "Gene Symbol" in res["input_type"]:
                case_legend = "UPPERCASE = Exonic Sequence | lowercase = Intronic Sequence"
            elif "rs-Number" in res["input_type"]:
                case_legend = "UPPERCASE = Target SNP Base | lowercase = Flanking Sequence"
            else:
                case_legend = "UPPERCASE = Target Region | lowercase = Flanking Sequence"

            st.caption(
                f"**Legend:** Blue brackets **`[`** **`]`** = Target Region | {case_legend} | "
                " <span style='color: #E53935; font-weight: bold; text-decoration: underline;'>Red Underlined</span> = Flagged SNP. "
                " *Hover over a flagged SNP for rsID/MAF details, or click to open NCBI dbSNP.*", 
                unsafe_allow_html=True
            )

            html_content = "".join(f["html_list"])

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

            st.subheader(f"📋 Primer3 Copy-Paste Output ({lbl})")
            st.text_area(
                f"Primer3 Template Sequence (Fragment {sub_idx+1}):",
                value=f["primer3_text"],
                height=150,
                key=f"p3_text_{sub_idx}"
            )

            st.download_button(
                label=f"💾 Download Fragment {sub_idx+1} (.txt)",
                data=f"SEQUENCE_ID=chr{res['chrom']}:{f['t_start']}-{f['t_end']}_frag{sub_idx+1}_hg38\nSEQUENCE_TEMPLATE={f['primer3_text']}\n",
                file_name=f"sniper3_fragment_{sub_idx+1}_chr{res['chrom']}_{f['t_start']}_{f['t_end']}_hg38.txt",
                mime="text/plain",
                key=f"dl_btn_{sub_idx}"
            )