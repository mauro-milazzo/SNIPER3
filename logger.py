import datetime
import streamlit as st
import gspread


def log_search_query(
    status: str,
    input_type: str,
    query_label: str,
    chrom: str = None,
    target_start: int = None,
    target_end: int = None,
    flank_size: int = None,
    target_padding: int = None,
    enable_snp_filtering: bool = None,
    snp_source: str = None,
    maf_threshold_pct: float = None,
    min_allele_n: int = None,
    enable_segmentation: bool = None,
    split_threshold: int = None,
    min_overlap: int = None,
    error_msg: str = None
):
    """
    Logs search parameters, options, and execution status/errors to Google Sheets.
    Executes ONLY when 'IS_PRODUCTION' is set to True in Streamlit Secrets.
    """
    # 1. Environment check: Run only when IS_PRODUCTION is enabled
    try:
        is_production = st.secrets.get("IS_PRODUCTION", False)
    except Exception:
        is_production = False

    if not is_production:
        return

    try:
        # 2. Retrieve credentials and parse private key formatting
        creds_dict = dict(st.secrets["gcp_service_account"])
        if "private_key" in creds_dict:
            creds_dict["private_key"] = creds_dict["private_key"].replace("\\n", "\n")

        # 3. Authenticate with Google Sheets
        gc = gspread.service_account_from_dict(creds_dict)
        sheet_id = st.secrets.get("GOOGLE_SHEET_ID")
        sheet = gc.open_by_key(sheet_id).sheet1

        # 4. Format row data safely
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        log_row = [
            timestamp,
            str(status),
            str(input_type) if input_type else "",
            str(query_label) if query_label else "",
            str(chrom) if chrom else "",
            int(target_start) if target_start else "",
            int(target_end) if target_end else "",
            int(flank_size) if flank_size is not None else "",
            int(target_padding) if target_padding is not None else "",
            bool(enable_snp_filtering) if enable_snp_filtering is not None else "",
            str(snp_source) if snp_source else "None",
            f"{maf_threshold_pct:.4f}%" if maf_threshold_pct is not None else "",
            int(min_allele_n) if min_allele_n is not None else "",
            bool(enable_segmentation) if enable_segmentation is not None else "",
            int(split_threshold) if split_threshold is not None else "",
            int(min_overlap) if min_overlap is not None else "",
            str(error_msg) if error_msg else ""
        ]

        # 5. Append row to sheet
        sheet.append_row(log_row)

    except Exception as e:
        # Fail silently to avoid breaking user interface on logging issues
        print(f"[Logger Exception] Failed to log query: {e}")