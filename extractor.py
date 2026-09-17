import os
import re
import pdfplumber
import pandas as pd
import numpy as np

def format_chemical_formula(formula_str):
    if not formula_str:
        return ''
    f_str = str(formula_str).replace('\n', ' ').strip()
    words = f_str.split()
    letters = [w for w in words if w.isalpha()]
    numbers = [w for w in words if w.isdigit()]
    
    if letters and numbers and len(letters) + len(numbers) == len(words):
        res = ""
        for i, l in enumerate(letters):
            res += l
            if i < len(numbers):
                res += numbers[i]
        return res
    return f_str

class UniversalPDFExtractor:
    def __init__(self, pdf_path, agency=None, analysis_type=None):
        self.pdf_path = pdf_path
        self.filename = os.path.basename(pdf_path)
        self._pdf = None
        self.agency = agency or self._identify_agency()
        self.analysis_type = analysis_type or self._identify_analysis_type()

    def open_pdf(self):
        if self._pdf is None:
            self._pdf = pdfplumber.open(self.pdf_path)
        return self._pdf

    def close(self):
        if self._pdf is not None:
            self._pdf.close()
            self._pdf = None

    def _identify_agency(self):
        try:
            pdf = self.open_pdf()
            if len(pdf.pages) > 0:
                first_page = (pdf.pages[0].extract_text() or "").lower()
                for name in ['eurofins', 'namsa', 'wuxi', 'ul']:
                    if name in first_page:
                        return name.capitalize()
                
                # Check Jordi
                if 'jordi' in first_page:
                    return 'Jordi'
                    
            # Check footer or body text on page 2 or page 1 for Jordi
            pages_to_check = pdf.pages[:min(5, len(pdf.pages))]
            for page in pages_to_check:
                text = (page.extract_text() or "").lower()
                if 'jordi' in text or 'support@jordilabs.com' in text:
                    return 'Jordi'
        except Exception as e:
            print(f"Error identifying agency for {self.filename}: {e}")
        return 'Unknown'

    def _identify_analysis_type(self):
        if self.agency not in ['Namsa', 'Wuxi']:
            return "Unknown"
            
        try:
            pdf = self.open_pdf()
            pages_to_check = pdf.pages[:min(2, len(pdf.pages))]
            for page in pages_to_check:
                text = (page.extract_text() or "").lower()
                if 'exhaustive' in text:
                    return 'Exhaustive'
                elif 'exaggerated' in text:
                    return 'Exaggerated'
                elif 'simulated' in text:
                    return 'Simulated'
        except Exception as e:
            print(f"Error identifying analysis type for {self.filename}: {e}")
            
        return "Unknown"

    def find_captions(self, max_pages=150, pages_subset=None):
        found = []
        try:
            import fitz
            doc = fitz.open(self.pdf_path)
            # ----------------- EUROFINS -----------------
            if self.agency == 'Eurofins':
                for page_idx in range(min(80, len(doc))):
                    if pages_subset is not None and page_idx not in pages_subset:
                        continue
                    text = doc.load_page(page_idx).get_text("text") or ""
                    if "INDEX OF TABLES" in text.upper():
                        for line in text.split('\n'):
                            m = re.search(r'(Table\s+[A-Z0-9.-]+.*?)(?:[\s\.\u2026_\-]{1,})(\d+)$', line.strip(), re.I)
                            if m:
                                cap_text = m.group(1).strip()
                                page_num = int(m.group(2))
                                
                                # Validation logic for Eurofins
                                is_valid = False
                                reason = "Does not match E&L analytical method, solvent, or result format"
                                
                                # Method regex
                                method_match = bool(re.search(r'GC/MS|LC/MS|LC/UV|ICP/MS', cap_text, re.I))
                                # Solvent regex
                                solvent_match = bool(re.search(r'Water|Ethanol|Isopropyl\s+alcohol', cap_text, re.I))
                                # Results regex
                                results_match = bool(re.search(r'Results', cap_text, re.I))
                                # Exclusions regex
                                has_exclusions = bool(re.search(r'suitability|limits|conditions|internal\s+standard|extraction\s+summary|summary\s+of\s+extraction', cap_text, re.I))
                                
                                if method_match and solvent_match and results_match and not has_exclusions:
                                    is_valid = True
                                    reason = "Valid E&L Analytical Results Table"
                                else:
                                    if has_exclusions:
                                        reason = "Excluded (suitability, limits, conditions, or summary table)"
                                    elif not results_match:
                                        reason = "Reference or method configuration table (no 'Results')"
                                    elif not method_match:
                                        reason = "Missing valid analytical method"
                                    elif not solvent_match:
                                        reason = "Missing valid extraction solvent"
                                
                                found.append({
                                    'page': page_num,
                                    'caption': cap_text,
                                    'is_valid': is_valid,
                                    'validation_reason': reason
                                })
                        if found:
                            break
                return found

            # ----------------- WUXI, NAMSA, UL, JORDI -----------------
            table_regex = re.compile(r'^(Table)\s+([A-Z]?\d+(?:\s+\(cont\.?\)|\b)?)(?:\s*[:\.-]|\s+|$)(.*)$', re.I)
            
            for i in range(min(len(doc), max_pages)):
                if pages_subset is not None and i not in pages_subset:
                    continue
                text = doc.load_page(i).get_text("text") or ""
                lines = [line.strip() for line in text.split('\n') if line.strip()]
                
                for line_idx, line in enumerate(lines):
                    m = table_regex.match(line)
                    if m:
                        tag = m.group(1)
                        tbl_id_part = m.group(2).strip()
                        rest = m.group(3).strip()
                        
                        rest_clean = re.sub(r'^[:\.-]\s*', '', rest).strip()
                        has_letters = bool(re.search(r'[a-zA-Z]', rest_clean))
                        
                        caption_parts = [f"Table {tbl_id_part}"]
                        if has_letters:
                            caption_parts.append(rest_clean)
                        else:
                            # Peek at exactly next non-empty line
                            if line_idx + 1 < len(lines):
                                next_line = lines[line_idx + 1]
                                is_invalid = (
                                    re.match(r'^(Table|Page|Analytical|Report|Eurofins|LIMS|support@)\b', next_line, re.I) or
                                    re.match(r'^\d+$', next_line) or
                                    len(next_line) > 120 or
                                    '|' in next_line or
                                    next_line.startswith('Standards ') or
                                    next_line.startswith('Device Name')
                                )
                                if not is_invalid:
                                    caption_parts.append(next_line)
                                    
                        caption_text = " ".join(caption_parts)
                        caption_text = re.sub(r'\s+', ' ', caption_text).strip()
                        
                        # Standardize format: Table X: Description
                        m_format = re.match(r'^Table\s+([A-Z]?\d+(?:\s+\(cont\.?\))?)(?:\s*[:\.-]\s*|\s+)(.+)$', caption_text, re.I)
                        if m_format:
                            tid = m_format.group(1)
                            desc = m_format.group(2).strip()
                            caption_text = f"Table {tid}: {desc}"
                            
                        # Skip footnote/sentence lines captured as captions
                        footnote_patterns = [
                            r'^Table\s+\d+:\s+All detected',
                            r'^Table\s+\d+:\s+The\s+',
                            r'^Table\s+\d+:\s+Note:\s*',
                        ]
                        is_footnote = any(re.match(p, caption_text, re.I) for p in footnote_patterns)
                        if is_footnote:
                            continue

                        # Validation logic per agency
                        is_valid = True
                        reason = "Valid E&L Analytical Table"

                        if self.agency == 'Ul':
                            is_valid = bool(re.search(r'concentration', caption_text, re.I))
                            if is_valid:
                                reason = "Valid E&L Analytical Results Table (Substance Concentrations)"
                            else:
                                reason = "Method or parameters configuration table (no 'CONCENTRATIONS')"

                        elif self.agency in ['Namsa', 'Wuxi']:
                            # Valid: must have 'Results' + recognized method; no setup/reference exclusions
                            has_results = bool(re.search(r'\bResults\b', caption_text, re.I))
                            has_method  = bool(re.search(
                                r'GC.?MS|LC.?MS|LC.?UV|ICP.?MS|QTOF|Headspace|Direct Inject',
                                caption_text, re.I
                            ))
                            has_exclusion = bool(re.search(
                                r'Parameters|Suitability|Observations?|\bSummary\b|DBT Selection|AET|Extraction Param',
                                caption_text, re.I
                            ))
                            if has_results and has_method and not has_exclusion:
                                is_valid = True
                                reason = "Valid E&L Analytical Results Table"
                            else:
                                is_valid = False
                                if has_exclusion:
                                    reason = "Excluded (parameters, suitability, observations, or summary table)"
                                elif not has_method:
                                    reason = "No recognized analytical method in caption"
                                elif not has_results:
                                    reason = "Caption does not describe an analytical results table"
                                else:
                                    reason = "Does not match E&L results table criteria"

                        elif self.agency == 'Jordi':
                            # Jordi valid tables: "Summary of QTOF-GCMS/LCMS Results", "HS-GCMS Results", "ICPMS Results"
                            # NOTE: Do NOT exclude 'Summary' — Jordi uses "Summary of [METHOD] Results" as valid captions
                            has_results = bool(re.search(r'\bResults\b', caption_text, re.I))
                            has_method  = bool(re.search(
                                r'QTOF|GC.?MS|LC.?MS|ICP.?MS|HS.?GC|ICPMS|GCMS|LCMS',
                                caption_text, re.I
                            ))
                            # Jordi exclusions: only genuine setup tables (no 'Summary' here)
                            has_exclusion = bool(re.search(
                                r'Parameters|Suitability|Observations?|DBT Selection|AET|Extraction Param'
                                r'|Identification of Test|List of Acronyms|Study Design|Device Characteristics'
                                r'|Gravimetric|Spike Study|Molecular Weight|Weight Percent|Standard and Solvent'
                                r'|Scope of ISO|PDMS Quantitation',
                                caption_text, re.I
                            ))
                            if has_results and has_method and not has_exclusion:
                                is_valid = True
                                reason = "Valid E&L Analytical Results Table"
                            else:
                                is_valid = False
                                if has_exclusion:
                                    reason = "Excluded (reference or setup table)"
                                elif not has_method:
                                    reason = "No recognized analytical method in caption"
                                elif not has_results:
                                    reason = "Caption does not describe an analytical results table"
                                else:
                                    reason = "Does not match E&L results table criteria"

                        if not any(f['page'] == i+1 and f['caption'] == caption_text for f in found):
                            if len(caption_text) > 5:
                                if not any(caption_text.lower().startswith(x) for x in ['table of contents', 'table contents']):
                                    found.append({
                                        'page': i+1,
                                        'caption': caption_text,
                                        'is_valid': is_valid,
                                        'validation_reason': reason
                                    })
            doc.close()
            return found
        except Exception as e:
            print(f"Error finding captions for {self.filename}: {e}")
            return []

    def _normalize_text(self, text):
        if not isinstance(text, str):
            return ""
        text = re.sub(r'Page \d+ of \d+|Analytical Report|Eurofins Number|Test Code:.*|LIMS Sample Number:.*|support@jordilabs\.com', '', text, flags=re.I)
        return re.sub(r'\s+', ' ', text).strip()

    def extract_table(self, page_num: int, caption: str, caption_bbox: tuple = None) -> pd.DataFrame:
        pdf = self.open_pdf()
        page = pdf.pages[page_num - 1]
        
        # If we know where the caption is, we can filter to only tables physically below it
        if caption_bbox:
            caption_y0 = caption_bbox[1]
            # Use find_tables to get tables with their bounding boxes
            found_tables = page.find_tables()
            valid_tables = [t for t in found_tables if t.bbox[1] >= caption_y0 - 15]
            
            if valid_tables:
                # Sort from top to bottom
                valid_tables.sort(key=lambda t: t.bbox[1])
                # We extract ONLY the first table right below the caption!
                # If a table is split into pieces, this will only grab the first piece,
                # but Namsa tables are usually contained within a single bounding box.
                t_obj = valid_tables[0]
                tables = [t_obj.extract()]
            else:
                tables = page.extract_tables()
        else:
            tables = page.extract_tables()
            
        if not tables:
            return None
                
        caption_idx_on_page = 0
        if len(tables) > 1 and caption:
            m = re.search(r'(Table\s+[A-Z0-9.-]+)', caption, re.IGNORECASE)
            if m:
                table_id = m.group(1).upper().rstrip('.')
                
                # 1. Check if the table ID is literally inside the first few rows of the table
                found_inside = False
                for i, t in enumerate(tables):
                    df = pd.DataFrame(t)
                    head_text = " ".join(df.head(3).astype(str).fillna('').values.flatten()).upper()
                    if table_id in head_text:
                        caption_idx_on_page = i
                        found_inside = True
                        break
                
                # 2. If not found inside, map based on the order of table IDs found in the page text
                if not found_inside:
                    text = page.extract_text() or ""
                    page_table_ids = []
                    
                    # Find all table IDs in order of appearance
                    for match in re.finditer(r'Table\s+([A-Z0-9.-]+)', text, re.IGNORECASE):
                        tid = f"TABLE {match.group(1).upper().rstrip('.')}"
                        if tid not in page_table_ids:
                            page_table_ids.append(tid)
                    
                    if table_id in page_table_ids:
                        caption_idx_on_page = page_table_ids.index(table_id)

        # ----------------- EUROFINS EXTRACTION -----------------
        if self.agency == 'Eurofins':
            # For Eurofins, tables are bordered so extract_tables() is usually 1:1 with captions
            if caption_idx_on_page < len(tables):
                tables = [tables[caption_idx_on_page]]

            
            extracted_rows = []
            current_config = "Unknown"
            
            for t in tables:
                df = pd.DataFrame(t)
                num_cols = df.shape[1]
                
                for idx, row in df.iterrows():
                    clean_row = [str(val).strip() if pd.notnull(val) else "" for val in row.tolist()]
                    row_text_upper = " ".join(clean_row).upper()
                    
                    # Skip header rows
                    if "RETENTION" in row_text_upper or "TENTATIVE IDENTIFICATION" in row_text_upper or "ELEMENT EST. CONC" in row_text_upper:
                        continue
                        
                    # Check for configuration row
                    non_empty_indices = [i for i, val in enumerate(clean_row) if val]
                    if len(non_empty_indices) == 1 and non_empty_indices[0] == 0:
                        val = clean_row[0]
                        if any(x in val for x in ["SU ", "MPU ", "Configuration"]):
                            current_config = val
                            continue
                            
                    if not any(clean_row):
                        continue
                        
                    extracted_rows.append({
                        "Configuration": current_config,
                        "RawData": clean_row,
                        "ColsCount": num_cols
                    })
                    
            if not extracted_rows:
                return None
                
            # Map rows based on columns count
            final_data = []
            for r in extracted_rows:
                raw = [re.sub(r'\n', ' ', v).strip() if v else v for v in r["RawData"]]
                cfg = r["Configuration"]
                cols_cnt = r["ColsCount"]
                
                if cols_cnt == 6: # GC/MS
                    final_data.append({
                        "Sample Configuration": cfg,
                        "Retention Time (min)": raw[0],
                        "Compound / Identification": raw[1],
                        "Est Conc (ug/mL)": raw[2],
                        "Est Conc (ug/cm2)": raw[3],
                        "Est Conc (ug/device)": raw[4],
                        "Confidence (%)": raw[5] if len(raw) > 5 else ""
                    })
                elif cols_cnt == 5: # LC/MS/UV
                    final_data.append({
                        "Sample Configuration": cfg,
                        "Retention Time (min)": raw[0],
                        "Compound / Identification": raw[1],
                        "Est Conc (ug/mL)": raw[2],
                        "Est Conc (ug/cm2)": raw[3],
                        "Est Conc (ug/g)": raw[4] if len(raw) > 4 else ""
                    })
                elif cols_cnt == 4: # ICP/MS
                    final_data.append({
                        "Sample Configuration": cfg,
                        "Element": raw[0],
                        "Est Conc (ug/mL)": raw[1],
                        "Est Conc (ug/cm2)": raw[2],
                        "Est Conc (ug/device)": raw[3] if len(raw) > 3 else ""
                    })
                else:
                    row_dict = {"Sample Configuration": cfg}
                    for c_idx, val in enumerate(raw):
                        row_dict[f"Col_{c_idx}"] = val
                    final_data.append(row_dict)
                    
            return pd.DataFrame(final_data)

        # ----------------- NAMSA & WUXI EXTRACTION -----------------
        elif self.agency in ['Namsa', 'Wuxi']:
            print(f"DEBUG: Namsa extraction started. Num tables={len(tables)}")
            result_df = None
            
            # Namsa tables are borderless, so pdfplumber returns arbitrary chunks.
            # First, find all chunks that look like the start of a valid Namsa table
            valid_chunks = []
            for idx_t, t in enumerate(tables):
                df = pd.DataFrame(t)
                if df.shape[1] >= 6: # Main tables have 8-10 columns
                    row0_str = " ".join(df.iloc[0].astype(str).fillna('').values).upper()
                    row1_str = " ".join(df.iloc[1].astype(str).fillna('').values).upper() if df.shape[0] > 1 else ""
                    if "CASRN" in row0_str or "COMPOUND" in row0_str or "CASRN" in row1_str:
                        valid_chunks.append(df)
                elif df.shape[1] >= 3: # Elemental tables have fewer columns
                    row0_str = " ".join(df.iloc[0].astype(str).fillna('').values).upper()
                    if "ELEMENT" in row0_str:
                        valid_chunks.append(df)
                        
            if valid_chunks:
                # Select the chunk that corresponds to the Nth caption on this page
                idx_to_use = caption_idx_on_page if caption_idx_on_page < len(valid_chunks) else -1
                result_df = valid_chunks[idx_to_use]
            
            if result_df is None:
                # fallback just in case
                if tables:
                    result_df = pd.DataFrame(tables[0])
                else:
                    return None
                    
            # Find header boundary
            header_rows_count = 0
            for idx in range(min(4, result_df.shape[0])):
                cells_upper = [str(c).upper().strip() for c in result_df.iloc[idx].fillna('').values]
                row_str = " ".join(cells_upper)
                if any(x in row_str for x in ["RETENTION TIME", "CASRN", "ID CATEGORY", "UG/DEVICE", "UG/ML"]) or any(c == "ELEMENT" for c in cells_upper):
                    header_rows_count = idx + 1
                    
            if header_rows_count == 0 and result_df.shape[0] > 0:
                row0_has_numbers = any(re.search(r'\d', str(c)) for c in result_df.iloc[0].fillna('').values)
                if not row0_has_numbers:
                    header_rows_count = 1
                    
            parsed_rows = []
            active_id_cat = ""
            active_casrn = ""
            active_compound = ""
            active_formula = ""
            active_surrogate = ""
            active_mol_wt = ""
            active_element = ""
            
            for idx in range(header_rows_count, result_df.shape[0]):
                row = result_df.iloc[idx].tolist()
                clean_row = [str(val).strip() if pd.notnull(val) and str(val).strip() != 'nan' else '' for val in row]
                
                if not any(clean_row):
                    continue
                    
                if result_df.shape[1] < 6:
                    # Elemental table parsing
                    element = clean_row[0]
                    if element: active_element = element
                    else: element = active_element
                    
                    conc = ""
                    for val in clean_row[1:]:
                        if val and re.search(r'[\d.]+', val):
                            conc = val
                            break
                            
                    if not element or not conc:
                        continue
                        
                    parsed_rows.append({
                        "Element": element,
                        "Est Conc (ug/device)": conc
                    })
                    continue
                    
                id_cat = clean_row[0]
                casrn = clean_row[1]
                compound = clean_row[2]
                formula = clean_row[3]
                
                if id_cat: active_id_cat = id_cat
                if casrn: active_casrn = casrn
                if compound: active_compound = compound
                if formula: active_formula = formula
                
                mol_wt = clean_row[4] if len(clean_row) > 4 else ""
                if mol_wt: active_mol_wt = mol_wt
                
                ret_time = ""
                for col_idx in [5, 6, 7]:
                    if col_idx < len(clean_row) and clean_row[col_idx] and not re.match(r'^[a-zA-Z]', clean_row[col_idx]):
                        ret_time = clean_row[col_idx]
                        break
                        
                ug_device = ""
                if len(clean_row) > 8 and clean_row[8]:
                    ug_device = clean_row[8]
                elif len(clean_row) > 7 and clean_row[7]:
                    ug_device = clean_row[7]
                    
                surrogate = ""
                if len(clean_row) > 9 and clean_row[9]:
                    surrogate = clean_row[9]
                    active_surrogate = surrogate
                else:
                    surrogate = active_surrogate
                    
                # Clean newline characters from text fields
                clean_compound = re.sub(r'\n', ' ', active_compound).strip()
                clean_formula = re.sub(r'\n', ' ', active_formula).strip()
                clean_surrogate = re.sub(r'\n', ' ', active_surrogate).strip()
                    
                parsed_rows.append({
                    "ID Category": active_id_cat,
                    "CASRN": active_casrn,
                    "Compound": clean_compound,
                    "Formula": clean_formula,
                    "Molecular Weight": active_mol_wt,
                    "Retention Time (min)": ret_time,
                    "Concentration (ug/device)": ug_device,
                    "Surrogate Used": clean_surrogate
                })
                
            print(f"DEBUG: parsed_rows length = {len(parsed_rows)}")
            if not parsed_rows:
                print(f"DEBUG: Returning None because parsed_rows is empty!")
                return None
                
            return pd.DataFrame(parsed_rows)

        # ----------------- UL EXTRACTION -----------------
        elif self.agency == 'Ul':
            text = page.extract_text() or ""
            
            if "Sample Name:" in text and "External Standard Report" in text:
                lines = text.split('\n')
                parsed_rows = []
                start_parsing = False
                
                for line in lines:
                    if "-------|" in line or "-------|--|------|" in line:
                        start_parsing = True
                        continue
                    if start_parsing:
                        if "Totals" in line or "====" in line or not line.strip():
                            if "Totals" in line:
                                break
                            continue
                            
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            ret_time = parts[0]
                            sig = parts[1]
                            
                            last_part = parts[-1]
                            compound_name = ""
                            cas = ""
                            
                            if '","' in last_part:
                                subparts = last_part.split('","')
                                compound_name = subparts[0].replace('"', '')
                                cas = subparts[1].replace('"', '')
                            elif ',' in last_part:
                                subparts = last_part.split(',')
                                compound_name = subparts[0].replace('"', '')
                                cas = subparts[1].replace('"', '')
                            else:
                                compound_name = last_part
                                cas = "N/A"
                                
                            type_val = "-"
                            area_val = "-"
                            amt_area = "-"
                            amount = "-"
                            
                            if len(parts) >= 7:
                                type_val = parts[2]
                                area_val = parts[3]
                                amt_area = parts[4]
                                amount = parts[5]
                            elif len(parts) == 6:
                                type_val = parts[2]
                                area_val = parts[3]
                                amt_area = parts[4]
                                
                            parsed_rows.append({
                                "Retention Time (min)": ret_time,
                                "Signal": sig,
                                "Type": type_val,
                                "Area": area_val,
                                "Amt/Area": amt_area,
                                "Amount (ug/M^3)": amount,
                                "Compound": compound_name,
                                "CASRN": cas
                            })
                return pd.DataFrame(parsed_rows) if parsed_rows else None
                
            for t in tables:
                df = pd.DataFrame(t)
                if df.shape[1] >= 8 and "PEAK_NUM" in str(df.iloc[0,0]).upper():
                    cols = [str(col).strip().replace('\n', ' ') for col in df.iloc[0]]
                    data_df = df.iloc[1:].copy()
                    data_df.columns = cols
                    return data_df
                    
                row0_str = " ".join(df.iloc[0].astype(str).fillna('').values).upper()
                if "CONCENTRATIONS OF" in row0_str:
                    parsed_rows = []
                    for idx in range(5, df.shape[0]):
                        row = df.iloc[idx].tolist()
                        clean_row = [str(val).strip() if pd.notnull(val) and str(val).strip() != 'nan' else '' for val in row]
                        if not any(clean_row):
                            continue
                        parsed_rows.append({
                            "CASRN": clean_row[0].replace('\n', ' '),
                            "Compound": clean_row[1].replace('†', '').replace('‡', '').replace('\n', ' ').strip(),
                            "Initial Operation (ug/m3)": clean_row[3] if len(clean_row) > 3 else "",
                            "24 hr (ug/m3)": clean_row[4] if len(clean_row) > 4 else "",
                            "72 hr (ug/m3)": clean_row[5] if len(clean_row) > 5 else "",
                            "168 hr (ug/m3)": clean_row[6] if len(clean_row) > 6 else ""
                        })
                    return pd.DataFrame(parsed_rows) if parsed_rows else None
                    
            return None

        # ----------------- JORDI EXTRACTION -----------------
        elif self.agency == 'Jordi':
            # Find the main data table (largest table with relevant keywords)
            extracted_df = None
            valid_chunks = []
            for t in tables:
                df = pd.DataFrame(t)
                if df.shape[0] < 2 or df.shape[1] < 3:
                    continue
                all_text = " ".join(df.astype(str).fillna('').values.flatten()).upper()
                if any(k in all_text for k in ["RETENTION", "RT", "CAS", "IDENTIFICATION", "ELEMENT", "PROPOSED"]):
                    valid_chunks.append(df)
                    
            if valid_chunks:
                idx_to_use = caption_idx_on_page if caption_idx_on_page < len(valid_chunks) else -1
                extracted_df = valid_chunks[idx_to_use]
                
            if extracted_df is None:
                if len(tables) > 0:
                    extracted_df = pd.DataFrame(tables[0])
                else:
                    return None
            
            # --- Identify table type from caption ---
            caption_upper = caption.upper()
            is_icpms = "ICPMS" in caption_upper or "ICP-MS" in caption_upper
            is_hsgcms = "HS-GCMS" in caption_upper or "HS-GC" in caption_upper
            
            # --- Find the first data row (skip title + header rows) ---
            first_data_idx = 0
            for idx in range(min(8, extracted_df.shape[0])):
                row_vals = extracted_df.iloc[idx].tolist()
                clean_vals = [str(v).strip() if pd.notnull(v) else '' for v in row_vals]
                row_text = " ".join(clean_vals).upper()
                
                # Skip title rows ("Table 17", "Summary of QTOF-...")
                if re.search(r'^TABLE\s+\d', row_text.strip()):
                    continue
                if "SUMMARY OF" in row_text or "RESULTS" in row_text:
                    if not any(c.replace('.','',1).isdigit() for c in clean_vals if c):
                        continue
                
                # Skip header rows (contain keywords but no numeric RT/element data)
                if any(k in row_text for k in ["RT", "RETENTION", "PROPOSED", "ELEMENT", "FORMULA", 
                                                 "CONFIDENCE", "MASS PER", "REPLICATE", "µG", "REP1", "REP 1",
                                                 "(MIN)", "DEVICE", "QUANTITATION", "METHOD BLANK",
                                                 "MEASURED", "CORRECTED", "CONCENTRATION", "NG/ML"]):
                    continue
                
                # If we reach here, this looks like a data row
                # Verify it has at least some non-empty content
                non_empty = [v for v in clean_vals if v and v.upper() != 'NONE']
                if len(non_empty) >= 2:
                    first_data_idx = idx
                    break
            
            # --- Extract data rows ---
            parsed_rows = []
            for idx in range(first_data_idx, extracted_df.shape[0]):
                row_vals = extracted_df.iloc[idx].tolist()
                clean_vals = [str(v).strip() if pd.notnull(v) else '' for v in row_vals]
                row_text = " ".join(clean_vals)
                
                # Skip footnote rows (start with superscript digits, "X:", calibration info)
                if re.match(r'^[¹²³⁴⁵\d]', clean_vals[0]) and len(clean_vals[0]) > 20:
                    continue
                if clean_vals[0].startswith('X:') or clean_vals[0].startswith('*'):
                    continue
                if 'calibration' in row_text.lower() or 'retention times are from' in row_text.lower():
                    continue
                
                # Skip empty rows
                non_empty = [v for v in clean_vals if v]
                if len(non_empty) < 2:
                    continue
                    
                # Skip table header rows that mistakenly get parsed
                if 'rt' in str(clean_vals[0]).lower() or 'min' in str(clean_vals[0]).lower():
                    continue
                if 'element' in str(clean_vals[0]).lower():
                    continue
                if len(clean_vals) > 2 and 'proposed' in str(clean_vals[2]).lower():
                    continue
                
                # --- Parse based on table type ---
                if is_icpms:
                    # ICPMS: Element | Method Blank | Rep1 Measured | Rep1 Corrected | Rep2... | Rep3...
                    element = clean_vals[0] if clean_vals[0] else ''
                    if not element:
                        continue
                    row_dict = {"Element": re.sub(r'\n', ' ', element).strip()}
                    # Collect remaining numeric values
                    nums = [v for v in clean_vals[1:] if v]
                    for ni, nv in enumerate(nums):
                        row_dict[f"Value_{ni+1}"] = re.sub(r'\n', ' ', nv).strip()
                    parsed_rows.append(row_dict)
                    
                # We MUST use non_empty to avoid index shifting due to empty PDF columns
                non_empty_vals = [v for v in clean_vals if v and str(v).upper() != 'NONE']
                
                if is_hsgcms:
                    # HS-GCMS: RT | Proposed ID | Formula | Confidence | CAS | Detected flags
                    if len(non_empty_vals) < 3:
                        continue
                    rt = non_empty_vals[0]
                    proposed_id = non_empty_vals[1] if len(non_empty_vals) > 1 else ''
                    formula = format_chemical_formula(non_empty_vals[2] if len(non_empty_vals) > 2 else '')
                    confidence = re.sub(r'[\s\n]+', '', non_empty_vals[3] if len(non_empty_vals) > 3 else '')
                    cas = non_empty_vals[4] if len(non_empty_vals) > 4 else ''
                    # Collect detection flags for replicates
                    det_flags = [v for v in non_empty_vals[5:] if v]
                    
                    parsed_rows.append({
                        "RT (min)": re.sub(r'\n', ' ', rt).strip(),
                        "Proposed Identification": re.sub(r'\n', ' ', proposed_id).strip(),
                        "Formula": format_chemical_formula(re.sub(r'\n', ' ', formula).strip()),
                        "Confidence Level": re.sub(r'\n', ' ', confidence).strip(),
                        "CAS": re.sub(r'\n', ' ', cas).strip(),
                        "Detection": ", ".join(det_flags)
                    })
                    
                else:
                    # GCMS / LCMS: RT | Proposed ID | Formula | Confidence | CAS | Rep1 | ... | Rep2 | ... | Rep3 | ... | Quant Std
                    if len(non_empty_vals) < 3:
                        continue
                    rt = non_empty_vals[0]
                    proposed_id = non_empty_vals[1] if len(non_empty_vals) > 1 else ''
                    formula = format_chemical_formula(non_empty_vals[2] if len(non_empty_vals) > 2 else '')
                    confidence = re.sub(r'[\s\n]+', '', non_empty_vals[3] if len(non_empty_vals) > 3 else '')
                    cas = non_empty_vals[4] if len(non_empty_vals) > 4 else ''
                    
                    # Extract numeric values (replicate concentrations)
                    numeric_vals = []
                    for v in non_empty_vals[5:]:
                        if v:
                            clean_v = re.sub(r'\n', ' ', v).strip()
                            # Only match if the entire string consists of valid numeric characters or specific flags
                            if re.match(r'^[\d.<>NDBQL\s,]+$', clean_v) or clean_v in ['ND', 'BQL']:
                                numeric_vals.append(clean_v)
                    
                    # Last text value(s) are typically the Quantitated By and Quantitation Standard
                    # Because they are text at the end after numeric replicates, we can find them.
                    quant_by = ''
                    quant_std = ''
                    text_vals_at_end = []
                    for v in reversed(non_empty_vals[5:]):
                        if v and re.search(r'[a-zA-Z]', str(v)):
                            # Ensure we don't grab something that is actually just a numeric replicate like "ND" or "< LOD"
                            # We'll just grab the strings at the end
                            text_vals_at_end.append(re.sub(r'\n', ' ', str(v)).strip())
                        else:
                            if text_vals_at_end:
                                # We hit numbers, so we got all trailing text
                                break
                    
                    # Re-reverse so it's left-to-right
                    text_vals_at_end.reverse()
                    
                    if len(text_vals_at_end) == 1:
                        quant_std = text_vals_at_end[0]
                    elif len(text_vals_at_end) >= 2:
                        # Typically second to last is Quantitated By, last is Quantitation Standard
                        quant_by = text_vals_at_end[-2]
                        quant_std = text_vals_at_end[-1]
                    
                    row_dict = {
                        "RT (min)": re.sub(r'\n', ' ', rt).strip(),
                        "Proposed Identification": re.sub(r'\n', ' ', proposed_id).strip(),
                        "Formula": re.sub(r'\n', ' ', formula).strip(),
                        "Confidence Level": confidence,
                        "CAS": cas
                    }
                    
                    # Add replicate values
                    for ri, rv in enumerate(numeric_vals):
                        row_dict[f"Rep{ri+1} (ug/device)"] = rv
                    row_dict["Quantitated by"] = quant_by
                    row_dict["Quantitation Standard"] = quant_std
                    parsed_rows.append(row_dict)
            
            if not parsed_rows:
                return None
                
            return pd.DataFrame(parsed_rows)
            
        return None


def _extract_table_worker(temp_pdf, page, caption_text, caption_bbox, agency, analysis_type):
    import traceback
    # from extractor import UniversalPDFExtractor  # Already in this file
    extractor = None
    try:
        extractor = UniversalPDFExtractor(temp_pdf, agency=agency, analysis_type=analysis_type)
        df = extractor.extract_table(page, caption_text, caption_bbox=caption_bbox)
        return df, None
    except Exception as e:
        return None, str(e)
    finally:
        if extractor is not None:
            extractor.close()
