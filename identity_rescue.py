import pandas as pd
import difflib
import re

def rescue_unresolved_compounds(df: pd.DataFrame, smiles_col: str = 'Standardized SMILES', name_col: str = 'Cleaned Compound Name') -> pd.DataFrame:
    """
    Tier 2 Surrogate Mapping
    Compares unresolved compounds to resolved compounds in the same batch.
    If textual similarity > 85%, passes them through Chemistry Gates.
    """
    df = df.copy()
    
    if smiles_col not in df.columns or name_col not in df.columns:
        return df

    if 'Identity Status' not in df.columns:
        df['Identity Status'] = "Resolved"
        
    # Mark initially unresolved
    df.loc[df[smiles_col].isna() | (df[smiles_col] == ""), 'Identity Status'] = "Unresolved"

    resolved_mask = df['Identity Status'] == "Resolved"
    unresolved_mask = df['Identity Status'] == "Unresolved"
    
    resolved_df = df[resolved_mask]
    if resolved_df.empty:
        return df
        
    resolved_dict = {}
    for idx, row in resolved_df.iterrows():
        name = str(row.get(name_col, '')).strip()
        if name:
            resolved_dict[name] = {
                'smiles': row.get(smiles_col, ''),
                'orig_smiles': row.get('Original SMILES', ''),
                'inchikey': row.get('InChIKey', ''),
                'inchi': row.get('InChI', ''),
                'resolved_name': row.get('Resolved Name', '')
            }

    if not resolved_dict:
        return df
        
    for idx, row in df[unresolved_mask].iterrows():
        unres_name = str(row.get(name_col, '')).strip()
        if not unres_name:
            continue
            
        best_match = None
        best_ratio = 0.0
        
        for res_name, res_data in resolved_dict.items():
            ratio = difflib.SequenceMatcher(None, unres_name.lower(), res_name.lower()).ratio()
            
            # Substrings (like "Irganox 1010 related compound") will have low ratio, so we override ratio
            is_substring = res_name.lower() in unres_name.lower() and len(res_name) > 5
            
            if ratio > best_ratio or is_substring:
                if is_substring and ratio < 0.85:
                    # Give substring matches a synthetic high score so they get picked up
                    ratio = 0.95 
                
                if ratio > best_ratio:
                    best_ratio = ratio
                    best_match = res_name
                
        if best_match and best_ratio >= 0.85:
            reason = pass_chemistry_gates(unres_name, best_match, best_ratio)
            if reason:
                res_data = resolved_dict[best_match]
                df.at[idx, smiles_col] = res_data['smiles']
                df.at[idx, 'Original SMILES'] = res_data['orig_smiles']
                df.at[idx, 'InChIKey'] = res_data['inchikey']
                df.at[idx, 'InChI'] = res_data['inchi']
                df.at[idx, 'Resolved Name'] = res_data['resolved_name']
                df.at[idx, 'Identity Status'] = f"{reason} of '{best_match}'"

    # Remaining Unresolved are flagged as Tier 3
    df.loc[(df[smiles_col].isna() | (df[smiles_col] == "")) & (df['Identity Status'] == "Unresolved"), 'Identity Status'] = "Tier 3: Class Fallback Required (No Structure)"

    return df


def pass_chemistry_gates(name_unres: str, name_res: str, ratio: float):
    nu = name_unres.lower()
    nr = name_res.lower()
    
    # Gate A: Oligomer (n=X)
    oligo_pattern = r'\(n=\d+\)'
    base_u = re.sub(oligo_pattern, '', nu).strip()
    base_r = re.sub(oligo_pattern, '', nr).strip()
    if base_u == base_r and base_u:
        return "Surrogate SMILES (Homologue)"
        
    # Gate B: Related Compound / Isomer
    suffix_pattern = r'(?i)\s*(related compound|isomer of|derivative|isomer)\s*'
    base_u2 = re.sub(suffix_pattern, '', nu).strip()
    base_r2 = re.sub(suffix_pattern, '', nr).strip()
    if base_u2 == base_r2 and base_u2:
        return "Surrogate SMILES (Related Compound/Isomer)"
        
    # Gate C: Strict Typo Rule (Requires >= 85% match)
    if ratio >= 0.85:
        # Prevent numbers from changing (e.g. 1-hexanol vs 2-hexanol)
        digits_u = re.sub(r'\D', '', nu)
        digits_r = re.sub(r'\D', '', nr)
        if digits_u != digits_r:
            return False
            
        # Prevent critical carbon-chain prefixes from shifting
        prefixes = ['meth', 'eth', 'prop', 'but', 'pent', 'hex', 'hept', 'oct', 'non', 'dec']
        for p in prefixes:
            if (p in nu) != (p in nr):
                return False
                
        return "Surrogate SMILES (Typo Corrected)"
        
    return False
