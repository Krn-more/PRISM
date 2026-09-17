import requests
import os
import time
from rdkit import Chem
from rdkit.Chem import inchi

def get_chembl_data(smiles: str) -> dict:
    """
    Queries the ChEMBL API for a given SMILES string (via InChIKey).
    Returns a dictionary of biological targets and max phase.
    """
    if not smiles:
        return {"Max_Phase": None, "Known_Targets": []}
        
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {"Max_Phase": None, "Known_Targets": []}
        
    try:
        inchi_key = inchi.MolToInchiKey(mol)
        if not inchi_key:
            return {"Max_Phase": None, "Known_Targets": []}
            
        # 1. Get Molecule ChEMBL ID
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule?molecule_structures__standard_inchi_key={inchi_key}&format=json"
        response = requests.get(url, timeout=5)
        
        if response.status_code != 200:
            return {"Max_Phase": None, "Known_Targets": []}
            
        data = response.json()
        if not data.get("molecules"):
            return {"Max_Phase": None, "Known_Targets": []}
            
        molecule = data["molecules"][0]
        chembl_id = molecule.get("molecule_chembl_id")
        max_phase = molecule.get("max_phase", 0)
        
        # 2. Get Mechanism of Action (Targets)
        mech_url = f"https://www.ebi.ac.uk/chembl/api/data/mechanism?molecule_chembl_id={chembl_id}&format=json"
        mech_response = requests.get(mech_url, timeout=5)
        
        targets = []
        if mech_response.status_code == 200:
            mech_data = mech_response.json()
            for mech in mech_data.get("mechanisms", []):
                target_name = mech.get("target_chembl_id", "Unknown Target")
                action = mech.get("action_type", "Unknown Action")
                targets.append(f"{target_name} ({action})")
                
        return {
            "Max_Phase": max_phase,
            "Known_Targets": list(set(targets))
        }
    except Exception as e:
        return {"Max_Phase": None, "Known_Targets": [], "Error": str(e)}

def fetch_chembl_for_dataframe(df, smiles_col='Standardized SMILES'):
    """
    Applies ChEMBL queries across a dataframe.
    """
    max_phases = []
    targets = []
    
    deadline = time.monotonic() + float(os.environ.get("PHASE_B_STAGE_DEADLINE_SECONDS", "300"))
    for smi in df[smiles_col]:
        if time.monotonic() >= deadline:
            remaining = len(df) - len(max_phases)
            max_phases.extend([None] * remaining)
            targets.extend(["Unavailable: StageDeadlineExceeded"] * remaining)
            break
        res = get_chembl_data(smi)
        max_phases.append(res.get("Max_Phase"))
        targets_str = ", ".join(res.get("Known_Targets", []))
        targets.append(targets_str if targets_str else "None")
        
    df['ChEMBL_Max_Phase'] = max_phases
    df['ChEMBL_Targets'] = targets
    return df
