import pandas as pd
import json
import os
import sys

# Append current directory to path to allow imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from main import assess_chemical
from schemas import AssessmentRequest
from database import SessionLocal

# Read input from the Wuxi max concentration excel file
import time
import math

namsa_file = r"C:\Users\kiran\Documents\Projects\xchem\output\Wuxi_ER 2245638 v00 - Attachment F - 1068837 E&L Analytical Report (2) (2)_max_concentration.xlsx"
df = pd.read_excel(namsa_file)

# Build inputs list from the dataframe
inputs = []
for _, row in df.iterrows():
    compound_name = str(row['Cleaned Compound Name']).strip()
    if compound_name == 'nan' or not compound_name:
        compound_name = str(row['Compound Name']).strip()
        
    cas_number = str(row['CASRN']).strip()
    if cas_number == 'nan':
        cas_number = None
        
    formula = str(row['Formula']).strip()
    if formula == 'nan':
        formula = None
        
    conc = row['Concentration']
    try:
        conc_float = float(conc)
    except:
        conc_float = None

    if compound_name and compound_name != 'nan':
        inputs.append({
            "Compound Name": compound_name,
            "CAS Number": cas_number,
            "Molecular Formula": formula,
            "Concentration": conc_float
        })

print(f"Loaded {len(inputs)} records from Namsa file.")

db = SessionLocal()

all_validation = []
all_identity = []
all_structure = []
all_fgs = []
all_moiety_class = []
all_analogs = []
all_final = []

print("Running pipeline for test compounds...")
for idx, item in enumerate(inputs):
    safe_name = item['Compound Name'].encode('ascii', 'ignore').decode('ascii')
    print(f"Processing: {safe_name}")
    req = AssessmentRequest(
        compound_name=item["Compound Name"],
        cas_number=item.get("CAS Number"),
        molecular_formula=item.get("Molecular Formula"),
        concentration=item.get("Concentration")
    )
    
    try:
        # Be polite to PubChem API
        time.sleep(0.5)
        res = assess_chemical(req, db)
        res_dict = res.model_dump()
        
        # Step 1: Validation (Implicitly valid if no exception)
        all_validation.append({"Compound Name": item["Compound Name"], "Validation Status": "Valid"})
        
        # Step 2: Identity
        id_info = res_dict["identity"]
        id_info["Input Name"] = item["Compound Name"]
        all_identity.append(id_info)
        
        # Step 3: Structure
        struct = res_dict["structure"]
        struct["Compound Name"] = item["Compound Name"]
        all_structure.append(struct)
        
        # Step 4: Functional Groups
        for fg in res_dict["functional_groups"]:
            fg_copy = dict(fg)
            fg_copy["Compound Name"] = item["Compound Name"]
            fg_copy["matched_atoms"] = str(fg_copy["matched_atoms"])
            all_fgs.append(fg_copy)
            
        # Step 5 & 6: Moiety & Class
        all_moiety_class.append({
            "Compound Name": item["Compound Name"],
            "Parent Moiety": res_dict["parent_moiety"],
            "Chemical Class": res_dict["chemical_class"]
        })
        
        # Step 8, 9, 10: Analogs
        for an in res_dict["analogs"]:
            an_copy = dict(an)
            an_copy["Target Compound"] = item["Compound Name"]
            all_analogs.append(an_copy)
            
        # Step 11, 12: Final
        all_final.append({
            "Compound Name": item["Compound Name"],
            "Confidence Score": res_dict["confidence"],
            "Rationale": res_dict["rationale"]
        })
        
    except Exception as e:
        all_validation.append({"Compound Name": item["Compound Name"], "Validation Status": f"Invalid: {str(e)}"})

db.close()

# Save to Excel
out_dir = "stepwise_outputs"
os.makedirs(out_dir, exist_ok=True)

pd.DataFrame(all_validation).to_excel(os.path.join(out_dir, "step_01_validation.xlsx"), index=False)
pd.DataFrame(all_identity).to_excel(os.path.join(out_dir, "step_02_identity.xlsx"), index=False)
pd.DataFrame(all_structure).to_excel(os.path.join(out_dir, "step_03_structure.xlsx"), index=False)
if all_fgs:
    pd.DataFrame(all_fgs).to_excel(os.path.join(out_dir, "step_04_functional_groups.xlsx"), index=False)
pd.DataFrame(all_moiety_class).to_excel(os.path.join(out_dir, "step_05_06_moiety_class.xlsx"), index=False)
if all_analogs:
    pd.DataFrame(all_analogs).to_excel(os.path.join(out_dir, "step_08_09_10_analogs.xlsx"), index=False)
pd.DataFrame(all_final).to_excel(os.path.join(out_dir, "step_11_12_final.xlsx"), index=False)

print("Saved all stepwise Excel outputs to 'stepwise_outputs' directory.")
