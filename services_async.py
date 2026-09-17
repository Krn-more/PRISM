import aiohttp
import asyncio
import os
import re
from urllib.parse import quote
from schemas import AssessmentRequest
from opsin_identity import OPSIN_WEB_SERVICE, opsin_identity_from_payload, opsin_name_is_eligible
import time

class RateLimiter:
    def __init__(self, rate: float):
        self.interval = 1.0 / rate
        self.lock = asyncio.Lock()
        self.last_time = 0.0

    async def wait(self):
        async with self.lock:
            now = time.monotonic()
            elapsed = now - self.last_time
            if elapsed < self.interval:
                await asyncio.sleep(self.interval - elapsed)
            self.last_time = time.monotonic()

async def resolve_identity_async(req: AssessmentRequest, session: aiohttp.ClientSession, rate_limiter: RateLimiter = None):
    search_term = req.cas_number if req.cas_number else req.compound_name
    
    fast_fail_keywords = ["poly(", "copolymer", "related", "degradant", "adduct", "dimer", "trimer", "unknown", "unidentified"]
    if any(k in search_term.lower() for k in fast_fail_keywords):
        return {
            "status": "Unresolved",
            "source": "Fast-Fail Rules",
            "name": search_term,
            "iupac_name": "",
            "cid": "",
            "synonyms": [],
            "cas": req.cas_number or "", 
            "smiles": "",
            "inchi": "",
            "inchikey": "",
            "formula": "",
            "molecular_weight": 0.0
        }

    search_terms = []
    if req.cas_number: search_terms.append((req.cas_number.strip(), "name"))
    if req.compound_name: search_terms.append((req.compound_name.strip(), "name"))
    if not search_terms: search_terms.append((search_term, "name"))

    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
    
    for term, search_type in search_terms:
        for attempt in range(3):
            try:
                pubchem_url = f"{base_url}/{search_type}/{term}/property/Title,CanonicalSMILES,IsomericSMILES,InChI,InChIKey,MolecularFormula,MolecularWeight,IUPACName/JSON"
                if rate_limiter: await rate_limiter.wait()
                async with session.get(pubchem_url, ssl=False, timeout=10) as res:
                    if res.status == 503 or res.status == 504:
                        await asyncio.sleep(1.5 * (attempt + 1))
                        continue
                    if res.status == 404:
                        break # Not found, move to next search term
                    if res.status == 200:
                        data = await res.json()
                    else:
                        break # Other error
                        
                    props = data['PropertyTable']['Properties'][0]
                    
                    cid = str(props.get('CID', ''))
                    title = props.get('Title', '')
                    iupac = props.get('IUPACName', '')
                    smiles = props.get('CanonicalSMILES') or props.get('IsomericSMILES') or props.get('SMILES') or props.get('ConnectivitySMILES') or ''
                    
                    if not smiles:
                        try:
                            cactus_url = f"https://cactus.nci.nih.gov/chemical/structure/{quote(term)}/smiles"
                            async with session.get(cactus_url, ssl=False, timeout=3) as c_res:
                                if c_res.status == 200:
                                    text = await c_res.text()
                                    if text: smiles = text.strip().split('\\n')[0]
                        except:
                            pass

                    inchi = props.get('InChI') or ''
                    inchikey = props.get('InChIKey') or ''
                    formula = props.get('MolecularFormula') or ''
                    mw = float(props.get('MolecularWeight') or 0.0)
                    
                    synonyms = []
                    syn_url = f"{base_url}/{search_type}/{term}/synonyms/JSON"
                    for syn_attempt in range(2):
                        try:
                            if rate_limiter: await rate_limiter.wait()
                            async with session.get(syn_url, ssl=False, timeout=5) as syn_res:
                                if syn_res.status == 503:
                                    await asyncio.sleep(1)
                                    continue
                                if syn_res.status == 200:
                                    syn_data = await syn_res.json()
                                    if 'InformationList' in syn_data:
                                        synonyms = syn_data['InformationList']['Information'][0].get('Synonym', [])[:5]
                                    break
                        except:
                            pass
                        
                    status = "Exact Match" if req.compound_name and (req.compound_name.lower() in title.lower() or title.lower() in req.compound_name.lower()) else "Partial Match"
                    
                    return {
                        "status": status,
                        "source": "PubChem",
                        "name": title,
                        "iupac_name": iupac,
                        "cid": cid,
                        "synonyms": synonyms,
                        "cas": req.cas_number or "", 
                        "smiles": smiles,
                        "inchi": inchi,
                        "inchikey": inchikey,
                        "formula": formula,
                        "molecular_weight": mw
                    }
            except Exception:
                pass

    for term, search_type in search_terms:
        try:
            cactus_url = f"https://cactus.nci.nih.gov/chemical/structure/{quote(term)}/smiles"
            async with session.get(cactus_url, ssl=False, timeout=4) as c_res:
                if c_res.status == 200:
                    text = await c_res.text()
                    if text:
                        smiles = text.strip().split('\\n')[0]
                        return {
                            "status": "Partial Match",
                            "source": "NCI/Cactus",
                            "name": req.compound_name or term,
                            "iupac_name": "",
                            "cid": "",
                            "synonyms": [],
                            "cas": req.cas_number or "", 
                            "smiles": smiles,
                            "inchi": "",
                            "inchikey": "",
                            "formula": "",
                            "molecular_weight": 0.0
                        }
        except Exception:
            pass

    epa_api_key = os.getenv("EPA_API_KEY")
    if epa_api_key:
        try:
            import requests
            epa_url = f"https://comptox.epa.gov/ctx-api/chemical/search/equal/{requests.utils.quote(search_term, safe='')}"
            headers = {"x-api-key": epa_api_key, "Accept": "application/json"}
            async with session.get(epa_url, headers=headers, ssl=False, timeout=2) as epa_res:
                if epa_res.status == 200:
                    epa_data = await epa_res.json()
                    if epa_data and isinstance(epa_data, list):
                        match = epa_data[0]
                        return {
                            "status": "Partial Match",
                            "source": "EPA CompTox",
                            "name": match.get("preferredName", search_term),
                            "iupac_name": match.get("iupacName", ""),
                            "cid": "",
                            "synonyms": [],
                            "cas": match.get("casrn", req.cas_number or ""), 
                            "smiles": match.get("smiles", ""),
                            "inchi": match.get("inchiString", ""),
                            "inchikey": match.get("inchiKey", ""),
                            "formula": match.get("formula", ""),
                            "molecular_weight": float(match.get("mass", 0.0))
                        }
        except:
            pass

    # OPSIN is the final network fallback for a complete systematic name. It
    # is intentionally not attempted for generic or polymer-like annotations.
    eligible, _reason = opsin_name_is_eligible(req.compound_name)
    if eligible:
        try:
            url = f"{OPSIN_WEB_SERVICE}/{quote(req.compound_name.strip(), safe='')}.json"
            async with session.get(url, ssl=False, timeout=5) as opsin_res:
                if opsin_res.status == 200:
                    parsed = opsin_identity_from_payload(req.compound_name.strip(), await opsin_res.json())
                    if parsed is not None:
                        parsed["cas"] = req.cas_number or ""
                        return parsed
        except Exception:
            pass

    return {
        "status": "Unresolved",
        "source": "None",
        "name": req.compound_name or search_term, 
        "iupac_name": "",
        "cid": "",
        "synonyms": [], 
        "cas": req.cas_number or "", 
        "smiles": "", 
        "inchi": "", 
        "inchikey": "", 
        "formula": "", 
        "molecular_weight": 0.0
    }

async def analog_search_and_score_async(target_smiles, target_fp, target_parent_moiety, target_fgs, db_session, session: aiohttp.ClientSession, rate_limiter: RateLimiter = None):
    analogs = []
    if not target_smiles or not target_fp:
        return analogs
        
    from urllib.parse import quote
    from rdkit import DataStructs
    from services import generate_fingerprint, detect_parent_moiety, detect_functional_groups, get_structure_image_base64
    
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/{quote(target_smiles)}/property/Title,CanonicalSMILES,IsomericSMILES/JSON?Threshold=60&MaxRecords=10"
    
    for attempt in range(3):
        try:
            if rate_limiter: await rate_limiter.wait()
            async with session.get(url, ssl=False, timeout=5) as res:
                if res.status == 503 or res.status == 504:
                    await asyncio.sleep(1.5 * (attempt + 1))
                    continue
                if res.status != 200:
                    break
                data = await res.json()
                if 'PropertyTable' not in data:
                    break
                    
                props = data['PropertyTable']['Properties']
                for prop in props:
                    analog_smiles = prop.get('CanonicalSMILES', prop.get('IsomericSMILES', prop.get('SMILES', '')))
                    if not analog_smiles or analog_smiles == target_smiles:
                        continue
                        
                    analog_fp = generate_fingerprint(analog_smiles)
                    if not analog_fp: continue
                        
                    tanimoto = DataStructs.TanimotoSimilarity(target_fp, analog_fp)
                    
                    classification = "Weak Analog"
                    if tanimoto >= 0.85: classification = "Very Strong Analog"
                    elif tanimoto >= 0.70: classification = "Strong Analog"
                    elif tanimoto >= 0.50: classification = "Moderate Analog"
                        
                    analogs.append({
                        "name": prop.get('Title', 'Unknown'),
                        "smiles": analog_smiles,
                        "similarity_score": round(tanimoto * 100, 1),
                        "confidence_tier": classification,
                        "image_base64": get_structure_image_base64(analog_smiles)
                    })
                break
        except:
            break
            
    analogs.sort(key=lambda x: x['similarity_score'], reverse=True)
    return analogs[:10]
