import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
import os
import re
import requests
import json
import csv
import io
import subprocess
import tempfile
import shutil
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from fastapi import HTTPException
from schemas import AssessmentRequest, AssessmentResponse, FunctionalGroup, AnalogInfo
from models import AssessmentRecord
from opsin_identity import OPSIN_WEB_SERVICE, opsin_identity_from_payload, opsin_name_is_eligible

_NO_PROXY_SESSION = requests.Session()
_NO_PROXY_SESSION.trust_env = False

def _http_get(*args, **kwargs):
    kwargs.setdefault("verify", False)
    return _NO_PROXY_SESSION.get(*args, **kwargs)

def _http_post(*args, **kwargs):
    kwargs.setdefault("verify", False)
    return _NO_PROXY_SESSION.post(*args, **kwargs)


def _resolve_with_opsin(name: object, cas_number: object = "") -> dict | None:
    """Use OPSIN only for an eligible systematic-name fallback.

    OPSIN parses nomenclature; it is not an identity registry. Returned
    structures therefore remain explicitly review-recommended.
    """
    eligible, _reason = opsin_name_is_eligible(name)
    if not eligible:
        return None
    try:
        from urllib.parse import quote

        url = f"{OPSIN_WEB_SERVICE}/{quote(str(name).strip(), safe='')}.json"
        response = _http_get(url, timeout=5)
        if response.status_code != 200:
            return None
        result = opsin_identity_from_payload(str(name).strip(), response.json())
        if result is not None:
            result["cas"] = str(cas_number or "").strip()
        return result
    except Exception:
        return None


def _parse_first_uri_list_line(text: str) -> str:
    for line in (text or "").splitlines():
        line = line.strip()
        if line.startswith("http://") or line.startswith("https://"):
            return line
    return ""


def _extract_text_result(text: str, ttc_unit: str = "µg/kg bw/day") -> dict:
    text = text or ""
    class_match = re.search(r"\bClass\s*(I{1,3}|IV|V)\b", text, re.IGNORECASE)
    ttc_match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*µg\s*/?\s*kg\s*bw\s*/?\s*day", text, re.IGNORECASE)
    if not ttc_match:
        ttc_match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*ug\s*/?\s*kg\s*bw\s*/?\s*day", text, re.IGNORECASE)
    return {
        "class": f"Class {class_match.group(1).upper()}" if class_match else "",
        "ttc_value": float(ttc_match.group(1)) if ttc_match else None,
        "ttc_unit": ttc_unit if ttc_match else "",
    }


def _sdf_from_smiles(smiles: str, name: str = "Query") -> bytes:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES for Toxtree dataset generation")
    mol.SetProp("_Name", name or "Query")
    mol.SetProp("SMILES", smiles)
    block = Chem.MolToMolBlock(mol)
    sdf = f"{block}\n>  <SMILES>\n{smiles}\n\n$$$$\n"
    return sdf.encode("utf-8")


def _csv_from_smiles(smiles: str, name: str = "Query") -> str:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=["NAME", "SMILES"])
    writer.writeheader()
    writer.writerow({"NAME": name or "Query", "SMILES": smiles})
    return buffer.getvalue()


def _detect_java_executable() -> str:
    java_path = os.getenv("TOXTREE_JAVA_PATH", "").strip()
    if java_path:
        return java_path
    return shutil.which("java") or ""


def _parse_toxtree_output_file(output_path: str, variant_key: str) -> dict:
    if not os.path.exists(output_path):
        return {"status": "Unavailable", "message": "Toxtree output file not found"}

    try:
        with open(output_path, "r", encoding="utf-8", errors="ignore") as fh:
            content = fh.read()
    except Exception as exc:
        return {"status": "Unavailable", "message": f"Unable to read Toxtree output: {exc}"}

    parsed_csv = None
    try:
        parsed_csv = list(csv.DictReader(io.StringIO(content)))
    except Exception:
        parsed_csv = None

    row = parsed_csv[0] if parsed_csv else {}
    normalized = {str(k).strip(): v for k, v in row.items()} if row else {}
    raw_text = content[:5000]

    class_candidates = [
        "Cramer Class",
        "Cramer_Class",
        "toxTree.tree.cramer.CramerRules",
        "cramer2.CramerRulesWithExtensions",
        "toxtree.tree.cramer3.RevisedCramerDecisionTree",
    ]
    path_candidates = [
        "Cramer Tree Result",
        "CramerTreeResult",
        "toxTree.tree.cramer.CramerTreeResult",
    ]
    ttc_candidates = [
        "TTC Value",
        "TTC",
        "Kroes TTC",
        "Cramer TTC",
    ]

    def _first_nonempty(keys):
        for key in keys:
            value = normalized.get(key)
            if value not in (None, "", "nan", "NaN"):
                return value
        return ""

    cramer_class = _first_nonempty(class_candidates)
    cramer_path = _first_nonempty(path_candidates)
    ttc_value = ""
    ttc_unit = "µg/kg bw/day"
    for key in ttc_candidates:
        value = normalized.get(key)
        if value not in (None, "", "nan", "NaN"):
            ttc_value = value
            break
    if not ttc_value:
        class_from_text = re.search(r"\bClass\s*(I{1,3}|IV|V)\b", raw_text, re.IGNORECASE)
        if class_from_text:
            cramer_class = f"Class {class_from_text.group(1).upper()}"
    if not cramer_path:
        cramer_path = normalized.get("path", "") or normalized.get("Tree Result", "") or ""

    return {
        "status": "Calculated" if cramer_class or ttc_value else "Unparsed",
        "algorithm": variant_key,
        "cramer_class": cramer_class,
        "ttc_value": ttc_value,
        "ttc_unit": ttc_unit,
        "cramer_path": cramer_path,
        "raw_response": raw_text,
        "message": ""
    }


def run_local_toxtree_variant(standardized_smiles: str, variant_key: str, module_class: str, compound_name: str = "Query") -> dict:
    """
    Run Toxtree locally in headless mode using a configured JAR path.

    This avoids any external credentials. The caller must set:
    - TOXTREE_JAR_PATH
    - optionally TOXTREE_JAVA_PATH
    """
    jar_path = os.getenv("TOXTREE_JAR_PATH", "").strip()
    java_path = _detect_java_executable()
    if not jar_path:
        return {
            "status": "Unavailable",
            "algorithm": module_class,
            "message": "TOXTREE_JAR_PATH is not configured"
        }
    if not os.path.exists(jar_path):
        return {
            "status": "Unavailable",
            "algorithm": module_class,
            "message": f"Toxtree JAR not found at {jar_path}"
        }
    if not java_path:
        return {
            "status": "Unavailable",
            "algorithm": module_class,
            "message": "Java runtime not found. Set TOXTREE_JAVA_PATH or install Java on PATH."
        }

    input_csv = None
    output_csv = None
    try:
        with tempfile.TemporaryDirectory(prefix="toxtree_") as tmpdir:
            input_csv = os.path.join(tmpdir, "input.csv")
            output_csv = os.path.join(tmpdir, f"output_{variant_key}.csv")
            with open(input_csv, "w", encoding="utf-8", newline="") as fh:
                fh.write(_csv_from_smiles(standardized_smiles, compound_name))

            cmd = [
                java_path,
                "-jar",
                jar_path,
                "-n",
                "-i",
                input_csv,
                "-o",
                output_csv,
                "-m",
                module_class
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
            if proc.returncode != 0 and not os.path.exists(output_csv):
                return {
                    "status": "Unavailable",
                    "algorithm": module_class,
                    "message": (proc.stderr or proc.stdout or f"Toxtree exited with code {proc.returncode}").strip()
                }

            parsed = _parse_toxtree_output_file(output_csv, variant_key)
            parsed["algorithm_uri"] = module_class
            parsed["java"] = java_path
            parsed["jar"] = jar_path
            parsed["stdout"] = (proc.stdout or "")[:3000]
            parsed["stderr"] = (proc.stderr or "")[:3000]
            return parsed
    except subprocess.TimeoutExpired:
        return {
            "status": "Unavailable",
            "algorithm": module_class,
            "message": "Toxtree timed out"
        }
    except Exception as exc:
        return {
            "status": "Unavailable",
            "algorithm": module_class,
            "message": str(exc)
        }


def query_ambit_compound(search_term: str, base_url: str = None, auth=None) -> dict:
    """
    Query an AMBIT host using its OpenTox query endpoints.

    This is a lightweight helper for self-hosted AMBIT deployments. It is
    intentionally permissive and returns the raw response text along with the
    endpoint used so the caller can decide how to interpret the search result.
    """
    base_url = (base_url or os.getenv("TOXTREE_AMBIT_BASE_URL", "")).strip().rstrip("/")
    if not base_url:
        return {"status": "Unavailable", "message": "AMBIT base URL is not configured"}

    endpoint_candidates = [
        f"{base_url}/query/compound/search/all",
        f"{base_url}/query/compound/search/name",
        f"{base_url}/query/compound/search/cas",
    ]
    for endpoint in endpoint_candidates:
        try:
            res = _http_get(endpoint, params={"search": search_term}, headers={"Accept": "application/json, text/uri-list, text/plain"}, timeout=30, auth=auth)
            if res.ok:
                return {
                    "status": "Calculated",
                    "endpoint": endpoint,
                    "search_term": search_term,
                    "content_type": res.headers.get("Content-Type", ""),
                    "raw_response": res.text[:5000],
                }
        except Exception as exc:
            last_error = str(exc)
    return {
        "status": "Unavailable",
        "search_term": search_term,
        "message": locals().get("last_error", "No AMBIT query endpoint returned a result")
    }


def evaluate_toxtree_cramer_variants(standardized_smiles: str, compound_name: str = "Query") -> dict:
    """
    Evaluate an AMBIT/OpenTox host for the available Toxtree Cramer modules.

    The host is configurable so the app can run against either a public demo or
    a self-hosted AMBIT instance. We attempt all three variants and return an
    explicit unavailable status when the host does not expose one.
    """
    mode = os.getenv("TOXTREE_EXECUTION_MODE", "").strip().lower()
    jar_path = os.getenv("TOXTREE_JAR_PATH", "").strip()
    base_url = os.getenv("TOXTREE_AMBIT_BASE_URL", "").strip().rstrip("/")
    if not mode:
        mode = "local" if jar_path else "ambit_rest"
    if mode == "local" and not jar_path:
        mode = "ambit_rest" if base_url else "local"
    algorithm_map = {
        "original": "toxtreecramer",
        "extensions": "toxtreecramer2",
        "revised": "toxtreecramer3",
    }

    results = {
        "source": jar_path if mode == "local" else base_url,
        "mode": mode,
        "available_algorithms": [],
        "variants": {},
        "status": "Unavailable",
        "message": ""
    }

    if mode == "local":
        results["message"] = "Local Toxtree runner selected"
        variant_modules = {
            "original": "toxTree.tree.cramer.CramerRules",
            "extensions": "cramer2.CramerRulesWithExtensions",
            "revised": "toxtree.tree.cramer3.RevisedCramerDecisionTree",
        }
        if compound_name and compound_name != "Query":
            results["query"] = {"status": "Not requested"}
        else:
            results["query"] = {"status": "Not requested"}
        for variant_key, module_class in variant_modules.items():
            results["variants"][variant_key] = run_local_toxtree_variant(
                standardized_smiles=standardized_smiles,
                variant_key=variant_key,
                module_class=module_class,
                compound_name=compound_name
            )
        results["available_algorithms"] = list(variant_modules.values())
    else:
        auth_user = os.getenv("TOXTREE_BASIC_AUTH_USER", "").strip()
        auth_pass = os.getenv("TOXTREE_BASIC_AUTH_PASS", "").strip()
        auth = (auth_user, auth_pass) if auth_user and auth_pass else None

        if compound_name and compound_name != "Query":
            results["query"] = query_ambit_compound(compound_name, base_url=base_url, auth=auth)
        else:
            results["query"] = {"status": "Not requested"}

        try:
            alg_res = _http_get(f"{base_url}/algorithm", headers={"Accept": "text/uri-list"}, timeout=20, auth=auth)
            if alg_res.ok and alg_res.text:
                available_algorithms = alg_res.text or ""
                results["available_algorithms"] = available_algorithms.splitlines()
            else:
                results["available_algorithms"] = [
                    f"{base_url}/algorithm/toxtreecramer",
                    f"{base_url}/algorithm/toxtreecramer2",
                    f"{base_url}/algorithm/toxtreecramer3",
                ]
                results["message"] = f"Algorithm registry unavailable or empty ({getattr(alg_res, 'status_code', 'n/a')})"
        except Exception as exc:
            results["available_algorithms"] = [
                f"{base_url}/algorithm/toxtreecramer",
                f"{base_url}/algorithm/toxtreecramer2",
                f"{base_url}/algorithm/toxtreecramer3",
            ]
            results["message"] = f"Unable to retrieve algorithm registry: {exc}"

        try:
            sdf_bytes = _sdf_from_smiles(standardized_smiles, compound_name)
            dataset_res = _http_post(
                f"{base_url}/dataset",
                data=sdf_bytes,
                headers={
                    "Content-Type": "chemical/x-mdl-sdfile",
                    "Accept": "text/uri-list"
                },
                timeout=30,
                auth=auth
            )
            dataset_uri = dataset_res.headers.get("Location", "").strip() or _parse_first_uri_list_line(dataset_res.text)
            if not dataset_uri and dataset_res.ok and dataset_res.text.strip().startswith("http"):
                dataset_uri = dataset_res.text.strip().splitlines()[0].strip()
            if not dataset_uri:
                if dataset_res.status_code == 403:
                    results["message"] = (
                        "Dataset upload was rejected by the AMBIT host (403). "
                        "This host appears to require authentication for POST operations."
                    )
                else:
                    results["message"] = f"Dataset upload did not return a dataset URI ({dataset_res.status_code})"
                return results
        except Exception as exc:
            results["message"] = f"Dataset upload failed: {exc}"
            return results

        results["dataset_uri"] = dataset_uri

        def _run_variant(variant_key: str, algorithm_id: str):
            try:
                res = _http_post(
                    f"{base_url}/algorithm/{algorithm_id}",
                    data={"dataset_uri": dataset_uri},
                    headers={"Accept": "text/uri-list"},
                    timeout=60,
                    auth=auth
                )
                candidate_uri = res.headers.get("Location", "").strip() or _parse_first_uri_list_line(res.text)
                raw_text = res.text or ""
                if not candidate_uri:
                    candidate_uri = _parse_first_uri_list_line(raw_text)
                if candidate_uri and "/model/" in candidate_uri:
                    pred_res = _http_post(
                        candidate_uri,
                        data={"dataset_uri": dataset_uri},
                        headers={"Accept": "text/csv, text/plain, application/json, text/uri-list"},
                        timeout=60,
                        auth=auth
                    )
                    raw_text = pred_res.text or raw_text
                    if not candidate_uri and pred_res.headers.get("Location"):
                        candidate_uri = pred_res.headers.get("Location", "").strip()
                parsed = _extract_text_result(raw_text)
                algorithm_uri = f"{base_url}/algorithm/{algorithm_id}"
                variant_status = "Calculated" if parsed.get("class") or parsed.get("ttc_value") is not None else "Available"
                return {
                    "status": variant_status,
                    "algorithm": algorithm_id,
                    "algorithm_uri": algorithm_uri,
                    "prediction_uri": candidate_uri,
                    "cramer_class": parsed.get("class") or "",
                    "ttc_value": parsed.get("ttc_value"),
                    "ttc_unit": parsed.get("ttc_unit") or "µg/kg bw/day",
                    "raw_response": raw_text[:5000],
                    "message": ""
                }
            except Exception as exc:
                return {
                    "status": "Unavailable",
                    "algorithm": algorithm_id,
                    "algorithm_uri": f"{base_url}/algorithm/{algorithm_id}",
                    "message": str(exc)
                }

        for variant_key, algorithm_id in algorithm_map.items():
            results["variants"][variant_key] = _run_variant(variant_key, algorithm_id)

    primary = results["variants"].get("original") or {}
    if primary.get("status") in ("Calculated", "Unparsed") and (primary.get("cramer_class") or primary.get("ttc_value") is not None):
        results["status"] = "Calculated"
    elif results["variants"].get("extensions", {}).get("status") in ("Calculated", "Unparsed"):
        results["status"] = "Calculated"
    elif results["variants"].get("revised", {}).get("status") in ("Calculated", "Unparsed"):
        results["status"] = "Calculated"

    return results

def validate_input(req: AssessmentRequest):
    """ STEP 1: Input Validation """
    direct_identity_supplied = bool(req.original_smiles or req.standardized_smiles)

    if not req.compound_name and not req.cas_number and not direct_identity_supplied:
        raise HTTPException(status_code=400, detail="Must provide either compound_name or cas_number")
    
    if req.cas_number:
        # Basic CAS validation (e.g., 50-00-0)
        if not re.match(r'^\d{2,7}-\d{2}-\d$', req.cas_number):
            raise HTTPException(status_code=400, detail="Invalid CAS Number format")

    if req.molecular_formula:
        # Simple Empirical Formula regex
        if not re.match(r'^[A-Z][a-z]?\d*([A-Z][a-z]?\d*)*$', req.molecular_formula):
            raise HTTPException(status_code=400, detail="Invalid Molecular Formula format")

    return "Valid"

def resolve_identity(req: AssessmentRequest, db=None):
    """ STEP 2: Chemical Identity Resolution """
    search_term = req.cas_number if req.cas_number else req.compound_name
    
    # Fast-fail for complex polymers or obvious fragments that won't resolve in simple APIs
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

    # Try both CAS and Name as search terms
    search_terms = []
    if req.cas_number:
        search_terms.append((req.cas_number.strip(), "name"))  # CAS often resolves via name endpoint
    if req.compound_name:
        search_terms.append((req.compound_name.strip(), "name"))
        
    if not search_terms:
        search_terms.append((search_term, "name"))

    # Attempt 1: PubChem PUG REST API
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
    
    for term, search_type in search_terms:
        try:
            pubchem_url = f"{base_url}/{search_type}/{term}/property/Title,CanonicalSMILES,IsomericSMILES,InChI,InChIKey,MolecularFormula,MolecularWeight,IUPACName/JSON"
            res = _http_get(pubchem_url, timeout=10) # Increased timeout to 10s
            if res.status_code == 200:
                data = res.json()
                props = data['PropertyTable']['Properties'][0]
                
                cid = str(props.get('CID', ''))
                title = props.get('Title', '')
                iupac = props.get('IUPACName', '')
                
                # Use robust fallback for missing or null properties
                smiles = props.get('CanonicalSMILES') or props.get('IsomericSMILES') or props.get('SMILES') or props.get('ConnectivitySMILES') or ''
                
                # Fallback 1: If PubChem succeeds but has no SMILES, try NCI/Cactus
                if not smiles:
                    try:
                        from urllib.parse import quote
                        cactus_url = f"https://cactus.nci.nih.gov/chemical/structure/{quote(term)}/smiles"
                        c_res = _http_get(cactus_url, timeout=3)
                        if c_res.status_code == 200 and c_res.text:
                            smiles = c_res.text.strip().split('\n')[0]
                    except:
                        pass

                inchi = props.get('InChI') or ''
                inchikey = props.get('InChIKey') or ''
                formula = props.get('MolecularFormula') or ''
                mw = float(props.get('MolecularWeight') or 0.0)
                
                synonyms = []
                syn_url = f"{base_url}/{search_type}/{term}/synonyms/JSON"
                syn_res = _http_get(syn_url, timeout=5)
                if syn_res.status_code == 200:
                    syn_data = syn_res.json()
                    if 'InformationList' in syn_data:
                        synonyms = syn_data['InformationList']['Information'][0].get('Synonym', [])[:5]
                
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
        except Exception as e:
            pass

    # Attempt 2: NCI/CADD Cactus API (Restored with strict timeout)
    from urllib.parse import quote
    for term, search_type in search_terms:
        try:
            cactus_url = f"https://cactus.nci.nih.gov/chemical/structure/{quote(term)}/smiles"
            c_res = _http_get(cactus_url, timeout=4)
            if c_res.status_code == 200 and c_res.text:
                smiles = c_res.text.strip().split('\n')[0]
                
                # We can also attempt to get IUPAC name and Formula from Cactus if desired,
                # but to save time, we return just the SMILES to proceed to cheminformatics steps.
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
    # Attempt 3: EPA CompTox API
    import os
    epa_api_key = os.getenv("EPA_API_KEY")
    if epa_api_key:
        try:
            import requests
            epa_url = f"https://comptox.epa.gov/ctx-api/chemical/search/equal/{requests.utils.quote(search_term, safe='')}"
            headers = {"x-api-key": epa_api_key, "Accept": "application/json"}
            epa_res = _http_get(epa_url, headers=headers, timeout=2)
            if epa_res.status_code == 200:
                epa_data = epa_res.json()
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

    # Attempt 4: Local prior assessment fallback
    if db is not None:
        try:
            q = db.query(AssessmentRecord)
            local_match = None

            if req.cas_number:
                local_match = (
                    q.filter(AssessmentRecord.input_cas_number == req.cas_number.strip())
                     .order_by(AssessmentRecord.created_at.desc())
                     .first()
                )

            if local_match is None and req.compound_name:
                local_match = (
                    q.filter(
                        AssessmentRecord.input_compound_name.isnot(None),
                        AssessmentRecord.input_compound_name.ilike(req.compound_name.strip())
                    )
                    .order_by(AssessmentRecord.created_at.desc())
                    .first()
                )

            if local_match and (local_match.resolved_smiles or local_match.standardized_smiles):
                return {
                    "status": local_match.resolution_status or "Exact Match",
                    "source": f"Local Cache ({local_match.source_database or 'Prior Assessment'})",
                    "name": local_match.resolved_name or req.compound_name or search_term,
                    "iupac_name": "",
                    "cid": "",
                    "synonyms": [],
                    "cas": local_match.resolved_cas or req.cas_number or "",
                    "smiles": local_match.resolved_smiles or local_match.standardized_smiles or "",
                    "inchi": local_match.resolved_inchi or "",
                    "inchikey": local_match.resolved_inchikey or "",
                    "formula": local_match.resolved_formula or req.molecular_formula or "",
                    "molecular_weight": float(local_match.resolved_mw or 0.0),
                }
        except Exception:
            pass

    # Attempt 5: OPSIN systematic-name parser. This runs only after registry
    # and local-cache resolution fail, and bypasses generic analytical names.
    opsin_identity = _resolve_with_opsin(req.compound_name, req.cas_number)
    if opsin_identity is not None:
        return opsin_identity

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

def standardize_structure(smiles: str):
    """ STEP 3: Structure Standardization (Layer 2) """
    if not smiles:
        return ""
    
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return ""
    
    try:
        from rdkit.Chem.MolStandardize import rdMolStandardize
        
        # 1. Disconnect organometallics
        mol = rdMolStandardize.MetalDisconnector().Disconnect(mol)
        
        # 2. Normalize (correct drawing errors, standard functional groups)
        mol = rdMolStandardize.Normalize(mol)
        
        # 3. Reionize (standardize charges)
        mol = rdMolStandardize.Reionize(mol)
        
        # 4. Remove salts (largest fragment)
        remover = rdMolStandardize.LargestFragmentChooser()
        mol = remover.choose(mol)
        
        # 5. Charge Neutralization
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)
        
        # 6. Tautomer Normalization
        te = rdMolStandardize.TautomerEnumerator()
        mol = te.Canonicalize(mol)
    except Exception as e:
        # Fallback to basic if MolStandardize fails
        frags = Chem.GetMolFrags(mol, asMols=True)
        if frags:
            mol = max(frags, default=mol, key=lambda m: m.GetNumAtoms())
            
    # Canonicalize
    return Chem.MolToSmiles(mol, canonical=True)

def get_structure_image_base64(smiles: str) -> str:
    """ Generate a Base64 PNG image of the structure from SMILES using RDKit """
    if not smiles: return ""
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return ""
    
    from rdkit.Chem import Draw
    import base64
    from io import BytesIO
    
    # White background with explicit size
    img = Draw.MolToImage(mol, size=(300, 300), bg_color=(255, 255, 255))
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    return base64.b64encode(buffered.getvalue()).decode("utf-8")

# Common functional group SMARTS
FG_PATTERNS = {
    "Hydroxyl": "[OX2H]",
    "Phenol": "[OX2H][c]",
    "Alcohol": "[OX2H][C;!$(C=O)]",
    "Ether": "[OD2]([#6])[#6]",
    "Aldehyde": "[CX3H1](=O)[#6]",
    "Ketone": "[#6][CX3](=O)[#6]",
    "Carboxylic Acid": "[CX3](=O)[OX2H1]",
    "Ester": "[#6][CX3](=O)[OX2H0][#6]",
    "Amide": "[NX3][CX3](=[OX1])[#6]",
    "Anhydride": "[CX3](=[OX1])[OX2][CX3](=[OX1])",
    "Amine": "[NX3;H2,H1,H0;!$(NC=O)]",
    "Nitro": "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]",
    "Nitrile": "[CX2]#[NX1]",
    "Sulfoxide": "[$([#16X3](=[OX1])([#6])[#6])]",
    "Sulfone": "[$([#16X4](=[OX1])(=[OX1])([#6])[#6])]",
    "Phosphate": "[PX4](=O)([OX2H,OX2H0])([OX2H,OX2H0])[OX2H,OX2H0]",
    "Siloxane": "[Si][OX2][Si]",
    "Halogen": "[#9,#17,#35,#53]",
    "Aromatic Ring": "a",
    "Aliphatic Chain": "[CX4,CX3!$(C=[O,N,S]),CX2!$(C#[N])]"
}

def detect_functional_groups(smiles: str):
    """ STEP 4: Functional Group Detection """
    mol = Chem.MolFromSmiles(smiles)
    fgs = []
    if not mol:
        return fgs
        
    for fg_name, smarts in FG_PATTERNS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern:
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                # Add each match or group them. Here we group by functional group type
                atoms_matched = list(set([atom for match in matches for atom in match]))
                fgs.append({
                    "functional_group": fg_name,
                    "matched_atoms": atoms_matched,
                    "SMARTS_pattern": smarts
                })
    return fgs

def detect_parent_moiety(smiles: str):
    """ STEP 5: Parent Moiety Detection """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "Unknown"
        
    # Use Murcko Scaffold as the parent moiety
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    scaffold_smiles = Chem.MolToSmiles(scaffold)
    
    if not scaffold_smiles:
        # Check if it's aliphatic chain
        if all(atom.GetIsAromatic() == False for atom in mol.GetAtoms()):
            return "Aliphatic Chain"
        return "Unknown Scaffold"
        
    # Heuristic mapping for common scaffolds
    if scaffold_smiles == "c1ccccc1":
        return "Benzene"
    
    return scaffold_smiles # Fallback to SMILES of the scaffold

def assign_chemical_class(parent_moiety: str, fgs: list):
    """ STEP 6: Chemical Class Assignment """
    fg_names = [fg['functional_group'] for fg in fgs]
    
    is_aromatic = "Aromatic Ring" in fg_names
    
    if "Carboxylic Acid" in fg_names:
        return "Aromatic Acid" if is_aromatic else "Aliphatic Acid"
    elif "Ester" in fg_names:
        return "Aromatic Ester" if is_aromatic else "Aliphatic Ester"
    elif "Alcohol" in fg_names:
        return "Aromatic Alcohol" if is_aromatic else "Aliphatic Alcohol"
    elif "Amine" in fg_names:
        return "Aromatic Amine" if is_aromatic else "Aliphatic Amine"
    elif "Siloxane" in fg_names:
        return "Siloxane Compound"
    elif "Phenol" in fg_names:
        return "Phenolic Compound"
    
    if is_aromatic:
        return "Aromatic Compound (Miscellaneous)"
    return "Aliphatic Compound (Miscellaneous)"

def generate_fingerprint(smiles: str):
    """ STEP 7: Molecular Fingerprint Generation """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    # Morgan Fingerprint Radius = 2 (ECFP4)
    from rdkit.Chem import rdFingerprintGenerator
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp = mfpgen.GetFingerprint(mol)
    return fp

def analog_search_and_score(target_smiles, target_fp, target_parent_moiety, target_fgs, db_session):
    """ STEPS 8, 9, 10: Analog Search, Tanimoto, Interpretation """
    analogs = []
    if not target_smiles or not target_fp:
        return analogs
        
    import requests
    from urllib.parse import quote
    
    # STEP 8: Search similarity database (Using PubChem Fast Similarity Search)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/{quote(target_smiles)}/property/Title,CanonicalSMILES,IsomericSMILES/JSON?Threshold=60&MaxRecords=10"
    try:
        res = _http_get(url, timeout=3)
    except Exception:
        return analogs
    
    if res.status_code != 200:
        return analogs
        
    data = res.json()
    if 'PropertyTable' not in data:
        return analogs
        
    props = data['PropertyTable']['Properties']
    target_fg_names = set(fg['functional_group'] for fg in target_fgs)
    
    for prop in props:
        analog_smiles = prop.get('CanonicalSMILES', prop.get('IsomericSMILES', prop.get('SMILES', '')))
        if not analog_smiles or analog_smiles == target_smiles:
            continue
            
        analog_fp = generate_fingerprint(analog_smiles)
        if not analog_fp: continue
            
        # STEP 9: Tanimoto Similarity Calculation
        tanimoto = DataStructs.TanimotoSimilarity(target_fp, analog_fp)
        
        # Determine analog properties for Step 10
        analog_parent = detect_parent_moiety(analog_smiles)
        analog_fgs = detect_functional_groups(analog_smiles)
        analog_fg_names = set(fg['functional_group'] for fg in analog_fgs)
        
        # STEP 10: Similarity Interpretation
        classification = "Weak Analog"
        if tanimoto >= 0.85:
            classification = "Very Strong Analog"
        elif tanimoto >= 0.70:
            classification = "Strong Analog"
        elif tanimoto >= 0.50:
            classification = "Moderate Analog"
            
        analogs.append({
            "name": prop.get('Title', 'Unknown'),
            "smiles": analog_smiles,
            "similarity_score": round(tanimoto * 100, 1),
            "confidence_tier": classification,
            "image_base64": get_structure_image_base64(analog_smiles)
        })
        
    # Sort descending by Tanimoto
    analogs.sort(key=lambda x: x['similarity_score'], reverse=True)
    return analogs[:10]

def determine_confidence(identity_status, standardized_smiles, parent_moiety, fgs):
    """ STEP 11: Confidence Scoring """
    if identity_status == "Exact Match" and standardized_smiles and parent_moiety != "Unknown" and len(fgs) > 0:
        return "High"
    elif identity_status == "Unresolved" or not standardized_smiles:
        return "Low"
    else:
        return "Medium"

def generate_rationale(identity_status, chemical_class, tanimoto_highest=0.0):
    """ Generate a textual rationale for the grouping/classification """
    rationale = f"Classification assigned as {chemical_class} based on the dominant functional group and scaffold context. "
    rationale += f"Identity resolution was {identity_status}. "
    return rationale


import pandas as pd
import time
def process_batch_dataframe(df: pd.DataFrame, db_session) -> pd.DataFrame:
    # Iterate rows and apply assess_chemical logic
    pass
