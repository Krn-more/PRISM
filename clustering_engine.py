import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem, MACCSkeys, rdMolDescriptors, rdFingerprintGenerator
from rdkit import DataStructs
from rdkit.Avalon import pyAvalonTools
from rdkit.ML.Cluster import Butina
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.cluster import DBSCAN, AgglomerativeClustering, SpectralClustering
from sklearn.metrics import silhouette_score
import uuid
import warnings
import re
warnings.filterwarnings('ignore')

import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from clustering_policy import load_clustering_policy
from clustering_class_layer import (
    annotate_cluster_class_evidence,
    build_class_distance_matrix,
    load_class_layer_rules,
)
from clustering_descriptor_layer import (
    annotate_descriptor_layer_evidence,
    build_descriptor_distance_matrix,
    load_descriptor_layer_rules,
)
from clustering_topology_layer import (
    annotate_topology_layer_evidence,
    build_topology_distance_matrix,
    load_topology_layer_rules,
)
from clustering_ionisation_layer import (
    annotate_ionisation_layer_evidence,
    build_ionisation_distance_matrix,
    load_ionisation_layer_rules,
)
from clustering_toxicophore_layer import (
    annotate_toxicophore_layer_evidence,
    build_toxicophore_distance_matrix,
    load_toxicophore_layer_rules,
)
from clustering_chemont_layer import (
    annotate_chemont_layer_evidence,
    build_chemont_distance_matrix,
    load_chemont_layer_rules,
)


def _chunk_worker_count(chunks) -> int:
    """Use bounded in-process workers so Phase B is safe in the frozen Windows app."""
    return max(1, min(4, len(chunks)))

def _trace_phase_b(trace: list[dict] | None, step: str, status: str, df: pd.DataFrame | None = None, note: str | None = None):
    if trace is None:
        return
    payload = {"step": step, "status": status}
    if df is not None:
        payload["rows"] = int(len(df))
        payload["columns"] = int(len(df.columns))
        payload["has_cluster_id"] = "Cluster ID" in df.columns
    if note:
        payload["note"] = note
    trace.append(payload)
    print(f"[Phase B] {step} [{status}]"
          f"{' rows=' + str(payload['rows']) if 'rows' in payload else ''}"
          f"{' cols=' + str(payload['columns']) if 'columns' in payload else ''}"
          f"{' | ' + note if note else ''}", flush=True)

def _process_chunk_adme(chunk_smiles):
    mw_list, logp_list, tpsa_list, min_estate_list, max_estate_list = [], [], [], [], []
    asph_list, pmi1_list, pmi2_list, pmi3_list, rad_gyr_list = [], [], [], [], []
    from rdkit.Chem import Descriptors, AllChem, Descriptors3D
    
    for smiles in chunk_smiles:
        mol = Chem.MolFromSmiles(str(smiles)) if pd.notnull(smiles) and str(smiles).strip() else None
        if mol:
            mw_list.append(round(Descriptors.ExactMolWt(mol), 2))
            logp_list.append(round(Descriptors.MolLogP(mol), 2))
            tpsa_list.append(round(Descriptors.TPSA(mol), 2))
            min_estate_list.append(round(Descriptors.MinEStateIndex(mol), 2))
            max_estate_list.append(round(Descriptors.MaxEStateIndex(mol), 2))
            try:
                mol3d = Chem.AddHs(mol)
                ret = AllChem.EmbedMolecule(mol3d, maxAttempts=10, randomSeed=42)
                if ret != -1:
                    asph_list.append(round(Descriptors3D.Asphericity(mol3d), 3))
                    pmi1_list.append(round(Descriptors3D.PMI1(mol3d), 3))
                    pmi2_list.append(round(Descriptors3D.PMI2(mol3d), 3))
                    pmi3_list.append(round(Descriptors3D.PMI3(mol3d), 3))
                    rad_gyr_list.append(round(Descriptors3D.RadiusOfGyration(mol3d), 3))
                else:
                    raise ValueError("Embedding failed")
            except:
                asph_list.append(0.0)
                pmi1_list.append(0.0)
                pmi2_list.append(0.0)
                pmi3_list.append(0.0)
                rad_gyr_list.append(0.0)
        else:
            mw_list.append(None); logp_list.append(None); tpsa_list.append(None)
            min_estate_list.append(None); max_estate_list.append(None)
            asph_list.append(None); pmi1_list.append(None); pmi2_list.append(None)
            pmi3_list.append(None); rad_gyr_list.append(None)
            
    return mw_list, logp_list, tpsa_list, min_estate_list, max_estate_list, asph_list, pmi1_list, pmi2_list, pmi3_list, rad_gyr_list

def calculate_adme_properties(df: pd.DataFrame, smiles_col: str = 'Standardized SMILES'):
    smiles_list = df[smiles_col].tolist()
    cores = os.cpu_count() or 4
    chunk_size = max(1, len(smiles_list) // cores)
    chunks = [smiles_list[i:i + chunk_size] for i in range(0, len(smiles_list), chunk_size)]
    
    results = [None] * len(chunks)
    # Process workers can fail to start from the one-file PyInstaller executable
    # and turn a recoverable Phase B calculation into HTTP 500. RDKit operations
    # are safe to run in these bounded in-process workers.
    with ThreadPoolExecutor(max_workers=_chunk_worker_count(chunks)) as executor:
        futures = {executor.submit(_process_chunk_adme, chunk): i for i, chunk in enumerate(chunks)}
        for future in as_completed(futures):
            results[futures[future]] = future.result()
            
    mw, logp, tpsa, min_e, max_e, asph, pmi1, pmi2, pmi3, rgyr = [], [], [], [], [], [], [], [], [], []
    for r in results:
        mw.extend(r[0]); logp.extend(r[1]); tpsa.extend(r[2]); min_e.extend(r[3]); max_e.extend(r[4])
        asph.extend(r[5]); pmi1.extend(r[6]); pmi2.extend(r[7]); pmi3.extend(r[8]); rgyr.extend(r[9])
            
    df['MW'] = mw
    df['LogP'] = logp
    df['TPSA'] = tpsa
    df['Min_EState'] = min_e
    df['Max_EState'] = max_e
    
    df['3D_Asphericity'] = asph
    df['3D_PMI1'] = pmi1
    df['3D_PMI2'] = pmi2
    df['3D_PMI3'] = pmi3
    df['3D_RadiusOfGyration'] = rgyr
    
    return df

def flag_alerts(df: pd.DataFrame, smiles_col: str = 'Standardized SMILES'):
    # Keep this legacy helper aligned with structural_evidence so that the
    # single alert and the complete alert profile cannot disagree.
    from structural_evidence import _retain_alert
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    catalog = FilterCatalog(params)
    
    alerts_list = []
    for smiles in df[smiles_col]:
        mol = Chem.MolFromSmiles(str(smiles)) if pd.notnull(smiles) and str(smiles).strip() else None
        if mol:
            matches = [match for match in catalog.GetMatches(mol) if _retain_alert(match.GetDescription(), mol)]
            alerts_list.append(matches[0].GetDescription() if matches else "None")
        else:
            alerts_list.append("N/A")
            
    df['Alerts'] = alerts_list
    return df

def extract_scaffolds(df: pd.DataFrame, smiles_col: str = 'Standardized SMILES'):
    scaffolds = []
    for smiles in df[smiles_col]:
        mol = Chem.MolFromSmiles(str(smiles)) if pd.notnull(smiles) and str(smiles).strip() else None
        if mol:
            try:
                core = MurckoScaffold.GetScaffoldForMol(mol)
                scaffold_smiles = Chem.MolToSmiles(core)
                scaffolds.append(scaffold_smiles if scaffold_smiles else "No_Scaffold")
            except:
                scaffolds.append("Error")
        else:
            scaffolds.append("Invalid")
    df['Scaffold'] = scaffolds
    return df

def _process_chunk_fp(chunk_smiles):
    from rdkit.Chem import rdFingerprintGenerator
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    ap_gen = rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048)
    tt_gen = rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=2048)
    
    fmorgan, fmaccs, fap, ftopo, fav = [], [], [], [], []
    for smi in chunk_smiles:
        mol = Chem.MolFromSmiles(str(smi)) if pd.notnull(smi) and str(smi).strip() else None
        if mol:
            fmorgan.append(morgan_gen.GetFingerprintAsNumPy(mol))
            fmaccs.append(np.array(MACCSkeys.GenMACCSKeys(mol)))
            fap.append(ap_gen.GetFingerprintAsNumPy(mol))
            ftopo.append(tt_gen.GetFingerprintAsNumPy(mol))
            av_fp = pyAvalonTools.GetAvalonFP(mol)
            arr_av = np.zeros((1,), dtype=int)
            DataStructs.ConvertToNumpyArray(av_fp, arr_av)
            fav.append(arr_av)
        else:
            fmorgan.append(np.zeros(1024))
            fmaccs.append(np.zeros(167))
            fap.append(np.zeros(2048))
            ftopo.append(np.zeros(2048))
            fav.append(np.zeros(512))
    return fmorgan, fmaccs, fap, ftopo, fav

def compute_distance_matrix(
    df: pd.DataFrame,
    smiles_col: str,
    class_layer_rules: dict | None = None,
    descriptor_layer_rules: dict | None = None,
    topology_layer_rules: dict | None = None,
    ionisation_layer_rules: dict | None = None,
    toxicophore_layer_rules: dict | None = None,
    chemont_layer_rules: dict | None = None,
):
    n = len(df)
    if n < 2:
        # Keep the return contract identical to the multi-row path so callers
        # can safely unpack metadata even for a one-row smoke test.
        return np.zeros((n, n)), np.zeros((n, 10)), np.array([], dtype=object), {
            "base_weight": 1.0,
            "class_layer_weight": 0.0,
            "descriptor_layer_weight": 0.0,
            "topology_layer_weight": 0.0,
            "ionisation_layer_weight": 0.0,
            "toxicophore_layer_weight": 0.0,
            "chemont_layer_weight": 0.0,
            "class_layer_status": "Not assessable",
            "descriptor_layer_status": "Not assessable",
            "topology_layer_status": "Not assessable",
            "ionisation_layer_status": "Not assessable",
            "toxicophore_layer_status": "Not assessable",
            "chemont_layer_status": "Not assessable",
        }

    def _series_from_column(column_name: str, default: object) -> pd.Series:
        if column_name not in df.columns:
            return pd.Series([default] * n, index=df.index)
        raw = df.loc[:, column_name]
        if isinstance(raw, pd.DataFrame):
            raw = raw.iloc[:, 0]
        if not isinstance(raw, pd.Series):
            raw = pd.Series(raw, index=df.index)
        return raw.reindex(df.index, fill_value=default)

    def _numeric_series(column_name: str, default: float = 0.0) -> np.ndarray:
        series = _series_from_column(column_name, default)
        series = pd.to_numeric(series, errors="coerce").fillna(default)
        return series.to_numpy(dtype=float)

    scaffolds = _series_from_column('Scaffold', "").astype(str).fillna("").values
    logp = _numeric_series('LogP')
    mw = _numeric_series('MW')
    tpsa = _numeric_series('TPSA')
    min_es = _numeric_series('Min_EState')
    max_es = _numeric_series('Max_EState')
    asph = _numeric_series('3D_Asphericity')
    pmi1 = _numeric_series('3D_PMI1')
    pmi2 = _numeric_series('3D_PMI2')
    pmi3 = _numeric_series('3D_PMI3')
    rad_gyr = _numeric_series('3D_RadiusOfGyration')
    alerts = _series_from_column('Alerts', "").astype(str).fillna("").values
    
    smiles_list = df[smiles_col].tolist()
    cores = os.cpu_count() or 4
    chunk_size = max(1, len(smiles_list) // cores)
    chunks = [smiles_list[i:i + chunk_size] for i in range(0, len(smiles_list), chunk_size)]
    
    results = [None] * len(chunks)
    with ThreadPoolExecutor(max_workers=_chunk_worker_count(chunks)) as executor:
        futures = {executor.submit(_process_chunk_fp, chunk): i for i, chunk in enumerate(chunks)}
        for future in as_completed(futures):
            results[futures[future]] = future.result()
            
    fps_morgan, fps_maccs, fps_atompair, fps_topo, fps_avalon = [], [], [], [], []
    for r in results:
        fps_morgan.extend(r[0])
        fps_maccs.extend(r[1])
        fps_atompair.extend(r[2])
        fps_topo.extend(r[3])
        fps_avalon.extend(r[4])
        
    fps_morgan = np.array(fps_morgan)
    fps_maccs = np.array(fps_maccs)
    fps_atompair = np.array(fps_atompair)
    fps_topo = np.array(fps_topo)
    fps_avalon = np.array(fps_avalon)
    
    adme_features = np.column_stack((mw, logp, tpsa, min_es, max_es, asph, pmi1, pmi2, pmi3, rad_gyr))
    stds = adme_features.std(axis=0)
    stds[stds == 0] = 1e-8
    adme_norm = (adme_features - adme_features.mean(axis=0)) / stds
    
    adme_diff = adme_norm[:, np.newaxis, :] - adme_norm[np.newaxis, :, :]
    adme_dist = np.sqrt(np.sum(adme_diff**2, axis=-1))
    max_adme = np.max(adme_dist)
    if max_adme > 0: adme_dist = adme_dist / max_adme
    
    alerts_obj = np.array(alerts.tolist(), dtype=object)
    tox_dist = (alerts_obj[:, np.newaxis] != alerts_obj[np.newaxis, :]).astype(float)
    
    scaffolds_obj = np.array(scaffolds.tolist(), dtype=object)
    scaffold_dist = (scaffolds_obj[:, np.newaxis] != scaffolds_obj[np.newaxis, :]).astype(float)
    
    def tanimoto_dist(fp_matrix):
        dot_prod = np.dot(fp_matrix, fp_matrix.T)
        norm_sq = np.sum(fp_matrix, axis=1)
        denom = norm_sq[:, np.newaxis] + norm_sq[np.newaxis, :] - dot_prod
        denom[denom == 0] = 1
        return 1.0 - (dot_prod / denom)

    morgan_dist = tanimoto_dist(fps_morgan)
    maccs_dist = tanimoto_dist(fps_maccs)
    ap_dist = tanimoto_dist(fps_atompair)
    tt_dist = tanimoto_dist(fps_topo)
    av_dist = tanimoto_dist(fps_avalon)
    
    combined_fp_dist = (morgan_dist + maccs_dist + ap_dist + tt_dist + av_dist) / 5.0

    class_dist = np.zeros((n, n), dtype=float)
    class_weight = 0.0
    if class_layer_rules is not None:
        try:
            class_dist, _ = build_class_distance_matrix(df, class_layer_rules)
            class_weight = float(class_layer_rules.get("soft_layer_weight", 0.0))
        except Exception as e:
            print("Chemical class distance layer unavailable:", e)
            class_dist = np.zeros((n, n), dtype=float)
            class_weight = 0.0

    descriptor_dist = np.zeros((n, n), dtype=float)
    descriptor_weight = 0.0
    descriptor_metadata = {"status": "Not assessable", "available_columns": [], "coverage_ratio": 0.0}
    if descriptor_layer_rules is not None:
        try:
            descriptor_dist, descriptor_metadata = build_descriptor_distance_matrix(df, descriptor_layer_rules)
            descriptor_weight = float(descriptor_layer_rules.get("soft_layer_weight", 0.0))
        except Exception as e:
            print("Descriptor distance layer unavailable:", e)
            descriptor_dist = np.zeros((n, n), dtype=float)
            descriptor_weight = 0.0

    topology_dist = np.zeros((n, n), dtype=float)
    topology_weight = 0.0
    topology_metadata = {"status": "Not assessable", "available_columns": [], "coverage_ratio": 0.0}
    if topology_layer_rules is not None:
        try:
            topology_dist, topology_signatures = build_topology_distance_matrix(df, topology_layer_rules)
            topology_weight = float(topology_layer_rules.get("soft_layer_weight", 0.0))
            topology_metadata = {
                "status": "Calculated",
                "available_columns": topology_layer_rules.get("topology_columns", []),
                "coverage_ratio": 1.0 if topology_signatures else 0.0,
            }
        except Exception as e:
            print("Topology distance layer unavailable:", e)
            topology_dist = np.zeros((n, n), dtype=float)
            topology_weight = 0.0

    ionisation_dist = np.zeros((n, n), dtype=float)
    ionisation_weight = 0.0
    ionisation_metadata = {"status": "Not assessable", "available_columns": [], "coverage_ratio": 0.0}
    if ionisation_layer_rules is not None:
        try:
            ionisation_dist, ionisation_signatures = build_ionisation_distance_matrix(df, ionisation_layer_rules)
            ionisation_weight = float(ionisation_layer_rules.get("soft_layer_weight", 0.0))
            ionisation_metadata = {
                "status": "Calculated",
                "available_columns": ionisation_layer_rules.get("ionisation_columns", []),
                "coverage_ratio": 1.0 if ionisation_signatures else 0.0,
            }
        except Exception as e:
            print("Ionisation distance layer unavailable:", e)
            ionisation_dist = np.zeros((n, n), dtype=float)
            ionisation_weight = 0.0

    toxicophore_dist = np.zeros((n, n), dtype=float)
    toxicophore_weight = 0.0
    toxicophore_metadata = {"status": "Not assessable", "available_columns": [], "coverage_ratio": 0.0}
    if toxicophore_layer_rules is not None:
        try:
            toxicophore_dist, toxicophore_signatures = build_toxicophore_distance_matrix(df, toxicophore_layer_rules)
            toxicophore_weight = float(toxicophore_layer_rules.get("soft_layer_weight", 0.0))
            toxicophore_metadata = {
                "status": "Calculated",
                "available_columns": toxicophore_layer_rules.get("toxicophore_columns", []),
                "coverage_ratio": 1.0 if toxicophore_signatures else 0.0,
            }
        except Exception as e:
            print("Toxicophore distance layer unavailable:", e)
            toxicophore_dist = np.zeros((n, n), dtype=float)
            toxicophore_weight = 0.0

    chemont_dist = np.zeros((n, n), dtype=float)
    chemont_weight = 0.0
    chemont_metadata = {"status": "Not assessable", "available_columns": [], "coverage_ratio": 0.0}
    if chemont_layer_rules is not None:
        try:
            chemont_dist, chemont_signatures = build_chemont_distance_matrix(df, chemont_layer_rules)
            chemont_weight = float(chemont_layer_rules.get("soft_layer_weight", 0.0))
            chemont_metadata = {
                "status": "Calculated",
                "available_columns": chemont_layer_rules.get("hierarchy_columns", []),
                "coverage_ratio": 1.0 if chemont_signatures else 0.0,
            }
        except Exception as e:
            print("ChemOnt distance layer unavailable:", e)
            chemont_dist = np.zeros((n, n), dtype=float)
            chemont_weight = 0.0

    base_dist = (0.30 * scaffold_dist + 
                 0.25 * combined_fp_dist + 
                 0.25 * adme_dist + 
                 0.20 * tox_dist)

    soft_total = class_weight + descriptor_weight + topology_weight + ionisation_weight + toxicophore_weight + chemont_weight
    base_weight = max(0.0, 1.0 - soft_total)
    if soft_total > 0:
        dist_matrix = np.clip(
            base_weight * base_dist
            + class_weight * class_dist
            + descriptor_weight * descriptor_dist
            + topology_weight * topology_dist,
            0.0,
            1.0,
        )
        if ionisation_weight > 0:
            dist_matrix = np.clip(dist_matrix + ionisation_weight * ionisation_dist, 0.0, 1.0)
        if toxicophore_weight > 0:
            dist_matrix = np.clip(dist_matrix + toxicophore_weight * toxicophore_dist, 0.0, 1.0)
        if chemont_weight > 0:
            dist_matrix = np.clip(dist_matrix + chemont_weight * chemont_dist, 0.0, 1.0)
    else:
        dist_matrix = base_dist

    matrix_metadata = {
        "base_weight": round(base_weight, 3),
        "class_layer_weight": round(class_weight, 3),
        "descriptor_layer_weight": round(descriptor_weight, 3),
        "topology_layer_weight": round(topology_weight, 3),
        "ionisation_layer_weight": round(ionisation_weight, 3),
        "toxicophore_layer_weight": round(toxicophore_weight, 3),
        "chemont_layer_weight": round(chemont_weight, 3),
        "class_layer_status": "Calculated" if class_layer_rules is not None else "Not assessable",
        "descriptor_layer_status": descriptor_metadata.get("status", "Not assessable"),
        "descriptor_layer_available_columns": descriptor_metadata.get("available_columns", []),
        "descriptor_layer_coverage_ratio": descriptor_metadata.get("coverage_ratio", 0.0),
        "topology_layer_status": topology_metadata.get("status", "Not assessable"),
        "topology_layer_available_columns": topology_metadata.get("available_columns", []),
        "topology_layer_coverage_ratio": topology_metadata.get("coverage_ratio", 0.0),
        "ionisation_layer_status": ionisation_metadata.get("status", "Not assessable"),
        "ionisation_layer_available_columns": ionisation_metadata.get("available_columns", []),
        "ionisation_layer_coverage_ratio": ionisation_metadata.get("coverage_ratio", 0.0),
        "toxicophore_layer_status": toxicophore_metadata.get("status", "Not assessable"),
        "toxicophore_layer_available_columns": toxicophore_metadata.get("available_columns", []),
        "toxicophore_layer_coverage_ratio": toxicophore_metadata.get("coverage_ratio", 0.0),
        "chemont_layer_status": chemont_metadata.get("status", "Not assessable"),
        "chemont_layer_available_columns": chemont_metadata.get("available_columns", []),
        "chemont_layer_coverage_ratio": chemont_metadata.get("coverage_ratio", 0.0),
    }
                   
    return dist_matrix, adme_norm, alerts_obj, matrix_metadata

def run_consensus_clustering(
    df: pd.DataFrame,
    smiles_col: str = 'Standardized SMILES',
    cutoff: float = 0.4,
    class_layer_rules: dict | None = None,
    descriptor_layer_rules: dict | None = None,
    topology_layer_rules: dict | None = None,
    ionisation_layer_rules: dict | None = None,
    toxicophore_layer_rules: dict | None = None,
    chemont_layer_rules: dict | None = None,
):
    n = len(df)

    def _has_valid_structure(value: object) -> bool:
        if value is None or (isinstance(value, float) and np.isnan(value)):
            return False
        text = str(value).strip()
        if not text:
            return False
        try:
            return Chem.MolFromSmiles(text) is not None
        except Exception:
            return False

    # Rows without a parseable standardized structure are not eligible for
    # similarity clustering. Process the eligible population independently so
    # invalid rows cannot influence cluster membership, centroids, or domain
    # decisions, then return explicit review-only statuses for those rows.
    valid_mask = df[smiles_col].map(_has_valid_structure) if smiles_col in df.columns else pd.Series(False, index=df.index)
    if n and not bool(valid_mask.all()):
        valid_indices = list(df.index[valid_mask])
        invalid_indices = list(df.index[~valid_mask])
        if valid_indices:
            eligible = df.loc[valid_indices].copy()
            processed = run_consensus_clustering(
                eligible,
                smiles_col=smiles_col,
                cutoff=cutoff,
                class_layer_rules=class_layer_rules,
                descriptor_layer_rules=descriptor_layer_rules,
                topology_layer_rules=topology_layer_rules,
                ionisation_layer_rules=ionisation_layer_rules,
                toxicophore_layer_rules=toxicophore_layer_rules,
                chemont_layer_rules=chemont_layer_rules,
            )
            for column in processed.columns:
                if column not in df.columns:
                    df[column] = None
                # Excel import can preserve an existing enrichment column as
                # pandas StringDtype.  B1a then legitimately produces numeric
                # values (for example Cluster Size) or structured evidence for
                # the eligible rows.  Widen only that destination column
                # before copy-back instead of coercing the calculated value to
                # text or failing the whole Phase B run.
                destination = df.loc[:, column]
                if isinstance(destination, pd.DataFrame):
                    raise ValueError(f"Duplicate enrichment column encountered: {column}")
                if isinstance(destination.dtype, pd.StringDtype):
                    df[column] = df[column].astype(object)
                df.loc[valid_indices, column] = processed[column].values
            df.attrs.update(getattr(processed, "attrs", {}))
        else:
            # Keep the expected cluster/domain columns present even when the
            # entire input population lacks a valid structure.
            for column in ("Cluster ID", "Cluster Size", "Cluster_Status", "Domain", "Decision"):
                df[column] = None
        df.loc[invalid_indices, "Cluster ID"] = "Not assessable"
        df.loc[invalid_indices, "Cluster Size"] = "Not assessable"
        df.loc[invalid_indices, "Cluster_Status"] = "Review Required"
        df.loc[invalid_indices, "Domain"] = "Not assessable"
        df.loc[invalid_indices, "Decision"] = "REVIEW"
        return df

    
    df['Cluster ID'] = [None] * n
    df['Cluster Size'] = [None] * n
    df['Cluster_Status'] = [""] * n
    df['Domain'] = [""] * n
    df['Decision'] = [""] * n
    
    if n == 0:
        return df
    if n == 1:
        df.at[df.index[0], 'Cluster ID'] = "Cluster_001"
        df.at[df.index[0], 'Cluster Size'] = 1
        df.at[df.index[0], 'Cluster_Status'] = "Valid Group (Size 1)"
        df.at[df.index[0], 'Domain'] = "INSIDE DOMAIN"
        df.at[df.index[0], 'Decision'] = "PASS"
        return df
        
    distance_result = compute_distance_matrix(
        df,
        smiles_col,
        class_layer_rules=class_layer_rules,
        descriptor_layer_rules=descriptor_layer_rules,
        topology_layer_rules=topology_layer_rules,
        ionisation_layer_rules=ionisation_layer_rules,
        toxicophore_layer_rules=toxicophore_layer_rules,
        chemont_layer_rules=chemont_layer_rules,
    )
    # Backward-compatible unpacking for injected/legacy distance providers
    # that return the original three-value contract.
    if len(distance_result) == 3:
        dist_matrix, features_norm, alerts_array = distance_result
        matrix_metadata = {}
    else:
        dist_matrix, features_norm, alerts_array, matrix_metadata = distance_result
    
    dists = []
    for i in range(1, n):
        for j in range(i):
            dists.append(dist_matrix[i, j])
    
    butina_clusters = Butina.ClusterData(dists, n, cutoff, isDistData=True)
    b_labels = np.zeros(n, dtype=int)
    for c_id, c_indices in enumerate(butina_clusters):
        for idx in c_indices: b_labels[idx] = c_id
        
    db = DBSCAN(eps=cutoff, min_samples=2, metric='precomputed')
    db_labels = db.fit_predict(dist_matrix)
    
    n_clus = max(2, n // 5) if n > 5 else 2
    if n > 2:
        agg = AgglomerativeClustering(n_clusters=n_clus, metric='precomputed', linkage='average')
        agg_labels = agg.fit_predict(dist_matrix)
    else:
        agg_labels = np.zeros(n, dtype=int)
        
    if n > 2:
        aff_matrix = np.exp(-dist_matrix ** 2 / (2.0 * (cutoff/2) ** 2))
        spec = SpectralClustering(n_clusters=n_clus, affinity='precomputed', random_state=42)
        spec_labels = spec.fit_predict(aff_matrix)
    else:
        spec_labels = np.zeros(n, dtype=int)
        
    co_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            match_count = 0
            if b_labels[i] == b_labels[j]: match_count += 1
            if db_labels[i] == db_labels[j]: match_count += 1
            if agg_labels[i] == agg_labels[j]: match_count += 1
            if spec_labels[i] == spec_labels[j]: match_count += 1
            co_matrix[i, j] = match_count / 4.0
            
    consensus_dist = 1.0 - co_matrix
    
    final_dists = []
    for i in range(1, n):
        for j in range(i):
            final_dists.append(consensus_dist[i, j])
            
    final_clusters = Butina.ClusterData(final_dists, n, 0.5, isDistData=True)
    
    final_labels = np.zeros(n, dtype=int)
    for c_id, c_indices in enumerate(final_clusters):
        for idx in c_indices: final_labels[idx] = c_id
        
    global_silhouette = -1.0
    if len(set(final_labels)) > 1:
        try:
            global_silhouette = silhouette_score(dist_matrix, final_labels, metric='precomputed')
        except:
            pass
            
    cluster_counter = 1
    for c_id, c_indices in enumerate(final_clusters):
        real_indices = [df.index[i] for i in c_indices]
        c_size = len(real_indices)
        
        domain_flags = {}
        cluster_alerts = set()
        
        if c_size > 1:
            c_features = features_norm[list(c_indices)]
            centroid = c_features.mean(axis=0)
            threshold = 2.0 * np.sqrt(features_norm.shape[1])
            
            for i, idx in zip(c_indices, real_indices):
                dist_to_centroid = np.linalg.norm(features_norm[i] - centroid)
                if dist_to_centroid > threshold:
                    domain_flags[idx] = "OUTSIDE DOMAIN"
                else:
                    domain_flags[idx] = "INSIDE DOMAIN"
                cluster_alerts.add(alerts_array[i])
        else:
            for i, idx in zip(c_indices, real_indices):
                domain_flags[idx] = "INSIDE DOMAIN (Size 1)"
                cluster_alerts.add(alerts_array[i])
                
        c_stat = "Valid Group"
        if c_size > 2 and global_silhouette != -1.0:
            c_stat = f"Valid Group (Silhouette: {global_silhouette:.2f})"
            if global_silhouette < 0.2:
                c_stat = f"Unstable Cluster (Silhouette: {global_silhouette:.2f})"
                
        for i, idx in zip(c_indices, real_indices):
            df.at[idx, 'Cluster ID'] = f"Cluster_{cluster_counter:03d}"
            df.at[idx, 'Cluster Size'] = c_size
            df.at[idx, 'Cluster_Status'] = c_stat
            df.at[idx, 'Domain'] = domain_flags[idx]
            
            decision = "PASS"
            if "Unstable" in c_stat or domain_flags[idx] == "OUTSIDE DOMAIN":
                decision = "FAIL"
            elif len(cluster_alerts) > 1:
                decision = "WARNING"
                
            df.at[idx, 'Decision'] = decision
            
        cluster_counter += 1

    if class_layer_rules is not None:
        try:
            df = annotate_cluster_class_evidence(df, class_layer_rules)
        except Exception as e:
            print("Chemical class cluster annotation unavailable:", e)

    if descriptor_layer_rules is not None:
        try:
            df = annotate_descriptor_layer_evidence(df, descriptor_layer_rules)
        except Exception as e:
            print("Descriptor cluster annotation unavailable:", e)

    if topology_layer_rules is not None:
        try:
            df = annotate_topology_layer_evidence(df, topology_layer_rules)
        except Exception as e:
            print("Topology cluster annotation unavailable:", e)

    if ionisation_layer_rules is not None:
        try:
            df = annotate_ionisation_layer_evidence(df, ionisation_layer_rules)
        except Exception as e:
            print("Ionisation cluster annotation unavailable:", e)

    if toxicophore_layer_rules is not None:
        try:
            df = annotate_toxicophore_layer_evidence(df, toxicophore_layer_rules)
        except Exception as e:
            print("Toxicophore cluster annotation unavailable:", e)

    if chemont_layer_rules is not None:
        try:
            df = annotate_chemont_layer_evidence(df, chemont_layer_rules)
        except Exception as e:
            print("ChemOnt cluster annotation unavailable:", e)

    df.attrs["distance_layer_metadata"] = matrix_metadata

    return df


def _safe_text(value: object) -> str:
    """Return a stable, non-placeholder text value for audit fields."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    text = str(value).strip()
    return "" if text.casefold() in {"", "nan", "none", "not available", "not assessable"} else text


def _cluster_compatibility_family(row: pd.Series, mol) -> tuple[str, str, str]:
    """Assign a deterministic pre-clustering compatibility domain.

    This is intentionally a chemical compatibility gate, not a toxicological
    conclusion.  It prevents structurally incompatible domains from reaching
    the same molecular-similarity partition.
    """
    smiles = _safe_text(row.get("Standardized SMILES"))
    alerts = _safe_text(row.get("All_Structural_Alerts")) or _safe_text(row.get("Alerts"))
    primary_group = _safe_text(row.get("Primary Functional Group"))
    corrected_class = _safe_text(row.get("Corrected Chemical Class"))
    has_silicon = "Si" in smiles
    has_peroxide = "peroxide" in (alerts + " " + primary_group + " " + corrected_class).casefold()
    formal_charge = pd.to_numeric(pd.Series([row.get("Formal_Charge")]), errors="coerce").fillna(0).iloc[0]

    if has_peroxide:
        return (
            "Reactive peroxide",
            "Reactive peroxide / separate review domain",
            "Peroxide functionality requires separate reactivity review and is not merged with non-peroxide analogues.",
        )
    if formal_charge != 0:
        return (
            "Charged compound",
            "Persistent charge / separate review domain",
            "Persistent formal charge is a compatibility boundary for the primary neutral-structure clustering space.",
        )
    if has_silicon:
        ring_count = mol.GetRingInfo().NumRings() if mol is not None else 0
        silicon_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 14) if mol is not None else 0
        if ring_count and silicon_atoms >= 4:
            ring_units = _silicon_ring_units(mol)
            return (
                f"Macrocyclic siloxane n={ring_units}",
                f"Siloxane / macrocyclic n={ring_units}",
                f"Cyclic Si-O backbone with {ring_units} silicon ring unit(s); only like-sized macrocyclic siloxanes enter the same primary clustering space.",
            )
        if ring_count:
            return (
                "Cyclic organosilicon",
                "Organosilicon / cyclic",
                "Cyclic organosilicon topology is evaluated separately from linear and macrocyclic siloxanes.",
            )
        if "O[Si]" in smiles or "[Si]O" in smiles:
            return (
                "Linear or oligomeric siloxane",
                "Siloxane / linear or oligomeric",
                "Acyclic Si-O backbone; evaluated separately from cyclic siloxanes.",
            )
        return (
            "Other organosilicon compound",
            "Organosilicon / other",
            "Organosilicon structure without a compatible siloxane backbone; retained in a separate review domain.",
        )

    if primary_group:
        return (
            f"Organic {primary_group}",
            f"Organic / {primary_group}",
            f"Primary functional-group compatibility gate: {primary_group}.",
        )
    if corrected_class:
        return (
            f"Organic {corrected_class}",
            f"Organic / {corrected_class}",
            f"Corrected chemical-class compatibility gate: {corrected_class}.",
        )
    return (
        "Organic structure pending functional-group review",
        "Organic / broad provisional",
        "No validated functional-group family was available; membership requires SME review.",
    )


def _silicon_ring_units(mol) -> int:
    """Return the largest count of silicon atoms in one perceived ring."""
    if mol is None:
        return 0
    counts = [
        sum(1 for atom_index in ring if mol.GetAtomWithIdx(atom_index).GetAtomicNum() == 14)
        for ring in mol.GetRingInfo().AtomRings()
    ]
    return max(counts, default=0)


def _cluster_letter(position: int) -> str:
    """Excel-style deterministic suffix: 0 -> A, 25 -> Z, 26 -> AA."""
    value = position + 1
    letters = ""
    while value:
        value, remainder = divmod(value - 1, 26)
        letters = chr(65 + remainder) + letters
    return letters


def _morgan_fingerprint(mol):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    return generator.GetFingerprint(mol)


PROPERTY_COMPATIBILITY_POLICY_VERSION = "1.0.0"
PROPERTY_COMPATIBILITY_MAX_MW_RATIO = 1.50
PROPERTY_COMPATIBILITY_MAX_LOGP_DELTA = 1.75
PROPERTY_COMPATIBILITY_MAX_TPSA_RATIO = 1.75
SUBSTANTIVE_SECONDARY_FEATURES = {
    "alcohol", "phenol", "amine", "carboxylic acid", "sulfonic acid",
    "peroxide", "epoxide", "phosphate ester", "phosphine oxide",
}


def _numeric_value(row: pd.Series, column: str) -> float | None:
    value = pd.to_numeric(pd.Series([row.get(column)]), errors="coerce").iloc[0]
    return None if pd.isna(value) else float(value)


def _substantive_secondary_features(row: pd.Series) -> set[str]:
    """Return material secondary functional groups relevant to analogue use."""
    values = (_safe_text(row.get("Secondary Functional Groups")) + ";" + _safe_text(row.get("Detected Feature Profile"))).split(";")
    return {
        value.strip().casefold()
        for value in values
        if value.strip().casefold() in SUBSTANTIVE_SECONDARY_FEATURES
    }


def _amide_unit_count(mol) -> int:
    """Count amide/lactam carbonyl units for conservative oligomer separation."""
    if mol is None:
        return 0
    pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    return len(mol.GetSubstructMatches(pattern)) if pattern is not None else 0


def _pair_cluster_compatible(left: pd.Series, right: pd.Series) -> tuple[bool, str]:
    """Require chemistry and property comparability in addition to similarity.

    Morgan fingerprints can saturate for repeated units.  This gate prevents a
    short molecule and a much larger oligomer, or an ester and hydroxy ester,
    from becoming a read-across candidate merely because their local fragments
    are similar.
    """
    left_mol, right_mol = left.get("_mol"), right.get("_mol")
    left_amides, right_amides = _amide_unit_count(left_mol), _amide_unit_count(right_mol)
    if left_amides and right_amides and left_amides != right_amides:
        return False, f"Amide/lactam repeat-unit count differs ({left_amides} vs {right_amides})."

    left_features, right_features = _substantive_secondary_features(left), _substantive_secondary_features(right)
    if left_features != right_features:
        return False, "Material secondary functional-group profile differs (" + ", ".join(sorted(left_features or {"none"})) + " vs " + ", ".join(sorted(right_features or {"none"})) + ")."

    left_mw, right_mw = _numeric_value(left, "MW"), _numeric_value(right, "MW")
    if left_mw and right_mw:
        ratio = max(left_mw, right_mw) / min(left_mw, right_mw)
        if ratio > PROPERTY_COMPATIBILITY_MAX_MW_RATIO:
            return False, f"Molecular-weight ratio {ratio:.2f} exceeds {PROPERTY_COMPATIBILITY_MAX_MW_RATIO:.2f}."

    left_logp, right_logp = _numeric_value(left, "LogP"), _numeric_value(right, "LogP")
    if left_logp is not None and right_logp is not None:
        delta = abs(left_logp - right_logp)
        if delta > PROPERTY_COMPATIBILITY_MAX_LOGP_DELTA:
            return False, f"LogP difference {delta:.2f} exceeds {PROPERTY_COMPATIBILITY_MAX_LOGP_DELTA:.2f}."

    left_tpsa, right_tpsa = _numeric_value(left, "TPSA"), _numeric_value(right, "TPSA")
    if left_tpsa and right_tpsa:
        ratio = max(left_tpsa, right_tpsa) / min(left_tpsa, right_tpsa)
        if ratio > PROPERTY_COMPATIBILITY_MAX_TPSA_RATIO:
            return False, f"TPSA ratio {ratio:.2f} exceeds {PROPERTY_COMPATIBILITY_MAX_TPSA_RATIO:.2f}."
    return True, "Passed functional-group, molecular-weight, LogP, and TPSA compatibility checks."


def _partition_property_compatible(members: tuple, member_rows: pd.DataFrame) -> tuple[bool, str]:
    rows = member_rows.iloc[list(members)]
    reasons = []
    for left in range(1, len(rows)):
        for right in range(left):
            compatible, reason = _pair_cluster_compatible(rows.iloc[left], rows.iloc[right])
            if not compatible:
                return False, reason
            reasons.append(reason)
    return True, reasons[0] if reasons else "Not assessable: one unique structure."


def _siloxane_identity_consistency(row: pd.Series, mol) -> tuple[str, str]:
    """Surface unresolved source-name/structure chain-length conflicts."""
    name = _safe_text(row.get("Compound Name"))
    if "siloxane" not in name.casefold():
        return "Not applicable", ""
    match = re.search(r"\(\s*n\s*=\s*0*(\d+)\s*\)", name, flags=re.IGNORECASE)
    if match is None or mol is None:
        return "Not assessable", "No explicit siloxane chain-length declaration was available for comparison."
    declared = int(match.group(1))
    resolved = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 14)
    if declared != resolved:
        return "Review required", f"Source name declares siloxane n={declared}, while the resolved standardized structure contains {resolved} silicon atom(s)."
    return "Consistent", f"Source name declares siloxane n={declared}, matching {resolved} silicon atom(s) in the resolved standardized structure."


def _pairwise_similarity_summary(fingerprints: list) -> tuple[float | None, float | None]:
    if len(fingerprints) < 2:
        return None, None
    values = [
        DataStructs.TanimotoSimilarity(fingerprints[left], fingerprints[right])
        for left in range(1, len(fingerprints))
        for right in range(left)
    ]
    return round(float(min(values)), 3), round(float(np.median(values)), 3)


def _complete_linkage_refine(partitions: tuple, fingerprints: list, member_rows: pd.DataFrame, similarity_cutoff: float) -> list[tuple]:
    """Split groups unless every pair passes similarity and chemistry gates."""
    refined: list[tuple] = []
    for partition in partitions:
        members = tuple(sorted(partition))
        member_fps = [fingerprints[index] for index in members]
        minimum, _ = _pairwise_similarity_summary(member_fps)
        compatible, _ = _partition_property_compatible(members, member_rows)
        if len(members) < 2 or (minimum is not None and minimum >= similarity_cutoff and compatible):
            refined.append(members)
            continue
        local_distances = np.zeros((len(member_fps), len(member_fps)))
        for left in range(1, len(member_fps)):
            for right in range(left):
                pair_compatible, _ = _pair_cluster_compatible(member_rows.iloc[members[left]], member_rows.iloc[members[right]])
                distance = 1.0 - DataStructs.TanimotoSimilarity(member_fps[left], member_fps[right])
                # An incompatible pair must never be joined by complete linkage.
                local_distances[left, right] = local_distances[right, left] = distance if pair_compatible else 1.01
        labels = AgglomerativeClustering(
            n_clusters=None,
            distance_threshold=1.0 - similarity_cutoff,
            metric="precomputed",
            linkage="complete",
        ).fit_predict(local_distances)
        for label in sorted(set(labels)):
            refined.append(tuple(members[position] for position, value in enumerate(labels) if value == label))
    return refined


def run_revalidated_clustering(
    df: pd.DataFrame,
    smiles_col: str = "Standardized SMILES",
    similarity_cutoff: float = 0.60,
) -> pd.DataFrame:
    """Run the transparent, deduplicated Phase B clustering model.

    Each unique standardized structure is clustered once with Morgan/ECFP4
    Tanimoto similarity.  Original rows are then mapped back to the resulting
    subcluster so aliases, homologues, and duplicate source records remain
    fully visible without repeatedly influencing membership.
    """
    result = df.copy()
    # PDF/table consolidation and recovered Excel workbooks can retain
    # duplicate source index labels.  They are not business identifiers and
    # cannot safely participate in pandas alignment during clustering.  Keep
    # every row while normalising the internal working index.
    if not result.index.is_unique:
        result = result.reset_index(drop=True)
    # Reader-facing terminology: this field identifies the broad structural
    # compatibility space, not a taxonomic chemical-family classification.
    # Accept older saved workbooks while ensuring new exports use one name.
    if "Cluster Family" in result.columns:
        if "Compatibility Group" not in result.columns:
            result = result.rename(columns={"Cluster Family": "Compatibility Group"})
        else:
            result = result.drop(columns=["Cluster Family"])
    required_columns = [
        "Cluster ID", "Compatibility Group", "Subcluster ID", "Cluster Size",
        "Cluster Unique Structure Count", "Cluster_Status", "Domain", "Decision",
        "Compatibility Gate", "Compatibility Gate Reason", "Duplicate Structure Count",
        "Nearest_Neighbour_Tanimoto", "Cluster_Min_Pairwise_Tanimoto",
        "Cluster_Median_Pairwise_Tanimoto", "Cluster_Consensus_Method",
        "Cluster_Membership_Parameters", "Cluster_Membership_Rationale",
        "Cluster Property Compatibility", "Cluster Property Compatibility Reason",
        "Identity Structure Consistency", "Identity Structure Consistency Reason",
    ]
    for column in required_columns:
        # Excel imports can retain prior enrichment columns as strict pandas
        # StringDtype.  The revalidated model writes numeric counts and
        # similarities to some of those columns, so widen only these
        # calculated destinations before assignment.
        if column in result.columns and isinstance(result[column].dtype, pd.StringDtype):
            result[column] = result[column].astype(object)
        if column not in result.columns:
            result[column] = ""
        # New blank columns can also be inferred as the strict string dtype in
        # recent pandas versions.  Every field in this calculated output set
        # may receive a numeric value, text, or an explicit not-assessable
        # marker, so object is the correct transport dtype until export.
        result[column] = result[column].astype(object)

    if result.empty or smiles_col not in result.columns:
        return result

    def valid_molecule(value: object):
        text = _safe_text(value)
        return Chem.MolFromSmiles(text) if text else None

    mols = result[smiles_col].map(valid_molecule)
    valid_mask = mols.map(lambda mol: mol is not None)
    invalid_indices = list(result.index[~valid_mask])
    if invalid_indices:
        result.loc[invalid_indices, "Cluster ID"] = "Not assessable"
        result.loc[invalid_indices, "Subcluster ID"] = "Not assessable"
        result.loc[invalid_indices, "Cluster_Status"] = "Review Required"
        result.loc[invalid_indices, "Domain"] = "Not assessable"
        result.loc[invalid_indices, "Decision"] = "REVIEW"
        result.loc[invalid_indices, "Compatibility Gate"] = "Identity unresolved"
        result.loc[invalid_indices, "Compatibility Gate Reason"] = "No valid standardized structure is available for structural clustering."

    valid = result.loc[valid_mask].copy()
    if valid.empty:
        return result
    valid["_mol"] = mols.loc[valid.index]
    valid["_structure_key"] = valid[smiles_col].map(_safe_text)
    identity_checks = valid.apply(lambda row: _siloxane_identity_consistency(row, row["_mol"]), axis=1, result_type="expand")
    identity_checks.columns = ["Identity Structure Consistency", "Identity Structure Consistency Reason"]
    for column in identity_checks.columns:
        valid[column] = identity_checks[column]
    gates = valid.apply(lambda row: _cluster_compatibility_family(row, row["_mol"]), axis=1, result_type="expand")
    gates.columns = ["_family", "_gate", "_gate_reason"]
    # ``join`` routes through pandas' merge machinery and can fail when a
    # recovered Excel payload has an unusual index/metadata shape.  Gates are
    # calculated row-for-row from ``valid``, so positional assignment is the
    # correct operation and avoids any merge-driven row duplication.
    if len(gates) != len(valid):
        raise ValueError("Compatibility-gate calculation did not return one result per valid structure.")
    for column in gates.columns:
        valid[column] = gates[column].to_numpy()

    # Structure-level representation: duplicate source rows become aliases,
    # never independent neighbours or votes in the clustering algorithm.
    unique = valid.sort_index(kind="stable").drop_duplicates("_structure_key", keep="first").copy()
    unique["_fp"] = unique["_mol"].map(_morgan_fingerprint)
    family_order = sorted(unique["_family"].unique(), key=lambda value: str(value).casefold())
    assignments: dict[str, dict] = {}
    # Cluster ID is the global, reviewer-facing sequence.  A compatibility
    # family may contain multiple independently reviewable clusters, so the
    # family remains separate audit context rather than leaking letter-suffixed
    # subcluster IDs into the primary workbook view.
    cluster_number = 0

    for family_number, family in enumerate(family_order, start=1):
        family_id = f"Family_{family_number:03d}"
        family_unique = unique[unique["_family"] == family].copy()
        fingerprints = family_unique["_fp"].tolist()
        count = len(family_unique)
        if count == 1:
            partitions = [(0,)]
        else:
            distances = [
                1.0 - DataStructs.TanimotoSimilarity(fingerprints[left], fingerprints[right])
                for left in range(1, count)
                for right in range(left)
            ]
            # Butina is the transparent primary assignment method.  The
            # threshold is explicit: distance <= 0.40 (similarity >= 0.60).
            partitions = Butina.ClusterData(distances, count, 1.0 - similarity_cutoff, isDistData=True, reordering=True)
            partitions = _complete_linkage_refine(partitions, fingerprints, family_unique, similarity_cutoff)

        ordered_partitions = sorted(
            [tuple(sorted(group)) for group in partitions],
            key=lambda group: min(group),
        )
        for subgroup_position, members in enumerate(ordered_partitions):
            cluster_number += 1
            cluster_id = f"Cluster_{cluster_number:03d}"
            member_rows = family_unique.iloc[list(members)]
            member_fps = member_rows["_fp"].tolist()
            minimum, median = _pairwise_similarity_summary(member_fps)
            structure_keys = member_rows["_structure_key"].tolist()
            source_members = valid[valid["_structure_key"].isin(structure_keys)]
            source_count = len(source_members)
            unique_count = len(structure_keys)
            duplicate_count = source_count - unique_count
            row_gate = member_rows["_gate"].iloc[0]
            gate_reason = member_rows["_gate_reason"].iloc[0]
            property_compatible, property_reason = _partition_property_compatible(members, family_unique)
            review_required = (
                unique_count < 2
                or "separate review domain" in row_gate.casefold()
                or median is None
                or median < similarity_cutoff
            )
            status = "SME Review Required" if review_required else "Structurally coherent pending SME review"
            rationale = (
                f"{row_gate}. {gate_reason} Primary assignment used deduplicated Morgan/ECFP4 "
                f"fingerprints with Tanimoto similarity >= {similarity_cutoff:.2f}. "
                f"{source_count} source row(s) represent {unique_count} unique structure(s); "
                f"{duplicate_count} duplicate/alias row(s) did not influence clustering repeatedly. "
                f"Property compatibility: {property_reason}"
            )
            for structure_key in structure_keys:
                assignments[structure_key] = {
                    "Cluster ID": cluster_id,
                    "Compatibility Group": family_id,
                    # Retained internally for compatibility with existing
                    # evidence functions; it is intentionally not exported.
                    "Subcluster ID": cluster_id,
                    "Cluster Size": source_count,
                    "Cluster Unique Structure Count": unique_count,
                    "Cluster_Status": status,
                    "Domain": "REVIEW REQUIRED" if review_required else "STRUCTURALLY COMPATIBLE",
                    "Decision": "REVIEW",
                    "Compatibility Gate": row_gate,
                    "Compatibility Gate Reason": gate_reason,
                    "Duplicate Structure Count": duplicate_count,
                    "Cluster_Min_Pairwise_Tanimoto": "Not assessable" if minimum is None else minimum,
                    "Cluster_Median_Pairwise_Tanimoto": "Not assessable" if median is None else median,
                    "Cluster_Consensus_Method": "Primary Morgan/ECFP4 Tanimoto Butina; complete-linkage diameter refinement; compatibility-gated; duplicate structures collapsed",
                    "Cluster_Membership_Parameters": f"Morgan radius=2; fpSize=1024; Tanimoto >= {similarity_cutoff:.2f}; Butina reordering=True; complete-linkage max distance <= {1.0 - similarity_cutoff:.2f}",
                    "Cluster_Membership_Rationale": rationale,
                    "Cluster Property Compatibility": "Passed" if property_compatible else "Failed",
                    "Cluster Property Compatibility Reason": property_reason,
                }

    for index, row in valid.iterrows():
        assignment = assignments[row["_structure_key"]]
        for column, value in assignment.items():
            result.at[index, column] = value
        result.at[index, "Identity Structure Consistency"] = row["Identity Structure Consistency"]
        result.at[index, "Identity Structure Consistency Reason"] = row["Identity Structure Consistency Reason"]

    # Nearest-neighbour evidence is calculated from unique structures in each
    # subcluster, then replicated to aliases.  It is not derived from a row's
    # own duplicate structure.
    for subcluster_id, group in result.loc[valid.index].groupby("Subcluster ID", sort=False):
        unique_group = group.drop_duplicates(smiles_col, keep="first")
        if len(unique_group) < 2:
            result.loc[group.index, "Nearest_Neighbour_Tanimoto"] = "Not assessable"
            continue
        fps = [_morgan_fingerprint(Chem.MolFromSmiles(_safe_text(smiles))) for smiles in unique_group[smiles_col]]
        nearest = []
        for current, fingerprint in enumerate(fps):
            nearest.append(max(DataStructs.TanimotoSimilarity(fingerprint, other) for other_index, other in enumerate(fps) if other_index != current))
        nearest_by_structure = dict(zip(unique_group[smiles_col].map(_safe_text), [round(value, 3) for value in nearest]))
        result.loc[group.index, "Nearest_Neighbour_Tanimoto"] = group[smiles_col].map(lambda value: nearest_by_structure.get(_safe_text(value), "Not assessable"))

    result.attrs["clustering_method"] = "revalidated_structure_deduplicated_morgan_butina_v1"
    result.attrs["clustering_similarity_cutoff"] = similarity_cutoff
    return result

def _evaluate_cramer_ttc(row: pd.Series) -> dict:
    """
    Transparent Cramer/TTC proxy used by the current workbook workflow.

    The exact formal Cramer variant still needs to be frozen against the
    reference Excel logic. Until that is confirmed, we keep the existing
    coarse decision structure but make the branching explicit and auditable.
    TTC is expressed in ug/kg bw/day.
    """
    mw = row.get('MW', 0) or 0
    alerts = str(row.get('Alerts', ''))
    scaffolds = str(row.get('Scaffold', ''))

    path = []
    evidence = []
    cramer = "Class III"
    ttc = "1.5 µg/kg bw/day"
    ttc_value = 1.5
    rule_version = "heuristic-v1"

    has_alerts = alerts not in ("None", "N/A", "", "nan", "NaN")
    has_scaffold = scaffolds not in ("No_Scaffold", "None", "N/A", "", "nan", "NaN")

    path.append(f"Structural alerts present: {'yes' if has_alerts else 'no'}")
    evidence.append(f"Alerts={alerts}")
    if has_alerts:
        cramer = "Class III"
        ttc = "1.5 µg/kg bw/day"
        ttc_value = 1.5
        path.append("Branch: alerts -> Class III")
        return {
            "Cramer_Class": cramer,
            "TTC_Limit": ttc,
            "TTC_Value": ttc_value,
            "TTC_Unit": "µg/kg bw/day",
            "Cramer_Rule_Version": rule_version,
            "Cramer_Path": " | ".join(path),
            "Cramer_Evidence": " | ".join(evidence),
        }

    path.append(f"Scaffold present: {'yes' if has_scaffold else 'no'}")
    evidence.append(f"Scaffold={scaffolds}")
    path.append(f"Molecular weight: {mw}")

    if mw < 200 and not has_scaffold:
        cramer = "Class I"
        ttc = "30 µg/kg bw/day"
        ttc_value = 30
        path.append("Branch: MW < 200 and no scaffold -> Class I")
    elif mw < 400 and has_scaffold:
        cramer = "Class II"
        ttc = "9 µg/kg bw/day"
        ttc_value = 9
        path.append("Branch: MW < 400 and scaffold present -> Class II")
    else:
        cramer = "Class III"
        ttc = "1.5 µg/kg bw/day"
        ttc_value = 1.5
        path.append("Branch: fallback -> Class III")

    return {
        "Cramer_Class": cramer,
        "TTC_Limit": ttc,
        "TTC_Value": ttc_value,
        "TTC_Unit": "µg/kg bw/day",
        "Cramer_Rule_Version": rule_version,
        "Cramer_Path": " | ".join(path),
        "Cramer_Evidence": " | ".join(evidence),
    }


def assign_cramer_and_ttc(df: pd.DataFrame):
    cramer_classes = []
    ttcs = []
    ttc_values = []
    ttc_units = []
    rule_versions = []
    paths = []
    evidence_list = []

    for _, row in df.iterrows():
        result = _evaluate_cramer_ttc(row)
        cramer_classes.append(result["Cramer_Class"])
        ttcs.append(result["TTC_Limit"])
        ttc_values.append(result["TTC_Value"])
        ttc_units.append(result["TTC_Unit"])
        rule_versions.append(result["Cramer_Rule_Version"])
        paths.append(result["Cramer_Path"])
        evidence_list.append(result["Cramer_Evidence"])

    df['Cramer_Class'] = cramer_classes
    df['TTC_Limit'] = ttcs
    df['TTC_Value'] = ttc_values
    df['TTC_Unit'] = ttc_units
    df['Cramer_Rule_Version'] = rule_versions
    df['Cramer_Path'] = paths
    df['Cramer_Evidence'] = evidence_list
    return df

def process_clustering_and_adme(df: pd.DataFrame, smiles_col: str = 'Standardized SMILES', debug_trace: list[dict] | None = None) -> pd.DataFrame:
    df = df.copy()
    if smiles_col not in df.columns or df.empty:
        return df
    df = process_clustering_and_adme_core(df, smiles_col=smiles_col, debug_trace=debug_trace)
    df = process_clustering_and_adme_optional(df, smiles_col=smiles_col, debug_trace=debug_trace)
    return df


def _load_phase_b_policies():
    try:
        clustering_policy = load_clustering_policy()
    except Exception as e:
        print("Clustering policy unavailable:", e)
        clustering_policy = None
    try:
        class_layer_rules = load_class_layer_rules()
    except Exception as e:
        print("Chemical class layer policy unavailable:", e)
        class_layer_rules = None
    try:
        descriptor_layer_rules = load_descriptor_layer_rules()
    except Exception as e:
        print("Descriptor layer policy unavailable:", e)
        descriptor_layer_rules = None
    try:
        topology_layer_rules = load_topology_layer_rules()
    except Exception as e:
        print("Topology layer policy unavailable:", e)
        topology_layer_rules = None
    try:
        ionisation_layer_rules = load_ionisation_layer_rules()
    except Exception as e:
        print("Ionisation layer policy unavailable:", e)
        ionisation_layer_rules = None
    try:
        toxicophore_layer_rules = load_toxicophore_layer_rules()
    except Exception as e:
        print("Toxicophore layer policy unavailable:", e)
        toxicophore_layer_rules = None
    try:
        chemont_layer_rules = load_chemont_layer_rules()
    except Exception as e:
        print("ChemOnt layer policy unavailable:", e)
        chemont_layer_rules = None
    # In deterministic-local mode, external taxonomy is optional evidence and
    # must not reserve any distance weight when the service is skipped.
    if os.environ.get("DETERMINISTIC_LOCAL_MODE", "false").strip().casefold() in {"1", "true", "yes", "on"}:
        if chemont_layer_rules is not None:
            chemont_layer_rules = dict(chemont_layer_rules)
            chemont_layer_rules["soft_layer_weight"] = 0.0
            chemont_layer_rules["mode"] = "disabled_local"
    return {
        "clustering_policy": clustering_policy,
        "class_layer_rules": class_layer_rules,
        "descriptor_layer_rules": descriptor_layer_rules,
        "topology_layer_rules": topology_layer_rules,
        "ionisation_layer_rules": ionisation_layer_rules,
        "toxicophore_layer_rules": toxicophore_layer_rules,
        "chemont_layer_rules": chemont_layer_rules,
    }


def process_clustering_and_adme_core(df: pd.DataFrame, smiles_col: str = 'Standardized SMILES', debug_trace: list[dict] | None = None) -> pd.DataFrame:
    """Deterministic core enrichment used by Phase B1a."""
    df = df.copy()
    if smiles_col not in df.columns or df.empty:
        return df
    _trace_phase_b(debug_trace, "start", "ok", df, f"smiles_col={smiles_col}")

    policies = _load_phase_b_policies()
    clustering_policy = policies["clustering_policy"]
    class_layer_rules = policies["class_layer_rules"]
    descriptor_layer_rules = policies["descriptor_layer_rules"]
    topology_layer_rules = policies["topology_layer_rules"]
    ionisation_layer_rules = policies["ionisation_layer_rules"]
    toxicophore_layer_rules = policies["toxicophore_layer_rules"]
    chemont_layer_rules = policies["chemont_layer_rules"]
    _trace_phase_b(debug_trace, "policy_load", "ok", df, "Loaded review-only clustering policies and layer rules")

    from structural_evidence import add_structural_evidence
    try:
        df = add_structural_evidence(df, smiles_col=smiles_col)
        _trace_phase_b(debug_trace, "structural_evidence", "ok", df, "Descriptor and 3D evidence attached")
    except Exception as e:
        _trace_phase_b(debug_trace, "structural_evidence", "error", df, str(e))
        raise

    try:
        df = flag_alerts(df, smiles_col=smiles_col)
        _trace_phase_b(debug_trace, "structural_alerts", "ok", df, "Alerts flagged")
    except Exception as e:
        _trace_phase_b(debug_trace, "structural_alerts", "error", df, str(e))
        raise

    try:
        df = extract_scaffolds(df, smiles_col=smiles_col)
        _trace_phase_b(debug_trace, "scaffold_extraction", "ok", df, "Scaffolds extracted")
    except Exception as e:
        _trace_phase_b(debug_trace, "scaffold_extraction", "error", df, str(e))
        raise

    try:
        df = run_revalidated_clustering(df, smiles_col=smiles_col)
        _trace_phase_b(
            debug_trace,
            "revalidated_clustering",
            "ok",
            df,
            "Deduplicated, compatibility-gated Morgan/Butina subclusters assigned",
        )
    except Exception as e:
        _trace_phase_b(debug_trace, "consensus_clustering", "error", df, str(e))
        raise

    from domain_evidence import add_domain_evidence
    try:
        df = add_domain_evidence(df, smiles_col=smiles_col)
        _trace_phase_b(debug_trace, "domain_evidence", "ok", df, "Domain evidence attached")
    except Exception as e:
        _trace_phase_b(debug_trace, "domain_evidence", "error", df, str(e))
        raise

    from uncertainty_register import add_uncertainty_register
    try:
        df, uncertainty_register = add_uncertainty_register(df)
        _trace_phase_b(debug_trace, "uncertainty_register", "ok", df, "Member and cluster uncertainty register created")
    except Exception as e:
        _trace_phase_b(debug_trace, "uncertainty_register", "error", df, str(e))
        raise

    try:
        from cluster_reasoning import add_cluster_membership_evidence
        df = add_cluster_membership_evidence(df, smiles_col=smiles_col)
        _trace_phase_b(debug_trace, "cluster_membership_evidence", "ok", df, "Compound-specific deterministic cluster evidence attached")
    except Exception as e:
        _trace_phase_b(debug_trace, "cluster_membership_evidence", "error", df, str(e))
        raise

    try:
        from evidence_ledger import attach_evidence_ledger
        df = attach_evidence_ledger(df)
        endpoint_evidence_ledger = df.attrs.get("endpoint_evidence_ledger")
        _trace_phase_b(debug_trace, "endpoint_ledger", "ok", df, "Endpoint evidence ledger attached")
    except Exception as e:
        print("Endpoint evidence ledger unavailable:", e)
        endpoint_evidence_ledger = None
        _trace_phase_b(debug_trace, "endpoint_ledger", "error", df, str(e))

    df.attrs["uncertainty_register"] = uncertainty_register
    if endpoint_evidence_ledger is not None:
        df.attrs["endpoint_evidence_ledger"] = endpoint_evidence_ledger
    if clustering_policy is not None:
        df.attrs["clustering_policy"] = clustering_policy
    if class_layer_rules is not None:
        df.attrs["class_layer_policy"] = class_layer_rules
    if descriptor_layer_rules is not None:
        df.attrs["descriptor_layer_policy"] = descriptor_layer_rules
    if topology_layer_rules is not None:
        df.attrs["topology_layer_policy"] = topology_layer_rules
    if ionisation_layer_rules is not None:
        df.attrs["ionisation_layer_policy"] = ionisation_layer_rules
    if toxicophore_layer_rules is not None:
        df.attrs["toxicophore_layer_policy"] = toxicophore_layer_rules
    if chemont_layer_rules is not None:
        df.attrs["chemont_layer_policy"] = chemont_layer_rules
    _trace_phase_b(debug_trace, "finish_core", "ok", df, "Phase B1a core pipeline complete")
    return df


def process_clustering_and_adme_optional(df: pd.DataFrame, smiles_col: str = 'Standardized SMILES', debug_trace: list[dict] | None = None) -> pd.DataFrame:
    """Optional enrichment stage used by Phase B1b."""
    df = df.copy()
    if smiles_col not in df.columns or df.empty:
        return df
    _trace_phase_b(debug_trace, "optional_start", "ok", df, f"smiles_col={smiles_col}")

    def _series_has_values(frame: pd.DataFrame, column: str) -> bool:
        if column not in frame.columns:
            return False
        series = frame[column]
        if isinstance(series, pd.DataFrame):
            series = series.iloc[:, 0]
        cleaned = series.apply(lambda value: "" if pd.isna(value) else str(value).strip())
        cleaned = cleaned.replace({"nan": "", "None": "", "NoneType": ""})
        return cleaned.ne("").any()

    core_ready_columns = ["Cluster ID", "Domain", "Decision"]
    if not all(_series_has_values(df, column) for column in core_ready_columns):
        df = process_clustering_and_adme_core(df, smiles_col=smiles_col, debug_trace=debug_trace)

    core_columns_to_preserve = [
        "Cluster ID",
        "Cluster Size",
        "Cluster_Status",
        "Domain",
        "Decision",
        "Domain_Status",
        "Domain_Reasons",
        "Domain_Rule_Version",
        "Nearest_Neighbour_Tanimoto",
        "Cluster_Min_Pairwise_Tanimoto",
        "Cluster_Median_Pairwise_Tanimoto",
        "Cluster_Scaffold_Coverage",
        "Cluster_Representative_Scaffold",
        "Property_Outlier_Flags",
        "Ionisation_Consistency",
        "Toxicophore_Consistency",
        "Uncertainty_Summary",
        "Uncertainty_Rule_Version",
        "Classification Status",
        "Classification Rule ID",
        "Classification Rule Version",
        "Manual Review Flag",
        "Manual Review Reason",
        "Final Classification Record",
        "Taxonomy Path",
        "Taxonomy Hierarchy",
    ]
    preserved_core = {
        column: df[column].copy()
        for column in core_columns_to_preserve
        if column in df.columns
    }

    try:
        from database import SessionLocal
        from chemont_adapter import enrich_with_chemont
        session = SessionLocal()
        try:
            df = enrich_with_chemont(df, session)
        finally:
            session.close()
        chemont_classifications = df.attrs.get("chemont_classifications")
        _trace_phase_b(debug_trace, "chemont_enrichment", "ok", df, "ChemOnt adapter completed")
    except Exception as e:
        print("ChemOnt enrichment unavailable:", e)
        df["ChemOnt_Retrieval_Status"] = "Unavailable: local adapter error"
        chemont_classifications = None
        _trace_phase_b(debug_trace, "chemont_enrichment", "error", df, str(e))

    try:
        from database import SessionLocal
        from epa_ctx_adapter import enrich_with_epa_ctx
        session = SessionLocal()
        try:
            df = enrich_with_epa_ctx(df, session)
        finally:
            session.close()
        from evidence_ledger import attach_evidence_ledger
        df = attach_evidence_ledger(df)
        endpoint_evidence_ledger = df.attrs.get("endpoint_evidence_ledger")
        _trace_phase_b(debug_trace, "epa_ctx", "ok", df, "EPA CTX enrichment refreshed local ledger")
    except Exception as e:
        print("EPA CTX enrichment unavailable:", e)
        df["EPA_CTX_Status"] = "Unavailable: local adapter error"
        df["EPA_CTX_DTXSID"] = ""
        df["EPA_CTX_Source_Version"] = ""
        endpoint_evidence_ledger = df.attrs.get("endpoint_evidence_ledger")
        _trace_phase_b(debug_trace, "epa_ctx", "error", df, str(e))

    try:
        from external_api import fetch_chembl_for_dataframe
        df = fetch_chembl_for_dataframe(df, smiles_col)
        _trace_phase_b(debug_trace, "chembl", "ok", df, "ChEMBL lookup completed")
    except Exception as e:
        print("ChEMBL fetch failed:", e)
        df['ChEMBL_Max_Phase'] = df.get('ChEMBL_Max_Phase', None)
        df['ChEMBL_Targets'] = df.get('ChEMBL_Targets', 'Unavailable')
        _trace_phase_b(debug_trace, "chembl", "error", df, str(e))

    try:
        from cluster_evidence_matrix import build_cluster_evidence_matrix
        uncertainty_register = df.attrs.get("uncertainty_register")
        cluster_evidence_matrix = build_cluster_evidence_matrix(
            df,
            uncertainty_register=uncertainty_register,
            endpoint_evidence_ledger=endpoint_evidence_ledger,
        )
        df.attrs["cluster_evidence_matrix"] = cluster_evidence_matrix
        _trace_phase_b(debug_trace, "cluster_evidence_matrix", "ok", df, "Cluster evidence matrix built")
    except Exception as e:
        print("Cluster evidence matrix unavailable:", e)
        cluster_evidence_matrix = None
        _trace_phase_b(debug_trace, "cluster_evidence_matrix", "error", df, str(e))

    try:
        from xai_engine import generate_batch_cluster_reasoning
        df = generate_batch_cluster_reasoning(df, cluster_col='Cluster ID')
        _trace_phase_b(debug_trace, "xai", "ok", df, "Grounded cluster interpretations generated or restored from cache")
    except Exception as e:
        print("XAI generation failed:", e)
        df['XAI_Status'] = 'Unavailable: batch reasoning failed'
        df['Read_Across_Justification'] = 'Deterministic membership evidence remains the authoritative explanation.'
        _trace_phase_b(debug_trace, "xai", "error", df, str(e))

    try:
        from graph_exporter import generate_cypher_queries
        graph_out = os.path.join(os.path.dirname(__file__), "knowledge_graph.cypher")
        generate_cypher_queries(df, output_file=graph_out)
        _trace_phase_b(debug_trace, "graph_export", "ok", df, "Cypher knowledge graph export written")
    except Exception as e:
        print("Graph export failed:", e)
        _trace_phase_b(debug_trace, "graph_export", "error", df, str(e))

    if cluster_evidence_matrix is not None:
        df.attrs["cluster_evidence_matrix"] = cluster_evidence_matrix
    if endpoint_evidence_ledger is not None:
        df.attrs["endpoint_evidence_ledger"] = endpoint_evidence_ledger
    if chemont_classifications is not None:
        df.attrs["chemont_classifications"] = chemont_classifications
    for column, values in preserved_core.items():
        if len(values) == len(df):
            df[column] = values
    _trace_phase_b(debug_trace, "finish_optional", "ok", df, "Phase B1b optional pipeline complete")
    return df


def process_clustering_and_adme_optional_stage(
    df: pd.DataFrame,
    stage: str,
    smiles_col: str = 'Standardized SMILES',
    debug_trace: list[dict] | None = None,
) -> pd.DataFrame:
    """Run exactly one bounded B1b substage.

    The PDF workflow serializes the optional enrichments so a slow external
    source cannot hold all of B1b open.  Core cluster/domain fields are
    preserved between requests; the caller may therefore send the previous
    JSON response back for the next substage.
    """
    stage_key = str(stage or '').strip().casefold().replace('-', '_')
    aliases = {
        'b1b1': 'chemont', 'b1b_1': 'chemont', 'chemont': 'chemont',
        'b1b2': 'epa_ctx', 'b1b_2': 'epa_ctx', 'epa': 'epa_ctx', 'epa_ctx': 'epa_ctx',
        'b1b3': 'chembl', 'b1b_3': 'chembl', 'chembl': 'chembl',
        'b1b4': 'finalize', 'b1b_4': 'finalize', 'finalize': 'finalize',
    }
    stage_key = aliases.get(stage_key, stage_key)
    if stage_key not in {'chemont', 'epa_ctx', 'chembl', 'finalize'}:
        raise ValueError(f"Unknown B1b substage: {stage}")

    result = df.copy()
    if smiles_col not in result.columns or result.empty:
        return result
    _trace_phase_b(debug_trace, f"{stage_key}_start", "ok", result, f"B1b substage={stage_key}")

    def _series_has_values(frame: pd.DataFrame, column: str) -> bool:
        if column not in frame.columns:
            return False
        series = frame[column]
        if isinstance(series, pd.DataFrame):
            series = series.iloc[:, 0]
        cleaned = series.apply(lambda value: "" if pd.isna(value) else str(value).strip())
        return cleaned.replace({"nan": "", "None": "", "NoneType": ""}).ne("").any()

    if not all(_series_has_values(result, column) for column in ("Cluster ID", "Domain", "Decision")):
        result = process_clustering_and_adme_core(result, smiles_col=smiles_col, debug_trace=debug_trace)

    core_columns = [
        "Cluster ID", "Cluster Size", "Cluster_Status", "Domain", "Decision",
        "Domain_Status", "Domain_Reasons", "Domain_Rule_Version",
        "Nearest_Neighbour_Tanimoto", "Cluster_Min_Pairwise_Tanimoto",
        "Cluster_Median_Pairwise_Tanimoto", "Cluster_Scaffold_Coverage",
        "Cluster_Representative_Scaffold", "Property_Outlier_Flags",
        "Ionisation_Consistency", "Toxicophore_Consistency", "Uncertainty_Summary",
        "Uncertainty_Rule_Version", "Classification Status", "Classification Rule ID",
        "Classification Rule Version", "Manual Review Flag", "Manual Review Reason",
        "Final Classification Record", "Taxonomy Path", "Taxonomy Hierarchy",
    ]
    preserved = {column: result[column].copy() for column in core_columns if column in result.columns}

    if stage_key == 'chemont':
        local_mode = os.environ.get("DETERMINISTIC_LOCAL_MODE", "false").strip().casefold() in {"1", "true", "yes", "on"}
        if local_mode:
            for column in ("ChemOnt_Kingdom", "ChemOnt_Superclass", "ChemOnt_Class", "ChemOnt_Subclass", "ChemOnt_Direct_Parent", "ChemOnt_Molecular_Framework", "ChemOnt_Substituents", "ChemOnt_Classification_Version"):
                if column not in result.columns:
                    result[column] = ""
            result["ChemOnt_Source"] = "ClassyFire ChemOnt API"
            result["ChemOnt_Retrieval_Status"] = "Not requested: deterministic local mode"
            result.attrs["chemont_classifications"] = result[["InChIKey", "Standardized SMILES", "ChemOnt_Kingdom", "ChemOnt_Superclass", "ChemOnt_Class", "ChemOnt_Subclass", "ChemOnt_Direct_Parent", "ChemOnt_Molecular_Framework", "ChemOnt_Substituents", "ChemOnt_Classification_Version", "ChemOnt_Source", "ChemOnt_Retrieval_Status"]].copy()
            _trace_phase_b(debug_trace, "chemont_enrichment", "skipped", result, "ChemOnt disabled in deterministic local mode")
        else:
            from database import SessionLocal
            from chemont_adapter import enrich_with_chemont
            session = SessionLocal()
            try:
                result = enrich_with_chemont(result, session)
            finally:
                session.close()
            _trace_phase_b(debug_trace, "chemont_enrichment", "ok", result, "ChemOnt substage completed")
    elif stage_key == 'epa_ctx':
        local_mode = os.environ.get("DETERMINISTIC_LOCAL_MODE", "false").strip().casefold() in {"1", "true", "yes", "on"}
        if local_mode:
            if "EPA_CTX_Status" not in result.columns:
                result["EPA_CTX_Status"] = "Not requested: deterministic local mode"
            else:
                result["EPA_CTX_Status"] = "Not requested: deterministic local mode"
            for column in ("EPA_CTX_DTXSID", "EPA_CTX_Source_Version"):
                if column not in result.columns:
                    result[column] = ""
            _trace_phase_b(debug_trace, "epa_ctx", "skipped", result, "EPA CTX disabled in deterministic local mode")
        else:
            from database import SessionLocal
            from epa_ctx_adapter import enrich_with_epa_ctx
            session = SessionLocal()
            try:
                result = enrich_with_epa_ctx(result, session)
            finally:
                session.close()
            _trace_phase_b(debug_trace, "epa_ctx", "ok", result, "EPA CTX substage completed")
    elif stage_key == 'chembl':
        from external_api import fetch_chembl_for_dataframe
        result = fetch_chembl_for_dataframe(result, smiles_col)
        _trace_phase_b(debug_trace, "chembl", "ok", result, "ChEMBL substage completed")
    else:
        from cluster_evidence_matrix import build_cluster_evidence_matrix
        from evidence_ledger import attach_evidence_ledger
        endpoint_ledger = result.attrs.get("endpoint_evidence_ledger")
        result = attach_evidence_ledger(result)
        endpoint_ledger = result.attrs.get("endpoint_evidence_ledger", endpoint_ledger)
        matrix = build_cluster_evidence_matrix(
            result,
            uncertainty_register=result.attrs.get("uncertainty_register"),
            endpoint_evidence_ledger=endpoint_ledger,
        )
        result.attrs["cluster_evidence_matrix"] = matrix
        try:
            from xai_engine import generate_batch_cluster_reasoning
            # The user-approved AI batch runs only after B1b-4 has frozen the
            # deterministic evidence used by the final workbook.
            result = generate_batch_cluster_reasoning(result, cluster_col='Cluster ID', enabled=True)
            _trace_phase_b(debug_trace, "xai", "ok", result, "Grounded cluster interpretations generated or restored from cache")
        except Exception as exc:
            result['XAI_Status'] = 'Unavailable: batch reasoning failed'
            result['Read_Across_Justification'] = 'Deterministic membership evidence remains the authoritative explanation.'
            _trace_phase_b(debug_trace, "xai", "error", result, str(exc))
        try:
            from graph_exporter import generate_cypher_queries
            graph_out = os.path.join(os.path.dirname(__file__), "knowledge_graph.cypher")
            generate_cypher_queries(result, output_file=graph_out)
            _trace_phase_b(debug_trace, "graph_export", "ok", result, "Cypher knowledge graph export written")
        except Exception as exc:
            _trace_phase_b(debug_trace, "graph_export", "error", result, str(exc))

    for column, values in preserved.items():
        if len(values) == len(result):
            result[column] = values
    _trace_phase_b(debug_trace, f"{stage_key}_finish", "ok", result, f"B1b substage {stage_key} complete")
    return result
