from pydantic import BaseModel, Field
from typing import List, Optional, Dict, Any

# --- Input Schemas ---
class AssessmentRequest(BaseModel):
    compound_name: Optional[str] = None
    cas_number: Optional[str] = None
    molecular_formula: Optional[str] = None
    concentration: Optional[float] = None
    concentration_unit: Optional[str] = None
    resolved_name: Optional[str] = None
    original_smiles: Optional[str] = None
    standardized_smiles: Optional[str] = None
    inchikey: Optional[str] = None
    inchi: Optional[str] = None

# --- Output Schemas ---
class IdentityInfo(BaseModel):
    name: str
    iupac_name: str
    cid: str
    synonyms: List[str]
    cas: str
    smiles: str
    inchi: str
    inchikey: str
    formula: str
    molecular_weight: float
    source: str
    status: str

class StructureInfo(BaseModel):
    original_smiles: str
    standardized_smiles: str
    image_base64: str

class FunctionalGroup(BaseModel):
    functional_group: str
    matched_atoms: List[int]
    SMARTS_pattern: str

class AnalogInfo(BaseModel):
    chemical_name: str
    cas: str
    smiles: str
    tanimoto_score: float
    classification: str # Very Strong Analog, Strong Analog, Moderate Analog, Weak Analog

class AnalogMatch(BaseModel):
    name: str
    smiles: str
    similarity_score: float
    confidence_tier: str
    image_base64: str

class AssessmentResponse(BaseModel):
    identity: Dict[str, Any]
    structure: Dict[str, Any]
    functional_groups: List[FunctionalGroup]
    primary_functional_group: str = ""
    secondary_functional_groups: List[FunctionalGroup] = Field(default_factory=list)
    parent_moiety: str
    chemical_class: str
    corrected_chemical_class: str = "Unclassified"
    final_classification_record: str = ""
    manual_review_flag: bool = False
    manual_review_reason: str = ""
    classification_input_source: str = ""
    classification_structure_smiles: str = ""
    classification_standardization_version: str = ""
    classification_stereochemistry_status: str = ""
    classification_status: str = ""
    classification_rule_id: str = ""
    classification_rule_version: str = ""
    classification_scope: str = ""
    classification_review_recommended: bool = False
    classification_review_recommendation: str = ""
    classification_suggested_action: str = ""
    detected_feature_profile: List[str] = Field(default_factory=list)
    taxonomy_path: str = ""
    taxonomy_path_steps: List[str] = Field(default_factory=list)
    matched_categories: List[str] = Field(default_factory=list)
    direct_parent: str = ""
    taxonomy_dictionary_version: str = ""
    taxonomy_hierarchy: Dict[str, str] = Field(default_factory=dict)
    ontology_compatibility: Dict[str, str] = Field(default_factory=dict)
    topology_profile: Dict[str, Any] = Field(default_factory=dict)
    feature_profile: Dict[str, Any] = Field(default_factory=dict)
    structural_evidence: Dict[str, Any] = Field(default_factory=dict)
    scaffold: str = ""
    cramer_class: str = ""
    ttc_limit: str = ""
    ttc_value: float = 0.0
    ttc_unit: str = ""
    cramer_rule_version: str = ""
    cramer_path: str = ""
    cramer_evidence: str = ""
    toxtree_results: Dict[str, Any] = Field(default_factory=dict)
    chemont_classification: Dict[str, Any] = Field(default_factory=dict)
    cluster_context: Dict[str, Any] = Field(default_factory=dict)
    final_classification_details: Dict[str, Any] = Field(default_factory=dict)
    classification_rationale: str = ""
    analogs: List[AnalogMatch]
    confidence: str
    rationale: str
