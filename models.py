from sqlalchemy import Column, Integer, String, Float, DateTime, JSON, ForeignKey, Text, UniqueConstraint, Boolean
from sqlalchemy.sql import func
from database import Base

class AssessmentRecord(Base):
    __tablename__ = "assessment_records"

    id = Column(Integer, primary_key=True, index=True)
    
    # Input data
    input_compound_name = Column(String, index=True)
    input_cas_number = Column(String, index=True)
    input_molecular_formula = Column(String)
    input_concentration = Column(Float)
    input_concentration_unit = Column(String)
    
    # Resolution (Identity)
    resolved_name = Column(String)
    resolved_cas = Column(String)
    resolved_smiles = Column(String)
    resolved_inchi = Column(String)
    resolved_inchikey = Column(String)
    resolved_formula = Column(String)
    resolved_mw = Column(Float)
    resolution_status = Column(String)  # Exact Match, Partial Match, Ambiguous Match, Unresolved
    source_database = Column(String)
    
    # Structure Standardization
    standardized_smiles = Column(String)
    
    # Analysis
    functional_groups = Column(JSON)  # List of identified FGs
    parent_moiety = Column(String)
    chemical_class = Column(String)
    morgan_fingerprint = Column(String)  # We might store it as a string representation or bytes
    
    # Similarity
    analogs = Column(JSON) # List of dictionaries with Tanimoto scores
    
    # Rationale and Confidence
    grouping_rationale = Column(String)
    confidence_score = Column(String) # High, Medium, Low
    
    # Audit
    created_at = Column(DateTime(timezone=True), server_default=func.now())

class ChemicalDatabase(Base):
    """
    Internal database for analog search (Step 8)
    """
    __tablename__ = "chemical_database"
    
    id = Column(Integer, primary_key=True, index=True)
    cas_number = Column(String, unique=True, index=True)
    name = Column(String, index=True)
    smiles = Column(String)
    standardized_smiles = Column(String)
    morgan_fingerprint = Column(String) # serialized binary or hex string
    parent_moiety = Column(String)
    functional_groups = Column(JSON)


class EndpointEvidenceRecord(Base):
    """Immutable Stage 4 endpoint-evidence provenance record."""
    __tablename__ = "endpoint_evidence_records"
    __table_args__ = (UniqueConstraint("raw_record_hash", name="uq_endpoint_evidence_raw_hash"),)

    id = Column(Integer, primary_key=True, index=True)
    evidence_id = Column(String, unique=True, index=True, nullable=False)
    identity_key = Column(String, index=True, nullable=False)
    compound_mapping_method = Column(String, nullable=False)
    endpoint = Column(String, nullable=False)
    result = Column(String, nullable=False)
    study_type = Column(String, nullable=False)
    species = Column(String, nullable=False)
    route = Column(String, nullable=False)
    dose = Column(String)
    dose_unit = Column(String)
    duration = Column(String)
    source = Column(String, nullable=False)
    source_record = Column(String, nullable=False)
    source_url = Column(String)
    source_database_version = Column(String)
    reliability = Column(String, nullable=False)
    evidence_status = Column(String, nullable=False)
    raw_record_hash = Column(String, nullable=False)
    vocabulary_version = Column(String, nullable=False)
    retrieved_at = Column(String)
    imported_at = Column(String, nullable=False)
    imported_by = Column(String, nullable=False)
    reviewer_status = Column(String, nullable=False)
    supersedes_evidence_id = Column(String)
    notes = Column(Text)


class EpaCtxCacheEntry(Base):
    """Stage 5 local cache of read-only EPA CTX responses."""
    __tablename__ = "epa_ctx_cache_entries"
    __table_args__ = (UniqueConstraint("dtxsid", "source_version", "resource_path", name="uq_epa_ctx_cache_key"),)

    id = Column(Integer, primary_key=True, index=True)
    dtxsid = Column(String, nullable=False, index=True)
    source_version = Column(String, nullable=False)
    resource_path = Column(String, nullable=False)
    retrieved_at = Column(String, nullable=False)
    http_status = Column(Integer, nullable=False)
    response_hash = Column(String, nullable=False)
    response_json = Column(JSON, nullable=False)


class ChemontClassificationCacheEntry(Base):
    """Stage 7 cache of ClassyFire/ChemOnt structural classifications.

    This contains only a public structure identifier and returned structural
    taxonomy.  It is not endpoint evidence and does not express a safety
    conclusion.
    """
    __tablename__ = "chemont_classification_cache_entries"
    __table_args__ = (UniqueConstraint("inchikey", "source", name="uq_chemont_cache_key"),)

    id = Column(Integer, primary_key=True, index=True)
    inchikey = Column(String, nullable=False, index=True)
    standardized_smiles = Column(String)
    source = Column(String, nullable=False)
    classification_version = Column(String)
    retrieved_at = Column(String, nullable=False)
    http_status = Column(Integer, nullable=False)
    response_hash = Column(String, nullable=False)
    response_json = Column(JSON, nullable=False)


class StructureClassificationRecord(Base):
    """Immutable audit row for the local structure classification decision."""
    __tablename__ = "structure_classification_records"

    id = Column(Integer, primary_key=True, index=True)
    identity_key = Column(String, index=True, nullable=False)
    input_smiles = Column(String)
    standardized_smiles = Column(String)
    classification_input_source = Column(String, nullable=False)
    corrected_chemical_class = Column(String, nullable=False)
    final_classification_record = Column(Text, nullable=False)
    manual_review_flag = Column(Boolean, nullable=False)
    manual_review_reason = Column(Text)
    classification_status = Column(String, nullable=False)
    classification_rule_id = Column(String, nullable=False)
    classification_rule_version = Column(String, nullable=False)
    detected_feature_profile = Column(JSON, nullable=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())
