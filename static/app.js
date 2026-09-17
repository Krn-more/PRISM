function showGlobalLoader(text = "Processing...") {
    document.getElementById('global-loader-text').innerText = text;
    document.getElementById('global-loader').style.display = 'flex';
}
function hideGlobalLoader() {
    document.getElementById('global-loader').style.display = 'none';
}

function switchTab(tabId) {
    document.querySelectorAll('.tab-content').forEach(el => el.style.display = 'none');
    document.querySelectorAll('.tab-btn').forEach(el => el.classList.remove('active'));
    document.getElementById('tab-' + tabId).style.display = 'block';
    const activeButton = document.querySelector(`.tab-btn[onclick="switchTab('${tabId}')"]`);
    if (activeButton) activeButton.classList.add('active');
    const layout = document.getElementById('assessmentLayout');
    if (layout) layout.classList.toggle('pdf-mode', tabId === 'pdf');
    if (tabId === 'pdf') updatePdfStepper();
}

function updatePdfStepper() {
    const stepIds = ['pdf-step0', 'pdf-step1', 'pdf-step234', 'pdf-step5', 'pdf-step55', 'pdf-step6'];
    const activeIndex = stepIds.findIndex(id => {
        const el = document.getElementById(id);
        return el && getComputedStyle(el).display !== 'none';
    });
    if (activeIndex < 0) return;
    const items = document.querySelectorAll('.pdf-stepper-item');
    items.forEach((item, index) => {
        item.classList.toggle('active', index === activeIndex);
        item.classList.toggle('complete', index < activeIndex);
        if (index === activeIndex) item.setAttribute('aria-current', 'step');
        else item.removeAttribute('aria-current');
    });
    const state = document.getElementById('pdfWorkflowState');
    if (state) state.textContent = activeIndex === 5 ? 'Complete' : `Step ${activeIndex + 1} of 6`;
}
// Initialize lucide icons
lucide.createIcons();

document.getElementById('assessmentForm').addEventListener('submit', async (e) => {
    e.preventDefault();
    
    // UI states
    const btn = document.getElementById('analyzeBtn');
    const resultsArea = document.getElementById('resultsArea');
    const loader = document.getElementById('loader');
    const resultsContent = document.getElementById('resultsContent');
    
    btn.disabled = true;
    btn.innerHTML = '<i data-lucide="loader" class="spinner-inline"></i> Processing...';
    lucide.createIcons();
    
    resultsArea.style.display = 'block';
    loader.style.display = 'flex';
    resultsContent.style.display = 'none';

    // Gather data
    const payload = {
        compound_name: document.getElementById('compound_name').value.trim(),
        cas_number: document.getElementById('cas_number').value.trim(),
        molecular_formula: document.getElementById('formula').value.trim(),
        concentration: parseFloat(document.getElementById('concentration').value) || 0,
        concentration_unit: document.getElementById('unit').value,
        resolved_name: document.getElementById('resolved_name').value.trim(),
        original_smiles: document.getElementById('original_smiles').value.trim(),
        standardized_smiles: document.getElementById('standardized_smiles').value.trim(),
        inchikey: document.getElementById('inchikey').value.trim(),
        inchi: document.getElementById('inchi').value.trim()
    };

    try {
        const response = await fetch('/assess', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payload)
        });

        if (!response.ok) {
            throw new Error(`Error: ${response.status} ${response.statusText}`);
        }

        const data = await response.json();
        
        // Populate UI
        document.getElementById('res-name').innerText = data.identity.name || payload.compound_name;
        
        const confBadge = document.getElementById('res-confidence');
        confBadge.innerText = data.confidence + ' Confidence';
        confBadge.style.color = data.confidence === 'High' ? 'var(--success)' : (data.confidence === 'Medium' ? 'var(--warning)' : 'var(--danger)');
        confBadge.style.backgroundColor = data.confidence === 'High' ? 'rgba(16, 185, 129, 0.1)' : (data.confidence === 'Medium' ? 'rgba(245, 158, 11, 0.1)' : 'rgba(239, 68, 68, 0.1)');
        confBadge.style.borderColor = data.confidence === 'High' ? 'rgba(16, 185, 129, 0.2)' : (data.confidence === 'Medium' ? 'rgba(245, 158, 11, 0.2)' : 'rgba(239, 68, 68, 0.2)');

        // Step 1: Identity & Structure
        let synHtml = data.identity.synonyms && data.identity.synonyms.length > 0 ? 
            `<div class="data-item" style="grid-column: 1 / -1;"><span class="data-label">Synonyms</span><span class="data-val" style="color: #94a3b8;">${data.identity.synonyms.join('; ')}</span></div>` : '';

        let imgHtml = data.structure.image_base64 ? 
            `<img src="data:image/png;base64,${data.structure.image_base64}" style="width: 200px; height: 200px; object-fit: contain; border-radius: 8px; box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1);" />` : 
            `<div style="width: 200px; height: 200px; background: #1e293b; border-radius: 8px; display:flex; align-items:center; justify-content:center; color:#64748b;">No Image</div>`;

        document.getElementById('identity-data').innerHTML = `
            <div style="grid-column: 1 / -1; display: flex; gap: 24px; align-items: flex-start;">
                <div style="flex-shrink: 0; background: #fff; padding: 10px; border-radius: 12px; box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1);">
                    ${imgHtml}
                </div>
                <div style="flex-grow: 1; display: grid; grid-template-columns: repeat(2, 1fr); gap: 1.5rem;">
                    <div class="data-item"><span class="data-label">Source Database</span><span class="data-val" style="color: #60a5fa;">${data.identity.source} ${data.identity.cid ? '(CID: ' + data.identity.cid + ')' : ''}</span></div>
                    <div class="data-item"><span class="data-label">CASRN</span><span class="data-val">${data.identity.cas || '-'}</span></div>
                    <div class="data-item" style="grid-column: 1 / -1;"><span class="data-label">IUPAC Name</span><span class="data-val">${data.identity.iupac_name || '-'}</span></div>
                    <div class="data-item"><span class="data-label">Formula</span><span class="data-val">${data.identity.formula || '-'}</span></div>
                    <div class="data-item"><span class="data-label">Mol. Weight</span><span class="data-val">${data.identity.molecular_weight ? data.identity.molecular_weight + ' g/mol' : '-'}</span></div>
                    ${synHtml}
                    <div class="data-item" style="grid-column: 1 / -1;"><span class="data-label">Original SMILES</span><span class="data-val" style="word-break: break-all; color: #cbd5e1;">${data.structure.original_smiles || '-'}</span></div>
                    <div class="data-item" style="grid-column: 1 / -1;"><span class="data-label">Standardized SMILES</span><span class="data-val" style="word-break: break-all;">${data.structure.standardized_smiles || '-'}</span></div>
                    <div class="data-item" style="grid-column: 1 / -1;"><span class="data-label">InChIKey</span><span class="data-val" style="color: #cbd5e1;">${data.identity.inchikey || '-'}</span></div>
                    <div class="data-item" style="grid-column: 1 / -1;"><span class="data-label">InChI</span><span class="data-val" style="word-break: break-all; color: #94a3b8;">${data.identity.inchi || '-'}</span></div>
                </div>
            </div>
        `;

        // Step 2: Topology & Functional Groups
        const featureProfileText = (data.detected_feature_profile && data.detected_feature_profile.length)
            ? data.detected_feature_profile.join('; ')
            : 'None';
        const taxonomyHierarchy = data.taxonomy_hierarchy || {};
        const ontologyCompatibility = data.ontology_compatibility || {};
        const topologyProfile = data.topology_profile || {};
        const structuralEvidence = data.structural_evidence || {};
        const primaryFunctionalGroup = data.primary_functional_group || data.corrected_chemical_class || 'Not available';
        const secondaryFunctionalGroups = Array.isArray(data.secondary_functional_groups)
            ? data.secondary_functional_groups
                .map(fg => fg.functional_group || fg.Functional_Group || fg.name || '')
                .filter(Boolean)
            : [];
        const taxonomyFields = [
            ['Parent Class', taxonomyHierarchy['Parent Class'] || (data.taxonomy_path_steps && data.taxonomy_path_steps[0]) || 'Not available'],
            ['Functional Group', taxonomyHierarchy['Functional Group'] || (data.taxonomy_path_steps && data.taxonomy_path_steps[1]) || data.corrected_chemical_class || 'Not available'],
            ['Subclass', taxonomyHierarchy['Subclass'] || (data.taxonomy_path_steps && data.taxonomy_path_steps[2]) || data.corrected_chemical_class || 'Not available'],
            ['Structural Type', taxonomyHierarchy['Structural Type'] || (data.taxonomy_path_steps && data.taxonomy_path_steps[3]) || data.corrected_chemical_class || 'Not available']
        ];
        const ontologyFields = [
            ['Mapping Status', ontologyCompatibility.mapping_status || 'Not mapped'],
            ['Ontology Version', ontologyCompatibility.ontology_mapping_version || 'Not available'],
            ['Kingdom', ontologyCompatibility.kingdom || 'Not available'],
            ['Superclass', ontologyCompatibility.superclass || 'Not available'],
            ['Class', ontologyCompatibility.class || 'Not available'],
            ['Subclass', ontologyCompatibility.subclass || 'Not available'],
            ['Direct Parent', ontologyCompatibility.direct_parent || 'Not available'],
            ['Mapping Evidence', ontologyCompatibility.mapping_evidence || 'No approved local ontology mapping for this structural rule', 2]
        ];
        const topologyFields = [
            ['Topology Class', topologyProfile['Topology Class'] || 'Not available'],
            ['Topology Modifiers', topologyProfile['Topology Modifiers'] || 'None'],
            ['Ring System', topologyProfile['Ring System'] || 'Not available'],
            ['Ring Profile', topologyProfile['Ring Profile'] || (data.final_classification_details && data.final_classification_details['Ring Profile']) || 'Not available']
        ];
        const descriptorFields = [
            ['MW', structuralEvidence.MW ?? data.identity?.molecular_weight ?? 'Not available'],
            ['LogP', structuralEvidence.LogP ?? 'Not available'],
            ['TPSA', structuralEvidence.TPSA ?? 'Not available'],
            ['Min EState', structuralEvidence.Min_EState ?? 'Not available'],
            ['Max EState', structuralEvidence.Max_EState ?? 'Not available'],
            ['HBD', structuralEvidence.HBD ?? 'Not available'],
            ['HBA', structuralEvidence.HBA ?? 'Not available'],
            ['Rotatable Bonds', structuralEvidence.Rotatable_Bonds ?? 'Not available'],
            ['Formal Charge', structuralEvidence.Formal_Charge ?? 'Not available'],
            ['Ring Count', structuralEvidence.Ring_Count ?? 'Not available'],
            ['Aromatic Ring Count', structuralEvidence.Aromatic_Ring_Count ?? 'Not available'],
            ['Fraction CSP3', structuralEvidence.Fraction_CSP3 ?? 'Not available'],
            ['Heavy Atom Count', structuralEvidence.Heavy_Atom_Count ?? 'Not available'],
            ['Molecular Refractivity', structuralEvidence.Molecular_Refractivity ?? 'Not available'],
            ['Ionisation Indicator', structuralEvidence.Ionisation_Indicator || 'Not available']
        ];
        const descriptor3DFields = [
            ['3D Asphericity', structuralEvidence['3D_Asphericity'] ?? 'Not available'],
            ['3D PMI1', structuralEvidence['3D_PMI1'] ?? 'Not available'],
            ['3D PMI2', structuralEvidence['3D_PMI2'] ?? 'Not available'],
            ['3D PMI3', structuralEvidence['3D_PMI3'] ?? 'Not available'],
            ['3D Radius of Gyration', structuralEvidence['3D_RadiusOfGyration'] ?? 'Not available']
        ];
        const alertFields = [
            ['Structural Evidence Status', structuralEvidence.Structural_Evidence_Status || 'Not available'],
            ['Structural Evidence Version', structuralEvidence.Structural_Evidence_Version || 'Not available'],
            ['Toxicophore Profile', structuralEvidence.Toxicophore_Profile || 'None'],
            ['All Structural Alerts', structuralEvidence.All_Structural_Alerts || 'None'],
            ['Alerts', structuralEvidence.Alerts || 'None'],
            ['Scaffold', data.scaffold || 'Not available']
        ];
        const formatMirrorValue = (value, fallback = 'Not available') => {
            if (value === null || value === undefined || value === '') return fallback;
            if (Array.isArray(value)) return value.length ? value.join('; ') : fallback;
            if (typeof value === 'object') {
                const entries = Object.entries(value).filter(([, v]) => v !== null && v !== undefined && v !== '');
                return entries.length ? entries.map(([k, v]) => `${k}: ${v}`).join(' | ') : fallback;
            }
            return String(value);
        };
        const taxonomyHierarchyText = formatMirrorValue(taxonomyHierarchy);
        const topologyProfileText = formatMirrorValue(topologyProfile);
        const secondaryFunctionalGroupText = secondaryFunctionalGroups.length ? secondaryFunctionalGroups.join('; ') : 'Not available';
        const taxonomyPathText = data.taxonomy_path || (
            Array.isArray(data.taxonomy_path_steps) && data.taxonomy_path_steps.length
                ? data.taxonomy_path_steps.join(' → ')
                : 'Not available'
        );
        const workflowStep = (title, subtitle, color, fields, columns = 2) => {
            const rowsPerColumn = Math.ceil(fields.length / columns);
            const columnSlices = Array.from({ length: columns }, (_, idx) => fields.slice(idx * rowsPerColumn, (idx + 1) * rowsPerColumn));
            return `
            <section style="padding: 14px; background: rgba(2, 6, 23, 0.72); border: 1px solid #334155; border-radius: 14px;">
                <div style="display: flex; align-items: flex-start; justify-content: space-between; gap: 10px; margin-bottom: 10px; flex-wrap: wrap;">
                    <div>
                        <div style="font-size: 0.76rem; font-weight: 800; letter-spacing: 0.04em; color: ${color}; text-transform: uppercase;">${title}</div>
                        ${subtitle ? `<div style="color: #94a3b8; font-size: 0.86rem; line-height: 1.45; margin-top: 4px;">${subtitle}</div>` : ''}
                    </div>
                </div>
                <div style="display: grid; grid-template-columns: repeat(${columns}, minmax(0, 1fr)); gap: 10px;">
                    ${columnSlices.map(col => `
                        <div style="display: flex; flex-direction: column; gap: 8px;">
                            ${col.map(([label, value]) => `
                                <div style="padding: 10px 12px; background: rgba(15,23,42,0.72); border: 1px solid #334155; border-radius: 10px;">
                                    <div style="font-size: 0.72rem; letter-spacing: 0.03em; text-transform: uppercase; color: #94a3b8; margin-bottom: 4px;">${label}</div>
                                    <div style="color: #e2e8f0; line-height: 1.45; word-break: break-word; font-size: 0.92rem; font-weight: 600;">${formatMirrorValue(value)}</div>
                                </div>
                            `).join('')}
                        </div>
                    `).join('')}
                </div>
            </section>
        `};
        const fixedAuditTemplate = (title, subtitle, color, fields, columns = 4) => `
            <section style="padding: 14px; background: rgba(2, 6, 23, 0.72); border: 1px solid #334155; border-radius: 14px;">
                <div style="display: flex; align-items: flex-start; justify-content: space-between; gap: 10px; margin-bottom: 10px; flex-wrap: wrap;">
                    <div>
                        <div style="font-size: 0.76rem; font-weight: 800; letter-spacing: 0.04em; color: ${color}; text-transform: uppercase;">${title}</div>
                        ${subtitle ? `<div style="color: #94a3b8; font-size: 0.85rem; line-height: 1.4; margin-top: 4px;">${subtitle}</div>` : ''}
                    </div>
                </div>
                <div style="display: grid; grid-template-columns: repeat(${columns}, minmax(0, 1fr)); gap: 8px;">
                    ${fields.map(([label, value, span = 1]) => `
                        <div style="grid-column: span ${span}; padding: 10px 12px; background: rgba(15,23,42,0.72); border: 1px solid #334155; border-radius: 10px; min-width: 0;">
                            <div style="font-size: 0.7rem; letter-spacing: 0.04em; text-transform: uppercase; color: #94a3b8; margin-bottom: 4px;">${label}</div>
                            <div style="color: #e2e8f0; line-height: 1.45; word-break: break-word; font-size: 0.92rem; font-weight: 600;">${formatMirrorValue(value)}</div>
                        </div>
                    `).join('')}
                </div>
            </section>
        `;
        const compactAuditTemplate = (title, subtitle, color, fields, columns = 3) => `
            <section style="padding: 12px; background: rgba(2, 6, 23, 0.72); border: 1px solid #334155; border-radius: 14px;">
                <div style="display: flex; align-items: flex-start; justify-content: space-between; gap: 8px; margin-bottom: 8px; flex-wrap: wrap;">
                    <div>
                        <div style="font-size: 0.74rem; font-weight: 800; letter-spacing: 0.04em; color: ${color}; text-transform: uppercase;">${title}</div>
                        ${subtitle ? `<div style="color: #94a3b8; font-size: 0.82rem; line-height: 1.35; margin-top: 3px;">${subtitle}</div>` : ''}
                    </div>
                </div>
                <div style="display: grid; grid-template-columns: repeat(${columns}, minmax(0, 1fr)); gap: 7px;">
                    ${fields.map(([label, value, span = 1]) => `
                        <div style="grid-column: span ${span}; padding: 8px 10px; background: rgba(15,23,42,0.72); border: 1px solid #334155; border-radius: 10px; min-width: 0;">
                            <div style="font-size: 0.66rem; letter-spacing: 0.04em; text-transform: uppercase; color: #94a3b8; margin-bottom: 3px;">${label}</div>
                            <div style="color: #e2e8f0; line-height: 1.35; word-break: break-word; font-size: 0.88rem; font-weight: 600;">${formatMirrorValue(value)}</div>
                        </div>
                    `).join('')}
                </div>
            </section>
        `;
        let topoHtml = `
            <div style="display: flex; flex-wrap: wrap; gap: 8px; margin-bottom: 14px;">
                <div class="tag" style="background: rgba(99,102,241,0.15); border-color: rgba(99,102,241,0.3); color: #c4b5fd;">PRISM Classification: ${data.corrected_chemical_class || 'Unclassified'}</div>
                <div class="tag" style="background: rgba(14,165,233,0.15); border-color: rgba(14,165,233,0.3); color: #7dd3fc;">Rule: ${data.classification_rule_id || '-'}</div>
                <div class="tag" style="background: rgba(168,85,247,0.15); border-color: rgba(168,85,247,0.3); color: #d8b4fe;">Version: ${data.classification_rule_version || '-'}</div>
                <div class="tag" style="background: rgba(34,197,94,0.15); border-color: rgba(34,197,94,0.3); color: #86efac;">Status: ${data.classification_status || '-'}</div>
                <div class="tag" style="background: rgba(34,197,94,0.15); border-color: rgba(34,197,94,0.3); color: #86efac;">Source: ${data.classification_input_source || '-'}</div>
                ${data.manual_review_flag ? '<div class="tag" style="background: rgba(239,68,68,0.15); border-color: rgba(239,68,68,0.3); color: #fca5a5;">Manual Review</div>' : ''}
            </div>
            <div style="display: flex; flex-direction: column; gap: 12px; margin-top: 8px;">
                ${fixedAuditTemplate(
                    'Step 3 · PRISM classification decision',
                    'Deterministic decision from the standardized chemical structure.',
                    '#c4b5fd',
                    [
                        ['PRISM Classification', data.corrected_chemical_class || 'Unclassified'],
                        ['Classification Status', data.classification_status || '-'],
                        ['Classification Scope', data.classification_scope || '-'],
                        ['Classification Rule ID', data.classification_rule_id || '-'],
                        ['Classification Rule Version', data.classification_rule_version || '-'],
                        ['Manual Review Flag', data.manual_review_flag ? 'Yes' : 'No'],
                        ['Manual Review Reason', data.manual_review_reason || 'None'],
                        ['Detected Feature Profile', featureProfileText],
                        ['Confidence', data.confidence || '-'],
                        ['Tanimoto Analogs Found', Array.isArray(data.analogs) ? data.analogs.length : 0]
                    ],
                    2
                )}
                ${fixedAuditTemplate(
                    'Step 4 · Functional group and topology interpretation',
                    'Fixed hierarchy and topology audit template.',
                    '#a5b4fc',
                    [
                        ['Parent Class', taxonomyHierarchy['Parent Class'] || (data.taxonomy_path_steps && data.taxonomy_path_steps[0]) || 'Not available'],
                        ['Primary Functional Group', primaryFunctionalGroup],
                        ['Secondary Functional Groups', secondaryFunctionalGroupText],
                        ['Subclass', taxonomyHierarchy['Subclass'] || (data.taxonomy_path_steps && data.taxonomy_path_steps[2]) || data.corrected_chemical_class || 'Not available'],
                    ],
                    4
                )}
                <section style="padding: 12px 14px; background: rgba(2, 6, 23, 0.72); border: 1px solid #334155; border-radius: 14px; margin-top: 10px;">
                    <div style="color: #e2e8f0; line-height: 1.55; word-break: break-word; font-size: 0.92rem; font-weight: 600;">${taxonomyPathText}</div>
                </section>
                ${fixedAuditTemplate(
                    'Ontology compatibility (supplementary)',
                    'Approved local mapping for selected rules. It does not replace the local structural decision or affect clustering.',
                    '#67e8f9',
                    ontologyFields,
                    2
                )}
                ${fixedAuditTemplate(
                    'Final Classification Record',
                    'Compact audit closure for validation and traceability.',
                    '#f8fafc',
                    [
                        ['Decision', (data.final_classification_details && data.final_classification_details['Primary Class']) || data.corrected_chemical_class || '-'],
                        ['Rule / Source', `${(data.final_classification_details && data.final_classification_details['Rule ID']) || data.classification_rule_id || '-'} · ${(data.final_classification_details && data.final_classification_details['Rule Source']) || data.classification_input_source || '-'}`],
                        ['Review', (data.final_classification_details && data.final_classification_details['Review']) || (data.manual_review_flag ? 'Manual review required' : 'No manual review required')],
                        ['Evidence', (data.final_classification_details && data.final_classification_details['Feature Evidence']) || featureProfileText || '-'],
                        ['Topology', (data.final_classification_details && data.final_classification_details['Topology Profile']) || (topologyProfile['Topology Class'] || '-')],
                        ['Ring Profile', (data.final_classification_details && data.final_classification_details['Ring Profile']) || '-']
                    ],
                    2
                )}
                ${workflowStep(
                    'Step 5 · Physicochemical descriptor generation',
                    'Workbook-style descriptor values used for comparison and validation.',
                    '#93c5fd',
                    descriptorFields,
                    2
                )}
                ${fixedAuditTemplate(
                    'Step 5B · 3D descriptor detail',
                    'Dedicated 3D descriptor panel for direct Excel parity checks.',
                    '#7dd3fc',
                    descriptor3DFields,
                    2
                )}
                ${workflowStep(
                    'Step 6 · Toxicophore and structural alert review',
                    'Hazard-relevant substructures, scaffold context, and screening limits.',
                    '#fbbf24',
                    alertFields,
                    2
                )}
                <div style="padding: 12px 14px; background: rgba(2, 6, 23, 0.82); border: 1px solid #334155; border-radius: 10px; color: #cbd5e1; line-height: 1.55;">
                    Validation note: fixed audit template for Steps 3–6. We validate the decision, hierarchy, descriptors, and alerts here first, then copy the same field structure and labels into the Excel workflow.
                </div>
            </div>`;
        document.getElementById('topology-data').innerHTML = topoHtml;

        // Transition views
        loader.style.display = 'none';
        resultsContent.style.display = 'block';

    } catch (err) {
        alert("Pipeline Error: " + err.message);
        loader.style.display = 'none';
    } finally {
        btn.disabled = false;
        btn.innerHTML = '<i data-lucide="zap"></i> Run Assessment Pipeline';
        lucide.createIcons();
    }
});

let currentFilename = null;
let currentExtractedData = null;

function resetPdfView() {
    document.getElementById('pdf-step1').style.display = 'none';
    document.getElementById('pdf-step234').style.display = 'none';
    document.getElementById('pdf-step5').style.display = 'none';
    document.getElementById('pdf-step6').style.display = 'none';
    document.getElementById('pdf-step0').style.display = 'block';
    currentFilename = null;
    currentExtractedData = null;
    currentResolvedData = null;
    phaseBRunId = null;
    phaseB1aResolvedData = null;
    phaseB1bResolvedData = null;
    phaseB1b1ResolvedData = null;
    phaseB1b2ResolvedData = null;
    phaseB1b3ResolvedData = null;
    phaseB1b4ResolvedData = null;
    finalExcelBase64 = null;
    finalFilename = null;
    finalCompounds = null;
    setPhaseBButtonsState({
        b1aDisabled: true,
        b1bDisabled: true,
        b1b1Disabled: true,
        b1b2Disabled: true,
        b1b3Disabled: true,
        b1b4Disabled: true,
        b2Disabled: true,
        b2Label: '<i data-lucide="download"></i> Create & Download Workbook',
        statusText: 'Phase A is ready. Enrich and classify the structures, then finalize the analysis before creating the workbook.',
        b1aState: 'pending',
        b1b1State: 'pending',
        b1b2State: 'pending',
        b1b3State: 'pending',
        b1b4State: 'pending',
        b2State: 'pending'
    });
    document.getElementById('pdf_file').value = '';
    try { sessionStorage.removeItem('pdfWorkflowStarted'); } catch (e) {}
    try { localStorage.removeItem('phaseBRunId'); } catch (e) {}
    updatePdfStepper();
}

// Global cancels
document.getElementById('pdfCancelBtn1')?.addEventListener('click', resetPdfView);
document.getElementById('pdfCancelBtn2')?.addEventListener('click', resetPdfView);

document.getElementById('pdfForm').addEventListener('submit', async (e) => {
    e.preventDefault();
    const btn = document.getElementById('pdfBtn');
    const fileInput = document.getElementById('pdf_file');
    if(!fileInput.files[0]) return;
    btn.disabled = true;
    btn.innerHTML = 'Ingesting PDF...';
    showGlobalLoader('Ingesting PDF...');
    const formData = new FormData();
    formData.append('file', fileInput.files[0]);
    try {
        const response = await fetch('/assess/pdf/step1', {method: 'POST', body: formData});
        if (!response.ok) throw new Error('Error ' + response.status);
        const data = await response.json();
        if(data.error) throw new Error(data.error);
        
        currentFilename = data.filename;
        try { sessionStorage.setItem('pdfWorkflowStarted', 'true'); } catch (e) {}
        
        // Populate UI
        document.getElementById('pdf-agency').innerText = data.agency || 'Unknown';
        document.getElementById('pdf-analysis-type').innerText = data.analysis_type || 'Unknown';
        
        // Show step 1
        document.getElementById('pdf-step0').style.display = 'none';
        document.getElementById('pdf-step1').style.display = 'block';
        updatePdfStepper();
        
    } catch (err) { alert(err.message); } 
    finally { hideGlobalLoader(); btn.disabled = false; btn.innerHTML = '<i data-lucide="cpu"></i> Auto-Extract'; lucide.createIcons(); }
});

document.getElementById('pdfContinue1').addEventListener('click', async () => {
    if (!currentFilename) return;
    const btn = document.getElementById('pdfContinue1');
    btn.disabled = true;
    btn.innerHTML = 'Extracting Tables (Please wait)...';
    showGlobalLoader('Extracting Tables...');
    
    const formData = new FormData();
    formData.append('filename', currentFilename);
    // Fast mode is an application default; the UI no longer exposes a toggle.
    const useParallel = true;
    formData.append('use_parallel', useParallel);
    
    try {
        const response = await fetch('/assess/pdf/step234', {method: 'POST', body: formData});
        if (!response.ok) throw new Error('Error ' + response.status);
        
        document.getElementById('pdf-stream-progress').style.display = 'block';
        const reader = response.body.getReader();
        const decoder = new TextDecoder();
        let buffer = '';
        let extractedData = null;
        
        while (true) {
            const { done, value } = await reader.read();
            if (done) break;
            
            buffer += decoder.decode(value, { stream: true });
            let lines = buffer.split('\n');
            buffer = lines.pop(); 
            
            for (let line of lines) {
                if (!line.trim()) continue;
                let msg;
                try {
                    msg = JSON.parse(line);
                } catch (e) {
                    if (e.message !== "Unexpected end of JSON input") {
                        console.error("Parse error on line:", line, e);
                    }
                    continue;
                }
                
                if (msg.type === 'progress') {
                    document.getElementById('pdf-stream-text').innerText = msg.message;
                    document.getElementById('global-loader-text').innerText = msg.message;
                } else if (msg.type === 'result') {
                    extractedData = msg.data;
                } else if (msg.type === 'error') {
                    throw new Error(msg.message);
                }
            }
        }
        
        document.getElementById('pdf-stream-progress').style.display = 'none';
        
        if (!extractedData) {
            throw new Error("Stream closed without sending result data");
        }
        
        const data = extractedData;
        
        // Populate UI
        document.getElementById('pdf-tables-summary').innerText = `Found ${data.total_tables} total tables. ${data.valid_tables} valid tables parsed, ${data.invalid_tables} invalid tables skipped.`;
        document.getElementById('pdf-raw-compounds').innerText = data.raw_compounds_count;
        
        // Render Captions
        window.currentCaptions = data.captions || [];
        renderCaptions();
        
        // Show step 234
        document.getElementById('pdf-step1').style.display = 'none';
        document.getElementById('pdf-step234').style.display = 'block';
        updatePdfStepper();
        
    } catch (err) { alert(err.message); } 
    finally { 
        hideGlobalLoader();
        btn.disabled = false; 
        btn.innerHTML = 'Continue to Table Extraction'; 
    }
});

document.getElementById('pdfContinue2').addEventListener('click', async () => {
    if (!currentFilename) return;
    const btn = document.getElementById('pdfContinue2');
    btn.disabled = true;
    btn.innerHTML = 'Deduplicating...';
    showGlobalLoader('Deduplicating Compounds...');
    
    const formData = new FormData();
    formData.append('filename', currentFilename);
    
    try {
        const response = await fetch('/assess/pdf/step5', {method: 'POST', body: formData});
        if (!response.ok) throw new Error('Error ' + response.status);
        const data = await response.json();
        if(data.error) throw new Error(data.error);
        
        currentExtractedData = data;
        currentExtractedData.filename = currentFilename; // ensure we pass filename to next step
        
        // Populate UI
        document.getElementById('pdf-unique-compounds').innerText = data.unique_compounds_count;
        
        // Show step 5
        document.getElementById('pdf-step234').style.display = 'none';
        document.getElementById('pdf-step5').style.display = 'block';
        updatePdfStepper();
        
    } catch (err) { alert(err.message); } 
    finally { 
        hideGlobalLoader();
        btn.disabled = false; 
        btn.innerHTML = 'Continue to Deduplication'; 
    }
});

let finalExcelBase64 = null;
let finalFilename = null;
let finalCompounds = null;
let phaseAExcelBase64 = null;
let phaseAFilename = null;
let phaseBRunId = null;

let currentResolvedData = null;
let phaseB1aResolvedData = null;
let phaseB1bResolvedData = null;
let phaseB1b1ResolvedData = null;
let phaseB1b2ResolvedData = null;
let phaseB1b3ResolvedData = null;
let phaseB1b4ResolvedData = null;

async function recoverPhaseBRun() {
    // Recover only runs started in this browser session. This prevents a
    // stale prior run from replacing the initial upload step on a new visit.
    let sessionStarted = false;
    try { sessionStarted = sessionStorage.getItem('pdfWorkflowStarted') === 'true'; } catch (e) {}
    if (!sessionStarted) return;
    let savedId = null;
    try { savedId = localStorage.getItem('phaseBRunId'); } catch (e) {}
    if (!savedId) return;
    try {
        const response = await fetch(`/assess/pdf/phase_b_state/${encodeURIComponent(savedId)}`);
        if (!response.ok) throw new Error('state unavailable');
        const state = await response.json();
        const completedStages = new Set(state.completed_stages || []);
        // Do not reopen the export panel for an identity-only (Phase A) run.
        // A fresh visit should begin at PDF Step 1; only an active Phase B
        // checkpoint is eligible for browser-refresh recovery.
        const phaseBStarted = ['B1A', 'B1B', 'B1B-finalize', 'B1B-chemont', 'B1B-epa_ctx', 'B1B-chembl', 'B2']
            .some(stage => completedStages.has(stage));
        if (!phaseBStarted) {
            try { sessionStorage.removeItem('pdfWorkflowStarted'); } catch (ignore) {}
            try { localStorage.removeItem('phaseBRunId'); } catch (ignore) {}
            return;
        }
        phaseBRunId = state.run_id;
        currentResolvedData = state;
        ['pdf-step0', 'pdf-step1', 'pdf-step234', 'pdf-step5', 'pdf-step6'].forEach(id => {
            const el = document.getElementById(id);
            if (el) el.style.display = 'none';
        });
        document.getElementById('pdf-step55').style.display = 'block';
        updatePdfStepper();
        const done = completedStages;
        const finalized = done.has('B1B-finalize') || done.has('B1B') || done.has('B1B-chembl');
        setPhaseBButtonsState({
            b1aDisabled: done.has('B1A'), b1b1Disabled: !done.has('B1A') || done.has('B1B-chemont'),
            b1b2Disabled: !done.has('B1B-chemont') || done.has('B1B-epa_ctx'),
            b1b3Disabled: !done.has('B1B-epa_ctx') || done.has('B1B-chembl'),
            b1b4Disabled: !done.has('B1B-chembl') || done.has('B1B-finalize'),
            b2Disabled: !done.has('B1A') && !finalized,
            statusText: `Recovered Phase B run (${state.compounds?.length || 0} rows). Continue with the next available checkpoint.`
        });
    } catch (e) {
        try { localStorage.removeItem('phaseBRunId'); } catch (ignore) {}
    }
}

function setPhaseBCheckpointState({ b1a = 'pending', b1b = 'pending', b1b1 = 'pending', b1b2 = 'pending', b1b3 = 'pending', b1b4 = 'pending', b2 = 'pending' } = {}) {
    const b1aEl = document.getElementById('phaseB1aCheckpoint');
    const b1bEl = document.getElementById('phaseB1bCheckpoint');
    const b1b1El = document.getElementById('phaseB1b1Checkpoint');
    const b1b2El = document.getElementById('phaseB1b2Checkpoint');
    const b1b3El = document.getElementById('phaseB1b3Checkpoint');
    const b1b4El = document.getElementById('phaseB1b4Checkpoint');
    const b2El = document.getElementById('phaseB2Checkpoint');
    const summaryEl = document.getElementById('phaseBCheckpointSummary');
    const progressEl = document.getElementById('phaseBProgressBar');
    const styles = {
        pending: { text: 'Pending', color: '#fbbf24' },
        running: { text: 'Running', color: '#38bdf8' },
        complete: { text: 'Complete', color: '#10b981' },
        skipped: { text: 'Skipped', color: '#94a3b8' }
    };
    const apply = (el, state) => {
        if (!el) return;
        const s = styles[state] || styles.pending;
        el.innerText = s.text;
        el.style.color = s.color;
    };
    apply(b1aEl, b1a);
    apply(b1bEl, b1b);
    apply(b1b1El, b1b1);
    apply(b1b2El, b1b2);
    apply(b1b3El, b1b3);
    apply(b1b4El, b1b4);
    apply(b2El, b2);

    if (summaryEl) {
        if (b2 === 'complete') {
            summaryEl.innerText = 'Analysis and workbook export are complete.';
        } else if (b2 === 'running') {
            summaryEl.innerText = 'Creating and downloading the final workbook.';
        } else if (b1b === 'complete') {
            summaryEl.innerText = 'Analysis is complete. The final workbook is ready to create and download.';
        } else if (b1a === 'complete') {
            summaryEl.innerText = 'Enrichment is complete. Finalize the analysis before creating the workbook.';
        } else if (b1b === 'running') {
            summaryEl.innerText = 'Finalizing analysis. Workbook creation remains locked until this step completes.';
        } else if (b1a === 'running') {
            summaryEl.innerText = 'Enriching and classifying structures. Final analysis remains locked until this step completes.';
        } else {
            summaryEl.innerText = 'Start by enriching and classifying the resolved structures.';
        }
    }
    if (progressEl) {
        const width = b2 === 'complete' ? 100 : b2 === 'running' ? 84 : b1b === 'complete' ? 67 : b1b === 'running' ? 50 : b1a === 'complete' ? 34 : b1a === 'running' ? 17 : 0;
        progressEl.style.width = `${width}%`;
    }
}

function setPhaseBButtonsState({ b1aDisabled, b1bDisabled, b1b1Disabled, b1b2Disabled, b1b3Disabled, b1b4Disabled, b2Disabled, b1aLabel, b1bLabel, b1b1Label, b1b2Label, b1b3Label, b1b4Label, b2Label, statusText, b1aState, b1bState, b1b1State, b1b2State, b1b3State, b1b4State, b2State } = {}) {
    const b1aBtn = document.getElementById('pdfRunB1aBtn');
    const b1bBtn = document.getElementById('pdfRunB1bBtn');
    const b1b1Btn = document.getElementById('pdfRunB1b1Btn');
    const b1b2Btn = document.getElementById('pdfRunB1b2Btn');
    const b1b3Btn = document.getElementById('pdfRunB1b3Btn');
    const b1b4Btn = document.getElementById('pdfRunB1b4Btn');
    const b2Btn = document.getElementById('pdfRunB2Btn');
    const b2LockedNote = document.getElementById('pdfB2LockedNote');
    const statusEl = document.getElementById('phaseBStatusText');
    const applyVisibleButton = (button, disabled, label, state, readyLabel, completeLabel) => {
        if (!button) return;
        if (typeof disabled === 'boolean') button.disabled = disabled;
        button.classList.toggle('is-complete', state === 'complete');
        if (state === 'complete') button.innerHTML = `<i data-lucide="check"></i> ${completeLabel}`;
        else if (state === 'pending') button.innerHTML = `<i data-lucide="play"></i> ${readyLabel}`;
        else if (label) button.innerHTML = label;
    };
    applyVisibleButton(b1aBtn, b1aDisabled, b1aLabel, b1aState, 'Enrich & Classify', 'Enrichment Complete');
    applyVisibleButton(b1bBtn, b1bDisabled, b1bLabel, b1bState, 'Finalize Analysis', 'Analysis Complete');
    [[b1b1Btn, b1b1Disabled, b1b1Label], [b1b2Btn, b1b2Disabled, b1b2Label], [b1b3Btn, b1b3Disabled, b1b3Label], [b1b4Btn, b1b4Disabled, b1b4Label]].forEach(([btn, disabled, label]) => {
        if (btn && typeof disabled === 'boolean') btn.disabled = disabled;
        if (btn && label) btn.innerHTML = label;
    });
    if (b2Btn) {
        if (typeof b2Disabled === 'boolean') b2Btn.disabled = b2Disabled;
        const exportAvailable = b1bState === 'complete' || b2State === 'running' || b2State === 'complete';
        b2Btn.style.display = exportAvailable ? 'inline-flex' : 'none';
        if (b2LockedNote) b2LockedNote.style.display = exportAvailable ? 'none' : 'flex';
        b2Btn.classList.toggle('is-complete', b2State === 'complete');
        if (b2State === 'complete') b2Btn.innerHTML = '<i data-lucide="check"></i> Workbook Downloaded';
        else if (b2State === 'pending') b2Btn.innerHTML = '<i data-lucide="download"></i> Create & Download Workbook';
        else if (b2Label) b2Btn.innerHTML = b2Label;
    }
    if (statusEl && statusText) statusEl.innerText = statusText;
    if (typeof b1aState !== 'undefined' || typeof b1bState !== 'undefined' || typeof b2State !== 'undefined') {
        setPhaseBCheckpointState({
            b1a: b1aState || 'pending',
            b1b: b1bState || 'pending',
            b1b1: b1b1State || 'pending',
            b1b2: b1b2State || 'pending',
            b1b3: b1b3State || 'pending',
            b1b4: b1b4State || 'pending',
            b2: b2State || 'pending'
        });
    }
    lucide.createIcons();
}

async function postDeterministicPhaseB(payload, stageLabel) {
    const postPhaseB = async (payload, stageLabel) => {
        const response = await fetch('/assess/pdf/run_deterministic_enrichment', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payload)
        });
        if (!response.ok) {
            let errText = await response.text();
            try {
                const errObj = JSON.parse(errText);
                errText = errObj.detail || errText;
                if (errObj.trace) errText += "\n\n" + errObj.trace;
                if (errObj.traceback) errText += "\n\n" + errObj.traceback;
                if (Array.isArray(errObj.phase_b_trace)) {
                    errText += "\n\nPhase trace:\n" + JSON.stringify(errObj.phase_b_trace, null, 2);
                }
            } catch (e) {}
            throw new Error(stageLabel + ' failed with HTTP ' + response.status + ': ' + errText);
        }
        const result = await response.json();
        if (result.run_id) {
            phaseBRunId = result.run_id;
            try { localStorage.setItem('phaseBRunId', result.run_id); } catch (e) {}
        }
        return result;
    };
    return await postPhaseB(payload, stageLabel);
}

async function downloadPhaseAExcel() {
    if (!currentResolvedData) return null;

    showGlobalLoader('Phase A: generating Excel snapshot...');
    const btn = document.getElementById('pdfDownloadPhaseATopBtn');
    if (btn) {
        btn.disabled = true;
        btn.innerHTML = '<i data-lucide="loader-circle"></i> Downloading Phase A...';
    }
    try {
        const response = await fetch('/assess/pdf/run_deterministic_enrichment', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                filename: currentResolvedData.filename || 'report.pdf',
                compounds: currentResolvedData.compounds || [],
                phase_b_stage: 'A',
                run_id: phaseBRunId || currentResolvedData.run_id || null
            })
        });
        if (!response.ok) {
            let errText = await response.text();
            try {
                const errObj = JSON.parse(errText);
                errText = errObj.detail || errText;
                if (errObj.trace) errText += "\n\n" + errObj.trace;
                if (errObj.traceback) errText += "\n\n" + errObj.traceback;
                if (Array.isArray(errObj.phase_b_trace)) {
                    errText += "\n\nPhase trace:\n" + JSON.stringify(errObj.phase_b_trace, null, 2);
                }
            } catch (e) {}
            throw new Error('Phase A snapshot failed with HTTP ' + response.status + ': ' + errText);
        }
        const data = await response.json();
        phaseAExcelBase64 = data.excel_base64;
        phaseAFilename = data.filename;
        const byteCharacters = atob(phaseAExcelBase64);
        const byteNumbers = new Array(byteCharacters.length);
        for (let i = 0; i < byteCharacters.length; i++) {
            byteNumbers[i] = byteCharacters.charCodeAt(i);
        }
        const byteArray = new Uint8Array(byteNumbers);
        const blob = new Blob([byteArray], {type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'});
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = phaseAFilename || 'phase_a_snapshot.xlsx';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        window.URL.revokeObjectURL(url);
        const statusEl = document.getElementById('phaseBStatusText');
        if (statusEl) statusEl.innerText = 'Phase A snapshot downloaded. Start enrichment and classification when ready.';
        return data;
    } finally {
        if (btn) {
            btn.disabled = false;
            btn.innerHTML = '<i data-lucide="file-down"></i> Download Phase A Excel';
        }
        hideGlobalLoader();
        lucide.createIcons();
    }
}

async function runPhaseB1a() {
    if (!currentResolvedData) return null;

    showGlobalLoader('Phase B1a: deterministic core enrichment...');
    setPhaseBButtonsState({
        b1aDisabled: true,
        b1bDisabled: false,
        b1b1Disabled: true,
        b1b2Disabled: true,
        b1b3Disabled: true,
        b1b4Disabled: true,
        b2Disabled: true,
        b1aLabel: '<i data-lucide="loader-circle"></i> Enriching & Classifying...',
        b1bLabel: '<i data-lucide="play"></i> Finalize Analysis',
        b2Label: '<i data-lucide="download"></i> Create & Download Workbook',
        statusText: 'Enriching and classifying the resolved structures.',
        b1aState: 'running',
        b1bState: 'pending',
        b2State: 'pending'
    });
    let b1Data;
    try {
        b1Data = await postDeterministicPhaseB({...currentResolvedData, phase_b_stage: 'B1A', run_id: phaseBRunId || currentResolvedData.run_id || null}, 'Phase B1a');
    } catch (error) {
        setPhaseBButtonsState({ b1aDisabled: false, b1bDisabled: true, b1b1Disabled: true, b1b2Disabled: true, b1b3Disabled: true, b1b4Disabled: true, b2Disabled: true, statusText: 'Phase B1a failed. Retry this checkpoint.', b1aState: 'pending', b1bState: 'pending', b2State: 'pending' });
        hideGlobalLoader();
        throw error;
    }
    phaseBRunId = b1Data.run_id || phaseBRunId;
    phaseB1aResolvedData = {
        ...currentResolvedData,
        compounds: b1Data.compounds || [],
        phase_b_trace: b1Data.phase_b_trace || []
        ,run_id: phaseBRunId
    };
    currentResolvedData = phaseB1aResolvedData;
    setPhaseBButtonsState({
        b1aDisabled: true,
        b1bDisabled: false,
        b1b1Disabled: true,
        b1b2Disabled: true,
        b1b3Disabled: true,
        b1b4Disabled: true,
        b2Disabled: true,
        b1aLabel: '<i data-lucide="play"></i> Run Phase B1a',
        b1bLabel: '<i data-lucide="play"></i> Run Phase B1b',
        b2Label: '<i data-lucide="download"></i> Export Final Workbook',
        statusText: 'Enrichment complete. Finalize the analysis before creating the workbook.',
        b1aState: 'complete',
        b1bState: 'pending',
        b1b1State: 'pending',
        b1b2State: 'pending',
        b1b3State: 'pending',
        b1b4State: 'pending',
        b2State: 'pending'
    });
    return b1Data;
}

async function runPhaseB1b() {
    const sourceData = phaseB1aResolvedData || currentResolvedData;
    if (!sourceData) return null;

    showGlobalLoader('Phase B1b: finalizing local evidence...');
    setPhaseBButtonsState({
        b1aDisabled: true,
        b1bDisabled: true,
        b2Disabled: true,
        b1aLabel: '<i data-lucide="play"></i> Run Phase B1a',
        b1bLabel: '<i data-lucide="loader-circle"></i> Finalizing Analysis...',
        b2Label: '<i data-lucide="download"></i> Create & Download Workbook',
        statusText: 'Finalizing local evidence before export.',
        b1aState: 'complete',
        b1bState: 'running',
        b2State: 'pending'
    });
    const b1bData = await postDeterministicPhaseB({...sourceData, phase_b_stage: 'B1B4', run_id: phaseBRunId || sourceData.run_id || null}, 'Phase B1b');
    phaseBRunId = b1bData.run_id || phaseBRunId;
    phaseB1bResolvedData = {
        ...sourceData,
        compounds: b1bData.compounds || [],
        phase_b_trace: b1bData.phase_b_trace || []
        ,run_id: phaseBRunId
    };
    currentResolvedData = phaseB1bResolvedData;
    setPhaseBButtonsState({
        b1aDisabled: true,
        b1bDisabled: true,
        b2Disabled: false,
        b1aLabel: '<i data-lucide="play"></i> Run Phase B1a',
        b1bLabel: '<i data-lucide="play"></i> Run Phase B1b',
        b2Label: '<i data-lucide="download"></i> Create & Download Workbook',
        statusText: 'Analysis complete. Create and download the final workbook.',
        b1aState: 'complete',
        b1bState: 'complete',
        b2State: 'pending'
    });
    return b1bData;
}

async function runPhaseB1bSubstage(stage, label, sourceData, nextStore) {
    if (!sourceData) return null;
    showGlobalLoader(`${label}: processing...`);
    const stageKey = stage.toLowerCase();
    const state = { b1aState: 'complete', b1b1State: 'pending', b1b2State: 'pending', b1b3State: 'pending', b1b4State: 'pending', b2State: 'pending' };
    state[`${stageKey}State`] = 'running';
    setPhaseBButtonsState({ b1aDisabled: true, b1b1Disabled: true, b1b2Disabled: true, b1b3Disabled: true, b1b4Disabled: true, b2Disabled: true, statusText: `${label} is running.`, ...state });
    let response;
    try {
        response = await postDeterministicPhaseB({ ...sourceData, phase_b_stage: stage, run_id: phaseBRunId || sourceData.run_id || null }, label);
        phaseBRunId = response.run_id || phaseBRunId;
        const nextData = { ...sourceData, compounds: response.compounds || [], phase_b_trace: response.phase_b_trace || [], run_id: phaseBRunId };
        nextStore(nextData);
        currentResolvedData = nextData;
    } catch (error) {
        setPhaseBButtonsState({
            b1aDisabled: true,
            b1b1Disabled: stageKey !== 'b1b1',
            b1b2Disabled: stageKey !== 'b1b2',
            b1b3Disabled: stageKey !== 'b1b3',
            b1b4Disabled: stageKey !== 'b1b4',
            b2Disabled: true,
            statusText: `${label} failed. Retry this checkpoint.`,
            b1aState: 'complete', b1b1State: 'pending', b1b2State: 'pending', b1b3State: 'pending', b1b4State: 'pending', b2State: 'pending'
        });
        hideGlobalLoader();
        throw error;
    }
    const completeState = { b1aState: 'complete', b1b1State: 'pending', b1b2State: 'pending', b1b3State: 'pending', b1b4State: 'pending', b2State: 'pending' };
    completeState[`${stageKey}State`] = 'complete';
    const order = { b1b1: 1, b1b2: 2, b1b3: 3, b1b4: 4 }[stageKey];
    ['b1b1', 'b1b2', 'b1b3', 'b1b4'].forEach((key, index) => { if (index < order) completeState[`${key}State`] = 'complete'; });
    setPhaseBButtonsState({
        b1aDisabled: true,
        b1b1Disabled: order >= 1,
        b1b2Disabled: order < 1,
        b1b3Disabled: order < 2,
        b1b4Disabled: order < 3,
        b2Disabled: order < 4,
        statusText: `${label} complete. Continue with the next checkpoint.`,
        ...completeState
    });
    return response;
}

async function runPhaseB2() {
    const sourceData = phaseB1b4ResolvedData || phaseB1b3ResolvedData || phaseB1b2ResolvedData || phaseB1b1ResolvedData || phaseB1bResolvedData || phaseB1aResolvedData || currentResolvedData;
    if (!sourceData) return null;
    showGlobalLoader('Phase B2: workbook export...');
    setPhaseBButtonsState({
        b1aDisabled: true,
        b1bDisabled: true,
        b2Disabled: true,
        b1aLabel: '<i data-lucide="play"></i> Run Phase B1a',
        b1bLabel: '<i data-lucide="play"></i> Run Phase B1b',
        b2Label: '<i data-lucide="loader-circle"></i> Creating Workbook...',
        statusText: 'Creating and downloading the final workbook.',
        b1aState: phaseB1aResolvedData ? 'complete' : 'pending',
        b1bState: phaseB1bResolvedData ? 'complete' : 'pending',
        b1b1State: phaseB1b1ResolvedData ? 'complete' : 'pending',
        b1b2State: phaseB1b2ResolvedData ? 'complete' : 'pending',
        b1b3State: phaseB1b3ResolvedData ? 'complete' : 'pending',
        b1b4State: phaseB1b4ResolvedData ? 'complete' : 'pending',
        b2State: 'running'
    });
    const b2Payload = {
        filename: sourceData.filename || currentResolvedData?.filename || 'report.pdf',
        compounds: sourceData.compounds || currentResolvedData?.compounds || [],
        phase_b_stage: 'B2',
        run_id: phaseBRunId || sourceData.run_id || currentResolvedData?.run_id || null
    };
    const data = await postDeterministicPhaseB(b2Payload, 'Phase B2');
    finalExcelBase64 = data.excel_base64;
    finalFilename = data.filename;
    finalCompounds = data.compounds;
    document.getElementById('downloadExcelBtn').click();
    document.getElementById('pdf-step55').style.display = 'none';
    document.getElementById('pdf-step6').style.display = 'block';
    updatePdfStepper();
    setPhaseBButtonsState({
        b1aDisabled: false,
        b1bDisabled: false,
        b2Disabled: true,
        b1aLabel: '<i data-lucide="play"></i> Run Phase B1a',
        b1bLabel: '<i data-lucide="play"></i> Run Phase B1b',
        b2Label: '<i data-lucide="download"></i> Export Final Workbook',
        statusText: 'Workbook export complete.',
        b1aState: 'complete',
        b1bState: 'complete',
        b1b1State: 'complete',
        b1b2State: 'complete',
        b1b3State: 'complete',
        b1b4State: 'complete',
        b2State: 'complete'
    });
    return data;
}

document.getElementById('pdfProceedBtn').addEventListener('click', async () => {
    if (!currentExtractedData) return;
    
    const btn = document.getElementById('pdfProceedBtn');
    btn.disabled = true;
    btn.innerHTML = 'Resolving Identities...';
    showGlobalLoader('Step 6: Resolving Chemical Structures...');
    
    // Fast mode is always enabled for identity resolution.
    currentExtractedData.use_parallel = true;
    
    try {
        const response = await fetch('/assess/pdf/resolve_identity', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(currentExtractedData)
        });
        if (!response.ok) {
            let errText = await response.text();
            try {
                const errObj = JSON.parse(errText);
                errText = errObj.detail || errText;
                if(errObj.trace) errText += "\n\n" + errObj.trace;
            } catch(e) {}
            throw new Error('Error ' + response.status + ': ' + errText);
        }
        
        document.getElementById('pdf-stream-progress').style.display = 'block';
        const reader = response.body.getReader();
        const decoder = new TextDecoder();
        let buffer = '';
        
        while (true) {
            const { done, value } = await reader.read();
            if (done) break;
            
            buffer += decoder.decode(value, { stream: true });
            let lines = buffer.split('\n');
            buffer = lines.pop(); 
            
            for (let line of lines) {
                if (!line.trim()) continue;
                let msg;
                try {
                    msg = JSON.parse(line);
                } catch (e) {
                    if (e.message !== "Unexpected end of JSON input") {
                        console.error("Parse error on line:", line, e);
                    }
                    continue;
                }
                
                if (msg.type === 'progress') {
                    document.getElementById('pdf-stream-text').innerText = msg.message;
                    document.getElementById('global-loader-text').innerText = msg.message;
                } else if (msg.type === 'result') {
                    currentResolvedData = msg;
                    phaseBRunId = msg.run_id || phaseBRunId;
                } else if (msg.type === 'error') {
                    throw new Error(msg.message + (msg.trace ? '\n' + msg.trace : ''));
                }
            }
        }
        
        if (!currentResolvedData) {
            throw new Error("Stream closed without sending result data");
        }
        
        // Hide stream progress UI once done
        document.getElementById('pdf-stream-progress').style.display = 'none';
        
        // Show the workbook-ready stage after Phase A. Phase B is now a
        // deterministic enrichment/export step, not AI.
        document.getElementById('pdf-step5').style.display = 'none';
        document.getElementById('pdf-step55').style.display = 'block';
        setPhaseBButtonsState({
            b1aDisabled: false,
            b1bDisabled: true,
            b1b1Disabled: true,
            b1b2Disabled: true,
            b1b3Disabled: true,
            b1b4Disabled: true,
            b2Disabled: true,
            b1aLabel: '<i data-lucide="play"></i> Run Phase B1a',
            b1bLabel: '<i data-lucide="play"></i> Run Phase B1b',
            b2Label: '<i data-lucide="download"></i> Export Final Workbook',
            statusText: 'Phase A is complete. Download the snapshot or run B1a, then export the final workbook with B2. AI reasoning remains off for this step.',
            b1aState: 'pending',
            b1bState: 'pending',
            b1b1State: 'pending',
            b1b2State: 'pending',
            b1b3State: 'pending',
            b1b4State: 'pending',
            b2State: 'pending'
        });
        
        let compoundsCount = currentResolvedData.compounds ? currentResolvedData.compounds.length : 0;
        document.getElementById('resolutionStatsText').innerText = `Successfully analyzed ${compoundsCount} extracted compounds across the database.`;
        const chemontPreview = (currentResolvedData.compounds || []).find(c => c.ChemOnt_Retrieval_Status);
        const chemontSummaryEl = document.getElementById('pdf-chemont-summary');
        if (chemontSummaryEl) {
            if (chemontPreview) {
                const chemontClass = chemontPreview.ChemOnt_Class || 'Not available';
                const chemontSubclass = chemontPreview.ChemOnt_Subclass || 'Not available';
                const chemontStatus = chemontPreview.ChemOnt_Retrieval_Status || 'Not available';
                chemontSummaryEl.innerText = `ChemOnt Class: ${chemontClass} | ChemOnt Subclass: ${chemontSubclass} | ChemOnt Retrieval Status: ${chemontStatus}`;
            } else {
                chemontSummaryEl.innerText = 'ChemOnt preview unavailable in the resolved rows.';
            }
        }

    } catch (err) { alert(err.message); } 
    finally { 
        hideGlobalLoader();
        btn.disabled = false; 
        btn.innerHTML = '<i data-lucide="download"></i> Run Workbook Export'; 
        lucide.createIcons(); 
    }
});

document.getElementById('pdfRunB1aBtn').addEventListener('click', async () => {
    if (!currentResolvedData) return;
    const btn = document.getElementById('pdfRunB1aBtn');
    btn.disabled = true;
    try {
        await runPhaseB1a();
    } catch (err) {
        alert(err.message);
    } finally {
        hideGlobalLoader();
        lucide.createIcons();
    }
});

document.getElementById('pdfRunB1bBtn').addEventListener('click', async () => {
    if (!phaseB1aResolvedData && !currentResolvedData) return;
    const btn = document.getElementById('pdfRunB1bBtn');
    btn.disabled = true;
    try { await runPhaseB1b(); }
    catch (err) { alert(err.message); }
    finally { hideGlobalLoader(); lucide.createIcons(); }
});

async function bindPhaseB1bSubstage(buttonId, stage, label, sourceGetter, storeSetter) {
    const btn = document.getElementById(buttonId);
    if (!btn) return;
    btn.addEventListener('click', async () => {
        const source = sourceGetter();
        if (!source) return;
        btn.disabled = true;
        try {
            await runPhaseB1bSubstage(stage, label, source, storeSetter);
        } catch (err) {
            alert(err.message);
        } finally {
            hideGlobalLoader();
            lucide.createIcons();
        }
    });
}

bindPhaseB1bSubstage('pdfRunB1b1Btn', 'B1B1', 'Phase B1b-1 (ChemOnt)', () => phaseB1aResolvedData || currentResolvedData, data => { phaseB1b1ResolvedData = data; });
bindPhaseB1bSubstage('pdfRunB1b2Btn', 'B1B2', 'Phase B1b-2 (EPA CTX)', () => phaseB1b1ResolvedData, data => { phaseB1b2ResolvedData = data; });
bindPhaseB1bSubstage('pdfRunB1b3Btn', 'B1B3', 'Phase B1b-3 (ChEMBL)', () => phaseB1b2ResolvedData, data => { phaseB1b3ResolvedData = data; });
bindPhaseB1bSubstage('pdfRunB1b4Btn', 'B1B4', 'Phase B1b-4 (Finalize evidence)', () => phaseB1b3ResolvedData, data => { phaseB1b4ResolvedData = data; phaseB1bResolvedData = data; });

async function handlePhaseADownload() {
    if (!currentResolvedData) return;
    const btn = document.getElementById('pdfDownloadPhaseATopBtn');
    btn.disabled = true;
    try {
        await downloadPhaseAExcel();
    } catch (err) {
        alert(err.message);
    } finally {
        btn.disabled = false;
        lucide.createIcons();
    }
}
document.getElementById('pdfDownloadPhaseATopBtn').addEventListener('click', handlePhaseADownload);

document.getElementById('pdfRunB2Btn').addEventListener('click', async () => {
    if (!phaseB1aResolvedData && !phaseB1bResolvedData) return;
    const btn = document.getElementById('pdfRunB2Btn');
    btn.disabled = true;
    try {
        await runPhaseB2();
    } catch (err) {
        setPhaseBButtonsState({
            b1aDisabled: true, b1bDisabled: true, b1b1Disabled: true,
            b1b2Disabled: true, b1b3Disabled: true, b1b4Disabled: true,
            b2Disabled: false,
            statusText: 'Phase B2 failed. Retry workbook export.',
            b1aState: 'complete', b1bState: 'complete', b2State: 'pending'
        });
        alert(err.message);
    } finally {
        hideGlobalLoader();
        lucide.createIcons();
    }
});

// PDF button logic removed as per user request

document.getElementById('downloadExcelBtn').addEventListener('click', () => {
    if (!finalExcelBase64) return;
    const byteCharacters = atob(finalExcelBase64);
    const byteNumbers = new Array(byteCharacters.length);
    for (let i = 0; i < byteCharacters.length; i++) {
        byteNumbers[i] = byteCharacters.charCodeAt(i);
    }
    const byteArray = new Uint8Array(byteNumbers);
    const blob = new Blob([byteArray], {type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'});
    
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    let outName = finalFilename;
    if (outName.toLowerCase().endsWith('.pdf')) {
        outName = outName.substring(0, outName.length - 4) + '.xlsx';
    } else if (!outName.toLowerCase().endsWith('.xlsx')) {
        outName = outName + '.xlsx';
    }
    a.download = outName;
    document.body.appendChild(a);
    a.click();
    a.remove();
    window.URL.revokeObjectURL(url);
});

// Table rendering fully removed

function renderCaptions() {
    const container = document.getElementById("pdf-captions-list");
    if(!container) return;
    const showValidOnly = document.getElementById("validTablesOnlyToggle").checked;
    
    if(!window.currentCaptions || window.currentCaptions.length === 0) {
        container.innerHTML = "<div>No tables found.</div>";
        return;
    }
    
    let html = "";
    window.currentCaptions.forEach(c => {
        if(showValidOnly && !c.is_valid) return;
        
        let icon = c.is_valid ? `<i data-lucide="check" style="color: #10b981; width:14px;"></i>` : `<i data-lucide="x" style="color: #ef4444; width:14px;"></i>`;
        html += `<div style="padding: 6px 0; border-bottom: 1px solid #1e293b; display: flex; align-items:flex-start; gap:8px;">
            <div style="flex-shrink:0; margin-top:2px;">${icon}</div>
            <div>
                <span style="color: ${c.is_valid ? "#cbd5e1" : "#64748b"};">Page ${c.page}: ${c.caption}</span>
            </div>
        </div>`;
    });
    container.innerHTML = html;
    lucide.createIcons();
}

document.getElementById("validTablesOnlyToggle").addEventListener("change", renderCaptions);

// Start every new page visit at the PDF upload step rather than silently
// restoring a stale run from browser storage.
resetPdfView();
