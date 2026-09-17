import os
import pandas as pd
from io import BytesIO
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from datetime import datetime

def generate_pdf_report(df: pd.DataFrame, filename: str) -> BytesIO:
    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter, rightMargin=30, leftMargin=30, topMargin=30, bottomMargin=30)
    
    styles = getSampleStyleSheet()
    title_style = styles['Heading1']
    subtitle_style = styles['Heading2']
    normal_style = styles['Normal']
    
    # Custom styles
    rationale_style = ParagraphStyle(
        'Rationale',
        parent=styles['Normal'],
        fontSize=8,
        leading=10,
        textColor=colors.darkslategray
    )
    
    elements = []
    
    # Title Page / Header
    elements.append(Paragraph("Enterprise Chemical Engine - Safety Dossier", title_style))
    elements.append(Spacer(1, 12))
    
    gen_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    elements.append(Paragraph(f"<b>Generated On:</b> {gen_time}", normal_style))
    elements.append(Paragraph(f"<b>Source File:</b> {filename}", normal_style))
    elements.append(Paragraph(f"<b>Total Compounds Assessed:</b> {len(df)}", normal_style))
    elements.append(Spacer(1, 24))
    
    # Overview
    elements.append(Paragraph("Executive Summary", subtitle_style))
    pass_count = len(df[df['Decision'].str.contains('PASS', na=False)])
    warn_count = len(df[df['Decision'].str.contains('WARNING', na=False)])
    fail_count = len(df[df['Decision'].str.contains('FAIL', na=False)])
    
    summary_text = f"Out of {len(df)} compounds, {pass_count} passed applicability domain checks, {warn_count} received structural warnings, and {fail_count} failed assessment or were out of domain."
    elements.append(Paragraph(summary_text, normal_style))
    elements.append(Spacer(1, 24))
    
    # Group by Cluster
    if 'Cluster ID' in df.columns:
        clusters = df.groupby('Cluster ID')
        for c_id, group in clusters:
            elements.append(Paragraph(f"Cluster: {c_id}", subtitle_style))
            
            # Table Data
            data = [['Compound', 'TTC Limit', 'Decision', 'Rationale']]
            
            for _, row in group.iterrows():
                name = str(row.get('Resolved Name', row.get('Cleaned Compound Name', 'Unknown')))
                ttc = str(row.get('TTC_Limit', 'N/A'))
                decision = str(row.get('Decision', 'N/A'))
                rationale_text = str(row.get('Read_Across_Justification', 'None provided'))
                
                # Wrap text for table
                p_name = Paragraph(name, normal_style)
                p_rationale = Paragraph(rationale_text, rationale_style)
                
                data.append([p_name, ttc, decision, p_rationale])
                
            t = Table(data, colWidths=[120, 60, 60, 300])
            t.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1e293b')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('VALIGN', (0, 0), (-1, -1), 'TOP'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), 10),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                ('BACKGROUND', (0, 1), (-1, -1), colors.HexColor('#f8fafc')),
                ('GRID', (0, 0), (-1, -1), 1, colors.HexColor('#cbd5e1')),
            ]))
            
            elements.append(t)
            elements.append(Spacer(1, 24))
            
    doc.build(elements)
    buffer.seek(0)
    return buffer
