"""
Automated Bioinformatics PDF Report Generator using ReportLab.
Compiles run metadata, QC statistics, candidate 25 bp deletions,
and locus alignment diagrams into a publication-quality PDF report.
"""

from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any, Optional
import pandas as pd

from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Preformatted, KeepTogether, HRFlowable
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle


def generate_pdf_report(
    output_pdf_path: str,
    target_bp: int,
    summary_stats: Dict[str, Any],
    deletions_df: pd.DataFrame,
    ascii_alignment: Optional[str] = None,
    is_admin: bool = True
) -> str:
    """Generate a publication-grade PDF summary report for target indel detection."""
    doc = SimpleDocTemplate(
        output_pdf_path,
        pagesize=letter,
        rightMargin=40,
        leftMargin=40,
        topMargin=40,
        bottomMargin=40
    )

    styles = getSampleStyleSheet()

    # Custom styles
    title_style = ParagraphStyle(
        'DocTitle',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=20,
        leading=24,
        textColor=colors.HexColor('#0f172a'),
        spaceAfter=6
    )
    subtitle_style = ParagraphStyle(
        'DocSubtitle',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=11,
        leading=15,
        textColor=colors.HexColor('#475569'),
        spaceAfter=15
    )
    h2_style = ParagraphStyle(
        'SectionH2',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=13,
        leading=17,
        textColor=colors.HexColor('#1e293b'),
        spaceBefore=12,
        spaceAfter=6
    )
    body_style = ParagraphStyle(
        'DocBody',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=9,
        leading=13,
        textColor=colors.HexColor('#334155')
    )
    code_style = ParagraphStyle(
        'AsciiCode',
        parent=styles['Normal'],
        fontName='Courier',
        fontSize=7.5,
        leading=10,
        textColor=colors.HexColor('#0f172a')
    )

    story = []

    # 1. Header Banner
    story.append(Paragraph("GENOMIC INDEL ANALYSIS REPORT", title_style))
    story.append(Paragraph(
        f"Automated NGS Detection & Breakpoint Profiling — Target: <b>{target_bp} bp Deletions</b>",
        subtitle_style
    ))
    story.append(HRFlowable(width="100%", thickness=1.5, color=colors.HexColor('#0284c7'), spaceAfter=15))

    # 2. Run Metadata & Configuration Table
    story.append(Paragraph("Pipeline Configuration & Execution Summary", h2_style))
    
    meta_data = [
        [
            Paragraph("<b>Timestamp:</b>", body_style),
            Paragraph(datetime.now().strftime("%Y-%m-%d %H:%M:%S UTC"), body_style),
            Paragraph("<b>Target Indel:</b>", body_style),
            Paragraph(f"{target_bp} bp Deletion", body_style),
        ],
        [
            Paragraph("<b>Sample Reference:</b>", body_style),
            Paragraph(str(summary_stats.get("reference", "N/A")), body_style),
            Paragraph("<b>Variant Caller:</b>", body_style),
            Paragraph(str(summary_stats.get("caller", "bcftools / GATK4")), body_style),
        ],
        [
            Paragraph("<b>Forward Reads (R1):</b>", body_style),
            Paragraph(str(summary_stats.get("r1_sample", "N/A")), body_style),
            Paragraph("<b>QC Filter:</b>", body_style),
            Paragraph("fastp (Q >= 30, Auto-PE)", body_style),
        ],
        [
            Paragraph("<b>Reverse Reads (R2):</b>", body_style),
            Paragraph(str(summary_stats.get("r2_sample", "N/A")), body_style),
            Paragraph("<b>Alignment:</b>", body_style),
            Paragraph("bwa mem | samtools sort", body_style),
        ],
    ]

    meta_table = Table(meta_data, colWidths=[1.4 * inch, 2.3 * inch, 1.4 * inch, 2.2 * inch])
    meta_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, -1), colors.HexColor('#f8fafc')),
        ('BOX', (0, 0), (-1, -1), 0.5, colors.HexColor('#cbd5e1')),
        ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#e2e8f0')),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('TOPPADDING', (0, 0), (-1, -1), 4),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
    ]))
    story.append(meta_table)
    story.append(Spacer(1, 12))

    # 3. Key Metrics Cards Table
    story.append(Paragraph("Key Biological & QC Metrics", h2_style))
    
    total_deletions = len(deletions_df)
    mean_depth = f"{deletions_df['Read Depth (DP)'].mean():.1f}x" if not deletions_df.empty and 'Read Depth (DP)' in deletions_df.columns else "N/A"
    mean_qual = f"{deletions_df['Quality Score (QUAL)'].mean():.1f}" if not deletions_df.empty and 'Quality Score (QUAL)' in deletions_df.columns else "N/A"
    total_raw_variants = summary_stats.get("total_raw_variants", total_deletions)

    kpi_data = [
        [
            Paragraph("<b>Total Target Deletions</b>", body_style),
            Paragraph("<b>Raw Variants Detected</b>", body_style),
            Paragraph("<b>Mean Read Depth</b>", body_style),
            Paragraph("<b>Mean Quality Score</b>", body_style),
        ],
        [
            Paragraph(f"<font size=14 color='#0284c7'><b>{total_deletions}</b></font>", body_style),
            Paragraph(f"<font size=14 color='#0f172a'><b>{total_raw_variants}</b></font>", body_style),
            Paragraph(f"<font size=14 color='#16a34a'><b>{mean_depth}</b></font>", body_style),
            Paragraph(f"<font size=14 color='#9333ea'><b>{mean_qual}</b></font>", body_style),
        ]
    ]

    kpi_table = Table(kpi_data, colWidths=[1.825 * inch] * 4)
    kpi_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#f1f5f9')),
        ('BACKGROUND', (0, 1), (-1, 1), colors.HexColor('#ffffff')),
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
        ('BOX', (0, 0), (-1, -1), 1, colors.HexColor('#cbd5e1')),
        ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#e2e8f0')),
        ('TOPPADDING', (0, 0), (-1, -1), 6),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
    ]))
    story.append(kpi_table)
    story.append(Spacer(1, 14))

    # 4. Detected Target Deletions Table
    story.append(Paragraph(f"Detected {target_bp} bp Deletion Candidates", h2_style))

    if deletions_df.empty:
        story.append(Paragraph("<i>No matching deletions found with the specified criteria.</i>", body_style))
    else:
        table_cols = ['Chromosome', 'Position', 'Reference Allele', 'Alternate Allele', 'Deletion Size (bp)', 'Genotype (HET/HOM)', 'Quality Score (QUAL)', 'Read Depth (DP)']
        
        # Prepare header
        header_row = [Paragraph(f"<b>{c}</b>", body_style) for c in table_cols]
        rows = [header_row]

        # Prepare rows (truncate alleles if too long for table width)
        for _, row in deletions_df.iterrows():
            ref_disp = str(row.get('Reference Allele', ''))
            alt_disp = str(row.get('Alternate Allele', ''))
            if len(ref_disp) > 12:
                ref_disp = ref_disp[:6] + ".." + ref_disp[-4:]
            if len(alt_disp) > 8:
                alt_disp = alt_disp[:4] + ".."

            rows.append([
                Paragraph(str(row.get('Chromosome', '')), body_style),
                Paragraph(f"{int(row.get('Position', 0)):,}", body_style),
                Paragraph(f"<code>{ref_disp}</code>", body_style),
                Paragraph(f"<code>{alt_disp}</code>", body_style),
                Paragraph(f"<b>{row.get('Deletion Size (bp)', '')} bp</b>", body_style),
                Paragraph(str(row.get('Genotype (HET/HOM)', '')), body_style),
                Paragraph(str(row.get('Quality Score (QUAL)', '')), body_style),
                Paragraph(f"{row.get('Read Depth (DP)', '')}x", body_style),
            ])

        del_table = Table(rows, colWidths=[0.85*inch, 0.85*inch, 1.25*inch, 0.85*inch, 1.0*inch, 1.1*inch, 0.7*inch, 0.7*inch])
        del_table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#0f172a')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.HexColor('#ffffff'), colors.HexColor('#f8fafc')]),
            ('BOX', (0, 0), (-1, -1), 0.5, colors.HexColor('#94a3b8')),
            ('INNERGRID', (0, 0), (-1, -1), 0.5, colors.HexColor('#e2e8f0')),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
        ]))
        story.append(del_table)

    story.append(Spacer(1, 14))

    # 5. Genomic Alignment Snapshot (ASCII Diagram)
    if ascii_alignment:
        story.append(KeepTogether([
            Paragraph("Genomic Breakpoint Sequence Alignment", h2_style),
            Spacer(1, 4),
            Preformatted(ascii_alignment, code_style)
        ]))

    story.append(Spacer(1, 16))
    if is_admin:
        story.append(Paragraph(
            "<font color='#059669'><b>[CONFIDENTIAL] Certified Molecular Diagnostic Report</b></font> &mdash; "
            "Analyzed by Lead Bioinformatician. Verified for precision clinical reporting.",
            subtitle_style
        ))
    else:
        story.append(Paragraph(
            "<font color='#0284c7'><b>[DEMO PREVIEW] Public Benchmark Run</b></font> &mdash; "
            "To commission full-scale certified clinical cohort analysis, contact: bioinformatics.services@contact.com",
            subtitle_style
        ))

    # Build PDF document
    doc.build(story)
    return output_pdf_path
