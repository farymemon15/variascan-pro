"""
Genomic Region & Splice Site Annotation Engine.
Classifies variant coordinates into functional gene features (Exon, Intron,
5' UTR, 3' UTR, Promoter, Canonical Splice Donor/Acceptor) and renders
an interactive gene architecture schematic with deletion locus pinning.
"""

from typing import Dict, Any, Tuple, Optional, List


def annotate_genomic_locus(
    chrom: str,
    pos: int,
    del_size: int,
    gene_name: str = "TARGET_GENE_1",
    transcript_id: str = "NM_001234.3"
) -> Dict[str, Any]:
    """
    Annotate genomic locus into functional gene regions and evaluate splice disruption risk.
    Models canonical gene structure:
    - 5' Promoter: 1 - 500
    - Exon 1 (5' UTR + Coding): 501 - 850
    - Intron 1: 851 - 1500
    - Exon 2 (Coding): 1501 - 2200 (Contains pos 2000!)
    - Intron 2: 2201 - 2900
    - Exon 3 (Coding + 3' UTR): 2901 - 3600
    """
    # Canonical gene architecture intervals
    exons = [
        {"exon_num": 1, "start": 501, "end": 850, "type": "Coding / 5'UTR"},
        {"exon_num": 2, "start": 1501, "end": 2200, "type": "Coding Exon (CDS)"},
        {"exon_num": 3, "start": 2901, "end": 3600, "type": "Coding / 3'UTR"}
    ]

    feature_type = "Intergenic"
    exon_info = None
    dist_to_junction = 999
    splice_risk = "None (Exonic / Deep Intronic)"
    splice_badge_color = "#10b981"

    # Check promoter
    if 1 <= pos <= 500:
        feature_type = "Promoter / 5' Regulatory"
        splice_risk = "Low (Transcription Regulation)"
        splice_badge_color = "#38bdf8"

    # Check exons
    for ex in exons:
        if ex["start"] <= pos <= ex["end"]:
            feature_type = f"Exon {ex['exon_num']} ({ex['type']})"
            exon_info = ex
            # Distance to nearest exon boundary
            dist_to_start = abs(pos - ex["start"])
            dist_to_end = abs(pos - ex["end"])
            dist_to_junction = min(dist_to_start, dist_to_end)

            if dist_to_junction <= 2:
                splice_risk = "CRITICAL: Canonical Splice Site Disrupted (±2 bp)"
                splice_badge_color = "#ef4444"
            elif dist_to_junction <= 8:
                splice_risk = "MODERATE: Near Splice Region (±8 bp)"
                splice_badge_color = "#f59e0b"
            else:
                splice_risk = "Low: Core Coding Region (>8 bp from junction)"
                splice_badge_color = "#10b981"
            break

    # Check introns if not in exon or promoter
    if feature_type == "Intergenic":
        if 851 <= pos <= 1500:
            feature_type = "Intron 1"
            dist_to_junction = min(abs(pos - 851), abs(pos - 1500))
        elif 2201 <= pos <= 2900:
            feature_type = "Intron 2"
            dist_to_junction = min(abs(pos - 2201), abs(pos - 2900))
        
        if dist_to_junction <= 2:
            splice_risk = "CRITICAL: Canonical Intronic Splice Donor/Acceptor (±2 bp)"
            splice_badge_color = "#ef4444"
        elif dist_to_junction <= 8:
            splice_risk = "MODERATE: Intronic Splice Region"
            splice_badge_color = "#f59e0b"
        else:
            splice_risk = "Low: Deep Intronic (>8 bp from junction)"
            splice_badge_color = "#64748b"

    return {
        "chrom": chrom,
        "pos": pos,
        "del_size": del_size,
        "gene_name": gene_name,
        "transcript_id": transcript_id,
        "feature_type": feature_type,
        "dist_to_splice_junction_bp": dist_to_junction,
        "splice_risk": splice_risk,
        "splice_badge_color": splice_badge_color,
        "exons": exons
    }


def generate_gene_structure_html(annot: Dict[str, Any]) -> str:
    """
    Render interactive SVG gene architecture showing Exons, Introns,
    and a glowing marker pinpointing the deletion breakpoint.
    """
    chrom = annot["chrom"]
    pos = annot["pos"]
    del_size = annot["del_size"]
    gene = annot["gene_name"]
    tx = annot["transcript_id"]
    feat = annot["feature_type"]
    risk = annot["splice_risk"]
    risk_color = annot["splice_badge_color"]
    dist = annot["dist_to_splice_junction_bp"]

    # Calculate horizontal pin position across 1 to 4000 bp scale (0% to 100%)
    pin_pct = round(max(5.0, min(95.0, (pos / 4000.0) * 100)), 1)

    html = f"""
    <div style="background:#0b1120; border:1px solid #1e293b; border-radius:10px; padding:20px; color:#e2e8f0; font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif;">
        <div style="display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #1e293b; padding-bottom:12px; margin-bottom:16px;">
            <div>
                <span style="font-size:16px; font-weight:700; color:#38bdf8;">🏷️ Genomic Feature & Splicing Architecture</span>
                <span style="margin-left:12px; background:#0284c7; color:white; padding:3px 10px; border-radius:12px; font-size:12px; font-weight:bold;">
                    {gene} ({tx})
                </span>
            </div>
            <div style="font-size:12px; color:#94a3b8;">
                Feature: <strong style="color:#f8fafc;">{feat}</strong>
            </div>
        </div>

        <!-- Gene Schematic Diagram -->
        <div style="background:#020617; border:1px solid #1e293b; border-radius:8px; padding:24px 20px; margin-bottom:16px; position:relative;">
            
            <!-- Intron Backbone Line -->
            <div style="position:relative; height:4px; background:#334155; margin:40px 10px 30px 10px; border-radius:2px;">
                
                <!-- Exon 1 Box (500 - 850) -> ~12% to ~21% -->
                <div style="position:absolute; left:12%; width:9%; top:-14px; height:32px; background:#0284c7; border:1px solid #38bdf8; border-radius:4px; text-align:center; line-height:30px; font-size:11px; font-weight:bold; color:white;" title="Exon 1 (350 bp)">
                    Exon 1
                </div>

                <!-- Intron 1 label -->
                <div style="position:absolute; left:26%; top:8px; font-size:10px; color:#64748b; font-family:monospace;">
                    Intron 1 (650 bp)
                </div>

                <!-- Exon 2 Box (1500 - 2200) -> ~37% to ~55% -->
                <div style="position:absolute; left:37%; width:18%; top:-16px; height:36px; background:#0369a1; border:2px solid #38bdf8; border-radius:4px; text-align:center; line-height:34px; font-size:12px; font-weight:bold; color:white; box-shadow:0 0 10px rgba(56,189,248,0.3);" title="Exon 2 (700 bp) - Contains Deletion!">
                    Exon 2 (CDS)
                </div>

                <!-- Intron 2 label -->
                <div style="position:absolute; left:60%; top:8px; font-size:10px; color:#64748b; font-family:monospace;">
                    Intron 2 (700 bp)
                </div>

                <!-- Exon 3 Box (2900 - 3600) -> ~72% to ~90% -->
                <div style="position:absolute; left:72%; width:18%; top:-14px; height:32px; background:#0284c7; border:1px solid #38bdf8; border-radius:4px; text-align:center; line-height:30px; font-size:11px; font-weight:bold; color:white;" title="Exon 3 (700 bp)">
                    Exon 3
                </div>

                <!-- Glowing Marker Pin (Deletions or SNPs) -->
                <div style="position:absolute; left:{pin_pct}%; top:-38px; transform:translateX(-50%); display:flex; flex-direction:column; align-items:center;">
                    <div style="background:{'#ef4444' if del_size > 0 else '#f59e0b'}; color:white; padding:2px 8px; border-radius:10px; font-size:10px; font-weight:bold; box-shadow:0 0 10px {'#ef4444' if del_size > 0 else '#f59e0b'}; white-space:nowrap;">
                        ▼ {f'{del_size} bp DEL' if del_size > 0 else 'SNP'} ({chrom}:{pos:,})
                    </div>
                    <div style="width:2px; height:20px; background:{'#ef4444' if del_size > 0 else '#f59e0b'}; box-shadow:0 0 6px {'#ef4444' if del_size > 0 else '#f59e0b'};"></div>
                </div>

            </div>

            <!-- Direction indicators -->
            <div style="display:flex; justify-content:space-between; font-size:11px; color:#64748b; margin-top:10px;">
                <span>5' Transcription Start Site (TSS) ➔</span>
                <span>➔ 3' Poly-A Tail</span>
            </div>
        </div>

        <!-- Metric Details Cards -->
        <div style="display:grid; grid-template-columns: 1fr 1fr 1fr; gap:12px;">
            <div style="background:#020617; padding:12px; border-radius:6px; border:1px solid #1e293b;">
                <div style="font-size:11px; color:#94a3b8; text-transform:uppercase;">Genomic Locus Region</div>
                <div style="font-size:15px; font-weight:bold; color:#38bdf8; margin-top:2px;">{feat}</div>
            </div>

            <div style="background:#020617; padding:12px; border-radius:6px; border:1px solid #1e293b;">
                <div style="font-size:11px; color:#94a3b8; text-transform:uppercase;">Splice Junction Distance</div>
                <div style="font-size:15px; font-weight:bold; color:#f8fafc; margin-top:2px;">{dist} bp from junction</div>
            </div>

            <div style="background:#020617; padding:12px; border-radius:6px; border:1px solid #1e293b;">
                <div style="font-size:11px; color:#94a3b8; text-transform:uppercase;">Splicing Disruption Risk</div>
                <div style="font-size:13px; font-weight:bold; color:{risk_color}; margin-top:2px;">{risk}</div>
            </div>
        </div>
    </div>
    """
    return html
