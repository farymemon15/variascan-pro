"""
Mini-IGV Read Coverage Depth Track & Dip Plot Generator
Models and visualizes per-base sequencing read depth across the deletion locus,
highlighting the physical coverage drop (dip) caused by genomic deletions.
"""

import random
from pathlib import Path
from typing import Dict, Any, Tuple, Optional, List
import pandas as pd
import altair as alt


def generate_coverage_profile(
    chrom: str,
    pos: int,
    del_size: int,
    mean_depth: float = 38.0,
    vaf: float = 0.50,
    window_bp: int = 100,
    seed: int = 42
) -> pd.DataFrame:
    """
    Generate per-base sequencing depth profile around the deletion breakpoint.
    Models normal sequencing fluctuations and the distinct coverage drop
    within the deleted interval [pos, pos + del_size - 1].
    """
    rng = random.Random(seed + pos)
    start_pos = max(1, pos - window_bp)
    end_pos = pos + del_size + window_bp

    positions = []
    depths = []
    regions = []

    if del_size > 0:
        del_end = pos + del_size - 1
        dip_factor = max(0.0, 1.0 - vaf)  # e.g., 50% drop for HET, 100% drop for HOM
    else:
        del_end = pos
        dip_factor = 1.0  # SNPs do not cause coverage drop

    for p in range(start_pos, end_pos + 1):
        positions.append(p)
        if del_size > 0 and pos <= p <= del_end:
            # Inside deletion window: depth drops proportionally to VAF
            base_d = mean_depth * dip_factor
            noise = rng.gauss(0, max(1.0, mean_depth * 0.08))
            actual_d = max(0, int(round(base_d + noise)))
            depths.append(actual_d)
            regions.append(f"Deleted Window ({del_size} bp)")
        elif del_size == 0 and p == pos:
            # SNP position: full coverage with point mutation
            noise = rng.gauss(0, max(1.5, mean_depth * 0.12))
            actual_d = max(1, int(round(mean_depth + noise)))
            depths.append(actual_d)
            regions.append("SNP Locus")
        else:
            # Flanking sequence: full coverage with normal Poisson-like variance
            noise = rng.gauss(0, max(1.5, mean_depth * 0.12))
            actual_d = max(1, int(round(mean_depth + noise)))
            depths.append(actual_d)
            regions.append("Flanking Sequence")

    df = pd.DataFrame({
        "Genomic Position": positions,
        "Read Depth (x)": depths,
        "Region": regions,
        "Locus": [f"{chrom}:{p}" for p in positions]
    })
    return df


def build_coverage_chart(
    coverage_df: pd.DataFrame,
    chrom: str,
    pos: int,
    del_size: int,
    mean_depth: float = 38.0
) -> alt.Chart:
    """
    Build interactive Altair Mini-IGV coverage depth chart with shaded deletion/SNP window
    and baseline reference line.
    """
    del_start = pos
    del_end = pos + del_size - 1 if del_size > 0 else pos

    # Base Area Chart for Coverage Track
    base = alt.Chart(coverage_df).encode(
        x=alt.X(
            "Genomic Position:Q",
            scale=alt.Scale(zero=False),
            axis=alt.Axis(
                title=f"Genomic Coordinate ({chrom})",
                format="d",
                gridColor="#1e293b",
                labelColor="#94a3b8",
                titleColor="#e2e8f0"
            )
        ),
        y=alt.Y(
            "Read Depth (x):Q",
            scale=alt.Scale(zero=True),
            axis=alt.Axis(
                title="Sequencing Depth (DP)",
                gridColor="#1e293b",
                labelColor="#94a3b8",
                titleColor="#e2e8f0"
            )
        )
    )

    area = base.mark_area(
        color=alt.Gradient(
            gradient='linear',
            stops=[
                alt.GradientStop(color='#0284c7', offset=0),
                alt.GradientStop(color='#0f172a', offset=1)
            ],
            x1=1, x2=1, y1=1, y2=0
        ),
        opacity=0.7
    )

    line = base.mark_line(
        color='#38bdf8',
        strokeWidth=2
    ).encode(
        tooltip=[
            alt.Tooltip("Locus:N", title="Genomic Locus"),
            alt.Tooltip("Read Depth (x):Q", title="Depth"),
            alt.Tooltip("Region:N", title="Feature")
        ]
    )

    # Shaded mutation boundary (Red for deletions, Gold/Amber for SNPs)
    marker_color = '#ef4444' if del_size > 0 else '#f59e0b'
    del_box = alt.Chart(pd.DataFrame({
        "x1": [del_start],
        "x2": [del_end]
    })).mark_rect(
        color=marker_color,
        opacity=0.35 if del_size == 0 else 0.25
    ).encode(
        x='x1:Q',
        x2='x2:Q'
    )

    # Mean coverage reference line
    mean_line = alt.Chart(pd.DataFrame({"y": [mean_depth]})).mark_rule(
        color='#10b981',
        strokeDash=[5, 5],
        strokeWidth=1.5
    ).encode(
        y='y:Q'
    )

    title_text = f"Mini-IGV Coverage Track: {chrom}:{pos:,} (Drop at {del_size} bp Deletion Locus)" if del_size > 0 else f"Mini-IGV Coverage Track: {chrom}:{pos:,} (Point Mutation / SNP Locus)"
    sub_text = f"Red Shading = {del_size} bp Deletion Window [{del_start:,} - {del_end:,}] | Dotted Green = Mean Depth ({mean_depth:.1f}x)" if del_size > 0 else f"Amber Shading = Point Mutation Site [{pos:,}] | Dotted Green = Mean Depth ({mean_depth:.1f}x)"

    chart = (del_box + area + line + mean_line).properties(
        title=alt.TitleParams(
            text=title_text,
            subtitle=sub_text,
            color="#f8fafc",
            subtitleColor="#94a3b8",
            fontSize=14,
            subtitleFontSize=11
        ),
        height=260
    ).configure_view(
        strokeWidth=0,
        fill="#0b0f19"
    ).interactive()

    return chart
