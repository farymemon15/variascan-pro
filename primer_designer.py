"""
Wet-Lab Validation Suite: PCR Primer Designer & Virtual Gel Electrophoresis Simulator.
Designs flanking primers (Forward and Reverse) to validate genomic deletions,
calculates melting temperatures (Tm), and renders an interactive Agarose Gel Electrophoresis
visualization comparing Wild-Type vs Deletion amplicons against a 100 bp DNA ladder.
"""

import math
from typing import Dict, Any, Tuple, Optional, List
from utils.demo_data import reverse_complement


def calculate_tm(seq: str) -> float:
    """
    Calculate primer melting temperature using nearest-neighbor / standard formula.
    Tm = 64.9 + 41 * (G+C - 16.4) / length
    """
    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    length = len(seq)
    if length == 0:
        return 0.0
    if length < 14:
        return (seq.count('A') + seq.count('T')) * 2 + gc_count * 4
    return round(64.9 + 41 * (gc_count - 16.4) / length, 1)


def calculate_gc(seq: str) -> float:
    """Calculate GC percentage of a nucleotide sequence."""
    seq = seq.upper()
    if not seq:
        return 0.0
    gc = seq.count('G') + seq.count('C')
    return round((gc / len(seq)) * 100, 1)


def design_pcr_primers(
    chrom: str,
    pos: int,
    del_size: int,
    upstream_seq: str,
    downstream_seq: str,
    flank_dist: int = 150,
    primer_len: int = 21
) -> Dict[str, Any]:
    """
    Design flanking PCR primers to validate the deletion by gel electrophoresis or Sanger sequencing.
    Places Forward primer in the upstream flanking region and Reverse primer in downstream flanking.
    """
    # Forward Primer (5' -> 3') from upstream region
    if len(upstream_seq) >= flank_dist + primer_len:
        fwd_start = len(upstream_seq) - flank_dist
        fwd_seq = upstream_seq[fwd_start:fwd_start + primer_len]
        fwd_offset = flank_dist
    elif len(upstream_seq) >= primer_len:
        fwd_seq = upstream_seq[:primer_len]
        fwd_offset = len(upstream_seq)
    else:
        # Fallback sequence
        fwd_seq = "AGTCGATCGTAGCTAGCTAGC"[:primer_len]
        fwd_offset = 150

    # Reverse Primer (5' -> 3') from downstream region (reverse complemented)
    if len(downstream_seq) >= flank_dist + primer_len:
        rev_raw = downstream_seq[flank_dist:flank_dist + primer_len]
        rev_seq = reverse_complement(rev_raw)
        rev_offset = flank_dist + primer_len
    elif len(downstream_seq) >= primer_len:
        rev_raw = downstream_seq[-primer_len:]
        rev_seq = reverse_complement(rev_raw)
        rev_offset = len(downstream_seq)
    else:
        rev_seq = "CTAGCTACGATCGACTAGCTA"[:primer_len]
        rev_offset = 150

    fwd_tm = calculate_tm(fwd_seq)
    fwd_gc = calculate_gc(fwd_seq)

    rev_tm = calculate_tm(rev_seq)
    rev_gc = calculate_gc(rev_seq)

    # Calculate expected amplicon sizes
    # WT Amplicon spans: fwd_offset + del_size + rev_offset
    wt_amplicon_len = fwd_offset + del_size + rev_offset
    mut_amplicon_len = wt_amplicon_len - del_size

    # Recommended annealing temperature: Ta = min(Tm) - 3°C
    anneal_temp = round(min(fwd_tm, rev_tm) - 3.0, 1)

    return {
        "chrom": chrom,
        "pos": pos,
        "del_size": del_size,
        "fwd_primer": {
            "name": f"Del_{del_size}bp_Fwd",
            "seq": fwd_seq,
            "len": len(fwd_seq),
            "tm": fwd_tm,
            "gc": fwd_gc,
            "position": f"{chrom}:{pos - fwd_offset}"
        },
        "rev_primer": {
            "name": f"Del_{del_size}bp_Rev",
            "seq": rev_seq,
            "len": len(rev_seq),
            "tm": rev_tm,
            "gc": rev_gc,
            "position": f"{chrom}:{pos + del_size + rev_offset}"
        },
        "wt_amplicon_bp": wt_amplicon_len,
        "mut_amplicon_bp": mut_amplicon_len,
        "size_diff_bp": del_size,
        "annealing_temp_c": anneal_temp,
        "recommended_gel_agarose_pct": 2.5 if del_size < 30 else 2.0
    }


def generate_gel_html(
    wt_size: int,
    mut_size: int,
    del_size: int,
    genotype: str = "HET"
) -> str:
    """
    Render a virtual UV transilluminator Agarose Gel Electrophoresis visualization.
    Uses log-molecular-weight migration physics to position DNA bands on the gel.
    """
    # Gel migration coordinate calculation (log10 scale)
    # Higher bp = slower migration (closer to top/wells)
    # Lower bp = faster migration (closer to bottom)
    min_bp = 80
    max_bp = 1000

    def get_y_percent(bp: int) -> float:
        """Map DNA fragment size to vertical gel position % (0% = top wells, 100% = bottom)."""
        log_min = math.log10(min_bp)
        log_max = math.log10(max_bp)
        log_val = math.log10(max(min_bp, min(max_bp, bp)))
        # Fraction migrated
        frac = (log_max - log_val) / (log_max - log_min)
        # Scale between 12% (wells) and 88% (gel bottom)
        return round(14 + frac * 72, 1)

    # 100 bp DNA Ladder markers
    ladder_markers = [1000, 800, 600, 500, 400, 300, 200, 100]
    ladder_bands_html = []
    for bp in ladder_markers:
        y = get_y_percent(bp)
        is_500 = (bp == 500)
        thickness = "3px" if is_500 else "2px"
        color = "#e0f2fe" if is_500 else "#7dd3fc"
        shadow = "0 0 6px #38bdf8" if is_500 else "0 0 3px #0284c7"
        ladder_bands_html.append(
            f'<div style="position:absolute; top:{y}%; left:12px; right:12px; height:{thickness}; '
            f'background:{color}; box-shadow:{shadow}; border-radius:1px;" title="{bp} bp"></div>'
            f'<span style="position:absolute; top:{y-2}%; left:-48px; font-size:10px; color:#64748b; font-family:monospace;">{bp}bp</span>'
        )

    # Lane 2: Wild-Type Control (only WT band)
    wt_y = get_y_percent(wt_size)
    lane2_bands = (
        f'<div style="position:absolute; top:{wt_y}%; left:12px; right:12px; height:2.5px; '
        f'background:#38bdf8; box-shadow:0 0 6px #0ea5e9; border-radius:1px;" title="WT Band: {wt_size} bp"></div>'
    )

    # Lane 3: Sample (HET = WT + Mutant, HOM = Mutant only)
    mut_y = get_y_percent(mut_size)
    is_het = "HET" in genotype or "0/1" in genotype

    if is_het:
        lane3_bands = (
            f'<div style="position:absolute; top:{wt_y}%; left:12px; right:12px; height:2px; '
            f'background:#38bdf8; box-shadow:0 0 5px #0ea5e9; border-radius:1px;" title="WT Allele Band: {wt_size} bp"></div>'
            f'<div style="position:absolute; top:{mut_y}%; left:12px; right:12px; height:2px; '
            f'background:#f43f5e; box-shadow:0 0 7px #ef4444; border-radius:1px;" title="Deletion Allele Band: {mut_size} bp (-{del_size} bp)"></div>'
        )
    else:
        lane3_bands = (
            f'<div style="position:absolute; top:{mut_y}%; left:12px; right:12px; height:2.5px; '
            f'background:#f43f5e; box-shadow:0 0 7px #ef4444; border-radius:1px;" title="Homozygous Deletion Band: {mut_size} bp (-{del_size} bp)"></div>'
        )

    ladder_str = "".join(ladder_bands_html)

    html = f"""
    <div style="background:#020617; border:1px solid #1e293b; border-radius:10px; padding:20px; color:#e2e8f0; font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif;">
        <div style="display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #1e293b; padding-bottom:10px; margin-bottom:16px;">
            <div>
                <span style="font-size:16px; font-weight:700; color:#38bdf8;">🧪 Virtual Agarose Gel Electrophoresis (2.5% EtBr)</span>
                <span style="margin-left:10px; font-size:12px; color:#94a3b8;">Resolving {del_size} bp deletion band shift</span>
            </div>
            <div style="font-size:11px; color:#94a3b8;">
                WT Band: <strong style="color:#38bdf8;">{wt_size} bp</strong> | Mutant Band: <strong style="color:#f43f5e;">{mut_size} bp</strong>
            </div>
        </div>

        <div style="display:flex; justify-content:center; padding:10px 0;">
            <div style="position:relative; width:460px; height:320px; background:#030712; border:2px solid #0f172a; border-radius:8px; box-shadow:inset 0 0 25px rgba(0,180,255,0.08); margin-left:50px;">
                
                <!-- Wells Header -->
                <div style="display:flex; justify-content:space-around; padding-top:8px; border-bottom:1px dashed #1e293b;">
                    <div style="text-align:center; width:90px;">
                        <span style="font-size:11px; font-weight:bold; color:#94a3b8;">Lane 1<br><span style="color:#64748b; font-size:10px;">100bp Ladder</span></span>
                        <div style="width:50px; height:6px; background:#000; border:1px solid #334155; margin:6px auto 0 auto; border-radius:2px;"></div>
                    </div>
                    <div style="text-align:center; width:90px;">
                        <span style="font-size:11px; font-weight:bold; color:#38bdf8;">Lane 2<br><span style="color:#64748b; font-size:10px;">Wild-Type</span></span>
                        <div style="width:50px; height:6px; background:#000; border:1px solid #334155; margin:6px auto 0 auto; border-radius:2px;"></div>
                    </div>
                    <div style="text-align:center; width:90px;">
                        <span style="font-size:11px; font-weight:bold; color:#f43f5e;">Lane 3<br><span style="color:#64748b; font-size:10px;">Sample ({genotype})</span></span>
                        <div style="width:50px; height:6px; background:#000; border:1px solid #334155; margin:6px auto 0 auto; border-radius:2px;"></div>
                    </div>
                </div>

                <!-- Lane Columns -->
                <div style="position:absolute; top:45px; bottom:10px; left:20px; width:100px;">
                    {ladder_str}
                </div>
                <div style="position:absolute; top:45px; bottom:10px; left:165px; width:100px;">
                    {lane2_bands}
                    <span style="position:absolute; top:{wt_y-2}%; right:-40px; font-size:10px; color:#38bdf8; font-family:monospace;">{wt_size}bp</span>
                </div>
                <div style="position:absolute; top:45px; bottom:10px; left:305px; width:100px;">
                    {lane3_bands}
                    <span style="position:absolute; top:{mut_y-2}%; right:-40px; font-size:10px; color:#f43f5e; font-family:monospace;">{mut_size}bp</span>
                </div>
            </div>
        </div>

        <div style="display:flex; justify-content:space-between; align-items:center; margin-top:14px; font-size:12px; color:#94a3b8; border-top:1px solid #1e293b; padding-top:10px;">
            <div>
                <span>Gel Matrix: <strong>2.5% High-Resolution Agarose</strong> in 1x TAE Buffer @ 100V</span>
            </div>
            <div>
                <span style="color:#38bdf8;">■ WT Band ({wt_size} bp)</span>
                <span style="margin-left:14px; color:#f43f5e;">■ Deletion Band ({mut_size} bp, Δ{del_size} bp)</span>
            </div>
        </div>
    </div>
    """
    return html
