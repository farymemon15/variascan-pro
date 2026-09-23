"""
Genomic Mini-Viewer & Sequence Alignment Visualizer
Renders ASCII and nucleotide-colored HTML representations of deletion loci with flanking sequences.
"""

from pathlib import Path
from typing import Optional, Tuple


NUCLEOTIDE_COLORS = {
    'A': '#10B981',  # Emerald Green
    'C': '#3B82F6',  # Bright Blue
    'G': '#F59E0B',  # Amber Yellow
    'T': '#EF4444',  # Crimson Red
    'N': '#94A3B8',  # Slate Grey
    '-': '#64748B',  # Gap Slate
}


def extract_flanking_sequence(
    fasta_path: str,
    chrom: str,
    pos: int,
    ref_len: int,
    flank_bp: int = 25
) -> Tuple[str, str]:
    """
    Extract upstream and downstream flanking sequences from a FASTA file.
    pos is 1-indexed.
    """
    p = Path(fasta_path)
    if not p.exists():
        return ("", "")

    try:
        current_chrom = None
        seq_parts = []
        with open(p, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    if current_chrom == chrom:
                        break
                    current_chrom = name
                elif current_chrom == chrom:
                    seq_parts.append(line.upper())

        full_seq = "".join(seq_parts)
        if not full_seq:
            return ("", "")

        # 1-indexed pos conversion
        zero_pos = pos - 1
        upstream_start = max(0, zero_pos - flank_bp)
        upstream = full_seq[upstream_start:zero_pos]

        downstream_start = zero_pos + ref_len
        downstream_end = min(len(full_seq), downstream_start + flank_bp)
        downstream = full_seq[downstream_start:downstream_end]

        return (upstream, downstream)
    except Exception:
        return ("", "")


def generate_ascii_alignment(
    chrom: str,
    pos: int,
    ref_allele: str,
    alt_allele: str,
    upstream: str = "",
    downstream: str = ""
) -> str:
    """Generate clean ASCII alignment for the deletion or SNP locus."""
    # Check if SNP
    if len(ref_allele) == 1 and len(alt_allele) == 1 and ref_allele != alt_allele:
        purines = {'A', 'G'}
        pyrimidines = {'C', 'T'}
        is_ti = (ref_allele in purines and alt_allele in purines) or (ref_allele in pyrimidines and alt_allele in pyrimidines)
        snp_class = "Transition (Ti)" if is_ti else "Transversion (Tv)"
        up_str = f"...{upstream}" if upstream else "..."
        down_str = f"{downstream}..." if downstream else "..."
        ref_line = f"REF (5'->3'): {up_str} [{ref_allele}] {down_str}"
        alt_line = f"ALT (5'->3'): {up_str} [{alt_allele}] {down_str}"
        up_matches = "|" * len(up_str)
        sub_mark = "*"
        down_matches = "|" * len(down_str)
        match_line = f"ALIGNMENT:    {up_matches}  {sub_mark}  {down_matches}"
        header = f"""╔══════════════════════════════════════════════════════════════════════════════╗
║ GENOMIC SNP LOCUS: {chrom}:{pos:<16} | TYPE: {ref_allele}➔{alt_allele} ({snp_class}){' ' * max(0, 18 - len(snp_class))}║
╚══════════════════════════════════════════════════════════════════════════════╝"""
        footer = f"""
Legend:
  • Reference Allele: [{ref_allele}]
  • Substituted Allele: [{alt_allele}] ({snp_class})
  • Flanking Context: 5' ({len(upstream)} bp) and 3' ({len(downstream)} bp)
"""
        return f"{header}\n\n{ref_line}\n{match_line}\n{alt_line}\n{footer}"

    del_size = len(ref_allele) - len(alt_allele)
    if del_size <= 0:
        return f"Variant at {chrom}:{pos} is an insertion or complex variant ({ref_allele} -> {alt_allele})."

    anchor = ref_allele[0]
    deleted_bases = ref_allele[1:]
    gap = "-" * len(deleted_bases)

    up_str = f"...{upstream}" if upstream else "..."
    down_str = f"{downstream}..." if downstream else "..."

    ref_line = f"REF (5'->3'): {up_str} [{anchor}] [{deleted_bases}] {down_str}"
    alt_line = f"ALT (5'->3'): {up_str} [{anchor}] [{gap}] {down_str}"

    up_matches = "|" * len(up_str)
    anchor_match = "|"
    del_matches = " " * len(deleted_bases)
    down_matches = "|" * len(down_str)
    match_line = f"ALIGNMENT:    {up_matches}  {anchor_match}   {del_matches}  {down_matches}"

    header = f"""╔══════════════════════════════════════════════════════════════════════════════╗
║ GENOMIC DELETION LOCUS: {chrom}:{pos:<12} | SIZE: {del_size:>2} bp DELETION               ║
╚══════════════════════════════════════════════════════════════════════════════╝"""

    footer = f"""
Legend:
  • Anchor Base: [{anchor}] (Preserved nucleotide at breakpoint)
  • Deleted Target: [{deleted_bases}] ({len(deleted_bases)} bp removed)
  • Alternate Allele: [{gap}] (Gap in sequence)
  • Flanking Context: 5' ({len(upstream)} bp) and 3' ({len(downstream)} bp)
"""

    return f"{header}\n\n{ref_line}\n{match_line}\n{alt_line}\n{footer}"


def generate_html_viewer(
    chrom: str,
    pos: int,
    ref_allele: str,
    alt_allele: str,
    upstream: str = "",
    downstream: str = "",
    genotype: str = "HET",
    qual: float = 99.0,
    dp: int = 35
) -> str:
    """Render a styled HTML nucleotide alignment widget with color-coded bases."""
    is_snp = (len(ref_allele) == 1 and len(alt_allele) == 1 and ref_allele != alt_allele)
    del_size = len(ref_allele) - len(alt_allele) if len(ref_allele) > len(alt_allele) else 0

    def base_badge(b: str, is_deleted: bool = False, is_gap: bool = False, is_sub: bool = False) -> str:
        color = NUCLEOTIDE_COLORS.get(b, '#94A3B8')
        if is_deleted:
            return (
                f'<span style="display:inline-block; margin:1px; padding:2px 5px; '
                f'background-color:#7f1d1d; color:#fca5a5; font-weight:bold; font-family:monospace; '
                f'border:1px solid #ef4444; border-radius:3px; text-decoration:line-through;" '
                f'title="Deleted Base: {b}">{b}</span>'
            )
        elif is_gap:
            return (
                f'<span style="display:inline-block; margin:1px; padding:2px 5px; '
                f'background-color:#1e293b; color:#94a3b8; font-weight:bold; font-family:monospace; '
                f'border:1px dashed #475569; border-radius:3px;" title="Deletion Gap">-</span>'
            )
        elif is_sub:
            return (
                f'<span style="display:inline-block; margin:1px; padding:2px 7px; '
                f'background-color:#78350f; color:#fde68a; font-weight:bold; font-family:monospace; '
                f'border:1.5px solid #f59e0b; border-radius:4px; box-shadow:0 0 10px rgba(245, 158, 11, 0.4);" '
                f'title="SNP Substitution: {b}">{b}</span>'
            )
        else:
            return (
                f'<span style="display:inline-block; margin:1px; padding:2px 5px; '
                f'background-color:#0f172a; color:{color}; font-weight:bold; font-family:monospace; '
                f'border:1px solid {color}44; border-radius:3px;" title="{b}">{b}</span>'
            )

    up_badges = "".join(base_badge(b) for b in upstream)
    down_badges = "".join(base_badge(b) for b in downstream)

    if is_snp:
        badge_label = f"SNP: {ref_allele}➔{alt_allele}"
        badge_bg = "#d97706"
        ref_middle = f'<span style="background:#0369a1; padding:2px 6px; border-radius:4px; font-weight:bold; border:1px solid #38bdf8;">{ref_allele}</span>'
        alt_middle = base_badge(alt_allele, is_sub=True)
        legend_detail = f"Point Mutation: <code style='color:#fde68a; background:#78350f; padding:2px 6px; border-radius:4px;'>{ref_allele} ➔ {alt_allele}</code>"
    else:
        badge_label = f"{del_size} bp DELETION"
        badge_bg = "#dc2626"
        anchor = ref_allele[0] if ref_allele else ""
        deleted_bases = ref_allele[1:] if len(ref_allele) > 1 else ""
        anchor_ref = base_badge(anchor)
        anchor_alt = base_badge(anchor)
        del_badges = "".join(base_badge(b, is_deleted=True) for b in deleted_bases)
        gap_badges = "".join(base_badge('-', is_gap=True) for _ in deleted_bases)
        ref_middle = f'{anchor_ref} <span style="background:#450a0a; padding:3px; border-radius:4px; margin:0 2px;">{del_badges}</span>'
        alt_middle = f'{anchor_alt} <span style="background:#1e1e38; padding:3px; border-radius:4px; margin:0 2px;">{gap_badges}</span>'
        legend_detail = f"Deleted Nucleotides: <code style='color:#f87171; background:#450a0a; padding:2px 6px; border-radius:4px; font-weight:bold;'>{deleted_bases}</code>"

    html = f"""
    <div style="background:#0b1329; border:1px solid #1e293b; border-radius:10px; padding:18px; margin:12px 0; color:#e2e8f0; font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif;">
        <div style="display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #1e293b; padding-bottom:10px; margin-bottom:14px;">
            <div>
                <span style="font-size:18px; font-weight:700; color:#38bdf8;">🧬 {chrom}:{pos:,}</span>
                <span style="margin-left:12px; background:{badge_bg}; color:#fee2e2; padding:3px 10px; border-radius:12px; font-size:12px; font-weight:600;">
                    {badge_label}
                </span>
                <span style="margin-left:8px; background:#3b82f6; color:#eff6ff; padding:3px 10px; border-radius:12px; font-size:12px; font-weight:600;">
                    {genotype}
                </span>
            </div>
            <div style="font-size:13px; color:#94a3b8;">
                <span style="margin-right:12px;">Depth: <strong style="color:#f8fafc;">{dp}x</strong></span>
                <span>QUAL: <strong style="color:#f8fafc;">{qual}</strong></span>
            </div>
        </div>

        <!-- Reference Allele -->
        <div style="margin-bottom:12px;">
            <div style="font-size:12px; font-weight:600; color:#94a3b8; margin-bottom:4px; text-transform:uppercase; letter-spacing:0.5px;">
                Reference Sequence (WT)
            </div>
            <div style="background:#020617; padding:10px; border-radius:6px; overflow-x:auto; white-space:nowrap; border:1px solid #1e293b;">
                <span style="color:#64748b; font-size:11px; margin-right:8px;">5'</span>
                {up_badges}
                <span style="border-left:2px solid #38bdf8; margin:0 4px;"></span>
                {ref_middle}
                <span style="border-left:2px solid #38bdf8; margin:0 4px;"></span>
                {down_badges}
                <span style="color:#64748b; font-size:11px; margin-left:8px;">3'</span>
            </div>
        </div>

        <!-- Alternate Allele -->
        <div style="margin-bottom:8px;">
            <div style="font-size:12px; font-weight:600; color:#94a3b8; margin-bottom:4px; text-transform:uppercase; letter-spacing:0.5px;">
                Alternate Sequence ({badge_label})
            </div>
            <div style="background:#020617; padding:10px; border-radius:6px; overflow-x:auto; white-space:nowrap; border:1px solid #1e293b;">
                <span style="color:#64748b; font-size:11px; margin-right:8px;">5'</span>
                {up_badges}
                <span style="border-left:2px solid #38bdf8; margin:0 4px;"></span>
                {alt_middle}
                <span style="border-left:2px solid #38bdf8; margin:0 4px;"></span>
                {down_badges}
                <span style="color:#64748b; font-size:11px; margin-left:8px;">3'</span>
            </div>
        </div>

        <!-- Legend and Summary -->
        <div style="display:flex; justify-content:space-between; align-items:center; margin-top:14px; font-size:12px; color:#64748b;">
            <div>
                {legend_detail}
            </div>
            <div style="display:flex; gap:10px;">
                <span><span style="color:#10b981;">■</span> A</span>
                <span><span style="color:#3b82f6;">■</span> C</span>
                <span><span style="color:#f59e0b;">■</span> G</span>
                <span><span style="color:#ef4444;">■</span> T</span>
                <span><span style="color:#78350f;">■</span> SNP Sub</span>
            </div>
        </div>
    </div>
    """
    return html

