"""
Functional Impact & Protein Translation Predictor for Genomic Indels
Analyzes reading frame changes (Frameshift vs In-Frame), translates
nucleotide sequences into amino acid peptides, and detects premature stop codons.
"""

from typing import Dict, Any, Tuple, Optional, List


# Standard NCBI Genetic Code
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

AMINO_ACID_NAMES = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    '*': 'Stop', 'X': 'Unknown'
}

AMINO_ACID_PROPERTIES = {
    'A': 'nonpolar', 'V': 'nonpolar', 'I': 'nonpolar', 'L': 'nonpolar',
    'M': 'nonpolar', 'F': 'nonpolar', 'Y': 'polar', 'W': 'nonpolar',
    'S': 'polar', 'T': 'polar', 'N': 'polar', 'Q': 'polar', 'C': 'polar',
    'G': 'nonpolar', 'P': 'special', 'R': 'positive', 'K': 'positive',
    'H': 'positive', 'D': 'negative', 'E': 'negative', '*': 'stop'
}


def translate_sequence(dna_seq: str) -> Tuple[str, List[str]]:
    """Translate DNA sequence into amino acid string and codon list."""
    dna_upper = dna_seq.upper().replace('U', 'T')
    amino_acids = []
    codons = []
    
    for i in range(0, len(dna_upper) - 2, 3):
        codon = dna_upper[i:i+3]
        if len(codon) == 3:
            aa = CODON_TABLE.get(codon, 'X')
            amino_acids.append(aa)
            codons.append(codon)
            
    return "".join(amino_acids), codons


def analyze_functional_impact(
    del_size: int,
    ref_allele: str,
    alt_allele: str,
    upstream_seq: str = "",
    downstream_seq: str = "",
    reading_frame_offset: int = 0
) -> Dict[str, Any]:
    """
    Analyze the functional and coding impact of a deletion or SNP.
    Determines Frameshift vs In-Frame for Indels, or Missense/Synonymous/Nonsense for SNPs,
    and models peptide translation alterations.
    """
    is_snp = (len(ref_allele) == 1 and len(alt_allele) == 1 and ref_allele != alt_allele)

    if is_snp:
        # Construct CDS window for SNP
        wt_dna = upstream_seq + ref_allele + downstream_seq
        mut_dna = upstream_seq + alt_allele + downstream_seq
        wt_protein, wt_codons = translate_sequence(wt_dna)
        mut_protein, mut_codons = translate_sequence(mut_dna)

        mut_first_stop = mut_protein.find('*')
        wt_first_stop = wt_protein.find('*')

        mut_codon_idx = len(upstream_seq) // 3
        wt_aa = wt_protein[mut_codon_idx] if mut_codon_idx < len(wt_protein) else "?"
        mut_aa = mut_protein[mut_codon_idx] if mut_codon_idx < len(mut_protein) else "?"

        has_premature_stop = False
        if wt_aa == mut_aa:
            consequence = f"Synonymous (Silent) Mutation ({wt_aa}➔{mut_aa})"
            impact_level = "LOW (Silent / Preserved)"
            impact_color = "#10b981"  # Emerald
            description = f"Point mutation {ref_allele}➔{alt_allele} preserves the translated amino acid ({wt_aa}). Reading frame is intact."
        elif mut_aa == '*':
            consequence = f"Nonsense Mutation (p.{wt_aa}{mut_codon_idx+1}*)"
            impact_level = "HIGH (Premature Truncation)"
            impact_color = "#ef4444"  # Red
            has_premature_stop = True
            description = f"Point mutation {ref_allele}➔{alt_allele} creates a premature stop codon at residue {mut_codon_idx+1}, leading to truncated protein."
        else:
            wt_name = AMINO_ACID_NAMES.get(wt_aa, wt_aa)
            mut_name = AMINO_ACID_NAMES.get(mut_aa, mut_aa)
            consequence = f"Missense Mutation (p.{wt_aa}{mut_codon_idx+1}{mut_aa})"
            impact_level = "MODERATE (Amino Acid Change)"
            impact_color = "#f59e0b"  # Amber
            description = f"Point mutation {ref_allele}➔{alt_allele} causes amino acid substitution {wt_aa} ({wt_name}) ➔ {mut_aa} ({mut_name}) at codon {mut_codon_idx+1}."

        return {
            "del_size": 0,
            "is_snp": True,
            "is_frameshift": False,
            "shift_offset": 0,
            "consequence": consequence,
            "impact_level": impact_level,
            "impact_color": impact_color,
            "description": description,
            "wt_protein": wt_protein,
            "mut_protein": mut_protein,
            "wt_codons": wt_codons,
            "mut_codons": mut_codons,
            "has_premature_stop": has_premature_stop,
            "premature_stop_pos": mut_first_stop + 1 if mut_first_stop != -1 else None
        }

    is_frameshift = (del_size % 3 != 0)
    shift_offset = del_size % 3

    if is_frameshift:
        consequence = "Frameshift Variant"
        impact_level = "HIGH (Loss of Function)"
        impact_color = "#ef4444"  # Red
        description = (
            f"Deletion length ({del_size} bp) is not a multiple of 3 (remainder: {shift_offset}). "
            f"This disrupts the downstream triplet reading frame, scrambling all subsequent amino acids "
            f"and typically resulting in a premature termination codon (nonsense truncation)."
        )
    else:
        num_aa_lost = del_size // 3
        consequence = f"In-Frame Deletion ({num_aa_lost} aa lost)"
        impact_level = "MODERATE"
        impact_color = "#f59e0b"  # Amber
        description = (
            f"Deletion length ({del_size} bp) is an exact multiple of 3. "
            f"Removes {num_aa_lost} intact amino acid(s) while preserving the downstream reading frame."
        )

    # Reconstruct CDS window for translation
    anchor = ref_allele[0] if ref_allele else ""
    deleted_bases = ref_allele[1:] if len(ref_allele) > 1 else ""

    # Adjust upstream to align with triplet reading frame
    pad_len = reading_frame_offset % 3
    active_upstream = upstream_seq[pad_len:] if len(upstream_seq) > pad_len else upstream_seq

    wt_dna = active_upstream + anchor + deleted_bases + downstream_seq
    mut_dna = active_upstream + anchor + downstream_seq

    wt_protein, wt_codons = translate_sequence(wt_dna)
    mut_protein, mut_codons = translate_sequence(mut_dna)

    # Detect premature stop codons in mutant protein
    wt_first_stop = wt_protein.find('*')
    mut_first_stop = mut_protein.find('*')
    
    has_premature_stop = False
    if mut_first_stop != -1 and (wt_first_stop == -1 or mut_first_stop < wt_first_stop):
        has_premature_stop = True

    return {
        "del_size": del_size,
        "is_snp": False,
        "is_frameshift": is_frameshift,
        "shift_offset": shift_offset,
        "consequence": consequence,
        "impact_level": impact_level,
        "impact_color": impact_color,
        "description": description,
        "wt_protein": wt_protein,
        "mut_protein": mut_protein,
        "wt_codons": wt_codons,
        "mut_codons": mut_codons,
        "has_premature_stop": has_premature_stop,
        "premature_stop_pos": mut_first_stop + 1 if mut_first_stop != -1 else None
    }


def generate_protein_viewer_html(impact_data: Dict[str, Any]) -> str:
    """Generate interactive visual peptide alignment comparing WT vs Mutant protein."""
    wt_prot = impact_data.get("wt_protein", "")
    mut_prot = impact_data.get("mut_protein", "")
    is_fs = impact_data.get("is_frameshift", True)
    impact_level = impact_data.get("impact_level", "HIGH")
    consequence = impact_data.get("consequence", "Frameshift")

    def aa_badge(aa: str, is_mutated: bool = False) -> str:
        prop = AMINO_ACID_PROPERTIES.get(aa, 'nonpolar')
        if aa == '*':
            # Stop codon
            return (
                '<span style="display:inline-block; margin:2px; padding:3px 8px; '
                'background:#7f1d1d; color:#fee2e2; border:2px solid #ef4444; border-radius:6px; '
                'font-weight:bold; font-family:monospace;" title="STOP Codon (Termination)">🛑 STOP</span>'
            )
        
        bg_map = {
            'nonpolar': ('#1e293b', '#94a3b8'),
            'polar': ('#0c4a6e', '#38bdf8'),
            'positive': ('#14532d', '#4ade80'),
            'negative': ('#701a75', '#f472b6'),
            'special': ('#78350f', '#fbbf24')
        }
        bg, text_col = bg_map.get(prop, ('#1e293b', '#94a3b8'))

        if is_mutated:
            border = "2px solid #ef4444"
            bg = "#450a0a"
            text_col = "#fca5a5"
        else:
            border = "1px solid #334155"

        name = AMINO_ACID_NAMES.get(aa, aa)
        return (
            f'<span style="display:inline-block; margin:2px; padding:3px 7px; '
            f'background:{bg}; color:{text_col}; border:{border}; border-radius:6px; '
            f'font-weight:bold; font-family:monospace; font-size:12px;" '
            f'title="{name} ({aa}) - {prop}">{aa}</span>'
        )

    wt_badges = "".join(aa_badge(aa) for aa in wt_prot[:40])
    
    # Identify divergent region in mutant
    mut_badges_list = []
    diverged = False
    for i, aa in enumerate(mut_prot[:40]):
        if i < len(wt_prot) and aa == wt_prot[i] and not diverged:
            mut_badges_list.append(aa_badge(aa, is_mutated=False))
        else:
            diverged = True
            mut_badges_list.append(aa_badge(aa, is_mutated=True))
    mut_badges = "".join(mut_badges_list)

    status_badge_bg = "#ef4444" if is_fs else "#f59e0b"

    html = f"""
    <div style="background:#0f172a; border:1px solid #1e293b; border-radius:10px; padding:18px; margin:12px 0; color:#e2e8f0; font-family:-apple-system,BlinkMacSystemFont,Segoe UI,Roboto,sans-serif;">
        <div style="display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #1e293b; padding-bottom:12px; margin-bottom:14px;">
            <div>
                <span style="font-size:16px; font-weight:700; color:#f8fafc;">🧬 Protein Impact & Translation Analysis</span>
                <span style="margin-left:12px; background:{status_badge_bg}; color:#ffffff; padding:4px 12px; border-radius:12px; font-size:12px; font-weight:700;">
                    {consequence.upper()}
                </span>
            </div>
            <div style="font-size:12px; font-weight:600; color:#94a3b8;">
                Predicted Severity: <strong style="color:{status_badge_bg};">{impact_level}</strong>
            </div>
        </div>

        <p style="font-size:13px; color:#94a3b8; line-height:1.5; margin-bottom:16px;">
            {impact_data.get('description', '')}
        </p>

        <!-- Wild-Type Protein -->
        <div style="margin-bottom:14px;">
            <div style="font-size:11px; font-weight:700; color:#38bdf8; text-transform:uppercase; letter-spacing:0.5px; margin-bottom:6px;">
                Wild-Type Peptide Sequence (WT Frame):
            </div>
            <div style="background:#020617; padding:10px; border-radius:6px; overflow-x:auto; white-space:nowrap; border:1px solid #1e293b;">
                <span style="color:#64748b; font-size:11px; margin-right:8px;">N-term</span>
                {wt_badges}
                <span style="color:#64748b; font-size:11px; margin-left:8px;">C-term</span>
            </div>
        </div>

        <!-- Mutant Protein -->
        <div style="margin-bottom:12px;">
            <div style="font-size:11px; font-weight:700; color:#f87171; text-transform:uppercase; letter-spacing:0.5px; margin-bottom:6px;">
                Mutant Translated Peptide (With Indel Shift):
            </div>
            <div style="background:#020617; padding:10px; border-radius:6px; overflow-x:auto; white-space:nowrap; border:1px solid #1e293b;">
                <span style="color:#64748b; font-size:11px; margin-right:8px;">N-term</span>
                {mut_badges}
                <span style="color:#64748b; font-size:11px; margin-left:8px;">C-term</span>
            </div>
        </div>

        <div style="display:flex; justify-content:space-between; align-items:center; margin-top:14px; font-size:11px; color:#64748b;">
            <div>
                <span style="border-left:3px solid #ef4444; padding-left:6px; color:#fca5a5; font-weight:bold;">Red Borders</span> = Scrambled / Altered Frameshift Codons
            </div>
            <div style="display:flex; gap:8px;">
                <span><span style="color:#38bdf8;">■</span> Polar</span>
                <span><span style="color:#4ade80;">■</span> Positive</span>
                <span><span style="color:#f472b6;">■</span> Negative</span>
                <span><span style="color:#94a3b8;">■</span> Nonpolar</span>
                <span><span style="color:#ef4444;">■</span> Stop</span>
            </div>
        </div>
    </div>
    """
    return html
