"""
Unit tests for Functional Impact & Protein Translation Predictor
"""

import unittest
from utils.impact_predictor import translate_sequence, analyze_functional_impact, generate_protein_viewer_html


class TestImpactPredictor(unittest.TestCase):

    def test_translation(self):
        """Test basic standard genetic code translation."""
        # ATG (M) GCT (A) TAA (*)
        seq = "ATGGCTTAA"
        protein, codons = translate_sequence(seq)
        self.assertEqual(protein, "MA*")
        self.assertEqual(codons, ["ATG", "GCT", "TAA"])

    def test_frameshift_detection(self):
        """Test that 25 bp is flagged as frameshift and 24 bp is flagged as in-frame."""
        # 25 bp deletion -> Frameshift
        impact_25 = analyze_functional_impact(
            del_size=25,
            ref_allele="CGCTAGCTAGCTAGCTAGCTAGCTAG",
            alt_allele="C",
            upstream_seq="ATGATGATG",
            downstream_seq="CCGCCGCCG"
        )
        self.assertTrue(impact_25["is_frameshift"])
        self.assertEqual(impact_25["shift_offset"], 1)
        self.assertIn("Frameshift", impact_25["consequence"])
        self.assertEqual(impact_25["impact_level"], "HIGH (Loss of Function)")

        # 24 bp deletion -> In-frame (8 amino acids lost)
        impact_24 = analyze_functional_impact(
            del_size=24,
            ref_allele="CGCTAGCTAGCTAGCTAGCTAGCTA",
            alt_allele="C",
            upstream_seq="ATGATGATG",
            downstream_seq="CCGCCGCCG"
        )
        self.assertFalse(impact_24["is_frameshift"])
        self.assertEqual(impact_24["shift_offset"], 0)
        self.assertIn("In-Frame", impact_24["consequence"])
        self.assertEqual(impact_24["impact_level"], "MODERATE")

    def test_protein_html_generation(self):
        """Test HTML rendering of peptide sequence alignment."""
        impact = analyze_functional_impact(
            del_size=25,
            ref_allele="CGCTAGCTAGCTAGCTAGCTAGCTAG",
            alt_allele="C",
            upstream_seq="ATGATGATG",
            downstream_seq="CCGCCGCCG"
        )
        html = generate_protein_viewer_html(impact)
        self.assertIn("Protein Impact & Translation Analysis", html)
        self.assertIn("FRAMESHIFT", html)
        self.assertIn("N-term", html)
        self.assertIn("C-term", html)


if __name__ == "__main__":
    unittest.main()
