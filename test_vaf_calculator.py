"""
Unit tests for Variant Allele Frequency (VAF) & CRISPR Editing Calculator
"""

import unittest
from utils.vaf_calculator import parse_allelic_depth, calculate_wilson_ci, calculate_vaf_metrics, generate_vaf_html


class TestVAFCalculator(unittest.TestCase):

    def test_parse_allelic_depth(self):
        """Test parsing of AD string 'ref,alt'."""
        r, a, tot = parse_allelic_depth("20,20", 40)
        self.assertEqual(r, 20)
        self.assertEqual(a, 20)
        self.assertEqual(tot, 40)

        # Missing AD fallback
        r2, a2, tot2 = parse_allelic_depth("N/A", 30)
        self.assertEqual(tot2, 30)

    def test_vaf_metrics_and_classification(self):
        """Test VAF percentage and biological classifications."""
        # 50% VAF -> Heterozygous
        m_het = calculate_vaf_metrics("20,20", 40)
        self.assertEqual(m_het["vaf_pct"], 50.0)
        self.assertIn("Heterozygous", m_het["classification"])
        self.assertIn("50.0% Editing Rate", m_het["crispr_efficiency_score"])

        # 95% VAF -> Homozygous
        m_hom = calculate_vaf_metrics("2,38", 40)
        self.assertEqual(m_hom["vaf_pct"], 95.0)
        self.assertIn("Homozygous", m_hom["classification"])

        # 25% VAF -> Subclonal / CRISPR
        m_sub = calculate_vaf_metrics("30,10", 40)
        self.assertEqual(m_sub["vaf_pct"], 25.0)
        self.assertIn("CRISPR Edited Subclone", m_sub["classification"])

    def test_wilson_ci(self):
        """Test Wilson binomial confidence interval calculation."""
        low, high = calculate_wilson_ci(20, 40)
        self.assertGreater(low, 30.0)
        self.assertLess(high, 70.0)

    def test_vaf_html(self):
        """Test HTML rendering of VAF gauge."""
        metrics = calculate_vaf_metrics("19,19", 38)
        html = generate_vaf_html(metrics)
        self.assertIn("Variant Allele Frequency", html)
        self.assertIn("50.0%", html)
        self.assertIn("Supporting Deletion", html)


if __name__ == "__main__":
    unittest.main()
