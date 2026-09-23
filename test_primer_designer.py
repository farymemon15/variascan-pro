"""
Unit tests for PCR Primer Designer & Virtual Gel Simulator
"""

import unittest
from utils.primer_designer import calculate_tm, calculate_gc, design_pcr_primers, generate_gel_html


class TestPrimerDesigner(unittest.TestCase):

    def test_tm_and_gc(self):
        """Test Tm and GC calculation."""
        seq = "ATGCGATCGATCGATCGATC"  # 20 bp
        tm = calculate_tm(seq)
        gc = calculate_gc(seq)
        self.assertGreater(tm, 45.0)
        self.assertLess(tm, 70.0)
        self.assertEqual(gc, 50.0)

    def test_primer_design(self):
        """Test flanking primer generation and amplicon sizing."""
        upstream = "A" * 200 + "GCATGCATGCATGCATGCAT" + "A" * 50
        downstream = "C" * 50 + "TGCATGCATGCATGCATGCA" + "C" * 200
        del_size = 25

        res = design_pcr_primers(
            chrom="chr1",
            pos=2000,
            del_size=del_size,
            upstream_seq=upstream,
            downstream_seq=downstream,
            flank_dist=150,
            primer_len=20
        )

        self.assertIn("fwd_primer", res)
        self.assertIn("rev_primer", res)
        self.assertEqual(res["size_diff_bp"], 25)
        # Expected mutant band should be exactly 25 bp smaller than WT
        self.assertEqual(res["wt_amplicon_bp"] - res["mut_amplicon_bp"], 25)

    def test_gel_html_rendering(self):
        """Test virtual agarose gel HTML output."""
        html = generate_gel_html(wt_size=350, mut_size=325, del_size=25, genotype="HET")
        self.assertIn("Virtual Agarose Gel Electrophoresis", html)
        self.assertIn("100bp Ladder", html)
        self.assertIn("350bp", html)
        self.assertIn("325bp", html)


if __name__ == "__main__":
    unittest.main()
