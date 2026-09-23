"""
Unit tests for Genomic Region & Splice Site Annotation Engine
"""

import unittest
from utils.region_annotator import annotate_genomic_locus, generate_gene_structure_html


class TestRegionAnnotator(unittest.TestCase):

    def test_exon_annotation(self):
        """Test that position 2000 is correctly mapped to Exon 2 (CDS)."""
        res = annotate_genomic_locus("chr1", 2000, 25)
        self.assertIn("Exon 2", res["feature_type"])
        self.assertIn("Coding Exon", res["feature_type"])
        self.assertIn("Low", res["splice_risk"])

    def test_splice_junction_proximity(self):
        """Test that coordinates near exon boundaries trigger splice disruption alerts."""
        # Exon 2 starts at 1501. Position 1502 is within 2 bp (canonical junction)
        res_splice = annotate_genomic_locus("chr1", 1502, 25)
        self.assertIn("CRITICAL", res_splice["splice_risk"])

        # Intron position
        res_intron = annotate_genomic_locus("chr1", 1100, 25)
        self.assertIn("Intron 1", res_intron["feature_type"])

    def test_gene_structure_html(self):
        """Test gene structure schematic rendering."""
        annot = annotate_genomic_locus("chr1", 2000, 25)
        html = generate_gene_structure_html(annot)
        self.assertIn("Genomic Feature & Splicing Architecture", html)
        self.assertIn("Exon 1", html)
        self.assertIn("Exon 2 (CDS)", html)
        self.assertIn("Exon 3", html)
        self.assertIn("25 bp DEL", html)


if __name__ == "__main__":
    unittest.main()
