"""
Unit tests for Mini-IGV Coverage Tracker
"""

import unittest
from utils.coverage_tracker import generate_coverage_profile, build_coverage_chart


class TestCoverageTracker(unittest.TestCase):

    def test_coverage_profile_generation(self):
        """Test per-base coverage dataframe and depth drop inside deletion window."""
        pos = 2000
        del_size = 25
        mean_depth = 40.0
        df = generate_coverage_profile(
            chrom="chr1",
            pos=pos,
            del_size=del_size,
            mean_depth=mean_depth,
            vaf=0.50,
            window_bp=50
        )

        self.assertFalse(df.empty)
        self.assertIn("Genomic Position", df.columns)
        self.assertIn("Read Depth (x)", df.columns)

        # Inside deletion region vs flanking region
        inside_del = df[df["Genomic Position"].between(pos, pos + del_size - 1)]
        outside_del = df[~df["Genomic Position"].between(pos, pos + del_size - 1)]

        self.assertEqual(len(inside_del), del_size)
        
        # Mean depth inside deletion window should be significantly lower than flanking
        mean_inside = inside_del["Read Depth (x)"].mean()
        mean_outside = outside_del["Read Depth (x)"].mean()
        self.assertLess(mean_inside, mean_outside)

    def test_chart_construction(self):
        """Test Altair chart builds properly."""
        df = generate_coverage_profile(
            chrom="chr1",
            pos=2000,
            del_size=25,
            mean_depth=35.0
        )
        chart = build_coverage_chart(df, "chr1", 2000, 25, 35.0)
        self.assertIsNotNone(chart)


if __name__ == "__main__":
    unittest.main()
