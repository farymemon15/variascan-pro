"""
Unit Tests for NCBI SRA & ENA Public Genomic Data Fetcher
Tests query generation, ENA URL resolution, curated gene fallbacks, and subsampling.
"""

import unittest
from pathlib import Path

from utils.ncbi_fetcher import (
    build_ncbi_query,
    get_ena_fastq_urls,
    search_sra_by_gene,
    download_fastq_subsample,
    CURATED_GENE_BENCHMARKS
)


class TestNCBIFetcher(unittest.TestCase):

    def test_build_ncbi_query_basic(self):
        """Test query construction with basic gene name."""
        q = build_ncbi_query("MYBPC3")
        self.assertIn("MYBPC3[All Fields]", q)
        self.assertIn('"Homo sapiens"[Organism]', q)
        self.assertIn('"illumina"[Platform]', q)
        self.assertIn('"paired"[Layout]', q)

    def test_build_ncbi_query_custom_criteria(self):
        """Test query construction with strategy, organism, and phenotype."""
        q = build_ncbi_query(
            gene="BRCA1",
            organism="Mus musculus",
            strategy="RNA-Seq",
            phenotype="Tumor"
        )
        self.assertIn("BRCA1[All Fields]", q)
        self.assertIn('"Mus musculus"[Organism]', q)
        self.assertIn('"rna seq"[Strategy]', q)
        self.assertIn('"Tumor"', q)

    def test_get_ena_fastq_urls_structure(self):
        """Test ENA URL generation produces valid HTTP/FTP FASTQ paths."""
        urls = get_ena_fastq_urls("SRR30316136")
        self.assertIn("r1_url", urls)
        self.assertIn("r2_url", urls)
        self.assertTrue(urls["r1_url"].startswith("http"))
        self.assertIn("SRR30316136", urls["r1_url"])
        self.assertIn("_1.fastq", urls["r1_url"])
        self.assertIn("_2.fastq", urls["r2_url"])

    def test_search_sra_curated_fallback(self):
        """Test that curated benchmarks return instant structured runs."""
        runs = search_sra_by_gene("MYBPC3", max_results=2)
        self.assertGreaterEqual(len(runs), 1)
        r0 = runs[0]
        self.assertIn("run_accession", r0)
        self.assertIn("r1_url", r0)
        self.assertIn("r2_url", r0)
        self.assertEqual(r0["layout"], "PAIRED")
        self.assertTrue(r0["run_accession"].startswith("SRR"))

    def test_curated_benchmarks_presence(self):
        """Verify standard genes (MYBPC3, TP53, BRCA1) have curated presets."""
        self.assertIn("MYBPC3", CURATED_GENE_BENCHMARKS)
        self.assertIn("TP53", CURATED_GENE_BENCHMARKS)
        self.assertIn("BRCA1", CURATED_GENE_BENCHMARKS)

    def test_download_subsample_invalid_url(self):
        """Verify invalid URL returns False without raising unhandled exception."""
        res = download_fastq_subsample("", "some_invalid_path.fastq")
        self.assertFalse(res)


if __name__ == "__main__":
    unittest.main()
