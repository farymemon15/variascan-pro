"""
Automated Test Suite for Genomic 25 bp Indel Pipeline
Tests demo data generation, VCF parsing, exact 25 bp indel filtering,
alignment visualization, PDF report generation, and input validation.
"""

import os
import shutil
import tempfile
import unittest
from pathlib import Path
import pandas as pd

from utils.demo_data import generate_demo_dataset, reverse_complement
from utils.vcf_parser import VCFParser, VCFRecord
from utils.visualizer import generate_ascii_alignment, generate_html_viewer, extract_flanking_sequence
from utils.report_generator import generate_pdf_report
from pipeline import IndelPipeline, PipelineConfig


class TestGenomicPipeline(unittest.TestCase):

    def setUp(self):
        self.test_dir = tempfile.mkdtemp(prefix="test_indel_")

    def tearDown(self):
        shutil.rmtree(self.test_dir, ignore_errors=True)

    def test_demo_data_generator(self):
        """Verify synthetic data generation creates valid fasta, fai, and fastq files."""
        meta = generate_demo_dataset(
            output_dir=self.test_dir,
            target_deletion_bp=25,
            genome_length=2000,
            read_length=100,
            read_depth=20,
            compress_fastq=False
        )

        ref_path = Path(meta["reference_path"])
        fai_path = Path(meta["fai_path"])
        r1_path = Path(meta["r1_path"])
        r2_path = Path(meta["r2_path"])

        self.assertTrue(ref_path.exists(), "Reference FASTA should exist")
        self.assertTrue(fai_path.exists(), "FASTA index (.fai) should exist")
        self.assertTrue(r1_path.exists(), "R1 FASTQ should exist")
        self.assertTrue(r2_path.exists(), "R2 FASTQ should exist")

        # Check .fai format
        fai_content = fai_path.read_text().strip().split("\t")
        self.assertEqual(fai_content[0], "chr1")
        self.assertEqual(int(fai_content[1]), 2000)

        # Check FASTQ line counts match
        r1_lines = r1_path.read_text().strip().split("\n")
        r2_lines = r2_path.read_text().strip().split("\n")
        self.assertEqual(len(r1_lines), len(r2_lines))
        self.assertEqual(len(r1_lines) % 4, 0, "FASTQ must have 4 lines per record")

        # Ground truth check
        truth = meta["truth_deletion"]
        self.assertEqual(truth["del_size"], 25)
        self.assertEqual(len(truth["deleted_sequence"]), 25)

    def test_vcf_parsing_and_exact_filtering(self):
        """Test strict isolation of 25 bp deletions: TYPE='indel' && len(REF)-len(ALT)==25."""
        vcf_content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t100\t.\tA\tG\t90.0\tPASS\tDP=30\tGT:DP:AD\t0/1:30:15,15
chr1\t200\t.\tCGTAGC\tC\t95.0\tPASS\tINDEL;DP=25\tGT:DP:AD\t0/1:25:12,13
chr1\t300\t.\tCGCTAGCTAGCTAGCTAGCTAGCTAG\tC\t99.5\tPASS\tINDEL;DP=40\tGT:DP:AD\t0/1:40:20,20
chr1\t400\t.\tC\tCGTAGCTAGCTA\t80.0\tPASS\tINDEL;DP=20\tGT:DP:AD\t1/1:20:0,20
"""
        vcf_path = Path(self.test_dir) / "test.vcf"
        vcf_path.write_text(vcf_content)

        parser = VCFParser(str(vcf_path))
        self.assertEqual(len(parser.records), 4)

        # Record 1: SNV (A -> G)
        self.assertFalse(parser.records[0].is_deletion)
        self.assertEqual(parser.records[0].deletion_size, 0)

        # Record 2: 5 bp deletion (CGTAGC -> C)
        self.assertTrue(parser.records[1].is_deletion)
        self.assertEqual(parser.records[1].deletion_size, 5)

        # Record 3: 25 bp deletion (len=26 -> len=1, diff=25)
        self.assertTrue(parser.records[2].is_deletion)
        self.assertEqual(parser.records[2].deletion_size, 25)
        self.assertEqual(parser.records[2].genotype, "HET (0/1)")
        self.assertEqual(parser.records[2].read_depth, 40)

        # Record 4: Insertion (C -> CGTAGCTAGCTA)
        self.assertTrue(parser.records[3].is_insertion)
        self.assertFalse(parser.records[3].is_deletion)

        # Test exact 25 bp filter
        del25 = parser.get_target_deletions(25)
        self.assertEqual(len(del25), 1)
        self.assertEqual(del25[0].pos, 300)
        self.assertEqual(del25[0].deletion_size, 25)

        # Test DataFrame conversion
        df = parser.to_dataframe(del25)
        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]["Position"], 300)
        self.assertEqual(df.iloc[0]["Deletion Size (bp)"], 25)

    def test_visualizer(self):
        """Test ASCII and HTML visualizer generators."""
        ref_allele = "CGCTAGCTAGCTAGCTAGCTAGCTAG"  # Anchor C + 25 bp
        alt_allele = "C"
        ascii_out = generate_ascii_alignment(
            chrom="chr1",
            pos=2000,
            ref_allele=ref_allele,
            alt_allele=alt_allele,
            upstream="AAAACCCCGGGG",
            downstream="TTTTCCCCGGGG"
        )
        self.assertIn("GENOMIC DELETION LOCUS: chr1:2000", ascii_out)
        self.assertIn("25 bp DELETION", ascii_out)
        self.assertIn("-------------------------", ascii_out)

        html_out = generate_html_viewer(
            chrom="chr1",
            pos=2000,
            ref_allele=ref_allele,
            alt_allele=alt_allele,
            upstream="AAAA",
            downstream="TTTT"
        )
        self.assertIn("25 bp DELETION", html_out)
        self.assertIn("chr1:2,000", html_out)

    def test_pdf_report_generation(self):
        """Test publication-quality PDF report creation with ReportLab."""
        pdf_path = Path(self.test_dir) / "test_report.pdf"
        sample_df = pd.DataFrame([{
            "Chromosome": "chr1",
            "Position": 2000,
            "Reference Allele": "CGCTAGCTAGCTAGCTAGCTAGCTA",
            "Alternate Allele": "C",
            "Deletion Size (bp)": 25,
            "Genotype (HET/HOM)": "HET (0/1)",
            "Quality Score (QUAL)": 99.9,
            "Read Depth (DP)": 38,
            "Allelic Depth (AD)": "19,19",
            "Filter": "PASS"
        }])

        summary_stats = {
            "reference": "reference.fasta",
            "r1_sample": "sample_R1.fastq",
            "r2_sample": "sample_R2.fastq",
            "caller": "bcftools",
            "total_raw_variants": 3
        }

        generated_pdf = generate_pdf_report(
            str(pdf_path),
            target_bp=25,
            summary_stats=summary_stats,
            deletions_df=sample_df,
            ascii_alignment="ASCII DIAGRAM"
        )

        self.assertTrue(Path(generated_pdf).exists())
        self.assertGreater(os.path.getsize(generated_pdf), 1000, "PDF file should have valid non-empty size")

    def test_pipeline_validation(self):
        """Test pipeline input validation rules."""
        config = PipelineConfig(
            reference_path="non_existent.fasta",
            r1_path="non_existent_R1.fastq",
            r2_path="non_existent_R2.fastq",
            output_dir=self.test_dir,
            target_deletion_bp=-5
        )
        pipeline = IndelPipeline(config)
        errors = pipeline.validate_inputs()
        self.assertTrue(any("Reference genome file does not exist" in e for e in errors))
        self.assertTrue(any("Forward reads file" in e for e in errors))
        self.assertTrue(any("Target deletion length must be positive" in e for e in errors))

    def test_end_to_end_pipeline_simulation(self):
        """Test end-to-end pipeline execution yielding events and generating all outputs."""
        meta = generate_demo_dataset(
            output_dir=os.path.join(self.test_dir, "demo_data"),
            target_deletion_bp=25,
            genome_length=2000,
            read_length=100,
            read_depth=20
        )

        pipeline_out = os.path.join(self.test_dir, "pipeline_out")
        config = PipelineConfig(
            reference_path=meta["reference_path"],
            r1_path=meta["r1_path"],
            r2_path=meta["r2_path"],
            output_dir=pipeline_out,
            target_deletion_bp=25,
            force_sim=True
        )

        pipeline = IndelPipeline(config)
        events = list(pipeline.execute())
        self.assertGreater(len(events), 5)

        # Verify output files
        csv_file = Path(pipeline_out) / "filtered_deletions.csv"
        vcf_file = Path(pipeline_out) / "filtered_deletions.vcf"
        pdf_file = Path(pipeline_out) / "filtered_deletions_report.pdf"

        self.assertTrue(csv_file.exists(), "filtered_deletions.csv must be generated")
        self.assertTrue(vcf_file.exists(), "filtered_deletions.vcf must be generated")
        self.assertTrue(pdf_file.exists(), "filtered_deletions_report.pdf must be generated")

        # Verify CSV content
        df = pd.read_csv(csv_file)
        self.assertEqual(len(df), 1)
        self.assertEqual(df.iloc[0]["Deletion Size (bp)"], 25)
        self.assertEqual(df.iloc[0]["Chromosome"], "chr1")

    def test_pipeline_snp_mode(self):
        """Test pipeline execution configured for SNPs only."""
        pipeline_out = os.path.join(self.test_dir, "snp_out")
        meta = generate_demo_dataset(self.test_dir, target_deletion_bp=25, seed=42)

        config = PipelineConfig(
            reference_path=meta["reference_path"],
            r1_path=meta["r1_path"],
            r2_path=meta["r2_path"],
            output_dir=pipeline_out,
            variant_type="snp",
            force_sim=True
        )

        pipeline = IndelPipeline(config)
        events = list(pipeline.execute())
        self.assertGreater(len(events), 5)

        csv_file = Path(pipeline_out) / "filtered_deletions.csv"
        self.assertTrue(csv_file.exists())
        df = pd.read_csv(csv_file)
        self.assertEqual(len(df), 2)
        # Verify that both detected variants are SNPs
        for _, row in df.iterrows():
            self.assertEqual(row["Variant Type"], "SNP")
            self.assertEqual(row["Deletion Size (bp)"], 0)

    def test_pipeline_all_variants_mode(self):
        """Test pipeline execution configured for all variants (SNPs + Indels)."""
        pipeline_out = os.path.join(self.test_dir, "all_out")
        meta = generate_demo_dataset(self.test_dir, target_deletion_bp=25, seed=42)

        config = PipelineConfig(
            reference_path=meta["reference_path"],
            r1_path=meta["r1_path"],
            r2_path=meta["r2_path"],
            output_dir=pipeline_out,
            variant_type="all",
            force_sim=True
        )

        pipeline = IndelPipeline(config)
        events = list(pipeline.execute())
        self.assertGreater(len(events), 5)

        csv_file = Path(pipeline_out) / "filtered_deletions.csv"
        self.assertTrue(csv_file.exists())
        df = pd.read_csv(csv_file)
        # Should contain both SNPs and the target deletion (3 records: 1 deletion, 2 SNPs)
        self.assertEqual(len(df), 3)
        variant_types = set(df["Variant Type"].unique())
        self.assertIn("SNP", variant_types)
        self.assertIn("Deletion", variant_types)

    def test_fastq_bundle_zip_packaging(self):
        """Verify FASTQ files can be bundled into a valid ZIP archive for system download."""
        import io
        import zipfile
        meta = generate_demo_dataset(self.test_dir, target_deletion_bp=25, seed=42)

        zip_buf = io.BytesIO()
        with zipfile.ZipFile(zip_buf, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.write(meta["r1_path"], arcname=Path(meta["r1_path"]).name)
            zf.write(meta["r2_path"], arcname=Path(meta["r2_path"]).name)
            zf.write(meta["reference_path"], arcname=Path(meta["reference_path"]).name)

        zip_bytes = zip_buf.getvalue()
        self.assertGreater(len(zip_bytes), 1000)

        # Re-read zip and check filenames
        with zipfile.ZipFile(io.BytesIO(zip_bytes), "r") as zf:
            names = zf.namelist()
            self.assertIn("demo_sample_R1.fastq", names)
            self.assertIn("demo_sample_R2.fastq", names)
            self.assertIn("reference.fasta", names)

    def test_pipeline_results_contain_fastq_paths(self):
        """Verify pipeline execution result dictionary contains fastq and reference paths."""
        pipeline_out = os.path.join(self.test_dir, "paths_out")
        meta = generate_demo_dataset(self.test_dir, target_deletion_bp=25, seed=42)

        config = PipelineConfig(
            reference_path=meta["reference_path"],
            r1_path=meta["r1_path"],
            r2_path=meta["r2_path"],
            output_dir=pipeline_out,
            force_sim=True
        )

        pipeline = IndelPipeline(config)
        events = list(pipeline.execute())
        self.assertTrue(len(events) > 0)
        res = pipeline._build_results(start_time=0, is_sim=True)

        self.assertIn("r1_path", res)
        self.assertIn("r2_path", res)
        self.assertIn("reference_path", res)
        self.assertTrue(os.path.exists(res["r1_path"]))
        self.assertTrue(os.path.exists(res["r2_path"]))


if __name__ == "__main__":
    unittest.main()
