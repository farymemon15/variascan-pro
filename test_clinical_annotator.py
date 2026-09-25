import unittest
from utils.clinical_annotator import (
    generate_hgvs_notation,
    get_clinvar_annotation,
    get_gnomad_population_frequencies,
    classify_acmg_pathogenicity,
    design_sanger_validation_protocol,
    compute_alignment_qc_metrics,
    generate_bam_igv_pileup_data,
    annotate_variant_clinically
)


class TestClinicalAnnotator(unittest.TestCase):

    def test_hgvs_notation_deletion(self):
        # 25 bp deletion: anchor + 25 deleted bases = 26 bp ref
        ref = "C" + "A"*25
        alt = "C"
        hgvs_c, hgvs_p, consequence = generate_hgvs_notation(
            chrom="chr11", pos=47333092, ref=ref, alt=alt, gene_symbol="MYBPC3"
        )
        self.assertIn("c.3662_", hgvs_c)
        self.assertIn("del25", hgvs_c)
        self.assertIn("fs*", hgvs_p)
        self.assertIn("Frameshift", consequence)

    def test_hgvs_notation_snp(self):
        # Stop-Gain SNP (C>T)
        hgvs_c, hgvs_p, consequence = generate_hgvs_notation(
            chrom="chr11", pos=47333092, ref="C", alt="T", gene_symbol="MYBPC3"
        )
        self.assertIn("c.3662C>T", hgvs_c)
        self.assertIn("Ter", hgvs_p)
        self.assertEqual(consequence, "Stop-Gain (Nonsense)")

        # Synonymous SNP (A>T)
        hgvs_c2, hgvs_p2, consequence2 = generate_hgvs_notation(
            chrom="chr11", pos=47333092, ref="A", alt="T", gene_symbol="MYBPC3"
        )
        self.assertIn("=", hgvs_p2)
        self.assertEqual(consequence2, "Synonymous (Silent)")

    def test_clinvar_annotation(self):
        # MYBPC3 25 bp deletion
        res = get_clinvar_annotation("MYBPC3", 47333092, "C" + "A"*25, "C", "Frameshift Deletion")
        self.assertEqual(res["clinvar_id"], "42123")
        self.assertEqual(res["clinical_significance"], "Pathogenic")
        self.assertIn("Hypertrophic Cardiomyopathy", res["phenotype"])

    def test_gnomad_frequencies(self):
        gnomad = get_gnomad_population_frequencies("MYBPC3", "Frameshift Deletion", is_founder=True)
        self.assertAlmostEqual(gnomad["global_af"], 0.00038, places=5)
        self.assertEqual(gnomad["popmax_population"], "South Asian (SAS)")
        self.assertGreater(gnomad["popmax_af"], 0.01)

    def test_acmg_pathogenicity(self):
        # Frameshift in LoF gene with low AF
        acmg = classify_acmg_pathogenicity(
            gene_symbol="MYBPC3",
            consequence="Frameshift Deletion",
            gnomad_af=0.00001,
            clinvar_sig="Pathogenic",
            del_size=25
        )
        self.assertEqual(acmg["acmg_tier"], "Pathogenic (Class 5)")
        self.assertIn("PVS1", acmg["criteria_summary"])
        self.assertIn("PM2", acmg["criteria_summary"])

    def test_sanger_primers(self):
        sanger = design_sanger_validation_protocol(
            chrom="chr11", pos=47333092, ref="C" + "A"*25, alt="C"
        )
        self.assertEqual(len(sanger["forward_primer"]), 20)
        self.assertEqual(len(sanger["reverse_primer"]), 20)
        self.assertGreaterEqual(sanger["forward_tm"], 50.0)
        self.assertEqual(sanger["amplicon_delta_bp"], 25)

    def test_alignment_qc(self):
        qc = compute_alignment_qc_metrics(read_depth=40, qual=99.0, allelic_depth="20,20")
        self.assertEqual(qc["vaf_pct"], 50.0)
        self.assertEqual(qc["mapq"], 60)
        self.assertEqual(qc["mean_base_qual"], 38)
        self.assertIn("%", qc["wilson_ci_95"])

    def test_bam_igv_pileup(self):
        bam = generate_bam_igv_pileup_data(chrom="chr11", pos=47333092, ref="C" + "A"*25, alt="C", read_depth=10)
        self.assertEqual(len(bam["reads"]), 10)
        self.assertGreater(bam["forward_count"], 0)
        self.assertGreater(bam["reverse_count"], 0)

    def test_master_annotate_variant_clinically(self):
        res = annotate_variant_clinically(
            chrom="chr11",
            pos=47333092,
            ref="C" + "A"*25,
            alt="C",
            qual=99.0,
            read_depth=38,
            allelic_depth="19,19",
            ref_build="GRCh38",
            target_gene="MYBPC3"
        )
        self.assertEqual(res["ref_genome_build"], "GRCh38")
        self.assertIn("chr11:47,333,092", res["genomic_coordinate"])
        self.assertIn("MYBPC3", res["gene_symbol"])
        self.assertIn("Exon 33", res["exon_intron"])
        self.assertEqual(res["acmg"]["acmg_tier"], "Pathogenic (Class 5)")
        self.assertEqual(res["clinvar"]["clinical_significance"], "Pathogenic")


if __name__ == "__main__":
    unittest.main()
