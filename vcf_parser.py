"""
Robust VCF Parser and Indel Extractor for Genomic Deletion Analysis.
Extracts variant metrics, calculates exact deletion sizes, determines zygosity,
and structures data for UI visualization and CSV/JSON/PDF exports.
"""

import gzip
import io
import os
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import pandas as pd

from utils.region_annotator import annotate_genomic_locus


class VCFRecord:
    def __init__(
        self,
        chrom: str,
        pos: int,
        id_: str,
        ref: str,
        alt: str,
        qual: float,
        filter_: str,
        info: Dict[str, Any],
        format_keys: List[str],
        sample_values: List[str],
        raw_line: str = ""
    ):
        self.chrom = chrom
        self.pos = pos
        self.id = id_
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter_
        self.info = info
        self.format_keys = format_keys
        self.sample_values = sample_values
        self.raw_line = raw_line

        # Derived calculations
        self.is_deletion = len(self.ref) > len(self.alt)
        self.deletion_size = len(self.ref) - len(self.alt) if self.is_deletion else 0
        self.is_insertion = len(self.alt) > len(self.ref)
        self.insertion_size = len(self.alt) - len(self.ref) if self.is_insertion else 0
        self.is_snp = (len(self.ref) == 1 and len(self.alt) == 1 and self.ref != self.alt)
        self.is_indel = (self.is_deletion or self.is_insertion)
        
        # Parse sample genotype and depths
        self.format_dict = dict(zip(self.format_keys, self.sample_values)) if format_keys and sample_values else {}
        self.gt_raw = self.format_dict.get("GT", "./.")
        self.genotype = self._parse_genotype(self.gt_raw)
        
        # Read Depth (DP)
        dp_val = self.format_dict.get("DP")
        if dp_val is None or dp_val == ".":
            dp_val = self.info.get("DP", 0)
        try:
            self.read_depth = int(dp_val)
        except (ValueError, TypeError):
            self.read_depth = 0

        # Allelic Depth (AD)
        self.allelic_depth = self.format_dict.get("AD", "N/A")

        # Functional coding impact and variant classification
        if self.is_snp:
            self.variant_type = "SNP"
            purines = {'A', 'G'}
            pyrimidines = {'C', 'T'}
            if (self.ref in purines and self.alt in purines) or (self.ref in pyrimidines and self.alt in pyrimidines):
                self.snp_type = "Transition (Ti)"
            else:
                self.snp_type = "Transversion (Tv)"
            self.size_display = f"{self.ref}➔{self.alt}"
            self.is_frameshift = False
            self.consequence = f"SNV: {self.ref}➔{self.alt} ({self.snp_type})"
        elif self.is_deletion:
            self.variant_type = "Deletion"
            self.snp_type = "N/A"
            self.size_display = f"-{self.deletion_size} bp"
            self.is_frameshift = (self.deletion_size % 3 != 0)
            self.consequence = "⚠️ Frameshift" if self.is_frameshift else f"ℹ️ In-Frame (-{self.deletion_size//3} aa)"
        elif self.is_insertion:
            self.variant_type = "Insertion"
            self.snp_type = "N/A"
            self.size_display = f"+{self.insertion_size} bp"
            self.is_frameshift = (self.insertion_size % 3 != 0)
            self.consequence = "⚠️ Frameshift (Ins)" if self.is_frameshift else f"ℹ️ In-Frame (+{self.insertion_size//3} aa)"
        else:
            self.variant_type = "Complex"
            self.snp_type = "N/A"
            self.size_display = f"{len(self.ref)}>{len(self.alt)} bp"
            self.is_frameshift = False
            self.consequence = "Complex Indel"

        # Variant Allele Frequency (VAF %) calculation
        self.vaf_pct = 50.0  # Default estimate
        if isinstance(self.allelic_depth, str) and "," in self.allelic_depth:
            try:
                parts = self.allelic_depth.split(",")
                r_cnt = int(parts[0])
                a_cnt = int(parts[1])
                tot = r_cnt + a_cnt
                if tot > 0:
                    self.vaf_pct = round((a_cnt / tot) * 100.0, 1)
            except Exception:
                pass
        elif "HOM" in self.genotype:
            self.vaf_pct = 100.0
        elif "HET" in self.genotype:
            self.vaf_pct = 50.0

        # Genomic region annotation
        self.region_annotation = annotate_genomic_locus(self.chrom, self.pos, self.deletion_size)

    @staticmethod
    def _parse_genotype(gt_str: str) -> str:
        """Categorize GT string into HET, HOM, or REF."""
        gt_clean = gt_str.replace("|", "/")
        parts = gt_clean.split("/")
        if len(parts) >= 2:
            a1, a2 = parts[0], parts[1]
            if a1 == "." or a2 == ".":
                return "Unknown"
            if a1 == "0" and a2 == "0":
                return "REF (0/0)"
            if (a1 == "0" and a2 != "0") or (a1 != "0" and a2 == "0") or (a1 != a2):
                return f"HET ({gt_str})"
            if a1 == a2 and a1 != "0":
                return f"HOM ({gt_str})"
        return gt_str

    def to_dict(self) -> Dict[str, Any]:
        """Return clean dictionary representation matching prompt requirements."""
        return {
            "Chromosome": self.chrom,
            "Position": self.pos,
            "Variant Type": self.variant_type,
            "Reference Allele": self.ref,
            "Alternate Allele": self.alt,
            "Variant Size / Change": self.size_display,
            "Deletion Size (bp)": self.deletion_size,
            "Genomic Region": self.region_annotation["feature_type"],
            "Consequence": self.consequence,
            "VAF (%)": self.vaf_pct,
            "Genotype (HET/HOM)": self.genotype,
            "Quality Score (QUAL)": round(self.qual, 2) if self.qual is not None else 0.0,
            "Read Depth (DP)": self.read_depth,
            "Allelic Depth (AD)": self.allelic_depth,
            "Filter": self.filter
        }


class VCFParser:
    """Parser for VCF 4.2+ files (.vcf or .vcf.gz)."""

    def __init__(self, vcf_path: str):
        self.vcf_path = Path(vcf_path)
        self.header_lines: List[str] = []
        self.column_headers: List[str] = []
        self.records: List[VCFRecord] = []
        self._parse()

    def _get_file_handle(self):
        if str(self.vcf_path).endswith(".gz"):
            return gzip.open(self.vcf_path, "rt", encoding="utf-8", errors="replace")
        return open(self.vcf_path, "r", encoding="utf-8", errors="replace")

    def _parse(self):
        if not self.vcf_path.exists():
            raise FileNotFoundError(f"VCF file not found at: {self.vcf_path}")

        with self._get_file_handle() as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("##"):
                    self.header_lines.append(line)
                elif line.startswith("#CHROM"):
                    self.column_headers = line[1:].split("\t")
                else:
                    fields = line.split("\t")
                    if len(fields) < 8:
                        continue
                    
                    chrom = fields[0]
                    try:
                        pos = int(fields[1])
                    except ValueError:
                        continue
                    id_ = fields[2]
                    ref = fields[3].upper()
                    alt_alleles = fields[4].upper().split(",")
                    
                    try:
                        qual = float(fields[5]) if fields[5] != "." else 0.0
                    except ValueError:
                        qual = 0.0
                    
                    filter_ = fields[6]
                    info_raw = fields[7]
                    
                    # Parse info
                    info = {}
                    for item in info_raw.split(";"):
                        if "=" in item:
                            k, v = item.split("=", 1)
                            info[k] = v
                        else:
                            info[item] = True

                    format_keys = fields[8].split(":") if len(fields) > 8 else []
                    sample_vals = fields[9].split(":") if len(fields) > 9 else []

                    # Handle multiple alternate alleles
                    for alt in alt_alleles:
                        rec = VCFRecord(
                            chrom=chrom,
                            pos=pos,
                            id_=id_,
                            ref=ref,
                            alt=alt,
                            qual=qual,
                            filter_=filter_,
                            info=info,
                            format_keys=format_keys,
                            sample_values=sample_vals,
                            raw_line=line
                        )
                        self.records.append(rec)

    def get_target_deletions(self, target_bp: int = 25) -> List[VCFRecord]:
        """Strictly isolate deletions where length difference equals target_bp."""
        return [r for r in self.records if r.is_deletion and r.deletion_size == target_bp]

    def get_range_deletions(self, min_bp: int, max_bp: int) -> List[VCFRecord]:
        """Isolate deletions within length range [min_bp, max_bp]."""
        return [r for r in self.records if r.is_deletion and (min_bp <= r.deletion_size <= max_bp)]

    def get_all_deletions(self) -> List[VCFRecord]:
        """Return all detected deletions regardless of size."""
        return [r for r in self.records if r.is_deletion]

    def get_snps(self) -> List[VCFRecord]:
        """Return all detected Single Nucleotide Polymorphisms (SNPs)."""
        return [r for r in self.records if r.is_snp]

    def get_all_variants(self) -> List[VCFRecord]:
        """Return all detected variants (SNPs + Indels)."""
        return self.records

    def to_dataframe(self, records: Optional[List[VCFRecord]] = None) -> pd.DataFrame:
        """Convert records to pandas DataFrame."""
        recs = records if records is not None else self.records
        if not recs:
            return pd.DataFrame(columns=[
                "Chromosome", "Position", "Variant Type", "Reference Allele", "Alternate Allele",
                "Variant Size / Change", "Deletion Size (bp)", "Genomic Region", "Consequence",
                "VAF (%)", "Genotype (HET/HOM)", "Quality Score (QUAL)", "Read Depth (DP)",
                "Allelic Depth (AD)", "Filter"
            ])
        return pd.DataFrame([r.to_dict() for r in recs])

    def export_csv(self, output_path: str, records: Optional[List[VCFRecord]] = None) -> str:
        """Export records to CSV."""
        df = self.to_dataframe(records)
        df.to_csv(output_path, index=False)
        return output_path

    def export_vcf(self, output_path: str, records: Optional[List[VCFRecord]] = None) -> str:
        """Export records back to standard VCF format."""
        recs = records if records is not None else self.records
        with open(output_path, "w", encoding="utf-8") as f:
            for h in self.header_lines:
                f.write(f"{h}\n")
            if self.column_headers:
                f.write("#" + "\t".join(self.column_headers) + "\n")
            for r in recs:
                if r.raw_line:
                    f.write(f"{r.raw_line}\n")
                else:
                    # Synthesize line if raw_line was not preserved
                    f.write(f"{r.chrom}\t{r.pos}\t{r.id}\t{r.ref}\t{r.alt}\t{r.qual}\t{r.filter}\tDP={r.read_depth}\n")
        return output_path
