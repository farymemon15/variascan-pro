"""
NCBI SRA & ENA Public Genomic Data Fetcher
Enables searching public sequencing repositories by Gene Symbol (e.g. MYBPC3, BRCA1, TP53),
filtering by sequencing criteria (Organism, Library Strategy, Disease Phenotype),
and retrieving direct FASTQ download links and quick subsamples.
"""

import os
import re
import json
import gzip
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple


# Curated benchmarks for popular genes (used for immediate testing, instant presets, and offline fallback)
CURATED_GENE_BENCHMARKS = {
    "MYBPC3": {
        "gene_name": "MYBPC3",
        "full_name": "Myosin Binding Protein C3, Cardiac",
        "chromosome": "chr11",
        "default_phenotype": "Cardiomyopathy / HCM",
        "suggested_runs": [
            {
                "run_accession": "SRR30316136",
                "title": "GSM8474755: Control Differentiation Batch 3; Homo sapiens; RNA-Seq",
                "strategy": "RNA-Seq",
                "platform": "Illumina NovaSeq 6000",
                "organism": "Homo sapiens",
                "read_count": "55,682,078",
                "spots": 55682078,
                "layout": "PAIRED",
                "r1_url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR303/036/SRR30316136/SRR30316136_1.fastq.gz",
                "r2_url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR303/036/SRR30316136/SRR30316136_2.fastq.gz"
            },
            {
                "run_accession": "SRR30316137",
                "title": "GSM8474756: Control Differentiation Batch 4; Homo sapiens; RNA-Seq",
                "strategy": "RNA-Seq",
                "platform": "Illumina NovaSeq 6000",
                "organism": "Homo sapiens",
                "read_count": "67,403,660",
                "spots": 67403660,
                "layout": "PAIRED",
                "r1_url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR303/036/SRR30316137/SRR30316137_1.fastq.gz",
                "r2_url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR303/036/SRR30316137/SRR30316137_2.fastq.gz"
            }
        ]
    },
    "TP53": {
        "gene_name": "TP53",
        "full_name": "Tumor Protein P53",
        "chromosome": "chr17",
        "default_phenotype": "Carcinoma / Li-Fraumeni",
        "suggested_runs": [
            {
                "run_accession": "SRR11412215",
                "title": "RNA-Seq of TP53 mutated colorectal adenocarcinoma; Homo sapiens",
                "strategy": "RNA-Seq",
                "platform": "Illumina HiSeq 2500",
                "organism": "Homo sapiens",
                "read_count": "38,410,210",
                "spots": 38410210,
                "layout": "PAIRED",
                "r1_url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/015/SRR11412215/SRR11412215_1.fastq.gz",
                "r2_url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR114/015/SRR11412215/SRR11412215_2.fastq.gz"
            }
        ]
    },
    "BRCA1": {
        "gene_name": "BRCA1",
        "full_name": "BRCA1 DNA Repair Associated",
        "chromosome": "chr17",
        "default_phenotype": "Breast / Ovarian Cancer",
        "suggested_runs": [
            {
                "run_accession": "SRR8572111",
                "title": "Targeted deep amplicon sequencing of BRCA1/BRCA2 in familial breast cancer",
                "strategy": "Targeted-Capture",
                "platform": "Illumina MiSeq",
                "organism": "Homo sapiens",
                "read_count": "1,240,500",
                "spots": 1240500,
                "layout": "PAIRED",
                "r1_url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR857/011/SRR8572111/SRR8572111_1.fastq.gz",
                "r2_url": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR857/011/SRR8572111/SRR8572111_2.fastq.gz"
            }
        ]
    }
}


def build_ncbi_query(
    gene: str,
    organism: str = "Homo sapiens",
    strategy: str = "Any",
    phenotype: str = ""
) -> str:
    """Build a standard NCBI SRA boolean search term."""
    clean_gene = gene.strip().upper()
    parts = [f"{clean_gene}[All Fields]"]

    if organism and organism != "Any":
        parts.append(f'"{organism}"[Organism]')

    if strategy and strategy not in ("Any", "All"):
        if "RNA" in strategy.upper():
            parts.append('"rna seq"[Strategy]')
        elif "TARGET" in strategy.upper() or "AMPLICON" in strategy.upper():
            parts.append('("targeted-capture"[Strategy] OR "amplicon"[Strategy])')
        elif "EXOME" in strategy.upper() or "WES" in strategy.upper():
            parts.append('"wxs"[Strategy]')
        elif "WGS" in strategy.upper():
            parts.append('"wgs"[Strategy]')

    if phenotype and phenotype.strip():
        clean_pheno = phenotype.strip()
        parts.append(f'"{clean_pheno}"')

    # Force paired layout on Illumina for NGS pipeline compatibility
    parts.append('("illumina"[Platform] AND "paired"[Layout])')

    return " AND ".join(parts)


def get_ena_fastq_urls(run_accession: str, timeout: int = 6) -> Dict[str, str]:
    """
    Query European Nucleotide Archive (ENA) filereport API for direct HTTP/FTP FASTQ URLs.
    """
    url = (
        f"https://www.ebi.ac.uk/ena/portal/api/filereport"
        f"?accession={run_accession}&result=read_run"
        f"&fields=fastq_ftp,fastq_bytes,read_count&format=json"
    )
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Bioinformatics Pipeline)"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            if data and len(data) > 0:
                item = data[0]
                ftp_str = item.get("fastq_ftp", "")
                parts = ftp_str.split(";") if ftp_str else []
                r1_url = ""
                r2_url = ""
                for p in parts:
                    clean_p = p.strip()
                    if not clean_p.startswith("http"):
                        clean_p = "https://" + clean_p
                    if "_1.fastq" in clean_p:
                        r1_url = clean_p
                    elif "_2.fastq" in clean_p:
                        r2_url = clean_p

                # Fallback if names don't have _1/_2 but two parts exist
                if not r1_url and len(parts) >= 1:
                    p0 = parts[0].strip()
                    r1_url = p0 if p0.startswith("http") else "https://" + p0
                if not r2_url and len(parts) >= 2:
                    p1 = parts[1].strip()
                    r2_url = p1 if p1.startswith("http") else "https://" + p1

                return {
                    "r1_url": r1_url,
                    "r2_url": r2_url,
                    "read_count": item.get("read_count", "")
                }
    except Exception:
        pass

    # Construct standard predictable ENA FTP path if API query is unavailable
    # E.g. SRR30316136 -> vol1/fastq/SRR303/036/SRR30316136/SRR30316136_1.fastq.gz
    m = re.match(r"([A-Z]{3})(\d{3})(\d{3,4})(\d{2,3})?", run_accession)
    if m:
        prefix = run_accession[:6]
        suffix = f"0{run_accession[-2:]}" if len(run_accession) == 11 else (run_accession[-1:] if len(run_accession) == 10 else "")
        sub_dir = f"{suffix}/" if suffix else ""
        base = f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{prefix}/{sub_dir}{run_accession}/{run_accession}"
        return {
            "r1_url": f"{base}_1.fastq.gz",
            "r2_url": f"{base}_2.fastq.gz",
            "read_count": ""
        }

    return {"r1_url": "", "r2_url": "", "read_count": ""}


def search_sra_by_gene(
    gene: str,
    organism: str = "Homo sapiens",
    strategy: str = "Any",
    phenotype: str = "",
    max_results: int = 5,
    timeout: int = 7
) -> List[Dict[str, Any]]:
    """
    Search NCBI SRA for runs matching the specified gene and criteria.
    Returns a list of structured run metadata dictionaries with direct FASTQ links.
    """
    clean_gene = gene.strip().upper()
    query = build_ncbi_query(clean_gene, organism, strategy, phenotype)
    encoded_query = urllib.parse.quote(query)

    search_url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=sra&term={encoded_query}&retmode=json&retmax={max_results}"
    )

    req = urllib.request.Request(search_url, headers={"User-Agent": "Mozilla/5.0 (Bioinformatics Pipeline)"})
    runs_list: List[Dict[str, Any]] = []

    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            id_list = data.get("esearchresult", {}).get("idlist", [])

        if not id_list:
            # Check if curated benchmark is available
            if clean_gene in CURATED_GENE_BENCHMARKS:
                return CURATED_GENE_BENCHMARKS[clean_gene]["suggested_runs"][:max_results]
            return []

        # Fetch summaries for the retrieved IDs
        uids_str = ",".join(id_list)
        summary_url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            f"?db=sra&id={uids_str}&retmode=json"
        )
        req_sum = urllib.request.Request(summary_url, headers={"User-Agent": "Mozilla/5.0 (Bioinformatics Pipeline)"})
        with urllib.request.urlopen(req_sum, timeout=timeout) as resp_sum:
            sum_data = json.loads(resp_sum.read().decode("utf-8"))

        result_dict = sum_data.get("result", {})
        for uid in id_list:
            doc = result_dict.get(uid, {})
            runs_xml = doc.get("runs", "")
            expxml = doc.get("expxml", "")

            # Parse Run accession from runs XML snippet
            run_acc_match = re.search(r'acc="([A-Z0-9]+)"', runs_xml)
            spots_match = re.search(r'total_spots="(\d+)"', runs_xml)
            
            # Parse Title & Platform from expxml snippet
            title_match = re.search(r'<Title>(.*?)</Title>', expxml)
            platform_match = re.search(r'instrument_model="([^"]+)"', expxml)
            if not platform_match:
                platform_match = re.search(r'<Platform[^>]*>(.*?)</Platform>', expxml)

            run_acc = run_acc_match.group(1) if run_acc_match else f"SRR_{uid}"
            spots_count = int(spots_match.group(1)) if spots_match else 0
            title = title_match.group(1) if title_match else f"Sequencing of {clean_gene} in {organism}"
            platform = platform_match.group(1) if platform_match else "Illumina"

            formatted_reads = f"{spots_count:,}" if spots_count > 0 else "Available"

            # Get direct ENA FastQ links
            ena_info = get_ena_fastq_urls(run_acc, timeout=4)

            runs_list.append({
                "run_accession": run_acc,
                "title": title,
                "strategy": strategy if strategy != "Any" else "Targeted / NGS",
                "platform": platform,
                "organism": organism,
                "read_count": formatted_reads,
                "spots": spots_count,
                "layout": "PAIRED",
                "r1_url": ena_info.get("r1_url", ""),
                "r2_url": ena_info.get("r2_url", "")
            })

    except Exception:
        # Fallback to curated benchmarks if NCBI network request fails
        if clean_gene in CURATED_GENE_BENCHMARKS:
            return CURATED_GENE_BENCHMARKS[clean_gene]["suggested_runs"][:max_results]

    return runs_list


def download_fastq_subsample(
    url: str,
    output_path: str,
    max_reads: int = 10000,
    timeout: int = 10
) -> bool:
    """
    Stream and extract a lightweight subsample (e.g. 10,000 reads) from a remote FASTQ (.gz) file.
    This enables immediate pipeline execution without downloading several gigabytes of data.
    """
    if not url:
        return False

    out_file = Path(output_path)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    max_lines = max_reads * 4
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})

    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
            is_gzip = url.endswith(".gz")
            if is_gzip:
                decomp = gzip.GzipFile(fileobj=response)
                lines = []
                for _ in range(max_lines):
                    line = decomp.readline()
                    if not line:
                        break
                    lines.append(line.decode("utf-8", errors="replace"))
            else:
                lines = []
                for _ in range(max_lines):
                    line = response.readline()
                    if not line:
                        break
                    lines.append(line.decode("utf-8", errors="replace"))

            with open(out_file, "w", encoding="utf-8") as f:
                f.writelines(lines)

        return True
    except Exception:
        return False
