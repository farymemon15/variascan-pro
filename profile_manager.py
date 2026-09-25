"""
Profile & Freelance Account Manager for VariaScan Pro.
Enables Admin to dynamically update contact details, Upwork, Fiverr, LinkedIn,
GitHub, and WhatsApp links without code modifications.
Persists settings to admin_profile.json.
"""

import json
import os
from pathlib import Path
from typing import Dict, Any

PROFILE_FILE = Path(__file__).resolve().parent.parent / "admin_profile.json"

DEFAULT_PROFILE = {
    "name": "Bioinformatics & Genomic Data Specialist",
    "title": "Lead Computational Biologist / NGS Pipeline Architect",
    "email": "contact@biogenomics.pro",
    "email_active": True,
    "upwork_url": "",
    "upwork_active": False,
    "fiverr_url": "",
    "fiverr_active": False,
    "linkedin_url": "",
    "linkedin_active": False,
    "github_url": "",
    "github_active": False,
    "whatsapp": "",
    "whatsapp_active": False,
    "bio": "Specialized in end-to-end Next-Generation Sequencing (NGS) pipeline development, variant calling, CRISPR off-target scoring, and clinical genomic reporting.",
    "services": [
        "Custom Automated NGS Pipelines (BWA-MEM, GATK4, fastp, BCFtools)",
        "Precision Indel & SNV Annotation (ACMG/AMP 2015 Framework)",
        "3D Protein Impact & Codon Translation Modeling",
        "Clinical & Diagnostic PDF Report Automation",
        "Cloud & Docker Containerized Bioinformatics Deployments"
    ]
}


def load_profile() -> Dict[str, Any]:
    """Load profile from admin_profile.json or return defaults."""
    if PROFILE_FILE.exists():
        try:
            with open(PROFILE_FILE, "r", encoding="utf-8") as f:
                data = json.load(f)
                # Merge with default to guarantee all keys exist
                merged = DEFAULT_PROFILE.copy()
                merged.update(data)
                return merged
        except Exception:
            return DEFAULT_PROFILE.copy()
    return DEFAULT_PROFILE.copy()


def save_profile(profile_data: Dict[str, Any]) -> bool:
    """Save updated profile to admin_profile.json."""
    try:
        with open(PROFILE_FILE, "w", encoding="utf-8") as f:
            json.dump(profile_data, f, indent=2, ensure_ascii=False)
        return True
    except Exception as e:
        print(f"Error saving profile: {e}")
        return False
