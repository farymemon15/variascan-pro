# Dockerfile for Genomic 25bp Indel Detection Web Application
FROM continuumio/miniconda3:latest

LABEL maintainer="Bioinformatics Software Engineer"
LABEL description="Complete automated pipeline for detecting 25 bp genomic deletions from paired-end NGS data"

WORKDIR /app

# Install system utilities
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    procps \
    curl \
    git \
    openjdk-17-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# Copy conda environment definition
COPY environment.yml /app/environment.yml

# Create Conda environment with Bioconda tools (fastp, bwa, samtools, bcftools, gatk4)
RUN conda env create -f environment.yml && conda clean -afy

# Set environment PATH to use the conda environment
ENV PATH=/opt/conda/envs/indel25-pipeline/bin:$PATH
ENV CONDA_DEFAULT_ENV=indel25-pipeline

# Copy application source code
COPY . /app

# Expose default Streamlit port
EXPOSE 8501

# Healthcheck to ensure Streamlit server is running
HEALTHCHECK --interval=30s --timeout=10s --start-period=30s --retries=3 \
  CMD curl -f http://localhost:8501/_stcore/health || exit 1

# Entrypoint to run Streamlit app
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "indel25-pipeline", "streamlit", "run", "app.py", "--server.port=8501", "--server.address=0.0.0.0"]
