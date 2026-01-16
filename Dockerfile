FROM mambaorg/micromamba:1.5.6

ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PYTHONUNBUFFERED=1

# Install dependencies
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.11 \
    viennarna \
    biopython \
    pandas \
    numpy \
    tqdm && \
    micromamba clean -a -y

WORKDIR /app

COPY codonopt/ codonopt/

ENV PYTHONPATH=/app

# IMPORTANT: use micromamba run
ENTRYPOINT ["micromamba", "run", "-n", "base", "python", "-m", "codonopt.main"]

