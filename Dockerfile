
FROM mambaorg/micromamba:1.5.6

ENV MAMBA_DOCKERFILE_ACTIVATE=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app
COPY environment.yml /app/environment.yml

RUN micromamba env create -f /app/environment.yml && \
    micromamba clean --all -y

COPY codonopt/ codonopt/
ENV PYTHONPATH=/app

ENTRYPOINT ["/opt/conda/envs/codonopt/bin/python", "-m", "codonopt.cli"]

