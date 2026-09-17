FROM condaforge/miniforge3:latest

COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && \
    mamba clean --all --yes && \
    rm /tmp/environment.yml

ENV PATH="/opt/conda/envs/pooled-profiling/bin:${PATH}"
