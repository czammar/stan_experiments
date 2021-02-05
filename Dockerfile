# Having fun :)
ARG BASE_CONTAINER=jupyter/scipy-notebook
FROM $BASE_CONTAINER

LABEL maintainer="Cesar Zamora"

USER root

# install some utilities
RUN apt-get update && \
    apt-get install -y htop tree nano && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

USER $NB_UID

# Install some conda packages
RUN conda install --quiet --yes pystan && \
    conda clean --all -f -y && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"