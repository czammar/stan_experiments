# Start from a core stack version
FROM jupyter/datascience-notebook

ENV DEBCONF_NOWARNINGS yes

# Install from requirements.txt file
COPY --chown=${NB_UID}:${NB_GID} requirements.txt /tmp/
RUN pip install --requirement /tmp/requirements.txt && \
    fix-permissions $CONDA_DIR && \
    fix-permissions /home/$NB_USER

USER jovyan

#COPY packages.r ${PWD}
#RUN Rscript ./packages.r

RUN conda clean --all -y && rm -rf /home/jovyan/.cache/*

WORKDIR /home/jovyan/work