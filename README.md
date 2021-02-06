# pystan_experiments

Some fun with Stan and Bayesian stuffs

## Dockerbuild

The following is a `Dockerimage` ready to work with Pystan, based in [jupyter-docker-stacks](https://jupyter-docker-stacks.readthedocs.io/en/latest/using/selecting.html)

```
git clone https://github.com/czammar/stan_experiments.git
cd stan_experiments
docker build . -t my_pystan
```

## Run Jupyter Lab

```
# Flags enable jupyter lab and sudo super powers to install stuffs
docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -e GRANT_SUDO=yes --user root -v /Users/$USER:/home/jovyan/work my_pystan
```

Open a web browser on localhost:8888 to enable juptyer lab