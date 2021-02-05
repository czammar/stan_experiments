# pystan_experiments

Some fun with Stan and Bayesian stuffs

## Dockerbuild

```
cd ~/github/stan_experiments
docker build . -t pystan_experiments:1.0
```

## Run Jupyter Lab

```
docker run --rm -p 8888:8888 -e JUPYTER_ENABLE_LAB=yes -e GRANT_SUDO=yes --user root -v /Users/$USER:/home/jovyan/ pystan_experiments:1.0
```

Open a web browser on localhost:8888 to enable jupyer lab