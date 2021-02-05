# pystan_experiments

Some fun with Stan and Bayesian stuffs

## Dockerbuild

```
cd ~/github/stan_experiments
docker build . -t pystan_experiments
```

## Run Jupyter Lab

```
docker run --rm -p 8888:8888 -p 4040:4040 -e JUPYTER_ENABLE_LAB=yes -v  /Users/$USER/home/jovyan/ czammar/
```

Open a web browser on localhost:8888 to enable jupyer lab