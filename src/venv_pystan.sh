# Pystan environment configuration with pyenv-virtualenv
KERNEL_NAME=custom_pystan
PYTHON_VERSION=3.8.1

pyenv install "$PYTHON_VERSION"

pyenv virtualenv "$PYTHON_VERSION" "$KERNEL_NAME"

# customize for packages installation
pip install numpy scipy matplotlib seaborn openpyxl statsmodels scikit-learn pandas
pip install pystan arviz

# adding jupyter lab 
 pip install jupyterlab

