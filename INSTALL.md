# Installation on MacOS with M2 chip

## First environment: model creation with ModelSEEDpy

### 0. Create the environment and activate it

```{bash}
conda create -n tfm python=3.10 -y
```

### 1. Install SciPy

Follow [this answer](https://github.com/scikit-learn/scikit-learn/issues/19137#issuecomment-936169173).

```{bash}
brew install openblas lapack gfortran

conda install cython pybind11 pythran -y
conda install numpy scipy pandas -y
```

### 2. Install ModelSEEDpy

```{bash}
pip install --no-use-pep517 scikit-learn"==0.23.2"
pip install modelseedpy
```

### 3. Install other dependencies
```{bash}
pip install bioservices
```

### 4. For working interactively with Python

```{bash}
conda install jupyter ipykernel -y
python -m ipykernel install --user
```



## Second environment: community modeling with MICOM

### 0. Create the environment and activate it

```{bash}
conda create -n tfm-micom python=3.10 -y
```

### 1. Install dependencies

NOTE: reinstall gcc if already installed.

```{bash}
brew install cmake gcc
```

NOTE: the error "Could not find compiler set in environment variable CC: arm64-apple-darwin20.0.0-clang" is solved by setting the following environment variables:

```{bash}
export CC=gcc
export CXX=g++
```

### 2. Install MICOM
```{bash}
pip install micom
```

### 3. Install other dependencies
```{bash}
conda install openpyxl xlrd tabula-py plotly
conda install -c conda-forge odfpy
```

### 4. For working interactively with Python

```{bash}
conda install jupyter ipykernel -y
python -m ipykernel install --user
```
