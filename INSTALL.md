# Installation on MacOS with M2 chip

## First environment: model creation with ModelSEEDpy

### 0. Create the environment and activate it

```{bash}
conda clean --all -y
conda create -n biofoundry python=3.10 -y
conda activate biofoundry
conda update --all -y
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

### 3. Install RetroPath2-wrapper

```{bash}
conda clean --all -y
conda update --all -y
conda install -c conda-forge retropath2_wrapper -y
```

A reinstallation is needed:

```{bash}
conda clean --all -y
conda update --all -y
conda install -c conda-forge retropath2_wrapper -y
```

### 4. Install other dependencies

```{bash}
conda install -c conda-forge rdkit python-kaleido plotly -y
conda install -c conda-forge tabula-py xlrd openpyxl -y
pip install odfpy
```

### 5. For working interactively with Python

```{bash}
pip install jupyter ipykernel
python -m ipykernel install --user
```

### 6. Unit testing

```{bash}
pip install -U pytest
```


## Second environment: community modeling with MICOM

### 0. Create the environment and activate it

```{bash}
conda create -n biofoundry-micom python=3.10 -y
conda activate biofoundry-micom
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
conda install -c conda-forge python-kaleido plotly -y
conda install -c conda-forge tabula-py xlrd openpyxl -y
pip install odfpy
```

Fix Java not found:
```{bash}
conda install -c conda-forge openjdk=17.0.3 -y
```

### 4. For working interactively with Python

```{bash}
conda install jupyter ipykernel -y
python -m ipykernel install --user
```

### 6. Unit testing

```{bash}
pip install -U pytest
```
