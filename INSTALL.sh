# Installation on MacOS with M1/2 architectures

# SciPy dependencies
brew install openblas lapack gfortran

conda create -n tfm python=3.10 -y && conda activate tfm

# Follow this answer: https://github.com/scikit-learn/scikit-learn/issues/19137#issuecomment-936169173
conda install cython pybind11 pythran -y
conda install numpy scipy pandas -y

# Required version for ModelSEEDPy
pip install --no-use-pep517 scikit-learn"==0.23.2"
pip install modelseedpy

# For interactively working with Python
conda install jupyter ipykernel -y
python -m ipykernel install --user
