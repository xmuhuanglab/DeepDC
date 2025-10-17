#!/bin/bash

conda create -n cp38tf1 python=3.8 -y
conda init bash
source ~/.bashrc
conda activate cp38tf1
pip install pandas tqdm pyBigWig pyfaidx omegaconf
pip install dm-sonnet==1.11
pip install einops==0.8.1
pip install scikit-learn==1.3.2
pip install tf-slim==1.1.0
pip install xgboost==2.1.4
pip install transformers==4.26.1
pip install biopython==1.70
pip install nvidia-tensorflow==1.15.5+nv23.03 --no-deps
pip install torch==2.4.1
pip install torchvision
pip install jupyter ipykernel
python -m ipykernel install --user --name cp38tf1
