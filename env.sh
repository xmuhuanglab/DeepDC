#!/bin/bash

conda create -n cp38tf1 python=3.8 -y
conda init bash
source ~/.bashrc
conda activate cp38tf1
pip install pandas tqdm pyBigWig pyfaidx omegaconf
pip install dm-sonnet==1.11
pip install tf-slim==1.1.0
pip install scikit-learn==1.3.2
pip install xgboost==2.1.4
pip install biopython==1.70
pip install nvidia-pyindex
pip install nvidia-tensorflow
conda install -c nvidia nccl=2.8 -y
pip install torch==2.4.1
pip install torchvision
pip install einops==0.8.1
pip install transformers==4.26.1
pip install jupyter ipykernel
python -m ipykernel install --user --name cp38tf1

SCRIPT_DIR=$(dirname "$(realpath "$0")")
rm -f $SCRIPT_DIR/chopchop/config.json
echo "{" > $SCRIPT_DIR/chopchop/config.json
echo "  \"PATH\": {" >> $SCRIPT_DIR/chopchop/config.json
echo "    \"PRIMER3\": \"$SCRIPT_DIR/chopchop/primer3_core\"," >> $SCRIPT_DIR/chopchop/config.json
echo "    \"BOWTIE\": \"$SCRIPT_DIR/chopchop/bowtie/bowtie\"," >> $SCRIPT_DIR/chopchop/config.json
echo "    \"TWOBITTOFA\": \"$SCRIPT_DIR/chopchop/twoBitToFa\"," >> $SCRIPT_DIR/chopchop/config.json
echo "    \"TWOBIT_INDEX_DIR\": \"$SCRIPT_DIR/genome_ref\"," >> $SCRIPT_DIR/chopchop/config.json
echo "    \"BOWTIE_INDEX_DIR\": \"$SCRIPT_DIR/genome_ref\"," >> $SCRIPT_DIR/chopchop/config.json
echo "    \"ISOFORMS_INDEX_DIR\": \"\"," >> $SCRIPT_DIR/chopchop/config.json
echo "    \"ISOFORMS_MT_DIR\": \"\"," >> $SCRIPT_DIR/chopchop/config.json
echo "    \"GENE_TABLE_INDEX_DIR\": \"\"" >> $SCRIPT_DIR/chopchop/config.json
echo "  }," >> $SCRIPT_DIR/chopchop/config.json
echo "  \"THREADS\": 4" >> $SCRIPT_DIR/chopchop/config.json
echo "}"  >> $SCRIPT_DIR/chopchop/config.json
