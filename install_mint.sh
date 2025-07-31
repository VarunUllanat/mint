#!/bin/bash

conda env create --name mint --file=environment.yml

# Ensure pip operations run inside the 'mint' environment
conda run -n mint pip uninstall -y torch
conda run -n mint pip install torch
conda run -n mint pip install seaborn
conda run -n mint pip install -e .