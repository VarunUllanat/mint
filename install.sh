#!/bin/bash

echo "Checking if conda is installed..."
if [ -d "$HOME/miniforge3" ]; then
    echo "Skipping..."
    eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
else
    echo "Installing conda..."
    (
        cd $HOME
        curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
        bash Miniforge3-$(uname)-$(uname -m).sh -b
        eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
        conda init
        conda config --set auto_activate_base false
    )
    eval "$($HOME/miniforge3/bin/conda shell.bash hook)"
fi


echo "Checking if conda mint environment exists..."
if conda list --name mint &> /dev/null; then
    echo "Skipping..."
else
    echo "Creating conda mint environment..."
    conda env create -f environment.yml

    conda activate mint
    pip install -e .
    python -c "import mint; print('MINT installed successfully!')"
fi


echo "Checking if mint.ckpt exists..."
if [ -f "mint.ckpt" ]; then
    echo "Skipping..."
else
    echo "Downloading mint.ckpt..."
    wget https://huggingface.co/varunullanat2012/mint/resolve/main/mint.ckpt
fi


echo "Checking if HumanPPI dataset is ready..."
if [ -d "HumanPPI" ]; then
    echo "Skipping..."
else
    echo "Downloading HumanPPI dataset..."
    # Downloading HumanPPI from here: https://github.com/westlake-repl/SaProt?tab=readme-ov-file#downstream-tasks
    gdown https://drive.google.com/uc?id=1ahgj-IQTtv3Ib5iaiXO_ASh2hskEsvoX
    sudo apt install -y unp
    unp HumanPPI.tar.gz
fi


jupyter execute downstream/GeneralPPI/human-ppi/prepare_data.ipynb

