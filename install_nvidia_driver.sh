#!/bin/bash

echo "Checking if NVIDIA driver is installed..."
if dpkg -s cuda-drivers >/dev/null 2>&1; then
    echo "Skipping..."
else
    echo "Installing NVIDIA driver..."
    sudo apt update -qq
    sudo apt install -y gcc -qq

    export distro=ubuntu2404
    export arch=x86_64

    wget https://developer.download.nvidia.com/compute/cuda/repos/$distro/$arch/cuda-keyring_1.1-1_all.deb
    sudo dpkg -i cuda-keyring_1.1-1_all.deb
    rm cuda-keyring_1.1-1_all.deb
    sudo apt update -qq
    sudo apt install cuda-drivers -y -qq
    echo "Reboot your system to complete the NVIDIA driver installation!"
fi
