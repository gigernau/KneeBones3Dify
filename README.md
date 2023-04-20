# KneeBones3Dify

# Setup for Ubuntu Linux OS

## 1) Update : 
	sudo apt update && sudo apt upgrade
  
## 2) Install dependencies:
    sudo apt install build-essential zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libreadline-dev libffi-dev python-tk python3-tk tk-dev python3-pil python3-pil.imagetk

## 3) Install pip3 for Python3: 
	sudo apt install python3-pip  && python3 -m pip install --upgrade pip

## 4) Install Python3 modules : 
	python3 -m pip install -r requirements.txt
	
## 5) Install [CUDA](https://developer.nvidia.com/cuda-toolkit)

## 6) Install CuPy module ( e.g., for CUDA 11.1 )
	python3 -m pip install cupy-cuda11x

## 7) Compile smoothPatch code
	gcc -shared -o smoothPatch.so smoothPatch.cpp


# USAGE EXAMPLE

## 1) Set VISDOM environment in a separate shell to view images and plots in a browser
	python3 SegOscanSST.py


[![Say Thanks!](https://img.shields.io/badge/Say%20Thanks-!-1EAEDB.svg)](https://saythanks.io/to/gianluca.delucia)
