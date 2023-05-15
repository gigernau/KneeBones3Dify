# KneeBones3Dify
<img src="https://github.com/gigernau/KneeBones3Dify/blob/main/KneeBones3Dify_logo.png" align="center" height="250" width="250">


# Setup for Ubuntu Linux OS

## 1) Update : 
	sudo apt update
  
## 2) Install dependencies:
    sudo apt install build-essential zlib1g-dev libncurses5-dev libgdbm-dev libnss3-dev libssl-dev libreadline-dev libffi-dev python-tk python3-tk tk-dev python3-pil python3-pil.imagetk

## 3) Install pip3 for Python3: 
	sudo apt install python3-pip  && python3 -m pip install --upgrade pip

## 4) Install Python3 modules: 
	python3 -m pip install -r requirements.txt
	
## 5) Check Cuda version or Install [CUDA](https://developer.nvidia.com/cuda-toolkit)
	nvcc -V

## 6) Install CuPy module ( e.g., for CUDA 11.1 )
	python3 -m pip install cupy-cuda11x

## 7) Compile smoothPatch code:
### Linux
	gcc -shared -o smoothPatch.so smoothPatch.cpp
### Windows ( install https://winlibs.com/ )
	gcc -c smoothPatch.cpp
	gcc -shared -o smoothPatch.dll smoothPatch.o


# USAGE EXAMPLE
	python3 SegOscanSST.py


[![Say Thanks!](https://img.shields.io/badge/Say%20Thanks-!-1EAEDB.svg)](https://saythanks.io/to/gianluca.delucia)
