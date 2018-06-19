#!/bin/bash

mkdir 273
mkdir 298
mkdir 345

python prepare-for-fah-solvate.py

python prepare-for-fah-equilibrate.py 273
python prepare-for-fah-equilibrate.py 298
python prepare-for-fah-equilibrate.py 345
