#!/bin/bash
cp prioritize2.py prioritize.py

python freeze.py build

rm prioritize.py

cp prioritize3.py prioritize.py

python3 freeze.py build

rm prioritize.py
