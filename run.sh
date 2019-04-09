#!/bin/bash


# Script that returns a plot
echo 'running main script '
python3 main.py

echo 'generating the pdf'

pdflatex handin_1.tex
