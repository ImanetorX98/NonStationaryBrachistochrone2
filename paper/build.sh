#!/usr/bin/env bash
# Builds paper3.pdf and supplement_fronts.pdf. The two documents cross-reference each other (package xr-hyper),
# so each is compiled after the other has written its .aux file. Requires a TeX distribution with revtex4-2.
set -e
cd "$(dirname "$0")"
pdflatex -interaction=nonstopmode paper3.tex > /dev/null
pdflatex -interaction=nonstopmode supplement_fronts.tex > /dev/null
bibtex paper3 > /dev/null
for i in 1 2; do
  pdflatex -interaction=nonstopmode paper3.tex > /dev/null
  pdflatex -interaction=nonstopmode supplement_fronts.tex > /dev/null
done
grep -H -E "^!|(Reference|Citation).*undefined|undefined references" paper3.log supplement_fronts.log && { echo "errors or undefined references, see the .log files"; exit 1; }
echo "built paper3.pdf and supplement_fronts.pdf"
