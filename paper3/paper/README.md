# Manuscript and supplementary material

LaTeX sources of the paper (`paper3.tex`) and of its supplementary material (`supplement_fronts.tex`), with the
bibliography `refs3.bib`, and the two PDF files built from them (`paper3.pdf`, `supplement_fronts.pdf`). The figures are read from `../figures/`, where the scripts that draw them are.

```
./build.sh      # pdflatex and bibtex in the right order; requires revtex4-2
```

The two documents cross-reference each other with the package `xr-hyper` (labels of the supplement appear in the
paper with the prefix `S-`, and labels of the paper in the supplement with the prefix `M-`), so each is compiled
after the other has written its `.aux` file; `build.sh` does this. The map from the statements of the paper to the
scripts of this repository is in [`../RESULTS_MAP.md`](../RESULTS_MAP.md).
