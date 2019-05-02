#!/bin/bash

# Package up paper files for submission to ApJ, arXiv, etc. This 
# compiles all the TeX source files into a single file, removes the
# dependency on PythonTeX, generates all of the figures, and fixes
# the figure paths.
#
# Usage:
#   ./package_for_aas.sh /path/to/top/level/dir
#   ./package_for_aas.sh /path/to/top/level/dir --cover-letter

# Create target directory
DIR_NAME="aas-submission"
DEST_DIR="$1/$DIR_NAME"
mkdir -p $DEST_DIR

# Check for cover letter option
if [[ $* == *--cover-letter* ]]
then
    COVER_LETTER=true
else
    COVER_LETTER=false
fi

# Copy over all paper files
cp -r paper/sections $DEST_DIR
cp -r paper/data $DEST_DIR
cp -r paper/python $DEST_DIR
cp paper/aasjournal.bst $DEST_DIR
cp paper/aastex62.cls $DEST_DIR
cp paper/pythontex.tex $DEST_DIR
cp paper/paper.tex $DEST_DIR
cp paper/references.bib $DEST_DIR
cp paper/software.bib $DEST_DIR
cp fix_figure_paths.py $DEST_DIR
if $COVER_LETTER ; then cp cover-letter.md $DEST_DIR ; fi

# Flatten to a single file
cd $DEST_DIR
mv paper.tex _paper.tex
latexpand _paper.tex --empty-comments --keep-comments -o paper.tex
rm _paper.tex

# Run PythonTeX
pdflatex -synctex=1 -interaction=nonstopmode -file-line-error paper.tex
pythontex --interpreter python:python paper.tex
bibtex paper
pdflatex -synctex=1 -interaction=nonstopmode -file-line-error paper.tex
pdflatex -synctex=1 -interaction=nonstopmode -file-line-error paper.tex

# Remove Pythontex dependence and fix figure paths
depythontex paper.tex --overwrite -o paper.tex
python fix_figure_paths.py -i paper.tex -o paper.tex
cp figures/* .
rm -r figures

# Build PDF again without pythontex
pdflatex -synctex=1 -interaction=nonstopmode -file-line-error paper.tex
bibtex paper
pdflatex -synctex=1 -interaction=nonstopmode -file-line-error paper.tex
pdflatex -synctex=1 -interaction=nonstopmode -file-line-error paper.tex

# (optionally) Generate PDF of cover letter
if $COVER_LETTER
then
    pandoc -s cover-letter.md -o cover-letter.pdf --bibliography="references.bib"
fi

# Cleanup
rm -r sections
rm -r data
rm -r pythontex-files-paper
rm -r python
rm *.aux
rm *.bbl
rm *.blg
rm *.depytx
rm *.log
rm *.out
rm *.pytxcode
rm *.synctex.gz
rm *.py
rm pythontex.tex
if $COVER_LETTER ; then rm cover-letter.md ; fi

# Compress directory and cleanup
cd $1
tar cvzf "$DEST_DIR.tar.gz" $DIR_NAME
rm -r $DIR_NAME
