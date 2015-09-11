#!/bin/bash

# Script to make latex files in doc folder
# Only works with tex files and png images

if [ -e "CMakeLists.txt" ]
then
	echo "Do not run this in the main directory."
	echo "Run in the build directory, please."
	return;
fi

if [ ! -d "doc" ]; 
then
	mkdir doc
fi
cd doc

# These files are needed, quick check if they are available, otherwise
if [ -z $(kpsewhich algorithmicx.sty) ];
then
	wget http://mirror.math.ku.edu/tex-archive/macros/latex/contrib/algorithmicx/algorithmicx.sty
fi
if [ -z $(kpsewhich algpseudocode.sty) ];
then
	wget http://mirror.unl.edu/ctan/macros/latex/contrib/algorithmicx/algpseudocode.sty
fi
for f in ../../; do
	cp $f*.tex .
	cp $f*.png .
done
for f in ../../*/; do
	cp $f*.tex .
	cp $f*.png .
done

# will do what it can
for tex in *.tex; do
	pdflatex -interaction=nonstopmode $tex
done 
rm *.aux *.log *.png *.tex
cd ..

