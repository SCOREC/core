TARGET=parmaGhost
OUTPUT_DIR=build

# make pdf by default
all: ${TARGET}.pdf

# it doesn't really need the .dvi, but this way all the refs are right
%.pdf : %.dvi
	pdflatex -output-directory=${OUTPUT_DIR} $* 

${TARGET}.bbl: references.bib
# in case we don't already have a .aux file listing citations
# this should probably be a separate makefile target/dependency instead
# of doing it every time... but *shrug*
	pdflatex -output-directory=${OUTPUT_DIR} ${TARGET}.tex 
# get the citations out of the bibliography
	bibtex ${OUTPUT_DIR}/${TARGET} 
# do it again in case there are out-of-order cross-references
	@pdflatex -output-directory=${OUTPUT_DIR} ${TARGET}.tex 

${TARGET}.dvi: ${TARGET}.bbl ${TARGET}.tex
	@pdflatex -output-directory=${OUTPUT_DIR} ${TARGET}.tex 

# shortcut, so we can say "make ps"
ps: ${TARGET}.ps

${TARGET}.ps: ${TARGET}.dvi
	@dvips -t a4 ${TARGET}.dvi

clean:
	rm -f ${OUTPUT_DIR}/${TARGET}.a* ${OUTPUT_DIR}/${TARGET}.b* ${OUTPUT_DIR}/${TARGET}.log

reallyclean: clean
	rm -f ${OUTPUT_DIR}/${TARGET}.pdf

.PHONY : ps all clean reallyclean
