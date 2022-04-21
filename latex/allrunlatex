#!/bin/sh

MAIN='main'

pdflatex --shell-escape ${MAIN}  
#biber     ${MAIN}  
bibtex   ${MAIN}
pdflatex --shell-escape ${MAIN}  
pdflatex --shell-escape ${MAIN} 

exit 0 
