## This is the Makefile for the template paper. All you need to do is 
## edit "TARGETS" below to the name of the main .tex file for your paper.

TARGETS = draft

TEXFILES = $(wildcard *.tex)
FIGURES =$(wildcard figures/*.pdf figures/*.tex)
GENPDFS = $(patsubst %.svg,%.pdf,$(wildcard figures/*.svg))
PDFS = $(addsuffix .pdf,$(TARGETS))
BIBFILES=draft.bib


all: once

%.pdf: %.tex %.blg %.toc $(TEXFILES) $(FIGURES) $(BIBFILES)
	pdflatex $*.tex
	bibtex $*
	pdflatex $*.tex
	pdflatex $*.tex

%.pdf: %.svg
	inkscape --export-pdf=$@ $<

once:
	pdflatex -shell-escape draft.tex
	bibtex draft
	pdflatex -shell-escape draft.tex
	pdflatex -shell-escape draft.tex

%.blg: $(GENPDFS) $(BIBFILES)
	pdflatex $*.tex
	bibtex $*
	pdflatex $*.tex

%.toc: %.tex $(GENPDFS)
	pdflatex $*.tex

clean:
	/bin/rm -f $(PDFS) *.dvi *.aux *.out *.ps *~ *.log *.lot *.lof *.toc *.blg *.bbl url.sty $(GENPDFS)

FORCE:

