PFILES	= research

FILES	= research

DEPS = research.tex intro.tex amg.tex emin.tex work.tex
          
#USER_TEXINPUTS:=$(TEXINPUTS)
#export TEXINPUTS=$(USER_TEXINPUTS):.:./sty:

all:	pdf

pdf:	$(addsuffix .pdf, $(FILES) $(PFILES))

ps:	$(addsuffix .ps, $(FILES))

# Build one example using pdflatex
research.pdf:	$(DEPS)
	@ echo "TEXINPUTS = $(TEXINPUTS)"
	pdflatex $<
	#bibtex $(basename $<)
	#pdflatex $<
	pdflatex $<

%.pdf:	%.ps
	ps2pdf13 $< $@

%.ps:	%.dvi
	dvips -Ppdf -o $@ $<

%.dvi:	%.tex parallel.bib SANDreport.cls
	latex $<
	bibtex $(basename $<)
	latex $<
	latex $<

clean:
	@rm -f $(addsuffix .aux, $(FILES) $(PFILES)) $(addsuffix .bbl, $(FILES) $(PFILES))
	@rm -f $(addsuffix .blg, $(FILES) $(PFILES)) $(addsuffix .lof, $(FILES) $(PFILES))
	@rm -f $(addsuffix .log, $(FILES) $(PFILES)) $(addsuffix .lot, $(FILES) $(PFILES))
	@rm -f $(addsuffix .toc, $(FILES) $(PFILES)) $(addsuffix .out, $(FILES) $(PFILES))
	@rm -f texput.log
	@rm -f Mark*.aux
	@rm -f distribution.aux

realclean:	clean
	@rm -f $(addsuffix .pdf, $(FILES) $(PFILES)) $(addsuffix .ps, $(FILES) $(PFILES))
	@rm -f $(addsuffix .dvi, $(FILES) $(PFILES)) $(addsuffix .out, $(FILES) $(PFILES))
