# makeenv optional file for additional setting
# LATEXARGS additional flags to latex to be set in makeenv
# LATEXMK latex compilation defaulting to latexmk, can be changed in makeenv
# Example: insert LATEXARGS=-file-line-error into makeenv

-include makeenv
TRGT=main

LATEKMK?=latexmk


all:
	@$(LATEKMK) -pdf $(LATEXARGS) $(TRGT) && cp $(TRGT).pdf ~/win/Documents/

clean:
	@latexmk -c $(TRGT)

allclean: clean
	@latexmk -C $(TRGT)
	@rm -fv $(TRGT).bbl

b:
	bibtex qbf_proofs
