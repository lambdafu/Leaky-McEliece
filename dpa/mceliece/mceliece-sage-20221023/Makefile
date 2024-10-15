SAGELIBS=\
byterepr parameters \
echelon matgen irreducible fieldordering keygen \
encode fixedweight encap \
interpolator approximant goppa decode decap \
test-padding test-checksums

default: runtests

SAGEPY=$(SAGELIBS:=.py)
SAGEOUT=$(SAGELIBS:=.out)

runtests: controlbits.out $(SAGEOUT)

controlbits.out: controlbits.py
	python3 $< > $@

$(SAGEOUT): %.out: %.py sagelibs
	sage $< > $@

sagelibs: $(SAGEPY)

%.sage.py: %.sage
	sage -preparse $<

$(SAGEPY): %.py: %.sage.py
	cp $< $@
