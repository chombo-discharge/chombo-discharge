# -*- Mode: Makefile -*- 

## Rolls all the ccse stuff we need into a single library.
## Rolls ccse/facade into another library.

lib: facade_lib

facade_lib: ccse_all_lib
	(cd ../CCSE_all/facade; $(MAKE) $(MAKEFLAGS) lib)

ccse_all_lib: ../CCSE_all
	(cd ../CCSE_all; $(MAKE) $(MAKEFLAGS) lib)

../CCSE_all: .
	./make_ccse_all.sh

clean: ccse_all_clean facade_clean

ccse_all_clean:
	(cd ../CCSE_all; $(MAKE) $(MAKEFLAGS) clean)

facade_clean:
	(cd ../CCSE_all/facade; $(MAKE) $(MAKEFLAGS) clean)

realclean:
	rm -rf ../CCSE_all
