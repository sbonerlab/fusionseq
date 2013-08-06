SUBDIRS = src
ACLOCAL_AMFLAGS = -I m4

DIST_SUBDIRS = src src/optional src/cgi src/test
EXTRA_DIST = scripts default.fusionseqrc

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = fusionseq.pc

optional:
	cd src/optional && $(MAKE) $(AM_MAKEFLAGS)

cgi:
	cd src/cgi && $(MAKE) $(AM_MAKEFLAGS) cgi

test:
	cd src/test && $(MAKE) "CFLAGS=-g -DDEBUG " $(AM_MAKEFLAGS)

install-exec-local:
	if test -e src/optional/bp2wig; then cd src/optional && $(MAKE) install; fi
	if test -e src/test/bp2fasta; then cd src/test && $(MAKE) install; fi

install-data-hook:
	scripts/install_rc.sh

deploy:
	scripts/install_deploy.pl

debug:
	$(MAKE) all "CFLAGS=-g -DDEBUG " $(AM_MAKEFLAGS)
	$(MAKE) optional "CFLAGS=-g -DDEBUG " $(AM_MAKEFLAGS)

if USE_DOXYGEN

doxygen:
	cd doc && doxygen Doxyfile.in

clean-local:
	cd src/optional && $(MAKE) clean
	cd src/test && $(MAKE) clean
	cd src/cgi && $(MAKE) clean
	rm -rf $(top_srcdir)/doc/html

.PHONY: doxygen clean-local install-exec-local install-data-hook deploy

else

clean-local:
	cd src/optional && $(MAKE) clean
	cd src/test && $(MAKE) clean
	cd src/cgi && $(MAKE) clean

.PHONY: clean-local install-exec-local install-data-hook deploy

endif
