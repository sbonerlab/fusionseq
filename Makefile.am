SUBDIRS = src
ACLOCAL_AMFLAGS = -I m4


optional:
	cd src && $(MAKE) $(AM_MAKEFLAGS) optional

cgi:
	cd src && $(MAKE) $(AM_MAKEFLAGS) cgi

test:
	cd src && $(MAKE) $(AM_MAKEFLAGS) test

.PHONY: optional cgi test
