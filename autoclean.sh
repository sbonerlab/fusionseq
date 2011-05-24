#!/bin/bash

echo autoclean: cleaning root
rm -rf aclocal.m4 \
       autom4te.cache \
       autoscan.log \
       config.guess \
       config.h \
       config.h.in~ \
       config.log \
       config.status \
       config.sub \
       configure \
       depcomp \
       install-sh \
       libtool \
       ltmain.sh \
       missing \
       Makefile.in \
       Makefile \
       stamp-h1

echo autoclean: cleaning m4
cd m4 && rm -rf \
       libtool.m4 \
       lt~obsolete.m4 \
       ltoptions.m4 \
       ltsugar.m4 \
       ltversion.m4
cd ..

echo autoclean: cleaning src
cd src && rm -rf \
       Makefile.in \
       Makefile
echo autoclean: cleaning src/optional
cd optional && rm -rf \
       Makefile.in \
       Makefile
cd ..
echo autoclean: cleaning src/test
cd test && rm -rf \
       Makefile.in \
       Makefile
cd ..
echo autoclean: cleaning src/cgi
cd cgi && rm -rf \
       Makefile.in \
       Makefile
cd ..
cd ..
