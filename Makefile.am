EXTRA_DIST = m4 scripts default.fusionseqrc
ACLOCAL_AMFLAGS = -I m4
AM_MAKEFLAGS = --no-print-directory

AM_CPPFLAGS = \
	-I${top_srcdir}/src

PC_SED = \
	$(AM_V_GEN)$(MKDIR_P) $(dir $@) && $(SED) \
	-e 's,@VERSION\@,$(VERSION),g' \
	-e 's,@prefix\@,$(prefix),g' \
	-e 's,@exec_prefix\@,$(exec_prefix),g' \
	-e 's,@libdir\@,$(libdir),g' \
	-e 's,@includedir\@,$(includedir),g' \
	< $< > $@ || rm $@

%.pc: %.pc.in Makefile
	$(PC_SED)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = fusionseq.pc
EXTRA_DIST += fusionseq.pc

# -----------------------------------------------------------------------------
# Common library
# -----------------------------------------------------------------------------

lib_LTLIBRARIES = libfusionseq.la
libfusionseq_la_SOURCES = \
	src/bp.c \
	src/gfr.c \
	src/sqvCircos.c \
	src/sqvUtil.c \
	src/sqvWeb.c \
	src/util.c

libfusionseq_la_LIBADD = -lm -lbios 
libfusionseq_la_LDFLAGS = -version-info $(FUSIONSEQ_SO_VERSION)

include_HEADERS = \
	src/bp.h \
	src/gfr.h \
	src/sqvCircos.h \
	src/sqvUtil.h \
	src/sqvWeb.h \
	src/util.h \
	src/uthash.h 

# -----------------------------------------------------------------------------
# FusionSeq programs
# -----------------------------------------------------------------------------

bin_PROGRAMS = \
	geneFusions \
	gfrAnnotationConsistencyFilter \
	gfrLargeScaleHomologyFilter \
	gfrSmallScaleHomologyFilter \
	gfrAbnormalInsertSizeFilter \
	gfrAddInfo \
	gfrProximityFilter \
	bpFilter \
	bp2alignment \
	gfrPCRFilter \
	gfrBlackListFilter \
	gfr2bpJunctions \
	bowtie2bp \
	sam2bp \
	gfrConfidenceValues \
	validateBpJunctions \
	gfrExpressionConsistencyFilter \
	gfrMitochondrialFilter \
	gfrRepeatMaskerFilter \
	gfrCountPairTypes \
	gfrRibosomalFilter \
	gfrSpliceJunctionFilter \
	gfrClassify \
	gfrWhiteListFilter \
	gfrRandomPairingFilter

if BUILD_CGI

EXTRA_PROGRAMS = \
	geneFusions_cgi \
	showDetails_cgi \
	seqViz_cgi \
	findFusionPartner.cgi

endif

if BUILD_OPTIONAL

bin_PROGRAMS +=	\
	gfr2bed \
	gfr2fasta \
	gfr2gff \
	bp2wig \
	export2mrf \
	gfrConfidenceValueTranscript

endif

if BUILD_TESTS

bin_PROGRAMS += \
	gfSimulation \
	bowtiePairedFix \
	mrfPseudogeneOverlap \
	plotIntraDistribution \
	quantifierAddInfo \
	bp2fasta \
	bpClustering \
	bpValidate

endif

# -----------------------------------------------------------------------------
# CORE: Finding candidate fusions
# -----------------------------------------------------------------------------

geneFusions_SOURCES = src/geneFusions.c
geneFusions_LDADD = libfusionseq.la -lbios -lmrf -lm 

gfrExpressionConsistencyFilter_SOURCES = src/gfrExpressionConsistencyFilter.c
gfrExpressionConsistencyFilter_LDADD = libfusionseq.la -lbios -lm 

gfrAnnotationConsistencyFilter_SOURCES = src/gfrAnnotationConsistencyFilter.c
gfrAnnotationConsistencyFilter_LDADD = libfusionseq.la -lbios -lm

gfrRibosomalFilter_SOURCES = src/gfrRibosomalFilter.c
gfrRibosomalFilter_LDADD = libfusionseq.la -lbios

gfrMitochondrialFilter_SOURCES = src/gfrMitochondrialFilter.c
gfrMitochondrialFilter_LDADD = libfusionseq.la -lbios

gfrLargeScaleHomologyFilter_SOURCES = src/gfrLargeScaleHomologyFilter.c
gfrLargeScaleHomologyFilter_LDADD = libfusionseq.la -lbios

gfrSmallScaleHomologyFilter_SOURCES = src/gfrSmallScaleHomologyFilter.c
gfrSmallScaleHomologyFilter_LDADD = libfusionseq.la -lbios

gfrAbnormalInsertSizeFilter_SOURCES = src/gfrAbnormalInsertSizeFilter.c
gfrAbnormalInsertSizeFilter_LDADD = libfusionseq.la -lbios

gfrPCRFilter_SOURCES = src/gfrPCRFilter.c
gfrPCRFilter_LDADD = libfusionseq.la -lbios

gfrBlackListFilter_SOURCES = src/gfrBlackListFilter.c
gfrBlackListFilter_LDADD = libfusionseq.la -lbios

gfrProximityFilter_SOURCES = src/gfrProximityFilter.c
gfrProximityFilter_LDADD = libfusionseq.la -lbios

gfrAddInfo_SOURCES = src/gfrAddInfo.c
gfrAddInfo_LDADD = libfusionseq.la -lbios

gfrConfidenceValues_SOURCES = src/gfrConfidenceValues.c
gfrConfidenceValues_LDADD = libfusionseq.la -lbios -lmrf

gfrRepeatMaskerFilter_SOURCES = src/gfrRepeatMaskerFilter.c
gfrRepeatMaskerFilter_LDADD = libfusionseq.la -lbios -lm

gfrCountPairTypes_SOURCES = src/gfrCountPairTypes.c
gfrCountPairTypes_LDADD = libfusionseq.la -lbios -lm

gfrSpliceJunctionFilter_SOURCES = src/gfrSpliceJunctionFilter.c
gfrSpliceJunctionFilter_LDADD = libfusionseq.la -lbios

gfrClassify_SOURCES = src/gfrClassify.c
gfrClassify_LDADD = libfusionseq.la -lbios -lm

gfrWhiteListFilter_SOURCES = src/gfrWhiteListFilter.c
gfrWhiteListFilter_LDADD = libfusionseq.la -lbios

gfrRandomPairingFilter_SOURCES = src/gfrRandomPairingFilter.c
gfrRandomPairingFilter_LDADD = libfusionseq.la -lbios

# -----------------------------------------------------------------------------
# CORE: Identifying sequences of the junction
# -----------------------------------------------------------------------------

gfr2bpJunctions_SOURCES = src/gfr2bpJunctions.c
gfr2bpJunctions_LDADD = libfusionseq.la -lbios -lm

bowtie2bp_SOURCES = src/bowtie2bp.c
bowtie2bp_LDADD = libfusionseq.la -lbios

sam2bp_SOURCES = src/sam2bp.c
sam2bp_LDADD = libfusionseq.la -lbios

validateBpJunctions_SOURCES = src/validateBpJunctions.c
validateBpJunctions_LDADD = libfusionseq.la -lbios

bpFilter_SOURCES = src/bpFilter.cpp
bpFilter_CXXFLAGS = -D_REENTRANT -pthread -m64
bpFilter_LDADD = libfusionseq.la -lstdc++ -lbios -L@ROOTLIBDIR@ @ROOTGLIBS@ @ROOTLIBS@ -ldl

bp2alignment_SOURCES = src/bp2alignment.c
bp2alignment_LDADD = libfusionseq.la -lbios


# -----------------------------------------------------------------------------
# CGI
# -----------------------------------------------------------------------------

if BUILD_CGI

geneFusions_cgi_SOURCES = src/cgi/geneFusions_cgi.c
geneFusions_cgi_LDADD = libfusionseq.la -lbios

showDetails_cgi_SOURCES = src/cgi/showDetails_cgi.c
showDetails_cgi_LDADD = libfusionseq.la -lbios

seqViz_cgi_SOURCES = src/cgi/seqViz_cgi.c
seqViz_cgi_LDADD = libfusionseq.la -lbios -lmrf

findFusionPartner_cgi_SOURCES = src/cgi/findFusionPartner_cgi.c
findFusionPartner_cgi_LDADD = libfusionseq.la -lbios

endif

# -----------------------------------------------------------------------------
# OPTIONAL
# -----------------------------------------------------------------------------

if BUILD_OPTIONAL

# gfr2images_SOURCES = src/optional/gfr2images.c
# gfr2images_LDADD = libfusionseq.la -lbios -lgd -lm

gfr2bed_SOURCES = src/optional/gfr2bed.c
gfr2bed_LDADD = libfusionseq.la -lbios

gfr2fasta_SOURCES = src/optional/gfr2fasta.c
gfr2fasta_LDADD = libfusionseq.la -lbios

gfr2gff_SOURCES = src/optional/gfr2gff.c
gfr2gff_LDADD = libfusionseq.la -lbios

bp2wig_SOURCES = src/optional/bp2wig.c
bp2wig_LDADD = libfusionseq.la -lbios

export2mrf_SOURCES = src/optional/export2mrf.c
export2mrf_LDADD = libfusionseq.la -lbios

gfrConfidenceValueTranscript_SOURCES = src/optional/gfrConfidenceValueTranscript.c
gfrConfidenceValueTranscript_LDADD = libfusionseq.la -lbios

endif

# -----------------------------------------------------------------------------
# TEST
# -----------------------------------------------------------------------------

if BUILD_TESTS

gfSimulation_SOURCES = src/test/gfSimulation.c
gfSimulation_LDADD = libfusionseq.la -lbios -lm

bowtiePairedFix_SOURCES = src/test/bowtiePairedFix.c
bowtiePairedFix_LDADD = libfusionseq.la -lbios -lmrf -lm

mrfPseudogeneOverlap_SOURCES = src/test/mrfPseudogeneOverlap.c
mrfPseudogeneOverlap_LDADD = libfusionseq.la -lbios -lmrf -lm

quantifierAddInfo_SOURCES = src/test/quantifierAddInfo.c
quantifierAddInfo_LDADD = libfusionseq.la -lbios

plotIntraDistribution_SOURCES = src/test/plotIntraDistribution.cpp
plotIntraDistribution_CXXFLAGS = -D_REENTRANT -pthread -m64
plotIntraDistribution_LDADD = libfusionseq.la -lstdc++ -lbios -L@ROOTLIBDIR@ @ROOTGLIBS@ @ROOTLIBS@ -ldl

bp2fasta_SOURCES = src/test/bp2fasta.c
bp2fasta_LDADD = libfusionseq.la -lbios -lm

bpClustering_SOURCES = src/test/bpClustering.c
bpClustering_LDADD = libfusionseq.la -lbios -lm

bpValidate_SOURCES = src/test/bpValidate.c
bpValidate_LDADD = libfusionseq.la -lbios -lm

endif

#------------------------------------------------------------------------------
# Misc
#------------------------------------------------------------------------------

install-data-hook:
	scripts/install_rc.sh

deploy:
	scripts/install_deploy.pl

if USE_DOXYGEN

doxygen:
	cd doc && doxygen Doxyfile.in

endif
