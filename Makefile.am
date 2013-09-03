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
# Convenience library
# -----------------------------------------------------------------------------

noinst_LTLIBRARIES = src/libfusionseq.la
src_libfusionseq_la_SOURCES = \
	src/bp.c \
	src/gfr.c \
	src/util.c

# -----------------------------------------------------------------------------
# FusionSeq programs
# -----------------------------------------------------------------------------

bin_PROGRAMS = \
	src/geneFusions \
	src/gfrAnnotationConsistencyFilter \
	src/gfrLargeScaleHomologyFilter \
	src/gfrSmallScaleHomologyFilter \
	src/gfrAbnormalInsertSizeFilter \
	src/gfrAddInfo \
	src/gfrProximityFilter \
	src/bpFilter \
	src/bp2alignment \
	src/gfrPCRFilter \
	src/gfrBlackListFilter \
	src/gfr2bpJunctions \
	src/bowtie2bp \
	src/sam2bp \
	src/gfrConfidenceValues \
	src/validateBpJunctions \
	src/gfrExpressionConsistencyFilter \
	src/gfrMitochondrialFilter \
	src/gfrRepeatMaskerFilter \
	src/gfrCountPairTypes \
	src/gfrRibosomalFilter \
	src/gfrSpliceJunctionFilter \
	src/gfrClassify \
	src/gfrWhiteListFilter \
	src/gfrRandomPairingFilter

if BUILD_CGI

noinst_PROGRAMS = \
	src/cgi/geneFusions_cgi \
	src/cgi/showDetails_cgi \
	src/cgi/seqViz_cgi \
	src/cgi/findFusionPartner.cgi

endif

if BUILD_OPTIONAL

bin_PROGRAMS +=	\
	src/optional/gfr2bed \
	src/optional/gfr2fasta \
	src/optional/gfr2gff \
	src/optional/bp2wig \
	src/optional/export2mrf \
	src/optional/gfrConfidenceValueTranscript

endif

if BUILD_TESTS

bin_PROGRAMS += \
	src/test/gfSimulation \
	src/test/bowtiePairedFix \
	src/test/mrfPseudogeneOverlap \
	src/test/plotIntraDistribution \
	src/test/quantifierAddInfo \
	src/test/bp2fasta \
	src/test/bpClustering \
	src/test/bpValidate

endif

# -----------------------------------------------------------------------------
# CORE: Finding candidate fusions
# -----------------------------------------------------------------------------

src_geneFusions_SOURCES = src/geneFusions.c
src_geneFusions_LDADD = src/libfusionseq.la -lbios -lmrf -lm 

src_gfrExpressionConsistencyFilter_SOURCES = src/gfrExpressionConsistencyFilter.c
src_gfrExpressionConsistencyFilter_LDADD = src/libfusionseq.la -lbios -lm 

src_gfrAnnotationConsistencyFilter_SOURCES = src/gfrAnnotationConsistencyFilter.c
src_gfrAnnotationConsistencyFilter_LDADD = src/libfusionseq.la -lbios -lm

src_gfrRibosomalFilter_SOURCES = src/gfrRibosomalFilter.c
src_gfrRibosomalFilter_LDADD = src/libfusionseq.la -lbios

src_gfrMitochondrialFilter_SOURCES = src/gfrMitochondrialFilter.c
src_gfrMitochondrialFilter_LDADD = src/libfusionseq.la -lbios

src_gfrLargeScaleHomologyFilter_SOURCES = src/gfrLargeScaleHomologyFilter.c
src_gfrLargeScaleHomologyFilter_LDADD = src/libfusionseq.la -lbios

src_gfrSmallScaleHomologyFilter_SOURCES = src/gfrSmallScaleHomologyFilter.c
src_gfrSmallScaleHomologyFilter_LDADD = src/libfusionseq.la -lbios

src_gfrAbnormalInsertSizeFilter_SOURCES = src/gfrAbnormalInsertSizeFilter.c
src_gfrAbnormalInsertSizeFilter_LDADD = src/libfusionseq.la -lbios

src_gfrPCRFilter_SOURCES = src/gfrPCRFilter.c
src_gfrPCRFilter_LDADD = src/libfusionseq.la -lbios

src_gfrBlackListFilter_SOURCES = src/gfrBlackListFilter.c
src_gfrBlackListFilter_LDADD = src/libfusionseq.la -lbios

src_gfrProximityFilter_SOURCES = src/gfrProximityFilter.c
src_gfrProximityFilter_LDADD = src/libfusionseq.la -lbios

src_gfrAddInfo_SOURCES = src/gfrAddInfo.c
src_gfrAddInfo_LDADD = src/libfusionseq.la -lbios

src_gfrConfidenceValues_SOURCES = src/gfrConfidenceValues.c
src_gfrConfidenceValues_LDADD = src/libfusionseq.la -lbios -lmrf

src_gfrRepeatMaskerFilter_SOURCES = src/gfrRepeatMaskerFilter.c
src_gfrRepeatMaskerFilter_LDADD = src/libfusionseq.la -lbios -lm

src_gfrCountPairTypes_SOURCES = src/gfrCountPairTypes.c
src_gfrCountPairTypes_LDADD = src/libfusionseq.la -lbios -lm

src_gfrSpliceJunctionFilter_SOURCES = src/gfrSpliceJunctionFilter.c
src_gfrSpliceJunctionFilter_LDADD = src/libfusionseq.la -lbios

src_gfrClassify_SOURCES = src/gfrClassify.c
src_gfrClassify_LDADD = src/libfusionseq.la -lbios -lm

src_gfrWhiteListFilter_SOURCES = src/gfrWhiteListFilter.c
src_gfrWhiteListFilter_LDADD = src/libfusionseq.la -lbios

src_gfrRandomPairingFilter_SOURCES = src/gfrRandomPairingFilter.c
src_gfrRandomPairingFilter_LDADD = src/libfusionseq.la -lbios

# -----------------------------------------------------------------------------
# CORE: Identifying sequences of the junction
# -----------------------------------------------------------------------------

src_gfr2bpJunctions_SOURCES = src/gfr2bpJunctions.c
src_gfr2bpJunctions_LDADD = src/libfusionseq.la -lbios -lm

src_bowtie2bp_SOURCES = src/bowtie2bp.c
src_bowtie2bp_LDADD = src/libfusionseq.la -lbios

src_sam2bp_SOURCES = src/sam2bp.c
src_sam2bp_LDADD = src/libfusionseq.la -lbios

src_validateBpJunctions_SOURCES = src/validateBpJunctions.c
src_validateBpJunctions_LDADD = src/libfusionseq.la -lbios

src_bpFilter_SOURCES = src/bpFilter.cpp
src_bpFilter_CXXFLAGS = -D_REENTRANT -pthread -m64
src_bpFilter_LDADD = src/libfusionseq.la -lstdc++ -lbios -L@ROOTLIBDIR@ @ROOTGLIBS@ @ROOTLIBS@ -ldl

src_bp2alignment_SOURCES = src/bp2alignment.c
src_bp2alignment_LDADD = src/libfusionseq.la -lbios


# -----------------------------------------------------------------------------
# CGI
# -----------------------------------------------------------------------------

if BUILD_CGI

src_cgi_geneFusions_cgi_SOURCES = src/cgi/geneFusions_cgi.c
src_cgi_geneFusions_cgi_LDADD = src/libfusionseq.la -lbios

src_cgi_showDetails_cgi_SOURCES = src/cgi/showDetails_cgi.c
src_cgi_showDetails_cgi_LDADD = src/libfusionseq.la -lbios

src_cgi_seqViz_cgi_SOURCES = \
	src/cgi/sqvUtil.c \
	src/cgi/sqvWeb.c \
	src/cgi/sqvCircos.c \
	src/cgi/seqViz_cgi.c
src_cgi_seqViz_cgi_LDADD = src/libfusionseq.la -lbios -lmrf

src_cgi_findFusionPartner_cgi_SOURCES = src/cgi/findFusionPartner_cgi.c
src_cgi_findFusionPartner_cgi_LDADD = src/libfusionseq.la -lbios

endif

# -----------------------------------------------------------------------------
# OPTIONAL
# -----------------------------------------------------------------------------

if BUILD_OPTIONAL

# src_optional_gfr2images_SOURCES = src/optional/gfr2images.c
# src_optional_gfr2images_LDADD = src/libfusionseq.la -lbios -lgd -lm

src_optional_gfr2bed_SOURCES = src/optional/gfr2bed.c
src_optional_gfr2bed_LDADD = src/libfusionseq.la -lbios

src_optional_gfr2fasta_SOURCES = src/optional/gfr2fasta.c
src_optional_gfr2fasta_LDADD = src/libfusionseq.la -lbios

src_optional_gfr2gff_SOURCES = src/optional/gfr2gff.c
src_optional_gfr2gff_LDADD = src/libfusionseq.la -lbios

src_optional_bp2wig_SOURCES = src/optional/bp2wig.c
src_optional_bp2wig_LDADD = src/libfusionseq.la -lbios

src_optional_export2mrf_SOURCES = src/optional/export2mrf.c
src_optional_export2mrf_LDADD = src/libfusionseq.la -lbios

src_optional_gfrConfidenceValueTranscript_SOURCES = src/optional/gfrConfidenceValueTranscript.c
src_optional_gfrConfidenceValueTranscript_LDADD = src/libfusionseq.la -lbios

endif

# -----------------------------------------------------------------------------
# TEST
# -----------------------------------------------------------------------------

if BUILD_TESTS

src_test_gfSimulation_SOURCES = src/test/gfSimulation.c
src_test_gfSimulation_LDADD = src/libfusionseq.la -lbios -lm

src_test_bowtiePairedFix_SOURCES = src/test/bowtiePairedFix.c
src_test_bowtiePairedFix_LDADD = src/libfusionseq.la -lbios -lmrf -lm

src_test_mrfPseudogeneOverlap_SOURCES = src/test/mrfPseudogeneOverlap.c
src_test_mrfPseudogeneOverlap_LDADD = src/libfusionseq.la -lbios -lmrf -lm

src_test_quantifierAddInfo_SOURCES = src/test/quantifierAddInfo.c
src_test_quantifierAddInfo_LDADD = src/libfusionseq.la -lbios

src_test_plotIntraDistribution_SOURCES = src/test/plotIntraDistribution.cpp
src_test_plotIntraDistribution_CXXFLAGS = -D_REENTRANT -pthread -m64
src_test_plotIntraDistribution_LDADD = src/libfusionseq.la -lstdc++ -lbios -L@ROOTLIBDIR@ @ROOTGLIBS@ @ROOTLIBS@ -ldl

src_test_bp2fasta_SOURCES = src/test/bp2fasta.c
src_test_bp2fasta_LDADD = src/libfusionseq.la -lbios -lm

src_test_bpClustering_SOURCES = src/test/bpClustering.c
src_test_bpClustering_LDADD = src/libfusionseq.la -lbios -lm

src_test_bpValidate_SOURCES = src/test/bpValidate.c
src_test_bpValidate_LDADD = src/libfusionseq.la -lbios -lm

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
