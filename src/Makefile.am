# -----------------------------------------------------------------------------
# Common library
# -----------------------------------------------------------------------------

lib_LTLIBRARIES = libfusionseq.la
libfusionseq_la_SOURCES = bp.c \
			  gfr.c \
			  sqvCircos.c \
			  sqvUtil.c \
			  sqvWeb.c \
			  util.c

libfusionseq_la_LIBADD = -lm -lbios

include_HEADERS = bp.h \
		  gfr.h \
		  sqvCircos.h \
		  sqvUtil.h \
		  sqvWeb.h \
		  util.h

# -----------------------------------------------------------------------------
# FusionSeq programs
# -----------------------------------------------------------------------------

bin_PROGRAMS = geneFusions \
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
	       gfrConfidenceValues \
	       validateBpJunctions \
	       gfrExpressionConsistencyFilter \
	       gfrMitochondrialFilter \
	       gfrRepeatMaskerFilter \
	       gfrCountPairTypes \
	       gfrRibosomalFilter \
	       gfrSpliceJunctionFilter \
	       gfrClassify \
	       gfrWhiteListFilter

EXTRA_PROGRAMS = gfr2images \
		 gfr2bed \
		 gfr2fasta \
		 gfr2gff \
		 bp2wig \
		 plotIntraDistribution \
		 export2mrf \
		 gfrConfidenceValueTranscript \
		 geneFusions_cgi \
		 showDetails_cgi \
		 seqViz_cgi \
		 findFusionPartner.cgi \
		 gfSimulation \
		 mrfPseudogeneOverlap \
		 quantifierAddInfo \
		 bp2fasta \
		 bpClustering \
		 bpValidate

optional: gfr2images \
	  gfr2bed \
	  gfr2fasta \
	  gfr2gff \
	  bp2wig \
	  plotIntraDistribution \
	  export2mrf \
	  gfrConfidenceValueTranscript

cgi: geneFusions_cgi \
     showDetails_cgi \
     seqViz_cgi \
     findFusionPartner.cgi

test: gfSimulation \
      mrfPseudogeneOverlap \
      quantifierAddInfo \
      bp2fasta \
      bpClustering \
      bpValidate

CLEANFILES = $(EXTRA_PROGRAMS)

.PHONY: optional cgi test

# -----------------------------------------------------------------------------
# CORE: Finding candidate fusions
# -----------------------------------------------------------------------------

geneFusions_SOURCES = geneFusions.c
geneFusions_LDADD = libfusionseq.la -lbios -lmrf -lm

gfrExpressionConsistencyFilter_SOURCES = gfrExpressionConsistencyFilter.c
gfrExpressionConsistencyFilter_LDADD = libfusionseq.la -lbios -lm

gfrAnnotationConsistencyFilter_SOURCES = gfrAnnotationConsistencyFilter.c
gfrAnnotationConsistencyFilter_LDADD = libfusionseq.la -lbios -lm

gfrRibosomalFilter_SOURCES = gfrRibosomalFilter.c
gfrRibosomalFilter_LDADD = libfusionseq.la -lbios

gfrMitochondrialFilter_SOURCES = gfrMitochondrialFilter.c
gfrMitochondrialFilter_LDADD = libfusionseq.la -lbios

gfrLargeScaleHomologyFilter_SOURCES = gfrLargeScaleHomologyFilter.c
gfrLargeScaleHomologyFilter_LDADD = libfusionseq.la -lbios

gfrSmallScaleHomologyFilter_SOURCES = gfrSmallScaleHomologyFilter.c
gfrSmallScaleHomologyFilter_LDADD = libfusionseq.la -lbios

gfrAbnormalInsertSizeFilter_SOURCES = gfrAbnormalInsertSizeFilter.c
gfrAbnormalInsertSizeFilter_LDADD = libfusionseq.la -lbios

gfrPCRFilter_SOURCES = gfrPCRFilter.c
gfrPCRFilter_LDADD = libfusionseq.la -lbios

gfrBlackListFilter_SOURCES = gfrBlackListFilter.c
gfrBlackListFilter_LDADD = libfusionseq.la -lbios

gfrProximityFilter_SOURCES = gfrProximityFilter.c
gfrProximityFilter_LDADD = libfusionseq.la -lbios

gfrAddInfo_SOURCES = gfrAddInfo.c
gfrAddInfo_LDADD = libfusionseq.la -lbios

gfrConfidenceValues_SOURCES = gfrConfidenceValues.c
gfrConfidenceValues_LDADD = libfusionseq.la -lbios -lmrf

gfrRepeatMaskerFilter_SOURCES = gfrRepeatMaskerFilter.c
gfrRepeatMaskerFilter_LDADD = libfusionseq.la -lbios -lm

gfrCountPairTypes_SOURCES = gfrCountPairTypes.c
gfrCountPairTypes_LDADD = libfusionseq.la -lbios -lm

gfrSpliceJunctionFilter_SOURCES = gfrSpliceJunctionFilter.c
gfrSpliceJunctionFilter_LDADD = libfusionseq.la -lbios

gfrClassify_SOURCES = gfrClassify.c
gfrClassify_LDADD = libfusionseq.la -lbios -lm

gfrWhiteListFilter_SOURCES = gfrWhiteListFilter.c
gfrWhiteListFilter_LDADD = libfusionseq.la -lbios

# -----------------------------------------------------------------------------
# CORE: Identifying sequences of the junction
# -----------------------------------------------------------------------------

gfr2bpJunctions_SOURCES = gfr2bpJunctions.c
gfr2bpJunctions_LDADD = libfusionseq.la -lbios -lm

bowtie2bp_SOURCES = bowtie2bp.c
bowtie2bp_LDADD = libfusionseq.la -lbios

validateBpJunctions_SOURCES = validateBpJunctions.c
validateBpJunctions_LDADD = libfusionseq.la -lbios

bpFilter_SOURCES = bpFilter.cpp
bpFilter_CXXFLAGS = -D_REENTRANT -pthread -m64
bpFilter_LDADD = libfusionseq.la -lstdc++ -lbios -L@ROOTLIBDIR@ @ROOTGLIBS@  @ROOTLIBS@ -ldl

bp2alignment_SOURCES = bp2alignment.c
bp2alignment_LDADD = libfusionseq.la -lbios

# -----------------------------------------------------------------------------
# OPTIONAL
# -----------------------------------------------------------------------------

plotIntraDistribution_SOURCES = plotIntraDistribution.cpp
plotIntraDistribution_CXXFLAGS = -D_REENTRANT -pthread -m64
plotIntraDistribution_LDADD = libfusionseq.la -lstdc++ -lbios -L@ROOTLIBDIR@ @ROOTGLIBS@ @ROOTLIBS@ -ldl

gfr2images_SOURCES = gfr2images.c
gfr2images_LDADD = libfusionseq.la -lbios -lgd -lm

gfr2bed_SOURCES = gfr2bed.c
gfr2bed_LDADD = libfusionseq.la -lbios

gfr2fasta_SOURCES = gfr2fasta.c
gfr2fasta_LDADD = libfusionseq.la -lbios

gfr2gff_SOURCES = gfr2gff.c
gfr2gff_LDADD = libfusionseq.la -lbios

bp2wig_SOURCES = bp2wig.c
bp2wig_LDADD = libfusionseq.la -lbios

export2mrf_SOURCES = export2mrf.c
export2mrf_LDADD = libfusionseq.la -lbios

gfrConfidenceValueTranscript_SOURCES = gfrConfidenceValueTranscript.c
gfrConfidenceValueTranscript_LDADD = libfusionseq.la -lbios

# -----------------------------------------------------------------------------
# CGI
# -----------------------------------------------------------------------------

geneFusions_cgi_SOURCES = geneFusions_cgi.c
geneFusions_cgi_LDADD = libfusionseq.la -lbios

showDetails_cgi_SOURCES = showDetails_cgi.c
showDetails_cgi_LDADD = libfusionseq.la -lbios

seqViz_cgi_SOURCES = seqViz_cgi.c
seqViz_cgi_LDADD = libfusionseq.la -lbios -lmrf

findFusionPartner_cgi_SOURCES = findFusionPartner_cgi.c
findFusionPartner_cgi_LDADD = libfusionseq.la -lbios

# -----------------------------------------------------------------------------
# TEST
# -----------------------------------------------------------------------------

gfSimulation_SOURCES = gfSimulation.c
gfSimulation_LDADD = libfusionseq.la -lbios -lm

mrfPseudogeneOverlap_SOURCES = mrfPseudogeneOverlap.c
mrfPseudogeneOverlap_LDADD = libfusionseq.la -lbios -lmrf -lm

quantifierAddInfo_SOURCES = quantifierAddInfo.c
quantifierAddInfo_LDADD = libfusionseq.la -lbios

bp2fasta_SOURCES = bp2fasta.c
bp2fasta_LDADD = libfusionseq.la -lbios -lm

bpClustering_SOURCES = bpClustering.c
bpClustering_LDADD = libfusionseq.la -lbios -lm

bpValidate_SOURCES = bpValidate.c
bpValidate_LDADD = libfusionseq.la -lbios -lm

