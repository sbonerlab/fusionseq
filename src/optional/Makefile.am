bin_PROGRAMS = gfr2images \
		 gfr2bed \
		 gfr2fasta \
		 gfr2gff \
		 bp2wig \
		 export2mrf \
		 gfrConfidenceValueTranscript

AM_CPPFLAGS = -I ../

# -----------------------------------------------------------------------------
# OPTIONAL
# -----------------------------------------------------------------------------

gfr2images_SOURCES = gfr2images.c
gfr2images_LDADD = ../libfusionseq.la @OTHERLDFLAGS@ -lbios -lgd -lm

gfr2bed_SOURCES = gfr2bed.c
gfr2bed_LDADD = ../libfusionseq.la @OTHERLDFLAGS@ -lbios

gfr2fasta_SOURCES = gfr2fasta.c
gfr2fasta_LDADD = ../libfusionseq.la @OTHERLDFLAGS@ -lbios

gfr2gff_SOURCES = gfr2gff.c
gfr2gff_LDADD = ../libfusionseq.la @OTHERLDFLAGS@ -lbios

bp2wig_SOURCES = bp2wig.c
bp2wig_LDADD = ../libfusionseq.la @OTHERLDFLAGS@ -lbios

export2mrf_SOURCES = export2mrf.c
export2mrf_LDADD = ../libfusionseq.la @OTHERLDFLAGS@ -lbios

gfrConfidenceValueTranscript_SOURCES = gfrConfidenceValueTranscript.c
gfrConfidenceValueTranscript_LDADD = ../libfusionseq.la @OTHERLDFLAGS@ -lbios
