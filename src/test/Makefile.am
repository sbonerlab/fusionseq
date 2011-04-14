bin_PROGRAMS =  gfSimulation \
		mrfPseudogeneOverlap \
		 quantifierAddInfo \
		 bp2fasta \
		 bpClustering \
		 bpValidate

AM_CPPFLAGS = -I ../

# -----------------------------------------------------------------------------
# TEST
# -----------------------------------------------------------------------------

gfSimulation_SOURCES = gfSimulation.c
gfSimulation_LDADD = ../libfusionseq.la -lbios -lm

mrfPseudogeneOverlap_SOURCES = mrfPseudogeneOverlap.c
mrfPseudogeneOverlap_LDADD = ../libfusionseq.la -lbios -lmrf -lm

quantifierAddInfo_SOURCES = quantifierAddInfo.c
quantifierAddInfo_LDADD = ../libfusionseq.la -lbios

bp2fasta_SOURCES = bp2fasta.c
bp2fasta_LDADD = ../libfusionseq.la -lbios -lm

bpClustering_SOURCES = bpClustering.c
bpClustering_LDADD = ../libfusionseq.la -lbios -lm

bpValidate_SOURCES = bpValidate.c
bpValidate_LDADD = ../libfusionseq.la -lbios -lm
