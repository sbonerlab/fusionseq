#include <stdlib.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>

#include "util.h"
#include "gfr.h"

/**
   @file gfrAddInfo.c
   @brief It includes additional information about the fusion transcript candidates such as gene symbols and gene description. 
   @details It includes additional information about the fusion transcript candidates such as gene symbols and gene description. This is a pre-requisite for gfrBlackListFilter and gfrAnnotationConsistencyFilter.
   
   @author Andrea Sboner  (andrea.sboner.w [at] gmail.com).  
   @version 0.8
   @date 2013.09.10
   @remarks WARNings will be output to stdout to summarize the filter results.
   @pre An external file that includes all descriptive information about the annotation set. The format of this file should follow kgXref.txt (from UCSC). Indeed, we use kgXref.txt for human, however, this could be modified by the user.
   @pre A valid GFR file as input, including stdin.
 */

int main (int argc, char *argv[])
{
	GfrEntry *currGE;
	Array kgXrefs;
	Stringa buffer;
	int count;

	config *Conf;
	
	if ((Conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL) {
	  die("%s:\tCannot find .fusionseqrc: %s", argv[0], getenv("FUSIONSEQ_CONFPATH"));
	  return EXIT_FAILURE;
	}
	if( confp_get( Conf, "ANNOTATION_DIR")==NULL ) {
	  die("%s:\tCannot find ANNOTATION_DIR in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
	  return EXIT_FAILURE;
	}
	if( confp_get( Conf, "KNOWN_GENE_XREF_FILENAME")==NULL ) {
	  die("%s:\tCannot find KNOWN_GENE_XREF_FILENAME in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
	  return EXIT_FAILURE;
	}

	buffer = stringCreate (100);
	stringPrintf (buffer,"%s/%s",
		      confp_get(Conf, "ANNOTATION_DIR"),
		      confp_get(Conf, "KNOWN_GENE_XREF_FILENAME"));

	kgXrefs = util_readKnownGeneXrefs (string (buffer));
	arraySort (kgXrefs,(ARRAYORDERF)sortKgXrefsByTranscriptName);
	stringDestroy (buffer);

	count = 0;
	gfr_init ("-");
	gfr_addNewColumnType (GFR_COLUMN_NAME_GENE_SYMBOL_TRANSCRIPT1);
	gfr_addNewColumnType (GFR_COLUMN_NAME_GENE_SYMBOL_TRANSCRIPT2);
	gfr_addNewColumnType (GFR_COLUMN_NAME_DESCRIPTION_TRANSCRIPT1);
	gfr_addNewColumnType (GFR_COLUMN_NAME_DESCRIPTION_TRANSCRIPT2);
	puts (gfr_writeHeader ());
	
	while (currGE = gfr_nextEntry ()){
		transcript2geneSymbolAndGeneDescription (kgXrefs,currGE->nameTranscript1,&currGE->geneSymbolTranscript1,&currGE->descriptionTranscript1);
		transcript2geneSymbolAndGeneDescription (kgXrefs,currGE->nameTranscript2,&currGE->geneSymbolTranscript2,&currGE->descriptionTranscript2);
		puts (gfr_writeGfrEntry (currGE));
		count++;
	}
	gfr_deInit ();
	warn ("%s_numGfrEntries: %d",argv[0],count);
	confp_close(Conf);

	return EXIT_SUCCESS;
}

