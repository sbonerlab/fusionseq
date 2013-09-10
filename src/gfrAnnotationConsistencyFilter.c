#include <bios/log.h>
#include <bios/format.h>

#include "gfr.h"

/**
   @file gfrAnnotationConsistencyFilter.c
   @brief It removes candidates involving genes with specific description, such as ribosomal, pseudogenes, etc. 
   @details It removes candidates involving genes with specific text in the gene description, such as ribosomal, pseudogenes, etc
   
   @author Andrea Sboner  (andrea.sboner.w [at] gmail.com).  
   @version 0.8
   @date 2013.09.10
   @remarks WARNings will be output to stdout to summarize the filter results.
   @prestring string the element to remove, ex. pseudogene.
   @pre A valid GFR file as input, including stdin.
 */

int main (int argc, char *argv[])
{
	GfrEntry *currGE;
	int count;
	int countRemoved;
 
	if (argc != 2) {
		usage ("%s <string>",argv[0]);
	}
	count = 0;
	countRemoved = 0;
	gfr_init ("-");
	puts (gfr_writeHeader ());
	while (currGE = gfr_nextEntry ()){
		if (currGE->descriptionTranscript1 == NULL ||
				currGE->descriptionTranscript2 == NULL) {
			die ("Transcript description is missing");
		}
		if (strCaseStr (currGE->descriptionTranscript1,argv[1]) ||
		    strCaseStr (currGE->descriptionTranscript2,argv[1])) {
			countRemoved++;
			continue;
		}
		puts (gfr_writeGfrEntry (currGE));
		count++;
	}
	gfr_deInit ();
	warn ("%s_string: %s",argv[0],argv[1]);
	warn ("%s_numRemoved: %d",argv[0],countRemoved);
	warn ("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}

