#include <bios/log.h>
#include <bios/format.h>
#include "gfr.h"

/**
   @file gfrAbnormalInsertSizeFilter.c
   @brief It removes candidates with an insert-size bigger than the normal insert-size.
   @details It removes candidates with an insert-size bigger than the normal insert-size. The fusion candidate insert-size is computed on the minimal fusion transcript fragment.
   
   @author Andrea Sboner  (andrea.sboner.w [at] gmail.com).  
   @version 0.8
   @date 2013.09.10
   @remarks WARNings will be output to stdout to summarize the filter results.
   @pre A valid GFR file as input, including stdin.
 */

int main (int argc, char *argv[])
{
	GfrEntry *currGE;
	double pvalueCutOff;
	int count;
	int countRemoved;

	if (argc != 2) {
		usage ("%s <pvalueCutOff>",argv[0]);
	}
	count = 0;
	countRemoved = 0;
	pvalueCutOff = atof (argv[1]);
	gfr_init ("-");
	puts (gfr_writeHeader ());
	while (currGE = gfr_nextEntry ()){
		if ((MAX(currGE->pValueAB,currGE->pValueBA) < pvalueCutOff) && 
		    (currGE->pValueAB != -1) ) {
			countRemoved++;
			continue;
		} 
		puts (gfr_writeGfrEntry (currGE));
		count++;
	}
	gfr_deInit ();
	warn ("%s_pvalueCutOff: %f",argv[0],pvalueCutOff);
	warn ("%s_numRemoved: %d",argv[0],countRemoved);
	warn ("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}

