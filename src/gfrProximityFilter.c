#include <bios/log.h>
#include <bios/format.h>

#include "gfr.h"

/**
   @file gfrProximityFilter.c
   @brief It removes 'cis' candidates that are likely due to mis-annotation of the 5' or 3' ends of the genes.
   @details It removes 'cis' candidates that are likely due to mis-annotation of the 5' or 3' ends of the genes.
   
   @author Andrea Sboner  (andrea.sboner.w [at] gmail.com).  
   @version 0.8
   @date 2013.09.10
   @remarks WARNings will be output to stdout to summarize the filter results.
   @pre [in] offset the minimum distance (in nucleotides) between the two genes to keep the candidate
   @pre A valid GFR file as input, including stdin.
 */

int main (int argc, char *argv[])
{
	GfrEntry *currGE;
	int offset;
	int count;
	int countRemoved;

	if (argc != 2) {
		usage ("%s <offset>",argv[0]);
	}
	count = 0;
	countRemoved = 0;
	offset = atoi (argv[1]);
	gfr_init ("-");
	puts (gfr_writeHeader ());
	while (currGE = gfr_nextEntry ()) {
		if (strEqual (currGE->fusionType,"cis") && 
		    currGE->strandTranscript1 != currGE->strandTranscript2 &&
		    (currGE->startTranscript2 - currGE->endTranscript1) < offset) {
			countRemoved++;
			continue;
		}
		puts (gfr_writeGfrEntry (currGE));
		count++;
	}
	gfr_deInit ();
	warn ("%s_offset: %d",argv[0],offset);
	warn ("%s_numRemoved: %d",argv[0],countRemoved);
	warn ("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}

