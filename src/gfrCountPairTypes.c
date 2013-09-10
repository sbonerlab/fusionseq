#include <bios/log.h>
#include <bios/format.h>

#include "gfr.h"

/**
   @file gfrCountPairTypes.c
   @brief Adding information regarding the inter-transcripts reads.
   @details It counts the number of inter-transcripts reads for each type of connection, i.e. exon-exon, exon-intron, etc.
   
   @author Andrea Sboner  (andrea.sboner.w [at] gmail.com).  
   @version 0.8
   @date 2013.09.10
   @remarks WARNings will be output to stdout to summarize the filter results.
   @pre A valid GFR file as input, including stdin.
 */


static int sortGfrInterReads (GfrInterRead *a, GfrInterRead *b)
{
	int diff;
  
	diff = a->pairType - b->pairType;
	if (diff != 0) {
		return diff;
	}
	diff = a->number1 - b->number1;
	if (diff != 0) {
		return diff;
	}
	return a->number2 - b->number2;
}



static void obtainPairCounts (GfrEntry *currGE)
{
	GfrPairCount *currPC;
	GfrInterRead *currGIR,*nextGIR;
	int i,j;

	currGE->pairCounts = arrayCreate (100,GfrPairCount);
	arraySort (currGE->interReads,(ARRAYORDERF)sortGfrInterReads);
	i = 0;
	while (i < arrayMax (currGE->interReads)) {
		currGIR = arrp (currGE->interReads,i,GfrInterRead);
		currPC = arrayp (currGE->pairCounts,arrayMax (currGE->pairCounts),GfrPairCount);
		currPC->number1 = currGIR->number1;
		currPC->number2 = currGIR->number2;
		currPC->pairType = currGIR->pairType;
		currPC->count = 1;
		j = i + 1;
		while (j < arrayMax (currGE->interReads)) {
			nextGIR = arrp (currGE->interReads,j,GfrInterRead);
			if (currGIR->pairType == nextGIR->pairType && currGIR->number1==nextGIR->number1 && currGIR->number2==nextGIR->number2) {
				currPC->count++;
			}
			else {
				break;
			}
			j++;
		}
		i = j;
	}
}



int main (int argc, char *argv[])
{
	GfrEntry *currGE;
	int count;

	gfr_init ("-");
	gfr_addNewColumnType (GFR_COLUMN_NAME_PAIR_COUNT);
	puts (gfr_writeHeader ());
	count = 0;
	while (currGE = gfr_nextEntry ()){
		obtainPairCounts (currGE);
		puts (gfr_writeGfrEntry (currGE)); fflush (stdout);
		count++;
	}
	gfr_deInit ();
	warn ("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}


