#include <stdio.h>

#include <bios/log.h>
#include <bios/format.h>

#include "gfr.h"



int main (int argc, char *argv[])
{
	GfrEntry *currGE;
	GfrInterRead *currGIR;
	int i;
	Stringa buffer;
	FILE *fp; //,*fp2;
	int count;

	count = 0;
	buffer = stringCreate (100);
	gfr_init ("-");
	puts (gfr_writeHeader ());
	while (currGE = gfr_nextEntry ()) {
		stringPrintf (buffer,"%s.bed",currGE->id);
		fp = fopen (string (buffer),"w");
	      
		if (fp == NULL ) { 
			die ("Unable to open BED files");
		}
		fprintf (fp,"browser full knownGene\n");
		fprintf (fp,"track name=\"Inter paired-end reads: %s\" visibility=2\n",currGE->id);
	
		for (i = 0; i < arrayMax (currGE->interReads); i++) {
			currGIR = arrp (currGE->interReads,i,GfrInterRead);
			fprintf (fp,"%s\t%d\t%d\t%s_R1_%d\n",currGE->chromosomeTranscript1,currGIR->readStart1,currGIR->readEnd1, currGE->id, i);
			fprintf (fp,"%s\t%d\t%d\t%s_R2_%d\n",currGE->chromosomeTranscript2,currGIR->readStart2,currGIR->readEnd2, currGE->id, i);
		}
		fclose (fp);
	
		puts (gfr_writeGfrEntry (currGE));
		count++;
	}
	gfr_deInit ();
	stringDestroy (buffer);
	warn ("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}

