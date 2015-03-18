#include <bios/log.h>
#include <bios/format.h>
#include <bios/intervalFind.h>
#include <bios/confp.h>
#include "gfr.h"

/**
   @file gfrSequenceComplexityFilter.c
   @brief Filter to remove artifacts due to low complexity of the reads.
   @details It removes candidates with reads with low complexity such as CCCCTTTTTCAAAAAAACAAAAAAAAAAAAAAAACACACAAAACAAAA. .
   
   @author Andrea Sboner  (andrea.sboner.w [at] gmail.com).  
   @version 0.8
   @date 2013.10.30
   @remarks WARNings will be output to stdout to summarize the filter results.
   @pre A valid GFR file as input, including stdin.
 */

static config *Conf = NULL; /**< Pointer to configuration file .fusionseqrc  */

static float getNumInter( GfrInterRead* currInter, int readLength ) { // computes the correct number of the inters by considering split reads on splice junctions.
  float numInter=0.0;
  if( (currInter->readEnd1 - currInter->readStart1 + 1) != readLength &
      (currInter->readEnd2 - currInter->readStart2 + 1) != readLength  ) {
    numInter += 0.25;
  } else if ( (currInter->readEnd1 - currInter->readStart1 + 1) != readLength  |
	      (currInter->readEnd2 - currInter->readStart2 + 1) != readLength ) {
    numInter += 0.5;
  } else {
    numInter += 1.0;
  }
  return( numInter );
}

int getNucleotideOverlap( int start, int end, Interval* currInterval ) {
  int k;
  int overlap=0;
  for( k=0; k<arrayMax( currInterval->subIntervals); k++ ) {
    SubInterval* currSubInterval = arrp( currInterval->subIntervals, k, SubInterval);
    overlap+=positiveRangeIntersection( start, end, currSubInterval->start, currSubInterval->end);
  }
  return overlap;
}

int main (int argc, char *argv[]) {
	GfrEntry *currGE;
	GfrInterRead *currGIR;
	int count,countRemoved, countChanges, readLength;
	int i, j;
	float numberOfInters;
	char* currRead;
	Stringa buffer = stringCreate(100);
	float minNumInterReads;

	if ((Conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL) {
	  die("%s:\tCannot find .fusionseqrc: %s", argv[0], getenv("FUSIONSEQ_CONFPATH"));
	  return EXIT_FAILURE;
	}
	if (argc != 2) {
	  usage ("%s <minNumInterReads>",argv[0]);
	}
	minNumInterReads = atof (argv[1]);
	count = 0;
	countRemoved = 0;
	gfr_init ("-");
	puts (gfr_writeHeader ());
	while (currGE = gfr_nextEntry ()){
	  for( i=0; i<arrayMax(currGE->readsTranscript1); i++ ) {
	    currRead = arru( currGE->readsTranscript1, i, char* );
	    countChanges = 0;
	    readLength = strlen( currRead );
	    for( j=0; j<(readLength-1); j++) {
	      if( currRead[j] != currRead[j+1]) countChanges++;
	    }
	    if( countChanges < (int)(readLength/2) ) {
	      currGIR = arrp (currGE->interReads,i,GfrInterRead);
	      currGE->numInter-= getNumInter( currGIR, readLength );
	      currGIR->flag = 1;
	      continue;
	    }
	  }
	  for( i=0; i<arrayMax(currGE->readsTranscript2); i++ ) {
	    currRead = arru( currGE->readsTranscript2, i, char* );
	    countChanges = 0;
	    readLength = strlen( currRead );
	    for( j=0; j<(readLength-1); j++) {
	      if( currRead[j] != currRead[j+1]) countChanges++;
	    }
	    if( countChanges < (int)(readLength/2) ) {
	      currGIR = arrp (currGE->interReads,i,GfrInterRead);
	      currGE->numInter-= getNumInter( currGIR, readLength );
	      currGIR->flag = 1;
	      continue;
	    }
	  }
	  if (currGE->numInter < (float)minNumInterReads) { 
	    countRemoved++;
	    continue;
	  }
	  puts (gfr_writeGfrEntry (currGE));
	  count++;
	}
	gfr_deInit ();
	warn ("%s_numRemoved: %d",argv[0],countRemoved);
	warn ("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}
