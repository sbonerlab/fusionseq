#include <bios/log.h>
#include <bios/format.h>
#include <bios/intervalFind.h>
#include <bios/confp.h>
#include "gfr.h"

/**
   @file gfrPseudogenesFlter.c
   @brief Filter to remove artifacts due to mis-alignment to pseudogenes.
   @details It removes candidates with reads overlapping repetitive sequences. It looks at non-exonic reads and, if some overlap exists with repetitive regions, the reads are excluded and the number of inter-reads is updated accordingly. MAX_OVERLAP_ALLOWED will determine if an overlap triggers the removal of the read. If the remaining number of reads is below the threshold (minNumberOfReads), the fusion candidate is removed.
   
   @author Andrea Sboner  (andrea.sboner.w [at] gmail.com).  
   @version 0.8
   @date 2013.10.30
   @remarks WARNings will be output to stdout to summarize the filter results.
   @pre A valid GFR file as input, including stdin.
   @pre repeatMasker.interval The repetitive regions in interval format; typically from RepeatMasker, defined in .fusionseqrc
   @pre [in] minNumberInterReads An integer representing the minimum number of reads to keep the fusion candidate.
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
	int count,countRemoved;
	int i, j;
	Array intervals;
	GfrInterRead *currGIR;
	int totalOverlaps;
	float minNumInterReads;
	float numberOfInters;
	Stringa buffer = stringCreate(100);

	if ((Conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL) {
	  die("%s:\tCannot find .fusionseqrc: %s", argv[0], getenv("FUSIONSEQ_CONFPATH"));
	  return EXIT_FAILURE;
	}
	if( confp_get( Conf, "PSEUDOGENE_DIR")==NULL ) {
	  die("%s:\tCannot find PSEUDOGENE_DIR in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
	  return EXIT_FAILURE;
	}
	if( confp_get( Conf, "PSEUDOGENE_FILENAME")==NULL ) {
	  die("%s:\tCannot find PSEUDOGENE_FILENAME in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
	  return EXIT_FAILURE;
	}

	if (argc != 2) {
	  usage ("%s <minNumInterReads>",argv[0]);
	}
	stringPrintf( buffer, "%s/%s", confp_get( Conf, "PSEUDOGENE_DIR"), confp_get( Conf, "PSEUDOGENE_FILENAME") );
	intervalFind_addIntervalsToSearchSpace (string(buffer),0);
	stringDestroy(buffer); 
	minNumInterReads = atof (argv[1]);
	count = 0;
	countRemoved = 0;
	gfr_init ("-");
	puts (gfr_writeHeader ());
	while (currGE = gfr_nextEntry ()){
	  int readLength = strlen( arru( currGE->readsTranscript1, 0, Texta ) );
	  numberOfInters = (float) currGE->numInter;
	  for (i = 0; i < arrayMax (currGE->interReads); i++) {
	    currGIR = arrp (currGE->interReads,i,GfrInterRead);
	    //if (currGIR->pairType == GFR_PAIR_TYPE_EXONIC_EXONIC) {
	    // continue;
	    // }
	    totalOverlaps = 0;
	    intervals = intervalFind_getOverlappingIntervals (currGE->chromosomeTranscript1,currGIR->readStart1,currGIR->readEnd1);
	    for(j=0; j < arrayMax( intervals ); j++) {
	      Interval* currInterval = arru( intervals, j, Interval*);
	      totalOverlaps = getNucleotideOverlap ( currGIR->readStart1,currGIR->readEnd1, currInterval );
	    }
	    if ( totalOverlaps >  ( ((double)(readLength)) * strtod(confp_get(Conf, "MAX_OVERLAP_ALLOWED"), NULL) ) ) {
	      currGE->numInter-= getNumInter( currGIR, readLength );
	      currGIR->flag = 1;
	      continue;
	    }
	    intervals = intervalFind_getOverlappingIntervals (currGE->chromosomeTranscript2,currGIR->readStart2,currGIR->readEnd2);
	     for(j=0; j < arrayMax( intervals ); j++) {
	      Interval* currInterval = arru( intervals, j, Interval*);
	      totalOverlaps = getNucleotideOverlap ( currGIR->readStart2,currGIR->readEnd2, currInterval );
	    }
	    if ( totalOverlaps >  ( ((double)(readLength)) * strtod(confp_get(Conf, "MAX_OVERLAP_ALLOWED"), NULL) ) ) {
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
	warn ( "%s_interval: %s/%s", argv[0], confp_get( Conf, "PSEUDOGENE_DIR"), confp_get( Conf, "PSEUDOGENE_FILENAME") );
	warn ("%s_numRemoved: %d",argv[0],countRemoved);
	warn ("%s_numGfrEntries: %d",argv[0],count);
	return 0;
}
