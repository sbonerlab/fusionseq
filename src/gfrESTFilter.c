#include "log.h"
#include "format.h"
#include "gfr.h"
#include "intervalFind.h"

typedef struct {
  char* gene1;
  char* gene2;
} BLEntry;

static void findCoordinates( GfrEntry *gfrE, int *start1, int *end1, int *start2, int *end2 )
{
  GfrInterRead *gfrIR;
  int i;
  *start1 = arrp( gfrE->interReads, 0, GfrInterRead )->readStart1;
  *end1 = arrp( gfrE->interReads, 0, GfrInterRead )->readEnd1;
  *start2 = arrp( gfrE->interReads, 0, GfrInterRead )->readStart2;
  *end2 = arrp( gfrE->interReads, 0, GfrInterRead )->readEnd2;
  for( i = 1; i< arrayMax( gfrE->interReads); i++ ) {
    gfrIR = arrp( gfrE->interReads, i, GfrInterRead );
    if( gfrIR->readStart1 < *start1 )
      *start1 = gfrIR->readStart1;
    if( gfrIR->readStart2 < *start2 )
      *start2 = gfrIR->readStart2;
    if( gfrIR->readEnd1 > *end1 )
      *end1 = gfrIR->readEnd1;
    if( gfrIR->readEnd2 > *end2 )
      *end2 = gfrIR->readEnd2;
  }
}


int main (int argc, char *argv[])
{
  GfrEntry *currGE;
  int count;
  int countRemoved; 
  int i, j;
  int foundEST;
 
  if (argc != 2) {
    usage ("%s <EST.interval>",argv[0]);
  }  
  intervalFind_addIntervalsToSearchSpace( argv[1], 0);	

  // beginFiltering
  count = 0;
  countRemoved = 0;
  gfr_init ("-");
  puts (gfr_writeHeader ());
  while (currGE = gfr_nextEntry ()) { // reading the gfr
    foundEST = 0;
    if( strEqual( currGE->fusionType, "cis" ) ) {
      if( ! strEqual( currGE->chromosomeTranscript1, currGE->chromosomeTranscript2 ) )
	die("The two genes are not on the same chromosomes: %s - %s",  currGE->chromosomeTranscript1, currGE->chromosomeTranscript2 );
      int start1, end1, start2, end2;
      findCoordinates( currGE, &start1, &end1, &start2, &end2 );
      
      Array intervals1 = arrayCopy( intervalFind_getOverlappingIntervals( currGE->chromosomeTranscript1, start1, end1 ) ); 
      Array intervals2 = intervalFind_getOverlappingIntervals( currGE->chromosomeTranscript2, start2, end2 );
      for( i=0; i<arrayMax( intervals1 ); i++ ) {
	Interval* currInterval1 = arru( intervals1, i, Interval* );
	for( j=0; j<arrayMax ( intervals2 ); j++ ) {
	  Interval* currInterval2 = arru( intervals2, j, Interval* );
	  if( currInterval1==currInterval2 ) {
	    foundEST = 1;
	    i = arrayMax( intervals1 );
	    j = arrayMax( intervals2 );
	  }
	}
      }
      arrayDestroy( intervals1 );
      
    }
    if( foundEST )
      countRemoved++;
    else {
      puts (gfr_writeGfrEntry (currGE));
      count++;
    }
  }	           
  gfr_deInit ();
  warn ("%s_EST_data: %s",argv[0], argv[1]);
  warn ("%s_numRemoved: %d",argv[0], countRemoved);
  warn ("%s_numGfrEntries: %d",argv[0],count);
  return 0;
}

