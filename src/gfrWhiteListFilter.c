#include <stdio.h>

#include <bios/log.h>
#include <bios/format.h>
#include <bios/common.h>
#include <bios/linestream.h>
#include <bios/intervalFind.h>

#include "gfr.h"

typedef struct {
  char* gene1;
  char* gene2;
} WLGeneEntry;

static int sortWhiteListByName1 (WLGeneEntry *a, WLGeneEntry *b) 
{
  int res = strcmp ( a->gene1, b->gene1);
  if( res==0 ) res = strcmp ( a->gene2, b->gene2 );
  return (res); //(strcmp ( a->gene1, b->gene1));
}

int main (int argc, char *argv[])
{
  GfrEntry *currGE;
  WLGeneEntry *currWLGE;
  WLGeneEntry currGeneQuery;
  FILE *fp, *fTempCoordinates;
  char *line;
  int count;

  Stringa buffer=stringCreate(10);

  int index, i, j;
  WordIter w;
  Array whiteGeneList = arrayCreate(20, WLGeneEntry);

  if (argc != 2) {
    usage ("%s <whiteList.txt>",argv[0]);
  }  
  fp = fopen( argv[1], "r" );
  if( !fp )  die("Unable to open file: %s", argv[1]);

  fTempCoordinates= fopen( "tmp_coordinates.interval", "w");
  if( !fTempCoordinates )  die("Unable to open file: tmp_coordinates.interval");
  
  // reading whitelist file
  LineStream ls = ls_createFromFile( argv[1] );
  index=0;
  while( line = ls_nextLine(ls) ) {
    if( strStartsWith( line, "#") ) // comments
      continue;
    if( strStartsWith( line, "@" ) ) { // coordinates
      Interval* currInterval;
      SubInterval *subInt;
      stringPrintf(buffer, "pair_%d", index );
      AllocVar( currInterval );
      w = wordIterCreate( line, "@:-\t", 1);
      currInterval->source=0;
      currInterval->name=hlr_strdup( string( buffer ) );
      currInterval->chromosome = hlr_strdup( wordNext(w) ); // chr1
      currInterval->strand = '.';
      currInterval->start = atoi( wordNext(w) ); // start1
      currInterval->end   = atoi( wordNext(w) ); // end1
      currInterval->subIntervalCount = 1;
      currInterval->subIntervals = arrayCreate( 1, SubInterval);
      subInt = arrayp( currInterval->subIntervals, arrayMax( currInterval->subIntervals), SubInterval );
      subInt->start = currInterval->start;
      subInt->end = currInterval->end;
      fprintf( fTempCoordinates, "%s\n", intervalFind_writeInterval( currInterval ) ); 
      freeMem( currInterval );
      AllocVar( currInterval );
      //wordNext(w); // @
      currInterval->source=0;
      currInterval->name=hlr_strdup( string( buffer ) );
      currInterval->chromosome = hlr_strdup( wordNext(w) ); // chr2
      currInterval->strand = '.';
      currInterval->start = atoi( wordNext(w) ); // start2
      currInterval->end   = atoi( wordNext(w) ); // end2
      currInterval->subIntervalCount = 1;
      currInterval->subIntervals = arrayCreate( 1, SubInterval);
      subInt = arrayp( currInterval->subIntervals, arrayMax( currInterval->subIntervals), SubInterval );
      subInt->start = currInterval->start;
      subInt->end = currInterval->end; 
      fprintf( fTempCoordinates, "%s\n", intervalFind_writeInterval( currInterval ) );
      freeMem( currInterval );
      index++;
    } else { // genes symbols
      w = wordIterCreate( line, "\t", 1);
      currWLGE = arrayp( whiteGeneList, arrayMax(whiteGeneList), WLGeneEntry);
      currWLGE->gene1 = hlr_strdup ( wordNext(w) );
      currWLGE->gene2 = hlr_strdup ( wordNext(w) );   
    }
    wordIterDestroy(w);
  }
  stringDestroy( buffer );
  fclose(fp);
  fclose(fTempCoordinates);
  intervalFind_addIntervalsToSearchSpace( "tmp_coordinates.interval", 0);
  arraySort( whiteGeneList, (ARRAYORDERF) sortWhiteListByName1);

  // beginFiltering
  count = 0;
  gfr_init ("-");
  puts (gfr_writeHeader ());
  while (currGE = gfr_nextEntry ()) { // reading the gfr
    // creating a new query to the gene white list
    currGeneQuery.gene1 = currGE->geneSymbolTranscript1;
    currGeneQuery.gene2 = currGE->geneSymbolTranscript2;
    // searching against read_1/read_2
    int res = arrayFind( whiteGeneList, &currGeneQuery, &index, (ARRAYORDERF) sortWhiteListByName1);  
    if( res ) { // found, write the instance to stdout, update the counts 
      puts (gfr_writeGfrEntry (currGE));
      count++;
      continue;
    } else { //not found: then searching against read_2/read_1
      currGeneQuery.gene1 = currGE->geneSymbolTranscript2;
      currGeneQuery.gene2 = currGE->geneSymbolTranscript1;
      res =  arrayFind( whiteGeneList, &currGeneQuery, &index, (ARRAYORDERF) sortWhiteListByName1 );
      
      if( res ) { // found, write the instance to stdout, update the counts
	puts (gfr_writeGfrEntry (currGE));
	count++;	
	continue;
      } 
    }
    // not found in the genes; search the coordinates
    Array overlap1 = arrayCopy(intervalFind_getOverlappingIntervals( currGE->chromosomeTranscript1, currGE->startTranscript1, currGE->endTranscript1 ));
    Array overlap2 = arrayCopy(intervalFind_getOverlappingIntervals( currGE->chromosomeTranscript2, currGE->startTranscript2, currGE->endTranscript2 ));
    if( (arrayMax(overlap1) >0 ) && (arrayMax( overlap2 ) > 0) ) {
      for( i=0; i<arrayMax( overlap1 ); i++ ) {
	for( j=0; j<arrayMax( overlap2 ); j++ ) {
	  Interval *currInterval1 = arru( overlap1, i, Interval*);
	  Interval *currInterval2 = arru( overlap2, j, Interval*);
	  if( strEqual( currInterval1->name, currInterval2->name ) ) {
	    puts( gfr_writeGfrEntry( currGE ));
	    count++;
	    i=arrayMax( overlap1 );
	    j=arrayMax( overlap2 );
	  }
	}
      }
    }
  }	           
  gfr_deInit ();
  arrayDestroy( whiteGeneList );
  hlr_system("rm tmp_coordinates.interval", 1);
  warn ("%s_WhiteListFilter: %s",argv[0], argv[1]);
  warn ("%s_numGfrEntries: %d",argv[0],count);
  return 0;
}

