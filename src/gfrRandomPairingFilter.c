#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <bios/log.h>
#include <bios/format.h>
#include <bios/intervalFind.h>
#include <bios/numUtil.h>
#include <bios/linestream.h>
#include <mrf/mrf.h>

#include <uthash.h>

struct geneFusionsCount {
  const char* geneID;
  int numFusions;
  UT_hash_handle hh;
};


#include "gfr.h"


/**
  @file gfrRandomPairingFilter.c 
  @brief Elimination of candidates with multiple partners.
  @details  Elimination of candidates with multiple partners likely due to random pairing during sample preparation. 
  @author Andrea Sboner (andrea.sboner.w [at] gmail.com).
  @version 0.8						
 
  @attention It requires a GFR file from stdin: @code $ gfrRandomPairingFilter < file.gfr @endcode

  @remarks It outputs to stdin and stderr. The stderr messages are mostly for logging purposes.
  @copyright GNU license: free for academic use
 */
#define MILLION 1000000.0


struct geneFusionsCount *myGeneFusionsCount = NULL;

void addGene( char* gene ) {
  struct geneFusionsCount *s;
  s = (struct geneFusionsCount*)malloc( sizeof( struct geneFusionsCount ) );
  s->geneID = hlr_strdup(gene);
  s->numFusions=1;
  HASH_ADD_KEYPTR( hh, myGeneFusionsCount, s->geneID, strlen(s->geneID), s );
}

struct geneFusionsCount *findGene( char* geneID ) 
{
  struct geneFusionsCount *s = NULL;
  HASH_FIND_STR( myGeneFusionsCount, geneID, s );
  return s;
}

int main (int argc, char *argv[])
{
  int i,j ;
  unsigned int countGFR, countRemoved=0;
  char *line;
  struct geneFusionsCount *currHash, *s, *tmp= NULL;
  Array gfrA = arrayCreate(20, GfrEntry); 

  // reading GFR file 
  gfr_init ( "-" );
  gfrA = gfr_parse();

  GfrEntry* gfrE;

  for( i=0; i<arrayMax(gfrA); i++) {
    gfrE = arrp( gfrA, i, GfrEntry );
    // adding/updating first transcript to the hash
    currHash = NULL;    
    currHash = findGene( gfrE->nameTranscript1 );
    if( currHash != NULL ) {
      // found; updating the count
      currHash->numFusions++;
    } else  {
      // not found; add element
      addGene( gfrE->nameTranscript1 );
    }

    // adding/updating second transcript to the hash
    currHash = NULL;    
    currHash = findGene( gfrE->nameTranscript2 );
    if( currHash != NULL ) {
      // found; updating the count
      currHash->numFusions++;
    } else  {
      // not found; add element
      addGene( gfrE->nameTranscript2 );
    }
  }
  printf("%s\n", gfr_writeHeader() );
  for( i=0; i<arrayMax(gfrA); i++) {
    GfrEntry* gfrE =  arrp( gfrA, i, GfrEntry );
    currHash = NULL;
    currHash = findGene( gfrE->nameTranscript1 );
    if( !currHash ) die("Error: transcript1 %s not found.", gfrE->nameTranscript1);
    if( currHash->numFusions < 4 ) {
      currHash = NULL;
      currHash = findGene( gfrE->nameTranscript2 );
      if( !currHash ) die("Error: transcript2 %s not found.", gfrE->nameTranscript2);
      if( currHash->numFusions < 4 ) {
	countGFR++;
	printf("%s\n" , gfr_writeGfrEntry(gfrE) );
      } else {
	countRemoved++;
      }
    }
  }
  gfr_deInit();
  HASH_ITER( hh, myGeneFusionsCount, s, tmp) {
    HASH_DEL( myGeneFusionsCount, s);
    free(s);
  }
  warn("%s_numRemoved: %d", argv[0], countRemoved);
  warn("%s_numGfrEntries: %d", argv[0], countGFR);
  arrayDestroy(gfrA);
  return 0;
}

