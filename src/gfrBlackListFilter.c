#include <stdio.h>

#include <bios/log.h>
#include <bios/format.h>
#include <bios/linestream.h>
#include <bios/confp.h>

#include "gfr.h"

/**
   @file gfrBlackListFilter.c
   @brief It removes candidates specified by the user in a file
   @details It removes candidates specified by the user in a file.
   
   @author Andrea Sboner  (andrea.sboner.w [at] gmail.com).  
   @version 0.8
   @date 2013.09.10
   @remarks WARNings will be output to stdout to summarize the filter results.
   @pre A valid GFR file as input, including stdin.
   @pre blacklist a tab delimited file with the two gene symbols to removed; defined in .fusionseqrc 
 */


typedef struct {
  char* gene1;
  char* gene2;
} BLEntry;

static int sortBlackListByName1 (BLEntry *a, BLEntry *b) 
{
  int res = strcmp ( a->gene1, b->gene1);
  if( res==0 ) res = strcmp ( a->gene2, b->gene2 );
  return (res); //(strcmp ( a->gene1, b->gene1));
}

int main (int argc, char *argv[])
{
  GfrEntry *currGE;
  BLEntry *currBLE;
  BLEntry currQuery;
  FILE *fp;
  char *line;
  int count;
  int countRemoved;
  
  int index;
  WordIter w;
  Array blackList = arrayCreate(20, BLEntry);
  config *Conf;

  if ((Conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL) {
    die("%s:\tCannot find .fusionseqrc: %s", argv[0], getenv("FUSIONSEQ_CONFPATH"));
    return EXIT_FAILURE;
  }
  if( confp_get( Conf, "ANNOTATION_DIR")==NULL ) {
    die("%s:\tCannot find ANNOTATION_DIR in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
    return EXIT_FAILURE;
  }
  if( confp_get( Conf, "BLACKLIST_FILENAME")==NULL ) {
    die("%s:\tCannot find BLACKLIST_FILENAME in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
    return EXIT_FAILURE;
  }
  Stringa buffer=stringCreate( 100 );
  stringPrintf( buffer, "%s/%s", confp_get( Conf, "ANNOTATION_DIR"), confp_get( Conf, "BLACKLIST_FILENAME") );
  fp = fopen( string( buffer ), "r" );
  stringDestroy( buffer );
  
  if( !fp )  die("Unable to open file: %s", confp_get( Conf, "BLACKLIST_FILENAME"));
  // reading blacklist file
  LineStream ls = ls_createFromFile( argv[1] );
  while( line = ls_nextLine(ls) ) {
    w = wordIterCreate( line, "\t", 1);
    currBLE = arrayp( blackList, arrayMax(blackList), BLEntry);
    currBLE->gene1 = hlr_strdup ( wordNext(w) );
    currBLE->gene2 = hlr_strdup ( wordNext(w) );    
    wordIterDestroy(w);
  }
  fclose(fp);
  arraySort( blackList, (ARRAYORDERF) sortBlackListByName1);

  // beginFiltering
  count = 0;
  countRemoved = 0;
  gfr_init ("-");
  puts (gfr_writeHeader ());
  while (currGE = gfr_nextEntry ()) { // reading the gfr
    // creating a new query to the black list
    currQuery.gene1 = currGE->geneSymbolTranscript1;
    currQuery.gene2 = currGE->geneSymbolTranscript2;
    if( strEqual( currQuery.gene1 , currQuery.gene2 ) ) {
	countRemoved++;
	continue;
      }
    // searching against read_1/read_2
    int res = arrayFind( blackList, &currQuery, 
			 &index,  (ARRAYORDERF) sortBlackListByName1);  
    
    if( !res ) { // not found, then searching against read_2/read_1
      currQuery.gene1 = currGE->geneSymbolTranscript2;
      currQuery.gene2 = currGE->geneSymbolTranscript1;
      
      res =  arrayFind( blackList, &currQuery, 
			&index, (ARRAYORDERF) sortBlackListByName1 );
      
      if( !res ) { // not found, write the instance to stdout, update the counts
	puts (gfr_writeGfrEntry (currGE));
	count++;	
      } else { // found: read2/read1
	countRemoved++;
      }	
    } else { //found: read1/read2
      countRemoved++;
    }
  }	           
  gfr_deInit ();
  arrayDestroy( blackList );
  warn ("%s_BlackListFilter: %s",argv[0], confp_get( Conf, "BLACKLIST_FILENAME"));
  warn ("%s_numRemoved: %d",argv[0],countRemoved);
  warn ("%s_numGfrEntries: %d",argv[0],count);
  confp_close( Conf);
  return 0;
}

