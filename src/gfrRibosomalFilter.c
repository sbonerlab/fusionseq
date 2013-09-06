#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/bowtieParser.h>
#include <bios/fasta.h>
#include <bios/intervalFind.h>
#include <bios/blatParser.h>
#include <bios/linestream.h>

#include "gfr.h"
#include "util.h"

int main (int argc, char *argv[])
{
  GfrEntry *currGE;
  int i,j,l;
  Stringa cmd;
  FILE *freads;
  Array gfrEntries;
  int ribosomalCount;
  int count;
  int countRemoved;
  int readSize1,readSize2,minReadSize;
  BlatQuery *blQ=NULL;
  
  config *conf;
  
  if ((conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL) {
    die("%s:\tCannot find .fusionseqrc", argv[0]);
    return EXIT_FAILURE;
  }
  if( confp_get( conf,"MAX_OVERLAP_ALLOWED")==NULL ) {
    die("%s:\tCannot find MAX_OVERLAP_ALLOWED in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
    return EXIT_FAILURE;
  }
  if( confp_get( conf,"MAX_FRACTION_HOMOLOGOUS")==NULL ) {
    die("%s:\tCannot find MAX_FRACTION_HOMOLOGOUS in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
    return EXIT_FAILURE;
  }
  if( confp_get( conf, "RIBOSOMAL_DIR")==NULL ) {
    die("%s:\tCannot find RIBOSOMAL_DIR in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
    return EXIT_FAILURE;
  }
  if( confp_get( conf, "RIBOSOMAL_FILENAME")==NULL ) {
    die("%s:\tCannot find RIBOSOMAL_FILENAME in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
    return EXIT_FAILURE;
  }
  if( confp_get( conf, "TMP_DIR")==NULL ) {
    die("%s:\tCannot find TMP_DIR in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
    return EXIT_FAILURE;
  }
  cmd = stringCreate (100);
  // initializing the gfServers
  stringPrintf( cmd, "gfServer status localhost 8080 2> /dev/null" );
  LineStream ls = ls_createFromPipe( string(cmd) );
  if( ls_nextLine( ls ) == NULL  ) { // not initialized
    ls_destroy_func( ls );
    stringPrintf( cmd , "gfServer -repMatch=100000 -tileSize=12 -canStop -log=%s/gfServer_ribosomal.log start localhost 8080 %s/%s &", confp_get(conf, "TMP_DIR"), confp_get(conf, "RIBOSOMAL_DIR"), confp_get(conf, "RIBOSOMAL_FILENAME"));
    hlr_system( string( cmd ), 0 );
    
    long int startTime = time(0);
    stringPrintf( cmd , "gfServer status localhost 8080 2> /dev/null");
    static unsigned short int initialized=0;
    while( !initialized && (time(0)-startTime)<600 ) {
      ls = ls_createFromPipe( string(cmd) );
      if( ls_nextLine( ls ) != NULL ) initialized=1; 
      ls_destroy_func( ls );
    }
    if( initialized==0 )  die("gfServer not initialized");
  } 
  gfr_init ("-");
  gfrEntries = arrayCreate( 100, GfrEntry );
  gfrEntries =  gfr_parse ();
  if (arrayMax (gfrEntries) == 0){
    puts (gfr_writeHeader ());
    gfr_deInit ();
    return 0;
  }
  
  
  count = 0;
  countRemoved = 0;
  puts (gfr_writeHeader ());
  j = 0;
  for (i = 0; i < arrayMax (gfrEntries); i++) {
    currGE = arrp (gfrEntries,i,GfrEntry);
    ribosomalCount = 0;
    
    // creating one fasta files with the reads
    Stringa readsFA = stringCreate( 100 ); 
    
    // creating the fasta files with the reads 
    stringPrintf( readsFA, "%s_reads.fa", currGE->id);
    freads = fopen ( string(readsFA) ,"w");
    if (freads == NULL) {
      die ("Unable to open file: %s",string (readsFA));
    }     
    if (arrayMax(currGE->readsTranscript1) != arrayMax(currGE->readsTranscript2))
      die("Error: different number of inter-transcript reads %d vs. %d", 
          arrayMax(currGE->readsTranscript1),
	  arrayMax( currGE->readsTranscript2));
    
    // writing the reads into file
    for (l = 0; l < arrayMax (currGE->readsTranscript1); l++) {      
      char* currRead1 = hlr_strdup( textItem (currGE->readsTranscript1,l)); // read1
      char* currRead2 = hlr_strdup( textItem (currGE->readsTranscript2,l)); // read2
      fprintf( freads, ">%d/1\n%s\n>%d/2\n%s\n", l+1, currRead1, l+1, currRead2 );
      readSize1 = strlen( currRead1 );
      readSize2 = strlen( currRead2 );
      
      if(readSize1 != readSize2 )
        die("The two reads have different lengths: 1:%d vs 2:%d", readSize1, readSize2);
      if( i == 0 )
	minReadSize = readSize1 ;
      else {
	if ( readSize1 < minReadSize )
	  minReadSize = readSize1 ;
      }
      hlr_free( currRead1 );
      hlr_free( currRead2 );
    }
    fclose( freads ); 
    freads=NULL;
    stringDestroy( readsFA );
    /*
    //blat of reads against the ribosomal genes
    // 2*stepSize + tileSize - 1 = min_num_nt_to_trigger_alignment
    // tileSize = 15
    // int stepSize = (int)((readSize1 * confp_get(conf, "MAX_OVERLAP_ALLOWED") + 1 + 11) * 0.5);
    stringPrintf(cmd, "blat -t=dna -q=dna -out=psl -fine -repMatch=1000000 -tileSize=15 %s/%s %s_reads.fa stdout", 
    confp_get(conf, "RIBOSOMAL_DIR"), 
    confp_get(conf, "RIBOSOMAL_FILENAME"), 
    currGE->id);
    
    */
    
    stringPrintf(cmd, "gfClient localhost 8080 / -t=dna -q=dna -minScore=%d -out=psl %s_reads.fa stdout" , minReadSize - 10 > 20 ? minReadSize - 10 : 20 , currGE->id);
    // reading the results of blast from Pipe
    blatParser_initFromPipe( string(cmd) );
    while( blQ = blatParser_nextQuery() ) {
      int nucleotideOverlap = getNucleotideOverlap ( blQ );
      if (nucleotideOverlap > (((double)readSize1) * strtod(confp_get(conf, "MAX_OVERLAP_ALLOWED"), NULL))) {
	ribosomalCount++;
      } 
    }
    blatParser_deInit();
    if (( (double)ribosomalCount / ( (double)currGE->numInter * 2.0 ) ) <= strtod(confp_get(conf, "MAX_FRACTION_HOMOLOGOUS"), NULL)) {       
      // writing the gfrEntry
      puts (gfr_writeGfrEntry (currGE));
      count++;
    } else {
      countRemoved++;
    }
    // removing temporary files
    stringPrintf (cmd,"rm -rf %s_reads.fa", currGE->id );
    hlr_system( string(cmd) , 0);      
  }
  
  gfr_deInit ();
  arrayDestroy ( gfrEntries );
  stringDestroy( cmd );
  warn ("%s_numRemoved: %d",argv[0],countRemoved);  
  warn ("%s_numGfrEntries: %d",argv[0],count);
 
  confp_close(conf);
  return 0;
}

