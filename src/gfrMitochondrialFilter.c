#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/blatParser.h>
#include <bios/linestream.h>
#include "gfr.h"


void writeFasta( GfrEntry* currGE, unsigned int *minReadSize ) {
  FILE *freads;
  int l, readSize1, readSize2;
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
	arrayMax( currGE->readsTranscript1),
	arrayMax( currGE->readsTranscript2) );
  
  // writing the reads into file
  for (l = 0; l < arrayMax (currGE->readsTranscript1); l++) {      
    char* currRead1 = hlr_strdup( textItem (currGE->readsTranscript1,l)); // read1
    char* currRead2 = hlr_strdup( textItem (currGE->readsTranscript2,l)); // read2
    fprintf( freads, ">%d/1\n%s\n>%d/2\n%s\n", l+1, currRead1, l+1, currRead2 );
    readSize1 = strlen( currRead1 );
    readSize2 = strlen( currRead2 );
    
    if(readSize1 != readSize2 )
      die("The two reads have different lengths: 1:%d vs 2:%d", readSize1, readSize2);
    if ( readSize1 < *minReadSize )
      *minReadSize = readSize1 ;
    hlr_free( currRead1 );
    hlr_free( currRead2 );
  }
  fclose( freads ); 
  freads=NULL;
  stringDestroy( readsFA );
}
  
int main (int argc, char *argv[])
{
  GfrEntry *currGE;
  int count;
  int countRemoved;
  int mitochondrialCount; 
  unsigned int minReadSize;
  int  i;
  Stringa cmd;
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
 if( confp_get( conf, "MITOCHONDRIAL_DIR")==NULL ) {
    die("%s:\tCannot find MITOCHONDRIAL_DIR in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
    return EXIT_FAILURE;
  }
  if( confp_get( conf, "MITOCHONDRIAL_FILENAME")==NULL ) {
    die("%s:\tCannot find MITOCHONDRIAL_FILENAME in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
    return EXIT_FAILURE;
  }
if( confp_get( conf, "TMP_DIR")==NULL ) {
    die("%s:\tCannot find TMP_DIR in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
    return EXIT_FAILURE;
  }
  count = 0;
  countRemoved = 0;
  
  cmd = stringCreate (100);
  // initializing the gfServers
  stringPrintf( cmd, "gfServer status localhost 8081 2> /dev/null" );
  LineStream ls = ls_createFromPipe( string(cmd) );
  if( ls_nextLine( ls ) == NULL  ) { // not initialized
    ls_destroy_func( ls );
    stringPrintf( cmd , "gfServer -repMatch=100000 -tileSize=12 -canStop -log=%s/gfServer_mitochondrial.log start localhost 8081 %s/%s  &", confp_get(conf, "TMP_DIR"), confp_get(conf, "MITOCHONDRIAL_DIR"), confp_get(conf,"MITOCHONDRIAL_FILENAME"));
    hlr_system( string( cmd ), 0 );
    long int startTime = time(0);
    stringPrintf( cmd , "gfServer status localhost 8081 2> /dev/null");
    int initialized=0;
    while( !initialized && (time(0)-startTime)<600 ) {
      ls = ls_createFromPipe( string(cmd) );
      if( ls_nextLine( ls ) != NULL ) initialized=1; 
      ls_destroy_func( ls );
    }
    if( initialized==0 )  {
      die("gfServer not initialized"); 
      return EXIT_FAILURE;
    }
  } 

 
  
  gfr_init ("-");
  puts (gfr_writeHeader ());
  while (currGE = gfr_nextEntry ()) {
    if (strEqual(currGE->chromosomeTranscript1, "chrM") || 
	strEqual(currGE->chromosomeTranscript2, "chrM")) {
      countRemoved++;
      continue;
    } else {
      mitochondrialCount = 0;
      minReadSize=1000;
      writeFasta( currGE, &minReadSize );
      stringPrintf(cmd, "gfClient localhost 8081 / -t=dna -q=dna -minScore=%d -out=psl %s_reads.fa stdout" , minReadSize - 10 > 20 ? minReadSize - 10 : 20 , currGE->id);
      // reading the results of blast from Pipe
      blatParser_initFromPipe( string(cmd) );
      while( blQ = blatParser_nextQuery() ) {
	int nucleotideOverlap = getNucleotideOverlap ( blQ );
	if (nucleotideOverlap > (((double)minReadSize) * strtod(confp_get(conf, "MAX_OVERLAP_ALLOWED"), NULL))) {
	  mitochondrialCount++;
	} 
      }
      blatParser_deInit();
      if (( (double)mitochondrialCount / ( (double)currGE->numInter * 2.0 ) ) <= strtod(confp_get(conf, "MAX_FRACTION_HOMOLOGOUS"), NULL)) {       
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
    
  }
  gfr_deInit ();
 
  stringDestroy( cmd );
  warn ("%s_numRemoved: %d",argv[0],countRemoved);
  warn ("%s_numGfrEntries: %d",argv[0],count);
  confp_close(conf);
  return 0;
}

 /*

 // stopping the gfServer
  //stringPrintf( cmd, "gfServer stop localhost %d", port);
  //hlr_system( string(cmd), 0 );
  //stringPrintf( cmd, "rm -rf gfServer_%d_mitochondrial.log", port);
  //hlr_system( string(cmd), 0 );
srand (  (unsigned int) getpid() ) ;
  port = 8080 + (int) (1000.0 * (random() / (RAND_MAX + 1.0 ) ) );
  stringPrintf( cmd , "gfServer -repMatch=100000 -tileSize=12 -canStop -log=gfServer_%d_mitochondrial.log start localhost %d %s/%s  &", port, port, "/home/asboner/FusionSeqData/human/hg19", "chrM.2bit");
  hlr_system( string( cmd ), 1 );*/
