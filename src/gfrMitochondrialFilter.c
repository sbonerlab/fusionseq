#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/blatParser.h>
//#include <bios/linestream.h>
#include <bios/confp.h>

#include "util.h"
#include "gfr.h"

/**
   @file gfrMitochondrialFilter.c
   @brief Filter to remove candidates involving mitochondrial genes.
   @details It removes candidates involving mitochondrial genes and also those with reads overlapping mitochondrial sequences.The MAX_OVERLAP_ALLOWED  determines if a read is sufficiently similar and MAX_FRACTION_HOMOLOGOUS determines the max number of reads that could have some homology to the mitochndrial chromosome
   
   @author Andrea Sboner  (andrea.sboner.w [at] gmail.com).  
   @version 0.8
   @date 2013.09.10

  @attention It requires a GFR file from stdin: @code $ gfrMitochondrial < file.gfr @endcode
 
   @remarks WARNings will be output to stdout to summarize the filter results.
   @remarks To speed up the computation, gfServer and gfClient (part of Blat suite) are used. If the server is not running, it will be initiated. The location of the tools must be defined in .fusionseqrc. Moreover, the mitochondrial reference is assigned to BLAT_GFSERVER_PORT + 2. 
 */
  
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
  config *conf = NULL; /**< Pointer to configuration file .fusionseqrc  */

  if ((conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL) {
    die("%s:\tCannot find .fusionseqrc: %s", argv[0], getenv("FUSIONSEQ_CONFPATH"));
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
 if( confp_get( conf, "BLAT_GFSERVER")==NULL ) {
    die("%s:\tCannot find BLAT_GFSERVER in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
    return EXIT_FAILURE;
  }
 if( confp_get( conf, "BLAT_GFCLIENT")==NULL ) {
    die("%s:\tCannot find BLAT_GFCLIENT in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
    return EXIT_FAILURE;
  }
if( confp_get( conf, "BLAT_GFSERVER_HOST")==NULL ) {
    die("%s:\tCannot find BLAT_GFSERVER_HOST in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
    return EXIT_FAILURE;
  }if( confp_get( conf, "BLAT_GFSERVER_PORT")==NULL ) {
    die("%s:\tCannot find BLAT_GFSERVER_PORT in the configuration file: %s)", argv[0], getenv("FUSIONSEQ_CONFPATH") );
    return EXIT_FAILURE;
  }

  count = 0;
  countRemoved = 0;
  
  cmd = stringCreate (100);
  // initializing the gfServers
  stringPrintf( cmd, "%s status %s %d &> /dev/null", confp_get( conf, "BLAT_GFSERVER"),  confp_get( conf, "BLAT_GFSERVER_HOST"), atoi(confp_get( conf, "BLAT_GFSERVER_PORT")) + 2);
  int ret = hlr_system( string(cmd), 1 );
   if( ret != 0 ) { // not initialized
    stringPrintf( cmd , "%s -repMatch=100000 -tileSize=12 -canStop -log=%s/gfServer_mitochondrial.log start %s %d %s/%s  &", confp_get( conf, "BLAT_GFSERVER"), confp_get( conf, "TMP_DIR"),  confp_get( conf, "BLAT_GFSERVER_HOST"), atoi(confp_get( conf, "BLAT_GFSERVER_PORT")) + 2, confp_get( conf, "MITOCHONDRIAL_DIR"), confp_get( conf,"MITOCHONDRIAL_FILENAME"));
    hlr_system( string( cmd ), 0 );
    long int startTime = time(0);
    stringPrintf( cmd , "%s status %s %d &> /dev/null", confp_get( conf, "BLAT_GFSERVER"), confp_get( conf, "BLAT_GFSERVER_HOST"), atoi(confp_get( conf, "BLAT_GFSERVER_PORT")) + 2);
    while( hlr_system( string(cmd), 1) && (time(0)-startTime)<600 ) ;
    if( hlr_system( string(cmd), 1 ) != 0 )  {
      die("gfServer for %s/%s not initialized: %s %s %s", confp_get( conf, "MITOCHONDRIAL_DIR"), confp_get( conf, "MITOCHONDRIAL_FILENAME"), confp_get( conf, "BLAT_GFSERVER"), confp_get( conf, "BLAT_GFSERVER_HOST"), confp_get( conf, "BLAT_GFSERVER_PORT")); 
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
      writeFasta( currGE, &minReadSize, confp_get( conf, "TMP_DIR") ); // in util.c
      stringPrintf(cmd, "cd %s;%s %s %d / -t=dna -q=dna -minScore=%d -out=psl %s_reads.fa %s.mito.psl &>/dev/null", confp_get( conf, "TMP_DIR"), confp_get( conf, "BLAT_GFCLIENT"), confp_get( conf, "BLAT_GFSERVER_HOST"), atoi(confp_get( conf, "BLAT_GFSERVER_PORT")) + 2, minReadSize - 5 > 20 ? minReadSize - 5 : 20 , currGE->id, currGE->id);
      int attempts=0;
      ret = hlr_system( string(cmd), 1 );
      while( hlr_system( string(cmd), 1 ) && attempts<50 ) attempts++;
      if( attempts == 50 ) {
	die("Cannot map the reads %s", string( cmd ));
	return EXIT_FAILURE;
      }

      // reading the results of blast from File
      stringPrintf(cmd,  "%s/%s.mito.psl", confp_get( conf, "TMP_DIR"), currGE->id);
      blatParser_initFromFile( string(cmd) );
      while( blQ = blatParser_nextQuery() ) {
	//warn("iter %d\tquery %s", iter, blQ->qName );iter++; 
	int nucleotideOverlap = getNucleotideOverlap ( blQ );
	if (nucleotideOverlap > (((double) minReadSize) * strtod(confp_get( conf, "MAX_OVERLAP_ALLOWED"), NULL))) {
	  char* value = strchr( blQ->qName,'/' );
	  if( value ) *value = '\0'; else die("Not a valid index in the blat query name:\t%s", blQ->qName );
	  int indexOfInter = atoi( blQ->qName ); // the following three lines should removed the read if writing the GFR entry
	  GfrInterRead *currGIR = arrp( currGE->interReads, indexOfInter, GfrInterRead );
	  currGIR->flag = 1;
	  mitochondrialCount++;
	} 
      }
      blatParser_deInit();
      if ( ( (double) mitochondrialCount / (double) ( arrayMax(currGE->readsTranscript1) + arrayMax(currGE->readsTranscript2) ) ) <= strtod(confp_get( conf, "MAX_FRACTION_HOMOLOGOUS"), NULL)) {   
	if( mitochondrialCount > 0 ) updateStats( currGE );
	// writing the gfrEntry
	puts (gfr_writeGfrEntry (currGE));
	count++;
      } else {
	countRemoved++;
      }
      // removing temporary files
      stringPrintf (cmd,"rm -rf %s/%s_reads.fa %s/%s.mito.psl", confp_get( conf, "TMP_DIR"),  currGE->id, confp_get( conf, "TMP_DIR"),  currGE->id );
      hlr_system( string(cmd) , 1);      
    } 
    
  }
  gfr_deInit ();
 
  stringDestroy( cmd );
  warn ("%s_numRemoved: %d",argv[0],countRemoved);
  warn ("%s_numGfrEntries: %d",argv[0],count);
  confp_close(conf);
  return 0;
}


/*int indexOfInter= atoi( *(strchr(blQ->qName,'/'))='\0' ); // this is to remove the reads, but we should keep it 
	  if( arrayRemoveD( currGE->interReads, indexOfInter ) != 0 ) 
	    die("Something went wrong when removing element[%d] from the array interReads(%d)", indexOfInter, arrayMax(currGE->interReads));
	  if( arrayRemoveD( currGE->readsTranscript1, indexOfInter ) != 0 ) 
	    die("Something went wrong when removing element[%d] from the array readsTranscript1(%d)", indexOfInter, arrayMax(currGE->readsTranscript1));
	  if( arrayRemoveD( currGE->readsTranscript2, indexOfInter ) != 0 ) 
	  die("Something went wrong when removing element[%d] from the array readsTranscript2(%d)", indexOfInter, arrayMax(currGE->readsTranscript2));*/
