#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>

#include <bios/confp.h>
#include <bios/log.h>
#include <bios/format.h>
#include <bios/bowtieParser.h>
#include <bios/fasta.h>
#include <bios/intervalFind.h>
#include <bios/blastParser.h>
#include <bios/linestream.h>

#include "gfr.h"
#include "util.h"


/**
   @file gfrSmallScaleHomologyFilter.c 
   @brief  It removes mismapping artifacts.
   @details It removes candidates with not supported by uniquely mapped reads. It reads can be mapped to other regions in the genome, the candidate is removed. MAX_OVERLAP_ALLOWED will determine if homologous regions are found. If the  number of homologous reads is above a user defined threshold (MAX_FRACTION_HOMOLOGOUS), the fusion candidate is removed.
   @author Andrea Sboner  (andrea.sboner.w [at] gmail.com).  
   @version 0.8
   @date 2013.09.10
   @remarks WARNings will be output to stdout to summarize the filter results.
   @remarks To speed up the computation, gfServer and gfClient (part of Blat suite) are used. If the server is not running, it will be initiated. The location of the tools must be defined in .fusionseqrc. Moreover, the human genome reference is assigned to  BLAT_GFSERVER_PORT. Note that gfrRibosomal filter uses gfServer on BLAT_GFSERVER_PORT + 1.
   @pre It requires 'blat' and 'fastx_collapser' to be used, defined in .fusionseqrc
   @pre A valid GFR file as input, including stdin.
   @pre humanReference.2bit A 2bit file with the sequences of the human genome, defined in .fusionseqrc
 */
/**
  \file gfrSmallScaleHomologyFilter.c 
  \brief Removal of mismapping artifacts.
  \details Filter to remove artifacts due to mismapping for small scale homology within the two genes, including issues with splice junctions.
  \author Andrea Sboner (andrea.sboner.w [at] gmail.com).
  \version 0.8
  \date 2012.08.22						
  \pre It requires 'blat' and 'fastx' to be on the path.
 */

/*static int sortBowtieQueriesBySequenceName (BowtieQuery *a, BowtieQuery *b)
{
  return strcmp (a->sequenceName,b->sequenceName);
  }*/



static int sortFastaSequencesByName (Seq *a, Seq *b) 
{
  return strcmp (a->name,b->name);
}



static char getTranscriptNumber (char *seqName)
{
  char *pos;

  pos = strchr (seqName,'|');
  return *(pos + 1);
}

int main (int argc, char *argv[])
{
  GfrEntry *currGE;
  int i,j,k,l, h,index;
  Stringa buffer,cmd,fnSequencesToAlign;
  FILE *fp;
  FILE *fp1;
  FILE *fp2;
  FILE *freads1;
  FILE *freads2;
  Array gfrEntries;
  BowtieQuery *currBQ,testBQ;
  BowtieEntry *currBE;
  Texta seqNames;
  int readSize1, readSize2, minReadSize;
  Array bowtieQueries;
  char transcriptNumber;
  int isHomologous,homologousCount;
  int count;
  int countRemoved;
  unsigned short int tooMany;
  BlatQuery *blQ;

  config *conf;

  if ((conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL) {
    die("%s:\tCannot find .fusionseqrc", argv[0]);
    return EXIT_FAILURE;
  } 
  if ( (confp_get( conf, "BLAT_TWO_BIT_TO_FA")) == NULL) {
    die("%s:\tCannot find BLAT_TWO_BIT_TO_FA in the configuration file: %s", argv[0], getenv("FUSIONSEQ_CONFPATH") );
    return EXIT_FAILURE;
  } 
  if ( (confp_get( conf,"BLAT_DATA_DIR")) == NULL) {
    die("%s:\tCannot find BLAT_DATA_DIR in the configuration file: %sc", argv[0], getenv("FUSIONSEQ_CONFPATH") );
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
  cmd = stringCreate (100);
  // initializing the gfServers
  stringPrintf( cmd, "%s status %s %s &> /dev/null", confp_get( conf, "BLAT_GFSERVER"), confp_get( conf, "BLAT_GFSERVER_HOST"), confp_get( conf, "BLAT_GFSERVER_PORT") );
  int ret = hlr_system( string(cmd), 1 );
  if( ret != 0 ) { // not initialized
    stringPrintf( cmd , "%s -repMatch=100000 -tileSize=12 -canStop -log=%s/gfServer_genome.log start %s %s %s/%s  &", confp_get( conf, "BLAT_GFSERVER"), confp_get(conf, "TMP_DIR"),confp_get( conf, "BLAT_GFSERVER_HOST"), confp_get( conf, "BLAT_GFSERVER_PORT"), confp_get(conf, "BLAT_DATA_DIR"), confp_get(conf, "BLAT_TWO_BIT_DATA_FILENAME"));
    hlr_system( string( cmd ), 0 );
    long int startTime = time(0);
    stringPrintf( cmd , "%s status %s %s &2> /dev/null", confp_get( conf, "BLAT_GFSERVER"), confp_get( conf, "BLAT_GFSERVER_HOST"), confp_get( conf, "BLAT_GFSERVER_PORT"));
    while( hlr_system( string(cmd), 1) && (time(0)-startTime)<600 ) ;
    if( hlr_system( string(cmd), 1 ) != 0 )  {
      die("gfServer for %s/%s not initialized: %s %s %s", confp_get(conf, "BLAT_DATA_DIR"), confp_get(conf, "BLAT_TWO_BIT_DATA_FILENAME"), confp_get( conf, "BLAT_GFSERVER"), confp_get( conf, "BLAT_GFSERVER_HOST"), confp_get( conf, "BLAT_GFSERVER_PORT")); 
      return EXIT_FAILURE;
    }
  } 
  // end initialization

  gfr_init ("-");
  gfrEntries =  gfr_parse ();
  if (arrayMax (gfrEntries) == 0){
    puts (gfr_writeHeader ());
    gfr_deInit ();
    return 0;
  }
  seqNames = textCreate (10000); 
  buffer = stringCreate (100);
  fnSequencesToAlign = stringCreate (100);
  count = 0;
  countRemoved = 0;
  puts (gfr_writeHeader ());
 
  for (i = 0; i < arrayMax (gfrEntries); i++) {
    currGE = arrp (gfrEntries,i,GfrEntry);
    homologousCount = 0;
    minReadSize=10000;
    // creating two fasta files with the two genes
    
    stringPrintf( cmd, "%s %s/%s -seq=%s -start=%d -end=%d %s/%s_transcript1.fa", confp_get(conf, "BLAT_TWO_BIT_TO_FA") , confp_get(conf, "BLAT_DATA_DIR"), confp_get(conf, "BLAT_TWO_BIT_DATA_FILENAME"), currGE->chromosomeTranscript1, currGE->startTranscript1, currGE->endTranscript1, confp_get(conf, "TMP_DIR"), currGE->id);
    hlr_system( string(cmd) , 0);   
    stringPrintf( cmd, "%s %s/%s -seq=%s -start=%d -end=%d %s/%s_transcript2.fa", confp_get(conf, "BLAT_TWO_BIT_TO_FA"),  confp_get(conf, "BLAT_DATA_DIR"), confp_get(conf, "BLAT_TWO_BIT_DATA_FILENAME"), currGE->chromosomeTranscript2, currGE->startTranscript2, currGE->endTranscript2, confp_get(conf, "TMP_DIR"), currGE->id);
    hlr_system( string(cmd) , 0);   
    
    Stringa fa1 = stringCreate( 100 ); 
    Stringa fa2 = stringCreate( 100 );
    
    // creating the two fasta files with the reads
    stringPrintf( fa1, "%s/%s_reads1.fa", confp_get(conf, "TMP_DIR"), currGE->id);
    if (!(freads1 = fopen ( string(fa1) ,"w"))) {
      die ("Unable to open file: %s",string (fa1));
    }   
    // writing the reads of the first end into file
    
    for (l = 0; l < arrayMax (currGE->readsTranscript1); l++) {
      char* currRead1 = hlr_strdup( textItem (currGE->readsTranscript1,l)); // read1
      readSize1 = strlen( currRead1 );
      if( readSize1 == 0 ) die("Read size cannot be zero: read1[ %s ]", currRead1);
      if( readSize1 < minReadSize ) minReadSize = readSize1;
      fprintf( freads1, ">%d\n%s\n", l, currRead1 );
      hlr_free( currRead1 );
    }
    fclose( freads1 );
    
    stringPrintf( fa2, "%s/%s_reads2.fa", confp_get(conf, "TMP_DIR"), currGE->id);
    if (!(freads2 = fopen ( string(fa2) ,"w"))) {
      die ("Unable to open file: %s",string (fa2));
    } 
    // writing the reads of the second end into file
    for (l = 0; l < arrayMax (currGE->readsTranscript2); l++) {
      char* currRead2 = hlr_strdup( textItem (currGE->readsTranscript2,l)); // read2
      readSize2 = strlen( currRead2 );
      if( readSize2 == 0 ) die("Read size cannot be zero: read2[ %s ]", currRead2);
      if( readSize2 < minReadSize ) minReadSize = readSize2;
      fprintf( freads2, ">%d\n%s\n", l, currRead2 );
      hlr_free( currRead2 );
    }
    fclose( freads2 );      
    
    // collapse the reads 2  ## requires the FASTX package
    stringPrintf( cmd, "%s -i %s/%s_reads2.fa -o %s/%s_reads2.collapsed.fa", confp_get(conf, "FASTX_COLLAPSER"), confp_get(conf, "TMP_DIR"), currGE->id, confp_get(conf, "TMP_DIR"), currGE->id  );
    hlr_system (string (cmd),0);
    
    //blat of reads2 against the first transcript
    stringPrintf( cmd, "%s -t=dna -out=psl -fine -tileSize=15 %s/%s_transcript1.fa %s/%s_reads2.collapsed.fa stdout",confp_get(conf, "BLAT_BLAT"), confp_get(conf, "TMP_DIR"), currGE->id, confp_get(conf, "TMP_DIR"), currGE->id );
    
    // reading the results of blast from Pipe
    blatParser_initFromPipe( string(cmd) );
    while( blQ = blatParser_nextQuery() ) {
      int nucleotideOverlap = getNucleotideOverlap ( blQ );
      if ( nucleotideOverlap > ( ((double)readSize2)* atof(confp_get(conf,"MAX_OVERLAP_ALLOWED"))) ) {
	char* value = strchr(blQ->qName,'-');
	homologousCount+=atoi(value+1);
      }
    }
    blatParser_deInit();
    
    // collapse the reads 1 ## requires the FASTX package on the path
    stringPrintf( cmd, "%s -i %s/%s_reads1.fa -o %s/%s_reads1.collapsed.fa", confp_get(conf, "FASTX_COLLAPSER"), confp_get(conf, "TMP_DIR"), currGE->id, confp_get(conf, "TMP_DIR"), currGE->id  );
    hlr_system (string (cmd),0);
    
    //blat of reads1 against the second transcript
    stringPrintf( cmd, "%s -t=dna -out=psl -fine -tileSize=15 %s/%s_transcript2.fa %s/%s_reads1.collapsed.fa stdout",confp_get(conf, "BLAT_BLAT"), confp_get(conf, "TMP_DIR"), currGE->id, confp_get(conf, "TMP_DIR"), currGE->id  );
    
    blatParser_initFromPipe( string(cmd) );
    while( blQ = blatParser_nextQuery() ) {		
      int nucleotideOverlap = getNucleotideOverlap ( blQ );
      if ( nucleotideOverlap > ( ((double)readSize1)* atof(confp_get(conf,"MAX_OVERLAP_ALLOWED"))) ) {
	char* value = strchr(blQ->qName,'-');
	homologousCount+=atoi(value+1);
      }
    }
    blatParser_deInit();
    stringPrintf (cmd,"cd %s;rm -rf %s_reads?.fa %s_reads?.collapsed.fa %s_transcript?.fa", confp_get(conf, "TMP_DIR"), currGE->id,currGE->id,currGE->id);
    hlr_system( string(cmd) , 0);      
    if (((double)homologousCount / (double)arrayMax(currGE->readsTranscript1)) <= atof(confp_get(conf, "MAX_FRACTION_HOMOLOGOUS")) ) { 
      homologousCount = 0;
      // there is no homology between the two genes, but what about the rest of the genome
      writeFasta( currGE, &minReadSize,  confp_get(conf, "TMP_DIR") );
      stringPrintf(cmd, "cd %s; %s %s %s / -t=dna -q=dna -minScore=%d -out=psl %s_reads.fa %s.smallhomology.psl &>/dev/null", confp_get(conf, "TMP_DIR"), confp_get( conf, "BLAT_GFCLIENT"), confp_get( conf, "BLAT_GFSERVER_HOST"), confp_get( conf, "BLAT_GFSERVER_PORT"), minReadSize - (int)(0.1 * minReadSize) > 20 ? minReadSize - (int) (0.1 * minReadSize) : 20 ,  currGE->id,  currGE->id);
      int attempts=0;
      ret = hlr_system( string(cmd), 1 );
      while( hlr_system( string(cmd), 1 ) && attempts<5000 ) attempts++;
      if( attempts == 5000 ) {
	die("Cannot map the reads %s", string( cmd ));
	return EXIT_FAILURE;
      }
      // reading the results of blast from File
      stringPrintf(cmd,  "%s/%s.smallhomology.psl", confp_get( conf, "TMP_DIR"), currGE->id);
      blatParser_initFromFile( string(cmd) );
      tooMany = 1;
      while( blQ = blatParser_nextQuery() ) {
	tooMany = 0;
	if( arrayMax( blQ->entries ) > 1 ) {
	  homologousCount+= arrayMax( blQ->entries ) - 1;
	  char* value = strchr( blQ->qName,'/' );
	  if( value ) *value = '\0'; else die("Not a valid index in the blat query name:\t%s", blQ->qName );
	  int indexOfInter = atoi( blQ->qName ); // the following three lines should removed the read if writing the GFR entry
	  GfrInterRead *currGIR = arrp( currGE->interReads, indexOfInter, GfrInterRead );
	  currGIR->flag = 1;
	}
      }
      blatParser_deInit();
      if (  tooMany == 1 || ( ( (double) homologousCount / (double) ( arrayMax(currGE->readsTranscript1) + arrayMax(currGE->readsTranscript2) ) )  > atof(confp_get(conf, "MAX_FRACTION_HOMOLOGOUS")) ) ) {
	countRemoved++;
	stringPrintf (cmd,"cd %s; rm -rf %s_reads*.fa %s_reads?.collapsed.fa %s_transcript?.fa %s.smallhomology.psl", confp_get(conf, "TMP_DIR"), currGE->id,currGE->id,currGE->id,currGE->id);
	hlr_system( string(cmd), 1 );
	continue;
      }
      // writing the gfrEntry, if everthing else didn't stop 
      if( homologousCount > 0 ) updateStats( currGE );
      puts (gfr_writeGfrEntry (currGE));
      count++;
      // removing temporary files
      stringPrintf (cmd,"cd %s;rm -rf %s_reads*.fa %s_reads?.collapsed.fa %s_transcript?.fa  %s.smallhomology.psl", confp_get(conf, "TMP_DIR"), currGE->id,currGE->id,currGE->id,currGE->id);
      hlr_system( string(cmd) , 1);      
    } else {
      countRemoved++;
    }
    
  }

  gfr_deInit ();

  stringDestroy (fnSequencesToAlign);
  stringDestroy (cmd);
  stringDestroy (buffer);
  warn ("%s_numRemoved: %d",argv[0],countRemoved);  
  warn ("%s_numGfrEntries: %d",argv[0],count);

  confp_close(conf);

  return EXIT_SUCCESS;
}


