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

#include "gfr.h"
#include "util.h"

/**
  \file gfrSmallScaleHomologyFilter.c 
  \brief Removal of mismapping artifacts.
  \details Filter to remove artifacts due to mismapping for small scale homology within the two genes, including issues with splice junctions.
  \author Andrea Sboner (andrea.sboner.w [at] gmail.com).
  \version 0.8
  \date 2012.08.22						
  \pre It requires 'blat' and 'fastx' to be on the path.
 */

static int sortBowtieQueriesBySequenceName (BowtieQuery *a, BowtieQuery *b)
{
  return strcmp (a->sequenceName,b->sequenceName);
}



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
  int i,j,k,l,index;
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
  int readSize1, readSize2;
  Array bowtieQueries;
  char transcriptNumber;
  int isHomologous,homologousCount;
  int count;
  int countRemoved;
  BlatQuery *blQ;

  config *conf;

  if ((conf = confp_open(getenv("FUSIONSEQ_CONFPATH"))) == NULL)
    return EXIT_FAILURE;
  gfr_init ("-");
  gfrEntries =  gfr_parse ();
  if (arrayMax (gfrEntries) == 0){
    puts (gfr_writeHeader ());
    gfr_deInit ();
    return 0;
  }
  seqNames = textCreate (10000); 
  buffer = stringCreate (100);
  cmd = stringCreate (100);
  fnSequencesToAlign = stringCreate (100);
  count = 0;
  countRemoved = 0;
  puts (gfr_writeHeader ());
  j = 0;
  for (i = 0; i < arrayMax (gfrEntries); i++) {
    currGE = arrp (gfrEntries,i,GfrEntry);
    homologousCount = 0;
  
    if (((double)homologousCount / arrayMax(currGE->readsTranscript1)) <= atof(confp_get(conf, "MAX_FRACTION_HOMOLOGOUS")) ) { 
      // creating two fasta files with the two genes
      //warn("%d %s", i, currGE->id);
      stringPrintf( cmd, "%s %s/%s -seq=%s -start=%d -end=%d %s_transcript1.fa", confp_get(conf, "BLAT_TWO_BIT_TO_FA") , confp_get(conf, "BLAT_DATA_DIR"), confp_get(conf, "BLAT_TWO_BIT_DATA_FILENAME"), currGE->chromosomeTranscript1, currGE->startTranscript1, currGE->endTranscript1, currGE->id);
      hlr_system( string(cmd) , 0);   
      stringPrintf( cmd, "%s %s/%s -seq=%s -start=%d -end=%d %s_transcript2.fa", confp_get(conf, "BLAT_TWO_BIT_TO_FA"),  confp_get(conf, "BLAT_DATA_DIR"), confp_get(conf, "BLAT_TWO_BIT_DATA_FILENAME"), currGE->chromosomeTranscript2, currGE->startTranscript2, currGE->endTranscript2, currGE->id);
      hlr_system( string(cmd) , 0);   

      Stringa fa1 = stringCreate( 100 ); 
      Stringa fa2 = stringCreate( 100 );

      // creating the two fasta files with the reads
      stringPrintf( fa1, "%s_reads1.fa", currGE->id);
      if (!(freads1 = fopen ( string(fa1) ,"w"))) {
	die ("Unable to open file: %s",string (fa1));
      }   
      // writing the reads of the first end into file
      for (l = 0; l < arrayMax (currGE->readsTranscript1); l++) {
	char* currRead1 = hlr_strdup( textItem (currGE->readsTranscript1,l)); // read1
	readSize1 = strlen( currRead1 );
	fprintf( freads1, ">%d\n%s\n", l, currRead1 );
	hlr_free( currRead1 );
      }
      fclose( freads1 );
      
      stringPrintf( fa2, "%s_reads2.fa", currGE->id);
      if (!(freads2 = fopen ( string(fa2) ,"w"))) {
	die ("Unable to open file: %s",string (fa2));
      } 
      // writing the reads of the second end into file
      for (l = 0; l < arrayMax (currGE->readsTranscript2); l++) {
	char* currRead2 = hlr_strdup( textItem (currGE->readsTranscript2,l)); // read2
	readSize2 = strlen( currRead2 );
	fprintf( freads2, ">%d\n%s\n", l, currRead2 );
	hlr_free( currRead2 );
      }
      fclose( freads2 );      
      
      // collapse the reads 2  ## requires the FASTX package on the path
      stringPrintf( cmd, "fastx_collapser -i %s_reads2.fa -o %s_reads2.collapsed.fa", currGE->id, currGE->id );
	hlr_system (string (cmd),0);

      //blat of reads2 against the first transcript
      stringPrintf( cmd, "blat -t=dna -out=psl -fine -tileSize=15 %s_transcript1.fa %s_reads2.collapsed.fa stdout", currGE->id, currGE->id );
  
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
      stringPrintf( cmd, "fastx_collapser -i %s_reads1.fa -o %s_reads1.collapsed.fa", currGE->id, currGE->id  );
      hlr_system (string (cmd),0);

      //blat of reads1 against the second transcript
      stringPrintf( cmd, "blat -t=dna -out=psl -fine -tileSize=15 %s_transcript2.fa %s_reads1.collapsed.fa stdout", currGE->id, currGE->id  );

      blatParser_initFromPipe( string(cmd) );
      while( blQ = blatParser_nextQuery() ) {		
	int nucleotideOverlap = getNucleotideOverlap ( blQ );
	if ( nucleotideOverlap > ( ((double)readSize1)* atof(confp_get(conf,"MAX_OVERLAP_ALLOWED"))) ) {
	  char* value = strchr(blQ->qName,'-');
	  homologousCount+=atoi(value+1);
	}
      }
      blatParser_deInit();

      if (((double)homologousCount / (double)arrayMax(currGE->readsTranscript1)) <= atof(confp_get(conf, "MAX_FRACTION_HOMOLOGOUS")) ) { 
	// writing the gfrEntry
	puts (gfr_writeGfrEntry (currGE));
	count++;
      } else {
	countRemoved++;
      }
      // removing temporary files
      stringPrintf (cmd,"rm -rf %s_reads?.fa %s_reads?.collapsed.fa %s_transcript?.fa", currGE->id,currGE->id,currGE->id);
      hlr_system( string(cmd) , 0);      
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

