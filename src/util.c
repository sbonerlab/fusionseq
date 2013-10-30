#include <bios/log.h>
#include <bios/format.h>
#include <bios/linestream.h>
#include <bios/bits.h>

#include "util.h"
#include "gfr.h"

int getNucleotideOverlap ( BlatQuery* blQ ) {
  int l;
  PslEntry* blE=NULL;     
  int qSize = arrp( blQ->entries, arrayMax(blQ->entries)-1, PslEntry)->qSize;
  int maxOverlap=0;
  //Bits* bitOverlap = bitAlloc( qSize  );	
  //bitClear( bitOverlap, qSize );
  //warn( "%d", arrayMax( blQ->entries ) );
  for( l=0; l < arrayMax ( blQ->entries ); l++ ) {
    blE = arrp( blQ->entries, l, PslEntry );
    if( blE->qSize != qSize ) die("Query size different from PslEntry query size: qSize(%d) - blE->qSize(%d)", qSize, blE->qSize);    
    if( ((blE->qEnd - blE->qStart)+1) > maxOverlap ) maxOverlap = (blE->qEnd - blE->qStart)+1;
    /*    if( (blE->qEnd > blE->qStart)  ) { // to ensure that at least two nucleotides matches and only significant hits are considered // 
      Bits* currBitEntry = bitAlloc( blE->qSize );
      bitClear ( currBitEntry, blE->qSize );
      if( (blE->qEnd - blE->qStart)> blE->qSize ) die("The query match is bigger than the read size: name(%s) - qSize(%d) - actualSize(%d)", blQ->qName, blE->qSize, (blE->qEnd - blE->qStart));
      bitSetRange( currBitEntry, (blE->qStart - 1), (blE->qEnd - blE->qStart));
      bitOr( bitOverlap, currBitEntry, qSize );
      bitFree( &currBitEntry );
    }//*/ 
  }
  //int overlap = 50; //bitCountRange( bitOverlap, 0, readSize);
  //bitFree( &bitOverlap);
  return maxOverlap;
}



void updateStats( GfrEntry* currGE ) {
  int i;
  float numInters=0.0;
  for( i = 0; i < arrayMax( currGE->interReads ); i++ ) {
    GfrInterRead* currGIR = arrp( currGE->interReads, i, GfrInterRead );
    if( currGIR->flag ) continue;
    char* read1 = arru( currGE->readsTranscript1, i, char* );
    char* read2 = arru( currGE->readsTranscript2, i, char* );
    if( read1==NULL || read2==NULL || (currGIR->readEnd1 - currGIR->readStart1 + 1) == 0 | 
	(currGIR->readEnd2 - currGIR->readStart2 + 1) == 0 ) 
      die("Something is wrong with the inter pairs: read1[%s](e%d-s%d +1 )=%d read2[%s](e%d-s%d + 1 )=%d", 
	  read1, currGIR->readEnd1, currGIR->readStart1, (currGIR->readEnd1 - currGIR->readStart1 + 1),
	  read2, currGIR->readEnd2, currGIR->readStart2, (currGIR->readEnd2 - currGIR->readStart2 + 1) );
    if( (currGIR->readEnd1 - currGIR->readStart1 + 1) != strlen(read1) &
	(currGIR->readEnd2 - currGIR->readStart2 + 1) != strlen(read2) ) {
      numInters += 0.25;
    } else if ( (currGIR->readEnd1 - currGIR->readStart1 + 1) != strlen(read1) |
		(currGIR->readEnd2 - currGIR->readStart2 + 1) != strlen(read2) ) {
      numInters += 0.5;
    } else {
      numInters += 1.0;
    }
  }
  currGE->numInter = numInters;
}


void writeFasta( GfrEntry* currGE, unsigned int *minReadSize, char* directory ) 
{
  FILE *freads;
  int l, readSize1, readSize2;
  // creating one fasta files with the reads
  Stringa readsFA = stringCreate( 100 ); 
  
  // creating the fasta files with the reads 
  stringPrintf( readsFA, "%s/%s_reads.fa", directory, currGE->id);
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
    fprintf( freads, ">%d/1\n%s\n>%d/2\n%s\n", l, currRead1, l, currRead2 );
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

Array util_readKnownGeneXrefs (char* fileName)
{
  WordIter w;
  LineStream ls;
  char *line,*pos;
  Array kgXrefs;
  KgXref *currKgXref;

  kgXrefs = arrayCreate (50000,KgXref);
  ls = ls_createFromFile (fileName);
  while (line = ls_nextLine (ls)) {
    if (line[0] == '\0') {
      continue;
    }
    currKgXref = arrayp (kgXrefs,arrayMax (kgXrefs),KgXref);
    w = wordIterCreate (line,"\t",0);
    currKgXref->transcriptName = hlr_strdup (wordNext (w));
    wordNext (w);
    currKgXref->swissProt = hlr_strdup (wordNext (w));
    currKgXref->uniprotId = hlr_strdup (wordNext (w));
    currKgXref->geneSymbol = hlr_strdup (wordNext (w));
    currKgXref->refseqId = hlr_strdup (wordNext (w));
    wordNext (w);
    currKgXref->refseqDescription = hlr_strdup (wordNext (w));
    wordIterDestroy (w);
    if (pos = strchr (currKgXref->uniprotId,'-')) {
      *pos = '\0';
    }
    if (pos = strchr (currKgXref->swissProt,'-')) {
      *pos = '\0';
    }
  }
  ls_destroy (ls);
  return kgXrefs;
}



Array util_readKnownGeneTreeFams (char* fileName) 
{
  LineStream ls;
  char* line;
  char* pos;
  Array kgTreeFams;
  KgTreeFam *currKgTreeFam;

  kgTreeFams = arrayCreate (30000,KgTreeFam);
  ls = ls_createFromFile (fileName);
  while (line = ls_nextLine (ls)) {
    if (line[0] == '\0') {
      continue;
    }
    if (pos = strchr (line,'\t')) {
      *pos = '\0';
      currKgTreeFam = arrayp (kgTreeFams,arrayMax (kgTreeFams),KgTreeFam);
      currKgTreeFam->transcriptName = hlr_strdup (line);
      currKgTreeFam->treeFamId = hlr_strdup (pos + 1);
    }
  }
  ls_destroy (ls);
  return kgTreeFams;
}



int sortKgXrefsByTranscriptName (KgXref *a, KgXref *b) 
{
  return strcmp (a->transcriptName,b->transcriptName);
}



static char* convert2string (Texta t)
{
  static Stringa buffer = NULL;
  int i;

  stringCreateClear (buffer,100);
  for (i = 0; i < arrayMax (t); i++) {
    stringAppendf (buffer,"%s%s",textItem (t,i),i < arrayMax (t) - 1 ? "|" : "");
  }
  return string (buffer);
}



void transcript2geneSymbolAndGeneDescription (Array kgXrefs, char *transcriptName, char** geneSymbol, char **description)
{
  Texta tokens;
  int i;
  KgXref testKX,*currKX;
  int index; 
  static Texta descriptions = NULL;
  static Texta geneSymbols = NULL;
  
  textCreateClear (descriptions,100);
  textCreateClear (geneSymbols,100);
  tokens = textFieldtokP (transcriptName,"|");
  for (i = 0; i < arrayMax (tokens); i++) {
    testKX.transcriptName = hlr_strdup (textItem (tokens,i));
    if (!arrayFind (kgXrefs,&testKX,&index,(ARRAYORDERF)sortKgXrefsByTranscriptName)) {
      warn ("Expected to find KgXref: %s",testKX.transcriptName);
      textAdd ( geneSymbols, testKX.transcriptName );
    } else {
      currKX = arrp (kgXrefs,index,KgXref);
      if (currKX->refseqDescription[0] != '\0') {
	textAdd (descriptions,currKX->refseqDescription);
      }
      if (currKX->geneSymbol[0] != '\0') {
	textAdd (geneSymbols,currKX->geneSymbol);
      } 
    }
    hlr_free (testKX.transcriptName);
  }
  textDestroy (tokens);
  textUniqKeepOrder (descriptions);
  textUniqKeepOrder (geneSymbols);
  *geneSymbol = hlr_strdup (convert2string (geneSymbols));
  *description = hlr_strdup (convert2string (descriptions));
}

