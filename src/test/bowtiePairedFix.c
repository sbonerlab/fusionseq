#include "log.h"
#include "format.h"
#include "mrf.h"
#include "linestream.h"
#include "stringUtil.h"
#include "common.h"
#include "seq.h"

int main (int argc, char *argv[])
{
  char* line;
  char* buffer1,*buffer2;
  char* prevID=NULL;
  char* currID=NULL;
  char *p;
  LineStream ls = ls_createFromFile ("-");
  while ( line = ls_nextLine(ls) ) {
    if( prevID==NULL ) {
      prevID = hlr_strdup( line );
      buffer1 = hlr_strdup( line );
      p = rindex( prevID, '/' );
      if( p == NULL ) die( "No '/' in the ID name as expected from bowtie output for paired end mapping.");
      *p='\0';
    } else {
      currID = hlr_strdup( line );
      buffer2 = hlr_strdup( line );
      p = rindex ( currID, '/' );
      if( p == NULL ) die( "No '/' in the ID name as expected from bowtie output for paired end mapping.");
      *p = '\0';
    } 
    if(  (currID != NULL) && (prevID!=NULL) ) {
      if( strEqual( prevID, currID ) ) {
	printf( "%s\n%s\n", buffer1, buffer2) ;
	hlr_free( prevID );
	hlr_free( currID );
	hlr_free( buffer1 );
	hlr_free( buffer2 );
	prevID=NULL;
	currID=NULL;
      } else {
	prevID = hlr_strdup(currID);
	buffer1 = hlr_strdup( buffer2);
	hlr_free( currID );
	currID=NULL;
      }
    }
  }
  ls_destroy( ls );
  return 0;
}    
