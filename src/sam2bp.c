#include <bios/log.h>
#include <bios/format.h>
#include <mrf/sam.h>



typedef struct {
  char *target;
  char *read;
  int position;
} BreakPoint;



typedef struct {
  char *tileCoordinate1;
  char *tileCoordinate2;
  Array breakPoints; // of type BreakPoint*
} SuperBreakPoint;



static int sortBreakPointsByTargetAndOffset (BreakPoint *a, BreakPoint *b) 
{
  int diff;

  diff = strcmp (a->target,b->target);
  if (diff != 0) {
    return diff;
  }
  return b->position - a->position;
}



static int sortSuperBreakPointsBySupport (SuperBreakPoint *a, SuperBreakPoint *b) 
{
  return arrayMax (b->breakPoints) - arrayMax (a->breakPoints);
}



int main (int argc, char *argv[])
{
  SamEntry *currSE;
  Array samEntries;
  Array breakPoints;
  BreakPoint *currBP,*nextBP;
  Array superBreakPoints;
  SuperBreakPoint *currSBP;
  char *targetCopy;
  char *pos;
  int i,j;

  samParser_initFromFile ("-");
  samEntries = samParser_getAllEntries ();
  samParser_deInit (); 
  breakPoints = arrayCreate (10000,BreakPoint);
  for (i = 0; i < arrayMax (samEntries); i++) {

    currSE = arrp (samEntries,0,SamEntry);
    currBP = arrayp (breakPoints,arrayMax (breakPoints),BreakPoint);
    currBP->target = hlr_strdup (currSE->rname);
    currBP->read = hlr_strdup (currSE->seq);
    currBP->position = currSE->pos;
  }
  arraySort (breakPoints,(ARRAYORDERF)sortBreakPointsByTargetAndOffset);
  targetCopy = NULL;
  superBreakPoints = arrayCreate (1000,SuperBreakPoint);
  i = 0;
  while (i < arrayMax (breakPoints)) {
    currBP = arrp (breakPoints,i,BreakPoint);
    currSBP = arrayp (superBreakPoints,arrayMax (superBreakPoints),SuperBreakPoint);
    currSBP->breakPoints = arrayCreate (100,BreakPoint*);
    strReplace (&targetCopy,currBP->target);
    pos = strchr (targetCopy,'|');
    if (pos == NULL) {
      die ("Unexpected target: %s",targetCopy);
    }
    *pos = '\0';
    currSBP->tileCoordinate1 = hlr_strdup (targetCopy);
    currSBP->tileCoordinate2 = hlr_strdup (pos + 1);
    array (currSBP->breakPoints,arrayMax (currSBP->breakPoints),BreakPoint*) = currBP;
    j = i + 1;
    while (j < arrayMax (breakPoints)) {
      nextBP = arrp (breakPoints,j,BreakPoint);
      if (strEqual (currBP->target,nextBP->target)) {
        array (currSBP->breakPoints,arrayMax (currSBP->breakPoints),BreakPoint*) = nextBP;
      }
      else {
        break;
      }
      j++;
    }
    i = j;
  }
  arraySort (superBreakPoints,(ARRAYORDERF)sortSuperBreakPointsBySupport);
  for (i = 0; i < arrayMax (superBreakPoints); i++) {
    currSBP = arrp (superBreakPoints,i,SuperBreakPoint);
    printf ("%s,%s,",currSBP->tileCoordinate1,currSBP->tileCoordinate2);
    for (j = 0; j < arrayMax (currSBP->breakPoints); j++) {
      currBP = arru (currSBP->breakPoints,j,BreakPoint*);
      printf ("%d:%s%s",currBP->position,currBP->read, j < arrayMax (currSBP->breakPoints) - 1 ? "|" : "\n");
    }
  }
  
  return 0;
}
