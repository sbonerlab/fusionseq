#ifndef DEF_GFR_H
#define DEF_GFR_H



#define GFR_COLUMN_TYPE_NUM_INTER 1
#define GFR_COLUMN_TYPE_INTER_MEAN_AB 2
#define GFR_COLUMN_TYPE_INTER_MEAN_BA 3
#define GFR_COLUMN_TYPE_PVALUE_AB 4
#define GFR_COLUMN_TYPE_PVALUE_BA 5
#define GFR_COLUMN_TYPE_NUM_INTRA1 6
#define GFR_COLUMN_TYPE_NUM_INTRA2 7
#define GFR_COLUMN_TYPE_FUSION_TYPE 8
#define GFR_COLUMN_TYPE_NAME_TRANSCRIPT1 9
#define GFR_COLUMN_TYPE_NUM_EXONS_TRANSCRIPT1 10
#define GFR_COLUMN_TYPE_EXON_COORDINATES_TRANSCRIPT1 11
#define GFR_COLUMN_TYPE_CHROMOSOME_TRANSCRIPT1 12
#define GFR_COLUMN_TYPE_STRAND_TRANSCRIPT1 13
#define GFR_COLUMN_TYPE_START_TRANSCRIPT1 14
#define GFR_COLUMN_TYPE_END_TRANSCRIPT1 15
#define GFR_COLUMN_TYPE_NAME_TRANSCRIPT2 16
#define GFR_COLUMN_TYPE_NUM_EXONS_TRANSCRIPT2 17
#define GFR_COLUMN_TYPE_EXON_COORDINATES_TRANSCRIPT2 18
#define GFR_COLUMN_TYPE_CHROMOSOME_TRANSCRIPT2 19
#define GFR_COLUMN_TYPE_STRAND_TRANSCRIPT2 20
#define GFR_COLUMN_TYPE_START_TRANSCRIPT2 21
#define GFR_COLUMN_TYPE_END_TRANSCRIPT2 22
#define GFR_COLUMN_TYPE_PAIR_COUNT 23
#define GFR_COLUMN_TYPE_INTER_READS 24
#define GFR_COLUMN_TYPE_ID 25
#define GFR_COLUMN_TYPE_READS_TRANSCRIPT1 26
#define GFR_COLUMN_TYPE_READS_TRANSCRIPT2 27
#define GFR_COLUMN_TYPE_GENE_SYMBOL_TRANSCRIPT1 28
#define GFR_COLUMN_TYPE_GENE_SYMBOL_TRANSCRIPT2 29
#define GFR_COLUMN_TYPE_DESCRIPTION_TRANSCRIPT1 30
#define GFR_COLUMN_TYPE_DESCRIPTION_TRANSCRIPT2 31
#define GFR_COLUMN_TYPE_SPER 32
#define GFR_COLUMN_TYPE_DASPER 33
#define GFR_COLUMN_TYPE_RESPER 34


#define GFR_COLUMN_NAME_NUM_INTER "numInter"
#define GFR_COLUMN_NAME_INTER_MEAN_AB "interMeanAB"
#define GFR_COLUMN_NAME_INTER_MEAN_BA "interMeanBA"
#define GFR_COLUMN_NAME_PVALUE_AB "pValueAB"
#define GFR_COLUMN_NAME_PVALUE_BA "pValueBA"
#define GFR_COLUMN_NAME_NUM_INTRA1 "numIntra1"
#define GFR_COLUMN_NAME_NUM_INTRA2 "numIntra2"
#define GFR_COLUMN_NAME_FUSION_TYPE "fusionType"
#define GFR_COLUMN_NAME_NAME_TRANSCRIPT1 "nameTranscript1"
#define GFR_COLUMN_NAME_NUM_EXONS_TRANSCRIPT1 "numExonsTranscript1"
#define GFR_COLUMN_NAME_EXON_COORDINATES_TRANSCRIPT1 "exonCoordinatesTranscript1"
#define GFR_COLUMN_NAME_CHROMOSOME_TRANSCRIPT1 "chromosomeTranscript1"
#define GFR_COLUMN_NAME_STRAND_TRANSCRIPT1 "strandTranscript1"
#define GFR_COLUMN_NAME_START_TRANSCRIPT1 "startTranscript1"
#define GFR_COLUMN_NAME_END_TRANSCRIPT1 "endTranscript1"
#define GFR_COLUMN_NAME_NAME_TRANSCRIPT2 "nameTranscript2"
#define GFR_COLUMN_NAME_NUM_EXONS_TRANSCRIPT2 "numExonsTranscript2"
#define GFR_COLUMN_NAME_EXON_COORDINATES_TRANSCRIPT2 "exonCoordinatesTranscript2"
#define GFR_COLUMN_NAME_CHROMOSOME_TRANSCRIPT2 "chromosomeTranscript2"
#define GFR_COLUMN_NAME_STRAND_TRANSCRIPT2 "strandTranscript2"
#define GFR_COLUMN_NAME_START_TRANSCRIPT2 "startTranscript2"
#define GFR_COLUMN_NAME_END_TRANSCRIPT2 "endTranscript2"
#define GFR_COLUMN_NAME_PAIR_COUNT "pairCount"
#define GFR_COLUMN_NAME_INTER_READS "interReads"
#define GFR_COLUMN_NAME_ID "id"
#define GFR_COLUMN_NAME_READS_TRANSCRIPT1 "readsTranscript1"
#define GFR_COLUMN_NAME_READS_TRANSCRIPT2 "readsTranscript2"
#define GFR_COLUMN_NAME_GENE_SYMBOL_TRANSCRIPT1 "geneSymbolTranscript1"
#define GFR_COLUMN_NAME_GENE_SYMBOL_TRANSCRIPT2 "geneSymbolTranscript2"
#define GFR_COLUMN_NAME_DESCRIPTION_TRANSCRIPT1 "descriptionTranscript1"
#define GFR_COLUMN_NAME_DESCRIPTION_TRANSCRIPT2 "descriptionTranscript2"
#define GFR_COLUMN_NAME_SPER "SPER" /**< SPER */
#define GFR_COLUMN_NAME_DASPER "DASPER" /**< DASPER */
#define GFR_COLUMN_NAME_RESPER "RESPER" /**< RESPER */

#define GFR_PAIR_TYPE_EXONIC_EXONIC 1 /**< exonic-exonic = 1 */
#define GFR_PAIR_TYPE_EXONIC_INTRONIC 2 /**< exonic-intronic = 2 */
#define GFR_PAIR_TYPE_EXONIC_JUNCTION 3 /**< exonic-boundary = 3 */
#define GFR_PAIR_TYPE_INTRONIC_EXONIC 4 /**< intronic-exonic = 4 */
#define GFR_PAIR_TYPE_INTRONIC_INTRONIC 5 /**< introni-intronic = 5 */
#define GFR_PAIR_TYPE_INTRONIC_JUNCTION 6 /**< intronic-boundary = 6 */
#define GFR_PAIR_TYPE_JUNCTION_JUNCTION 7  /**< boundary-boundary = 7 */
#define GFR_PAIR_TYPE_JUNCTION_EXONIC 8  /**< boundary-exonic = 8 */
#define GFR_PAIR_TYPE_JUNCTION_INTRONIC 9  /**< boundary-intronic = 9 */




/**
    @file gfr.h 
    @brief Gene Fusion Report (GFR) data structure and functions.
    @details This file defines the main entities of a GFR file. A GFR file is a required input for all the GFR filters. 
    @author Lukas Habegger
    @author Andrea Sboner
   
*/


/**
   Virtual exon connection summary structure
*/
typedef struct {
  int pairType;/**< virtual exon connection type */
  int number1;/**< virtual exon of transcript 1  */
  int number2;/**< virtual exon of transcript 1 */
  float count;/**< number of reads for the connection, accounting for spliced reads */
} GfrPairCount;


/**
   Data structure for inter-transcript reads, including the virtual-exon connection information.
 */
typedef struct {
  int readStart1;/**< genomic start location of read 1 */
  int readStart2;/**< genomic start location of read 2 */
  int readEnd1;/**< genomic end location of read 1 */
  int readEnd2;/**< genomic end location of read 2 */
  int pairType;/**< type of connection */
  int number1;/**< virtual exon of transcript 1  */
  int number2;/**< virtual exon of transcript 2 */
  int flag;/**< valid read @attention if flagged, then the read is not considered in the analysis, but still reported */
} GfrInterRead;


/**
   exon coordinates
*/
typedef struct {
  int start; /**< genomic start position */
  int end;/**< genomic end position */
} GfrExonCoordinate;

/**
   Data structure of the fusion transcript candidates
 */
typedef struct {
  double numInter; /**< number of inter-transcript reads supporting the fusion candidate */
  double interMeanAB;/**< median insert size of the minimal fusion transcript @attention it consider transcript 1 before transcript 2 */
  double interMeanBA;/**< median insert size of the minimal fusion transcript @attention it consider transcript 2 before transcript 1 */
  double pValueAB;/**< p-value of the median insert size @attention it consider transcript 1 before transcript 2 */
  double pValueBA;/**< p-value of the median insert size @attention it consider transcript 2 before transcript 1 */
  double numIntra1;/**< number of intra-transcript reads on transcript 1 */
  double numIntra2;/**< number of intra-transcript reads on transcript 2 */
  char *fusionType;/**< fusion type: inter, intra, read-through, or cis */
  char *nameTranscript1;/**< ID of transcript 1 */
  char *chromosomeTranscript1;/**< chromosome of transcript 1 */
  char strandTranscript1;/**< strand of transcript 1 */
  int numExonsTranscript1;/**< number of exons of transcript 1 */
  Array exonCoordinatesTranscript1;/**< coordinates of exons of transcript 1 @remark type GfrExonCoordinate */
  int startTranscript1;/**< genomic start position of transcript 1 */
  int endTranscript1;/**< genomic end position of transcript 1 */
  char *geneSymbolTranscript1;/**< gene symbol of transcript 1 @remark multiple symbols are "piped" together via '|'*/
  char *descriptionTranscript1;/**< description of transcript 1 @remark multiple descriptions are "piped" together via '|'*/
  char *nameTranscript2;/**< ID of transcript 2 */
  char *chromosomeTranscript2;/**< chromosome of transcript 2 */
  char strandTranscript2;/**< strand of transcript 2 */
  int numExonsTranscript2;/**< number of exons of transcript 2 */
  Array exonCoordinatesTranscript2;/**< coordinates of exons of transcript 2 @remark type GfrExonCoordinate */
  int startTranscript2;/**< genomic start position of transcript 2 */
  int endTranscript2;/**< genomic end position of transcript 2 */
  char *geneSymbolTranscript2;/**< gene symbol of transcript 2 @remark multiple symbols are "piped" together via '|'*/
  char *descriptionTranscript2;/**< description of transcript 2 @remark multiple descriptions are "piped" together via '|'*/
  Array interReads; /**< description of connection between the virtual exons @remark type GfrInterReads */
  Array pairCounts;/**< number of inter-transcript reads for each  virtual-exon connection @remark type GfrPairCount */
  char *id;/**< ID of the fusion candidate @remark format prefix_0000# */
  Texta readsTranscript1;/**< sequences of inter-transcript reads for transcript 1 @remark the sequences are piped together via '|' */
  Texta readsTranscript2;/**< sequences of inter-transcript reads for transcript 2 @remark the sequences are piped together via '|' */
  double SPER;/**< Supportive Paired-End Reads */
  double DASPER;/**< Difference between the Analytically computed SPER and the observed SPER */
  double RESPER;/**< Ratio between the Empirically computed SPER and the observed SPER */
} GfrEntry;


/** initialization (constructor) of the gfr module. */
extern int gfr_init (char* fileName /**< [in] pointer to the filename. @remark use "-" to denote stdin */);
/** de-initialization (desctructor) of the gfr module.  @pre the gfr module has been initialized with gfr_init().*/
extern void gfr_deInit (void);
/** add a column to the gfr module.  @pre the gfr module has been initialized with gfr_init().*/
extern void gfr_addNewColumnType (char* columnName /**< [in] string encoding the column name. */);
/** obtain a pointer to the next GfrEntry. @pre the gfr module has been initialized with gfr_init(). @param [out] GfrEntry* a pointer to a GfrEntry */
extern GfrEntry* gfr_nextEntry (void);
/** Retrieve all entries from a GFR file. @return an Array with all the gfr entries. @pre the gfr module has been initialized with gfr_init(). */
extern Array gfr_parse (void);
/** write the header of a GFR file.  @pre the gfr module has been initialized with gfr_init(). */
extern char* gfr_writeHeader (void);
/** write the gfr entry to a string.  @pre the gfr module has been initialized with gfr_init(). */
extern char* gfr_writeGfrEntry (GfrEntry *currEntry /**< [in] pointer to the current gfr entry.*/);



#endif
