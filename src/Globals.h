#ifndef _GLOBALS_H_
#define _GLOBALS_H_

#include "alphabet.h"
#include "seqset.h"
#include "log.h"
#include "ThresholdChecker.h"

#include "getopt_pp/getopt_pp.h"
#include "memoryPool/pool_alloc.h"
#include <string.h>
#include <stdlib.h>
#include <unordered_set>
#include <math.h>
#include <sys/stat.h>

#ifdef NDEBUG
const bool kDebug = false;
#else
const bool kDebug = true;
#endif  // DEBUG

// maximum number of mutations calibrated in conservation track
#define MAX_CONS_MUTATIONS 30
// maximum fraction of mutated residues
#define MAX_MUTATED_FRACTION 1
// maximum calibrated length of conservation track
#define MAX_CONSERVATION_LENGTH 15
#define PWM_LENGTH 30

// indices for A, C, G, T mutations in conservation track
static const int mutIndex[4] = {( int )pow( 2, 0 ), ( int )pow( 2, 4 ), ( int )pow( 2, 8 ), ( int )pow( 2, 12 )};

enum motif_type{
	ALL,
	PALINDROME,
	TANDEM,
	NOPALINDROME,
	NOTANDEM,
	FIVEMERS,
	NO_VALID_MOTIF_TYPE
};

enum merge_type{
	LOW,
	MEDIUM,
	HIGH,
	NO_VALID_MERGE_MODE
};
std::ostream& operator<<( std::ostream &os, const merge_type &m );

enum TerminusMode{
	NONE,
	POS,
	NEG,
	BOTH
};

// type of supplementary information to use
enum SuppInfMode{
	SUPP_NO,
	SUPP_DISOCONS,
	SUPP_NNET
};
std::ostream& operator<<( std::ostream &os, const SuppInfMode &v );

enum Types_t{
	CS_63,
	BLOSUM45_21,
	BLOSUM62_21,
	BLOSUM80_21,
	GONNET_21,
	UNDEF
};

class StateType{

public:

	Types_t type;
	StateType() : type( UNDEF ){}
	StateType( const Types_t &t ) : type( t ){}
	bool operator==( const Types_t &t ) const{
		return type == t;
	}
	bool operator!=( const Types_t &t ) const{
		return !( *this==t );
	}
	static std::string toString( const StateType &t ){
		switch( t.type ){
		case CS_63 : return std::string( "CS_63" );
		case BLOSUM45_21 : return std::string( "BLOSUM45_21" );
		case BLOSUM62_21 : return std::string( "BLOSUM62_21" );
		case BLOSUM80_21 : return std::string( "BLOSUM80_21" );
		case GONNET_21 : return std::string( "GONNET_21" );
		case UNDEF : return std::string( "UNDEF" );
		default : std::cerr << "Unknown StateType!" << std::endl; exit(1);
		}
	}
	static Types_t fromString( const std::string &s ){
		if( s == "CS_63" ) return CS_63;
		else if( s=="BLOSUM45_21" ) return BLOSUM45_21;
		else if( s=="BLOSUM62_21" ) return BLOSUM62_21;
		else if( s=="BLOSUM80_21" ) return BLOSUM80_21;
		else if( s=="GONNET_21" ) return GONNET_21;
		else{ assert( false ); return UNDEF; }
	}
	static int sizeOfAlphabet( const StateType &t ){
		switch( t.type ){
		case CS_63 : return 63;
		case BLOSUM45_21 : return 21;
		case BLOSUM62_21 : return 21;
		case BLOSUM80_21 : return 21;
		case GONNET_21 : return 21;
		case UNDEF : return std::numeric_limits<int>::min();
		default : std::cerr << "Unknown StateType!" << std::endl; exit(1);
		}
	}
};

typedef std::list<int, Pool_alloc<int> > motif_columns_type;

class Global{

public:

	Global( int nopt, char *options[] );
	~Global();

	static a_type 		A;						// alphabet
	static ss_type 		posSet;					// positive sequences
	static ss_type 		negSet;					// negative/background sequences

	static bool 		usePositionalProbs;
	static int 			backgroundOrder;

	static bool			positionalProbsRanking;

	static int			GAPS;					// number of gap combinations in start motifs to consider
	static merge_type	mergeMode;
	static double		gapOpening;
	static double		gapExtension;
	static int 			maxMultipleSequences;
	static int 			conservationLength;
	static int			maxPosSetSize;
	static double*		motifNbCorrection;
	static ThresholdChecker instanceThreshold;	// threshold to decide whether an instance is a motif
	static bool         removeHomology;			// remove perfectly conserved motifs

	static double		consCorrection;			// p-value correction for conservation p-value
	static double		overrepCorrection;		// p-value correction for overrepresentation p-value
	static double		consPvalWeight;			// p-value combination weighting
	static int			maxSeqCount;

	static int 			maxMotifsPerSequence;	// maximum number of motifs per sequence

	static bool 	    useRankPvalues;
	static bool			useAliFree;

	static bool			maximizeMotifLength;
	static bool			noRefinementPhase;

	static bool   		multipleOccurrence;		// multiple occurrence per sequence
	static bool   		oneOccurrence;			// one occurrence per sequence
	static bool			zeroOrOneOccurrence;	// zero or one occurrence per sequence
	static bool			ungappedOutput;			// output file in PWM directory has gapped output
	static bool			revcomp;				// search on reverse complements of positive sequences too
	static bool 		repeatFiltering;
	static bool			lowComplexityFilter;

	static char* 		startMotif;				// start motif (IUPAC pattern string) for motif discovery
	static char*		profFile;				// file with start profile (PWM) for motif discovery
	static int			startRegion;			// expected start position of region enriched for motif occurrences
	static int 			endRegion;				// expected end position of region enriched for motif occurrences
	static motif_type 	type;					// seed pattern types: ALL, FIVEMERS, PALINDROME, TANDEM, NOPALINDROME, NOTANDEM
	static seq_format	seqFormat;				// format of positive and negative/background sequence sets: FASTA, CLUSTALW

  	static float***		conservationProbs;		// conservation probability for a n-mer with mutations in n species
  	static float***		alignmentFreeProbs;		// conservation probability for a n-mer with mutations in n species

  	static double* 		posBg_log;				// logarithm of base frequencies in positive sequences
  	static double* 		posBg;					// base frequencies in positive sequences
  	static double* 		negBg_log;				// logarithm of base frequencies in negative/background sequences
  	static double*		negBg;					// base frequencies in negative/background sequences

    static double		pseudo;					// fraction of PWM pseudocounts
    static double		plusFrac;				// plus fraction in motif iteration

    static int			neff_pwm;				// effective number of different bases in single PWM columns
    static int 			neff_discrete;			// effective number of different bases in single IUPAC extensions

	static int 			downstream;				// distance between the anchor position and the end of positive sequences

	static char* 		outputDirectory;		// output directory for the results
	static char* 		tmpDirectory;			// temporary XXmotif directory
	static char*		name;					// positive sequence file name
	static char*		shortFileName;			// positive sequence file basename
	static char*		negFile;				// negative/background sequence file name
	static char* 		benchmarkFolder;		// directory for the benchmark results
	static char* 		pwmFolder;				// directory for PWM results of Bulyk benchmark

	static int			maxMotifLevel;			// maximum number of extensions to consider per level
	static double		minCoverage;			// minimum number of sequences with an instance of the motif

	static int			minMatchPositions;		// minimum number of non-wildcard motif positions
	static int			maxMatchPositions;		// maximum number of non-wildcard motif positions

	// homogeneous background BMM options
	static int order;							// model order
	static float pseudocountsFactor;			// prior strength
	static float countsOffset;					// counts offset

	/*
	 * BaMM!motif options
	 */

	static bool	em;

	/*
	 * Options to initialize a single BMM from file
	 */

	// file with binding sites of equal length
	// one binding site per line
	static char* bindingSiteFile;
	// length of binding sites
	static int bindingSiteLength;
	// file with BMM probabilities
	// omit .conds and .probs filename extensions
	static char* markovModelFile;
	// length of Markov model
	static int markovModelLength;

	/*
	 * Options to initialize one or more BMMs from XXmotif PWMs
	 */

	// minimum number of PWMs
	// options pValueThreshold and minOccurrence are ignored
	static int minModels;
	// maximum number of PWMs
	static int maxModels;
	// maximum p-value of PWMs
	// filter is not applied to the top minimum number of PWMs (minModels)
	static double pValueThreshold;
	// minimum fraction of sequences that contain the motif
	// filter is not applied to the top minimum number of PWMs (minModels)
	static float minOccurrence;
	// PWM ranks in XXmotif results
	// former options to initialize BMMs from PWMs are ignored
	static std::vector<int> nrModels;

	// set the prior probability for a positive sequence to contain a motif to
	//   the fraction of sequences that contain a binding site instance that
	//   XXmotif used to build the PWM
	// option only works for a single PWM
	static bool msq;

	/*
	 * Options for (inhomogeneous) motif BMMs
	 */

	// order
	static int modelOrder;
	// order-specific prior strength
	static std::vector<float> alpha;
	// calculate order-specific alphas according to beta x gamma^(k-1)
	static float beta;
	// calculate order-specific alphas according to beta x gamma^(k-1)
	static float gamma;
    // use position-specific alphas
    static bool positionSpecificAlphas;
    // calculate prior probabilities from lower-order probabilities instead of
    //   background frequencies of mononucleotides
	static bool interpolate;
	// add uniformly initialized positions to the left/right of initial BMMs
	static std::vector<int> addColumns;

	/*
	 * Options for the (homogeneous) background BMM
	 */

	// order
	static int modelOrderBg;
	// prior strength
	static float alphaBg;

	/*
	 * EM options
	 */

	// prior probability for a positive sequence to contain a motif
	static float q;
	static float qmax;
	// the EM algorithm is deemed to be converged when the sum over the absolute
	//   differences in BMM probabilities from successive EM rounds is smaller
	//   than epsilon
	static float epsilon;
	// the EM algorithm is deemed to be converged when the likelihood converges
	static bool likelihoodConvergence;
	// limit the number of EM iterations
	static int maxEMIterations;
    // initialize BMMs only
	static bool noExpectationMaximizationPhase;

	// optimize alphas
    static bool learnHyperParameter;
    // verbose printouts to debug alpha learning code
    static bool debugAlphalearning;

    // calculate BMM probabilities using pseudocounts from the previous EM
    //   iteration
	static bool lastCondsPseudoCounts;
	// calculate 0'th-order BMM probabilities using pseudocounts calculated from
	//   initial 0'th-order BMM probabilities
	static bool monoProbsPseudoCounts;
	// calculate 0'th-order BMM probabilities using pseudocounts calculated from
	//   initial 0'th-order BMM probabilities and alpha = N*q*alphaZeroFactor
	static float alphaZeroFactor;

	/*
	 * Options for weighting positive sequences
	 */

	// file with intensities for positive sequences (one intensity per line)
	//    used to weight sequences in the EM algorithm
	// the order of intensities must conform to the order of positive sequences
	// higher intensities produce higher sequence weights
	static char* sequenceIntsFile;
	// use intensities to initialize BMMs from weighted instances of XXmotif
	static bool initInts;
	// use intensity ranks instead of intensities to calculate weights
	static bool rankWeighting;
	// the quantile of intensities (or ranks) that defines the background
	//   intensity (rank) used to translate intensities (ranks) into weights
	// the weight of sequences with intensities (ranks) below (above) the
	//   background intensity (rank) is set to zero
	static float backgroundQuantile;
	// the intensity that defines the background intensity used to translate
	//   intensities into weights
	// the weight of sequences with intensities below the background intensity
	//   is set to zero
	static float backgroundIntensity;
	// the rank that defines the background rank used to translate ranks into
	//   weights
	// the weight of sequences with ranks above the background rank is set to
	//   zero
	static float backgroundRank;

	/*
	 * Options for weighting binding sites
	 */

	// file with intensities for binding site sequences (one per line) used to
	//   initialize BMMs from weighted binding sites
	// the order of intensities must conform to the order of sequences in the
	//   binding sites file
	// higher intensities produce higher weights
	static char* bindingSiteIntsFile;
	// use intensity ranks instead of intensities to calculate weights
	static bool bindingSiteRankWeighting;
	// the quantile of intensities (or ranks) that defines the background
	//   intensity (rank) used to translate intensities (ranks) into weights
	// the weight of binding sites with intensities (ranks) below (above) the
	//   background intensity (rank) is set to zero
	static float bindingSiteBackgroundQuantile;
	// the intensity that defines the background intensity used to translate
	//   intensities into weights
	// the weight of binding sites with intensities below the background
	//   intensity is set to zero
	static float bindingSiteBackgroundIntensity;
	// the rank that defines the background rank used to translate ranks into
	//   weights
	// the weight of binding sites with ranks above the background rank is set
	//   to zero
	static float bindingSiteBackgroundRank;

	/*
	 * Options to score sequences
	 */

	// score positive (training) sequences with optimized BMMs
	static bool testPosSequences;
	// score background (training) sequences with optimized BMMs
	static bool testNegSequences;
	// score test sequences with optimized BMMs
	// test sequences can be provided in a single or multiple FASTA files
	static std::vector<std::string> testSequenceFile;
	// score sequences with XXmotif PWMs
	static bool evaluatePWMs;
	// calculate log instead of log-odds scores
	static bool logProbs;

	/*
	 * Output options
	 */

	// write initialized BMM(s) to disk
	static bool saveInitModels;
	// write optimized BMM(s) to disk
	static bool saveModels;
	// write sequence likelihoods and positional odds scores to disk after each
	//   EM iteration
	static bool saveExpectationMaximizationLikelihoods;
	// write BBMs to disk after each EM iteration
	static bool saveExpectationMaximizationModels;
	// verbose terminal printouts
	static bool	verbose;

	/*
	 * BaMM!motif variables
	 */

	// momomer background frequencies
	static double* freqs;
	// calculate background probabilities for k-mers with gaps
	// gaps are mandatory in order to initialize from XXmotif PWMs
	static bool gaps;
	// list of test sequences to score
	static std::list<ss_type> testSet;

	// CS-BLAST variables
	static int			cswlen;
	static std::string	csprofiles;
	static int 			csbest;

	// MadonaPro variables
	static double		disoconsWeight;		// weight for p-value combination
	static TerminusMode termMode;			// append dollar to positive sequence termini?
	static float		dollarScaleFactor;  // factor to scale $ frequencies
	static StateType	type_of_states;		// type of profile states
	static int			maxIterations;		// maximum number of PWM iterations
	typedef std::unordered_set<std::string> tracked_t;
	static tracked_t trackedMotifs;			// all dimer/trimer substrings are tracked, i.e. [KR][DE]EL
	static bool isTracked( const std::string &s );
	static bool 		trackedElongation;
	static bool			trackedOnly;		// use only tracked motifs as initial seeds
	static double		extensionECut;		// initial seeds with E <= cut are extended
	static int			extensionMinCut;	// minimum number of extended seeds (if enough)
	static int			extensionMaxCut;	// maximum number of extended seeds
	static double		aaStateSigThresh;	// minimum score s[a] to consider a state relevant for amino acid s
	static double		aaSeqFreqThresh;	// minimum fraction of sequences with aa conserved for extension
	static bool			batch_mode;			// if running non-interactively (suppress progress indicators)

	static bool fixedPosition;
	static double finalFilterThreshold;		// e-value threshold for filtering final motifs

	static bool 		DEBUG;
	static std::string	argv0;				// name of executable as given on the command line

	static std::string nnetFilename;

	static char* String( const char *s );
	static void createDirectory( const char *s );
	static void deleteDirectory( const char *s );

private:

	bool readCommandLineOptions( int argc, char *argv[] );
	void printHelp();
};

inline char* Global::String( const char *s ){
	return strdup( s );
}

inline void Global::createDirectory( const char *s ){
	struct stat St;
	if( stat( s, &St ) != 0 ){
		char* command = ( char* )calloc( 1024, sizeof( char ) );
		sprintf( command, "mkdir %s", s );
		if( system( command ) != 0 ){
			fprintf( stderr, "Error: Output directory %s could not be created\n", s );
			exit( -1 );
		}
		free( command );
	}
}

inline void Global::deleteDirectory( const char *s ){
	struct stat St;
	if( stat( s, &St ) == 0 ){
		char* command = ( char* )calloc( 1024, sizeof( char ) );
		sprintf( command, "rm -rf %s", s );
		if( system( command ) != 0 ){
			fprintf( stderr, "Error: Temporary directory %s could not be deleted\n", s );
			exit( -1 );
		}
		free( command );
	}
}

namespace GetOpt
{
	template <> inline _Option::Result convert<char*>( const std::string& s, char*& d, std::ios::fmtflags ){
		_Option::Result ret = _Option::BadType;
		d = Global::String( s.c_str() );
		ret = _Option::OK;
		return ret;
	}

	template <> inline _Option::Result convert<motif_type>( const std::string& s, motif_type& d, std::ios::fmtflags ){
		_Option::Result ret = _Option::BadType;

		if( s.compare( "PALINDROME" ) == 0 ){
			d = PALINDROME; std::cerr << "Motif seeds: PALINDROME" << std::endl;
			ret = _Option::OK;
		}else if( s.compare( "TANDEM" ) == 0 ){
			d = TANDEM; std::cerr << "Motif seeds: TANDEM" << std::endl;
			ret = _Option::OK;
		}else if( s.compare( "NOPALINDROME" ) == 0 ){
			d = NOPALINDROME; std::cerr << "Motif seeds: NOPALINDROME" << std::endl;
			ret = _Option::OK;
		}else if( s.compare( "NOTANDEM" ) == 0 ){
			d = NOTANDEM; std::cerr << "Motif seeds: NOTANDEM" << std::endl;
			ret = _Option::OK;
		}else if( s.compare( "ALL" ) == 0 ){
			d = ALL; std::cerr << "Motif seeds: ALL" << std::endl;
			ret = _Option::OK;
		}else if( s.compare( "FIVEMERS" ) == 0 ){
			d = FIVEMERS;
			std::cerr << "Motif seeds: FIVEMERS" << std::endl;
			ret = _Option::OK;
		}else {
			d = NO_VALID_MOTIF_TYPE;
			std::cerr << "Error: Motif seeds \"" << s << "\" not accepted" << std::endl;
		}
		return ret;
	}

	template <> inline _Option::Result convert<merge_type>( const std::string& s, merge_type& d, std::ios::fmtflags ){
		_Option::Result ret = _Option::BadType;

		if( s.compare( "LOW" ) == 0 ){
			d = LOW;
			std::cerr << "Similarity threshold for merging motifs set to: LOW" << std::endl;
			ret = _Option::OK;
		}else if( s.compare( "MEDIUM" ) == 0 ){
			d = MEDIUM;
			std::cerr << "Similarity threshold for merging motifs set to: MEDIUM" << std::endl;
			ret = _Option::OK;
		}else if( s.compare( "HIGH" ) == 0 ){
			d = HIGH;
			std::cerr << "Similarity threshold for merging motifs set to: HIGH" << std::endl;
			ret = _Option::OK;
		}else{
			d = NO_VALID_MERGE_MODE;
			std::cerr << "Error: Similarity threshold for merging motifs \"" << s << "\" not accepted" << std::endl;
		}
		return ret;
	}

	template <> inline _Option::Result convert<seq_format>( const std::string& s, seq_format& d, std::ios::fmtflags ){
		_Option::Result ret = _Option::BadType;

		if( s.compare( "CLUSTALW" ) == 0 ){
			d = CLUSTALW;
			std::cout << "Format of input sequences: CLUSTALW" << std::endl;
			ret = _Option::OK;
		}else if( s.compare( "FASTA" ) == 0 ){
			d = FASTA;
			std::cout << "Format of input sequences: FASTA" << std::endl;
			ret = _Option::OK;
		}else if( s.compare( "MFASTA" ) == 0 ){
			d = MFASTA;
			std::cout << "Format of input sequences: multiple FASTA" << std::endl;
			ret = _Option::OK;
		}else if( s.compare( "CUSTOM" ) == 0 ){
			d = CUSTOM;
			std::cout << "Format of input sequences: customized multiple FASTA" << std::endl;
			ret = _Option::OK;
		}else{
			d = NO_VALID_SEQ_FORMAT;
			std::cerr << "Error: Sequence format \"" << s << "\" not accepted" << std::endl;
		}
		return ret;
	}
}

inline bool Global::isTracked( const std::string &s ){
   for( auto m = trackedMotifs.begin(); m != trackedMotifs.end(); ++m ){
      if ( s.compare( *m ) == 0 ){
        return true;
      }
  }
  return false;
}

#endif /* _GLOBALS_H_ */
