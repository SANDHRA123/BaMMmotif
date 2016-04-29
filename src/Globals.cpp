#include <fstream>
#include <limits.h>
#include <math.h>

#include "getopt_pp/getopt_pp.h"

#include "backgroundDistribution.h"
#include "Globals.h"
#include "NullModel.h" // homogeneous background BMM used in XXmotif
#include "em/hoUtils.h"
#include "em/hoNullModel.h" // homogeneous background BMM used in BaMM!motif
#include "refinementPhase/Motif.h"

using GetOpt::GetOpt_pp;
using GetOpt::GlobalOption;
using GetOpt::Option;
using GetOpt::OptionPresent;

a_type 		Global::A = NULL;
ss_type		Global::posSet = NULL;					// positive sequences
ss_type 	Global::negSet = NULL;					// negative/background sequences
bool		Global::usePositionalProbs = false;
bool		Global::positionalProbsRanking = false;

int 		Global::GAPS = 0;						// number of gap combinations in start motifs to consider
merge_type	Global::mergeMode;
double 		Global::gapOpening = 1;
double 		Global::gapExtension = 1;
int 		Global::maxMultipleSequences = 100;
int			Global::maxPosSetSize = INT_MAX;
ThresholdChecker	Global::instanceThreshold;

double 		Global::overrepCorrection=1.05;
double 		Global::consCorrection=1.00;
double 		Global::consPvalWeight=1.0/3;

int			Global::maxMotifsPerSequence = 1000;
int			Global::maxSeqCount = 1000;

bool 		Global::useAliFree = false;
bool 		Global::useRankPvalues = false;

bool   		Global::multipleOccurrence = false;		// multiple occurrence per sequence
bool   		Global::oneOccurrence = false;			// one occurrence per sequence
bool		Global::zeroOrOneOccurrence = false;	// zero or one occurrence per sequence
bool		Global::revcomp= false;					// search on reverse complements of positive sequences too
bool 		Global::repeatFiltering = false;
bool 		Global::lowComplexityFilter = false;
bool		Global::noRefinementPhase = false;

char* 		Global::startMotif = NULL;		// start motif (IUPAC pattern string) for motif discovery
char*		Global::profFile = NULL;		// file with start profile (PWM) for motif discovery **/
int			Global::startRegion = 0;		// expected start position of region enriched for motif occurrences
int			Global::endRegion = 1;			// expected start position of region enriched for motif occurrences
motif_type	Global::type = ALL;				// seed pattern types: ALL, FIVEMERS, PALINDROME, TANDEM, NOPALINDROME, NOTANDEM
seq_format	Global::seqFormat = FASTA;		// format of positive and negative/background sequence sets: FASTA, CLUSTALW

float***	Global::conservationProbs = NULL;
float***	Global::alignmentFreeProbs = NULL;
bool        Global::removeHomology = false;

double* 	Global::posBg_log = NULL;		// logarithm of base frequencies in positive sequences
double* 	Global::posBg = NULL;			// base frequencies in positive sequences
double* 	Global::negBg_log = NULL; 		// logarithm of base frequencies in negative/background sequences
double*		Global::negBg = NULL;			// base frequencies in negative/background sequences

double		Global::pseudo = 0.1;			// fraction of PWM pseudocounts
double		Global::plusFrac = 0;
int			Global::neff_pwm;				// effective number of different bases in single PWM columns
int 		Global::neff_discrete;			// effective number of different bases in single IUPAC extensions

int			Global::downstream = 0;			// distance between the anchor position and the end of positive sequences

char* 		Global::outputDirectory = NULL; // output directory for the results
char* 		Global::tmpDirectory = NULL;	// temporary XXmotif directory
char*		Global::name = NULL;			// positive sequence file name
char*		Global::shortFileName = NULL;	// positive sequence file basename
char*		Global::negFile = NULL;			// negative/background sequence file name
char* 		Global::benchmarkFolder = NULL; // directory for the benchmark results
char*		Global::pwmFolder = NULL;
int			Global::maxMotifLevel = 3;		// maximum number of extensions per level to consider
double		Global::minCoverage;

int			Global::minMatchPositions;
int			Global::maxMatchPositions;

bool		Global::maximizeMotifLength = true;

StateType 	Global::type_of_states;


// cs blast options
int			Global::cswlen;
std::string	Global::csprofiles;
int 		Global::csbest;

double      Global::disoconsWeight = 1.0/3;
float		Global::dollarScaleFactor = 3.0f;
TerminusMode Global::termMode = BOTH;
int			Global::maxIterations = std::numeric_limits<int>::max();
Global::tracked_t Global::trackedMotifs;
bool 		Global::trackedElongation = false;
bool		Global::trackedOnly = false;
double		Global::extensionECut;
int			Global::extensionMinCut;
int			Global::extensionMaxCut;
double		Global::aaStateSigThresh = 2.0;
double		Global::aaSeqFreqThresh = 0.75;
bool		Global::batch_mode = false;

bool 		Global::fixedPosition = false;
double 		Global::finalFilterThreshold = 1e2;

std::string Global::nnetFilename;

bool 		Global::DEBUG = false;
std::string Global::argv0;

// homogeneous background BMM options used in XXmotif
int Global::order = 2;						// model order
float Global::pseudocountsFactor = 10.0f;	// prior strength
float Global::countsOffset = 0.0f;			// counts offset

bool Global::em = true;

/*
 * Parameters to initialize models from file
 * */

/*
 * binding sites to initialize a single Markov model. Sequence lengths must not
 * differ and be provided line-by-line
 */
char* Global::bindingSiteFile = NULL;
/*
 * binding sites lengths
 */
int Global::bindingSiteLength = 30;

/*
 * Markov model file name (without ending) to initialize a single Markov model.
 * Files .conds and .probs need to be available
 */
char* Global::markovModelFile = NULL;
/*
 * Markov model size
 */
int Global::markovModelLength = 0;

/*
 * Parameters to initialize models from XXmotif results. Parameters
 * bindingSiteFile and markovModelFile must not be provided simultaneously
 * */

/*
 * number of one or more XXmotif models in the ranking used to initialize Markov
 * models. The remaining parameters available to choose models from XXmotif
 * results are ignored
 */
std::vector<int> Global::nrModels;
/*
 * min. number of XXmotif PWMs used to initialize Markov models. Independent
 * on pValueThreshold and minOccurrence
 */
int Global::minModels = 1;
/*
 * max. number of XXmotif models used to initialize Markov models
 */
int Global::maxModels = std::numeric_limits<int>::max();
/*
 * max. p-value of XXmotif models used to initialize Markov models. Not applied
 * to min. number of models
 */
double Global::pValueThreshold = 1.0;
/*
 * min. fraction of sequences containing a binding site instance. Not applied
 * to min. number of models
 */
float Global::minOccurrence = 0.05f;
/*
 * add columns to the left and right of XXmotif models used to initialize Markov
 * models
 */
std::vector<int> Global::addColumns( 2, 0 );
/*
 * use model-specific specificity factor by calculating the percentage of
 * positive sequences containing a corresponding binding site instance
 */
bool Global::msq = false;

/*
 * Markov model parameters
 * */

/*
 * Markov model order
 */
int Global::modelOrder = 2;
/*
 * order-specific prior strength. Order k fixes vector size to k+1.
 */
std::vector<float> Global::alpha( modelOrder+1, 1.0f );
/*
 * calculate order-specific alphas according to beta x gamma^(k-1)
 * (for k > 0).
 */
float Global::beta = 20.0f;
/*
 * calculate order-specific alphas according to beta x gamma^(k-1)
 * (for k > 0).
 */
float Global::gamma = 3.0f;

/*
 * learn hyper parameter alpha during interpolation of markov models
 */
bool Global::learnHyperParameter = true;
/*
 * flag in order to print output for debugging the alpha learning part only
 */
bool Global::debugAlphalearning = true;
/*
 *  whether or not to use position specific alhpa's
 */
bool Global::positionSpecificAlphas = false;

/*
 * Markov model pseudo-counts factor(s). Markov model order k fixes vector size
 * to k+1
 */
//std::vector<float> Global::eta( modelOrder+1, 90.0f );
/*
 * interpolate between higher- and lower-order probabilities
 */
bool Global::interpolate = true;

/*
 * Interpolated Markov background model parameters
 * */

/*
 * Background model order
 */
int Global::modelOrderBg = 2;
/*
 * Background model pseudo-counts factor
 */
float Global::alphaBg = 10.0f;

/*
 * EM parameters
 * */

/*
 * initialize Markov model but skip EM phase
 */
bool Global::noExpectationMaximizationPhase = false;
/*
 * specificity factor approximates the percentage of sequences contributing to
 * the Markov model
 */
float Global::q = 0.9f;
float Global::qmax = 0.99999f;
/*
 * EM convergence parameter
 */
float Global::epsilon = 0.001f;
/*
 * likelihood or max. order model parameter EM convergence
 */
bool Global::likelihoodConvergence = false;
/*
 * max. number of EM iterations
 */
int Global::maxEMIterations = std::numeric_limits<int>::max();
/*
 * update interpolated Markov model probabilities with last EM iteration's
 * pseudo-counts
 */
bool Global::lastCondsPseudoCounts = false;
/*
 * calculate 0th-order interpolated Markov model pseudo-counts from initial
 * 0th-order probabilities
 */
//bool Global::monoProbsPseudoCounts = false;
/*
 * calculate 0th-order interpolated Markov model pseudo-counts from initial
 * 0th-order probabilities using 0th-order pseudo-counts factor N * q *
 * alphaZeroFactor
 */
float Global::alphaZeroFactor = 5.0f;

/*
 * Weighting parameters
 * */

/*
 * intensity or significance values for positive sequences. The higher the
 * values the higher the weights
 */
char* Global::sequenceIntsFile = NULL;
/*
 * parameter to initialize models from XXmotif results by weighting instances
 * with corresponding sequence weigths. Option --sequenceIntsFile must be
 * provided simultaneously. Options --bindingSiteFile and --markovModelFile must
 * not be provided simultaneously
 */
bool Global::initInts = false;
/*
 * rank-based weighting. Defaults to intensity-based weighting
 */
bool Global::rankWeighting = false;
/*
 * quantile to estimate the background intensity value (or rank). Sequences
 * having their intensity value (rank) below (above) the background intensity
 * value (rank) get assigned to weight zero. Defaults to 0
 */
float Global::backgroundQuantile = 0.0f;
/*
 * background intensity value. Sequences having their intensity value below the
 * background intensity value get assigned to weight zero. Defaults to the min.
 * intensity value
 */
float Global::backgroundIntensity = std::numeric_limits<float>::min();
/*
 * background intensity rank. Sequences having their intensity rank above the
 * background intensity rank get assigned to weight zero. Defaults to the max.
 * rank
 */
float Global::backgroundRank = std::numeric_limits<float>::max();
/*
 * intensity or significance values for binding site sequences. The higher the
 * values the higher the weights. Parameter bindingSiteFile must be provided
 * simultaneously
 */
char* Global::bindingSiteIntsFile = NULL;
/*
 * binding site rank-based weighting. Defaults to intensity-based weighting
 */
bool Global::bindingSiteRankWeighting = false;
/*
 * quantile to estimate the background intensity value (or rank). Binding sites
 * having their intensity value (rank) below (above) the background intensity
 * value (rank) get assigned to weight zero. Defaults to 0
 */
float Global::bindingSiteBackgroundQuantile = 0.0f;
/*
 * background intensity value. Binding sites having their intensity value below
 * the background intensity value get assigned to weight zero. Defaults to the
 * min. binding site intensity value
 */
float Global::bindingSiteBackgroundIntensity = std::numeric_limits<float>::min();
/*
 * background intensity rank. Binding sites having their intensity rank above
 * the background intensity rank get assigned to weight zero. Defaults to the
 * max. binding site rank
 */
float Global::bindingSiteBackgroundRank =std::numeric_limits<float>::max();

/*
 * Scoring parameters
 * */

/*
 * evaluate model(s) on training sequences
 */
bool Global::testPosSequences = false;
/*
 * evaluate model(s) on background sequences
 */
bool Global::testNegSequences = false;
/*
 * evaluate model(s) on sequences in FASTA format. Specify one or more files.
 * Sequence lengths may differ
 */
std::vector<std::string> Global::testSequenceFile;
/*
 * evaluate PWM model(s) used to initialize Markov model(s) on test sequences
 */
bool Global::evaluatePWMs = false;
/*
 * calculate log probabilities instead of log likelihood ratios
 */
bool Global::logProbs = false;

/*
 * Output parameters
 * */

/*
 * save Markov models after initialization to file
 */
bool Global::saveInitModels = false;
/*
 * save Markov models after EM phase to file
 */
bool Global::saveModels = false;
/*
 * save EM iteration's sequence likelihoods and positional odds to file
 */
bool Global::saveExpectationMaximizationLikelihoods = false;
/*
 * save EM iteration's Markov models to file
 */
bool Global::saveExpectationMaximizationModels = false;
/*
 * verbose printouts
 */
bool Global::verbose = false;

/*
 * Internal parameters
 */

/*
 * monomer background frequencies
 */
double* Global::freqs = NULL;
/*
 * calculate background probabilities for k-mers with gaps. Gaps are mandatory
 * in order to initialize from XXmotif's models
 */
bool Global::gaps = true;
/*
 * list of sequences to evaluate model(s)
 */
std::list<ss_type> Global::testSet;

Global::Global( int argc, char *argv[] ){

	argv0 = std::string( argv[0] );

	if( argc < 2 ){
		 // BaMMmotif
		printHelp();
	}
	if( argc == 2 && ( strcmp( argv[1], "-h" ) == 0 || strcmp( argv[1], "--help"
			                                                 ) == 0 ) ){
		// BaMMmotif -h
		// BaMMmotif --help
		printHelp();
	}

	seq_format format = FASTA;
	for( long i = 2; i < argc; i++ ){
//		if( strcmp( argv[i], "--aa" ) == 0 ){
//			A = MkAlpha( "XACDEFGHIKLMNPQRSTVWY$" );
//			aa = true;
//		} else{
		if( strcmp( argv[i], "--terminus-mode" ) == 0 ){
			if( i+1 < argc && strcmp( argv[i+1], "NONE" ) == 0 ){
				termMode = NONE;
			} else if( i+1 < argc && strcmp( argv[i+1], "POS" ) == 0){
				termMode = POS;
			} else if( i+1 < argc && strcmp( argv[i+1], "NEG" ) == 0){
				termMode = NEG;
			} else if( i+1 < argc && strcmp( argv[i+1], "BOTH" ) == 0){
				termMode = BOTH;
			} else{
				fprintf( stderr, "Error: Illegal terminus mode\n" );
				exit( -1 );
			}
		} else if( strcmp( argv[i], "--maxPosSequences" ) == 0 ){
			if( i+1 == argc ){
				fprintf( stderr, "Error: Please provide a value for option --maxPosSequences\n" );
				exit( -1 );
			}
			maxPosSetSize = atoi( argv[i+1] );
		} else if( strcmp(argv[i], "--XX-format" ) == 0 ){
			if( i+1 == argc ){
				fprintf( stderr, "Error: Please provide a value for option --XX-format\n" );
				exit( -1 );
			}
			if( i+1 < argc && strcmp( argv[i+1], "CLUSTALW" ) == 0 ){
				format = CLUSTALW;
			}
			if( i+1 < argc && strcmp( argv[i+1], "MFASTA" ) == 0 ){
				format = MFASTA;
			}
			if( i+1 < argc && strcmp( argv[i+1], "CUSTOM" ) == 0 ){
				format = CUSTOM;
			}
		} else if( strcmp( argv[i], "--negSequenceSet" ) == 0 ){
			negFile = argv[i+1];
		} else if( strcmp( argv[i], "--lcf" ) == 0 ){
			lowComplexityFilter = true;
		}
	}

	struct stat sts;
	if( ( ( stat( argv[2], &sts ) ) == -1 ) || S_ISDIR( sts.st_mode ) ){
		fprintf( stderr, "Error: Positive sequence file %s not found\n", argv[2] );
		exit( -1 );
	}
	if( negFile != NULL && ( ( ( stat( negFile, &sts ) ) == -1 ) || S_ISDIR( sts.st_mode ) ) ){
		fprintf( stderr, "Error: Negative/background sequence file %s not found\n", negFile );
		exit( -1 );
	}

	termMode = NONE;

	if( A == NULL ){
		A = MkAlpha( "NACGT" );
	}

	// read in positive sequences
	posSet = readSeqSet( argv[2], A, format, maxPosSetSize, termMode==BOTH || termMode == POS );
	if( posSet->min_leng != posSet->max_leng ){
		fprintf( stderr, "Error: Please provide positive sequences of equal length\n" );
		exit( -1 );
	}
	calculateSequenceFeatures( posSet, A );

    posBg_log = ( double* )calloc( nAlpha(A)+1, sizeof( double ) );
    posBg     = ( double* )calloc( nAlpha(A)+1, sizeof( double ) );
    negBg_log = ( double* )calloc( nAlpha(A)+1, sizeof( double ) );
    negBg 	  = ( double* )calloc( nAlpha(A)+1, sizeof( double ) );

	startRegion = 1;
	endRegion = Global::posSet->max_leng;

	if( !( readCommandLineOptions( argc, argv ) ) ){
		printHelp();
	}
	if( negFile != NULL ){
		negSet = readSeqSet( negFile, A, format, INT_MAX );
	}
	if( negSet != NULL ){
		if( repeatFiltering ){
			filter_repeats( negSet, A );
		}
		if( lowComplexityFilter ){
			cerr << "Negative sequence set: ";
			filter_lowComplexity( negSet, A );
		}
		if( revcomp ){
			createRevcomp( negSet, A );
		}
		fillGapsInMultipleAlignment( negSet, A );
		filterMaxMultipleSequences( negSet, maxMultipleSequences );
		calculateSequenceFeatures( negSet, A );
	}

	if( repeatFiltering ){
		filter_repeats( posSet, A );
	}
	if( lowComplexityFilter ){
		cerr << "Positive sequence set: ";
		filter_lowComplexity( posSet, A );
	}
	if( revcomp ){
		createRevcomp( posSet, A );
	}
	fillGapsInMultipleAlignment( posSet, A );
	filterMaxMultipleSequences( posSet, maxMultipleSequences );

	checkSequenceSet(); // do some tests whether the negative set is a good choice
	setBackgroundDistribution(); // calculate trimer probabilities and conditional probabilities from negative/background sequences
}

void Global::printHelp(){

	bool developerHelp = false; // show developer-specific options

	printf( "\n" );
	printf( "SYNOPSIS\n" );
	printf( "      BaMMmotif DIRPATH FILEPATH [OPTIONS]\n\n" );
	printf( "DESCRIPTION\n" );
	printf( "      Bayesian Markov Model motif discovery software.\n\n" );
	printf( "      DIRPATH\n"
			"          Output directory for the results.\n\n" );
	printf( "      FILEPATH\n"
			"          FASTA file with positive sequences of equal length.\n\n" );

	printf("OPTIONS\n");
	printf( "  Sequence options\n" );
	printf( "      --negSequenceSet <FILEPATH>\n"
			"          FASTA file with negative/background sequences used to learn the\n"
			"          (homogeneous) background BMM. If not specified, the background BMM is\n"
			"          learned from the positive sequences.\n\n" );
	printf( "      --reverseComp\n"
			"          Search motifs on both strands (positive sequences and reverse\n"
			"          complements). This option is e.g. recommended when using sequences\n"
			"          derived from ChIP-seq experiments.\n\n" );
	if( developerHelp ){
		printf( "      --maxPosSequences <INTEGER> (*)\n"
				"          Maximum number of positive sequences to read in. The default is to\n"
				"          read in all sequences.\n\n");
	}

	printf( "  Options to initialize a single BMM from file\n" );
	printf( "      --bindingSiteFile <FILEPATH>\n"
			"          File with binding sites of equal length (one per line).\n\n" );
	printf( "      --markovModelFile <FILEPATH>\n"
			"          File with BMM probabilities as obtained from BaMM!motif (omit\n"
			"          filename extension).\n\n" );

	printf( "  Options to initialize one or more BMMs from XXmotif PWMs\n" );
	printf( "      --minPWMs <INTEGER>\n"
			"          Minimum number of PWMs. The options --maxPValue and --minOccurrence\n"
			"          are ignored. The default is 1.\n\n" );
	printf( "      --maxPWMs <INTEGER>\n"
			"          Maximum number of PWMs.\n\n" );
	printf( "      --maxPValue <FLOAT>\n"
			"          Maximum p-value of PWMs. This filter is not applied to the top\n"
			"          minimum number of PWMs (see --minPWMs). The default is 1.0.\n\n" );
	printf( "      --minOccurrence <FLOAT>\n"
			"          Minimum fraction of sequences that contain the motif. This filter is\n"
			"          not applied to the top minimum number of PWMs (see --minPWMs). The\n"
			"          default is 0.05.\n\n" );
	printf( "      --rankPWMs <INTEGER> [<INTEGER>...]\n"
			"          PWM ranks in XXmotif results. The former options to initialize BMMs\n"
			"          from PWMs are ignored.\n\n" );
	if( developerHelp ){
		printf( "      --msq (*)\n"
				"          Set the prior probability for a positive sequence to contain a motif\n"
				"          (see -q) to the fraction of sequences that contain a binding site\n"
				"          instance that XXmotif used to build the PWM. Since the prior\n"
				"          probability is a global parameter this option only works for a single\n"
				"          PWM.\n\n" );
	}

	printf( "  Options for (inhomogeneous) motif BMMs\n" );
	printf( "      -k <INTEGER>\n"
			"          Order. The default is 2.\n\n" );
	printf( "      -a|--alpha <FLOAT> [<FLOAT>...]\n"
			"          Order-specific prior strength. The default is 1.0 (for k = 0) and\n"
			"          20 x 3^(k-1) (for k > 0). The options -b and -g are ignored.\n\n" );
	printf( "      -b|--beta <FLOAT>\n"
			"          Calculate order-specific alphas according to beta x gamma^(k-1) (for\n"
			"          k > 0). The default is 20.0.\n\n" );
	printf( "      -g|--gamma <FLOAT>\n"
			"          Calculate order-specific alphas according to beta x gamma^(k-1) (for\n"
			"          k > 0). The default is 3.0.\n\n" );
	printf( "      --extend <INTEGER>{1,2}\n"
			"          Extend BMMs by adding uniformly initialized positions to the left\n"
			"          and/or right of initial BMMs. Invoking e.g. with --extend 0 2 adds\n"
			"          two positions to the right of initial BMMs. Invoking with --extend 2\n"
			"          adds two positions to both sides of initial BMMs. By default, BMMs\n"
			"          are not being extended.\n\n" );
	if( developerHelp ){
		printf( "      --nonBayesian (*)\n"
				"          Calculate prior probabilities from background frequencies of\n"
				"          mononucleotides instead of lower-order probabilities.\n\n" );
	}

	printf( "  Options for the (homogeneous) background BMM\n" );
	printf( "      -K <INTEGER>\n"
			"          Order. The default is 2.\n\n" );
	printf( "      -A|--Alpha <FLOAT>\n"
			"          Prior strength. The default is 10.0.\n\n" );

	printf( "  EM options.\n" );
	printf( "      -q <FLOAT>\n"
			"          Prior probability for a positive sequence to contain a motif. The\n"
			"          default is 0.9.\n\n" );
	printf( "      -e|--epsilon <FLOAT>\n"
			"          The EM algorithm is deemed to be converged when the sum over the\n"
			"          absolute differences in BMM probabilities from successive EM rounds\n"
			"          is smaller than epsilon. The default is 0.001.\n\n" );
	if( developerHelp ){
		printf( "      --maxEMIterations <INTEGER> (*)\n"
				"          Limit the number of EM iterations.\n\n" );
		printf( "      --noEM (*)\n"
				"          Initialize BMMs only.\n\n" );
	}

	printf( "  XXmotif options\n" );
	printf( "      --XX-ZOOPS\n"
			"          Use the zero-or-one-occurrence-per-sequence model (default).\n\n" );
	printf( "      --XX-MOPS\n"
			"          Use the multiple-occurrence-per-sequence model.\n\n" );
	printf( "      --XX-OOPS\n"
			"          Use the one-occurrence-per-sequence model.\n\n" );
	printf( "      --XX-seeds ALL|FIVEMERS|PALINDROME|TANDEM|NOPALINDROME|NOTANDEM\n"
			"          Define the nature of seed patterns. The default is to start using ALL\n"
			"          seed pattern variants.\n\n" );
	printf( "      --XX-gaps 0|1|2|3\n"
			"          Maximum number of gaps used for seed patterns. The default is 0.\n\n");
	printf( "      --XX-pseudoCounts <FLOAT>\n"
			"          Percentage of pseudocounts. The default is 10.0.\n\n" );
	printf( "      --XX-mergeMotifsThreshold LOW|MEDIUM|HIGH\n"
			"          Define the similarity threshold used to merge PWMs. The default is to\n"
			"          merge PWMs with LOW similarity in order to reduce runtime.\n\n");
	printf( "      --XX-maxPositions <INTEGER>\n"
			"          Limit the number of motif positions to reduce runtime. The default is\n"
			"          17.\n\n" );
	printf( "      --XX-noLengthOptimPWMs\n"
			"          Omit the length optimization of PWMs.\n\n");
	printf( "      --XX-K <INTEGER>\n"
			"          Order of the (homogeneous) background BMM. The default is either 2\n"
			"          (when learned on positive sequences) or 8 (when learned on background\n"
			"          sequences).\n\n" );
	printf( "      --XX-A <FLOAT>\n"
			"          Prior strength of the (homogeneous) background BMM. The default is\n"
			"          10.0.\n\n");
	printf( "      --XX-jumpStartPatternStage <STRING>\n"
			"          Jump-start pattern stage using an IUPAC pattern string.\n\n");
	printf( "      --XX-jumpStartPWMStage <FILEPATH>\n"
			"          Jump-start PWM stage reading in a PWM from file.\n\n" );
	if( developerHelp ){
		printf( "      --XX-minMatchPositions <INTEGER> (*)\n"
				"          Minimum number of non-wildcard motif positions. The default is 4.\n\n" );
		printf( "      --XX-maxMotifsPerSequence <INTEGER> (*)\n"
				"          Maximum number of motif occurrences per sequence.\n\n");
		printf( "      --XX-track <STRING> (*)\n"
				"          Track extensions and refinements of IUPAC pattern string.\n\n");
		printf( "      --XX-effectiveIUPACStates <INTEGER> (*)\n"
				"          Effective number of different states in single IUPAC extensions. The\n"
				"          default is 6.\n\n");
		printf( "      --XX-effectivePWMStates <INTEGER> (*)\n"
				"          Effective number of different states in single PWM columns. The\n"
				"          default is 10.\n\n");
		printf( "      --XX-gapOpening <NUMBER> (*)\n"
				"          Bit penalty for each gap opening.\n\n");
		printf( "      --XX-gapExtension <INTEGER> (*)\n"
				"          Bit penalty for each gap extension.\n\n");
	}
	printf( "      --XX-localization\n"
			"          Calculate p-values for positional clustering of motif occurrences in\n"
			"          positive sequences of equal length. Improves the sensitivity to find\n"
			"          weak but positioned motifs.\n\n" );
	printf( "      --XX-localizationRanking\n"
			"          Rank motifs according to localization statistics.\n\n" );
	printf( "      --XX-downstreamPositions <INTEGER>\n"
			"          Distance between the anchor position (e.g. the transcription start\n"
			"          site) and the last positive sequence nucleotide. Corrects motif\n"
			"          positions in result plots. The default is 0.\n\n");
	if( developerHelp ){
		printf( "      --XX-startPosEnrichedReg <INTEGER> (*)\n"
				"          Expected start position of region enriched for motif occurrences,\n"
				"          relative to anchor position (see --localization).\n\n");
		printf( "      --XX-endPosEnrichedReg <INTEGER> (*)\n"
				"          Expected end position of region enriched for motif occurrences,\n"
				"          relative to anchor position (see --localization).\n\n");
		printf( "      --XX-format FASTA|MFASTA (*)\n"
				"          Use conservation information from multiple sequence alignments in\n"
				"          FASTA format or provide positive sequences in FASTA format (default).\n\n");
		printf( "      --XX-maxMFASTASequences <INTEGER> (*)\n"
				"          Limit the number of sequences to use from multiple sequence\n"
				"          alignments. By default, all sequences are used.\n\n" );
		printf( "      --XX-debug (*)\n"
				"          Printout evolving PWMs during the PWM stage.\n\n");
	}
	printf( "      --XX-batch\n"
			"          Suppress progress bars.\n\n" );

	if( developerHelp ){
		printf( "  Options for weighting positive sequences (*)\n" );
		printf( "      --posSequenceIntensities <FILEPATH>\n"
				"          File with intensities for positive sequences (one per line) used to\n"
				"          weight sequences in the EM algorithm. The order of intensities must\n"
				"          conform to the order of positive sequences. Higher intensities\n"
				"          produce higher sequence weights.\n\n" );
		printf( "      --useIntensitiesToInitBMMs\n"
				"          Use intensities to initialize BMMs from weighted instances of XXmotif\n"
				"          PWMs.\n\n" );
		printf( "      --useRanks\n"
				"          Use intensity ranks instead of intensities to calculate weights.\n\n" );
		printf( "      --backgroundQuantile <FLOAT>\n"
				"          The quantile of intensities (or ranks) that defines the background\n"
				"          intensity (rank) used to translate intensities (ranks) into weights.\n"
				"          The weight of sequences with intensities (ranks) below (above) the\n"
				"          background intensity (rank) is set to zero. The default is 0.0.\n\n" );
		printf( "      --backgroundIntensity <FLOAT>\n"
				"          The intensity that defines the background intensity used to translate\n"
				"          intensities into weights. The weight of sequences with intensities\n"
				"          below the background intensity is set to zero. The default is the\n"
				"          minimum intensity of positive sequences.\n\n" );
		printf( "      --backgroundRank <INTEGER>\n"
				"          The rank that defines the background rank used to translate ranks\n"
				"          into weights. The weight of sequences with ranks above the background\n"
				"          rank is set to zero. The default is the maximum rank (i.e. the\n"
				"          number of positive sequences).\n\n" );

		printf( "  Options for weighting binding sites (*)\n" );
		printf( "      --bindingSiteIntensities <FILEPATH>\n"
				"          File with intensities for binding site sequences (one per line) used\n"
				"          to initialize BMMs from weighted binding sites. The order of\n"
				"          intensities must conform to the order of sequences in the binding\n"
				"          sites file. Higher intensities produce higher weights.\n\n" );
		printf( "      --useBindingSiteRanks\n"
				"          Use intensity ranks instead of intensities to calculate weights.\n\n" );
		printf( "      --bindingSiteBackgroundQuantile <FLOAT>\n"
				"          The quantile of intensities (or ranks) that defines the background\n"
				"          intensity (rank) used to translate intensities (ranks) into weights.\n"
				"          The weight of binding sites with intensities (ranks) below (above)\n"
				"          the background intensity (rank) is set to zero. The default is 0.0.\n\n" );
		printf( "      --bindingSiteBackgroundIntensity <FLOAT>\n"
				"          The intensity that defines the background intensity used to translate\n"
				"          intensities into weights. The weight of binding sites with\n"
				"          intensities below the background intensity is set to zero. The\n"
				"          default is the minimum intensity of binding sites.\n\n" );
		printf( "      --bindingSiteBackgroundRank <INTEGER>\n"
				"          The rank that defines the background rank used to translate ranks\n"
				"          into weights. The weight of binding sites with ranks above the\n"
				"          background rank is set to zero. The default is the maximum rank (i.e.\n"
				"          the number of binding sites).\n\n" );
	}

	printf( "  Options to score sequences\n" );
	printf( "      --scorePosSequenceSet\n"
			"          Score positive (training) sequences with optimized BMMs.\n\n" );
	printf( "      --scoreNegSequenceSet\n"
			"          Score background (training) sequences with optimized BMMs.\n\n" );
	printf( "      --scoreTestSequenceSet <FILEPATH> [<FILEPATH>...]\n"
			"          Score test sequences with optimized BMMs. Test sequences can be\n"
			"          provided in a single or multiple FASTA files.\n\n" );
	if( developerHelp ){
		printf( "      --evaluatePWMs (*)\n"
				"          Score sequences with XXmotif PWMs.\n\n" );
		printf( "      --calculateLogScores (*)\n"
				"          Calculate log instead of log-odds scores.\n\n" );
	}

	printf( "  Output options\n" );
	printf( "      --saveInitBMMs\n"
			"          Write initialized BMM(s) to disk.\n\n" );
	printf( "      --saveBMMs\n"
			"          Write optimized BMM(s) to disk.\n\n" );
	if( developerHelp ){
		printf( "      --saveEMLikelihoods (*)\n"
				"          Write sequence likelihoods and positional odds scores to disk after\n"
				"          each EM iteration.\n\n" );
		printf( "      --saveEMBMMs (*)\n"
				"          Write BBM(s) to disk after each EM iteration.\n\n" );
	}
	printf( "      --verbose\n"
			"          Verbose terminal printouts.\n\n" );
	printf( "      -h, --help\n"
			"          Printout this help.\n\n" );

	if( developerHelp ){
		printf( "      (*) Developer options\n\n" );
	}

	exit( -1 );
}

bool Global::readCommandLineOptions( int argc, char *argv[] ){

	/*
	 * 1. process flags & option-preceding arguments
	 *    * flags
	 *    * option-preceding arguments
	 * 2. process arguments without preceding option
	 *    * output directory
	 *    * input sequence file
	 * 3. process settings depending on arguments (in 2nd)
	 *
	 * order necessary to accurately use GlobalOption()
	 */

	GetOpt_pp ops( argc, argv );

	/*
	 *  mode setting
	 */

	if( ops >> OptionPresent( 'h', "help" ) ){
		return false;
	}
	ops >> OptionPresent( "XX-debug", DEBUG );
//	ops >> OptionPresent( "aa", aa );
	ops >> OptionPresent( "evaluatePWMs", evaluatePWMs );

	/*
	 * mode-specific default setting
	 */

	order = ( negFile == NULL ) ? 2 : 8; // background model order
	mergeMode = LOW;
	maxMatchPositions = 17;
	minMatchPositions = 5;

	/*
	 * process flags
	 */

	ops >> OptionPresent( "XX-batch", batch_mode );

	if( ops >> OptionPresent( "XX-noLengthOptimPWMs" ) ){
		maximizeMotifLength = false;
	}

	ops >> OptionPresent( "XX-MOPS", multipleOccurrence );
	ops >> OptionPresent( "XX-OOPS", oneOccurrence );
	ops >> OptionPresent( "XX-ZOOPS", zeroOrOneOccurrence );
	if( ( oneOccurrence && multipleOccurrence) || ( oneOccurrence && zeroOrOneOccurrence ) || ( zeroOrOneOccurrence && multipleOccurrence ) ){
		fprintf( stderr, "Error: Please choose at most one of the options --zoops, --mops, and --oops\n" );
		exit( -1 );
	}
	if( !( oneOccurrence ) && !( zeroOrOneOccurrence ) && !( multipleOccurrence ) ){
		zeroOrOneOccurrence = true;
	}

	ops >> OptionPresent( "reverseComp", revcomp );

	ops >> OptionPresent( "XX-localization", usePositionalProbs );
	if( usePositionalProbs && ( Global::posSet->max_leng != Global::posSet->min_leng ) ){
		fprintf( stderr, "Error: Option --localization can only be applied for positive sequences of equal length\n" );
		exit( -1 );
	}
	if( usePositionalProbs ){
		ops >> OptionPresent( "XX-localizationRanking", positionalProbsRanking );
	}

	if( em ){ // default
		ops >> OptionPresent( "msq", msq );
		if( ops >> OptionPresent( "nonBayesian" ) ){
			interpolate = false;
		}
		ops >> OptionPresent( "noEM", noExpectationMaximizationPhase );
		ops >> OptionPresent( "useIntensitiesToInitBMMs", initInts );
		ops >> OptionPresent( "saveInitBMMs", saveInitModels );
		ops >> OptionPresent( "saveBMMs", saveModels );
		ops >> OptionPresent( "verbose", verbose );
		ops >> OptionPresent( "learnAlpha", learnHyperParameter );
		ops >> OptionPresent( "posAlpha", positionSpecificAlphas );
		ops >> OptionPresent( "debugAlpha", debugAlphalearning );
	}

	if( em || evaluatePWMs ){
		ops >> OptionPresent( "scorePosSequenceSet", testPosSequences );
		ops >> OptionPresent( "scoreNegSequenceSet", testNegSequences );
		ops >> OptionPresent( "calculateLogScores", logProbs );
	}

	ops >> OptionPresent( "ali-free", useAliFree );
	ops >> OptionPresent( "extra-homology-filter", removeHomology );
	ops >> OptionPresent( "filtering", repeatFiltering );
	ops >> OptionPresent( "lcf", lowComplexityFilter );
	ops >> OptionPresent( "noRefinementPhase", noRefinementPhase );
	ops >> OptionPresent( "ranks", useRankPvalues );

	if( em ){
		ops >> OptionPresent( "saveEMLikelihoods",
							   saveExpectationMaximizationLikelihoods );
		ops >> OptionPresent( "saveEMBMMs",
							   saveExpectationMaximizationModels );
	}

	/*
	 * process option-preceding arguments
	 */

	ops >> Option( "negSequenceSet", negFile );
	ops >> Option( "XX-K", order );

	if( ops >> Option( "XX-pseudoCounts", pseudo ) ){
		pseudo /= 100;
	}

	ops >> Option( "XX-gaps", GAPS );

	ops >> Option( "XX-seeds", type );
	if( type == NO_VALID_MOTIF_TYPE ){
		return false;
	}

	ops >> Option( "XX-mergeMotifsThreshold", mergeMode );
	if( mergeMode == NO_VALID_MERGE_MODE ){
		return false;
	}

	if( ops >> Option( "XX-maxPositions", maxMatchPositions ) ){
		if( maxMatchPositions > 26 ){
			fprintf( stderr, "Error: The maximum number of motif positions in XXmotif is 26\n");
			exit( -1 );
		}
	}

	ops >> Option( "maxPosSequences", maxPosSetSize );

	std::string tr;
	if( ops >> Option( "XX-track", tr ) ){
	  trackedMotifs.insert( tr );
	}

	ops >> Option( "XX-format", seqFormat );
	if( seqFormat == NO_VALID_SEQ_FORMAT ){
		return false;
	}

	ops >> Option( "XX-maxMFASTASequences", maxMultipleSequences );

	ops >> Option( "XX-downstreamPositions", downstream );

	ops >> Option( "XX-jumpStartPatternStage", startMotif );
	ops >> Option( "XX-jumpStartPWMStage", profFile );

 	int offset = posSet->max_leng - downstream;
	if( ops >> Option( "XX-startPosEnrichedReg", startRegion ) ){
		startRegion += offset;
	} else{
		startRegion = 0;
	}
	if( ops >> Option( "XX-endPosEnrichedReg", endRegion ) ){
		Global::endRegion += offset;
	} else{
		endRegion = posSet->max_leng;
	}

	if( em ){

		ops >> Option( 'k', modelOrder );

		if( ops >> Option( "bindingSiteFile", bindingSiteFile ) ){
			if( ops >> OptionPresent( "bindingSiteLength" ) ){
				ops >> Option( "bindingSiteLength", bindingSiteLength );
			} else{

				/*
				 * determine the length of binding site sequences
				 */

				FILE* fp;
				if( ( fp = fopen( bindingSiteFile, "r" ) ) == NULL ){
			        fprintf( stderr, "Error: Cannot open binding sites file %s\n", bindingSiteFile );
			        exit( -1 );
				}

				int c;
				while( ( c = fgetc( fp ) ) != EOF && ( c == '\n' || c == '\r' ) ){
					; // skip leading blank lines
				}
				if( c == EOF ){
					fprintf( stderr, "Error: No binding site sequences in %s\n", bindingSiteFile );
					exit( -1 );
				}

				int L = -1;			// length of first binding site sequence
				int l = 1;			// length or current binding site sequence
				int ncounter = 0;	// counter for \n and \r

				while( ( c = fgetc( fp ) ) != EOF ){
					if( c == '\n' || c == '\r' ){
						if( ncounter > 0 ){
							continue;
						} else if( L == -1 ){
							// remember the length of the first binding site sequence
							L = l;
						} else if ( L != l ){
							fprintf( stderr, "Error: Please provide binding site sequences of equal length\n" );
							exit( -1 );
						}
						l = 0; // reset for next binding site sequence
						ncounter++;
					} else{
						l++;
						ncounter = 0; // reset counter for \n and \r
					}
				}

				fclose( fp );

				bindingSiteLength = L;
			}

			if( ops >> OptionPresent( "extend" ) ){
				addColumns.clear();
				ops >> Option( "extend", addColumns );
				if( addColumns.size() < 1 || addColumns.size() > 2 ){
					fprintf( stderr, "Error: Wrong format of option --extend\n" );
					exit( -1 );
				}
				if( addColumns.size() == 1 ){
					addColumns.resize( 2, addColumns.back() );
				}
				bindingSiteLength = addColumns.at( 0 ) + bindingSiteLength + addColumns.at( 1 );
			}

			if( bindingSiteLength > posSet->min_leng ){
				fprintf( stderr, "Error: Binding sites are longer than positive sequences\n" );
				exit( -1 );
			}

			if( !( modelOrder < bindingSiteLength ) ){
				modelOrder = bindingSiteLength - 1;
			}

			if( ops >> Option( "bindingSiteIntensities", bindingSiteIntsFile ) ){

				std::ifstream fs( bindingSiteIntsFile, std::ios_base::in );
				if( fs.fail() ){
			        fprintf( stderr, "Error: Cannot open file %s with binding site intensities\n", bindingSiteIntsFile );
					exit( -1 );
				}

				float weight;
				std::vector<float> weights;

				while( fs >> weight ){
					if( weight < 0.0f ){
						fprintf( stderr, "Error: Negative binding site intensities\n" );
						exit( -1 );
					}
					weights.push_back( weight );
				}
				fs.close();

				ops >> OptionPresent( "useBindingSiteRanks", bindingSiteRankWeighting );
				if( bindingSiteRankWeighting ){
					float N = static_cast<float>( weights.size() );
					ops >> Option( "bindingSiteBackgroundQuantile", bindingSiteBackgroundQuantile );
					if( !( ops >> Option( "bindingSiteBackgroundRank", bindingSiteBackgroundRank ) ) ){
						if( bindingSiteBackgroundQuantile > 1 - ( ( 1.0f / 2.0f ) / static_cast<float>( N ) ) ){
							bindingSiteBackgroundRank = 1.0f;
						} else{
							bindingSiteBackgroundRank = N - roundf( N * bindingSiteBackgroundQuantile );
						}
					}
				} else{

					ops >> Option( "bindingSiteBackgroundQuantile", bindingSiteBackgroundQuantile );
					if( !( ops >> Option( "bindingSiteBackgroundIntensity", bindingSiteBackgroundIntensity ) ) ){
						bindingSiteBackgroundIntensity = quantile( weights, bindingSiteBackgroundQuantile, 7 );
					}
				}
			}
		} else if( ops >> OptionPresent( "markovModelFile" ) ){
			ops >> Option( "markovModelFile", markovModelFile );

			/* determine model length and order */

			FILE* fp;
			std::stringstream str;

			str << markovModelFile << ".conds";

			if( ( fp = fopen( str.str().c_str(), "r" ) ) == NULL ){
		        fprintf( stderr, "Error: Cannot open file %s with BMM probabilities\n", str.str().c_str() );
		        exit(-1);
			}

			int c;
			while( ( c = fgetc( fp ) ) != EOF && ( c == '\n' || c == '\r' ) ){
				; // skip leading blank lines
			}
			if( c == EOF ){
				fprintf( stderr, "Error: Cannot find BMM probabilities in file %s\n", str.str().c_str() );
				exit( -1 );
			}

			int lines = 0; // lines without blank lines
			int positions = 0;

			int ncounter = 0; // counter for \n
			int rcounter = 0; // counter for \r
			while( ( c=fgetc( fp ) ) != EOF ){
				if( c == '\n' || c == 'r' ){
					if( c == 'n' ){
						if( ncounter == 0 && rcounter == 0 ){
							lines++;
						}
						ncounter++;
					} else{ // c == 'r'
						if( rcounter == 0 && ncounter == 0 ){
							lines++;
						}
						rcounter++;
					}
				} else{
					if( ncounter > 1 ){
						positions++;
					} else if( rcounter > 1 ){
						positions++;
					}
					ncounter = 0; // reset counter for \n
					rcounter = 0; // reset counter for \r
				}
			}
			fclose( fp );

			if( ncounter == 0 && rcounter == 0 ){
				lines++; // handle last line without newline in kate
			}

			markovModelLength = positions + 1;
			modelOrder = ( lines / markovModelLength ) - 1;

			if( !( modelOrder < markovModelLength ) ){
				exit( -1 );
			}

			if( ops >> OptionPresent( "extend" ) ){
				addColumns.clear();
				ops >> Option( "extend", addColumns );
				if( addColumns.size() < 1 || addColumns.size() > 2 ){
					fprintf( stderr, "Error: Wrong format of option --extend\n" );
					exit( -1 );
				}
				if( addColumns.size() == 1 ){
					addColumns.resize( 2, addColumns.back() );
				}
				markovModelLength = addColumns.at( 0 ) + markovModelLength + addColumns.at( 1 );
			}

			if( markovModelLength > posSet->min_leng ){
				fprintf( stderr, "Error: The BMM is longer than the positive sequences\n" );
				exit( -1 );
			}

		} else{

			if( !( modelOrder < PWM_LENGTH ) ){
				modelOrder = PWM_LENGTH - 1;
			}

			// number of one or more models in the ranking to pursue
			if( ops >> OptionPresent( "rankPWMs" ) ){
				ops >> Option( "rankPWMs", nrModels );
				std::sort( nrModels.begin(), nrModels.end() );
			}
			// min. number of models to pursue
			ops >> Option( "minPWMs", minModels );
			// max. number of models to pursue
			ops >> Option( "maxPWMs", maxModels );
			// max. p-value of models
			ops >> Option( "maxPValue", pValueThreshold );
			pValueThreshold = log( pValueThreshold );
			// min. percentage of sequences containing a binding site instance
			ops >> Option( "minOccurrence", minOccurrence );
			// add columns (left, right) to models
			if( ops >> OptionPresent( "extend" ) ){
				addColumns.clear();
				ops >> Option( "extend", addColumns );
				if( addColumns.size() < 1 || addColumns.size() > 2 ){
					fprintf( stderr, "Error: Wrong format of option --extend\n" );
					exit( -1 );
				}
				if( addColumns.size() == 1 ){
					addColumns.resize( 2, addColumns.back() );
				}
			}
		}

		if( ops >> Option( "posSequenceIntensities", sequenceIntsFile ) ){
			std::ifstream fs( sequenceIntsFile, std::ios_base::in );
			if( fs.fail() ){
				fprintf( stderr, "Error: Cannot open file %s with positive sequence intensities\n", Global::sequenceIntsFile );
				exit( -1 );
			}

			float weight;
			std::vector<float> weights;

			while( fs >> weight ){
				if( weight < 0.0f ){
					fprintf( stderr, "Error: negative positive sequence intensities\n" );
					exit( -1 );
				}
				weights.push_back( weight );
			}
			fs.close();

			if( Global::posSet->nent != static_cast<int>( weights.size() ) ){
				fprintf( stderr, "Error: Differing number of positive sequences and positive sequence intensities\n" );
				exit( -1 );
			}

			if( verbose ){
				std::cout << std::endl;
				std::cout << " _____________" << std::endl;
				std::cout << "|             |" << std::endl;
				std::cout << "| INTENSITIES |" << std::endl;
				std::cout << "|_____________|" << std::endl;
				std::cout << std::endl;
				for( std::vector<float>::iterator iter = weights.begin(); iter <
				     weights.end(); ++iter ){
					std::cout << *iter << " ";
				}
				std::cout << std::endl;
			}

			ops >> OptionPresent( "useRanks", rankWeighting );

			// calculate weights
			if( rankWeighting ){

				float N = static_cast<float>( weights.size() );

				ops >> Option( "backgroundQuantile", backgroundQuantile );

				if( !( ops >> Option( "backgroundRank", backgroundRank ) ) ){

					if( backgroundQuantile > 1 - ( ( 1.0f / 2.0f ) /
						static_cast<float>( N ) ) ){
						backgroundRank = 1.0f;
					} else{
						backgroundRank = N - roundf( N * backgroundQuantile );
					}
				}

				// calculate rank-based weights
				calculateWeights( weights, backgroundRank, true );
			} else{

				ops >> Option( "backgroundQuantile", backgroundQuantile );

				if( !( ops >> Option( "backgroundIntensity", backgroundIntensity
					) ) ){

					backgroundIntensity = quantile( weights, backgroundQuantile,
							                        7 );
				}
				// calculate intensity-based weights
				calculateWeights( weights, backgroundIntensity, false );
			}

			if( verbose ){
				std::cout << std::endl;
				std::cout << " _________" << std::endl;
				std::cout << "|         |" << std::endl;
				std::cout << "| WEIGHTS |" << std::endl;
				std::cout << "|_________|" << std::endl;
				std::cout << std::endl;
				for( std::vector<float>::iterator iter = weights.begin(); iter <
				     weights.end(); ++iter ){
					std::cout << *iter << " ";
				}
				std::cout << std::endl;
			}

			// assign sequence weights to sequences
			for( unsigned i=0; i < weights.size(); i++ ){
				Global::posSet->entity[i+1]->weight = weights.at(i);
			}
		} else{
			// assign default sequence weights to sequences
			for( int i=1; i <= Global::posSet->nent; i++ ){
				Global::posSet->entity[i]->weight = 1.0f;
			}
		}

		if( ops >> OptionPresent( 'q' ) ){
			ops >> Option( 'q', q );
			if( q <= 0 || q >= 1 ){
				fprintf( stderr, "Error: The value of option -q is restricted to ]0,1[\n" );
				exit( -1 );
			}
		}
		if( ops >> OptionPresent( 'a', "alpha" ) ){
			alpha.clear();
			ops >> Option( 'a', "alpha", alpha );
			if( static_cast<int>( alpha.size() ) != modelOrder+1 ){
				if( static_cast<int>( alpha.size() ) > modelOrder+1 ){
					alpha.resize( modelOrder+1 );
				} else{
					alpha.resize( modelOrder+1, alpha.back() );
				}
			}
			float b;
			if( ops >> OptionPresent( 'b', "beta" ) ){
				ops >> Option( 'b', "beta", b );
				fprintf( stderr, "Option -b is ignored.\n" );
			}
			float g;
			if( ops >> OptionPresent( 'g', "gamma" ) ){
				ops >> Option( 'g', "gamma", g );
				fprintf( stderr, "Option -g is ignored.\n" );
			}
		} else{
			if( ops >> OptionPresent( 'b', "beta" ) ){
				ops >> Option( 'b', "beta", beta );
			}
			if( ops >> OptionPresent( 'g', "gamma" ) ){
				ops >> Option( 'g', "gamma", gamma );
			}
			if( modelOrder > 0 ){
				for( unsigned int i=1; i <= alpha.size(); i++ ){
					alpha[i] = beta * powf( gamma,
							   static_cast<float>( modelOrder )-1.0f );
				}
			}
		}

		ops >> Option( 'K', modelOrderBg );
		ops >> Option( 'A', "Alpha", alphaBg );

	} else if( evaluatePWMs ){

		if( ops >> OptionPresent( "rankPWMs" ) ){
			ops >> Option( "rankPWMs", nrModels );
			std::sort( nrModels.begin(), nrModels.end() );
		}
		ops >> Option( "minPWMs", minModels );
		ops >> Option( "maxPWMs", maxModels );
		ops >> Option( "maxPValue", pValueThreshold );
		pValueThreshold = log( pValueThreshold );
		ops >> Option( "minOccurrence", minOccurrence );
	}

	ops >> Option( "benchmarkFolder", benchmarkFolder );
	ops >> Option( "counts-offset", countsOffset );

	ops >> Option( "final-filter", finalFilterThreshold );
	finalFilterThreshold = log( finalFilterThreshold );

	ops >> Option( "XX-gapExtension", gapExtension );
	ops >> Option( "XX-gapOpening", gapOpening );

	std::string thresh_init;
	if( ops >> Option( "instance-threshold", thresh_init ) ) {
		instanceThreshold = ThresholdChecker( thresh_init );
	} else{
		instanceThreshold = ThresholdChecker( "<1.0" );
	}

	ops >> Option( "maxIterations", maxIterations );
	ops >> Option( "maxMotifLevel", maxMotifLevel );
	ops >> Option( "XX-maxMotifsPerSequence", maxMotifsPerSequence );
	ops >> Option( "XX-minMatchPositions", minMatchPositions );

	int n_eff = -1;
	neff_discrete = 6;
	neff_pwm = 10;
	if( ops >> Option( "neff", n_eff ) ){
		neff_pwm = n_eff;
		neff_discrete = n_eff;
	}
	ops >> Option( "XX-effectiveIUPACStates", neff_discrete );
	ops >> Option( "XX-effectivePWMStates", neff_pwm );

	if( ops >> Option( "plusFrac", plusFrac ) ){
		plusFrac /= 100;
	}

	ops >> Option( "XX-A", pseudocountsFactor );
	ops >> Option( "termFreqScale", dollarScaleFactor );

	ops >> Option( "maxEMIterations", maxEMIterations );


	ops >> Option( "cswlen", cswlen, 13 );
	ops >> Option( "csbest", csbest, 0 );
	ops >> Option( "csprofiles", csprofiles );
	if( csbest != 0 ){
		if( csprofiles.length() == 0 ){
			fprintf( stderr, "Error: CS profile output file is missing (use --csprofiles).\n" );
			exit( -1 );
		}
	}

	double min_cov;
	if( ops >> Option( "min-coverage", min_cov ) ){
		if( min_cov <= 1.0 ){
			minCoverage = min_cov * posSet->nent;
		} else{
			minCoverage = min_cov;
		}
	} else{
		minCoverage = 0;
	}

	if( minCoverage > posSet->nent ){
		minCoverage = posSet->nent;
	}

	if( em ){
		ops >> Option( 'e', "epsilon", epsilon );
	}

	if( em || evaluatePWMs ){
		if( ops >> Option( "scoreTestSequenceSet", testSequenceFile ) ){
			for( unsigned int i=0; i < testSequenceFile.size(); i++ ){
				testSet.push_back( readSeqSet( const_cast<char*>( testSequenceFile.at(i).c_str() ), Global::A, FASTA, INT_MAX ) );
				if( revcomp ){
					createRevcomp( testSet.back(), A );
				}
			}
		}
	}

	/*
	 * process arguments without preceding option
	 * * output directory
	 * * input sequence file
	 */

	std::vector<std::string> args;
	ops >> GlobalOption( args );
	if( !( args.size() == 2 ) ){
//		for( unsigned int i = 0; i < args.size(); i++ ){
//			printf( "%s\n", args.at(i).c_str() );
//		}
		return false;
	}

	outputDirectory = String( args[0].c_str() );
	createDirectory( outputDirectory );
	pwmFolder = String( outputDirectory );

	tmpDirectory = ( char* )calloc( 1024, sizeof( char ) );
	sprintf( tmpDirectory, "%s/tmp", outputDirectory );
	createDirectory( tmpDirectory );

	int i = 0, start = 0, end = 0;
	name = String( args[1].c_str() );
	while( name[++i] != '\0' ){
		if( name[i] == '.' ){
			end = i - 1;
		}
	}
	while( --i != 0 && name[i] != '/' );
	if( i == 0 ){
		start = 0;
	} else{
		start = i + 1;
	}
	shortFileName = ( char* )malloc( ( end-start+2 ) * sizeof( char ) );
	for( i=start; i <= end; i++ ){
		shortFileName[i - start] = name[i];
	}
	shortFileName[i - start] = '\0';

	/*
	 * process settings depending on arguments
	 * * output directory
	 * * input sequence file
	 */

	ops >> Option( "write-pwm-file", pwmFolder );

	std::stringstream s;
	s << Global::outputDirectory << "/" << Global::shortFileName << ".mtf";

	if( em && verbose ){

		std::cout << std::endl;
		std::cout << " ____________________" << std::endl;
		std::cout << "|                    |" << std::endl;
		std::cout << "| BaMM!motif setting |" << std::endl;
		std::cout << "|____________________|" << std::endl;
		std::cout << std::endl;
		std::cout << "sequence file" << "\t\t\t\t" << name << std::endl;
		std::cout << std::endl;
		if( bindingSiteFile == NULL && markovModelFile == NULL ){
			if( !( nrModels.empty() ) ){
				std::cout << "ranking models\t\t\t\t";
				for( std::vector<int>::iterator iter = nrModels.begin(); iter <
				     nrModels.end(); ++iter ){
					std::cout << *iter << " ";
				}
				std::cout << std::endl;
			}
			std::cout << "min. models\t\t\t\t" << minModels << std::endl;
			std::cout << "max. models\t\t\t\t" << maxModels << std::endl;
			std::cout << "max. p-value\t\t\t\t" << exp( pValueThreshold )
					  << std::endl;
			std::cout << "min. occurrence\t\t\t\t" << minOccurrence
					  << std::endl;
			std::cout << "add columns\t\t\t\t";
			for( std::vector<int>::iterator iter = addColumns.begin(); iter <
			     addColumns.end(); ++iter ){
				std::cout << *iter << " ";
			}
			std::cout << std::endl;
			std::cout << "msq" << "\t\t\t\t\t" << ( msq ? "true" : "false" )
					  << std::endl;
		} else if( bindingSiteFile != NULL ){
			std::cout << "binding site file" << "\t\t\t" << bindingSiteFile
					  << std::endl;
			std::cout << "binding site lengths" << "\t\t\t" << bindingSiteLength
				      << std::endl;
		} else{
			std::cout << "Markov model file" << "\t\t\t" << markovModelFile
					  << std::endl;
			std::cout << "Markov model size" << "\t\t\t" << markovModelLength
					  << std::endl;
		}
		std::cout << std::endl;
		std::cout << "order" << "\t\t\t\t\t" << modelOrder << std::endl;
		std::cout << "pseudo-counts factor(s)" << "\t\t\t";
		for( std::vector<float>::iterator iter = alpha.begin(); iter <
		     alpha.end(); ++iter ){
			std::cout << *iter << " ";
		}
		std::cout << std::endl;
		std::cout << "interpolate\t\t\t\t" << ( interpolate ? "true" : "false" )
				  << std::endl;
		std::cout << std::endl;
		if( negFile != NULL ){
			std::cout << "background sequence file" << "\t\t" << negFile
					  << std::endl;
		}
		std::cout << "background order" << "\t\t\t" << modelOrderBg
				  << std::endl;
		std::cout << "background pseudo-counts factor" << "\t\t" << alphaBg
				  << std::endl;
		std::cout << std::endl;
		std::cout << "EM phase" << "\t\t\t\t"
				  << ( noExpectationMaximizationPhase ? "false" : "true" )
				  << std::endl;
		std::cout << "specificity factor q" << "\t\t\t" << q << std::endl;
		std::cout << "convergence e" << "\t\t\t\t" << epsilon << std::endl;
		std::cout << std::endl;
		if( sequenceIntsFile != NULL ){
			std::cout << "sequence intensities file" << "\t\t"
					  << sequenceIntsFile << std::endl;
			if( bindingSiteFile == NULL && markovModelFile == NULL ){
				std::cout << "use intensities to initialize models" << "\t"
						  << ( initInts ? "true" : "false" ) << std::endl;
			}
		}
		if( sequenceIntsFile ){
			std::cout << "sequence weighting" << "\t\t\t"
					  << ( rankWeighting ? "rank" : "intensity" )
					  << "-based weighting" << std::endl;
			if( rankWeighting ){
				std::cout << "background rank" << "\t\t\t\t" << backgroundRank
						  << std::endl;
			} else{
				std::cout << "background intensity" << "\t\t\t"
						  << backgroundIntensity << std::endl;
			}
		}
		std::cout << std::endl;
		if( bindingSiteFile != NULL && bindingSiteIntsFile != NULL ){
			std::cout << "binding site intensities file" << "\t\t"
					  << bindingSiteIntsFile << std::endl;
			std::cout << "binding site weighting" << "\t\t\t"
					  << ( bindingSiteRankWeighting ? "rank" : "intensity" )
					  << "-based weighting" << std::endl;
			if( bindingSiteRankWeighting ){
				std::cout << "binding site background rank" << "\t\t"
						  << bindingSiteBackgroundRank << std::endl;
			} else{
				std::cout << "binding site background intensity" << "\t"
						  << bindingSiteBackgroundIntensity << std::endl;
			}
		}

		std::cout << std::endl;
		std::cout << "score training sequences" << "\t\t" << ( testPosSequences ?
				     "true" : "false" ) << std::endl;
		std::cout << "score background sequences" << "\t\t"
				  << ( testNegSequences ? "true" : "false" ) << std::endl;
		if( !( testSet.empty() ) ){
			std::cout << "test sequence file" << "\t\t\t"
					  << testSequenceFile.at(0) << std::endl;
			if( testSet.size() > 1 ){
				for( unsigned int i=1; i < testSet.size(); i++ ){
					std::cout << "\t\t\t\t\t" << testSequenceFile.at(i)
							  << std::endl;
				}
			}
		}
		std::cout << "evaluate PWMs" << "\t\t\t\t" << ( evaluatePWMs ? "true" :
				     "false" ) << std::endl;
		std::cout << "calculate log likelihood ratios" << "\t\t" << ( logProbs ?
				     "false" : "true" ) << std::endl;
	}

	if( ops.options_remain() ){
		std::cerr << "Unknown options remaining..." << std::endl;
		return false;
	}

	return true;
}


Global::~Global(){
   	if( name != NULL ){
   		free(name);
   	}
   	if( tmpDirectory != NULL ){
   		deleteDirectory( tmpDirectory );
   		free( tmpDirectory );
   	}
   	if( outputDirectory != NULL ){
   		free( outputDirectory );
   	}
   	if( shortFileName != NULL ){
   		free( shortFileName );
   	}
   	if( negFile != NULL ){
   		free( negFile );
   	}
   	if( startMotif != NULL ){
   		free( startMotif );
   	}
   	if( profFile != NULL ){
   		free(profFile);
   	}
   	if( benchmarkFolder != NULL ){
   		free( benchmarkFolder );
   	}
   	if( pwmFolder != NULL ){
   		free( pwmFolder );
   	}

	ss_type set = posSet;
	if( negSet != NULL ){
		set = negSet;
	}

	double POW_2_16 = pow( 2, 16 );
	if( conservationProbs != NULL ){
		for( int i = 0; i <= POW_2_16; i++ ){
			if( conservationProbs[i] == NULL ){
				continue;
			}
			for( int j = 1; j < set->max_MultSeq; j++ ){
				free( conservationProbs[i][j] );
			}
			free( conservationProbs[i] );
		}
		free( conservationProbs );
	}

	if( alignmentFreeProbs != NULL ){
		for( int i = 0; i <= POW_2_16; i++ ){
			if( alignmentFreeProbs[i] == NULL ){
				continue;
			}
			for( int j = 1; j < set->max_MultSeq; j++ ){
				free( alignmentFreeProbs[i][j] );
			}
			free( alignmentFreeProbs[i] );
		}
		free( alignmentFreeProbs );
	}

	free( posBg_log );
   	free( negBg_log );
   	free( posBg );
   	free( negBg );

	NullModel::destruct();

   	if( negSet != NULL ){
   		NilSeqSet( negSet );
   	}

   	NilSeqSet( posSet );
   	NilAlpha( A );

   	if( em ){

   		hoNullModel::destruct();

   		if( bindingSiteFile != NULL ){
   			free( bindingSiteFile );
   		}
   		if( markovModelFile != NULL ){
   			free( markovModelFile );
   		}
   		if( freqs != NULL ){
   			free( freqs );
   		}
   		if( !( testSet.empty() ) ){
			std::list<ss_type>::const_iterator iter;
			for( iter=Global::testSet.begin(); iter != Global::testSet.end(); iter++ ){
				NilSeqSet( *iter );
			}
   		}
   	}
}

std::ostream& operator<<( std::ostream &os, const merge_type &m ){
  switch( m ){
    case LOW: os << "LOW"; break;
    case MEDIUM: os << "MEDIUM"; break;
    case HIGH: os << "HIGH"; break;
    case NO_VALID_MERGE_MODE: os << "NO_VALID_MERGE_MODE"; break;
    default: os << "UNKNOWN - THIS MUST NOT HAPPEN!"; exit( -1 );
  }
  return os;
}

std::ostream& operator<<( std::ostream &os, const SuppInfMode &v ){
	switch( v ){
	case SUPP_NO: os << "NO"; break;
	case SUPP_DISOCONS: os << "DISOCONS"; break;
	case SUPP_NNET: os << "NNET"; break;
	default: os << "UNKNOWN - THIS MUST NOT HAPPEN!"; exit( -1 );
	}
	return os;
}
