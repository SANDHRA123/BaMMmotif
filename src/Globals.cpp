#include <fstream>
#include <limits.h>
#include <math.h>

#include "refinementPhase/Motif.h"
#include "backgroundDistribution.h"
#include "getopt_pp/getopt_pp.h"
#include "Globals.h"
#include "NullModel.h"
#include "em/hoUtils.h"

#include "em/hoNullModel.h" /* Expectation Maximization(EM)-specific */

using GetOpt::GetOpt_pp;
using GetOpt::GlobalOption;
using GetOpt::Option;
using GetOpt::OptionPresent;

/* TODO
 * improve sequence & readability */

a_type 		Global::A = NULL;
ss_type		Global::posSet = NULL;			/** positive sequence set **/
ss_type 	Global::negSet = NULL;			/** negative sequence set **/
bool		Global::usePositionalProbs = false;
bool		Global::positionalProbsRanking = false;

int 		Global::GAPS = 0;				/** number of gap combinations allowed in the start motifs */
merge_type	Global::mergeMode; /* default is set below for proteins/dna */
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

bool   		Global::multipleOccurrence = false;	/** multiple occurrences should be found in one sequence **/
bool   		Global::oneOccurrence = false;	/** multiple occurrences should be found in one sequence **/
bool		Global::zeroOrOneOccurrence = false;
bool		Global::revcomp= false;			/** search on reverse complement of sequences **/
bool 		Global::repeatFiltering = false;
bool 		Global::lowComplexityFilter = false;
bool		Global::noRefinementPhase = false;

char* 		Global::startMotif = NULL;		/** start motif for motif discovery **/
char*		Global::profFile = NULL;		/** file with a the start profile for motif discovery **/
int			Global::startRegion = 0;		/** start of enriched region **/
int			Global::endRegion = 1;			/** end of enriched region **/
motif_type	Global::type = ALL;			/* type of motif: ALL, FIVEMERS, PALINDROME, TANDEM, NOPALINDROME, NOTANDEM*/
seq_format	Global::seqFormat = FASTA;      /* format of input sequences: FASTA, CLUSTALW */

float***	Global::conservationProbs = NULL;
float***	Global::alignmentFreeProbs = NULL;
bool        Global::removeHomology = false;

double* 	Global::posBg_log = NULL;		/* logarithm of distribution of positive set **/
double* 	Global::posBg = NULL;			/* background of positive set **/
double* 	Global::negBg_log = NULL; 		/* logarithm of distribution of negative set **/
double*		Global::negBg = NULL;			/* background of negative set **/

double		Global::pseudo = 0.1;			/** pseudocounts for pwm for finding new instances**/
double		Global::plusFrac = 0;
int			Global::neff_pwm;				/** effective number of different bases in one pwm column **/
int 		Global::neff_discrete;			/** effective number of different bases before pwm phase **/

int			Global::downstream = 0;			/** distance between the alignment point and the end of the input sequences **/

char* 		Global::outputDirectory = NULL; /** output Directory **/
char*		Global::name = NULL;			/** input file name **/
char*		Global::shortFileName = NULL;	/** input file name without path and .**/
char*		Global::negFile = NULL;			/** negative input file name **/
char* 		Global::benchmarkFolder = NULL; /** folder in which the benchmark results are stored **/
char*		Global::pwmFolder = NULL;
int			Global::maxMotifLevel = 3;		/** max number of extensions tried per level */
double		Global::minCoverage;

int			Global::minMatchPositions;
int			Global::maxMatchPositions;

bool		Global::maximizeMotifLength = true;

StateType 	Global::type_of_states;


/* cs blast */
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

/* background model order */
int Global::order = 2;
/* background model pseudo-counts factor */
float Global::pseudocountsFactor = 10.0f;
/* background model counts offset */
float Global::countsOffset = 0.0f;

/*
 * Expectiation Maximization (EM) parameters
 * */

/*
 * higher-order modeling and EM mode
 */
bool Global::em = false;

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
 * min. number of XXmotif models used to initialize Markov models. Independent
 * on pValueThreshold and minOccurrence
 */
int Global::minModels = 0;
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
 * min. percentage of sequences containing a binding site instance. Not applied
 * to min. number of models
 */
float Global::minOccurrence = 0.0f;
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
int Global::modelOrder = 0;
/*
 * Markov model pseudo-counts factor(s). Markov model order k fixes vector size
 * to k+1
 */
std::vector<float> Global::alpha( modelOrder+1, 10.0f );

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
bool Global::interpolate = false;

/*
 * Interpolated Markov background model parameters
 * */

/*
 * Background model order
 */
int Global::modelOrderBg = 0;
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
float Global::q = 0.1f;
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

/* TODO
 * improve readability & (remove) comment(s) */
Global::Global(int argc, char *argv[]){
	argv0 = std::string(argv[0]);

	if(argc < 2) printHelpOutput();
	seq_format format = FASTA;
	for(long i = 2; i < argc; i++){
		//if(strcmp(argv[i], "--aa") == 0){
		//	A = MkAlpha("XACDEFGHIKLMNPQRSTVWY$");
		//	aa = true;
		//}else
		if(strcmp(argv[i], "--terminus-mode") == 0){
			if(i+1 < argc && strcmp(argv[i+1], "NONE") == 0) termMode = NONE;
			else if(i+1 < argc && strcmp(argv[i+1], "POS") == 0) termMode = POS;
			else if(i+1 < argc && strcmp(argv[i+1], "NEG") == 0) termMode = NEG;
			else if(i+1 < argc && strcmp(argv[i+1], "BOTH") == 0) termMode = BOTH;
			else { fprintf(stderr, "Illegal terminus mode\n"); exit(1); }
		}else if(strcmp(argv[i], "--maxPosSetSize") == 0){
			if(i+1 == argc){ fprintf(stderr, "\n\nERROR: no value given for option --maxPosSetSize\n\n"); exit(-1); }
			maxPosSetSize = atoi(argv[i+1]);
		}else if(strcmp(argv[i], "--format") == 0){
			if(i+1 == argc){ fprintf(stderr, "\n\nERROR: no value given for option --format\n\n"); exit(-1); }
			if(i+1 < argc && strcmp(argv[i+1], "CLUSTALW") == 0) format = CLUSTALW;
			if(i+1 < argc && strcmp(argv[i+1], "MFASTA") == 0) format = MFASTA;
			if(i+1 < argc && strcmp(argv[i+1], "CUSTOM") == 0) format = CUSTOM;
		}else if(strcmp(argv[i], "--negSet") == 0){
			negFile = argv[i+1];
		}else if(strcmp(argv[i], "--lcf") == 0){
			lowComplexityFilter = true;
		}
	}

	struct stat sts;
	if (((stat (argv[2], &sts)) == -1) || S_ISDIR(sts.st_mode)){
		fprintf(stderr, "\n !!! SEQFILE %s does not exist !!! \n\n", argv[2]);
		exit(-1);
	}
	if (negFile != NULL && (((stat (negFile, &sts)) == -1) || S_ISDIR(sts.st_mode))){
		fprintf(stderr, "\n !!! NEGFILE %s does not exist !!! \n\n", negFile);
		exit(-1);
	}

	//if (!aa) termMode=NONE;
	termMode=NONE;

	if(A==NULL) A = MkAlpha("NACGT");

	posSet = readSeqSet(argv[2], A, format, maxPosSetSize, termMode==BOTH || termMode==POS);
	calculateSequenceFeatures(posSet, A);

    posBg_log = (double*)calloc(nAlpha(A)+1, sizeof(double));
    posBg     = (double*)calloc(nAlpha(A)+1, sizeof(double));
    negBg_log = (double*)calloc(nAlpha(A)+1, sizeof(double));
    negBg 	  = (double*)calloc(nAlpha(A)+1, sizeof(double));

	startRegion = 1;
	endRegion = Global::posSet->max_leng;

	if(!readCommandLineOptions(argc, argv))	printHelpOutput();
			if(negFile != NULL){
			negSet = readSeqSet(negFile,A,format, INT_MAX);
		}
		if(negSet != NULL ){
			if(repeatFiltering){ filter_repeats(negSet, A); }
			if(lowComplexityFilter){ cerr << "negSet: "; filter_lowComplexity(negSet, A); }
			if(revcomp){ createRevcomp(negSet, A); }
			fillGapsInMultipleAlignment(negSet, A);
			filterMaxMultipleSequences(negSet, maxMultipleSequences);
			calculateSequenceFeatures(negSet, A);
		}

		if(repeatFiltering){ filter_repeats(posSet, A); }
		if(lowComplexityFilter){ cerr << "posSet: "; filter_lowComplexity(posSet, A); }
		if(revcomp){ createRevcomp(posSet, A); }
		fillGapsInMultipleAlignment(posSet, A);
		filterMaxMultipleSequences(posSet, maxMultipleSequences);


	checkSequenceSet(); // do some tests whether the negative set is a good choice
	setBackgroundDistribution(); // get trinuc distribution and conditional trinuc distribution of the background set
}


/* TODO
 * improve readability in non-EM parts */
void Global::printHelpOutput(){

	bool developerHelp = true; /* use emHelp & emDeveloperHelp to enable EM-specific help */

	printf( "BaMM!motif version 1.0\n" );
	printf( "\n" );
	printf( "SYNOPSIS\n" );
	printf( "      BaMMmotif DIR FILE [OPTIONS]\n\n" );
	printf( "DESCRIPTION\n" );
	printf( "      Bayesian Markov model motif discovery software.\n\n" );
	printf( "      DIR\n"
			"          Output directory for all results.\n\n" );
	printf( "      FILE\n"
			"          FASTA file with positive sequences of equal length.\n\n" );

	printf("OPTIONS\n");
	printf( "      --negSequenceSet <FILE>\n"
			"          FASTA file with negative/background sequences.\n\n" );
	printf( "      --reverseComp\n"
			"          Search motifs on both positive sequences and their reverse\n"
			"          complements. This option is e.g. recommended when using positive\n"
			"          sequences derived from ChIP-seq experiments.\n\n" );

	printf( "  Options to initialize BMMs from file\n" );
	printf( "      --bindingSiteFile <FILE>\n"
			"          File with binding sites of equal length (one per line).\n\n" );
	printf( "      --markovModelFile <FILE>\n"
			"          File with BMM probabilities (<FILE>.conds and <FILE>.probs).\n\n" );

	printf( "  Options to initialize BMMs from XXmotif PWMs\n" );
	printf( "      --PWMRank <INTEGER> [<INTEGER>...]\n"
			"          Rank(s) of PWM(s) in XXmotif results. Ignores the remaining options\n"
			"          to initialize BMMs from PWMs.\n\n" );
	printf( "      --minPWMNumber <INTEGER>\n"
			"          Minimum number of PWMs. Ignores the options --maxPvalue and\n"
			"          --minOccurrence. The default is 1.\n\n" );
	printf( "      --maxPWMNumber <INTEGER>\n"
			"          Maximum number of PWMs.\n\n" );
	printf( "      --maxPValue <FLOAT>\n"
			"          Maximum p-value of PWMs. This filter is not applied to the top\n"
			"          minimum number of PWMs (--minPWMNumber). The default is 1.0.\n\n" );
	printf( "      --minOccurrence <FLOAT>\n"
			"          Minimum fraction of sequences that contain a PWM instance. This\n"
			"          filter is not applied to the top minimum number of PWMs\n"
			"          (--minPWMNumber). The default is 0.05.\n\n" );

	if( developerHelp ){
		printf( "  Options to initialize BMMs from XXmotif PWMs or binding site files\n" );
		printf( "      --msq (*)\n"
				"          Calculate BMM-specific q value from PWM-specific fraction of\n"
				"          sequences that contain a corresponding PWM instance.\n\n" );
	}

	printf( "  Options for (inhomogeneous) motif BMMs\n" );
	printf( "      -k <INTEGER>\n"
			"          Order. The default is 2.\n\n" );
	printf( "      -a|--alpha <FLOAT> [<FLOAT>...]\n"
			"          Order-specific prior strength. The default is 1 (for k = 0) and 20 x\n"
			"          3^(k-1) (for k > 0). Ignores the options -b/--beta and -g/--gamma.\n\n" );
	printf( "      -b|--beta <FLOAT>\n"
			"          Recalculate order-specific alphas according to beta x gamma^(k-1)\n"
			"          (for k > 0). The default is 20.\n\n" );
	printf( "      -g|--gamma <FLOAT>\n"
			"          Recalculate order-specific alphas according to beta x gamma^(k-1)\n"
			"          (for k > 0). The default is 3.\n\n" );
	printf( "      --extend <INTEGER>{1,2}\n"
			"          Extend BMMs by adding uniformly initialized positions to the left\n"
			"          and/or right of initial BMMs. Invoking e.g. with --extend 0 2 adds\n"
			"          two positions to the right of initial BMMs. By default, BMMs are not\n"
			"          extended.\n\n" );
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
			"          Prior probability for a sequence to contain a motif. The default is\n"
			"          0.9.\n\n" );
	printf( "      -e|--epsilon <FLOAT>\n"
			"          The EM algorithm is deemed to be converged when the sum over the\n"
			"          absolute differences in BMM parameters from successive EM rounds is\n"
			"          smaller than epsilon. The default is 0.001\n\n" );
	if( developerHelp ){
		printf( "      --maxEMIterations <INTEGER> (*)\n"
				"          Maximum number of EM iterations.\n\n" );
		printf( "      --noEM (*)\n"
				"          Initialize BMMs only.\n\n" );
	}

	printf( "  XXmotif options\n" );
	printf( "      --XX-ZOOPS\n"
			"          Use the zero-or-one-occurrence-per-sequence model. This is the\n"
			"          default.\n\n" );
	printf( "      --XX-MOPS\n"
			"          Use the multiple-occurrence-per-sequence model.\n\n" );
	printf( "      --XX-OOPS\n"
			"          Use the one-occurrence-per-sequence model.\n\n" );
	printf( "      --XX-K <INTEGER>\n"
			"          Order of the (homogeneous) background BMM. The default is 2 (when\n"
			"          learned on positive) or 8 (when learned on background sequences).\n\n" );
	printf( "      --XX-A|--XX-Alpha <FLOAT>\n"
			"          Prior strength. The default is 10.0.\n\n");
	printf( "      --XX-seeds ALL|FIVEMERS|PALINDROME|TANDEM|NOPALINDROME|NOTANDEM\n"
			"          Define the nature of seeds. The default is to start using ALL seed\n"
			"          variants.\n\n");
	printf( "      --XX-gaps 0|1|2|3\n"
			"          Maximum number of gaps used for start seeds. The default is 0.\n\n");
	printf( "      --XX-mergeMotifsThreshold LOW|MEDIUM|HIGH\n"
			"          Define the similarity threshold for merging PWMs in the refinement\n"
			"          phase. The default is to merge PWMs with LOW similarity in order to\n"
			"          reduce runtime.\n\n");
	if( developerHelp ){
		printf( "      --XX-minMatchPositions <INTEGER>\n"
				"          Minimum number of non-wildcard motif positions. The default is 4.\n\n" );
		printf( "      --XX-maxPositions <INTEGER>\n"
				"          Maximum number of motif positions. The default is 17.\n\n" );
		printf( "      --XX-maxOccurrencePerSequence\n"
				"          Maximum number of motif occurrences per sequence.\n\n");
		printf( "      --XX-maxNumberPosSequences <INTEGER>\n"
				"          Maximum number of sequences to use from positive sequences. The\n"
				"          default is to use all sequences.\n\n");
		printf( "      --XX-track <SEED>\n"
				"          Track extensions and refinements of IUPAC <SEED>.\n\n");
		printf( "      --effectiveIUPACStateNumber <INTEGER>\n"
				"          Effective number of different states in a single IUPAC extension. The\n"
				"          default is 6.\n\n");
		printf( "      --effectivePWMStateNumber <INTEGER>\n"
				"          Effective number of different states in a single PWM column. The\n"
				"          default is 10.\n\n");
		printf( "      --gapOpening <NUMBER>\n"
				"          Bit penalty for each gap opening\n\n");
		printf( "      --gapExtension <INTEGER>\n"
				"          Bit penalty for each gap extension.\n\n");
		printf( "      --cons-length <INTEGER>\n"
				"          Number of nucleotides used to calculate conservation p-values. The\n"
				"          default is 8.\n\n");
	}
	printf( "      --format FASTA|MFASTA\n"
			"          Use conservation information from multiple sequence\n"
			"          alignments in MFASTA format or provide positive sequences in FASTA\n"
			"          format (default).\n\n");
	printf( "      --maxMFASTASequences <INTEGER>\n"
			"          Maximum number of sequences to use from multiple sequence alignments.\n"
			"          By default, all sequences are used.\n\n" );
	printf( "      --localization\n"
			"          Calculate p-values for positional clustering of motif occurrences in\n"
			"          positive sequences of equal length. Improves the sensitivity to find\n"
			"          weak, positioned motifs.\n\n" );
	printf( "      --localizationRanking\n"
			"          Rank motifs according to localization statistics.\n\n" );
	printf( "      --downstreamPositions <INTEGER>\n"
			"          Distance between the anchor position (e.g. the transcription start\n"
			"          site) and the last positive sequence nucleotide. Corrects motif\n"
			"          positions in result plots. The default is 0.\n\n");
	printf( "      -m|--startMotif <MOTIF>\t\t\tStart motif (IUPAC characters)\n");
	printf( "      -p|--startProfile <FILE>\t\t\tprofile file\n");
	printf( "      --startRegion <NUMBER>\t\t\texpected start position for motif occurrences relative to anchor point (--localization)\n");
	printf( "      --endRegion <NUMBER>\t\t\texpected end position for motif occurrences relative to anchor point (--localization)\n");
	if( developerHelp ){
		printf( "      --XX-debug\n"
				"          Show PWMs during refinement phase.\n\n");
	}
	printf( "      --XX-batch\n"
			"          Suppress progress bars\n\n" );

	if( developerHelp ){
		printf( "  Options for weighting positive sequences(*)\n" );
		printf( "      --posSequenceIntensities <FILE>\n"
				"          File with intensities for positive sequences (one per line) used to\n"
				"          weight sequences in the EM algorithm. The order of intensities must\n"
				"          conform to the order of positive sequences. Higher intensities\n"
				"          produce higher sequence weights.\n\n" );
		printf( "      --useIntensitiesToInitBMMs\n"
				"          Use intensities to initialize BMMs from weighted instances of XXmotif\n"
				"          PWMs.\n\n" );
		printf( "      --useRanks\n"
				"          Use intensity ranks instead of intensities to calculate weights.\n\n" );
		printf( "      --quantileToBg <FLOAT>\n"
				"          The quantile of intensities (or ranks) that defines the background\n"
				"          intensity (rank) used to translate intensities (ranks) into weights.\n"
				"          The weight of sequences with intensities (ranks) below (above) the\n"
				"          background intensity (rank) is set to zero. The default is 0.0.\n\n" );
		printf( "      --intensityToBg <FLOAT>\n"
				"          The intensity that defines the background intensity used to translate\n"
				"          intensities into weights. The weight of sequences with intensities\n"
				"          below the background intensity is set to zero. The default is the\n"
				"          minimum intensity of positive sequences.\n\n" );
		printf( "      --rankToBg <INTEGER>\n"
				"          The rank that defines the background rank used to translate ranks\n"
				"          into weights. The weight of sequences with ranks above the background\n"
				"          rank is set to zero. The default is the maximum rank (i.e. the\n"
				"          number) of positive sequences.\n\n" );

		printf( "  Options for weighting binding sites (*)\n" );
		printf( "      --bindingSiteIntensities <FILE>\n"
				"          File with intensities for binding site sequences (one per line) used\n"
				"          to initialize BMMs from weighted binding sites. The order of\n"
				"          intensities must conform to the order of sequences in the binding\n"
				"          site file. Higher intensities produce higher weights.\n\n" );
		printf( "      --useBindingSiteRanks\n"
				"          Use intensity ranks instead of intensities to calculate weights.\n\n" );
		printf( "      --quantileToBindingSiteBg <FLOAT>\n"
				"          The quantile of intensities (or ranks) that defines the background\n"
				"          intensity (rank) used to translate intensities (ranks) into weights.\n"
				"          The weight of binding sites with intensities (ranks) below (above)\n"
				"          the background intensity (rank) is set to zero. The default is 0.0.\n\n" );
		printf( "      --intensityToBindingSiteBg <FLOAT>\n"
				"          The intensity that defines the background intensity used to translate\n"
				"          intensities into weights. The weight of binding sites with\n"
				"          intensities below the background intensity is set to zero. The\n"
				"          default is the minimum intensity of binding sites.\n\n" );
		printf( "      --rankToBindingSiteBg <INTEGER>\n"
				"          The rank that defines the background rank used to translate ranks\n"
				"          into weights. The weight of binding sites with ranks above the\n"
				"          background rank is set to zero. The default is the maximum rank (i.e.\n"
				"          the number) of binding sites.\n\n" );
	}

	printf( "  Options to score sequences\n" );
	printf( "      --scorePosSequenceSet\n"
			"          Score positive (training) sequences with optimized BMMs.\n\n" );
	printf( "      --scoreNegSequenceSet\n"
			"          Score background (training) sequences with optimized BMMs.\n\n" );
	printf( "      --scoreTestSequenceSet <FILE> [<FILE>...]\n"
			"          Score test sequences with optimized BMMs. Test sequences can be\n"
			"          provided in a single or multiple FASTA files.\n\n" );
	if( developerHelp ){
		printf( "      --evaluatePWMs (*)\n"
				"          Score sequences with XXmotif PWMs\n\n." );
		printf( "      --calculateLogScores (*)\n"
				"          Calculate log instead of log-odds scores.\n\n" );
	}

	printf( "  Output options\n" );
	printf( "      --saveInitBMMss\n"
			"          Write initialized BMM(s) to disk.\n\n" );
	printf( "      --saveBMMs\n"
			"          Write optimized BMM(s) to disk.\n\n" );
	printf( "      --verbose\n"
			"          Verbose console printouts.\n\n" );
	printf( "      -h, --help\n"
			"          Printout this help and exist.\n\n" );
	if( developerHelp ){
		printf( "      --saveEMLikelihoods (*)\n"
				"          Write sequence likelihoods and positional odds scores to disk after\n"
				"          each EM iteration.\n\n" );
		printf( "      --saveEMBMMs (*)\n"
				"          Write BBM(s) to disk after each EM iteration.\n\n" );
	}

	if( developerHelp ){
		printf( "      (*) Developer options\n" );
	}

	exit(-1);
}


bool Global::readCommandLineOptions( int argc, char *argv[] ){

	/* 1. process flags & option-preceding arguments:
	 *    * flags
	 *      * user flags
	 *        * --aa- and/or --em-specific
	 *      * developer flags
	 *        * --aa- and/or --em-specific
	 *    * option-preceding arguments
	 *      * user arguments
	 *        * --aa- and/or --em-specific
	 *      * developer arguments
	 *        * --aa- and/or --em-specific
	 * 2. process arguments without preceding option:
	 *    * output directory
	 *    * input sequence file
	 * 3. process settings depending on arguments (in 2nd)
	 *
	 * order necessary to accurately use GlobalOption( args ) */

	GetOpt_pp ops( argc, argv );

	/* mode setting */

	if( ops >> OptionPresent( 'h', "help" ) ){
		return false;
	}
	ops >> OptionPresent( "debug", DEBUG );
	//ops >> OptionPresent( "aa", aa );
	ops >> OptionPresent( "em", em );
	ops >> OptionPresent( "evaluatePWMs", evaluatePWMs );

	/* mode-dependent default setting */

	order = ( negFile == NULL ) ? 2 : 8; // background model order
	mergeMode = HIGH;
	maxMatchPositions = 17;
	minMatchPositions = 5;

	/*
	 * process flags
	 * * user
	 */

	ops >> OptionPresent( "batch", batch_mode );
	if( batch_mode ){
		printf( "Running in batch mode: no progress bars.\n" );
	}

	if( ops >> OptionPresent( "no-pwm-length-optimization" ) ){
		maximizeMotifLength = false;
	}

	ops >> OptionPresent( "mops", multipleOccurrence );
	ops >> OptionPresent( "oops", oneOccurrence );
	ops >> OptionPresent( "zoops", zeroOrOneOccurrence );
	if( ( oneOccurrence && multipleOccurrence) ||
		  ( oneOccurrence && zeroOrOneOccurrence ) ||
		  ( zeroOrOneOccurrence && multipleOccurrence ) ){
		fprintf( stderr, "One-occurrence-per-sequence model (--oops) and/or "
				"multiple-occurrence-per-sequence model (--mops) and/or "
				"zero-or-one-occurrence-per-sequence model (--zoops) cannot be "
				"used simultaneously.\n" );
		exit(-1);
	}
	if( !( oneOccurrence ) && !( zeroOrOneOccurrence )
		&& !( multipleOccurrence ) ){
		zeroOrOneOccurrence = true;
	}
	if( em && !( zeroOrOneOccurrence ) ){
		fprintf( stderr, "--em requires the zero-or-one-occurrence-per-sequence "
				 "(--zoops) mode.\n" );
		exit(-1);
	}

	ops >> OptionPresent( "revcomp", revcomp );

	ops >> OptionPresent( "localization", usePositionalProbs );
	if( usePositionalProbs &&
		( Global::posSet->max_leng != Global::posSet->min_leng ) ){
		fprintf( stderr, "Localization information is only applicable for input"
				" sequences having all the same length.\n");
		exit( -1 );
	}
	if( usePositionalProbs ){
		ops >> OptionPresent( "localization-ranking", positionalProbsRanking );
	}

	/*
	 * process flags: user
	 * --em-specific
	 */

	if( em ){
		ops >> OptionPresent( "msq", msq );
		ops >> OptionPresent( "interpolate", interpolate );
		ops >> OptionPresent( "noExpectationMaximizationPhase",
							   noExpectationMaximizationPhase );
		ops >> OptionPresent( "initInts", initInts );
		ops >> OptionPresent( "saveInitModels", saveInitModels );
		ops >> OptionPresent( "saveModels", saveModels );
		ops >> OptionPresent( "verbose", verbose );
		ops >> OptionPresent( "learnAlpha", learnHyperParameter );
		ops >> OptionPresent( "posAlpha", positionSpecificAlphas );
		ops >> OptionPresent( "debugAlpha", debugAlphalearning );


	}

	if( em || evaluatePWMs ){
		ops >> OptionPresent( "testPosSet", testPosSequences );
		ops >> OptionPresent( "testNegSet", testNegSequences );
		ops >> OptionPresent( "logProbs", logProbs );
	}

	/*
	 * process flags: developer
	 */

	ops >> OptionPresent( "ali-free", useAliFree );
	ops >> OptionPresent( "extra-homology-filter", removeHomology );
	ops >> OptionPresent( "filtering", repeatFiltering );
	ops >> OptionPresent( "lcf", lowComplexityFilter );
	ops >> OptionPresent( "noRefinementPhase", noRefinementPhase );
	ops >> OptionPresent( "ranks", useRankPvalues );

	/*
	 * process flags: developer
	 * --em-specific
	 */

	if( em ){
		ops >> OptionPresent( "saveExpectationMaximizationLikelihoods",
							   saveExpectationMaximizationLikelihoods );
		ops >> OptionPresent( "saveExpectationMaximizationModels",
							   saveExpectationMaximizationModels );
	}


	/*
	 * process option-preceding arguments
	 * * user
	 */

	ops >> Option( "negSet", negFile );
	ops >> Option( "background-model-order", order );

	if( ops >> Option( "pseudo", pseudo ) ){
		pseudo /= 100;
	}

	ops >> Option( 'g', "gaps", GAPS );

	ops >> Option( "type", type );
	if( type == NO_VALID_MOTIF_TYPE ){
		return false;
	}

	ops >> Option( "merge-motif-threshold", mergeMode );
	if( mergeMode == NO_VALID_MERGE_MODE ){
		return false;
	}

	if( ops >> Option( "max-match-positions", maxMatchPositions ) ){
		if( maxMatchPositions > 26 ){
			fprintf( stderr, "maxMatchPositions > 26 not possible.\n");
			exit( -1 );
		}
	}

	ops >> Option( "maxPosSetSize", maxPosSetSize );

	std::string tr;
	if( ops >> Option( "trackedMotif", tr ) ){
	  trackedMotifs.insert( tr );
	}

	ops >> Option( "format", seqFormat );
	if( seqFormat == NO_VALID_SEQ_FORMAT ){
		return false;
	}

	ops >> Option( "maxMultipleSequences", maxMultipleSequences );

	ops >> Option( "downstream", downstream );

	ops >> Option( 'm', "startMotif", startMotif );
	ops >> Option( 'p', "profileFile", profFile );

 	int offset = posSet->max_leng - downstream;
	if( ops >> Option( "startRegion", startRegion ) ){
		startRegion += offset;
	} else{
		startRegion = 0;
	}
	if( ops >> Option( "endRegion", endRegion ) ){
		Global::endRegion += offset;
	} else{
		endRegion = posSet->max_leng;
	}

	/*
	 * process option-preceding arguments
	 * * user
	 *   * --em-specific
	 */

	if( em ){

		ops >> Option( 'k', modelOrder );

		if( ops >> Option( "bindingSiteFile", bindingSiteFile ) ){
			if( ops >> OptionPresent( "bindingSiteLength" ) ){
				ops >> Option( "bindingSiteLength", bindingSiteLength );
			} else{

				/* determine binding site lengths */

				FILE* fp;
				if( ( fp = fopen( bindingSiteFile, "r" ) ) == NULL ){
			        fprintf( stderr, "Cannot open bindingSiteFile %s\n",
			        		 bindingSiteFile );
			        exit(-1);
				}

				int c;
				while( ( c=fgetc( fp ) ) == '\n' ){
					; // skip first blank lines
				}

				int l = 1; // character counter per line
				int L = -1; // character number of 1st line

				int ncounter = 0; // \n-counter
				while( ( c=fgetc( fp ) ) != EOF ){
					if( c == '\n' ){
						if( ncounter > 0 ){
							// skip blank lines
							continue;
						} else if( L == -1 ){
							// remember the 1st line's character number
							L = l;
						} else if ( L != l ){
							fprintf( stderr, "Binding site sequence lengths "
									 "must not differ\n" );
							exit( -1 );
						}
						l = 0; // reset character counter
						++ncounter;
					} else{ // next character
						++l;
						ncounter = 0; // reset \n-counter
					}
				}

				fclose( fp );

				bindingSiteLength = L;
			}

			if( ops >> OptionPresent( "addColumns" ) ){
				addColumns.clear();
				ops >> Option( "addColumns", addColumns );
				if( addColumns.size() < 1 || addColumns.size() > 2 ){
					fprintf( stderr, "--addColumns format error\n" );
					exit( -1 );
				}
				if( addColumns.size() == 1 ){
					addColumns.resize( 2, addColumns.back() );
				}
				bindingSiteLength = addColumns.at( 0 ) + bindingSiteLength +
						            addColumns.at( 1 );
			}

			if( bindingSiteLength > posSet->min_leng ){
				fprintf( stderr, "Binding site sequence lengths exeed positive "
						"sequence lengths\n" );
				exit( -1 );
			}

			if( !( modelOrder < bindingSiteLength ) ){
				modelOrder = bindingSiteLength - 1;
			}

			if( ops >> Option( "bindingSiteIntsFile", bindingSiteIntsFile ) ){

				std::ifstream fs( bindingSiteIntsFile, std::ios_base::in );
				if( fs.fail() ){
					fprintf( stderr, "Cannot open bindingSiteIntsFile %s\n",
							 Global::bindingSiteIntsFile );
					exit( -1 );
				}

				float weight;
				std::vector<float> weights;

				while( fs >> weight ){
					if( weight < 0.0f ){
						fprintf( stderr, "Use non-negative binding site "
								 "weights\n" );
						exit( -1 );
					}
					weights.push_back( weight );
				}
				fs.close();

				ops >> OptionPresent( "bindingSiteRankWeighting",
						              bindingSiteRankWeighting );

				if( bindingSiteRankWeighting ){

					float N = static_cast<float>( weights.size() );

					ops >> Option( "bindingSiteBackgroundQuantile",
							       bindingSiteBackgroundQuantile );

					if( !( ops >> Option( "bindingSiteBackgroundRank",
						bindingSiteBackgroundRank ) ) ){

						if( bindingSiteBackgroundQuantile > 1 - ( ( 1.0f / 2.0f
							) / static_cast<float>( N ) ) ){
							bindingSiteBackgroundRank = 1.0f;
						} else{
							bindingSiteBackgroundRank =
							N - roundf( N * bindingSiteBackgroundQuantile );
						}
					}
				} else{

					ops >> Option( "bindingSiteBackgroundQuantile",
							       bindingSiteBackgroundQuantile );

					if( !( ops >> Option( "bindingSiteBackgroundIntensity",
						bindingSiteBackgroundIntensity ) ) ){

						bindingSiteBackgroundIntensity =
						quantile( weights, bindingSiteBackgroundQuantile, 7 );
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
		        fprintf( stderr, "Cannot open markovModelFile %s\n",
		        		 str.str().c_str() );
		        exit(-1);
			}

			int c;
			while( ( c=fgetc( fp ) ) == '\n' ){
				; // skip first blank lines
			}

			int lines = 0; // lines without blank lines
			int positions = 0;

			int ncounter = 0; // \n-counter
			while( ( c=fgetc( fp ) ) != EOF ){
				if( c == '\n' ){
					++lines;
					++ncounter;
				} else{
					if( ncounter > 1 ){
						++positions;
						lines -= ( ncounter-1 ); // substract blank lines
					}
					ncounter = 0; // reset \n-counter
				}
			}
			if( ncounter > 1 ){ // ignore last blank lines
				lines -= ( ncounter-1 );
			}

			fclose( fp );

			markovModelLength = positions + 1;
			modelOrder = ( lines / markovModelLength ) - 1;

			if( !( modelOrder < markovModelLength ) ){
				fprintf( stderr, "Markov model format error\n" );
				exit( -1 );
			}

			if( ops >> OptionPresent( "addColumns" ) ){
				addColumns.clear();
				ops >> Option( "addColumns", addColumns );
				if( addColumns.size() < 1 || addColumns.size() > 2 ){
					fprintf( stderr, "--addColumns format error\n" );
					exit( -1 );
				}
				if( addColumns.size() == 1 ){
					addColumns.resize( 2, addColumns.back() );
				}
				markovModelLength = addColumns.at( 0 ) + markovModelLength +
						            addColumns.at( 1 );
			}

			if( markovModelLength > posSet->min_leng ){
				fprintf( stderr, "Markov model columns exeed positive sequence "
						"lengths\n" );
				exit( -1 );
			}

		} else{

			if( !( modelOrder < PWM_LENGTH ) ){
				modelOrder = PWM_LENGTH - 1;
			}

			// number of one or more models in the ranking to pursue
			if( ops >> OptionPresent( "nrModels" ) ){
				ops >> Option( "nrModels", nrModels );
				std::sort( nrModels.begin(), nrModels.end() );
			}
			// min. number of models to pursue
			ops >> Option( "minModels", minModels );
			// max. number of models to pursue
			ops >> Option( "maxModels", maxModels );
			// max. p-value of models
			ops >> Option( "maxPvalue", pValueThreshold );
			pValueThreshold = log( pValueThreshold );
			// min. percentage of sequences containing a binding site instance
			ops >> Option( "minOccurrence", minOccurrence );
			// add columns (left, right) to models
			if( ops >> OptionPresent( "addColumns" ) ){
				addColumns.clear();
				ops >> Option( "addColumns", addColumns );
				if( addColumns.size() < 1 || addColumns.size() > 2 ){
					fprintf( stderr, "--addColumns format error\n" );
					exit( -1 );
				}
				if( addColumns.size() == 1 ){
					addColumns.resize( 2, addColumns.back() );
				}
			}
		}

		if( ops >> Option( "sequenceIntsFile", sequenceIntsFile ) ){
			std::ifstream fs( sequenceIntsFile, std::ios_base::in );
			if( fs.fail() ){
				fprintf( stderr, "Cannot open sequenceIntsFile %s\n",
						 Global::sequenceIntsFile );
				exit( -1 );
			}

			float weight;
			std::vector<float> weights;

			while( fs >> weight ){
				if( weight < 0.0f ){
					fprintf( stderr, "Use non-negative sequence weights\n" );
					exit( -1 );
				}
				weights.push_back( weight );
			}
			fs.close();

			if( Global::posSet->nent != static_cast<int>( weights.size() ) ){
				fprintf( stderr, "The number of sequences and sequence weights "
						 "differs\n" );
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

			ops >> OptionPresent( "rankWeighting", rankWeighting );

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

		/* check specificity and pseudo-counts factor(s) */
		unsigned checkSum = 0;
		if( ops >> OptionPresent( 'q' ) ){
			ops >> Option( 'q', q );
			if( q <= 0 || q >= 1 ){
				fprintf( stderr,
						 "Specificity factor q restricted to 0 < q < 1.\n" );
				exit( -1 );
			}
			checkSum = checkSum + 4;
		}
		if( ops >> OptionPresent( 'a', "alpha" ) ){
			alpha.clear();
			checkSum = checkSum + 2;
		}
		switch( checkSum ){
		case 0:
			break;
		case 1:
			break;
		case 2:
			ops >> Option( 'a', "alpha", alpha );
			break;
		case 4:
			break;
		case 5:
			break;
		case 6:
			ops >> Option( 'a', "alpha", alpha );
			break;
		case 3: case 7:
			fprintf( stderr, "Specify either eta (--eta) or alpha (--alpha).\n"
					);
			exit( -1 );
		default:
			fprintf( stderr, "Specificity and pseudo-counts factor fatal error."
					"\n" );
			exit( -1 );
		}

		/* interpolated Markov model order k determines pseudo-counts factor
		 * vector size to be k+1
		 *
		 * remember: eta.size() equals alpha.size() */
		if( static_cast<int>( alpha.size() ) != modelOrder+1 ){
			if( static_cast<int>( alpha.size() ) > modelOrder+1 ){
				alpha.resize( modelOrder+1 );
			} else{
				alpha.resize( modelOrder+1, alpha.back() );
			}
		}

		ops >> Option( 'K', modelOrderBg );
		ops >> Option( 'A', "Alpha", alphaBg );

	} else if( evaluatePWMs ){

		// number of one or more models in the ranking to pursue
		if( ops >> OptionPresent( "nrModels" ) ){
			ops >> Option( "nrModels", nrModels );
			std::sort( nrModels.begin(), nrModels.end() );
		}
		// min. number of models to pursue
		ops >> Option( "minModels", minModels );
		// max. number of models to pursue
		ops >> Option( "maxModels", maxModels );
		// max. p-value of models
		ops >> Option( "maxPvalue", pValueThreshold );
		pValueThreshold = log( pValueThreshold );
		// min. percentage of sequences containing a binding site instance
		ops >> Option( "minOccurrence", minOccurrence );
	}

	/*
	 * process option-preceding arguments
	 * * developer
	 */

	ops >> Option( "benchmarkFolder", benchmarkFolder );
	ops >> Option( "counts-offset", countsOffset );

	ops >> Option( "final-filter", finalFilterThreshold );
	finalFilterThreshold = log( finalFilterThreshold );

	ops >> Option( "gapExtension", gapExtension );
	ops >> Option( "gapOpening", gapOpening );

	std::string thresh_init;
	if( ops >> Option( "instance-threshold", thresh_init ) ) {
		instanceThreshold = ThresholdChecker( thresh_init );
	} else{
		instanceThreshold = ThresholdChecker( "<1.0" );
	}

	ops >> Option( "maxIterations", maxIterations );
	ops >> Option( "maxMotifLevel", maxMotifLevel );
	ops >> Option( "maxSeqOcc", maxMotifsPerSequence );
	ops >> Option( "min-match-positions", minMatchPositions );

	/* neffs can be set identically (--neff) or seperately for discrete
	 * (--neff-states) and PWM (--neff-pwm) phase */
	int n_eff = -1;
	neff_discrete = 6;
	neff_pwm = 10;
	if( ops >> Option( "neff", n_eff ) ){
		neff_pwm = n_eff;
		neff_discrete = n_eff;
	}
	ops >> Option( "neff-states", neff_discrete );
	ops >> Option( "neff-pwm", neff_pwm );

	if( ops >> Option( "plusFrac", plusFrac ) ){
		plusFrac /= 100;
	}

	ops >> Option( "pseudocounts-factor", pseudocountsFactor );
	ops >> Option( "termFreqScale", dollarScaleFactor );

	/*
	 * process option-preceding arguments
	 * * developer
	 *   * --em-specific
	 */

	ops >> Option( "maxEMIterations", maxEMIterations );


	ops >> Option( "cswlen", cswlen, 13 );
	ops >> Option( "csbest", csbest, 0 );
	ops >> Option( "csprofiles", csprofiles );
	if( csbest != 0 ){
		if( csprofiles.length() == 0 ){
			fprintf( stderr, "Missing cs profile output file (--csprofiles).\n"
					);
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
		fprintf( stderr, "Minimum coverage > number of sequences. Does this "
				"makes sense?\n" );
		exit( -1 );
	}

	/*
	 * process option-preceding arguments
	 * * developer
	 *   * --em-specific
	 */

	if( em ){
		ops >> Option( "epsilon", epsilon );
	}

	if( em || evaluatePWMs ){
		if( ops >> Option( "testSet", testSequenceFile ) ){
			for( unsigned int i=0; i < testSequenceFile.size(); i++ ){
				testSet.push_back( readSeqSet( const_cast<char*>(
						           testSequenceFile.at(i).c_str() ), Global::A,
						           FASTA, INT_MAX ) );
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
		for( unsigned int i=0; i< args.size(); i++ ){
			printf( "%s\n", args.at(i).c_str() );
		}
		return false;
	}

	/* create output directory */
	outputDirectory = String( args[0].c_str() );
	createDirectory( outputDirectory );
	pwmFolder = String( outputDirectory );

	/* create temporary directory */
	char* tmpFolder = ( char* )calloc( 1024, sizeof( char ) );
	sprintf( tmpFolder, "%s/tmp", outputDirectory );
	createDirectory( tmpFolder );
	free( tmpFolder );

	/* extract file name from sequence file */
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

	/*
	 * process settings depending on arguments
	 * * --aa-specific
	 */

	std::stringstream s;
	s << Global::outputDirectory << "/" << Global::shortFileName << ".mtf";

	/*
	 * process settings depending on arguments
	 * * --em-specific
	 */

	if( em && verbose ){

		std::cout << std::endl;
		std::cout << " ____________________" << std::endl;
		std::cout << "|                    |" << std::endl;
		std::cout << "| EM MODE PARAMETERS |" << std::endl;
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
   	if(name != NULL)  free(name);
   	if(outputDirectory != NULL) free(outputDirectory);
   	if(shortFileName != NULL) free(shortFileName);
   	if(negFile != NULL) free(negFile);
   	if(startMotif != NULL) free(startMotif);
   	if(profFile != NULL) free(profFile);
   	if(benchmarkFolder != NULL) free(benchmarkFolder);
   	if(pwmFolder != NULL) free(pwmFolder);

	ss_type set = posSet;
	if(negSet != NULL) set = negSet;

	double POW_2_16 = pow(2,16);
	if(conservationProbs != NULL){
		for(int i=0; i<= POW_2_16; i++){
			if(conservationProbs[i] == NULL) continue;
			for(int j=1; j< set->max_MultSeq; j++){
				free(conservationProbs[i][j]);
			}
			free(conservationProbs[i]);
		}
		free(conservationProbs);
	}

	if(alignmentFreeProbs != NULL){
		for(int i=0; i<= POW_2_16; i++){
			if(alignmentFreeProbs[i] == NULL) continue;
			for(int j=1; j< set->max_MultSeq; j++){
				free(alignmentFreeProbs[i][j]);
			}
			free(alignmentFreeProbs[i]);
		}
		free(alignmentFreeProbs);
	}

	free(posBg_log);
   	free(negBg_log);
   	free(posBg);
   	free(negBg);

	NullModel::destruct();


   	if(negSet != NULL) NilSeqSet(negSet);

   	NilSeqSet(posSet);
   	NilAlpha(A);

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


std::ostream& operator<<(std::ostream &os, const merge_type &m) {
  switch(m) {
    case LOW: os << "LOW"; break;
    case MEDIUM: os << "MEDIUM"; break;
    case HIGH: os << "HIGH"; break;
    case NO_VALID_MERGE_MODE: os << "NO_VALID_MERGE_MODE"; break;
    default: os << "UNKNOWN - THIS MUST NOT HAPPEN!"; exit(1);
  }
  return os;
}


std::ostream& operator<<(std::ostream &os, const SuppInfMode &v) {
	switch (v) {
	case SUPP_NO: os << "NO"; break;
	case SUPP_DISOCONS: os << "DISOCONS"; break;
	case SUPP_NNET: os << "NNET"; break;
	default: os << "UNKNOWN - THIS MUST NOT HAPPEN!"; exit(1);
	}
	return os;
}

