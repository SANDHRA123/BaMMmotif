# BaMM!motif

SYNOPSIS

      BaMMmotif DIRPATH FILEPATH [OPTIONS]

DESCRIPTION

      Bayesian Markov Model motif discovery software.

      DIRPATH
          Output directory for the results.

      FILEPATH
          FASTA file with positive sequences of equal length.

OPTIONS

  Sequence options

      --negSequenceSet <FILEPATH>
          FASTA file with negative/background sequences used to learn the
          (homogeneous) background BMM. If not specified, the background BMM is
          learned from the positive sequences.

      --reverseComp
          Search motifs on both strands (positive sequences and reverse
          complements). This option is e.g. recommended when using sequences
          derived from ChIP-seq experiments.

  Options to initialize a single BMM from file

      --bindingSiteFile <FILEPATH>
          File with binding sites of equal length (one per line).

      --markovModelFile <FILEPATH>
          File with BMM probabilities as obtained from BaMM!motif (omit
          filename extension).

  Options to initialize one or more BMMs from XXmotif PWMs

      --minPWMs <INTEGER>
          Minimum number of PWMs. The options --maxPValue and --minOccurrence
          are ignored. The default is 1.

      --maxPWMs <INTEGER>
          Maximum number of PWMs.

      --maxPValue <FLOAT>
          Maximum p-value of PWMs. This filter is not applied to the top
          minimum number of PWMs (see --minPWMs). The default is 1.0.

      --minOccurrence <FLOAT>
          Minimum fraction of sequences that contain the motif. This filter is
          not applied to the top minimum number of PWMs (see --minPWMs). The
          default is 0.05.

      --rankPWMs <INTEGER> [<INTEGER>...]
          PWM ranks in XXmotif results. The former options to initialize BMMs
          from PWMs are ignored.

  Options for (inhomogeneous) motif BMMs

      -k <INTEGER>
          Order. The default is 2.

      -a|--alpha <FLOAT> [<FLOAT>...]
          Order-specific prior strength. The default is 1.0 (for k = 0) and
          20 x 3^(k-1) (for k > 0). The options -b and -g are ignored.

      -b|--beta <FLOAT>
          Calculate order-specific alphas according to beta x gamma^(k-1) (for
          k > 0). The default is 20.0.

      -g|--gamma <FLOAT>
          Calculate order-specific alphas according to beta x gamma^(k-1) (for
          k > 0). The default is 3.0.

      --extend <INTEGER>{1,2}
          Extend BMMs by adding uniformly initialized positions to the left
          and/or right of initial BMMs. Invoking e.g. with --extend 0 2 adds
          two positions to the right of initial BMMs. Invoking with --extend 2
          adds two positions to both sides of initial BMMs. By default, BMMs
          are not being extended.

  Options for the (homogeneous) background BMM

      -K <INTEGER>
          Order. The default is 2.

      -A|--Alpha <FLOAT>
          Prior strength. The default is 10.0.

  EM options

      -q <FLOAT>
          Prior probability for a positive sequence to contain a motif. The
          default is 0.9.

      -e|--epsilon <FLOAT>
          The EM algorithm is deemed to be converged when the sum over the
          absolute differences in BMM probabilities from successive EM rounds
          is smaller than epsilon. The default is 0.001.

  XXmotif options

      --XX-ZOOPS
          Use the zero-or-one-occurrence-per-sequence model (default).

      --XX-MOPS
          Use the multiple-occurrence-per-sequence model.

      --XX-OOPS
          Use the one-occurrence-per-sequence model.

      --XX-seeds ALL|FIVEMERS|PALINDROME|TANDEM|NOPALINDROME|NOTANDEM
          Define the nature of seed patterns. The default is to start using ALL
          seed pattern variants.

      --XX-gaps 0|1|2|3
          Maximum number of gaps used for seed patterns. The default is 0.

      --XX-pseudoCounts <FLOAT>
          Percentage of pseudocounts. The default is 10.0.

      --XX-mergeMotifsThreshold LOW|MEDIUM|HIGH
          Define the similarity threshold used to merge PWMs. The default is to
          merge PWMs with LOW similarity in order to reduce runtime.

      --XX-maxPositions <INTEGER>
          Limit the number of motif positions to reduce runtime. The default is
          17.

      --XX-noLengthOptimPWMs
          Omit the length optimization of PWMs.

      --XX-K <INTEGER>
          Order of the (homogeneous) background BMM. The default is either 2
          (when learned on positive sequences) or 8 (when learned on background
          sequences).

      --XX-A <FLOAT>
          Prior strength of the (homogeneous) background BMM. The default is
          10.0.

      --XX-jumpStartPatternStage <STRING>
          Jump-start pattern stage using an IUPAC pattern string.

      --XX-jumpStartPWMStage <FILEPATH>
          Jump-start PWM stage reading in a PWM from file.

      --XX-localization
          Calculate p-values for positional clustering of motif occurrences in
          positive sequences of equal length. Improves the sensitivity to find
          weak but positioned motifs.

      --XX-localizationRanking
          Rank motifs according to localization statistics.

      --XX-downstreamPositions <INTEGER>
          Distance between the anchor position (e.g. the transcription start
          site) and the last positive sequence nucleotide. Corrects motif
          positions in result plots. The default is 0.

      --XX-batch
          Suppress progress bars.

  Options to score sequences

      --scorePosSequenceSet
          Score positive (training) sequences with optimized BMMs.

      --scoreNegSequenceSet
          Score background (training) sequences with optimized BMMs.

      --scoreTestSequenceSet <FILEPATH> [<FILEPATH>...]
          Score test sequences with optimized BMMs. Test sequences can be
          provided in a single or multiple FASTA files.

  Output options

      --saveInitBMMs
          Write initialized BMM(s) to disk.

      --saveBMMs
          Write optimized BMM(s) to disk.

      --verbose
          Verbose terminal printouts.

      -h, --help
          Printout this help.
