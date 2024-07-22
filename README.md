Artefacts for Leaky McEliece Paper
==================================

Organization
------------

 - src/ - attack code
 - support-splitting/ - support splitting algorithm for the case of known Goppa points
 - sage/ - code to generate leaks for Classic McEliece toy parameter set
 - docker/ - code to generate all other leaks
 - results/ - scaffolding and output for the attack experiments
 - 3rdparty - third party support packages


Entrypoints
-----------

sage/synthetic-leak.sage
 - generate leak files for Classic McEliece toy parameter sets

(docker) leak-mceliece{348864,...}
 - generate leak files for Classic McEliece parameter set

(docker) leak-botan
 - generate leak files for Botan (any parameter set)

build/Stockfish
 - run the attack on a leak file (any implementation, any parameter set)

results/generate-artefacts.sh
results/generate-all.sh
 - generate a series of leak files and run the attack for a series of errors
 - note that the experiments are large, and you might need to copy the script
   and reduce the number of trials or parameter range

results/timing-measurements.sh
 - run a series of timing measurements

support-splitting/experiments.sage
 - run the support splitting algorithm on test matrices
 - verification of the results is slow, so you may want to disable it in experiments.sage


Build and Usage Instructions
----------------------------

    $ (cd sage/mceliece-sage-20221023 && make sagelibs)
    $ (cd docker && docker build -t mceliece .)
    $ (mkdir build && cd build && cmake ../src && make)

Example usage:

    $ ./Stockfish --debug --threads=8 recover classic-mceliece348864.leak

There is a --help option at the top level and for each subcommand. Currently there are only two subcommands:

Top-level options:
 - The number of parallel threads can be set with "--threads <NR>", and is initialized to the number of CPU threads.
 - Some additional output can be requested with "--debug".
 
./Stockfish setup <KEM>
 - This command pre-generates the codebook for the column candidates and saves it to a file.
 - Optional, as the `recover' subcommand will generate the codebook file on demand.
 - <KEM> must be provided and is in the form "botan51220" or "mceliece51220" for Botan and Classic McEliece, resp.
 - The codebook is saved as "cached-<KEM>.bin" in the current directory.

./Stockfish recover <LEAKFILE>
 - Run the attack on <LEAKFILE>, which must be provided.
 - The error probability can be set with "--error <FLOAT>", where the default for <FLOAT> is 0.3.
 - The leak error seed can be set with "--seed <UINT>", where 0 uses OS random.
   With a seed, the error generation and thus the attack is fully deterministic.
 - A report file (in JSON) can be requested by "--with-report".
   The report is saved to "report.json" in the current directory, or an alternative path by "--report file <FILENAME>".
 - For Classic McEliece, a histogram of the distances in the codebook can be requested by "--with-distances".
   The report is saved to "distances-report.csv" in the current directory, or an alternative path by "--distances-file <FILENAME>"
   Note that generating the full distances histogram for an attack is very CPU intensive and can take weeks of computation,
   so using many threads is recommended.

Generating leaks
----------------

### Format

The leak files are human readable text files that contain the following lines:

kem = 'mceliece348864'
 - The implementation and parameter set.
 - 'mceliece...' for Classic McEliece and 'botan...' for Botan.
 
seed = '0001020304050607...'
 - The random seed used for key generation. Length depends on implementation:
 - Classic McEliece: 48 Bytes (96 hex digits)
 - Botan: 32 Bytes (64 hex digits)

pk = '...'
 - The public key. Format depends on implementation.
 - Classic McEliece: Format defined by NIST submission, see src/classic.h::parse_pk().
 - Botan: For simplicity, we do not use the Botan ASN.1-based output format, but a simplified form:
   Raw public key matrix as hex string, see src/botan.h::parse_pk().
   See src/botan.h::parse:pk() for details.

sk = '...'
 - The secret key. Only used for validation. Format depends on implementation.
 - Classic McEliece: Format defined by NIST submission, see src/classic.h::ClassicSecret().
 - Botan: For simplicity, we do not use the Botan ASN.1-based output format, but a simplified form:
   Raw coefficients of the Goppa polynomial, followed by raw inverse permutation of L (before mapping to gray code).
   See src/botan.h::BotanSecret() for details.

leak = '00010001...' # The leak bits, as a sequence of 00 or 01 bytes.
 - Length is at least (mt)^2 - mt bytes.
 - Note that the leak my contain data from several key generation attempts (in case FAIL is returned).
   In this case, our attack only uses the leak data from the last, successful attempt.

mat = '...' # The matrix H used during public key generation.
 - Length is (mt)*(n/8) bytes.
 - This is used for differential power analysis. Due to restrictions on the system under test, we can not
   run the whole key generation on the device, so we run all steps up to gaussian elimination offline.
   The value of mat is then transfered to the device, where gaussian elimination is run under testing conditions.

ss = '...'
 - A sample secret message, Classic McEliece only.
 - Unused. It can be used to verify decryption, for example using the Sage implementation of Classic McEliece.

ct = '...'
 - A sample ciphertext, Classic McEliece only. Encryption of ss.
 - Unused. It can be used to verify decryption, for example using the Sage implementation of Classic McEliece.


### Classic McEliece, Toy Parameters

For Classic McEliece, toy parameters: The reference implementation does not
implement toy parameters, so there is a separate script sage/synthetic-leak.sage
based on the Classic McEliece sage implementation:

    $ sage sage/synthetic-leak.sage mceliece51220 sample-classic51220.leak
    $ sage sage/synthetic-leak.sage mceliece102450 sample-classic102450.leak

Note that key generation is deterministic. To create fresh keys, provide a different
96 character (=48 bytes) hex string on the command line, e.g.:

    $ docker run --rm -ti mceliece /leak-mceliece348864 C0FFEE118BADF00DC0FFEE118BADF00DC0FFEE118BADF00DC0FFEE118BADF00DC0FFEE118BADF00DC0FFEE118BADF00D > sample-classic348864-2.leak


### Classic McEliece, Regular Parameters

For Classic McEliece, regular parameters: We build the reference implementation
in a docker container, patched with a leak function that outputs the execution matrix
along with the secret and public key:

    $ docker run --rm -ti mceliece /leak-mceliece348864 > sample-classic348864.leak
    $ docker run --rm -ti mceliece /leak-mceliece460896 > sample-classic460896.leak
    $ docker run --rm -ti mceliece /leak-mceliece6960119 > sample-classic6960119.leak
    $ docker run --rm -ti mceliece /leak-mceliece6688128 > sample-classic6688128.leak
    $ docker run --rm -ti mceliece /leak-mceliece8192128 > sample-classic8192128.leak

Note that key generation is deterministic. To create fresh keys, provide a different
96 character (=48 bytes) hex string on the command line, e.g.:

    $ docker run --rm -ti mceliece /leak-mceliece348864 C0FFEE118BADF00DC0FFEE118BADF00DC0FFEE118BADF00DC0FFEE118BADF00DC0FFEE118BADF00DC0FFEE118BADF00D > sample-classic348864-2.leak


### Botan, Toy and Regular Parameters

For Botan, we build the botan CLI tool, patched with a leak function that outputs the
execution matrix along with the secret and public key. It works for all parameters,
toy and regular:

    $ docker run --rm -ti mceliece /leak-botan leak --n=512 --t=20 > sample-botan51220.leak
    $ docker run --rm -ti mceliece /leak-botan leak --n=3488 --t=64 > sample-botan348864.leak
    $ docker run --rm -ti mceliece /leak-botan leak --n=3488 --t=64 > sample-botan348864.leak
    $ docker run --rm -ti mceliece /leak-botan leak --n=3488 --t=64 > sample-botan348864.leak
    $ docker run --rm -ti mceliece /leak-botan leak --n=3488 --t=64 > sample-botan348864.leak

Note that key generation is deterministic. To create fresh keys, provide a different
64 character (=32 bytes) hex string on the command line, e.g.:

    $ docker run --rm -ti mceliece /leak-botan leak --n=512 --t=20 --seed=C0FFEE118BADF00DC0FFEE118BADF00DC0FFEE118BADF00DC0FFEE118BADF00D > sample-botan51220-2.leak

Experiments
-----------

### Threshold Probabilities

Threshold probabilities in Table 2 are proved by:
 - results/proof-thresholds.sage

When run, the file should output the threshold probabilities and not throw an assertion.

### Attack Experiments

For our experiments, we generated one leak file per experiment and ran the attack, outputting a result file.

The following script can be called with a KEM argument ("mceliece51220", "botan102450", etc.) to run a series of experiments
for one particular implementation and parameter set.
 - results/generate-artefacts.sage
   
The following parameters can be configured in the file:
 - START_TRIAL, END_TRIAL: The trials to run. Default is 0-99.
 - START_ERROR, END_ERROR, STEP_ERROR: The error probabilities to try (default: 0.00-0.50 or 0.30-0.50, depending on parameters).
 - STOCKFISH_THREADS: How many threads to allocate per experiment. Default is 1 for toy and 16 for regular parameters.
 - THREADS: How many experiments to run in parallel. Default is OS cpu count divided by STOCKFISH_THREADS.

Note that the total number of threads executed in parallel can be up to THREADS*STOCKFISH_THREADS.

The following script runs all experiments in sequence. Note that this can take a considerable amount of time (days to weeks),
even on very large hardware.
 - results/generate-all.sh

Afterwards, the report files are generated in subdirectories, one for each implementation/parameter set. The reports can be further
summarized with the following script:
 - results/generate-evaluation.sage

The output of the last script is a lot of CSV files. For each implementation/parameter set, there is:
 - results/error-rate-report-<KEM>.csv: threshold probability (linear interpolation), relative entropy, information leak
 - results/success-probability-estimated-<KEM>.csv: p(tau) according to our model
 - results/success-probability-report-<KEM>.csv: p(tau) according to our experiments

Because the leak and result.json files are very large (367 GB uncompressed), we do not include them in these artefacts. However,
for each implementation and parameter set, there is a pair of sample files for one leak and attack run:
 - results/samples/{*.leak,*.json}


### Distances Report

The option --with-distances also generates a .csv file with the inter-candidate distances in the candidate matrix.
This is very compute expensive. Using an AMD Epy 7763 with 256 threads this took:
 - distances report for 51220: real: 2m14 total: 505m37.650s
 - distances report for 102450: real: 145m4.394s total: 35013m14,973s 

The raw data that was output from running a successful attack on toy parameters can be found here:
 - results/distances-report-mceliece51220.csv
 - results/distances-report-mceliece102450.csv

Visualizations can be generated with:
 - results/distances-visualize.sage
 - results/distances-report-mceliece51220.png
 - results/distances-report-mceliece102450.png

Acknowledgements
----------------

The following 3rdparty libraries are used:

 - Parallelization: https://github.com/bshoshany/thread-pool
 - Fast Logging: https://github.com/gabime/spdlog
 - Command Line Parsing: https://github.com/CLIUtils/CLI11
 - JSON File Parsing: https://github.com/nlohmann/json 
