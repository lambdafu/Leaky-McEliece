import os
import subprocess
import sys

KEM=sys.argv[1]
if not KEM.startswith("mceliece") and not KEM.startswith("botan"):
    print("Usage: sage " + sys.argv[0] + " KEM")
    exit()

START_TRIAL=0
END_TRIAL=99
START_ERROR=0.0
if not KEM.endswith("51220") and not KEM.endswith("102450") and not KEM.endswith("348864"):
    START_ERROR=0.3
END_ERROR=0.5
STEP_ERROR=0.01
STOCKFISH_THREADS=1
THREADS=os.cpu_count()
if not KEM.endswith("51220") and not KEM.endswith("102450"):
    STOCKFISH_THREADS=16
    THREADS=int(os.cpu_count()/16)

SAGE="sage"
DOCKER="docker"
# Note that the synthetic leak already randomizes the seed.
LEAKGEN="../../sage/synthetic-leak.sage"
STOCKFISH="../../build/Stockfish"
DRY_RUN=False

os.environ["SAGE_NUM_THREADS"] = str(THREADS)

def run_cmd(args, **kwargs):
    print("Running: " + " ".join(args))
    if not DRY_RUN:
        return subprocess.run(args, **kwargs)

@parallel
def generate_one(kem, error, precision, trial):
    print("Running Job:", kem, error, trial)
    set_random_seed()
    leakfile = "{}-e{}-t{:03}.leak".format(kem, error, trial)
    reportfile = "{}-e{}-t{:03}.json".format(kem, error, trial)
    if kem == "mceliece51220" or kem == "mceliece102450":
        leakgen = run_cmd([SAGE, LEAKGEN, kem, leakfile])
    elif kem.startswith("mceliece"):
        byte_array = [randint(0, 255) for _ in range(48)]
        seed = ''.join('{:02x}'.format(byte) for byte in byte_array)
        with open(leakfile, "w") as fh:
            leakgen = run_cmd([DOCKER, "exec", "mceliece", "./leak-" + kem, seed], stdout=fh)
    else:
        byte_array = [randint(0, 255) for _ in range(32)]
        seed = ''.join('{:02x}'.format(byte) for byte in byte_array)
        with open(leakfile, "w") as fh:
            param_map = { "botan51220": (512, 20),
                          "botan102450": (1024, 50),
                          "botan348864": (3488, 64),
                          "botan460896": (4608, 96),
                          "botan6688128": (6688, 128),
                          "botan6960119": (6960, 119),
                          "botan8192128": (8192, 128) }
            pn, pt = param_map[kem]
            leakgen = run_cmd([DOCKER, "exec", "mceliece", "./leak-botan", "leak",
                               "--seed=" + seed, "--n=" + str(pn), "--t=" + str(pt)], stdout=fh)
    recover = run_cmd([STOCKFISH, "--threads=" + str(STOCKFISH_THREADS), "recover", "--error={}".format(error), "--with-report", "--report-file", reportfile, leakfile])
    return True

trials = range(START_TRIAL, END_TRIAL+1)
errors = [START_ERROR + i*STEP_ERROR for i in range(int((END_ERROR-START_ERROR)/STEP_ERROR) + 1)]
precision=6 #len(str(float(abs(STEP_ERROR)))) - 2
errors = [("{:01." + str(precision) + "f}").format(e) for e in errors]

jobs = []
for t in trials:
    for e in errors:
        jobs.append((KEM, e, 2, t))

# Ensure cache file exists.
run_cmd([STOCKFISH, "setup", KEM])

results = generate_one(jobs)
for f in results:
    pass

