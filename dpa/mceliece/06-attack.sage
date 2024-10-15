#!/bin/env sage

import numpy as np
import random
import time
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
import struct
import logging
from types import SimpleNamespace
import sys
import os
import subprocess
import shutil
import binascii

params = SimpleNamespace(t=128, m=13, n=129)
TRACEDIR="traces"
RESULTDIR="results"
TRACES=10
SAMPLES=287002932

# Some derived values: Number of rows, number of PoIs we expect
PK_NROWS = params.t * params.m

# This is the output of running 04-analysis.sage!
# This defines PERIOD_ZERO, PERIOD_SHIFT and model_leakpos.
load(os.path.join(RESULTDIR, "leakpos.sage"))

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("mceliece")
logger.setLevel(logging.DEBUG)

file_handler = logging.FileHandler(os.path.join(RESULTDIR, '06-attack.log'), mode='w')
formatter = logging.Formatter('%(asctime)s:%(levelname)s :%(name)s:%(message)s')
file_handler.setFormatter(formatter)
file_handler.setLevel(logging.DEBUG)
logger.addHandler(file_handler)

# Simple reader only.
class McElieceKey(object):
    def __init__(self, fname=None):
        with open(fname, 'r') as fh:
            self.output = fh.read()
        self.loads(self.output)

    def loads(self, kvpairs):
        result = {}
        lines = kvpairs.split('\n')
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            key, value = map(str.strip, line.split('=', 1))
            value = value.strip("'\"")
            if key not in ["kem"]:
                value = bytes.fromhex(value)
            setattr(self, key, value)

        # self.leak now contains the full leak, including data from failed generations. Trim that.
        nr_leakbits = (PK_NROWS - 1) * PK_NROWS
        self.leak = self.leak[-nr_leakbits:]

    def save(self, fname):
        with open(fname, "w") as fh:
            keys = ['kem', 'seed', 'pk', 'sk', 'ct', 'ss', 'leak']
            for key in keys:
                value = getattr(self, key)
                if key not in ["kem"]:
                    value = binascii.hexlify(value).decode().upper()
                fh.write("{} = '{}'\n".format(key, value))

class Traces(object):
    def __init__(self, name):
        self.name = name
        self.fname_base = os.path.join(TRACEDIR, name)
        self.fname = self.fname_base + ".npy"
        self.fname_nfo = self.fname_base + ".nfo"
        self.fname_avg = self.fname_base + "-ravg.npy"
        # welford's online m2: https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
        self.fname_wm2 = self.fname_base + "-rwm2.npy"
        self.count = int(np.load(self.fname))
        if os.path.exists(self.fname_nfo):
            self.key = McElieceKey(fname=self.fname_nfo)
        else:
            self.key = None
        self.avg = np.load(self.fname_avg)
        self.wm2 = np.load(self.fname_wm2)

    def get_trace(self, idx):
        fname_base = "{}-t{:04}".format(self.fname_base, idx)
        key = McElieceKey(fname_base + ".nfo")
        trace = np.load(fname_base + ".npy")
        return key, trace

def attack(name, tr, key, trace):
    leakpos = []
    for c in range(params.n):
        for r in range(PK_NROWS):
            if c == r:
                continue
            leakpos.append(model_leakpos(r, c))
    average = tr.avg[leakpos]
    victim = trace[leakpos]

    # We save this for calculation of SNR and the true error probability.
    groundtruth = np.frombuffer(key.leak[:len(leakpos)], dtype=np.uint8)

    # Binary mean threshold classifier
    leak = np.where(victim > average, 0, 1)
    key.leak = leak.astype(np.uint8).tobytes()
    fname_nfo = os.path.join(RESULTDIR, name + ".nfo")
    key.save(fname_nfo)

    # Calculate the error rate.
    errors = 0
    for i in range(len(leakpos)):
        if groundtruth[i] != leak[i]:
            errors = errors + 1
    tau = errors / len(leakpos)
    logger.info("leak bit error rate = %.4f", errors / len(leakpos))
    # true error rate considers diagonal
    logger.info("true error rate tau = %.4f", (params.n/2 + errors) / (params.n + len(leakpos)))
    
    # Calculate the SNR.
    class_mean = np.load(os.path.join(RESULTDIR, "class-mean.npy"))

    signal = np.zeros(len(leakpos))
    signal[groundtruth == 0] = class_mean[0, groundtruth == 0]
    signal[groundtruth == 1] = class_mean[1, groundtruth == 1]
    noise = victim - signal
    snr = 10 * np.log10(np.mean(signal**2) / np.mean(noise**2))
    logger.info("SNR: %.4f", snr)
    return True, tau, snr, fname_nfo

def main():
    tr = Traces("seedr")
    tv = Traces("seedv")

    count = tv.count
    successes = np.zeros(count)
    taus = np.zeros(count)
    snrs = np.zeros(count)

    for i in range(count):
        logger.info("attacking key %s/%i", tv.name, i)
        key, trace = tv.get_trace(i)
        successes[i], taus[i], snrs[i], fname_nfo = attack("victim-{:02}".format(i), tr, key, trace)
        # Now let's run the key recovery.
        fname_json = os.path.join(RESULTDIR, "victim-{:02}.json".format(i))
        fname_sfl = os.path.join(RESULTDIR, "victim-{:02}.sfl".format(i))
        command = [ "../../Leaky-McEliece/build/Stockfish", "recover",
                "--external-error", "--error=%.4f" % (taus[i],),
                "--with-report", "--report-file=" + fname_json,
                fname_nfo ]
        with open(fname_sfl, "w") as fh:
            result = subprocess.run(command, stdout=fh, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            logger.warning("stockfish failed with return code: %i", result.returncode)

    logger.info("SUMMARY")
    logger.info("success rate: %.4f (%i out of %i)", sum(successes) / count, sum(successes), count)
    logger.info("tau: %.4f +- %.4f", np.mean(taus), np.std(taus))
    logger.info("SNR: %.4f +- %.4f", np.mean(snrs), np.std(snrs))

main()

