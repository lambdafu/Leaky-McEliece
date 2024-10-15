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
import scipy.stats

params = SimpleNamespace(t=128, m=13, n=129)
TRACEDIR="traces"
RESULTDIR="results"
TRACES=1000
VICTIM_TRACES=10
SAMPLES=287002932

# Some derived values: Number of rows, number of PoIs we expect
PK_NROWS = params.t * params.m
POIS_DIAG = sum(range(PK_NROWS - 1, PK_NROWS - 1 - params.n, -1))
POIS_ZERO = (PK_NROWS - 1) * params.n
POIS = POIS_DIAG + POIS_ZERO


# This is the output of running 04-analysis.sage!
# This defines PERIOD_ZERO, PERIOD_SHIFT and model_leakpos.
load(os.path.join(RESULTDIR, "leakpos.sage"))

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("mceliece")
logger.setLevel(logging.DEBUG)

file_handler = logging.FileHandler(os.path.join(RESULTDIR, '05-verify-classifier.log'), mode='w')
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

    def analyse_groundtruth(self):

        # Let's solve the binary classification problem for the leakpos bits that are not missing.
        # Collect all sample values for 0 and 1 for each leakpos.
        leakpos = []
        for c in range(params.n):
            for r in range(PK_NROWS):
                if c == r:
                    continue
                leakpos.append(model_leakpos(r, c))

        # For development, set PARTIAL auf True to work on a part of the data for speed.
        PARTIAL=False
        if PARTIAL:
            logger.warning("development mode, not calculating on full set of data!")
            count = 100 
            leakpos = leakpos[:100]
        else:
            count = self.count

        leak_data = np.zeros((count, len(leakpos)), dtype=np.uint8)
        trace_data = np.zeros((count, len(leakpos)))
        logger.info("reading trace data")
        for t in tqdm(range(count)):
            key, trace = self.get_trace(t)
            leak_data[t] = np.frombuffer(key.leak[:len(leakpos)], dtype=np.uint8)
            trace_data[t] = trace[leakpos]

        nr_leakpos = len(leakpos)
        class_count = np.zeros((2, nr_leakpos))
        class_mean = np.zeros((2, nr_leakpos))
        class_std = np.zeros((2, nr_leakpos))

        snr = np.zeros(nr_leakpos)

        mask = [leak_data == 0, leak_data == 1]
        logger.info("generating statistics for %s", self.name)
        for idx in tqdm(range(len(leakpos))):
            class_count[1, idx] = np.sum(leak_data[:, idx])
            class_count[0, idx] = count - class_count[1, idx]

            for cls in [0, 1]:
                col = trace_data[:, idx][mask[cls][:, idx]]
                class_mean[cls, idx] = np.mean(col)
                class_std[cls, idx] = np.std(col)

            # SNR analysis
            signal = class_mean[0, idx]*(1-leak_data[:, idx]) + class_mean[1, idx]*leak_data[:, idx]
            noise = trace_data[:, idx] - signal
            snr[idx] = 10*np.log10(np.mean(signal**2) / np.mean(noise**2))

        for cls in [0, 1]:
            logger.info("statistics for class %i", cls)
            logger.info("class %i count: %.4f +- %.4f", cls, np.mean(class_count[cls]), np.std(class_count[cls]))
            logger.info("class %i mean: %.4f +- %.4f", cls, np.mean(class_mean[cls]), np.std(class_mean[cls]))
            logger.info("class %i std: %.4f +- %.4f", cls, np.mean(class_std[cls]), np.std(class_std[cls]))
        logger.info("SNR(dB): %.4f +- %.4f", np.mean(snr), np.std(snr))
       
        diff_of_means = abs(np.mean(class_mean[0]) - np.mean(class_mean[1]))
        sum_of_std = np.mean(class_std[0]) + np.mean(class_std[1])
        logger.info("diff of means is %.4f, sum of std is %.4f", diff_of_means, sum_of_std)
        if diff_of_means > sum_of_std:
            logger.info("diff of means is larger than sum of std, good!")
        else:
            logger.info("diff of means is smaller than sum of std, bad :(")
        zval = (diff_of_means / 2) / (sum_of_std / 2)
        cum_prob =  scipy.stats.norm.cdf(zval)
        logger.info("z-value: %.4f, estimated error rate: %.4f", zval, 1-cum_prob)

        # It's expensive to calcuate the class statistics at each position, so we store them in a file.
        np.save(os.path.join(RESULTDIR, "class-count.npy"), class_count)
        np.save(os.path.join(RESULTDIR, "class-mean.npy"), class_mean)
        np.save(os.path.join(RESULTDIR, "class-std.npy"), class_std)
        np.save(os.path.join(RESULTDIR, "snr.npy"), snr)
              
def main():
    tr = Traces("seedr")

    tr.analyse_groundtruth()

main()

