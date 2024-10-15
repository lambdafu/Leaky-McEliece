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

params = SimpleNamespace(t=128, m=13, n=129)
TRACEDIR="traces"
RESULTDIR="results"
TRACES=1000
VICTIM_TRACES=10
SAMPLES=287002932

# Set to true if you want to run slow, optional diagnosts
RUN_SLOW=True

# Some derived values: Number of rows, number of PoIs we expect
PK_NROWS = params.t * params.m
POIS_DIAG = sum(range(PK_NROWS - 1, PK_NROWS - 1 - params.n, -1))
POIS_ZERO = (PK_NROWS - 1) * params.n
POIS = POIS_DIAG + POIS_ZERO

# Special value to indicate missing pois in lists of indices.
MISSING = -1

# Some values we know from previous analysis. We can use them for automatic verification.
VERIFY = True
VERIFY_PERIOD_DIAG = 630
VERIFY_PERIOD_ZERO = 732
VERIFY_PERIOD_SHIFT = 33
def VERIFY_MODEL_LEAKPOS(r, c):
    if r == c:
        return -1
    PERIOD_SHIFT=33
    PERIOD_ZERO=732
    pos = -315.0*c^2 + 2264764.75*c + 1047823.0
    pos += -5.75*(c%8)
    pos += PERIOD_ZERO*r
    if r > c:
        pos += PERIOD_SHIFT - PERIOD_ZERO
    return int(pos)

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("mceliece")
logger.setLevel(logging.DEBUG)

file_handler = logging.FileHandler(os.path.join(RESULTDIR, '04-find-leak.log'), mode='w')
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

def report_traces(traces):
    for t in traces:
        logger.info("%s: %i traces", t.name, t.count)


class FixedVsRandom(object):
    def __init__(self, trace_fixed, trace_random):
        self.tf = trace_fixed
        self.tr = trace_random
        self.name = "FvR({},{})".format(self.tf.name, self.tr.name)
        self.logger = logging.getLogger(self.name)
        self.logger.setLevel(logging.DEBUG)
        self.logger.addHandler(file_handler)

        # Calculate the Welch t-statistic point-wise on our sample array.
        count = self.tf.count
        assert(count == self.tr.count)
        self.stat = np.abs((self.tf.avg - self.tr.avg)
                / sqrt((self.tf.wm2 + self.tr.wm2) / (count - 1) / count))

        # Diagnostics for the t-values.
        stat_max = np.max(self.stat)
        self.logger.info("max t-value: %.2f", stat_max)

        if RUN_SLOW:
            overcount = 1.1
            self.logger.info("Generating t-value approximation with overcount factor %.2f ...", overcount)
            perc = np.percentile(self.stat, 100*overcount*POIS/SAMPLES)
            stat_sorted = np.sort(self.stat)
            approx_tvalue = stat_sorted[-int(perc * SAMPLES)]
            self.approx_tvalue = approx_tvalue
            self.logger.info("PoIs percentile: %.4f, approx cutoff: %.2f",
                perc, approx_tvalue)
        else:
            self.approx_tvalue = 13

    # Find the periods of the for-loops for the diag and zero stage.
    # SLIP is the percentage of allowed difference to estimate.
    def find_loop_periods(self, slip=0.25):
        # The approximate period of the signal we are looking for.
        approx_period = ceil(SAMPLES / POIS)
        # We want to be within 25% of our approximate signal.
        period_min = (1 - slip) * approx_period
        period_max = (1 + slip) * approx_period
        self.logger.info("find loop periods between %i and %i samples", period_min, period_max)
        
        # We will do the FFT over a range of samples approximately covering the
        # first 8 columns of gaussian elimination.
        # We need plenty of samples here to avoid aliasing effects (side bands appearing instead of the actual signal).
        fftcount = approx_period * (PK_NROWS - 1) * 8

        # Just as a sanity check, we check that we have sufficient samples.
        # We want to be able to nail down the bucket for lowest frequency to half a sample.
        min_fft_samples = ceil(1/(1/period_max - 1/(period_max+0.5)))
        assert(fftcount > min_fft_samples)

        # Now do the actual FFT.
        sig = self.stat[:fftcount]
        fft = abs(np.fft.rfft(sig))
        fftfreq = np.fft.rfftfreq(int(fftcount))

        # This acts as a simple band filter, that removes a couple of problems:
        # 1. There are some high-amplitude signals at very low frequencies (including the DC offset)
        # 2. There are some high-amplitude signals at very high frequencies (for example with period 7)
        fftidx_lo = np.searchsorted(fftfreq, 1/period_max)
        fftidx_hi = np.searchsorted(fftfreq, 1/period_min)

        # Store a data table for the frequency graph. Note that we scale by dividing with length of signal.
        fft_graph = sorted(zip(1/(fftfreq[range(fftidx_lo, fftidx_hi)]),
            fft[fftidx_lo:fftidx_hi]))
        with open(os.path.join(RESULTDIR, self.tf.name + "-fft.csv"), "w") as fh:
            fh.write("x;a\n");
            for x, a in fft_graph:
                fh.write("%.4f;%.8f\n" % (x, a / fftcount))

        peakpos = np.argsort(fft)
        peakpos = peakpos[(peakpos >= fftidx_lo) & (peakpos <= fftidx_hi)]
        peakpos = np.flip(peakpos)
        # Values below threshold are excluded, as are ZONE samples before and after a peak.
        PEAK_THRESHOLD = fft.mean() + 9*fft.std()
        PEAK_ZONE = 10
        peaks = []
        self.logger.info("searching for peaks above %.2f with clearance %i either side", PEAK_THRESHOLD, PEAK_ZONE)
        def peak_period(peak):
            return floor(0.5+1/fftfreq[peak])
        while peakpos.size > 0:
            idx = peakpos[0]
            if fft[idx] < PEAK_THRESHOLD:
                break
            peaks.append(idx)
            remove_lo = np.searchsorted(fftfreq, 1/(peak_period(idx)+PEAK_ZONE))
            remove_hi = np.searchsorted(fftfreq, 1/(peak_period(idx)-PEAK_ZONE))
            peakpos = peakpos[(peakpos < remove_lo) | (peakpos > remove_hi)]
        for peak in peaks:
            self.logger.info("found signal with period %i and strength %.2f", peak_period(peak), fft[peak])

        result = sorted([peak_period(idx) for idx in peaks])
        if VERIFY:
            assert(len(result) == 2)
            PERIOD_DIAG, PERIOD_ZERO = result
            assert(PERIOD_DIAG == VERIFY_PERIOD_DIAG)
            assert(PERIOD_ZERO == VERIFY_PERIOD_ZERO) 

        return result

    def find_poi_candidates(self, threshold=13, clearance=5):
        # Values below threshold are excluded, as are ZONE samples before and after a peak.
        POI_THRESHOLD = threshold
        POI_ZONE = clearance

        poipos = np.where(self.stat > POI_THRESHOLD)[0]
        # Gymnastics to reorder poipos by the value in stat
        stat_values_at_poipos = self.stat[poipos]
        sorted_indices = np.argsort(stat_values_at_poipos)
        poipos = poipos[sorted_indices]
        poipos = np.flip(poipos)

        pois = []
        self.logger.info("searching for pois with t-values above %.2f and clearance %i", POI_THRESHOLD, 2*POI_ZONE+1)
        pois_skip = np.zeros(SAMPLES, dtype='bool')
        pbar = tqdm(range(poipos.size))
        for i in pbar:
            idx = poipos[i]
            if pois_skip[idx]:
                continue
            if self.stat[idx] < POI_THRESHOLD:
                break
            if len(pois) % 1000 == 0:
                pbar.set_description("strength {:.2f}".format(self.stat[idx]))
            pois.append(idx)
            r1 = max(0, idx - POI_ZONE)
            r2 = min(SAMPLES, idx + POI_ZONE)
            pois_skip[r1:r2+1] = True
        pois = np.array(pois)
        self.logger.info("found %i poi candidates", len(pois))
        return pois        

    def find_long_sequences(self, indices, PERIOD, ALLOWED_MISSING=1, MIN_LENGTH=7):
        # Ensure the indices are sorted
        indices = np.sort(indices)
        
        all_seqs = []
        longest_seq = []
        
        i = 0
        k = 0
        total_indices = len(indices)
        self.logger.info("searching for long sequences at least %i samples long", MIN_LENGTH)
        with tqdm(total=len(indices)) as pbar:
            while i < len(indices):
                current_seq = [indices[i]]
                next_value = indices[i] + PERIOD
            
                allowed_missing = ALLOWED_MISSING
                last_idx = 0
                for j in range(i + 1, len(indices)):
                    if indices[j] > next_value + 1:
                        allowed_missing -= 1
                        if allowed_missing < 0:
                            break
                        current_seq.append(MISSING)
                        next_value += PERIOD
                    if indices[j] == next_value - 1 or indices[j] == next_value or indices[j] == next_value+1:
                        current_seq.append(indices[j])
                        last_idx = j
                        allowed_missing = ALLOWED_MISSING
                        next_value += PERIOD
                    
                while current_seq[-1] == MISSING:
                    current_seq = current_seq[:-1]

                if len(current_seq) > len(longest_seq):
                    longest_seq = current_seq
                    if len(longest_seq) >= MIN_LENGTH:
                        i = last_idx
                        all_seqs.append(np.array(longest_seq))
                        longest_seq = []
                i = i + 1
                if k < 10:
                    k = k + 1
                pbar.update(int(i - pbar.n))
        assert(len(all_seqs) == 2 * params.n  - MIN_LENGTH)
        list_of_lens = sorted([len(x) for x in all_seqs])
        expected_list_of_lens = list(range(MIN_LENGTH, params.n)) + list(range(PK_NROWS - params.n, PK_NROWS)) 
        assert(list_of_lens == expected_list_of_lens)
        return all_seqs

    # Helper function used in cleanup_sequences
    def find_period_shift(self, PERIOD_ZERO, seqs_before_if, seqs_after_if):
        # We expect a small shift when the if condition is true.
        # Measure this shift here.
        before_and_after_deltas = []
        for col in range(params.n):
            seq1 = seqs_before_if[col]
            seq2 = seqs_after_if[col]
            seq1 = (seq1[seq1 != MISSING]) % PERIOD_ZERO
            seq2 = (seq2[seq2 != MISSING]) % PERIOD_ZERO
            if len(seq1) == 0 or len(seq2) == 0:
                continue
            before_and_after_deltas.append(int(round(seq2.mean() - seq1.mean()) % PERIOD_ZERO))
        all_deltas = list(set(before_and_after_deltas))
        if len(all_deltas) != 1:
            self.logger.warning("there is not a unique delta for if-branch. use majority vote instead?")
        PERIOD_SHIFT = all_deltas[0]
        self.logger.info("detected period shift for if-branch: %i", PERIOD_SHIFT)
        if VERIFY:
            assert(PERIOD_SHIFT == VERIFY_PERIOD_SHIFT)
        return PERIOD_SHIFT

    # seqs is now a complete set of leak sequences.
    def derive_model(self, PERIOD_ZERO, PERIOD_SHIFT, seqs_before_if, seqs_after_if):
        # Now let's model the location of the PoIs.
        # Note that this can not be a simple regression, because there is a second order effect
        # due to bit packing. 8 columns are packed in one byte, and accessing columns for the
        # check if the matrix is systematic takes time depending on col%8.
        # This is why we only take every 8th column here, and correct for that later.
        COL_TIMESTAMPS_STEP=8
        col_timestamps = []
        fit_col_timestamps = []
        for col in range(0, params.n):
            seq = seqs_before_if[col]
            if len(seq) == 0 or seq[0] == MISSING:
                continue
            col_timestamps.append([col, seq[0]])
            if col % COL_TIMESTAMPS_STEP == 0:
                fit_col_timestamps.append([col, seq[0]])

        col_timestamps = np.array(col_timestamps)
        fit_col_timestamps = np.array(fit_col_timestamps)

        # Intercept, slope, and acceleration
        var('ci', 'cs', 'ca')
        model(x) = ca*x**2 + cs*x + ci
        fit_result = find_fit(fit_col_timestamps, model)
        fit_values = {}
        for eq in fit_result:
            fit_values[eq.left()] = round(eq.rhs(), 4)
        # Result is: {ca: -315.0, ci: 1047823.0, cs: 2264764.75}
        first_approx = model(ca=fit_values[ca], cs=fit_values[cs], ci=fit_values[ci])

        # Now do a linear fit on col % COL_TIMESTAMP_STEP
        # Can't use 0-8 because our full before-seqs start at 5. So we shift from 8 down to 0.
        fit_diff = np.array([[col-COL_TIMESTAMPS_STEP, ts-first_approx(x=col)] for col,ts in col_timestamps])
        fit_diff=fit_diff[(fit_diff[:,0] >= 0) & (fit_diff[:,0] < COL_TIMESTAMPS_STEP)]
        fit_diff = fit_diff.astype(np.float64)
        var('mci', 'mcs')
        model_diff(x) = mcs*x + mci
        fit_result_diff = find_fit(fit_diff, model_diff)
        fit_values_diff = {}
        for eq in fit_result_diff:
            fit_values_diff[eq.left()] = round(eq.rhs(), 4)
        # Result is: {ca: -315.0, ci: 1047823.0, cs: 2264764.75}
        second_approx = model_diff(mci=fit_values_diff[mci], mcs=fit_values_diff[mcs])

        def model_start_ts_col(col):
            return Integer(first_approx(x=col) + second_approx(x=col%8))

        diff_to_model = [model_start_ts_col(c) for c in range(params.n - col_timestamps.shape[0], params.n)] - col_timestamps[:,1]
        if (diff_to_model != 0).any():
            self.logger.info("Model of start of zero stage columns differs for some values")

        # Find the PoI for zero stage processing column col and row.
        def model_leakpos(row, col):
            if row == col:
                return MISSING
            elif row < col:
                return model_start_ts_col(col) + row * PERIOD_ZERO
            else: # row > col
                return model_start_ts_col(col) + (row - 1) * PERIOD_ZERO + PERIOD_SHIFT

        model_file = ("# generated by 04-find-leak.sage\n"
                "MISSING=%i\n"
                "PERIOD_SHIFT=%i\n"
                "PERIOD_ZERO=%i\n"
                "def model_leakpos(r, c):\n"
                "    if r == c:\n"
                "        return MISSING\n"
                "    pos = %s\n"
                "    pos += %s\n"
                "    pos += PERIOD_ZERO*r\n"
                "    if r > c:\n"
                "        pos += PERIOD_SHIFT - PERIOD_ZERO\n"
                "    return int(pos)\n" % (
                MISSING, PERIOD_SHIFT, PERIOD_ZERO,
                repr(first_approx).replace('x', 'c'),
                repr(second_approx).replace('x', '(c%8)'))
        )
        self.logger.info("found model for leak positions:\n%s", model_file)

        if VERIFY:
            for r in range(PK_NROWS):
                for c in range(params.n):
                    assert(model_leakpos(r, c) == VERIFY_MODEL_LEAKPOS(r, c))

        with open(os.path.join(RESULTDIR, "leakpos.sage"), "w") as fh:
            fh.write(model_file)

        return model_leakpos

    def cleanup_sequences(self, seqs, PERIOD_ZERO):
        # We need two sequences for each column, a short one and a long one.
        MIN_LENGTH = 2 * params.n - len(seqs)

        # Now verify that we collected all the bits.
        idx = 0
        seqs_before_if = []
        seqs_after_if = []
        leakpos = []
        for col in range(params.n):
            seq = seqs[idx]
            if col < MIN_LENGTH:
                # For the first MIN_LENGTH columns, we do not expect to find the iterations before
                # the if-condition.
                seqs_before_if.append(np.array([MISSING]*col))
            else:
                # For later columns, we actually identify all iterations of the loop.
                assert(len(seq) == col)
                bits_missing = sum(seq == MISSING)
                if bits_missing > 0:
                    self.logger.info("Bits missing in column {} before if: {} {}".format(col, bits_missing, np.where(seq == MISSING)[0]))
                seqs_before_if.append(seq)
                idx = idx + 1
                seq = seqs[idx]
            assert(len(seq) == PK_NROWS - 1 - col)
            bits_missing = sum(seq == MISSING)
            if bits_missing > 0:
                self.logger.info("Bits missing in column {} after if: {} {}".format(col, bits_missing, np.where(seq == MISSING)[0]))
            seqs_after_if.append(seq)
            idx = idx + 1
            leakpos.extend(seqs_before_if[-1])
            leakpos.extend(seqs_after_if[-1])
        leakpos = np.array(leakpos)
        expected_leakpos_len = (PK_NROWS - 1) * params.n
        if (len(leakpos) != expected_leakpos_len):
            self.logger.info("warning: did not get expected number of leaked bits:", len(leakpos), expected_leakpos_len)
        assert(len(leakpos) == expected_leakpos_len)
        # We expect at least MIN_LENGTH*(MIN_LENGTH-1)/2 missing elements.
        nr_missing_expected = MIN_LENGTH*(MIN_LENGTH-1)/2
        nr_missing = sum(leakpos == MISSING)
        if nr_missing != nr_missing_expected:
            self.logger.info("warning: some bits were missing: %i", nr_missing - nr_missing_expected)

        PERIOD_SHIFT = self.find_period_shift(PERIOD_ZERO, seqs_before_if, seqs_after_if)
        model_leakpos = self.derive_model(PERIOD_ZERO, PERIOD_SHIFT, seqs_before_if, seqs_after_if)

        np.save(os.path.join(RESULTDIR, "leakpos-{}.npy".format(self.tf.name)), leakpos)

        return leakpos, model_leakpos

    def verify_model(self, leakpos, model_leakpos):
        predicted_leakpos = []
        for c in range(params.n):
            for r in range(PK_NROWS):
                if c == r:
                    continue
                predicted_leakpos.append(model_leakpos(r, c))
        predicted_leakpos = np.array(predicted_leakpos)
        non_missing_positions = np.where(leakpos != MISSING)[0]
        diff_to_prediction = leakpos[non_missing_positions] - predicted_leakpos[non_missing_positions]
        if (abs(diff_to_prediction) > 1).any():
            self.logger.warning("model differs a lot from reality")
        differing_positions = np.where(diff_to_prediction != 0)[0]
        all_differences = np.array([differing_positions, diff_to_prediction[differing_positions]]).transpose()
        if len(differing_positions) == 0:
            self.logger.info("model fits reality perfectly")
        else:
            self.logger.info("model differs from measurement in %i positions:\n%s", len(differing_positions), str(all_differences))

    def find_pois(self):
        periods = self.find_loop_periods()
        # The longer period is that of the zero loop.
        period_zero = periods[-1]
        candidates = self.find_poi_candidates(threshold=self.approx_tvalue)
        seqs = self.find_long_sequences(candidates, period_zero)
        leakpos, model_leakpos = self.cleanup_sequences(seqs, period_zero)
        self.verify_model(leakpos, model_leakpos)
        return leakpos, model_leakpos

def main():
    ta = Traces("seeda")
    tb = Traces("seedb")
    tr = Traces("seedr")
    tv = Traces("seedv")

    report_traces([ta, tb, tr, tv])

    # can also try tb instead of ta.
    analysis = FixedVsRandom(ta, tr)
    leakpos_a, model_leakpos_a = analysis.find_pois()
    
    analysis = FixedVsRandom(tb, tr)
    leakpos_b, model_leakpos_b = analysis.find_pois()

    if (leakpos_a == leakpos_b).all():
        self.logger.info("fixed vs random yields identical results for two fixed vectors")
    else:
        for i in range((PK_NROWS - 1) * params.n):
            if leakpos_a[i] != leakpos_b[i]:
                logger.info("leakpos differs: %s[%i] = %i vs %s[%i] = %i",
                    ta.name, i, leakpos_a[i],
                    tb.name, i, leakpos_b[i])

main()

