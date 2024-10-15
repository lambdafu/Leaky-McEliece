#!/bin/env python

import chipwhisperer as cw
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

params = SimpleNamespace(t=128, m=13, n=128+8)
TRACEDIR="traces"
TRACES=1000
VICTIM_TRACES=10
MOCK=False

logger = logging.getLogger("mceliece")
logger.setLevel(logging.DEBUG)

class McElieceKey(object):
    def __init__(self, seed=None, fname=None):
        self.SEEDLEN = 48
        if fname is not None:
            self.load(fname)
        else:
            self.set(seed)

    def set(self, seed):
        if seed is None:
            seed = np.random.bytes(self.SEEDLEN)
        if len(seed) != self.SEEDLEN:
            logger.error("seed wrong length")
            return False
        self.seed = seed
        seed_str = ''.join(format(byte, '02x') for byte in seed)
        cmd = ["docker", "run", "--rm", "-ti", "mceliece", "/leak-mceliece8192128", seed_str]
        try:
            result = subprocess.run(cmd, check=True, text=True, capture_output=True)
            output = result.stdout
        except subprocess.CalledProcessError as e:
            output = e.output
            err = e.stderr.strip()
            logger.error("docker failed with error: %s", err)
            raise RuntimeError(f"Docker command failed with error: {e.stderr.strip()}")
        self.output = output
        self.loads(output)
        return True

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

    def save(self, fname):
        with open(fname, 'w') as fh:
            fh.write(self.output)

    def load(self, fname):
        with open(fname, 'r') as fh:
            self.output = fh.read()
        self.loads(self.output)

class Scope(object):
    def __init__(self):
        self.PLATFORM = "CW308_STM32L5"
        self.SAMPLES = 287002932
        self.FREQUENCY = 10000000
        # We can upload 4 lines at once through the s command.
        self.STORE_NROWS = 4
        self.scope = None
        self.target = None

    def reset(self):
        if self.scope is None:
            logger.error("not connected")
            return
        logger.info("reset")
        self.scope.errors.clear()
        self.scope.io.nrst = 'low'
        time.sleep(0.25)
        self.scope.io.nrst = 'high_z'
        time.sleep(0.25)

    def connect(self):
        if self.scope is not None:
            logger.error("already connected")
            return
        self.scope = cw.scope()
        logger.info('fw version=%s', self.scope.fw_version)
        self.target = cw.target(self.scope)
        self.scope.default_setup()
        self.scope.gain.db += 9
        self.scope.clock.adc_mul = 1
        self.scope.adc.clip_errors_disabled = True
        self.scope.adc.stream_mode = True
        self.scope.adc.timeout = 140
        self.scope.adc.samples = self.SAMPLES
        freq_ratio = self.scope.clock.clkgen_freq / self.target.baud
        self.scope.clock.clkgen_freq = self.FREQUENCY
        self.target.baud = self.FREQUENCY / freq_ratio
        self.reset()

    def disconnect(self):
        if self.scope is None:
            logger.error("not connected")
            return
        self.scope.errors.clear()
        self.target.dis()
        self.scope.dis()
        self.target, self.scope = None, None

    # The mceliece firmware has a checksum function in C that corresponds to this checksum function
    # https://stackoverflow.com/questions/35205702/calculating-crc16-in-python
    @staticmethod
    def crc16(data):
        data = bytes(data)
        crc = 0xFFFF
        for i in range(len(data)):
            crc ^= data[i] << 8
            for j in range(0,8):
                if (crc & 0x8000) > 0:
                    crc =(crc << 1) ^ 0x1021
                else:
                    crc = crc << 1
        return crc & 0xFFFF

    @staticmethod
    def _mat_to_binmat(mat):
        H = np.zeros(shape=(params.t * params.m, params.n), dtype='bool')
        for j in range(params.m * params.t):
            for c in range(params.n // 8):
                for b in range(8):
                    if (mat[j, c] >> b) & 1:
                       H[j, c*8 + b] = 1
        return H

    @staticmethod
    def _binmat_to_mat(H):
        # Now calculate the bytes of the matrix.
        mat = np.zeros(shape=(params.m * params.t, params.n // 8), dtype=np.uint8)
        for j in range(params.m * params.t):
            for c in range(params.n // 8):
                for b in range(8):
                    if H[j, c*8 + b]:
                        mat[j,c] |= 1 << b
        return mat

    @staticmethod
    def _binmat_ge(H):
        H = np.copy(H)
        for j in range(params.t + 1):
            for i in range(j + 1, params.m * params.t):
                mask = H[j, j] ^ H[i, j]
                if mask:
                    H[j] ^= H[i]
            assert(H[j,j] == True)
            for i in range(params.m * params.t):
                if i != j:
                    mask = H[i, j]
                    if mask:
                        H[i] ^= H[j]
        return H

    @staticmethod
    def _mat_ge(mat):
        H = Scope._mat_to_binmat(mat)
        H = Scope._binmat_ge(H)
        return Scope._binmat_to_mat(H)

    # Run gauss. elim. and verify the result with the device checksum.
    def verify_final_checksum(self, mat, device_chksum):
        mat2 = Scope._mat_ge(mat)
        mat2chk = Scope.crc16(mat2)
        if mat2chk == device_chksum:
            return True

        logger.error("final checksum mismatch!")
        matd = self.download_matrix()
        np.savez("chksumerror",
                 mat=mat, matchk=Scope.crc16(mat),
                 mat2=mat2, mat2chk=mat2chk,
                 matd=matd, matdchk=Scope.crc16(matd),
                 device=device_chksum)
        return False

    # Upload the matrix, and read it back for verification if verify is True.
    def upload_matrix(self, mat, verify=False):
        with logging_redirect_tqdm():
            for r in tqdm(range(0, params.t * params.m, self.STORE_NROWS), desc="uploading matrix data"):
                rowdata = bytes(mat[r:r + self.STORE_NROWS])
                padding = b''
                wanted_len = self.STORE_NROWS * params.n // 8
                if len(rowdata) < wanted_len:
                    padding = b'\x00' * (wanted_len - len(rowdata))
                header = struct.pack("!H", r)
                msg = header + rowdata + padding
                self.target.simpleserial_write('s', msg)
                response = self.target.simpleserial_read('r', 2, timeout=250)
                written, = struct.unpack("!H", response)
                if written != len(rowdata):
                    logger.error("wrote %i bytes instead of %i", written, len(rowdata))

                # Read back and verify (this is an alternative to checksum verification)
                if verify:
                    self.target.simpleserial_write('l', header)
                    response = self.target.simpleserial_read('r', 2 + wanted_len, timeout=1000)
                    gotten, = struct.unpack("!H", response[:2])
                    if gotten != written:
                        logger.error("got %i bytes instead of %i", gotten, written)
                    if response[2:2+gotten] != rowdata:
                        logger.error("got wrong data on readback: %s instead of %s", repr(response[2:2+gotten]), repr(rowdata))
                    if response[2+gotten:] != padding:
                        logger.error("got wrong padding on readback: %s instead of %s", repr(response[2+gotten:]), repr(padding))

    # after uploading the matrix, checksum it
    def verify_checksum(self, mat):
        matchk=Scope.crc16(mat)
        chksum = struct.pack("!H", matchk)
        self.target.simpleserial_write('c', b'')
        target_chksum = self.target.simpleserial_read('r', 2, timeout=4000)
        if chksum != target_chksum:
            logger.error("chksum error: %s instead of %s", repr(target_chksum), repr(chksum))

    def download_matrix(self):
        devmat = b''
        with logging_redirect_tqdm():
            for r in tqdm(range(0, params.t * params.m, self.STORE_NROWS), desc="downloading matrix data"):
                header = struct.pack("!H", r)
                wanted_len = self.STORE_NROWS * params.n // 8
                self.target.simpleserial_write('l', header)
                response = self.target.simpleserial_read('r', 2 + wanted_len, timeout=1000)
                gotten, = struct.unpack("!H", response[:2])
                devmat += response[2:2+gotten]
        return devmat

    def capture_trace(self, mat):
        self.upload_matrix(mat)
        self.verify_checksum(mat)
        self.scope.arm()
        # Launch gaussian elimination process.
        self.target.simpleserial_write('g', b'')
        ret = self.scope.capture(poll_done=True)
        if ret:
            logger.error("timeout during capture")
            return None
        # Long timeout, see https://forum.newae.com/t/trigger-error-in-lab-2-1/3340/11
        target_chksum2 = self.target.simpleserial_read('r', 2, timeout=20000)
        target_chksum2, = struct.unpack("!H", target_chksum2)
        #logger.info("final checksum: %i", target_chksum2)
        passed = self.verify_final_checksum(mat, target_chksum2)
        if not passed:
            logger.error("final checksum verification failed")
            return None
        wave = self.scope.get_last_trace()
        return wave
    
    def my_capture_trace(self, key):
        # We only need to first t+1 columns of the matrix H, rounded up to multiple of 8.
        mat = np.zeros(shape=(params.m * params.t, params.n // 8), dtype=np.uint8)
        key_row_len = len(key.mat) // (params.t * params.m)
        row_len = params.n // 8
        for r in range(params.t * params.m):
            mat[r] = [x for x in key.mat[r * key_row_len:r * key_row_len + row_len]]
        return self.capture_trace(mat)

# Traces keeps a series of traces in an experiment.
# scope: A scope object.
# name: A name for the trace series.
# random_seed: if True, use a fresh random seed for each capture, else use the same (also random) seed
# save_all: if True, saves each trace in a dedicated file. otherwise, only mean, variance and selected traces are saved.
class Traces(object):
    def __init__(self, scope, name, fresh_keys=False, save_all=False):
        self.scope = scope
        self.name = name
        self.save_all = save_all
        self.fname_base = os.path.join(TRACEDIR, name)
        self.fname = self.fname_base + ".npy"
        self.fname_nfo = self.fname_base + ".nfo"
        self.fname_avg = self.fname_base + "-ravg.npy"
        # welford's online m2: https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
        self.fname_wm2 = self.fname_base + "-rwm2.npy"
        if os.path.exists(self.fname):
            self.load()
        else:
            logger.info("creating experiment file %s", self.fname)
            self.init(fresh_keys=fresh_keys)
            self.save(backup=False)

    def init(self, fresh_keys=False):
        self.count = 0
        # average, variance
        self.avg = np.zeros(scope.SAMPLES, dtype=np.float64)
        self.wm2 = np.zeros(scope.SAMPLES, dtype=np.float64)
        if not fresh_keys:
            self.key = McElieceKey()
        else:
            self.key = None

    def load(self):
        self.count = int(np.load(self.fname))
        if os.path.exists(self.fname_nfo):
            self.key = McElieceKey(fname=self.fname_nfo)
        else:
            self.key = None
        self.avg = np.load(self.fname_avg)
        self.wm2 = np.load(self.fname_wm2)

    def backup(self):
        shutil.copy(self.fname, self.fname + ".bak")
        shutil.copy(self.fname_avg, self.fname_avg + ".bak")
        shutil.copy(self.fname_wm2, self.fname_wm2 + ".bak")

    def save(self, backup=True):
        if backup:
            self.backup()
        np.save(self.fname, self.count)
        if self.key:
            self.key.save(self.fname_nfo)
        np.save(self.fname_avg, self.avg)
        np.save(self.fname_wm2, self.wm2)

    def add_trace(self, auto_save=True):
        logger.info("capturing %s/%i", self.name, self.count)
        fname_base = "{}-t{:04}".format(self.fname_base, self.count)
        fname = fname_base + ".npy"
        fname_nfo = fname_base + ".nfo"
        key = self.key
        if key is None:
            key = McElieceKey()
            key.save(fname_nfo)
        wave = self.scope.my_capture_trace(key)
        if len(wave) != self.scope.scope.adc.samples:
            logger.warning("wrong sample count in %s/%i: %i", self.name, self.count, len(wave))
        if self.scope.scope.adc.samples != self.scope.scope.adc.trig_count:
            logger.warning("trig_count mismatch in trace %s/%i: %i", self.name, self.count, self.scope.scope.adc.trig_count)
        if wave.min() < -0.499 or wave.max() > 0.499:
            logger.warning("clipping in trace %s/%i: %f, %f", self.name, self.count, wave.min(), wave.max())
        trunc_wave = wave.astype('float16')
        np.save(fname, trunc_wave)
        if not self.save_all and self.count >= 10:
            old_fname_base = "{}-t{:04}".format(self.fname_base, self.count - 10)
            old_fname = old_fname_base + ".npy"
            old_fname_nfo = fname_base + ".nfo"
            if os.path.exists(old_fname):
                os.remove(old_fname)
            if os.path.exists(old_fname_nfo):
                os.remove(old_fname_nfo)
        # Update count
        self.count += 1
        # Update avg and welford's m2.
        delta = wave - self.avg
        self.avg += delta / self.count
        delta2 = wave - self.avg
        self.wm2 += delta * delta2
        if auto_save:
            self.save()

    def capture(self, nr):
        if nr < self.count:
            logger.warning("skipping existing trace: %s/%i", self.name, nr)
            return
        if nr != self.count:
            logger.warning("refusing to generate non-continuous trace: %s/%i", self.name, nr)
            return
        self.add_trace()

class MockScope(object):
    def __init__(self):
        self.SAMPLES = 10
        adc = SimpleNamespace(samples=self.SAMPLES, trig_count=self.SAMPLES)
        self.scope = SimpleNamespace(adc=adc,errors=None)
    def my_capture_trace(self, key):
        wave = np.random.normal(0, 0.15, self.SAMPLES).astype(np.float32)
        # To make sure that resuming measurements work, we add a seed-based fixed value here.
        wave[0] = key.mat[0]
        wave[1] = key.mat[1]
        return wave
    def connect(self):
        pass
    def disconnect(self):
        pass
    
np.random.seed()
if MOCK:
    scope = MockScope()
else:
    scope = Scope()
scope.connect()

ta = Traces(scope, "seeda")
tb = Traces(scope, "seedb")
tr = Traces(scope, "seedr", fresh_keys=True, save_all=True)

traces = [ta, tb, tr]
for trace in traces:
    trace.save()
min_count = min([t.count for t in traces])
if min_count < TRACES:
    logger.warning("Resuming from %i", min_count)


with logging_redirect_tqdm():
    pbar = tqdm(range(min_count, TRACES), desc="capture traces")
    for i in pbar:
        # Use this to stop the process instead of CTRL-C to avoid USB hangups.        
        if os.path.exists("STOP"):
            logger.warning("stopping due to user request")
            break

        for t in traces:
            pbar.set_description(t.name)
            t.capture(i)

            if os.path.exists("STOP"):
                break

# Victim traces. We do not care about their average and mean, but rather want
# the individual traces.
tv = Traces(scope, "seedv", fresh_keys=True, save_all=True)
with logging_redirect_tqdm():
    pbar = tqdm(range(tv.count, VICTIM_TRACES), desc="capture victim traces")
    for i in pbar:
        # Use this to stop the process instead of CTRL-C to avoid USB hangups.        
        if os.path.exists("STOP"):
            logger.warning("stopping due to user request")
            break

        tv.capture(i)

logger.info(scope.scope.errors)
scope.disconnect()

