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
import json
from multiprocessing import Pool
load("support_splitting/support_splitting.sage")


params = SimpleNamespace(t=128, m=13, n=129)
RESULTDIR="results"
TRACES=10

# Some derived values: Number of rows, number of PoIs we expect
PK_NROWS = params.t * params.m

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger("mceliece")
logger.setLevel(logging.DEBUG)

file_handler = logging.FileHandler(os.path.join(RESULTDIR, '07-support-splitting.log'), mode='w')
formatter = logging.Formatter('%(asctime)s:%(levelname)s :%(name)s:%(message)s')
file_handler.setFormatter(formatter)
file_handler.setLevel(logging.DEBUG)
logger.addHandler(file_handler)

sys.path.append("mceliece-sage-20221023")
from parameters import parameters
from byterepr import to_privatekey, to_publickey, to_fieldelement
logger.info("loading McEliece parameters")
mc_params = parameters("mceliece8192128")


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

    def secretkey(self):
        return to_privatekey(self.sk, mc_params)

    def publickey(self):
        return to_publickey(self.pk, mc_params)

def recovered_points(i):
    fname_json = os.path.join(RESULTDIR, "victim-%02i.json" % i)
    with open(fname_json) as fh:
        recovery_data = json.load(fh)
    points = []
    for c in recovery_data["columns"]:
        alpha = struct.pack("<H", c["alpha"])
        beta = struct.pack("<H", c["g_of_alpha"])
        alpha = to_fieldelement(alpha, mc_params)
        beta = to_fieldelement(beta, mc_params)
        points.append((alpha, beta))
    return points

def makeColumn(input):
    g, a = input
    H_j = [0] * (mc_params.t * mc_params.m)
    a_int = a.integer_representation()
    a_bits = [(a_int >> i) & 1 for i in range(mc_params.m)]
    for i in range(mc_params.t):
        val = (a^i/g(a)).integer_representation()
        for b in range(mc_params.m):
            H_j[i*mc_params.m + b] = (val >> b) & 1
    return H_j

from itertools import repeat
from tqdm.contrib.concurrent import process_map
def goppa_code(g, alphas):
    #H = process_map(makeColumn, zip(repeat(g), alphas), max_workers=128)
    pool = Pool(256)
    H = pool.map(makeColumn, zip(repeat(g), alphas))
    H = matrix(GF(2), H).transpose()
    return H

def attack(i):
    fname_nfo = os.path.join(RESULTDIR, "victim-{:02}.nfo".format(i))
    key = McElieceKey(fname_nfo)
    sk = key.secretkey() # (delta,c,g,alpha,s)
    pk = key.publickey()
    points = recovered_points(i)

    # Verify the recovered points just to detect failures early.
    verify_g = sk[2]
    for i, (alpha, beta) in enumerate(points):
        verify_alpha = sk[3][i]
        verify_beta = verify_g(verify_alpha)
        if verify_alpha != alpha:
            logger.critical("column %i wrong alpha %s (expected %s)", i, repr(alpha), repr(verify_alpha))
        if verify_beta != beta:
            logger.critical("column %i wrong beta %s (expected %s)", i, repr(beta), repr(verify_beta))

    # Now we can generate the goppa polynomial and the standard goppa matrix
    logger.info("lagrange interpolation...")
    R.<x> = mc_params.Fq[]
    g = 0
    for i, (a_i, b_i) in enumerate(points):
        g_i = 1
        for j, (a_j, b_j) in enumerate(points):
            if i != j:
                g_i *= (x - a_j) / (a_i - a_j)
        g += b_i * g_i
    if g != verify_g:
        logger.critical("goppa polynomial differs:\ng_poly=%s\ng_vrfy=%s", g, verify_g)

    logger.info("generating parity check matrix of goppa code")
    alphas_unsorted = list(mc_params.Fq)
    H = goppa_code(g, alphas_unsorted)
    
    logger.info("bringing parity check matrix into systematic form")
    n = H.ncols()
    codim = H.nrows()
    H = H.echelon_form()
    H = H.transpose()

    permutation_sysForm = list(range(n))

    j = H.nrows()-1
    for i in range(codim):
        while H.submatrix(0,0,i+1,i+1).rank() < i+1:
            H.swap_rows(i,j)
            permutation_sysForm[j], permutation_sysForm[i] = permutation_sysForm[i], permutation_sysForm[j]
            j -= 1
            
    H = H.transpose().echelon_form()

    logger.info("starting support splitting")
    H_pk = identity_matrix(GF(2), codim).augment(pk)
    P = supportSplitting(H,H_pk)

    logger.info("sorting alphas")
    alphas = [0]*n
    for col in range(n):
        row = 0
        while P[row,col] == 0:
            row += 1
        alphas[col] = alphas_unsorted[permutation_sysForm[row]]

    # Verify correctness of solution
    for i in range(n):
        verify_alpha = sk[3][i]
        alpha = alphas[i]   

        if verify_alpha != alpha:
            logger.critical("column %i wrong alpha %s (expected %s)", i, repr(alpha), repr(verify_alpha))

    logger.info("recovered secret key succesfully. PASS")

def main():
    count = 10
    for i in range(count):
        logger.info("attacking key victim %i", i)
        attack(i)

main()

