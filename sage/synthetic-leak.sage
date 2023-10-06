import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "mceliece-sage-20221023"))

from parameters import parameters
from byterepr import from_privatekey, from_ciphertext, from_sessionkey, from_publickey, from_vector, from_fieldelement
from byterepr import to_privatekey, to_ciphertext, to_sessionkey, to_publickey, to_vector, to_fieldelement
from byterepr import from_fieldelement
from decap import decap_abstract
from keygen import keygen_abstract

def mceliece_leak(kem):
    params = parameters(kem, allowtestparams=True)
    publickey, secretkey = keygen_abstract(randombits, params)
    g=secretkey[2]
    alpha=secretkey[3]

    H = matrix(GF(2), params.t * params.m, params.n)
    for j, a in enumerate(alpha):
        a_int = a.integer_representation()
        a_bits = [(a_int >> i) & 1 for i in range(params.m)]
        for i in range(params.t):
            val = (a^i/g(a)).integer_representation()
            for b in range(params.m):
                H[i*params.m + b, j] = (val >> b) & 1

    leak = []
    P = copy(H)
    for j in range(params.m * params.t):
        for i in range(j + 1, params.m * params.t):
            mask = P[j, j] + P[i, j]
            P.add_multiple_of_row(j, i, mask)
        assert(P[j,j] == 1)
        for i in range(params.m * params.t):
            if i != j:
                mask = P[i, j]
                P.add_multiple_of_row(i, j, mask)
                leak.append(mask)
    P1 = P[:,:params.m*params.t]
    P2 = P[:, params.m*params.t:]
    public_key_correct = P1 == identity_matrix(GF(2), params.m*params.t), P2 == publickey
    assert(public_key_correct)
    leak_correct = len(leak) == (params.m*params.t)*(params.m*params.t-1)
    assert(leak_correct)

    sk = from_privatekey(secretkey, params)
    pk = from_publickey(publickey, params)
    return sk, pk, leak

import random
def randombits(n):
    return [random.getrandbits(1) for i in range(n)]

def generate_leak(seed, kem):
    seed = 0
    if seed == 0:
        set_random_seed()
        seed = initial_seed()
    set_random_seed(seed)
    if kem.startswith("mceliece"):
        sk, pk, leak = mceliece_leak(kem)
    else:
        sk, pk, leak = botan_leak(kem)
    return { 'kem': kem, 'seed': seed, 'sk': sk, 'pk': pk, 'leak': leak }

def as_hex(raw):
    return ''.join(['{:02X}'.format(x) for x in raw])

from collections import OrderedDict
def write_leak_file(fname, result):
    content = OrderedDict()
    content["kem"] = result['kem']
    content["comment"] = "synthetic leak generated with synthetic-leak.sage"
    content["seed"] = result['seed']
    content["sk"] = as_hex(result['sk'])
    content["pk"] = as_hex(result['pk'])
    content["leak"] = ''.join(["01" if x else "00" for x in result['leak']])

    with open(fname, "w") as fh:
        for k, v in content.items():
            fh.write("{} = '{}'\n".format(k, v))

#kem = "mceliece51220"
#kem = "mceliece102450"
kem = sys.argv[-2]
fname = sys.argv[-1]
#"test-" + kem + ".leak"
seed = 0
result = generate_leak(seed, kem)
write_leak_file(fname, result)
