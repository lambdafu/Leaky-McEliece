import encode
import fixedweight
import parameters

def encap_abstract(T,randombits,params):
  r'''
  Return the output of the abstract Classic McEliece Encap function
  on input T
  using the specified source of random bits
  for the specified parameters.

  "Abstract" means that this function does not include encodings
  of the inputs and outputs as byte strings.
  See encap() for the full function including encodings.

  INPUT:

  "T" - a matrix with mt rows and k columns

  "randombits" - a function that, on input r, returns a list of r bits

  "params" - a Classic McEliece parameter set
  '''

  assert isinstance(params,parameters.parameters)
  m = params.m
  n = params.n
  t = params.t
  k = params.k
  H = params.H

  assert T.nrows() == m*t
  assert T.ncols() == k

  # "Use FixedWeight to generate a vector e in F_2^n of weight t."
  e = fixedweight.fixedweight(randombits,params)
  assert len(e) == n
  assert sum(1 if ej else 0 for ej in e) == t

  if parameters.supportpc:
    # for testing against previous versions
    if params.pc:
      C0 = encode.encode(e,T,params)
      C1 = H((2,e))
      C = C0,C1
      K = H((1,e,C))
      return C,K

  # "Compute C = Encode(e,T)."
  C = encode.encode(e,T,params)

  # "Compute K = H(1,e,C)"
  K = H((1,e,C))

  # "Output ciphertext C and session key K."
  return C,K

import byterepr

def encap(T,randombytes,params):
  r'''
  Return the output of the Classic McEliece Encap function
  on input T
  using the specified source of random bytes
  for the specified parameters.

  This is the full function, including encodings
  of the inputs and outputs as byte strings.

  INPUT:

  "T" - a string of bytes representing a public key

  "randombytes" - a function that, on input r, returns an r-byte string

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)

  try:
    T = byterepr.to_publickey(T,params)
  except ValueError:
    # "Narrowly Decoded Classic McEliece rejects inputs
    #  (ciphertexts and public keys) where padding bits are nonzero;
    #  rejection means returning False."
    return False

  randombits = byterepr.randombits_from_randombytes(randombytes)
  C,K = encap_abstract(T,randombits,params)
  return byterepr.from_ciphertext(C,params),byterepr.from_sessionkey(K,params)

# ----- miscellaneous tests

def test1():
  import os
  import keygen

  def randombits(r):
    return [randrange(2) for j in range(r)]

  def randombytes(r):
    return os.urandom(r)

  for system in parameters.alltests:
    P = parameters.parameters(system,allowtestparams=True)
    n = P.n
    m = P.m
    t = P.t
    l = P.l

    print('encap_abstract %s' % system)
    sys.stdout.flush()

    T,privatekey = keygen.keygen_abstract(randombits,P)
    C,K = encap_abstract(T,randombits,P)
    if parameters.supportpc and P.pc:
      # for testing against previous versions
      C0,C1 = C
      assert len(C0) == m*t
      assert len(C1) == l
    else:
      assert len(C) == m*t
    assert len(K) == l
    
    print('encap %s' % system)
    sys.stdout.flush()

    T,privatekey = keygen.keygen(randombytes,P)
    C,K = encap(T,randombytes,P)
    if parameters.supportpc and P.pc:
      # for testing against previous versions
      assert len(C) == ceil(m*t/8)+ceil(l/8)
    else:
      assert len(C) == ceil(m*t/8)
    assert len(K) == ceil(l/8)
      
if __name__ == '__main__':
  test1()
