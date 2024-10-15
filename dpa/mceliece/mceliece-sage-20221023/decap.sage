import decode
import parameters

def decap_abstract(C,priv,params):
  r'''
  Return the output of the abstract Classic McEliece Decap function
  on input (C,priv)
  for the specified parameters.

  "Abstract" means that this function does not include encodings
  of the inputs and outputs as byte strings.
  See decap() for the full function including encodings.

  INPUT:

  "C" - a ciphertext

  "priv" - a private key returned by keygen.keygen()

  "params" - a Classic McEliece parameter set
  '''

  assert isinstance(params,parameters.parameters)
  m = params.m
  n = params.n
  t = params.t
  H = params.H

  # "Extract s in F_2^n and Gammaprime = (G,alpha'_0,alpha'_1,...,alpha'_{n-1}) from the private key."
  delta,c,g,alpha,s = priv
  alphaprime = alpha[:n]
  Gammaprime = tuple([g]+alphaprime)

  if parameters.supportpc:
    # for testing against previous versions
    if params.pc:
      C0,C1 = C
      C0 = list(C0)
      C1 = list(C1)
      e = decode.decode(C0,Gammaprime,params)
      if e == False: return H((0,s,C))
      C1prime = H((2,e))
      if C1prime != C1: return H((0,s,C))
      return H((1,e,C))

  # "takes as input a ciphertext C"
  C = list(C)
  assert len(C) == m*t

  # "Set b <- 1."
  b = 1

  # "Compute e <- Decode(C,Gammaprime)."
  e = decode.decode(C,Gammaprime,params)

  # "If e = False, set e <- s and b <- 0."
  if e == False:
    e = s
    b = 0

  # "Compute K = H(b,e,C)"
  K = H((b,e,C))

  # "Output session key K."
  return K

import byterepr

def decap(C,priv,params):
  r'''
  Return the output of the Classic McEliece Decap function
  on input (C,priv)
  for the specified parameters.

  This is the full function, including encodings
  of the inputs and outputs as byte strings.

  INPUT:

  "C" - a byte string that represents a ciphertext

  "priv" - a byte string that represents a private key

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)

  try:
    C = byterepr.to_ciphertext(C,params)
  except ValueError:
    # "Narrowly Decoded Classic McEliece rejects inputs
    #  (ciphertexts and public keys) where padding bits are nonzero;
    #  rejection means returning False."
    return False

  priv = byterepr.to_privatekey(priv,params)

  K = decap_abstract(C,priv,params)
  return byterepr.from_sessionkey(K,params)

# ----- miscellaneous tests

def test1():
  import keygen
  import encap
  import os

  def randombits(r):
    return [randrange(2) for j in range(r)]

  def randombytes(r):
    return os.urandom(r)

  for system in parameters.alltests:
    P = parameters.parameters(system,allowtestparams=True)

    print('decap_abstract %s' % system)
    sys.stdout.flush()

    T,priv = keygen.keygen_abstract(randombits,P)
    C,K = encap.encap_abstract(T,randombits,P)
    assert decap_abstract(C,priv,P) == K

    print('decap %s' % system)
    sys.stdout.flush()

    T,priv = keygen.keygen(randombytes,P)
    C,K = encap.encap(T,randombytes,P)
    assert decap(C,priv,P) == K

    # not tested here: decap on invalid ciphertexts

if __name__ == '__main__':
  test1()
