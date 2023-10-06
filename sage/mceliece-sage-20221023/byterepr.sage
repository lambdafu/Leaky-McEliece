import parameters
import controlbits

def from_vector(v):
  r'''
  Return the byte string representing bit vector v.

  INPUT:

  "v" - a list of bits
  '''

  v = [int(vj) for vj in v]
  assert all(vj in (0,1) for vj in v)

  r = len(v)
  expectedlen = ceil(r/8)

  if r%8 != 0:
    # "If r is not a multiple of 8
    #  then an r-bit vector v = (v_0,v_1,...,v_{r-1}) in F_2^r
    #  is zero-padded on the right to length between r+1 and r+7,
    #  whichever is a multiple of 8"
    padding = [0]*(8-(r%8))
    assert len(padding) >= 1
    assert len(padding) <= 7
    v = v+padding
    assert len(v)%8 == 0
    # "and then represented as above"
    r = len(v)

  assert r == len(v)
  assert r%8 == 0
  assert expectedlen == r/8

  # "the following sequence of r/8 bytes:
  #  (v_0+2v_1+4v_2+...+128v_7,
  #   v_8+2v_9+4v_10+...+128v_15,
  #   ...,
  #   v_{r-8}+2v_{r-7}+4v_{r-6}+...+128v_{r-1})"
  b = [sum(int(v[j*8+i])<<i for i in range(8)) for j in range(r/8)]
  assert len(b) == expectedlen
  return bytes(bytearray(b))

def to_vector(b,r,narrowly=True):
  r'''
  Return the list of r bits
  represented by byte string b.
  Raise an exception if b has the wrong length.

  INPUT:

  "b" - bytes

  "r" - a nonnegative integer

  "narrowly" (optional, default True) -
  raise a ValueError exception if there are padding bits with nonzero value
  '''
  assert parent(b) == bytes
  r = ZZ(r)
  assert r >= 0
  assert narrowly in (True,False)

  b = list(bytearray(b))
  assert len(b) == ceil(r/8)

  v = [1&(bi>>j) for bi in b for j in range(8)]
  assert len(v) == 8*len(b)

  if narrowly:
    # "Narrowly Decoded Classic McEliece rejects inputs ...
    #  where padding bits are nonzero"
    if any(vi != 0 for vi in v[r:]):
      raise ValueError('padding bits do not all have value 0')
  v = v[:r]

  assert all(vj in (0,1) for vj in v)
  assert len(v) == r
  return v

def randombits_from_randombytes(randombytes):
  def randombits(r):
    # parameter sets allowed in this implementation have r%8 == 0
    assert r%8 == 0
    v = randombytes(r//8)
    return to_vector(v,r)
  return randombits

def from_sessionkey(K,params):
  r'''
  Return the byte string representing session key K.

  INPUT:

  "K" - a list of l bits

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  l = params.l

  # "A session key K is an element of F_2^l.
  #  It is represented as a ceil(l/8)-byte string."
  K = list(K)
  assert len(K) == l
  b = from_vector(K)
  assert len(b) == ceil(l/8)
  return b

# not used in tests currently
def to_sessionkey(b,params):
  r'''
  Return the session key represented by byte string b.

  INPUT:

  "b" - a byte string

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  l = params.l
  assert len(b) == ceil(l/8)
  K = to_vector(b,l)
  assert len(K) == l
  return K

def from_ciphertext(C,params):
  r'''
  Return the byte string representing ciphertext C.

  INPUT:

  "C" - a list of m*t bits

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  m = params.m
  t = params.t

  if parameters.supportpc:
    # for testing against previous versions
    if params.pc:
      C0,C1 = C
      return from_vector(C0)+from_vector(C1)

  # "A ciphertext C is an element of F_2^{mt}.
  #  It is represented as a ceil(mt/8)-byte string."
  C = list(C)
  assert len(C) == m*t
  b = from_vector(C)
  assert len(b) == ceil(m*t/8)
  return b

def to_ciphertext(b,params):
  r'''
  Return the ciphertext represented by byte string b.

  INPUT:

  "b" - a byte string

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  m = params.m
  t = params.t

  if parameters.supportpc:
    # for testing against previous versions
    if params.pc:
      l = params.l
      b0,b1 = b[:ceil(m*t/8)],b[ceil(m*t/8):]
      return to_vector(b0,m*t),to_vector(b1,l)

  assert len(b) == ceil(m*t/8)
  C = to_vector(b,m*t)
  assert len(C) == m*t
  return C

def from_hashinput(x,params):
  r'''
  Return the byte string representing hash input x.

  INPUT:

  "x" - a hash input

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  n = params.n
  m = params.m
  t = params.t

  if parameters.supportpc:
    # for testing against previous versions
    if params.pc:
      if x[0] == 2:
        assert len(x) == 2
        return bytes(bytearray([2]))+from_vector(x[1])
      u,v,C = x
      u = int(u)
      v = list(v)
      return bytes(bytearray([u]))+from_vector(v)+from_ciphertext(C,params)

  # "There are two types of hash inputs: (1,v,C), and (0,v,C).
  #  Here v in F_2^n, and C is a ciphertext."
  u,v,C = x
  u = int(u)
  assert u in (0,1)
  v = list(v)
  assert len(v) == n

  # "The initial 0 or 1 is represented as a byte."
  result = bytes(bytearray([u]))

  # "The vector v is represented as the next ceil(n/8) bytes."
  result += from_vector(v)

  # "The ciphertext is represented as the next ceil(mt/8) bytes."
  result += from_ciphertext(C,params)

  # "All hash inputs thus begin with byte 0 or 1"
  assert bytearray(result)[0] in (0,1)

  assert len(result) == 1+ceil(n/8)+ceil(m*t/8)
  return result

def from_publickey(T,params):
  r'''
  Return the byte string representing public key T.

  INPUT:

  "T" - an mt x k matrix over F_2

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  m = params.m
  t = params.t
  k = params.k

  # "The public key T, which is an mt x k matrix,
  #  is represented in a row-major fashion.
  #  Each row of T is represented as a ceil(k/8)-byte string,
  #  and the public key is represented as
  #  the mt ceil(k/8)-byte concatenation of these strings."
  assert T.nrows() == m*t
  assert T.ncols() == k
  b = b''.join(from_vector(Ti) for Ti in T)
  assert len(b) == m*t*ceil(k/8)
  return b

def to_publickey(b,params):
  r'''
  Return the public key represented by byte string b.

  INPUT:

  "b" - a byte string

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  m = params.m
  t = params.t
  k = params.k
  assert len(b) == m*t*ceil(k/8)
  rows = [to_vector(b[j*ceil(k/8):(j+1)*ceil(k/8)],k) for j in range(m*t)]
  T = matrix(GF(2),rows)
  assert T.nrows() == m*t
  assert T.ncols() == k
  return T

def from_fieldelement(u,params):
  r'''
  Return the byte string representing field element u.

  INPUT:

  "u" - an element of F_q

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  m = params.m

  uint = u.integer_representation()
  ubits = [1&(uint>>j) for j in range(m)]
  b = from_vector(ubits)
  assert len(b) == ceil(m/8)
  return b

def to_fieldelement(b,params):
  r'''
  Return the field element represented by byte string b.

  INPUT:

  "b" - a byte string

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  m = params.m
  Fq = params.Fq
  assert len(b) == ceil(m/8)
  ubits = to_vector(b,m)
  u = Fq(ubits)
  assert from_fieldelement(u,params) == b
  return u

def from_poly(g,params):
  r'''
  Return the byte string representing monic degree-t polynomial g.

  INPUT:

  "g" - a monic degree-t polynomial in F_q[x]

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  m = params.m
  t = params.t

  # "The monic-irreducible degree-t polynomial
  #  g = g_0 + g_1 x + ... + g_{t-1} x^{t-1} + x^t
  #  is represented as t ceil(m/8) bytes,
  #  namely the concatenation of the representations
  #  of the field elements g_0,g_1,...,g_{t-1}."
  assert g.degree() == t
  assert g.is_monic()
  # this encoding doesn't care about irreducibility
  b = b''.join(from_fieldelement(g[i],params) for i in range(t))
  assert len(b) == t*ceil(m/8)
  return b

def to_poly(b,params):
  r'''
  Return the monic degree-t polynomial represented by byte string b.

  INPUT:

  "b" - a byte string

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  m = params.m
  t = params.t
  Fq = params.Fq
  Fqx.<x> = Fq[]
  assert len(b) == t*ceil(m/8)
  coeffs = [to_fieldelement(b[i*ceil(m/8):(i+1)*ceil(m/8)],params) for i in range(t)]
  g = Fqx(coeffs+[1])
  assert g.degree() == t
  assert g.is_monic()
  return g

def from_fieldordering(alpha,params):
  r'''
  Return the byte string representing field ordering alpha.

  INPUT:

  "alpha" - a list of q distinct elements of F_q

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  m = params.m
  q = params.q
  Fq = params.Fq

  # "a sequence (alpha_0,...,alpha_{q-1})
  #  of q distinct elements of F_q"
  alpha = list(alpha)
  assert len(alpha) == q
  assert len(set(alpha)) == q
  assert sorted(alpha) == sorted(Fq)

  # "Define pi as the permutation of {0,1,...,q-1}
  #  such that alpha_i = sum_{j=0}^{m-1} pi(i)_j * z^{m-1-j}
  #  for all i in {0,1,...,q-1}."
  alphaint = [alpha[i].integer_representation() for i in range(q)]
  pi = [sum(((alphaint[i]>>(m-1-j))&1)<<j for j in range(m)) for i in range(q)]
  
  # "The ordering (alpha_0,...,alpha_{q-1}) is represented
  #  as a sequence of (2m-1)2^(m-1) control bits
  #  for an in-place Benes network for pi. ...
  #  This document requirs that a permutation pi be converted to 
  #  specifically the control bits defined by controlbits in Figure 1."
  c = controlbits.controlbits(pi)
  assert len(c) == (2*m-1)<<(m-1)

  # "This vector is represented as ceil((2m-1)2^{m-4}) bytes as above."
  b = from_vector(c)
  assert len(b) == ceil((2*m-1)*2^(m-1)/8)

  return b

def to_fieldordering(b,params):
  r'''
  Return the field ordering represented by byte string b.

  INPUT:

  "b" - a byte string

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  m = params.m
  q = params.q
  Fq = params.Fq

  assert len(b) == ceil((2*m-1)*2^(m-1)/8)
  c = to_vector(b,(2*m-1)<<(m-1))
  pi = controlbits.permutation(c)
  alpha = [Fq([1&(pi[i]>>(m-1-j)) for j in range(m)]) for i in range(q)]

  assert len(alpha) == q
  assert len(set(alpha)) == q
  assert sorted(alpha) == sorted(Fq)
  return alpha

def from_columnselection(c,params):
  r'''
  Return the byte string representing column selection c.

  INPUT:

  "c" - a list of mu integers in increasing order between mt-mu and mt-mu+nu-1

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  m = params.m
  t = params.t
  mu = params.mu
  nu = params.nu

  # "a sequence c = (c_{mt-mu},...,c_{mt-1}) of mu integers
  #  in increasing order between mt-mu and mt-mu+nu-1."
  c = [int(cj) for cj in c]
  assert len(c) == mu
  assert c == sorted(set(c))
  assert all(cj >= m*t-mu for cj in c)
  assert all(cj <= m*t-mu+nu-1 for cj in c)

  # "represented as a ceil{nu/8}-byte string,
  #  the little-endian format of the integer
  #  sum_{i=0}^{mu-1} 2^{c_{mt-mu+i}-(mt-mu)}."
  cint = sum(2^(ci-(m*t-mu)) for ci in c)
  b = from_vector([1&(cint>>j) for j in range(nu)])
  assert len(b) == ceil(nu/8)

  # "However, for (mu,nu) = (0,0),
  #  the sequence is instead represented as
  #  the 8-byte string which is the little-endian format
  #  of 2^32-1, i.e., 4 bytes of value 255 followed by 4 bytes of value 0."
  if (mu,nu) == (0,0):
    assert len(b) == 0
    b = bytes(bytearray([255,255,255,255,0,0,0,0]))
    assert len(b) == 8

  return b

def to_columnselection(b,params):
  r'''
  Return the column selection represented by byte string b.

  INPUT:

  "b" - the byte string (255,255,255,255,0,0,0,0) if (mu,nu) = (0,0);
  otherwise a byte string of length ceil(nu/8)
  with exactly mu bits set

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  m = params.m
  t = params.t
  mu = params.mu
  nu = params.nu

  if (mu,nu) == (0,0):
    assert b == bytes(bytearray([255,255,255,255,0,0,0,0]))
    return []

  assert len(b) == ceil(nu/8)
  cvec = to_vector(b,nu)
  c = [m*t-mu+j for j in range(nu) if cvec[j] != 0]
  assert len(c) == mu
  assert c == sorted(set(c))
  assert all(cj >= m*t-mu for cj in c)
  assert all(cj <= m*t-mu+nu-1 for cj in c)
  return c

def from_privatekey(priv,params):
  r'''
  Return the byte string representing private key priv.

  INPUT:

  "priv" - a private key (delta,c,g,alpha,s)

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  m = params.m
  n = params.n
  t = params.t
  l = params.l
  mu = params.mu
  nu = params.nu

  # "A private key (delta,c,g,alpha,s) is represented
  #  as the concatenation of five parts:"
  delta,c,g,alpha,s = priv

  # "The ceil(l/8)-byte string representing delta in F_2^l."
  assert len(delta) == l
  b0 = from_vector(delta)
  assert len(b0) == ceil(l/8)

  # "The string representing the column selections c.
  #  This string has ceil(nu/8) bytes,
  #  or 8 bytes of (mu,nu) = (0,0)."
  b1 = from_columnselection(c,params)
  assert len(b1) == 8 if (mu,nu) == (0,0) else ceil(nu/8)

  # "The t ceil(m/8)-byte string representing the polynomial g."
  b2 = from_poly(g,params)
  assert len(b2) == t*ceil(m/8)

  # "The ceil((2m-1)2^(m-4)) bytes representing the field ordering alpha."
  b3 = from_fieldordering(alpha,params)
  assert len(b3) == ceil((2*m-1)*2^(m-4))

  # "The ceil(n/8)-byte string representing s in F_2^n."
  assert len(s) == n
  b4 = from_vector(s)
  assert len(b4) == ceil(n/8)

  b = b0+b1+b2+b3+b4
  return b

def to_privatekey(b,params):
  r'''
  Return the private key represented by byte string b.

  INPUT:

  "b" - a byte string

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  m = params.m
  n = params.n
  t = params.t
  l = params.l
  mu = params.mu
  nu = params.nu

  b0len = ceil(l/8)
  b1len = 8 if (mu,nu) == (0,0) else ceil(nu/8)
  b2len = t*ceil(m/8)
  b3len = ceil((2*m-1)*2^(m-4))
  b4len = ceil(n/8)

  assert len(b) == b0len+b1len+b2len+b3len+b4len
  b0,b = b[:b0len],b[b0len:]
  b1,b = b[:b1len],b[b1len:]
  b2,b = b[:b2len],b[b2len:]
  b3,b = b[:b3len],b[b3len:]
  b4,b = b[:b4len],b[b4len:]

  delta = to_vector(b0,l)
  c = to_columnselection(b1,params)
  g = to_poly(b2,params)
  alpha = to_fieldordering(b3,params)
  s = to_vector(b4,n)

  return delta,c,g,alpha,s

# ----- miscellaneous tests

def test_vector():
  assert from_vector([0,0,1,0]) == b'\4'
  assert from_vector([0,0,0,1,0,0,1,0, 1,0,1,0,0,1,1,0]) == b'He'

  assert to_vector(b'He',16) == [0,0,0,1,0,0,1,0, 1,0,1,0,0,1,1,0]
  assert to_vector(b'He',15) == [0,0,0,1,0,0,1,0, 1,0,1,0,0,1,1]
  assert to_vector(b'He',14,narrowly=False) == [0,0,0,1,0,0,1,0, 1,0,1,0,0,1]
  ok = False
  try:
    assert to_vector(b'He',14) == [0,0,0,1,0,0,1,0, 1,0,1,0,0,1]
  except:
    ok = True
  assert ok

  for r in range(100):
    print('bitvector %d' % r)
    sys.stdout.flush()
    v = [randrange(2) for i in range(r)]
    b = from_vector(v)
    assert len(b) == ceil(r/8)
    assert v == to_vector(b,r)
    assert v == to_vector(b,r,narrowly=False)

if __name__ == '__main__':
  test_vector()
