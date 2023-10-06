import parameters

def encode(e,T,params):
  r'''
  Return the output of the Classic McEliece Encode function
  on input (e,T) for the specified parameters.

  INPUT:

  "e" - a length-n weight-t vector over F_2

  "T" - an mt x n matrix over F_2

  "params" - a Classic McEliece parameter set
  '''

  assert isinstance(params,parameters.parameters)
  m = params.m
  n = params.n
  t = params.t
  k = params.k

  e = vector(GF(2),list(e))
  assert len(e) == n
  assert sum(1 for ej in e if ej != 0) == t

  assert T.nrows() == m*t
  assert T.ncols() == k
  assert T.base_ring() == GF(2)

  # "Define H = (I_{mt} | T)."
  H = identity_matrix(GF(2),m*t).augment(T)

  # "Compute and return C = He in F_2^{mt}."
  C = H*e
  assert len(C) == m*t
  return C

# ----- miscellaneous tests

def test1():
  for system in parameters.alltests:
    P = parameters.parameters(system,allowtestparams=True)
    m = P.m
    t = P.t
    n = P.n
    k = P.k

    print('encode %s' % system)
    sys.stdout.flush()

    for loop in range(10):
      T = matrix(GF(2),[[randrange(2) for j in range(k)] for i in range(m*t)])
      e = [0]*n
      while sum(e) < t: e[randrange(n)] = 1
  
      C = encode(e,T,P)
  
      assert len(C) == m*t
      assert C.base_ring() == GF(2)
      for i in range(m*t):
        assert C[i] == e[i]+sum(T[i,j]*e[j+m*t] for j in range(k))

if __name__ == '__main__':
  test1()
