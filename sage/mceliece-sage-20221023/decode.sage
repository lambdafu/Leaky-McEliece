import goppa
import parameters

def decode(C,Gammaprime,params):
  r'''
  Return the output of the Classic McEliece Decode function
  on input (C,Gammaprime) for the specified parameters.

  The output is either False or a length-n vector over F_2.

  INPUT:

  "Gammaprime" - a tuple (g,alphaprime[0],...,alphaprime[n-1])
  where g is a monic irreducible polynomial in F_q[x] of degree t
  and alphaprime[0],...,alphaprime[n-1] are distinct elements of F_q

  "params" - a Classic McEliece parameter set
  '''
  assert isinstance(params,parameters.parameters)
  m = params.m
  n = params.n
  t = params.t
  k = params.k
  Fq = params.Fq

  Gammaprime = tuple(Gammaprime)
  assert len(Gammaprime) == n+1
  g,alphaprime = Gammaprime[0],Gammaprime[1:]
  assert all(alphaj in Fq for alphaj in alphaprime)
  assert len(set(alphaprime)) == n
  assert g.base_ring() == Fq
  assert g.is_monic()
  assert g.is_irreducible()
  assert g.degree() == t

  assert len(C) == m*t

  # "Extend C to v = (C,0,...,0) in F_2^n by appending k zeros."
  v = list(C)+[0]*k
  v = list(vector(GF(2),v))

  # "Find the unique c in F_2^n such that (1) Hc = 0 and
  #  (2) c has Hamming distance <=t from v.
  #  If there is no such c, return False. ...
  #  Set e = v+c."
  e = goppa.goppa_errors(n,t,Fq,alphaprime,g,v)
  if e is None: return False
  c = vector(GF(2),[e[i]+v[i] for i in range(n)])

  # "If wt(e) = t and C = He, return e. Otherwise return False."
  if sum(1 if ej else 0 for ej in e) != t: return False

  kpoly = g.parent()
  xminusalpha = [kpoly([-alphaprime[i],1]) for i in range(n)]
  A = kpoly(prod(xminusalpha))
  if sum(c[i]*A//xminusalpha[i] for i in range(n))%g != 0: return False

  return e

def test1():
  import matgen
  import encode

  for system in parameters.alltests:
    P = parameters.parameters(system,allowtestparams=True)
    Fq = P.Fq
    Fqx.<x> = Fq[]
    m = P.m
    n = P.n
    t = P.t
    k = P.k
    mu = P.mu
    nu = P.nu

    tries = 0
    while True:
      print('matgen %s' % system)
      sys.stdout.flush()
      tries += 1
      while True:
        g = Fqx.random_element(t)
        if g.is_irreducible(): break
      g /= g.leading_coefficient()
      alpha = list(Fq)
      shuffle(alpha)
      alpha = alpha[:n]
      Gamma = tuple([g]+alpha)
      result = matgen.matgen(Gamma,P)

      if result == False: continue

      T = result[0]
      Gammaprime = result[-1]
      break

    print('successful matgen; tries: %d' % tries)
    sys.stdout.flush()

    # "T is an mt x k matrix over F_2"
    assert T.nrows() == m*t
    assert T.ncols() == k

    # "Gamma' has the form (g,alpha'_0,alpha'_1,...,alpha'_{n-1})"
    g,alphaprime = Gamma[0],Gamma[1:]
    assert len(alphaprime) == n

    # "g is a monic irreducible polynomial in F_q[x] of degree t"
    assert g.base_ring() == Fq
    assert g.is_monic()
    assert g.is_irreducible()
    assert g.degree() == t

    # "alpha'_0,alpha'_1,...,alpha'_{n-1} are distinct elements of F_q"
    assert all(alphaj in Fq for alphaj in alphaprime)
    assert len(set(alphaprime)) == n

    # ----- decoding a legitimate ciphertext

    e = [0]*n
    while sum(e) < t: e[randrange(n)] = 1
    C = encode.encode(e,T,P)

    # "If C = Encode(e,T) then Decode(C,Gammaprime) = e."
    assert list(decode(C,Gammaprime,P)) == list(e)
    print('successful decode')
    sys.stdout.flush()

    # ----- decoding a ciphertext for weight t-1

    if t > 0:
      while True:
        j = randrange(n)
        if not e[j]: continue
        edelta = vector(GF(2),[i==j for i in range(n)])
        Cmod = list(C)
        for i in range(m*t):
          Cmod[i] += edelta[i]+sum(T[i,j]*edelta[j+m*t] for j in range(k))
        break
      
      assert decode(Cmod,Gammaprime,P) == False
      print('successful anti-decode for reduced weight')
      sys.stdout.flush()

    # ----- decoding a random ciphertext

    # "If C does not have the form He for any weight-t vector e in F_2^n,
    #  then Decode(G,Gamma') = False."
    C = vector(GF(2),[randrange(2) for j in range(m*t)])
    assert decode(C,Gammaprime,P) == False
    print('successful anti-decode for random ciphertext')
    sys.stdout.flush()

if __name__ == '__main__':
  test1()
