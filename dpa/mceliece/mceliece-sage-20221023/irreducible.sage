import parameters

def irreducible(d,params):
  r'''
  Return the output of the Classic McEliece Irreducible function
  on input d for the specified parameters.

  INPUT:

  "d" - a list of sigma1*t bits

  "params" - a Classic McEliece parameter set
  '''

  assert isinstance(params,parameters.parameters)
  m = params.m
  t = params.t
  Fq = params.Fq
  z = Fq.gen()
  Fqt = params.Fqt
  y = Fqt.gen()
  sigma1 = params.sigma1

  d = list(d)
  assert len(d) == sigma1*t

  # "Define beta_j = sum_{i=0}^{m-1} d_{sigma_1 j+i} z^i 
  #  for each j in {0,1,...,t-1}."
  beta = [sum(d[sigma1*j+i]*z^i for i in range(m)) for j in range(t)]

  # "Define beta = beta_0 + beta_1 y + ... + beta_{t-1} y^{t-1} in F_q[y]/F(y)."
  beta = Fqt(sum(beta[j]*y^j for j in range(t)))

  # "Compute the minimal polynomial g of beta over F_q. ...
  #  Return g if g has degree t. Otherwise return False."
  betapow = [(beta^i).list() for i in range(t+1)]
  M = matrix(Fq,[[betapow[i][j] for j in range(t)] for i in range(t)])
  v = vector(Fq,[betapow[t][j] for j in range(t)])
  if not M.is_invertible(): return False
  gcoeffs = v*M.inverse()

  Fqx.<x> = Fq[]

  g = Fqx(list(gcoeffs)+[1])
  assert g(beta) == 0
  assert g.degree() == t
  return g

# ----- miscellaneous tests

def test1():
  for system in parameters.alltests:
    P = parameters.parameters(system,allowtestparams=True)
    Fq = P.Fq
    Fqt = P.Fqt
    Fqx.<x> = Fq[]
    m = P.m
    t = P.t
    sigma1 = P.sigma1

    numsuccess = 0
    for loop in range(10):
      d = [randrange(2) for j in range(sigma1*t)]
      g = irreducible(d,P)
      if g != False:
        numsuccess += 1
        assert g.is_monic()
        assert g.is_irreducible()
        assert g.degree() == t

    print('irreducible %s numsuccess %d' % (system,numsuccess))
    sys.stdout.flush()

if __name__ == '__main__':
  test1()
