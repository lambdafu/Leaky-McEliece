import parameters

def fixedweight(randombits,params):
  r'''
  Return the output of the Classic McEliece FixedWeight function
  using the specified source of random bits
  for the specified parameters.

  INPUT:

  "randombits" - a function that, on input r, returns a list of r bits

  "params" - a Classic McEliece parameter set
  '''

  assert isinstance(params,parameters.parameters)
  m = params.m
  q = params.q
  n = params.n
  t = params.t
  sigma1 = params.sigma1

  # "The integer tau is defined as t if n=q;
  #  as 2t if q/2 <= n < q; as 4t if q/4 <= n < q/2; etc."
  tau = t
  while tau*n < t*q: tau *= 2
  
  # "All of the selected parameter sets have q/2 <= n <= q,
  #  so tau in {t,2t}."
  assert q/2 <= n
  assert n <= q
  assert tau == t or tau == 2*t

  while True:
    # "Generate sigma_1 tau uniform random bits b_0,b_1,...,b_{sigma_1 tau - 1}."
    b = randombits(sigma1*tau)

    # "Define d_j = sum_{i=0}^{m-1} b_{sigma_1 j + i} 2^i for each j in {0,1,...,tau-1}."
    d = [sum(b[sigma1*j+i]<<i for i in range(m)) for j in range(tau)]

    # "Define a_0,a_1,...,a_{t-1} as the first t entries in d_0,d_1,...,d_{tau-1}
    #  in the range {0,1,...,n-1}.
    #  If there are fewer than t such entries, restart the algorithm."
    a = [dj for dj in d if dj<n]
    if len(a) < t: continue
    a = a[:t]

    # "If a_0,a_1,...,a_{t-1} are not all distinct, restart the algorithm."
    if len(set(a)) < t: continue

    # "Define e = (e_0,e_1,...,e_{n-1}) in F_2^n
    #  as the weight-t vector such that e_{a_i} = 1 for each i."
    e = [0]*n
    for ai in a: e[ai] = 1

    # "Return e."
    return e

# ----- miscellaneous tests

def test1():
  def randombits(r):
    return [randrange(2) for j in range(r)]

  for system in parameters.alltests:
    P = parameters.parameters(system,allowtestparams=True)
    n = P.n
    t = P.t

    print('fixedweight %s' % system)
    sys.stdout.flush()

    for loop in range(10):
      e = fixedweight(randombits,P)
      
      "outputs a vector e in F_2^n of weight t"
      assert len(e) == n
      assert sum(1 if ej else 0 for ej in e) == t

if __name__ == '__main__':
  test1()
