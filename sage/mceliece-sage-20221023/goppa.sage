# goppa_errors() copied from https://eprint.iacr.org/2022/473 Algorithm 6.2
# see Theorems 6.4 and 6.5 for proof

from interpolator import interpolator
from approximant import approximant

def goppa_errors(n,t,k,alpha,g,r):
  r'''
  Return the unique e in F_2^n
  with sum_i (r[i]-e[i]) A/(x-alpha[i]) in gk[x]
  and wt e <= t,
  or None if no such e exists.
  Here A = prod_i (x-alpha[i]).

  INPUTS:

  * "n" - a nonnegative integer

  * "t" - a nonnegative integer

  * "k" - a field containing F_2

  * "alpha" - a list of n distinct elements of k

  * "g" - a squarefree element of the polynomial ring k[x] with deg g = t and g(alpha[j]) nonzero for each j

  * "r" - a list of n elements of k
  '''

  alpha,r = list(alpha),list(r)
  assert k.is_field() and k.characteristic() == 2
  assert g.base_ring() == k and g.degree() == t and g.is_squarefree()
  assert len(alpha) == n and len(set(alpha)) == n and len(r) == n
  kpoly = g.parent()
  A = kpoly(prod(kpoly([-alpha[j],1]) for j in range(n)))
  Aprime = A.derivative()
  rtwist = [r[i]*Aprime(alpha[i])/g(alpha[i])^2 for i in range(n)]
  B = interpolator(n,k,alpha,rtwist)
  a,b = approximant(t,k,A,B)
  aprime = a.derivative()
  if a.divides(A):
    if a.divides(g^2*b-aprime):
      if a*B-b*A == 0 or (a*B-b*A).degree() < n-2*t+a.degree():
        return [k(a(alpha[j]) == 0) for j in range(n)]

# ---- miscellaneous tests
# copied from https://eprint.iacr.org/2022/473 Figure A.4

def test_smallrandom():
  for m in range(1,10):
    q = 2^m
    print('goppa_errors %d' % q)
    sys.stdout.flush()
    k = GF(q)
    kpoly.<x> = k[]
    for loop in range(100):
      while True:
        n = randrange(q+1)
        t = randrange(3+n//m)
        if t >= n: t = n
        a = list(k)
        shuffle(a)
        a = a[:n]
        g = kpoly([k.random_element() for j in range(t)]+[1])
        if g.is_squarefree(): 
          if all(g(aj) != 0 for aj in a):
            break
  
      assert g.degree() == t
      A = kpoly(prod(x-aj for aj in a))
      Aprime = A.derivative()
      for aj in a: assert Aprime(aj) != 0
  
      for known in True,False:
        if known:
          f = kpoly([k.random_element() for j in range(n-2*t)])
          r = [(f*g^2)(aj)/Aprime(aj) for aj in a]
          if randrange(2):
            e = [1]*t+[0]*(n-t)
          else:
            actualweight = randrange(t+1)
            e = [1]*actualweight+[0]*(n-actualweight)
          shuffle(e)
          assert len([ej for ej in e if ej != 0]) <= t
          for j in range(n): r[j] += e[j]
        else:
          e = 'unknown' # cut off data flow from previous iteration
          r = [k.random_element() for j in range(n)]
        e2 = goppa_errors(n,t,k,a,g,r)
        if e2 == None:
          assert not known
        else:
          assert len(e2) == n
          if known: assert e2 == e
          assert len([ej for ej in e2 if ej != 0]) <= t
          assert g.divides(sum((r[i]-e2[i])*A//(x-a[i]) for i in range(n)))

if __name__ == '__main__':
  test_smallrandom()
