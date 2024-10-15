# approximant() copied from https://eprint.iacr.org/2022/473 Algorithm 4.4
# see Theorem 4.1 for proof

# uses basic matrices and polys, not Sage's berlekamp_massey()

def approximant(t,k,A,B):
  r'''
  Return a,b in the polynomial ring k[x]
  with gcd{a,b} = 1,
  deg a <= t,
  deg b < t,
  and deg(aB-bA) < deg A - t.

  INPUT:

  "t" - a nonnegative integer

  "k" - a field
  
  "A" - an element of k[x]

  "B" - an element of k[x] with deg A > deg B
  '''

  assert t >= 0 and A.base_ring() == k and B.base_ring() == k
  kpoly,n = A.parent(),A.degree()
  assert n > B.degree()
  M = [   [ B[t+n-1-i-j] for i in range(t+1)]
        + [-A[t+n-1-i-j] for i in range(t)  ] for j in range(2*t)]
  M = matrix(k,2*t,2*t+1,M)
  ab = list(M.right_kernel().gens()[0])
  a,b = kpoly(ab[:t+1]),kpoly(ab[t+1:])
  d = gcd(a,b)
  return a//d,b//d

# ---- miscellaneous tests
# copied from https://eprint.iacr.org/2022/473 Figure A.2

def test_smallrandom():
  for q in range(100):
    q = ZZ(q)
    if not q.is_prime_power(): continue
    print('approximant %d' % q)
    sys.stdout.flush()
    k = GF(q)
    kpoly.<x> = k[]
    for loop in range(100):
      Adeg = randrange(100)
      A = kpoly([k.random_element() for j in range(Adeg)]+[1])
      if Adeg == 0:
        B = kpoly(0)
      else:
        Bdeg = randrange(Adeg)
        B = kpoly([k.random_element() for j in range(Bdeg+1)])
        # note that B could actually have lower degree
      t = randrange(Adeg+3)
      a,b = approximant(t,k,A,B)
      assert gcd(a,b) == 1
      assert a.degree() <= t
      assert b.degree() < t
      assert a != 0
      assert a*B-b*A == 0 or (a*B-b*A).degree() < A.degree()-t

if __name__ == '__main__':
  test_smallrandom()
