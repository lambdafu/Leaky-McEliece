F2 = GF(2)
F2z.<z> = F2[]
F2zy.<y> = F2z[]

f12 = z^12+z^3+1
f13 = z^13+z^4+z^3+z+1
F64 = y^64+y^3+y+z
F96 = y^96+y^10+y^9+y^6+1
F119 = y^119+y^8+1
F128 = y^128+y^7+y^2+y+1

selected = {
  # "mceliece348864": "m = 12, n = 3488, t = 64. Field polynomials f(z) = z^12+z^3+1 and F(y) = y^64+y^3+y+z."
  'mceliece348864':{'m':12,'n':3488,'t':64,'f':f12,'F':F64},
  # "mceliece348864f": "m = 12, n = 3488, t = 64. Field polynomials f(z) = z^12+z^3+1 and F(y) = y^64+y^3+y+z. Semi-systematic parameters (mu,nu) = (32,64)."
  'mceliece348864f':{'m':12,'n':3488,'t':64,'f':f12,'F':F64,'mu':32,'nu':64},
  # "mceliece460896": "m = 13, n = 4608, t = 96. Field polynomials f(z) = z^13+z^4+z^3+z+1 and F(y) = y^96+y^10+y^9+y^6+1."
  'mceliece460896':{'m':13,'n':4608,'t':96,'f':f13,'F':F96},
  # "mceliece460896f": "m = 13, n = 4608, t = 96. Field polynomials f(z) = z^13+z^4+z^3+z+1 and F(y) = y^96+y^10+y^9+y^6+1. Semi-systematic parameters (mu,nu) = (32,64)."
  'mceliece460896f':{'m':13,'n':4608,'t':96,'f':f13,'F':F96,'mu':32,'nu':64},
  # "mceliece6688128": "m = 13, n = 6688, t = 128. Field polynomials f(z) = z^13+z^4+z^3+z+1 and F(y) = y^128+y^7+y^2+y+1."
  'mceliece6688128':{'m':13,'n':6688,'t':128,'f':f13,'F':F128},
  # "mceliece6688128f": "m = 13, n = 6688, t = 128. Field polynomials f(z) = z^13+z^4+z^3+z+1 and F(y) = y^128+y^7+y^2+y+1. Semi-systematic parameters (mu,nu) = (32,64)."
  'mceliece6688128f':{'m':13,'n':6688,'t':128,'f':f13,'F':F128,'mu':32,'nu':64},
  # "mceliece6960119": "m = 13, n = 6960, t = 119. Field polynomials f(z) = z^13+z^4+z^3+z+1 and F(y) = y^119+y^8+1."
  'mceliece6960119':{'m':13,'n':6960,'t':119,'f':f13,'F':F119},
  # "mceliece6960119f": "m = 13, n = 6960, t = 119. Field polynomials f(z) = z^13+z^4+z^3+z+1 and F(y) = y^119+y^8+1. Semi-systematic parameters (mu,nu) = (32,64)."
  'mceliece6960119f':{'m':13,'n':6960,'t':119,'f':f13,'F':F119,'mu':32,'nu':64},
  # "mceliece8192128": "m = 13, n = 8192, t = 128. Field polynomials f(z) = z^13+z^4+z^3+z+1 and F(y) = y^128+y^7+y^2+y+1."
  'mceliece8192128':{'m':13,'n':8192,'t':128,'f':f13,'F':F128},
  # "mceliece8192128f": "m = 13, n = 8192, t = 128. Field polynomials f(z) = z^13+z^4+z^3+z+1 and F(y) = y^128+y^7+y^2+y+1. Semi-systematic parameters (mu,nu) = (32,64)."
  'mceliece8192128f':{'m':13,'n':8192,'t':128,'f':f13,'F':F128,'mu':32,'nu':64},
}

f9 = z^9+z^1+1 # for testing
f10 = z^10+z^3+1 # for testing
F20 = z^20+z^3+1 # for testing
F50 = y^50+y^5+y^2+z # for testing

testing = {
  'mceliece51220f':{'m':9,'n':512,'t':20,'f':f9,'F':F20,'mu':32,'nu':64},
  'mceliece51220':{'m':9,'n':512,'t':20,'f':f9,'F':F20},
  'mceliece102450f':{'m':10,'n':1024,'t':50,'f':f10,'F':F50,'mu':32,'nu':64},
  'mceliece102450':{'m':10,'n':1024,'t':50,'f':f10,'F':F50},
}

supportpc = True # for testing against previous versions
if supportpc:
  for which in selected,testing:
    for system in list(which):
      systempc = system+'pc'
      systempc = systempc.replace('fpc','pcf')
      assert systempc not in which
      which[systempc] = dict(which[system])
      which[systempc]['pc'] = True

alltests = list(testing)+list(selected)

import byterepr
import hashlib

class parameters:
  def __init__(self,system,checkirred=False,allowtestparams=False):
    r'''
    Return one of the selected Classic McEliece parameter sets,
    specifically the one named by "system".

    INPUT:

    "system" - a string such as 'mceliece6960119'

    "checkirred" (optional, default False) -
    check irreducibility of the f and F polynomials

    "allowtestparams" (optional, default False) -
    allow test parameter sets, not just selected parameter sets
    '''

    if allowtestparams and system in testing:
      S = testing[system]
    else:
      assert system in selected
      S = selected[system]

    m = S['m']
    n = S['n']
    t = S['t']
    f = S['f']
    F = S['F']

    if 'mu' in S and 'nu' in S:
      mu = S['mu']
      nu = S['nu']
    else:
      # "parameter sets that do not mention these parameters define them as (0,0) by default"
      mu = 0
      nu = 0

    q = 1<<m
    k = n-m*t

    if supportpc:
      self.pc = S.get('pc',False)

    # "positive integer m"
    assert parent(m) is ZZ
    assert m > 0
    self.m = m

    # "also defines a parameter q = 2^m"
    assert parent(q) is ZZ
    assert q == 2^m
    self.q = q

    # "positive integer n with n <= q"
    assert parent(n) is ZZ
    assert n > 0
    assert n <= q
    self.n = n

    # "positive integer t >= 2 with mt < n"
    assert parent(t) is ZZ
    assert t >= 2
    assert m*t < n
    self.t = t

    # "also defines a parameter k = n-mt"
    assert parent(k) is ZZ
    assert k == n-m*t
    self.k = k

    # "monic irreducible polynomial f(z) in F_2[z] of degree m"
    assert parent(f) is F2z
    assert f.is_monic()
    if checkirred: assert f.is_irreducible()
    assert f.degree() == m
    self.f = f

    # "defines a representation F2[z]/f(z) of the field Fq"
    Fq.<zinFq> = GF(q,name='zinFq',modulus=f,check_irreducible=checkirred)
    assert Fq.is_field()
    assert Fq.cardinality() == q
    self.Fq = Fq

    F = F.change_ring(Fq) # now F is in Fq[y] instead of F2[z][y]

    # "monic irreducible polynomial F(y) in F_q[y] of degree t"
    assert F.base_ring() == Fq
    assert F.is_monic()
    if checkirred: assert F.is_irreducible()
    assert F.degree() == t
    self.F = F

    # "defines a representation F_q[y]/F(y) of the field F_{q^t} = F_{2^{mt}}"
    Fqt.<yinFqt> = Fq.extension(F) # no easy way to pass checkirred to this
    assert Fqt.is_field()
    assert Fqt.cardinality() == q^t
    self.Fqt = Fqt

    # "integers nu >= mu >= 0 with nu <= k+mu"
    assert parent(mu) is ZZ
    assert parent(nu) is ZZ
    assert nu >= mu
    assert mu >= 0
    assert nu <= k+mu
    self.mu = mu
    self.nu = nu

    # ----- symmetric-cryptography parameters

    # "positive integer l"
    # selected symmetric-cryptography parameters:
    # "The integer l is 256"
    l = 256
    self.l = l

    # "integer sigma1 >= m"
    # selected symmetric-cryptography parameters:
    # "The integer sigma1 is 16"
    sigma1 = 16
    assert sigma1 >= m
    self.sigma1 = sigma1

    # "integer sigma2 >= 2m"
    # selected symmetric-cryptography parameters:
    # "The integer sigma2 is 32"
    sigma2 = 32
    assert sigma2 >= 2*m
    self.sigma2 = sigma2

    # extra parameter constraints for the H,G implementations below,
    # all satisfied by the selected parameter sets:
    assert n%8 == 0
    assert l%8 == 0
    assert sigma1%8 == 0
    assert sigma2%8 == 0

    def shake256(input,outlen):
      h = hashlib.shake_256()
      h.update(input)
      return h.digest(int(outlen))

    # "cryptographic hash function H that outputs l bits"
    def H(x):
      # selected symmetric-cryptography parameters:
      # "The l-bit string H(x) is defined as the first l bits of output of SHAKE256(x)."
      outlen = l
      assert outlen%8 == 0

      inbytes = byterepr.from_hashinput(x,self)

      if supportpc and self.pc:
        # for testing against previous versions
        assert list(bytearray(inbytes[0:1])) in ([0],[1],[2])
      else:
        # "All H inputs used in Classic McEliece begin with byte 0 or 1"
        assert list(bytearray(inbytes[0:1])) in ([0],[1])

      result = shake256(inbytes,outlen/8)
      return byterepr.to_vector(result,outlen)

    self.H = H

    # "pseudorandom bit generator G mapping a string of l bits
    #   to a string of n + sigma2 q + sigma1 t + l bits"
    def G(delta):
      # selected symmetric-cryptography parameters:
      # "The (n + sigma_2 q + sigma_1 t + l)-bit string G(delta)
      #  is defined as the first n + sigma_2 q + sigma_1 t + l bits
      #  of output of SHAKE256(64,delta).
      #  Here 64,delta means the 33-byte string
      #  that begins with byte 64 and continues with delta."

      inbytes = byterepr.from_vector(delta)
      inbytes = bytes(bytearray([64]))+inbytes
      outlen = n+sigma2*q+sigma1*t+l
      assert outlen%8 == 0

      # G input begins with byte 64
      assert list(bytearray(inbytes[0:1])) == [64]

      result = shake256(inbytes,outlen/8)
      return byterepr.to_vector(result,outlen)

    self.G = G

# ----- miscellaneous tests

def test_assertions():
  for system in selected:
    print('parameters %s' % system)
    sys.stdout.flush()
    P = parameters(system,checkirred=True)
    assert P.m == selected[system]['m']
    assert P.n == selected[system]['n']
    assert P.t == selected[system]['t']

if __name__ == '__main__':
  test_assertions()
