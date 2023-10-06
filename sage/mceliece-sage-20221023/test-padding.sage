import parameters
import keygen
import encap
import decap

def randombytes(r):
  return os.urandom(r)

systems = parameters.alltests
if len(sys.argv) > 1:
  systems = sys.argv[1:]

for system in systems:
  print(system)
  sys.stdout.flush()

  params = parameters.parameters(system,allowtestparams=True)

  pk,sk = keygen.keygen(randombytes,params)
  C,sessionkey = encap.encap(pk,randombytes,params)
  assert decap.decap(C,sk,params) == sessionkey

  k = params.k
  mt = params.m*params.t

  if k < 8*ceil(k/8):
    for row in range(mt):
      rowpos = row*8*ceil(k/8)
      for j in range(rowpos+k,rowpos+8*ceil(k/8)):
        print('pk padding bit',j)
        sys.stdout.flush()
        pk2 = bytearray(pk)
        pk2[j//8] |= 1<<(j%8)
        pk2 = bytes(pk2)
        assert pk2 != pk
        assert encap.encap(pk2,randombytes,params) == False

  for loop in range(100):
    row = randrange(mt)
    rowpos = row*8*ceil(k/8)
    j = randrange(rowpos,rowpos+k)
    print('pk real bit',j)
    sys.stdout.flush()
    pk2 = bytearray(pk)
    pk2[j//8] = ZZ(pk2[j//8]).__xor__(1<<(j%8))
    pk2 = bytes(pk2)
    assert pk2 != pk
    assert encap.encap(pk2,randombytes,params) != False

  if mt < 8*ceil(mt/8):
    for j in range(mt,8*ceil(mt/8)):
      print('C padding bit',j)
      sys.stdout.flush()
      C2 = bytearray(C)
      C2[j//8] |= 1<<(j%8)
      C2 = bytes(C2)
      assert C2 != C
      assert decap.decap(C2,sk,params) == False

  for loop in range(10):
    j = randrange(mt)
    print('C real bit',j)
    sys.stdout.flush()
    C2 = bytearray(C)
    C2[j//8] = ZZ(C2[j//8]).__xor__(1<<(j%8))
    C2 = bytes(C2)
    assert C2 != C
    assert decap.decap(C2,sk,params) not in (False,sessionkey)
