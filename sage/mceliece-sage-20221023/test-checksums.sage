import parameters
import encap
import decap
import keygen

def randomuint8():
  global randombytes_key
  global randombytes_buf
  global randombytes_pos
  if randombytes_pos == len(randombytes_buf):
    h = []

    def quarterround(a,b,c,d):
      a += b
      a &= 0xffffffff
      d = ZZ(d).__xor__(a)
      d = (d<<16)|(d>>16)
      d &= 0xffffffff
      c += d
      c &= 0xffffffff
      b = ZZ(b).__xor__(c)
      b = (b<<12)|(b>>20)
      b &= 0xffffffff
      a += b
      a &= 0xffffffff
      d = ZZ(d).__xor__(a)
      d = (d<<8)|(d>>24)
      d &= 0xffffffff
      c += d
      c &= 0xffffffff
      b = ZZ(b).__xor__(c)
      b = (b<<7)|(b>>25)
      b &= 0xffffffff
      return a,b,c,d

    for pos in range(12):
      x = [1634760805,857760878,2036477234,1797285236]
      k = randombytes_key
      for i in range(0,32,4):
        x += [k[i]+(k[i+1]<<8)+(k[i+2]<<16)+(k[i+3]<<24)]
      x += [pos,0,0,0]
      y = copy(x)
      for i in range(10):
        x[0],x[4],x[8],x[12] = quarterround(x[0],x[4],x[8],x[12])
        x[1],x[5],x[9],x[13] = quarterround(x[1],x[5],x[9],x[13])
        x[2],x[6],x[10],x[14] = quarterround(x[2],x[6],x[10],x[14])
        x[3],x[7],x[11],x[15] = quarterround(x[3],x[7],x[11],x[15])
        x[0],x[5],x[10],x[15] = quarterround(x[0],x[5],x[10],x[15])
        x[1],x[6],x[11],x[12] = quarterround(x[1],x[6],x[11],x[12])
        x[2],x[7],x[8],x[13] = quarterround(x[2],x[7],x[8],x[13])
        x[3],x[4],x[9],x[14] = quarterround(x[3],x[4],x[9],x[14])
      for i in range(16):
        z = x[i]+y[i]
        z &= 0xffffffff
        h += [z&255]; z >>= 8
        h += [z&255]; z >>= 8
        h += [z&255]; z >>= 8
        h += [z&255]; z >>= 8

    randombytes_key = h[:32]
    randombytes_buf = h[32:]
    randombytes_pos = 0
  c = randombytes_buf[randombytes_pos]
  randombytes_buf[randombytes_pos] = None
  randombytes_pos += 1
  return c

def randombytes(r):
  return bytes(bytearray([randomuint8() for j in range(r)]))

def L32(x,c):
  x &= 0xffffffff
  result = (x << c) | (x >> (32 - c))
  result &= 0xffffffff
  return result

def ld32(x):
  assert len(x) == 4
  u = x[3]
  u = (u<<8)|x[2]
  u = (u<<8)|x[1]
  return (u<<8)|x[0]

def st32(x):
  result = []
  for i in range(4):
    result += [255&x]
    x >>= 8
  return result

def core(block,k):
  assert len(block) == 16

  x = [None]*16

  sigma = b'expand 32-byte k'
  sigma = list(bytearray(sigma))

  for i in range(4):
    x[5*i] = ld32(sigma[4*i:4*i+4])
    x[1+i] = ld32(k[4*i:4*i+4])
    x[6+i] = ld32(block[4*i:4*i+4])
    x[11+i] = ld32(k[4*i+16:4*i+20])

  y = copy(x)

  for i in range(20):
    w = [None]*16
    for j in range(4):
      t = [None]*4
      for m in range(4):
        t[m] = ZZ(x[(5*j+4*m)%16])
      t[1] = t[1].__xor__(L32(t[0]+t[3], 7));
      t[2] = t[2].__xor__(L32(t[1]+t[0], 9));
      t[3] = t[3].__xor__(L32(t[2]+t[1],13));
      t[0] = t[0].__xor__(L32(t[3]+t[2],18));
      for m in range(4):
        w[4*j+((j+m)%4)] = t[m]
    for m in range(16):
      x[m] = w[m]

  result = []

  for i in range(16):
    z = 0xffffffff & (x[i]+y[i])
    result += st32(z)

  return result

def checksum(x):
  global checksum_state

  if type(x) == type(b'123'):
    x = list(bytearray(x))

  info = 'checksum %s' % ''.join('%02x'%xi for xi in x)

  while len(x) >= 16:
    checksum_state = core(x[:16],checksum_state)
    x = x[16:]

  info += ' %s' % ''.join('%02x'%ci for ci in checksum_state)

  block = copy(x) + [1] + [0]*(15-len(x))
  checksum_state[0] = checksum_state[0].__xor__(1)
  checksum_state = core(block,checksum_state)

  info += ' %s' % ''.join('%02x'%ci for ci in checksum_state)
  # print(info)
  
def salsa20(outlen,n,k):
  assert outlen >= 0
  assert len(n) == 8
  assert len(k) == 32

  result = []

  if outlen == 0: return result

  z = n + [0]*8

  while outlen >= 64:
    result += core(z,k)
    outlen -= 64

    for i in range(8,16):
      z[i] = 255&(z[i]+1)
      if z[i]: break

  if outlen > 0:
    result += core(z,k)[:outlen]

  return result

def testvector(outlen):
  k = b'generate inputs for test vectors'
  k = list(bytearray(k))
  result = salsa20(outlen,testvector_n,k)
  for i in range(8):
    testvector_n[i] = 255&(testvector_n[i]+1)
    if testvector_n[i]: break
  return result

def myrandom():
  x = testvector(8)
  return sum(x[i]<<(8*i) for i in range(8))

systems = parameters.alltests
if len(sys.argv) > 1:
  systems = sys.argv[1:]

for system in systems:
  randombytes_key = [0]*32
  randombytes_buf = [0]*736
  randombytes_pos = len(randombytes_buf)

  checksum_state = [0]*64

  testvector_n = [0]*8

  params = parameters.parameters(system,allowtestparams=True)

  result = system
  print(result)
  sys.stdout.flush()

  for loop in range(64):
    print(loop)
    sys.stdout.flush()

    pk,sk = keygen.keygen(randombytes,params)
    checksum(pk)
    checksum(sk)

    C,k = encap.encap(pk,randombytes,params)
    checksum(C)
    checksum(k)

    assert decap.decap(C,sk,params) == k
    checksum(k)
    
    for loop2 in range(3):
      Clen = len(C)
      C = list(bytearray(C))

      offset = 1 + (myrandom() % 255)
      pos = myrandom() % Clen
      C[pos] = 255&(C[pos]+offset)

      C = bytes(bytearray(C))

      k2 = decap.decap(C,sk,params)

      if k2 == False:
        checksum(C)
      else:
        checksum(k2)

    if loop in [7,63]:
      checksumhex = ''
      for i in range(32):
        checksumhex += '%x' % (15&(checksum_state[i]>>4))
        checksumhex += '%x' % (15&checksum_state[i])

      result += ' ' + checksumhex
      print(result)
      sys.stdout.flush()

  print(result)
  sys.stdout.flush()
