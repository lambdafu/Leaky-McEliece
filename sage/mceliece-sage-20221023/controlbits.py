#!/usr/bin/env python3

# copied from https://eprint.iacr.org/2020/1493
# see proofs there

def permutation(c):
  m = 1
  while (2*m-1)<<(m-1) < len(c): m += 1
  assert (2*m-1)<<(m-1) == len(c)

  n = 1<<m
  pi = list(range(n))
  for i in range(2*m-1):
    gap = 1<<min(i,2*m-2-i)
    for j in range(n//2):
      if c[i*n//2+j]:
        pos = (j%gap)+2*gap*(j//gap)
        pi[pos],pi[pos+gap] = pi[pos+gap],pi[pos]

  return pi

def composeinv(c,pi):
  return [y for x,y in sorted(zip(pi,c))]

def controlbits(pi):
  n = len(pi)
  m = 1
  while 1<<m < n: m += 1
  assert 1<<m == n

  if m == 1: return [pi[0]]
  p = [pi[x^1] for x in range(n)]
  q = [pi[x]^1 for x in range(n)]

  piinv = composeinv(range(n),pi)
  p,q = composeinv(p,q),composeinv(q,p)

  c = [min(x,p[x]) for x in range(n)]
  p,q = composeinv(p,q),composeinv(q,p)
  for i in range(1,m-1):
    cp,p,q = composeinv(c,q),composeinv(p,q),composeinv(q,p)
    c = [min(c[x],cp[x]) for x in range(n)]

  f = [c[2*j]%2 for j in range(n//2)]
  F = [x^f[x//2] for x in range(n)]
  Fpi = composeinv(F,piinv)
  l = [Fpi[2*k]%2 for k in range(n//2)]
  L = [y^l[y//2] for y in range(n)]
  M = composeinv(Fpi,L)
  subM = [[M[2*j+e]//2 for j in range(n//2)] for e in range(2)]
  subz = map(controlbits,subM)
  z = [s for s0s1 in zip(*subz) for s in s0s1]
  return f+z+l

# ----- miscellaneous tests

import sys

def test_onepermutation(pi):
  pi = list(pi)
  n = len(pi)
  if n < 2:
    raise Exception('testing only permutations of length at least 2')
  m = 1
  while 1<<m < n: m += 1
  if 1<<m != n:
    raise Exception('testing only permutations of power-of-2 length')
  assert sorted(pi) == list(range(n))

  c = controlbits(pi)
  assert pi == permutation(c)

def test_small():
  def doit(pi,pos):
    if pos == 1:
      test_onepermutation(pi)
    else:
      for i in range(pos):
        doit(pi,pos-1)
        if pos%2:
          pi[0],pi[pos-1] = pi[pos-1],pi[0]
        else:
          pi[i],pi[pos-1] = pi[pos-1],pi[i]

  for m in range(2,4):
    n = 1<<m
    print('controlbits scan %d %d' % (m,n))
    sys.stdout.flush()
    pi = list(range(n))
    doit(pi,n)

import random

def test_random():
  for m in range(1,14):
    n = 2**m
    print('controlbits random %d %d' % (m,n))
    sys.stdout.flush()
    for loop in range(10):
      r = [random.randrange(2**64)*n+j for j in range(n)]
      r.sort()
      pi = [x%n for x in r]
      test_onepermutation(pi)

if __name__ == '__main__':
  test_small()
  test_random()
