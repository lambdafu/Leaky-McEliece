import sys
load("support_splitting.sage")

params = [
  (512,20,9),
  (1024,50,10),
  (3488,64,12),
  (4608,96,13),
  (6688,128,13),
  (6960,119,13),
  (8192,128,13)
]

set_random_seed(0)

def main(n,k):

  #Construct random k x n matrix in systematic form G.  
  G = zero_matrix(GF(2),k,n)
  while G.submatrix(0,0,k,k).rank() < k:
      G = random_matrix(GF(2),k,n).echelon_form()

  G_perm = copy(G)

  #Apply random permutation from right to G_perm.
  G_perm = G_perm.transpose()
  ctr = 0
  while ctr < n or G_perm.submatrix(0,0,k,k).rank() < k:
      i = randint(0,n-1)
      j = randint(0,n-1)
      G_perm.swap_rows(i,j)
      ctr += 1
  G_perm = G_perm.transpose()

  #Apply random invertible matrix from left to G_perm.
  for _ in range(n):
      i = randint(0,k-1)
      j = randint(0,k-1)
      if i != j:
        G_perm.add_multiple_of_row(i,j,1)
  G_perm = G_perm.echelon_form()

  # Run support splitting on G and G_perm
  tt = cputime()
  perm = supportSplitting(G,G_perm)
  time = cputime(tt)
  
  # Verify correctness of solution
  # Can be pretty slow.
  assert (G*perm).image() == G_perm.image()
  
  return time

for n,t,m in params:
  k = n - t*m
  
  print(n,t,m)
  
  tMin = Infinity
  tMax = 0
  tAvg = 0
  
  for _ in range(10):
    time = main(n,k)
    
    if time < tMin:
      tMin = time
    if time > tMax:
      tMax = time
      
    tAvg += time / 10
    
    
  print("tAvg: %f\ttMin: %f\ttMax: %f" % (tAvg,tMin,tMax))
  print("")