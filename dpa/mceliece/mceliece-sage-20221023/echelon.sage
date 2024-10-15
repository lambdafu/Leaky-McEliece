def echelon(X):
  r'''
  Return X reduced to row-echelon form.

  INPUT:

  "X" - a matrix
  '''
  return copy(X.rref())

def echelon_positions(R):
  r'''
  If R is in reduced row-echelon form,
  return the list c such that
  len(c) is the rank of R and, for each j in {0,1,...,r-1},
  c[j] is the position of the leading 1 in row j of R.

  If R is not in reduced row-echelon form,
  raise an exception.
  
  INPUT:

  "R" - a matrix over a field
  '''

  r = R.rank()
  c = [min(b for b,Rab in enumerate(R[a]) if Rab != 0) for a in range(r)]

  # "c_0 < c_1 < ... < c_{r-1}"
  assert c == sorted(set(c))

  for a in range(r):
    # "row 0 of R begins with a 1 in column c_0"
    # "row 1 of R begins with a 1 in column c_1"
    # etc.
    assert R[a][:c[a]] == 0
    assert R[a,c[a]] == 1

  for a in range(r):
    for b in range(R.nrows()):
      if b != a:
        # "this is the only nonzero entry in column c_0"
        # "the only nonzero entry in column c_1"
        # etc.
        assert R[b,c[a]] == 0

  for a in range(r,R.nrows()):
    # "all subsequent rows of R are 0"
    assert R[a] == 0

  return c

def is_systematic(R):
  r'''
  Return True if R is in systematic form,
  else False.

  INPUT:

  "R" - a matrix in reduced row-echelon form over a field
  '''
  c = echelon_positions(R)
  r = len(c)

  # "R has exactly zero rows, i.e., there are no zero rows"
  if R.nrows() != r: return False

  # "c_i = i for 0 <= i < r"
  okc = all(c[i] == i for i in range(r))

  # "equivalent to simply saying c_{r-1} == r-1,
  #  except in the degenerate case r = 0"
  if r > 0: assert okc == (c[r-1] == r-1)

  if not okc: return False

  # "in other words, R has the form (I_r | T),
  #  where I is an rxr identity matrix."
  for a in range(r):
    for b in range(r):
      assert R[a,b] == (a == b)

  return True

def is_semi_systematic(R,mu,nu):
  r'''
  Return True if R is in systematic form,
  else False.

  INPUT:

  "R" - a matrix in reduced row-echelon form over a field

  "mu" - a nonnegative integer

  "nu" - an integer with nu >= mu
  '''
  c = echelon_positions(R)
  r = len(c)
  assert mu >= 0
  assert nu >= mu

  # "r >= mu"
  if r < mu: return False

  # "there are at least r-mu+nu columns"
  if R.ncols() < r-mu+nu: return False

  # "R has r rows"
  if R.nrows() != r: return False

  # "c_i = i for 0 <= i < r-mu"
  # c_i <= i-mu+nu for 0 <= i < r"
  okc = all(c[i] == i for i in range(r-mu)) and all(c[i] <= i-mu+nu for i in range(r))

  # "equivalent to simply c_{r-mu-1} = r-mu-1 and c_{r-1} <= r-mu+nu-1
  #  except in the degenerate case r = mu"
  if r > mu: assert okc == ((c[r-mu-1] == r-mu-1) and (c[r-1] <= r-mu+nu-1))

  return okc

# ----- miscellaneous tests

def test_manual():
  print('echelon misc')
  sys.stdout.flush()

  X = matrix(GF(2),[[1,1,0],[0,0,1]])
  assert X == echelon(X)
  assert echelon_positions(X) == [0,2]
  assert echelon(matrix(GF(2),[[0,0,1],[1,1,1]])) == X
  assert echelon(matrix(GF(2),[[1,1,1],[1,1,0]])) == X

  X = matrix(GF(2),[[1,0,1],[0,1,1]])
  assert echelon(X) == X
  assert echelon_positions(X) == [0,1]
  assert is_systematic(X)

  X = matrix(GF(2),[[1,1,0,1],[0,0,1,1]])
  assert echelon(X) == X
  assert echelon_positions(X) == [0,2]
  assert not is_semi_systematic(X,0,0)
  assert not is_semi_systematic(X,1,1)
  assert not is_semi_systematic(X,0,1)
  assert is_semi_systematic(X,1,2)
  assert not is_semi_systematic(X,0,2)
  assert is_semi_systematic(X,1,3)

def test_smallrandom():
  for q in range(100):
    q = ZZ(q)
    if not q.is_prime(): continue
    k = GF(q)

    nummatrices = 0
    numsystematic = 0
    numsemisystematic = 0
    for rows in range(10):
      for cols in range(10):
        nummatrices += 1
        X = matrix(k,[[randrange(q) for b in range(cols)] for a in range(rows)])
        R = echelon(X)
        assert R.row_space() == X.row_space()
        assert R.rank() == X.rank()
        c = echelon_positions(R)
        assert len(c) == R.rank()
        if len(c) == R.nrows() and (len(c) == 0 or c[-1] == len(c)-1):
          numsystematic += 1
          assert is_systematic(R)
          assert matrix([Ra[:len(c)] for Ra in R]) == identity_matrix(len(c))
        else:
          assert not is_systematic(R)
          assert matrix([Ra[:len(c)] for Ra in R]) != identity_matrix(len(c))
        if len(c) == R.nrows():
          if len(c) == 0:
            assert is_semi_systematic(R,0,0)
          else:
            assert is_semi_systematic(R,len(c),c[-1]+1)
    print('echelon field %d nummatrices %d numsystematic %d' % (q,nummatrices,numsystematic))
    sys.stdout.flush()

    for mu in range(5):
      for nu in range(mu,mu+5):
        for rows in range(mu,mu+5):
          if rows == 0: continue
          for cols in range(rows-mu+nu,rows-mu+nu+5):
            X = matrix(k,[[a==b if b<rows else randrange(q) for b in range(cols)] for a in range(rows)])
            perm = list(range(nu))
            shuffle(perm)
            perm = [b for b in range(rows-mu)]+[rows-mu+b for b in perm]+[b for b in range(rows-mu+nu,cols)]
            assert len(perm) == cols
            X = X.matrix_from_columns(perm)
            assert is_semi_systematic(echelon(X),mu,nu)

if __name__ == '__main__':
  test_manual()
  test_smallrandom()
