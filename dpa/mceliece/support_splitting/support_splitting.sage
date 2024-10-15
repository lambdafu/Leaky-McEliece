from random import randint

#Need Cython for diagonal standard form, because Sage matrix access is way too slow.
load("support_splitting/diagonal_standard_form.spyx")

"""
    Input:
        Generator matrices G, G_perm of some equivalent codes C, C_perm,
        have to be in systematic form.
    Output:
        Permutation matrix perm, satisfying C(G*perm) = C(G_perm).
"""
def supportSplitting(G, G_perm):
    table = SupportSplittingTable(G)
    table_perm = SupportSplittingTable(G_perm)
    
    hull = table.C.getHullOfPunctured()
    
    n = G.ncols()
    used_indices = []
    
    
    while len(table.discriminatedPositions) < n:
        #Construct set of indices J, at which we will puncture C and C^bot.
        candidates = [j for j in table.discriminatedPositions if j not in used_indices and j in hull.support]
        J = []

        #Algorithm works best, if both C_J and C^bot_J have trivial support.
        #We try to construct such a J.
        if len(candidates) > 0:
          
          for v in hull.B:
            j = 0
            foundPos = False

            while not foundPos and j < n:
              
              if v[j] == 1 and j in candidates:
                foundPos = True
                J.append(j)
                candidates.remove(j)
              
              j += 1
              
        if len(J) < hull.B.nrows():
          #Could not construct J, for which both C_J and C^bot_J have trivial support.
          
          #TODO: Can we do something more clever here?
          candidates = [j for j in table.discriminatedPositions if j not in used_indices]
          j = candidates[randint(0,len(candidates)-1)]
          J = [j]
            
        J_perm = table_perm.indices(table.rows(J))
        
        used_indices += J
      
        table.update(J)
        table_perm.update(J_perm)
        
    perm = zero_matrix(GF(2), n)
    
    for i in range(n):
        i_perm = table_perm.indices(table.rows([i]))[0]
        
        perm[i,i_perm] = 1
    
    return perm

class SupportSplittingTable:
    """
        Input:
            Generator matrix G of some code C,
            has to be in systematic form.
    """
    def __init__(self, G):
        self.C = Code(G)
        self.discriminatedPositions = []
        self.table = [ [] for _ in range(self.C.length) ]
        
        self.update()
    
    """
        Input:
            List J.
            
        Updates support splitting table by adding signature of code punctured at J.
    """
    def update(self, J = []):
        hull = self.C.getHullOfPunctured(J)
        hull_dual = self.C.getHullOfPunctured(J, dual = True)
        
        
        #TODO: Parallelize?
        for i in [i for i in range(self.C.length) if i not in self.discriminatedPositions]:
          sig = self.__signatures(i,hull,hull_dual)

          self.table[i].append( sig )
          
        reversedTable = {}
        for i in range(self.C.length):
          if i not in self.discriminatedPositions:
            key = repr(self.table[i])
            if key in reversedTable:
              reversedTable[key].append(i)
            else:
              reversedTable[key] = [i]
        
        for key in reversedTable:
          if len( reversedTable[key] ) == 1:
            self.discriminatedPositions.append( reversedTable[key][0] )
        
        
    """
        Input:
          Integer i,
          hull hull_1 of some code C,
          hull hull_2 of some code C'.
        Output:
          ( W(H(C_i)), W(H(C^bot_i)), W(H(C'_i), W(H((C'^bot_i)) )
    """
    def __signatures(self,i,hull_1,hull_2):
      S_1, S_2 = hull_1.signature(i)
      S_3, S_4 = hull_2.signature(i)
      
      makeUniqueHashable = lambda s : tuple( sorted( s.items() ) )
      
      return (
        makeUniqueHashable(S_1),
        makeUniqueHashable(S_2),
        makeUniqueHashable(S_3),
        makeUniqueHashable(S_4)
      )
    
    """
        Input:
            List J = [j_1,...,j_m].
        Output:
            List of containing the j_1,...,j_m-th rows of support splitting table.
    """
    def rows(self, J):
        return [ self.table[j] for j in J ]
    
    """
        Input:
            List L of rows r_1,...,r_m contained in support splitting table.
        Output:
            Indices of r_1,...,r_m.
    """
    def indices(self, L):
        return [ self.table.index(row) for row in L ]
        
class Code:
    """
        Input:
            Generator matrix G of some code C,
            has to be in systematic form.
    """
    def __init__(self, G):
        self.dim = G.nrows()
        self.length = G.ncols()
        self.codim = self.length - self.dim
        
        if G.submatrix(0,0,self.dim,self.dim) != identity_matrix(GF(2), self.dim):
            raise ValueError("SupportSplittingTable expects G to be in systematic form.")
        
        #Basis of code C.
        self.G = G
        
        A = self.G.submatrix(0,self.dim,self.dim,self.codim)
        A_t = A.transpose()
        
        #Basis of dual code C^bot.
        self.H = A_t.augment( identity_matrix(GF(2),self.codim) ) 

        #Basis of C + C^bot.
        self.M = self.G.stack(self.H)

        X = identity_matrix(GF(2),self.codim) + A_t * A
        E,U = diagonalStandardForm(X, identity_matrix(GF(2),self.codim))
        AU = A*U

        #Diagonal standard form of M
        self.D = block_matrix([
            [identity_matrix(GF(2),self.dim),A+A*E],
            [zero_matrix(GF(2),self.codim,self.dim),E]
        ])
        
        # S*M = D
        self.S = block_matrix([
            [identity_matrix(GF(2),self.dim)-AU*A_t,AU],
            [U*A_t,U]
        ])
    
    """
        Input:
            List J, boolean dual (optional).
        Output:
            Hull of code punctured at J.
            
            If dual = True, then hull of dual code punctured at J is returned.
    """
    def getHullOfPunctured(self, J=[], dual = False):
        G = self.G
        H = self.H
        
        if dual:
            H = self.__puncture(H,J)
        else:
            G = self.__puncture(G,J)
        
        M = G.stack(H)
        X = self.M+M
        
        D, S = diagonalStandardForm(self.D + self.S*X,self.S)
        
        #Basis of hull
        B = (D + identity_matrix(GF(2),self.length)).transpose()
        B = self.__puncture(B,J)
        B = self.__removeZeroRows(B)
        
        
        #Support of hull
        support = set()
        for v in B:
            for i,v_i in enumerate(v):
                if v_i == 1:
                    support.add(i)
        support = list(support)
        
        #Splitting unit vectors into code vectors and dual vectors.
        L = {}
        
        S_left = S.submatrix(0,0,self.length,self.dim)
        S_right = S.submatrix(0,self.dim,self.length,self.codim)
        
        #TODO: Speed-up possible?
        if dual:
          X = S_right*H
        else:
          X = S_left*G
            
        I_n = identity_matrix(GF(2), self.length) 
        
        for i in [i for i in range(self.length) if i not in support and i not in J]:
            e_i = I_n[i]
            
            if dual:
                x = X[i]
            else:
                x = X[i]
            
            y = x + e_i
            
            L[i] = (x,y)
        

        return Hull(B, support, L)
    
    """
        Input:
            Matrix B, list J.
        Output:
            B punctured at J.
    """
    def __puncture(self,B,J):
        B = copy(B)
        
        for j in J:
            for i in range(B.nrows()):
                B[i,j] = 0
        return B
    
    """
        Input:
            Matrix B.
        Output
            If B != 0, returns a copy of B, where all zero rows are removed.
            If B = 0, then a single zero row is returned.
    """
    def __removeZeroRows(self,B):
        B_ = []
        if B.is_zero():
            B_ = [ B[0] ]
        else:
            for b in B:
                if not b.is_zero():
                    B_.append(b)
        return matrix(B_)

class Hull:
    """
        Input:
            Basis B of some hull H(C),
            support of H(C),
            dict L of tuples L[i] = (x,y), s.t. x in C, y in C^bot, and x+y = e_i.
    """
    def __init__(self, B, support, L):
        self.B = B
        self.support = support
        self.L = L

        #TODO: It does not make sense to pass support and L as arguments. We should move the stuff from Code.getHullOfPunctured() here.
        
        self.space = list(B.image()) #Casting to list seems to make for-loop faster.
        
        self.WeightEnumerator = {}
        for v in self.space:
            self.__addToWeightEnumerator(self.WeightEnumerator, v)
            
    """
        Input:
            Integer i.
        Output:
            Weight locator of hulls H(C_i) and H(C^bot_i).
    """
    def signature(self,i):
        if i in self.support:
            shortened_hull_wt = [ v.hamming_weight() for v in self.space if v[i] != 0 ]
            
            W = copy(self.WeightEnumerator)
            for wt in shortened_hull_wt:
                W[wt] -= 1
            
            W_dual = W
            
        elif i in self.L:
            x,y = self.L[i]
            
            if x[i]==1:
                W = self.__getWeightEnumeratorOfSum(y)
                W_dual = self.WeightEnumerator
            else:
                W_dual = self.__getWeightEnumeratorOfSum(x)
                W = self.WeightEnumerator
        
        return(W,W_dual)
    
    """
          Input:
            Vector x.
          Output:
            Weight enumerator of H(C + GF(2) * x)
    """
    def __getWeightEnumeratorOfSum(self, x):
      W = copy(self.WeightEnumerator)
      for v in self.space:
        self.__addToWeightEnumerator(W,v+x)
      
      return W
    
    """
      Input:
        Weight enumerator W of some code C',
        vector v.
      
      Adds v in place to W.
    """
    def __addToWeightEnumerator(self, W, v):
      wt = v.hamming_weight()
      if wt in W:
        W[wt] += 1
      else:
        W[wt] = 1