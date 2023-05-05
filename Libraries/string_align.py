import numpy as np

def sequence_align(seq1, seq2, penmatrix, pengap):
    #print str(penmatrix)
    #print penmatrix
   # print pengap
    m = len(seq1)
    n = len(seq2)
    score = np.zeros([m+1, n+1], dtype=np.float)
    pointer = np.zeros([m+1, n+1], dtype=long)
    
    max_score = 0.0
    max_i = 0
    max_j = 0  

    # Precalculo
    for i in xrange(1, m+1):
        score[i, 0] = i*pengap
    for j in xrange(1, n+1):
        score[0, j] = j*pengap
    for i in xrange(1, m+1):
        for j in xrange(1, n+1):
            score_up = score[i-1, j] + pengap
            score_down = score[i, j-1] + pengap
            #print penmatrix(seq1[i-1],seq2[j-1])
            score_diagonal = score[i-1, j-1] + penmatrix(seq1[i-1],seq2[j-1])
            
            score[i, j] = max(0, score_up, score_down, score_diagonal)
            
            if score[i, j] == score_diagonal:
                pointer[i, j] = 3 # Venimos de la diagonal
            elif score[i, j] == score_down:
                pointer[i, j] = 2 # Venimos de la izquierda
            elif score[i, j] == score_up:
                pointer[i, j] = 1 # Venimos de arriba
            elif score[i, j] == 0:
                pointer[i, j] = 1 # No sigas...
                
            if score[i, j] >= max_score:
                max_i = i
                max_j = j
                max_score = score[i, j]
                
    # Ya podemos calcular el alineamiento
    i, j = max_i, max_j
    
    align1 = []
    align2 = []
    while pointer[i, j] != 0:
        if pointer[i, j] == 3:
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i = i-1
            j = j-1
        elif pointer[i, j] == 2:
            align1.append("-")
            align2.append(seq2[j-1])
            j = j-1
        elif pointer[i, j] == 1:
            align1.append(seq1[i-1])
            align2.append("-")
            i = i-1
        
    align1 = align1[::-1]
    align2 = align2[::-1] 
        
    while i > 0:
        align1[0:0] = [seq1[i-1]]
        i = i-1
    while j > 0:
        align2[0:0] = [seq2[j-1]]
        j = j-1
        
    MAX = max(len(align1), len(align2)) 
    while len(align1) < MAX:
        align1[0:0] = ["-"]
    while len(align2) < MAX:
        align2[0:0] = ["-"]  
    
    return align1, align2, max_score/max(len(seq1), len(seq2))
