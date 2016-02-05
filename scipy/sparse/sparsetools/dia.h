#ifndef __DIA_H__
#define __DIA_H__

#include <algorithm>


/*
 * Compute Y += A*X for DIA matrix A and dense vectors X,Y
 *
 *
 * Input Arguments:
 *   I  n_row            - number of rows in A
 *   I  n_col            - number of columns in A
 *   I  n_diags          - number of diagonals
 *   I  L                - length of each diagonal
 *   I  offsets[n_diags] - diagonal offsets 
 *   T  diags[n_diags,L] - nonzeros 
 *   T  Xx[n_col]        - input vector
 *
 * Output Arguments:
 *   T  Yx[n_row]        - output vector 
 *
 * Note:
 *   Output array Yx must be preallocated
 *   Negative offsets correspond to lower diagonals
 *   Positive offsets correspond to upper diagonals
 *
 */
template <class I, class P, class T>
void dia_matvec(const I n_row,
                const I n_col,
                const I n_diags,
                const I L,
	            const P offsets[], 
	            const T diags[], 
	            const T Xx[],
	                  T Yx[])
{
    for(P i = 0; i < n_diags; i++){
        const P k = offsets[i];  //diagonal offset

        const P i_start = std::max<P>(0,-k);
        const P j_start = std::max<P>(0, k);
        const P j_end   = std::min<P>(std::min<P>(n_row + k, n_col),L);

        const P N = j_end - j_start;  //number of elements to process

        const T * diag = diags + (npy_intp)i*L + j_start;
        const T * x = Xx + j_start;
              T * y = Yx + i_start;

        for(P n = 0; n < N; n++){
            y[n] += diag[n] * x[n]; 
        }
    }
}


#endif
