/* wrapper routines for basic sparse linear algebra 
 *
 * PCx 1.1 11/97
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Michael Wagner, Steve Wright.
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

#include <stdio.h>
#include "main.h"
#include "memory.h"


//  #define max(a,b) \
//    ({ __typeof__ (a) _a = (a); \
//        __typeof__ (b) _b = (b); \
//      _a > _b ? _a : _b; })
// #define abs(a) \
//    ({ __typeof__  (a) _a = (a);\
//      _a > 0 ? _a : - _a; }) 

/* copies a sparseMatrix type into an MMTtype, but doesn't fill in the three
 * index arrays for the transpose structure (this is done in a later routine
 * ) */

MMTtype        *copyMMT(A, NumRows, NumCols)
  sparseMatrix    A;
  int             NumRows, NumCols;
{
  int             i, Nonzeros;
  MMTtype        *Anew, *NewMMTtype();

  Nonzeros = A.pEndRow[NumCols - 1];
  Anew = NewMMTtype(NumRows, NumCols, Nonzeros);
  for (i = 0; i < NumCols; i++) {
    Anew->pBeginRow[i] = A.pBeginRow[i];
    Anew->pEndRow[i] = A.pEndRow[i];
  }
  for (i = 0; i < Nonzeros; i++) {
    Anew->Row[i] = A.Row[i];
    Anew->Value[i] = A.Value[i];
  }
  return Anew;
}

/* wrapper for RealSparseMatrixVectorProductPlusx() */
int             Axplusb(A, x, b)
  MMTtype        *A;
  double         *x, *b;

{
  return RealSparseMatrixVectorProductPlusx
    (A->Value, A->pBeginRow, A->pEndRow, A->Row, x, b,
     &(A->NumRows), &(A->NumCols));
}


/* wrapper for RealSparseMatrixVectorProduct */
int             Ax(A, x, b)
  MMTtype        *A;
  double         *x, *b;

{
  return RealSparseMatrixVectorProduct
    (A->Value, A->pBeginRow, A->pEndRow, A->Row, x, b,
     &(A->NumRows), &(A->NumCols));
}


/* The following two routines are for A of type sparseMatrix.  These are
 * usually called with LP->A as the argument, where LP is an LPtype data
 * structure.  */

int             SparseSaxpy(A, x, y)
  sparseMatrix    A;
  double         *x, *y;

{
  int             m, n;

  m = A.NumRows;
  n = A.NumCols;

  return RealSparseMatrixVectorProductPlusx
    (A.Value, A.pBeginRow, A.pEndRow, A.Row,
     x, y, &m, &n);
}


int             SparseSaxpyT(A, x, y)
  sparseMatrix    A;
  double         *x, *y;

{
  int             m, n;

  m = A.NumRows;
  n = A.NumCols;

  return RealSparseMatrixTransposeVectorProductPlusx
    (A.Value, A.pBeginRow, A.pEndRow, A.Row,
     x, y, &m, &n);
}


/* The following two routines are for A of type *MMTtype.  */

int             SparseSaxpyM(A, x, y)
  MMTtype        *A;
  double         *x, *y;

{
  int             m, n;

  m = A->NumRows;
  n = A->NumCols;

  return RealSparseMatrixVectorProductPlusx
    (A->Value, A->pBeginRow, A->pEndRow, A->Row,
     x, y, &m, &n);
}


int             SparseSaxpyTM(A, x, y)
  MMTtype        *A;
  double         *x, *y;

{
  int             m, n;

  m = A->NumRows;
  n = A->NumCols;

  return RealSparseMatrixTransposeVectorProductPlusx
    (A->Value, A->pBeginRow, A->pEndRow, A->Row,
     x, y, &m, &n);
}


double          TwoNorm2(x, n)
  double         *x;
  int            *n;

{
  double          temp;

  if (*n <= 0)
    return 0.0;
  NormTwoSquareRealDenseVector(x, n, &temp);
  return temp;
}


double          MaxNormVector(x, n)
  double         *x;
  int            n;

{
  double          temp;
  if (n <= 0)
    return 0.0;
  static double  *count;
  count = x + n;
  temp = absLRS(*x);
  int teste;
  teste = 1;
   for ( ; x < count ; x++ ){
      temp = maxLRS(temp,absLRS(*(x+1)));
      }      
  return temp;
}

/* MaxNormMatrix.c, v1.0, 08/26/91, Sanjay Mehrotra */

/* The function in this file computes a matrix vector product.  The matrix a
 * is assumed to be stored column-wise and the given vector q is assumed to
 * be dense. The result is written in vector aq.  The begining and end of
 * column i in a is stored in ipbra[i] and ipera[i].  The row indices
 * corresponding to each column in a are stored in ira[]. All the array
 * information is passed through pointers.  The pointers point to the first
 * element of the array. The subroutine assumes that the row/column indices
 * run from 1..m/n. As oppose to 0..m-1, which is the standard C.  All
 * changes required are done internally in a subroutine.  */

/* WARNING: ROW/COLUMN INDEX ZERO IS NOT ALLOWED. */

//           (A->Value, A->pBeginRow, A->pEndRow, A->Row, x, b,
//      &(A->NumRows), &(A->NumCols))
// double             
// MaxNormMatrix(a_p, ipbra_p, ipera_p, ira_p, q_p, aq_p,
//             nrow_p, ncol_p)
//      double         *a_p, *q_p, *aq_p;
//      int            *ipbra_p, *ipera_p, *ira_p;
//      int            *ncol_p, *nrow_p;
// {
//    int             error_code;
//    int            *ira_ps;
//    int             jbeg, jend, irow;
//    double         *q_ps, *aq_ps, *a_ps;
//    double          qval, maxnorm;
   
//    /* initialize *nrow_p elements of aq */
   
//    if (ZeroRealDenseVector(aq_p, nrow_p)) 
//       {
//    error_code = 2;
//    fprintf(stdout, "Error: RealSparseMatrixVectorProduct: Error\n");
//    fprintf(stdout, " in ZeroRealDenseVector initializing vector.\n");
//    return error_code;
//       }
//    /* indices start from 1 not zero so shift relevant pointers */
   
//    aq_ps = aq_p - 1;
//    a_ps = a_p - 1;
//    ira_ps = ira_p - 1;
   
//    /* now let us start the matrix vector product accumulate one column at a
//     * time.  */
//    maxnorm = 0.0;
//    q_ps = q_p + *ncol_p;
//    for (; q_p < q_ps; q_p++) 
//       {
//    qval = *q_p;
//    jbeg = *ipbra_p++;
//    jend = *ipera_p++;
//       for (; jbeg <= jend; jbeg++) 
//         {
//            irow = *(ira_ps + jbeg);
//            maxnorm = max(maxnorm, abs(*(a_ps + jbeg)));
//           *(aq_ps + irow) += *(a_ps + jbeg) * qval;
//         }
//         }
//    return maxnorm;
// }