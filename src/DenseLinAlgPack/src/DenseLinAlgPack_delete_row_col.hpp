// //////////////////////////////////////////////////////////////////
// delete_row_col.h

#ifndef DELETE_ROW_COL_H
#define DELETE_ROW_COL_H

#include "LinAlgPackTypes.h"

namespace LinAlgPack {

///
/** Delete a symmetric row and a column form a triangular matrix.
 *
 * If #M# is a lower triangular matrix then we partition it
 * as:
 \begin{verbatim}

   1 |\
     |  \
     |    \
     | M11  \
     |________\ _
  kd |_________|_|
     |         | |\
     |         | |  \
     |         | |    \
     |   M31   | | M33  \
   n |         | |        \
     ----------------------
     1         kd         n

 \end{verbatim}
 * In order to delete row #kd# and column #kd# the rectangular
 * matrix #M31# is moved up one row and the triangular matrix
 * #M33# is moved up one row and to the left one column.
 *
 * If #M# is an upper triangular matrix then we partition it
 * as:
 \begin{verbatim}

  1         kd      n
  -------------------- 1
  \        | |       |
    \  M11 | |  M13  |
      \    | |       |
        \  | |       |
          \|_|_______|
           |_|_______| kd
             \       |
               \ M33 |
                 \   |
                   \ | n
 
 \end{verbatim}
 *
 * In order to delete row #kd# and column #kd# the matrix
 * #M13# is moved one column to the left and the upper
 * triangular matrix #M33# is moved one row up and
 * on column to the left.
 *
 * Preconditions:\begin{itemize}
 * \item #M != NULL#
 * \item #M->rows() >= 1#
 * \item #1 <= kd <= M->rows()#
 * \end{itemize}
 */
void delete_row_col( size_type kd, tri_ele_gms* M );

} // end namespace LinAlgPack

#endif  // DELETE_ROW_COL_H
