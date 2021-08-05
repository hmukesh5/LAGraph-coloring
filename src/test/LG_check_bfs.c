//------------------------------------------------------------------------------
// LAGraph/src/test/LG_check_bfs: stand-alone test for BFS
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
//
// See additional acknowledgments in the LICENSE file,
// or contact permission@sei.cmu.edu for the full terms.

#define LAGRAPH_FREE_WORK                   \
{                                           \
    LAGraph_Free ((void **) &queue) ;       \
    LAGraph_Free ((void **) &level_check) ; \
    LAGraph_Free ((void **) &level_in) ;    \
    LAGraph_Free ((void **) &parent_in) ;   \
    LAGraph_Free ((void **) &visited) ;     \
    LAGraph_Free ((void **) &neighbors) ;   \
    GrB_free (&Row) ;                       \
}

#define LAGRAPH_FREE_ALL                    \
{                                           \
    LAGRAPH_FREE_WORK ;                     \
    LAGraph_Free (&Ap) ;                    \
    LAGraph_Free (&Aj) ;                    \
    LAGraph_Free (&Ax) ;                    \
}

#include "LG_internal.h"
#include "LG_test.h"

//------------------------------------------------------------------------------
// extract the contents of an int64 vector
//------------------------------------------------------------------------------

bool get_vector
(
    int64_t *x,
    GrB_Vector X,
    int64_t n
)
{
    for (int64_t i = 0 ; i < n ; i++)
    {
        int64_t t ;
        int info = GrB_Vector_extractElement_INT64 (&t, X, i) ;
        if (info == GrB_SUCCESS)
        {
            x [i] = t ;
        }
        else if (info == GrB_NO_VALUE)
        {
            x [i] = -1 ;
        }
        else
        {
            return (false) ;    // method failed
        }
    }
    return (true) ;             // success
}

//------------------------------------------------------------------------------
// test the results from a BFS
//------------------------------------------------------------------------------

int LG_check_bfs
(
    // input
    GrB_Vector Level,       // optional; may be NULL
    GrB_Vector Parent,      // optional; may be NULL
    LAGraph_Graph G,
    GrB_Index src,
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    GrB_Vector Row = NULL ;
    GrB_Index *Ap = NULL, *Aj = NULL ;
    void *Ax = NULL ;
    GrB_Index Ap_size, Aj_size, Ax_size, n, ncols ;
    bool iso, jumbled ;
    int64_t *queue = NULL, *level_in = NULL, *parent_in = NULL,
        *level_check = NULL, *neighbors = NULL ;
    bool *visited = NULL ;
    LG_CHECK (LAGraph_CheckGraph (G, msg), -1, "graph is invalid") ;
    GrB_TRY (GrB_Matrix_nrows (&n, G->A)) ;
    GrB_TRY (GrB_Matrix_ncols (&ncols, G->A)) ;
    LG_CHECK (n != ncols, -1001, "G->A must be square") ;

//  printf ("\nTEST BFS\n") ;
//  GrB_TRY (GxB_Matrix_fprint (G->A, "adjacency matrix", 2, stdout)) ;
//  if (Level != NULL)
//  {
//      GrB_TRY (GxB_Vector_fprint (Level, "level", 2, stdout)) ;
//  }
//  if (Parent != NULL)
//  {
//      GrB_TRY (GxB_Vector_fprint (Parent, "parent", 2, stdout)) ;
//  }

    //--------------------------------------------------------------------------
    // allocate workspace
    //--------------------------------------------------------------------------

    queue = LAGraph_Malloc (n, sizeof (int64_t)) ;
    level_check = LAGraph_Malloc (n, sizeof (int64_t)) ; 
    LG_CHECK (queue == NULL || level_check == NULL , -1003, "out of memory") ;

    //--------------------------------------------------------------------------
    // get the contents of the Level and Parent vectors
    //--------------------------------------------------------------------------

    if (Level != NULL)
    {
        level_in = LAGraph_Malloc (n, sizeof (int64_t)) ;
        LG_CHECK (level_in == NULL, -1003, "out of memory") ;
        LG_CHECK (!get_vector (level_in, Level, n), -1004, "invalid level") ;
    }

    if (Parent != NULL)
    {
        parent_in = LAGraph_Malloc (n, sizeof (int64_t)) ;
        LG_CHECK (parent_in == NULL, -1003, "out of memory") ;
        LG_CHECK (!get_vector (parent_in, Parent, n), -1005, "invalid parent") ;
    }

    //--------------------------------------------------------------------------
    // compute the level of each node
    //--------------------------------------------------------------------------

    queue [0] = src ;
    int64_t head = 0 ;
    int64_t tail = 1 ;
    visited = LAGraph_Calloc (n, sizeof (bool)) ;
    LG_CHECK (visited == NULL, -1003, "out of memory") ;
    visited [src] = true ;      // src is visited, and is level 0

    for (int64_t i = 0 ; i < n ; i++)
    {
        level_check [i] = -1 ;
    }
    level_check [src] = 0 ;

    GrB_TRY (GrB_Vector_new (&Row, GrB_BOOL, n)) ;
    neighbors = LAGraph_Malloc (n, sizeof (GrB_Index)) ;
    LG_CHECK (neighbors == NULL, -1003, "out of memory") ;

    while (head < tail)
    {
        // dequeue the node at the head of the queue
        int64_t u = queue [head++] ;

        // extract the indices of entries in A(u,:)
        GrB_TRY (GrB_Col_extract (Row, NULL, NULL, G->A, GrB_ALL, n, u,
            GrB_DESC_T0)) ;
        int64_t degree = n ;
        GrB_TRY (GrB_Vector_extractTuples_BOOL (neighbors, NULL, &degree,
            Row)) ;

        // traverse all entries in A(u,:)
        for (int64_t k = 0 ; k < degree ; k++)
        {
            // consider edge (u,v)
            int64_t v = neighbors [k] ;
            if (!visited [v])
            {
                // node v is not yet visited; set its level and add to the
                // end of the queue
                visited [v] = true ;
                level_check [v] = level_check [u] + 1 ;
                queue [tail++] = v ;
            }
        }
    }

    //--------------------------------------------------------------------------
    // check the level of each node
    //--------------------------------------------------------------------------

    if (level_in != NULL)
    {
        for (int64_t i = 0 ; i < n ; i++)
        {
            bool ok = (level_in [i] == level_check [i]) ;
            LG_CHECK (!ok, -1004, "invalid level") ;
        }
    }

    //--------------------------------------------------------------------------
    // check the parent of each node
    //--------------------------------------------------------------------------

    if (parent_in != NULL)
    {
        for (int64_t i = 0 ; i < n ; i++)
        {
            if (i == src)
            {
                // src node is its own parent
                bool ok = (parent_in [src] == src) && (visited [src]) ;
                LG_CHECK (!ok, -1005, "invalid parent") ;
            }
            else if (visited [i])
            {
                int64_t pi = parent_in [i] ;
                // ensure the parent pi is valid and has been visited
                bool ok = (pi >= 0 && pi < n) && visited [pi] ;
                LG_CHECK (!ok, -1005, "invalid parent") ;
                // ensure the edge (pi,i) exists
                bool x ;
                int info = GrB_Matrix_extractElement_BOOL (&x, G->A, pi, i) ;
                ok = (info == GrB_SUCCESS) ;
                LG_CHECK (!ok, -1005, "invalid parent") ;
                // ensure the parent's level is ok
                ok = (level_check [i] == level_check [pi] + 1) ;
                LG_CHECK (!ok, -1005, "invalid parent") ;
            }
        }
    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    LAGRAPH_FREE_WORK ;
    return (0) ;
}
