#include "LG_internal.h"
#include "LAGraphX.h"
// add this algorithm to LAGraphX.h

#define LG_FREE_WORK                \
    GrB_free (&local_color) ;       \
    GrB_free (&weight) ;            \
    GrB_free (&in_curr_subset) ;    \
    GrB_free (&max_weights) ;

#define LG_FREE_ALL                 \
    LG_FREE_WORK ;                  \
    GrB_free (&color) ;

int LAGraph_coloring_independent_set_optimized
(
    // output
    GrB_Vector *color,
    int *num_colors,

    // input
    LAGraph_Graph G,
    char *msg
)
{
    bool verbose = false;
    GrB_Vector local_color = NULL;
    GrB_Vector weight = NULL;
    GrB_Vector in_curr_subset = NULL;
    GrB_Vector max_weights = NULL;

    GrB_Index n;
    GRB_TRY (GrB_Matrix_nrows(&n, G->A)) ;

    GrB_Type Int = (n < UINT32_MAX) ? GrB_UINT32 : GrB_UINT64 ;

    /* initialize local copy of color to SPARSE vector */
    GRB_TRY(GrB_Vector_new(&local_color, Int, n));

    // lg_set_format_hint -> bitmap

    /* weights initialized randomly
    *  seed of 20 was chosen arbitrarily */   
    GRB_TRY(GrB_Vector_new(&weight, GrB_UINT64, n));
    GRB_TRY(GrB_assign (weight, NULL, NULL, 0, GrB_ALL, n, NULL));
    LAGraph_Random_Seed(weight, 20, msg);

    GRB_TRY(GrB_Vector_new(&in_curr_subset, GrB_BOOL, n));

    GRB_TRY(GrB_Vector_new(&max_weights, GrB_UINT64, n));

    /* algorithm start */
    int64_t curr_color;
    for (curr_color = 1; curr_color < n+1; curr_color++) {
        /* mxv - find maximum of all neighboring weights */

        // FIXME: try using a set of sparse candidate nodes, not yet colored
        GRB_TRY(GrB_mxv(max_weights, local_color, GrB_NULL,
            GrB_MAX_SECOND_SEMIRING_UINT64, G->A, weight, GrB_DESC_RSC));

        /* eWiseAdd - 1 if current weight > max neighboring weight */
        GRB_TRY(GrB_eWiseMult(in_curr_subset, GrB_NULL, GrB_NULL, GrB_GT_UINT64, weight, max_weights, GrB_NULL));
        /* select - select all entries in in_curr_subset that are true, and delete falses */
        GRB_TRY(GrB_select(in_curr_subset, GrB_NULL, GrB_NULL, GrB_VALUEEQ_BOOL, in_curr_subset, true, GrB_NULL));

        /* reduce - OR all entries in in_curr_subset - if false, break */
        // FIXME: just check nvals(in_curr_subset)
        bool subset_exists;
        GRB_TRY(GrB_reduce(&subset_exists, GrB_NULL, GrB_LOR_MONOID_BOOL, in_curr_subset, GrB_NULL));
        if (subset_exists == false) { break; }

        // FIXME: future: if in_curr_subset is empty, but nvals (local_color) < n, then BROKEN

        /* assign - write current color to C vector according to in_curr_subset mask */
        GRB_TRY(GrB_assign(local_color, in_curr_subset, GrB_NULL, curr_color, GrB_ALL, n, GrB_DESC_S));
        
        /* assign - write 0 to weight according to in_curr_subset mask */
        GRB_TRY(GrB_assign(weight, in_curr_subset, GrB_NULL, 0, GrB_ALL, n, GrB_DESC_S));
    }

    (*num_colors) = curr_color - 1;
    (*color) = local_color;
    local_color = NULL ;
    LG_FREE_ALL ;
    return (GrB_SUCCESS) ;
}
