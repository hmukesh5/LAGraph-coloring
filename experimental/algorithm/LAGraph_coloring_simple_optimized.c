#include "LG_internal.h"
#include "LAGraphX.h"
// add this algorithm to LAGraphX.h

int LAGraph_coloring_simple_optimized
(
    // output
    GrB_Vector *color, // lowercase "color"
    // add number of colors as output

    // input
    LAGraph_Graph G,
    char *msg
)
{
    bool verbose = true;

    GrB_Index n;
    GrB_Matrix_nrows(&n, G->A);

    /* initialize local copy of color to SPARSE vector */
    GrB_Vector local_color = NULL;
    GRB_TRY(GrB_Vector_new(&local_color, GrB_INT32, n));

    // lg_set_format_hint -> bitmap

    /* weights initialized randomly
    *  seed of 20 was chosen arbitrarily */   
    GrB_Vector weight = NULL;
    GRB_TRY(GrB_Vector_new(&weight, GrB_UINT64, n));
    GRB_TRY(GrB_assign (weight, NULL, NULL, 0, GrB_ALL, n, NULL));
    LAGraph_Random_Seed(weight, 20, msg);
    /* DEBUG PRINT */ if (verbose) { printf("\n[ DEBUG ] weight vector\n"); GxB_print(weight, 3); }

    GrB_Vector in_curr_subset = NULL;
    GRB_TRY(GrB_Vector_new(&in_curr_subset, GrB_BOOL, n)); // bool

    GrB_Vector max_weights = NULL;
    GRB_TRY(GrB_Vector_new(&max_weights, GrB_UINT64, n));

    GrB_Scalar reduced_scalar = NULL;
    GrB_Scalar_new(&reduced_scalar, GrB_INT32);

    GrB_Scalar color_scalar = NULL;
    GrB_Scalar_new(&color_scalar, GrB_INT32);
    
    /* algorithm start */
    for (int curr_color = 1; curr_color < n+1; curr_color++) {
        /* vxm - find maximum of all neighboring weights */
        GRB_TRY(GrB_mxv(max_weights, local_color, GrB_NULL, GrB_MAX_SECOND_SEMIRING_UINT64, G->A, weight, GrB_DESC_RSC));
        /* DEBUG PRINT */ if (verbose) { printf("ITERATION %d\n[DEBUG] max_weights\n", curr_color); GxB_print(max_weights, 3); }

        /* eWiseAdd - 1 if current weight > max neighboring weight */
        GRB_TRY(GrB_eWiseMult(in_curr_subset, GrB_NULL, GrB_NULL, GrB_GT_UINT64, weight, max_weights, GrB_NULL));
        /* DEBUG PRINT */ if (verbose) { printf("[DEBUG] in_curr_subset\n"); GxB_print(in_curr_subset, 3); }

        /* reduce - OR all entries in in_curr_subset - if false, break */
        bool subset_exists;
        GRB_TRY(GrB_reduce(&subset_exists, GrB_NULL, GrB_LOR_MONOID_BOOL, in_curr_subset, GrB_NULL));
        /* DEBUG PRINT */ if (verbose) { printf("[DEBUG] subset_exists: %d\n", subset_exists); }
        if (subset_exists == false) { break; }

        /* assign - write current color to C vector according to in_curr_subset mask */
        GRB_TRY(GrB_assign(local_color, in_curr_subset, GrB_NULL, curr_color, GrB_ALL, n, GrB_NULL));
        if(verbose) { printf("[DEBUG] color vector, current color: %d\n", curr_color); GxB_print(local_color, 3);}
        
        /* assign - write 0 to weight according to in_curr_subset mask */
        GRB_TRY(GrB_assign(weight, in_curr_subset, GrB_NULL, 0, GrB_ALL, n, GrB_NULL));
        /* DEBUG PRINT */ if (verbose) { printf("[DEBUG] new weight\n"); GxB_print(weight, 3);}
    }

    (*color) = local_color;

    return (GrB_SUCCESS) ;
}
