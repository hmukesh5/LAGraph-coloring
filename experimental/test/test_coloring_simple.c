// no idea about these 2 includes
#include <stdio.h>
#include <acutest.h>

// includes from LAGraph, also dunno
#include <LG_internal.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>

// error msg buffer
char msg[LAGRAPH_MSG_LEN];

LAGraph_Graph G = NULL;

#define LEN 512
char filename[LEN + 1];

void test_coloring_simple(void)
{
    /* defines whether my [ DEBUG ] messages are printed */
    bool verbose = false;

    /* required initialization (found from other test files) */
    LAGraph_Init(msg);
    LAGraph_Random_Init(msg);

    /* initializing A (matrix) and C (color vector) */
    GrB_Matrix A = NULL;
    GrB_Vector C = NULL;

    /* open matrix market file */
    snprintf(filename, LEN, LG_DATA_DIR "%s", "ldbc-undirected-example-unweighted.mtx");
    FILE *f = fopen(filename, "r");
    TEST_CHECK(f != NULL);
    OK(LAGraph_MMRead(&A, f, msg));
    OK(fclose(f));
    OK(LAGraph_New(&G, &A, LAGraph_ADJACENCY_UNDIRECTED, msg));
    TEST_CHECK(A == NULL); // A has been moved into G->A

    /* run the algorithm */
    /* DEBUG PRINT */ if (verbose) { printf("\n"); LAGraph_Matrix_Print(G->A, LAGraph_SHORT, stdout, msg); printf("[ DEBUG ] running alg"); }

    GRB_TRY (GxB_set (GxB_BURBLE, false)) ;

    double time = LAGraph_WallClockTime();
    LAGraph_coloring_simple_optimized(&C, G, msg);
    time = LAGraph_WallClockTime() - time;

    GRB_TRY (GxB_set (GxB_BURBLE, false)) ;

    printf("\nTook %g seconds\n", time);
    /* DEBUG PRINT */ if (verbose) { printf("[ DEBUG ] finished alg, color vector:\n"); LAGraph_Vector_Print(C, LAGraph_SHORT, stdout, msg); }

    /* extract graph in CSC form
    *  CSC form:
    *  Ap: start and end indices of Ai that represent a column of A
    *  Ai: values that represent the row indices of A where there is a value
    *       - in our case, these are the neighbors' IDs
    *  Ax: the values of that edge, stored in the same order as Ai
    *       - in our case, all values are 1
    *  note: CSC is same as CSR for undirected graphs
    *        maybe use the one that's more efficient
    *        ( prevent converting from one to other )
    */
    GrB_Index *Ap = NULL;
    GrB_Index *Ai = NULL;
    void *Ax = NULL;
    GrB_Index Ap_size, Ai_size, Ax_size;
    GRB_TRY(GxB_Matrix_unpack_CSC(G->A, &Ap, &Ai, &Ax, &Ap_size, &Ai_size, &Ax_size, NULL, NULL, NULL));
    
    /*  convert Ap_size from bytes to indices
    *   - make sure to loop only up to Ap_size - 1
    *     when checking [Ap_index] to [Ap_index + 1]
    */
    Ap_size = Ap_size / sizeof(GrB_Index);

    /* traverse through unpacked matrix and
    *  check current node's color against its neighbors
    *   - Ap_index: current node
    *   - Ai_index: a neighbor
    */
   /* DEBUG PRINT */ if (verbose) { printf("[ DEBUG ] traversing unpacked matrix\n"); }
    GrB_Index Ap_index;
    GrB_Index Ai_index;
    GrB_Index Ai_index_start, Ai_index_end;
    int current_color, neighbor_color;
    for (Ap_index = 0; Ap_index < Ap_size - 1; Ap_index++) {
        
        Ai_index_start = Ap[Ap_index];
        Ai_index_end = Ap[Ap_index + 1];

        GrB_Vector_extractElement(&current_color, C, Ap_index);

        for (Ai_index = Ai_index_start; Ai_index < Ai_index_end; Ai_index++) {
            GrB_Vector_extractElement(&neighbor_color, C, Ai[Ai_index]);

            /* DEBUG PRINT */ if (verbose) { printf("[ DEBUG ] current_color: %d, neighbor_color: %d\n", current_color, neighbor_color); }

            TEST_ASSERT(neighbor_color != current_color);
        }
    }

    /* DEBUG PRINT */ if (verbose) { printf("[ DEBUG ] finished traversing unpacked matrix\n"); }



    /* get number of colors - aka get max of vector C */
    int num_colors = 0;
    GrB_reduce(&num_colors, GrB_NULL, GrB_MAX_MONOID_INT32, C, NULL);
    printf("Number of Colors: %d\n", num_colors);



    /* clean up (don't understand this) */
    OK(LAGraph_Delete(&G, msg));
    LAGraph_Finalize(msg);
    LAGraph_Random_Finalize(msg);
}


/* used to run tests
*  always end on {NULL, NULL}
*/
TEST_LIST =
{
    {"coloring_simple", test_coloring_simple},
    {NULL, NULL}
};
