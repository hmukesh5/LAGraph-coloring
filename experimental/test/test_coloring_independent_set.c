// no idea about these 2 includes
#include <stdio.h>
#include <acutest.h>

// includes from LAGraph, also dunno
#include <LG_internal.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <LG_Xtest.h>
#include <LG_test.h>


char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;

#define LEN 512
char filename[LEN + 1];

const char* matrix_files[] = {
    "ldbc-undirected-example-unweighted.mtx",
}

void test_coloring_independent_set(void)
{
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
    GxB_set (GxB_BURBLE, false) ;

    int num_colors = 0;
    double time = LAGraph_WallClockTime();    
    LAGraph_coloring_independent_set_optimized(&C, &num_colors, G, msg);
    time = LAGraph_WallClockTime() - time;

    GxB_set (GxB_BURBLE, false) ;

    printf("\nTook %g seconds\n", time);
    printf("Final color vector:\n"); LAGraph_Vector_Print(C, LAGraph_SHORT, stdout, msg); }

    
    // ------------------------------------------------
    // check if coloring is valid
    // ------------------------------------------------

    OK (LG_check_coloring(G->A, C, msg));


    printf("Number of Colors: %d\n", num_colors);


    /* clean up (don't understand this) */
    OK(LAGraph_Delete(&G, msg));
    LAGraph_Finalize(msg);
    LAGraph_Random_Finalize(msg);
}

TEST_LIST =
{
    {"coloring_independent_set", test_coloring_independent_set},
    {NULL, NULL}
};
