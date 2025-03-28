//--------------------------------------------------------------------------
// error catching macros
//--------------------------------------------------------------------------
// free workspace (runs on success)
#define LG_FREE_WORK                \
{                                   \
    GrB_free (&weights) ;           \
    GrB_free (&max_neighbor_weights); \
    GrB_free (&empty) ;             \
    GrB_free (&candidates) ;        \
    GrB_free (&independent_set) ;   \
    GrB_free (&independent_set_neighbors) ; \
    GrB_free (&independent_set_neighbors_colors) ; \
    GrB_free (&MIS_candidates) ;    \
}
// free everything (runs on error)
#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK ;                  \
    GrB_free (&JP_coloring_copy);   \
}


//--------------------------------------------------------------------------
// error codes
//--------------------------------------------------------------------------
#define JP_COLORING_NON_SYMMETRIC -5501
#define JP_COLORING_STALLED       -5502

#define NUM_ALLOWED_STALLS 32

//--------------------------------------------------------------------------
// includes
//--------------------------------------------------------------------------
#include "LG_internal.h"    // header for internal lagraph
#include "LAGraphX.h"       // header for experimental lagraph


//--------------------------------------------------------------------------
// algorithm
//--------------------------------------------------------------------------
int LAGraph_coloring_JP
(
    // outputs
    GrB_Vector *JP_coloring,    // vector of natural numbers representing
                                // the colors of each node
    int *JP_num_colors,         // number of unique colors
    
    // inputs
    LAGraph_Graph G,            // input graph
    uint64_t seed,              // random number seed
    char *msg                   // LAGraph error msg

)
{   
    // clear error msg
    LG_CLEAR_MSG ;
    
    //--------------------------------------------------------------------------
    // define objects
    //--------------------------------------------------------------------------
    GrB_Vector JP_coloring_copy = NULL;             // local version of JP_coloring
    int JP_num_colors_copy = 0;                     // local version of JP_num_colors

    GrB_Vector candidates = NULL;                   // candidates    
    GrB_Vector MIS_candidates = NULL;               // MIS candidates
    GrB_Vector independent_set = NULL;              // independent set
    GrB_Vector weights = NULL;                      // random weights
    GrB_Vector max_neighbor_weights = NULL;         // maximum random weight of neighbors
    GrB_Vector independent_set_neighbors = NULL;    // neighbors of independent set
    GrB_Matrix independent_set_neighbors_colors = NULL; // colors of neighbors

    GrB_Vector empty = NULL;                        // empty vector for sparsification
    GrB_Matrix A = NULL;                            // adjacency matrix of G
    GrB_Index n;                                    // number of nodes

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;                      // check graph G
    LG_ASSERT (JP_coloring != NULL, GrB_NULL_POINTER) ;         // check output
    LG_ASSERT (JP_num_colors != NULL, GrB_NULL_POINTER) ;       // check output
    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED ||              // check symmetric
        (G->kind == LAGraph_ADJACENCY_DIRECTED &&
         G->is_symmetric_structure == LAGraph_TRUE))
    {
        A = G->A ;                                              // symmetric
    }
    else {
        LG_ASSERT_MSG (false, JP_COLORING_NON_SYMMETRIC, "G->A must be symmetric") ;
    }


    //--------------------------------------------------------------------------
    // initialize objects
    //--------------------------------------------------------------------------
    GRB_TRY (GrB_Matrix_nrows (&n, A)) ;
    GRB_TRY (GrB_Vector_new (&weights, GrB_UINT64, n)) ;
    GRB_TRY (GrB_assign (weights, NULL, NULL, 0, GrB_ALL, n, NULL)) ;
    LG_TRY (LAGraph_Random_Seed(weights, seed, msg)) ;
    GRB_TRY (GrB_Vector_new (&max_neighbor_weights, GrB_UINT64, n)) ;
    GRB_TRY (GrB_Vector_new (&empty, GrB_BOOL, n)) ;
    GRB_TRY (GrB_Vector_new (&independent_set, GrB_BOOL, n)) ;
    GRB_TRY (GrB_Vector_new (&independent_set_neighbors, GrB_BOOL, n)) ;
    GRB_TRY (GrB_Matrix_new (&independent_set_neighbors_colors, GrB_UINT64, n, n))
    GRB_TRY (GrB_Vector_new (&candidates, GrB_BOOL, n)) ;
    GRB_TRY (GrB_assign (candidates, GrB_NULL, GrB_NULL, true, GrB_ALL, n, GrB_NULL)) ;
    GRB_TRY (GrB_Vector_new (&MIS_candidates, GrB_BOOL, n)) ; 

    GRB_TRY (GrB_Vector_new (&JP_coloring_copy, GrB_UINT64, n)) ;   

    //--------------------------------------------------------------------------
    // optional - handle singletons and ignore_node
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // main algorithm
    //--------------------------------------------------------------------------
    GrB_Index num_candidates = 0;
    GRB_TRY (GrB_Vector_nvals (&num_candidates, candidates)) ;
    int64_t curr_color = 0;

    printf("--------- A \n");
    LAGraph_Matrix_Print(A, 2, stdout, msg);
    printf("--------- weights \n");
    LAGraph_Vector_Print(weights, 2, stdout, msg);

    while (num_candidates > 0) {
        // STEP 0: copy candidates to MIS_candidates
        // assign: copy + replace old values
        // FIXME: check if vector_dup is faster
        GRB_TRY (GrB_assign (MIS_candidates, GrB_NULL, GrB_NULL, candidates,
            GrB_ALL, n, GrB_DESC_R)) ;

        GrB_Index num_MIS_candidates = 0;        
        GRB_TRY (GrB_Vector_nvals (&num_MIS_candidates, MIS_candidates)) ;
        GrB_Index last_num_MIS_candidates = num_MIS_candidates;
        GrB_Index num_stalls = 0;

        // prints: check if everything is good        
        printf("--------- candidates \n");
        LAGraph_Vector_Print(candidates, 2, stdout, msg);

        while (num_MIS_candidates > 0) {
            
            // STEP 1: pick winners
            // mxv: find the maximum nearby weight and add to max_neighbor_weights
            //      - max_neighbor_weight will only calculate for candidates
            //      - colored neighbors not considered - weight will be empty
            //      - DESC_RS - delete old values (replace) + structural mask
            // eWiseAdd: winners are added to independent_set
            //           - add means weights without max_neighbor_weights fall thru
            //             this is expected - these nodes have no neighbors to check
            //           - DESC_R - delete old values
            // select: keep 1s and delete 0s -> turns into structural vector
            //
            // FIXME: add push vs pull        
            GRB_TRY(GrB_mxv(max_neighbor_weights, MIS_candidates, GrB_NULL,                 
                GrB_MAX_SECOND_SEMIRING_UINT64, A, weights, GrB_DESC_RS));
            GRB_TRY(GrB_eWiseAdd(independent_set, MIS_candidates, GrB_LOR,
                GrB_GT_UINT64, weights, max_neighbor_weights, GrB_DESC_S));
            GRB_TRY(GrB_select(independent_set, GrB_NULL, GrB_NULL, 
                GrB_VALUEEQ_BOOL, independent_set, true, GrB_NULL));

            
            // STEP 2: cleanup candidates, check quit condition + stall
            // assign: remove independent_set from candidates
            // mxv: find neighbors of independent set + replace vector
            // assign: remove independent_set from candidates
            // assign: remove weights for non-candidates for step 1
            // FIXME: add push vs pull            
            GRB_TRY (GrB_assign (MIS_candidates, independent_set, GrB_NULL, empty,
                GrB_ALL, n, GrB_DESC_S)) ;
            GRB_TRY (GrB_mxv (independent_set_neighbors, MIS_candidates, GrB_NULL,
                LAGraph_any_one_bool, A, independent_set, GrB_DESC_RS)) ;
            GRB_TRY (GrB_assign (MIS_candidates, independent_set_neighbors, GrB_NULL,
                empty, GrB_ALL, n, GrB_DESC_S)) ;

            // print matrices
            printf("--------- independent_set \n");
            LAGraph_Vector_Print(independent_set, 2, stdout, msg) ;
            printf("--------- MIS_candidates \n");
            LAGraph_Vector_Print(MIS_candidates, 2, stdout, msg) ;            
            
            
            // STEP 3: check quit condition
            // nvals + if: break if no more candidates
            // if: check if candidates is same as last iteration
            //     if so, increment stall count, break if too many
            //     and redo weights
            // save last num candidates
            GRB_TRY (GrB_Vector_nvals (&num_MIS_candidates, MIS_candidates)) ;
            if (num_MIS_candidates == 0) { break ; }
            if (num_MIS_candidates == last_num_MIS_candidates) {
                num_stalls++ ;
                LG_ASSERT_MSG (num_stalls <= NUM_ALLOWED_STALLS, JP_COLORING_STALLED, "MIS stalled") ;
                LG_TRY (LAGraph_Random_Next (weights, msg)) ;
            }
            last_num_MIS_candidates = num_MIS_candidates;

        }

        // run JP
        // at this point, independent_set is now maximal
        // mxm: find all colors of independent_set neighbors
        // convert into bitmap
        // GrB_mxv(indpendent_set_neighbors_colors, independent_set, GrB_NULL,
        //     GrB_SECOND_SEMIRING_UINT64, A, JP_coloring_copy, GrB_DESC_RS) ;
        
        // alternate approach: MIS
        // color independnet set
        GRB_TRY(GrB_assign(JP_coloring_copy, independent_set, GrB_NULL, curr_color, GrB_ALL, n, GrB_DESC_S)) ;
        curr_color ++;

        // prepare for next iteration
        // remove independent_set from candidates
        // reset independent_set
        GRB_TRY(GrB_assign(candidates, independent_set, GrB_NULL, empty, GrB_ALL, n, GrB_DESC_S)) ;
        GRB_TRY (GrB_Vector_nvals (&num_candidates, candidates)) ;
        GRB_TRY (GrB_assign (independent_set, GrB_NULL, GrB_NULL, empty, GrB_ALL, n, GrB_NULL)) ;
    }
    JP_num_colors_copy = curr_color - 1;
    
    //--------------------------------------------------------------------------
    // clean up and return outputs
    //--------------------------------------------------------------------------
    GRB_TRY (GrB_wait (JP_coloring_copy, GrB_MATERIALIZE)) ;
    (*JP_coloring) = JP_coloring_copy ;
    (*JP_num_colors) = JP_num_colors_copy ;
    LG_FREE_WORK ;
    return (GrB_SUCCESS) ;
}