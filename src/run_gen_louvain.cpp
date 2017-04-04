#include "stability_gen.h"
#include "louvain_gen.h"
#include "vector_partition.h"

#include <lemon/smart_graph.h>
#include <vector>

#include <Rcpp.h>
using namespace Rcpp;

template<typename G, typename E>
bool read_edgelist_weighted_from_data(NumericMatrix data, G &graph, E &weights) {
    // function to read in graph data from numerical [from, to, weight] list 'data'
    // and assemble a (templated) lemon graph with topology stored in 'graph' and Edgemap
    // with weights stored in 'weights'.

    // define Node class for convenience
    typedef typename G::Node Node;

    // node index counted
    int max_node_id_seen =-1;
    // largest dimension of the list == number of rows
    int num_l_dim = data.nrow();

    // loop over complete list
    for (int i = 0; i < num_l_dim; ++i) {

        // get node ids --- should be first and second colum of list
        int node1_id = data(i,0);
        int node2_id = data(i,1);

        // if we have not seen a node with such a high id before add all nodes in between
        // to keep id counts contiguous and in 1:1 correspondence within the lemon graph
        if (node1_id > max_node_id_seen) {
            int difference = node1_id - max_node_id_seen;
            for (int i=0; i<difference; ++i) {
                graph.addNode();
            }
            max_node_id_seen = node1_id;
        }
        if (node2_id > max_node_id_seen) {
            int difference = node2_id - max_node_id_seen;
            for (int i=0; i<difference; ++i) {
                graph.addNode();
            }
            max_node_id_seen = node2_id;
        }

        // read in list is undirected, edges should only be created once
        // if (node1_id > node2_id) {
        //     std::cout << "you provided a double sided list";
        //     std::cout << "this is not necessary -- some entries will be ignored";
        // }

        // read in weight and assign to weight map
        double weight = data(i,2);
        typename G::Edge edge = graph.addEdge(graph.nodeFromId(node1_id),
                                              graph.nodeFromId(node2_id));
        weights.set(edge, weight);
    }

    return true;
}

std::vector<std::vector<double> > read_null_model(NumericMatrix nullmod, int num_nodes) {
    // function to read in the null model vectors

    // read in number of nm vectors. Should be a multiple of 2
    int num_null_model_vectors = nullmod.ncol();
    if (num_null_model_vectors % 2 != 0){
        //TODO: return some kind of error!
    }

    unsigned int num_null_vec_entries = num_nodes;
    std::vector<std::vector<double> > null_vectors(num_null_model_vectors,std::vector<double>(num_null_vec_entries,0));


    for (int i= 0; i < num_null_model_vectors; i++){
        for(unsigned int j = 0; j < num_null_vec_entries; j++){
            // data should be read in from nullmod..
            // TODO: double check that this ordering is correct..
            null_vectors[i][j] = nullmod(j,i);

        } 
    }
    return null_vectors;
}

// Generic Louvain optimisation
// input: [graph] list of edges and weights; [nullmod] null model, stationary distribution; [time] markov time, resolution parameter; [M] number of iterations.
// output: vector of length M with stability values; matrix n x M, where n is number of nodes, containing M partition vectors.
//
// TODO: create specialised version when no null model is given / can be extracted from matrix!?

// [[Rcpp::export]]
List run_gen_louvain(NumericMatrix graph, NumericMatrix nullmod, double time, int M) {

	// initialise markov time from input
    double current_markov_time = time;

    // initialise graph and edge map
    lemon::SmartGraph input_graph;
    lemon::SmartGraph::EdgeMap<double> input_graph_weights(input_graph);

    // read in data of graph
    read_edgelist_weighted_from_data(graph, input_graph, input_graph_weights);
    int num_nodes = lemon::countNodes(input_graph);

    // read in data of null_model vectors
    std::vector<std::vector<double> > null_model = read_null_model(nullmod, num_nodes);
	
	// initialise output
	NumericVector out_stab(M);
	NumericMatrix out_mat(num_nodes, M);

	for (int i=0; i<M; ++i) {
		
		// set seed
		int seed = R::runif(-1e9,+1e9);
		srand(seed);
		
	    // define type for vector partition,initialise start partition, quality functions etc.
	    typedef clq::VectorPartition partition;
	    partition start_partition(num_nodes);
	    start_partition.initialise_as_singletons();

	    clq::find_linearised_generic_stability quality(current_markov_time);
	    clq::linearised_generic_stability_gain quality_gain(current_markov_time);
	    std::vector<partition> optimal_partitions;

	    // call Louvain
	    double stability = clq::find_optimal_partition_louvain_gen<partition>(
	                           input_graph, input_graph_weights, null_model, quality, quality_gain,
	                           start_partition, optimal_partitions, 1e-12);
    
	    // get optimal partition and store it in an R compatible data structure
	    partition best_partition = optimal_partitions.back();
	    NumericVector output_partition(num_nodes);
	    for (int node = 0; node < num_nodes; ++node) {
	        output_partition[node] = int(best_partition.find_set(node));
	    }
	
		// write the output
		out_stab(i) = stability;
		out_mat(_,i) = output_partition;
	}

    // output list
    List ReturnArg;
	ReturnArg["stability"] = out_stab;
	ReturnArg["partition"] = out_mat;
    return ReturnArg;
}
