
#include <sys/time.h>
#include <sys/resource.h>
#include <cassert>
#include <valarray>
#include <fstream>
#include <iostream>

//////////////////////REZA////////////////////////////////////////>>>>>>>>>>>>>>>>>>>>
//////////////////////////////////////////////////////////////////

#include <regex>
#include <iterator>
#include <set>
#include <string>

//************ADDED FOR SPR neihborhood************>>>>>>>>>>>>>>>>>>>>>>>>>>>
#include <cstdio>
#include <cstdlib>
#include <ctime>
//#include <string>
#include <cstring>
//#include <iostream>
//#include <fstream>
#include <sstream>
#include <climits>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <list>
#include <time.h>
#include "rspr.h"

#include "Forest.h"
#include "ClusterForest.h"
#include "LCA.h"
#include "ClusterInstance.h"
#include "UndoMachine.h"
#include "lgt.h"
#include "sparse_counts.h"
#include "node_glom.h"
#include <regex>


void split(const string &s, char delim, vector<string> &elems);
vector<string> split(const string &s, char delim);
set<string> find_non_common_taxa_set(const string &supertree, const string &source_tree);
string change_branch_lengths(const string &tree, int percentage_to_be_reweighted, int new_weight);
void write_line_to_File(string s, const char* file_name);
double calculate_rf_btwn_ST_n_source_tree(int argc, char** argv);
double calculate_total_rf(string source_trees[], set<string> non_shared_taxa[], string supertree, bool is_weighted);



void writeToFile(string s, const char* file_name);
void adjustTree(Node* tree);
int produce_all_spr_neighbors(Node* tree, int number_of_taxa);
int produce_n_percent_of_spr_neighbors(Node* myTree, int number_of_taxa, int percentage);
int number_of_taxa(string const tree);
void total_number_of_nodes(Node* node, int& total_nodes);

//////////////////Added for restrict_st()/////////////////////////////////////////
void split(const string &s, char delim, vector<string> &elems);
vector<string> split(const string &s, char delim);
set<string> find_non_common_taxa_set(const string &supertree, const string &source_tree);
void restrict_supertree(Node& supertree, set<string>& non_shared_taxon_set);
void find_node_to_be_removed(Node& tree, string taxon, Node* & taxan_to_be_removed);
std::string ReplaceAll(std::string str, const std::string& from, const std::string& to);
static unsigned int NUMBER_OF_SOURCE_TREES_ZZ;
static unsigned int NUMBER_OF_TAXA;

////////////////RFS implementation//////////////////////
#include <unordered_map>
void initialize_alpha_beta(Node* Q);
void set_cluster_size(Node* Q);
void find_int_labels_for_leaves_in_supertree(Node* root, int& current_label, unordered_map<string, int>& map);
void set_int_labels_for_leaves_in_source_tree(Node* root, unordered_map<string, int>& map);
void compute_lca_mapping_S_R(Node& S, Node& R);
void compute_lca_mapping_helper_1(vector<int>& cluster, Node* R, int& lca);
void compute_lca_mapping_helper_2(vector<int>& cluster, Node* R, Node*& lca, bool& found);
void set_cluster_and_cluster_size(Node* T);
void set_cluster_size_in_supertree(Node* T);
void reset_fields_to_initial_values(Node* T);
Node* apply_SPR_RS_algorithm_to_find_best_regraft_place(Node& T, Node& spr_on, Node* source_trees_array[], set<string> non_shared_taxon_set[], bool weighted);
void traverse_S_and_update_alpha_beta_in_Q(Node& T, Node& Q_in_restricted_st, Node& v_in_restricted_st, Node& S, Node& S_prime, bool weighted);
void find_b_in_lemma12(Node& S, Node& R); // i'm not using this now
void suppress_nodes_with_mapping_in_Rv(Node& S_prime, Node& v);
void find_best_regraft_place(Node& n, Node*& best_regraft_place, int& max);
void preorder_traversal(Node& n);
int find_best_node_to_prune_and_its_best_regraft_place(Node& T, Node* source_trees_array[], set<string> non_shared_taxon_set[], Node* & best_node_to_prune, Node* & best_node_to_regraft, bool weighted);
void find_F_T(Node& S, Node& T, int& num_clusters_not_in_T_prime);
void find_weighted_rf_dist(Node& S, Node& T_prime, int& dist);
void put_internal_nodes_in_vector(Node& n, vector<Node*>& nodes);
void put_all_nodes_in_vector(Node& n, vector<Node*>& nodes);
void find_f_u_with_regard_to_R(vector<int>& cluster, vector<int>& source_tree_leaf_set, Node * R, Node*& lca, bool & found, int& f_u);



void print_weighted_tree(Node& n);

template <typename T>
bool IsSubset(std::vector<T> A, std::vector<T> B)
{
	//std::sort(A.begin(), A.end());
	//std::sort(B.begin(), B.end());
	return std::includes(A.begin(), A.end(), B.begin(), B.end());
}


//////////////////////from ratchet:
vector<string> split(const string &s, char delim);
void split(const string &s, char delim, vector<string> &elems);
set<string> find_non_common_taxa_set(const string &supertree, const string &source_tree);




/////////////////////////////     main()     /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
int NUM_SPR_NGHBRS = 0;

int main(int argc, char** argv) {
	srand((unsigned)time(NULL));

	cout << "Please enter the pre-specified number of 'ratchet' iterations: " << endl;
	//pre-specified number of iterations if it didn't stop after this number of iterations
	int number_of_ratchet_iterations;
	cin >> number_of_ratchet_iterations;

	cout << "Please enter the percentage of clades to be re-weighted in Ratchet search (between 0 and 100): " << endl;
	int percentage_of_clades_to_be_reweighted ;
	cin >> percentage_of_clades_to_be_reweighted;

	cout << "Please enter the new weight to which clades be re-weighted in Ratchet search (between 0 and 10): " << endl;
	int ratchet_weight;
	cin >> ratchet_weight;

	cout << "********************************************************************" << endl;


	string init_supertree;
	ifstream init_sup("initial_supertree");
	if (init_sup.good())
	{
		string sLine;
		getline(init_sup, sLine);
		init_supertree = sLine;
	} else {
		cout << "Initial supertree file not found!!" << endl;
		return -1;
	}

	int number_of_source_trees = 0;
	ifstream myfile1(argv[1]);
	string l0;
	while (std::getline(myfile1, l0))
	{
		number_of_source_trees ++;
	}
	myfile1.close();

	::NUMBER_OF_SOURCE_TREES_ZZ = number_of_source_trees;

	string source_trees_newick[number_of_source_trees];
	Node* source_trees_root[number_of_source_trees];
	int cntr0 = 0;
	ifstream myfile2(argv[1]);
	string l1;  //hopefully size of trees are not larger than str.max_size()
	//cout << "Source Trees:\n" ;
	while (std::getline(myfile2, l1))
	{
		source_trees_newick[cntr0] = l1;
		source_trees_root[cntr0] = build_tree(l1);
		adjustTree(source_trees_root[cntr0]);
		//cout << source_trees_root[cntr0]->str_subtree() << endl;
		cntr0 ++;
	}
	myfile2.close();


	//non-shared taxa between each source tree and any supertree is the same in any iteration.
	//thus, I calculate it once here and just use it.
	set<string> non_shared_taxa_arr[number_of_source_trees];
	int cntr2 = 0;
	for (string t : source_trees_newick) {
		non_shared_taxa_arr[cntr2] = find_non_common_taxa_set(init_supertree, t);
		cntr2 ++;
	}


	///////////setting int labels for leaves, i.e. mapping names to integers from 1 to n
	Node* supertree = build_tree(init_supertree);
	adjustTree(supertree);
	//cout << "\n\nInitial Supertree:\n" ;
	//cout << supertree->str_subtree() << endl;
	set_cluster_and_cluster_size(supertree);

	unordered_map<string, int> int_label_map;
	int starting_label = 1;
	find_int_labels_for_leaves_in_supertree(supertree, starting_label, int_label_map);  //finding int labels for all taxa in supertree
	::NUMBER_OF_TAXA = starting_label;
	// set int labels for taxa in source trees
	for (int i = 0; i < number_of_source_trees; i++) {
		set_int_labels_for_leaves_in_source_tree(source_trees_root[i], int_label_map);
		//also, since source trees won't change, for each node in each source tree, lets find clusters only ones, and store it in DS
		set_cluster_and_cluster_size(source_trees_root[i]);
	}




	cout << "\nsource trees and supertree has been read, let's start running edge-ratchet algorithm using RFS" << endl;
	cout << "********************************************************************" << endl;



	//running time
	clock_t start_time_total, finish_time_total;
	start_time_total = clock();


	//just to keep track of best ST seen throughout the algorithm
	string the_best_supertree_seen = init_supertree;
	int the_best_rf_distance_seen = INT_MAX;


	////////////////////////////////////////////////
	//////////////////////ratchet///////////////////
	////////////////////////////////////////////////

	//ratchet search loop
	for (int ratchet_counter = 1; ratchet_counter < number_of_ratchet_iterations + 1; ratchet_counter++) {

		bool ratchet = true; //when true, use re-weighted input trees; when false, use unweighted trees.
		if (ratchet) {
			for (int i = 0; i < number_of_source_trees; ++i) {
				source_trees_root[i]->reweight_edges_in_source_tree(percentage_of_clades_to_be_reweighted, ratchet_weight);
				//cout << source_trees_root[i]->str_subtree() << endl;
			}
			//print_weighted_tree(*source_trees_root[0]);
			//cout <<  "\n\n";
		}


		int best_score_of_current_hill = INT_MAX;
		string best_supertree_of_current_hill;

		//search until reaching a local optimum
		int iteration = 0;
		while (true) {
			NUM_SPR_NGHBRS = 0;
			clock_t start_time_iter, finish_time_iter;
			start_time_iter = clock();

			iteration ++;

			Node* best_node_to_prune;
			Node* best_node_to_regraft;
			//cout << "\ninit  ST: " << supertree->str_subtree() << "\n";
			int current_score = find_best_node_to_prune_and_its_best_regraft_place(*supertree, source_trees_root, non_shared_taxa_arr, best_node_to_prune, best_node_to_regraft, ratchet);
			//cout << "T: " << supertree->str_subtree() << endl;
			//cout << "best node to  prune: " << best_node_to_prune->str_subtree() << endl;
			//cout << "best node to regrft: " << best_node_to_regraft->str_subtree() << endl;

			if (current_score < best_score_of_current_hill) {

				int which_sibling = 0;
				Node* old_sibling = best_node_to_prune->spr(best_node_to_regraft, which_sibling);
				set_cluster_and_cluster_size(supertree);


				supertree = best_node_to_prune->find_root();
				adjustTree(supertree);

				best_score_of_current_hill = current_score;
				best_supertree_of_current_hill = supertree->str_subtree();

				finish_time_iter = clock();
				float diff ((float)finish_time_iter - (float)start_time_iter);
				float iter_time = diff / CLOCKS_PER_SEC;
				cout << "Iter: " << iteration << ", num_spr_neighbours: " <<
				     NUM_SPR_NGHBRS << ", RF_dist: " << current_score << ", time(sec) : " << iter_time << "\n";
				//cout << "\n" << supertree->str_subtree() << "\n\n\n";

			}
			else { // local optimum

				if (ratchet) {
					cout << "\nratchet local opt (re-weighted) reached." << endl;
					cout << "---------------------------------------------------\n" << endl;
				} else {
					cout << "\nregular local opt (original weights) reached. ###" << endl;
				}

				//cout << "We have reached a local optimum which has better \(or equal\) WRF distance than all its spr-neighbors\n";
				//cout << "The best SuperTree found after " << iteration << " number of iterations is: " "\n";
				//cout << "And its RF distance is " << best_score_of_current_hill << "\n";
				//cout << "--------------------------------------------------------" << endl;



				if (!ratchet) { //we are at the end of one ratchet iteration
					if (best_score_of_current_hill < the_best_rf_distance_seen) { //keep track of best supertree seen so far
						cout << "##############Whoooooooop!! Better supertree seen####################" << endl;
						the_best_rf_distance_seen = best_score_of_current_hill;
						the_best_supertree_seen = best_supertree_of_current_hill;
					}

					cout << "=======================================end of " << ratchet_counter << "-th ratchet iter=========================================" << endl;
					cout << "========================================================================================================" << endl;

					break;  //end of second phase of ONE ratchet iteration

				} else { //now perform a regular branch swapping

					ratchet = false;
					iteration = 0;
				}

				
				best_score_of_current_hill = INT_MAX;	//this line so necessary!!
				//Because we are about to start a completely new hill climbing with new objective function.
				//without this line, after 2nd ratchet iter, no improvement will be made since the weighted
				//distance is always larger than un-weighted and the init ST from previous step is local opt
				//already.
			}

		}//end of one branch swapping search loop

	}//end of ratchet search loop


	finish_time_total = clock();
	float diff ((float)finish_time_total - (float)start_time_total);
	float seconds = diff / CLOCKS_PER_SEC;
	//cout << "The #of clock ticks of this iteration: " << diff << endl;
	//cout << "\n" << source_trees_root[0]->str_subtree() << endl;
	//cout << "init_supertree:\n" << init_supertree << endl;
	//cout << "The best tree found:\n" << the_best_supertree_seen << endl;

	cout << "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
	cout << "The best SuperTree found after " << number_of_ratchet_iterations << " number of ratchet iterations is: " << endl;
	cout << "And its RF distance is " << the_best_rf_distance_seen << endl;
	cout << "the running time is: " << seconds << " sec." << "\n\n";
	cout << the_best_supertree_seen << endl;
	cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;



	//don't forget to free memory!!
	supertree->delete_tree();
	for (int i = 0; i < number_of_source_trees; ++i) {
		source_trees_root[i]->delete_tree();
	}

	return 0;
}


//tetsting
void preorder_traversal(Node& n) {
	cout << n.str_subtree() << endl;
	cout << " prenum is : " << n.get_preorder_number() << endl;
	//cout << "lca mapping node is: " << n.get_lca_mapping().str_subtree() << "\n" << endl;
	//cout << n.get_alpha() << "-" << n.get_beta() << "= " << n.get_alpha() - n.get_beta() << "\n" << endl;
	//cout << n.get_edge_weight() << "\n" << endl;

	//just for testing:
	list<Node *>::iterator c;
	list<Node *> children = n.get_children();
	for (c = children.begin(); c != children.end(); c++) {
		preorder_traversal(**c);
	}
}

//tetsting
void print_weighted_tree(Node& n) {
	//cout << n.str_subtree() << endl;
	//cout << " prenum is : " << n.get_preorder_number() << endl;
	//cout << "lca mapping node is: " << n.get_lca_mapping().str_subtree() << "\n" << endl;
	//cout << n.get_alpha() << "-" << n.get_beta() << "= " << n.get_alpha() - n.get_beta() << "\n" << endl;
	//cout << n.get_edge_weight() << "\n" << endl;
	cout << n.str_subtree() << ": " << n.get_edge_weight() << "  ,";

	//just for testing:
	list<Node *>::iterator c;
	list<Node *> children = n.get_children();
	for (c = children.begin(); c != children.end(); c++) {
		print_weighted_tree(**c);
	}
}


//For each internal node of T to be pruned, find best regraft place
//alongside, keep track of best SPR neighbor seen (keep track of best node to be pruned and its best regraft place):
//the neighbour with smaller |F_T| is better:
//because RF_dist(S,T)=|I(T)| + |I(S)| - 2|F_T|
//returns min_F found
int find_best_node_to_prune_and_its_best_regraft_place(Node& T, Node* source_trees_array[], set<string> non_shared_taxon_arr[], Node* & best_node_to_prune, Node* & best_node_to_regraft, bool weighted) {

	int min_F = INT_MAX;
	//cout << "T: " << T.str_subtree() << endl;
	vector<Node*> internal_nodes;
	put_all_nodes_in_vector(T, internal_nodes);
	vector<Node*>::iterator iter, end;
	//any internal node except root may be pruned, iter plays the role of "v" in the algorithm
	for (iter = internal_nodes.begin(), end = internal_nodes.end(); iter != end;  iter++) { //for each internal node v in T to be pruned:

		if (!(*iter)->get_p()) { //root
			continue;
		}

		//cout << "---------------------------------" << endl;

		Node* best_regraft_place_for_current_v = apply_SPR_RS_algorithm_to_find_best_regraft_place(T, **iter, source_trees_array, non_shared_taxon_arr, weighted);

		//cout << "\n------T before: " << T.str_subtree() << endl;
		//cout << "-----spr_on(v): " << (*iter)->str_subtree() << endl;
		//cout << "--best regraft: " << best_regraft_place_for_current_v->str_subtree() << endl;

		int which_sibling = 0;
		Node* old_sibling = (*iter)->spr(best_regraft_place_for_current_v, which_sibling);
		adjustTree(&T);
		//cout << "------T after : \n" << T.str_subtree() << endl;
		//cout << "\n-----------------------------------------------------\n" ;


		//Now that we know for "v" to be pruned, what's the best regraft place,
		//we have a best SPR neighbour, T'. We need to know how good is this neighbour,
		//to ompare it to other best neighbours for different v's. Note, T is now T' (after spr move)
		//Thus, we have to restrict it again.
		vector<Node*> restricted_T_primes;
		for (int i = 0; i <::NUMBER_OF_SOURCE_TREES_ZZ; i++) {
			Node* current_sup_tree = build_tree(T.str_subtree()); //note spr() has been applied to T, and now T is actually T'
			adjustTree(current_sup_tree);
			current_sup_tree->copy_fields_for_supertree(T);	//DON"T forget this!!
			restrict_supertree(*current_sup_tree, non_shared_taxon_arr[i]);
			restricted_T_primes.push_back(current_sup_tree);
			//cout << "\n" << current_sup_tree->str_subtree() << endl;
		}

		int current_F = 0 ; //which (for unweighted case) is "num_clusters_not_in_T_prime", i.e. RF distance
		if (weighted) {
			for (int i = 0; i <::NUMBER_OF_SOURCE_TREES_ZZ; i++) {
				find_weighted_rf_dist(*source_trees_array[i], *restricted_T_primes[i], current_F);
			}
		} else {
			for (int i = 0; i <::NUMBER_OF_SOURCE_TREES_ZZ; i++) {
				find_F_T(*source_trees_array[i], *restricted_T_primes[i], current_F);
			}
		}

		//cout << "best F: " << min_F << ", and curent F: " << current_F << endl;

		if (current_F < min_F) {
			//cout << "\n\nbest F was : " << min_F << ">>>>>>>>>>>>>>>>>>>>>>>>>>> better SPR move with F: " << current_F << endl;
			//cout << "------T after : \n" << T.str_subtree() << endl;
			//cout << "\n" << T.str_subtree() << "\n\n";
			//cout << "--best regraft: " << best_regraft_place_for_current_v->str_subtree() <<  endl;

			min_F = current_F;
			best_node_to_prune = *iter;
			best_node_to_regraft = best_regraft_place_for_current_v;
		}
		//restore tree to consider next node to be pruned
		if (best_regraft_place_for_current_v == old_sibling) {
			cout << "-------------------OOOOOOOOOOOOOOOOOOOpsssssssss\n";
		}
		//else {
		(*iter)->spr(old_sibling, which_sibling);
		//reset_alpha_beta(T);
		adjustTree(&T);
		//}


		//prevent mem leak
		for (int i = 0; i <::NUMBER_OF_SOURCE_TREES_ZZ; i++) {
			restricted_T_primes[i]->delete_tree();
		}
	}

	return min_F;

}


//preorderly
void put_internal_nodes_in_vector(Node& n, vector<Node*>& nodes) {
	if (n.is_leaf()) {
		//not internal node
	} else {
		nodes.push_back(&n);

		list<Node *>::iterator c;
		list<Node *> children = n.get_children();
		for (c = children.begin(); c != children.end(); c++) {
			put_internal_nodes_in_vector(**c, nodes);
		}
	}
}


//preorderly
void put_all_nodes_in_vector(Node& n, vector<Node*>& nodes) {

	nodes.push_back(&n);

	list<Node *>::iterator c;
	list<Node *> children = n.get_children();
	for (c = children.begin(); c != children.end(); c++) {
		put_all_nodes_in_vector(**c, nodes);
	}

}


//Actually finds RF distance: finds the number of internal nodes whose corresponding lca does not exist in T
//Note the way I implemented find_best_node_to_prune_and_its_best_regraft_place(), I don't want to change T so that I have to do calculations (lca_mapping, ...) for each v to be pruned
//The node T_prime being passed to this function is best neghibour for current v to be pruned.
//THUS the lca_mappings are NOT valid anymore, and that's why I compute it again here
void find_F_T(Node& S, Node& T_prime, int& num_clusters_not_in_T_prime) {
	//if S curresponds to a trivial bipartition, i.e. bipartition with one leaf on one side, DO NOT count it. There are 2 cases:
	//1- it is a leaf
	//2- it is child of root, and it's sibling is a leaf
	//Also we don't count bipartition curresponding to roo ROOT
	if (S.is_leaf()) { //leaf
		return;
	}
	else {
		//preorderly
		list<Node *>::iterator c;
		list<Node *> children = S.get_children();
		for (c = children.begin(); c != children.end(); c++) {
			find_F_T(**c, T_prime, num_clusters_not_in_T_prime);
		}

		if (S.get_p() == NULL) { //root
			//do not count
		} else if (S.get_p()->get_p() == NULL && S.get_sibling()->is_leaf()) {
			//do not count
		} else {  //a non-trivial bipartition
			vector<int> cluster = S.get_cluster();
			Node* lca_mapping_in_T_prime;
			bool f = false;
			compute_lca_mapping_helper_2(cluster, &T_prime, lca_mapping_in_T_prime, f);

			if ( (cluster.size()) != (lca_mapping_in_T_prime->number_of_leaves()) ) { //f(u)=0 , here S is actually u --->NOTE you CAN'T use get_cluster_size() cuz it will return cluster size of "this" in T not T'
				//cout << "..........cluster corrsponding to the following does not exist in T': " << S.str_subtree() << endl;
				//cout << ".........................................cluster corrsponding to lca: " << lca_mapping_in_T_prime-> str_subtree() << endl;

				num_clusters_not_in_T_prime ++;
			}
		}
	}

}


//calculates the (weighted_RF_dist - num_bipartitions_in_T), that's because the "num_bipartitions_in_T" is constant when
//comparing RF_dist of two different T's and S
void find_weighted_rf_dist(Node& S, Node& T_prime, int& dist) {
	//if S curresponds to a trivial bipartition, i.e. bipartition with one leaf on one side, DO NOT count it. There are 2 cases:
	//1- it is a leaf
	//2- it is child of root, and it's sibling is a leaf
	//Also we don't count bipartition curresponding to roo ROOT
	if (S.is_leaf()) { //leaf
		return;
	}
	else {
		//preorderly
		list<Node *>::iterator c;
		list<Node *> children = S.get_children();
		for (c = children.begin(); c != children.end(); c++) {
			find_weighted_rf_dist(**c, T_prime, dist);
		}

		if (S.get_p() == NULL) { //root
			//do not count
		} else if (S.get_p()->get_p() == NULL && S.get_sibling()->is_leaf()) {
			//do not count
		} else {  //a non-trivial bipartition
			vector<int> cluster = S.get_cluster();
			Node* lca_mapping_in_T_prime;
			bool f = false;
			compute_lca_mapping_helper_2(cluster, &T_prime, lca_mapping_in_T_prime, f);

			//cout << ".....................Node u is: " << S.str_subtree() << endl;
			//cout << "a, i.e. lca_mapping_in_T_prime: " << lca_mapping_in_T_prime->str_subtree() << endl;

			//f(u)=0 , here S is actually u --->NOTE you CAN'T use get_cluster_size() cuz it will return cluster size of "this" in T not T'
			if ( (cluster.size()) != (lca_mapping_in_T_prime->number_of_leaves()) ) { //if u belongs to X (see your notes)
				//cout << "weight of this bipartition: " << S.get_edge_weight() << endl;
				dist += S.get_edge_weight();
			} else {  //if u belongs to Z (see your notes)
				dist += (S.get_edge_weight() - 2);
			}
		}
	}

}



//finds node for which "alpha - beta" is maximized
//ASSUMES T is the root of Q, i.e. only containing VALID regraft places for v under examination
void find_best_regraft_place(Node& T, Node*& best_regraft_place, int& max) {
	NUM_SPR_NGHBRS++;
	int best = T.get_alpha() - T.get_beta();
	//cout << "....alpha-beta....for: " << T.str_subtree() << ": " <<  T.get_alpha() << "-" << T.get_beta() << " = " << best << endl;
	if (best > max) {
		//cout << "is better: " << T.str_subtree() << endl;
		max = best;
		best_regraft_place = &T;
	}

	//after considering each node to be pruned, we update alpha beta values, BUT they should be reset to 0
	//I do it here cuz 1- I don't change T untill I find the best spr move (i.e. best node to be pruned and best place for it to be regrafted),
	//2-this is the last time I will need current alpha beta for the current v to be pruned
	T.set_alpha(0);
	T.set_beta(0);

	list<Node *>::iterator c;
	list<Node *> children = T.get_children();
	for (c = children.begin(); c != children.end(); c++) {
		find_best_regraft_place(**c, best_regraft_place, max);
	}
}


//given trees T, S and node v on T to be pruned, finds, in O(n), the best place to regraft v so that RF-dist(T,S) is minimized,
//i.e. it solves SPR(v) for T (and S as source tree)
//it changes the T to its best SPR neighbor, too.
//Returns best place to regraft
Node* apply_SPR_RS_algorithm_to_find_best_regraft_place(Node& T, Node& v, Node* source_trees_array[], set<string> non_shared_taxon_arr[], bool weighted) {

	//cout << "--------------v to be pruned: " << v.str_subtree() << " --------------\n";
	Node* best_regraft_place = NULL;
	Node* Q = NULL;
	Node* Q_in_restricted_st = NULL;
	Node* v_in_restricted_st = NULL;
	Node* old_sibling = NULL;
	Node* R = NULL;

	//1- preprocessing step:
	////////////////////////////////////////////////////////
	//we make tree R = SPR(v,rt(T)) from T with the following steps:
	//1- create a new node, R, and make rt(T) to be its child
	//2- Make the spr move as v->spr(rt(T),0), i.e. prune p(v) from tree and regraft it between rt(T) and R
	//depending on whether parent of v is T, preprocessin step (and hence step 2 of alg) should be implemented a little different

	//2- Calculating alpha and beta in Q:
	///////////////////////////////////////////////////////////
	//The algorithm traverses through S, and updates alpha and beta values in Q (which is the parameter T).
	//Note I pass T (rather than R) to the following function, since T repreents Q in the paper's notation, and we only care about alpha&beta vaues in Q.


	//BUT if v's parent is T (root), then handle it separately
	if (v.get_p() == &T) {
		//cout << "---v is child of T, and should be handled differently---" << endl;
		Q = v.get_sibling();

		for (int i = 0; i <::NUMBER_OF_SOURCE_TREES_ZZ; i++) {
			//cout << "current source tree:" << source_trees_array[i]->str_subtree() << endl;

			//There are two cases where SPR move won'y change RF distance between current ST and source_trees_array[i],
			//and thus we should ignore them:
			//////////////////////////case 1:
			//if the induced subtree at v, doesn't have any taxa in a source tree, S, then pruning and regrafting
			//v (at any place) will not change the RF distance, i.e.  RF_dist(T,S) = RF_dist(T',S)
			//////////////////////////case 2:
			//if the induced subtree at v contains all the taxa in source_trees_array[i].
			//since the clusters are sorted we can use std:include()

			vector<int> v_cluster = v.find_cluster_int_labels();
			vector<int> current_source_trre_root_cluster = source_trees_array[i]->get_cluster();
			vector<int> intersection;
			set_intersection(v_cluster.begin(), v_cluster.end(), current_source_trre_root_cluster.begin(), current_source_trre_root_cluster.end(), back_inserter(intersection));

			//case 1
			if (intersection.size() == 0 ) {
				//cout << "the induced subtree at v does not include any taxa in current source tree, thus we ignore it.\n";
				continue;
			}

			//case 2
			if (intersection.size() == current_source_trre_root_cluster.size()) {
				//cout << "the induced subtree at v includes all taxa in current source tree, thus we ignore it.\n";
				continue;
			}
			// In other words, if root of current source tree's LCA is v or a descendant of it.
			//Overlooking above case will cause 'double free' error when you try delete S_prime since S_prime is going to be a singletone node
			//which has already been free'ed once in suppress_nodes_with_mapping_in_Rv()


			//find restricted supertrees to current source tree
			Node* restricted_suptree = build_tree(T.str_subtree());
			restricted_suptree->copy_fields_for_supertree(T);	//prenum is also copied so no adjustTree() needed
			restrict_supertree(*restricted_suptree, non_shared_taxon_arr[i]);

			v_in_restricted_st = restricted_suptree->find_by_prenum(v.get_preorder_number());//corresponding node in restricted_st

			//if sibling is a leaf not present in current source tree,  then there is no valid place to regraft it, thus ignore
			if (! (Q_in_restricted_st =  v_in_restricted_st->get_sibling()) ) {
				continue;
			}

			compute_lca_mapping_S_R(*source_trees_array[i], *restricted_suptree);  //No need for R here, T will do the work
			set_cluster_size_in_supertree(restricted_suptree);

			//Find S' in lemma 12:
			Node* S_prime = build_tree(source_trees_array[i]->str_subtree() + ";"); //make a copy of S (not the best way to mak ea copy, I know, should be fixed)
			S_prime->copy_fields_for_source_tree(*source_trees_array[i]);
			suppress_nodes_with_mapping_in_Rv(*S_prime, *v_in_restricted_st);  //S' in lemma 12 is now  constructed from S

			traverse_S_and_update_alpha_beta_in_Q(*Q, *Q_in_restricted_st, *v_in_restricted_st, *source_trees_array[i], *S_prime, weighted);

			S_prime->delete_tree();  //prevent memory leak
			restricted_suptree->delete_tree();

		}//end of for-loop

	} else { //normal case
		//cout << "---normal case: i.e. v is not the child of root---" << endl;

		Q = &T; //here Q is going to be T
		R = new Node();
		R->add_child(&T);
		//cunstructing R which is SPR(v,rt(T))
		int which_sibling = 0;
		old_sibling = v.spr(&T, which_sibling);
		adjustTree(R);//DON'T FORGET THIS!!!

		for (int i = 0; i <::NUMBER_OF_SOURCE_TREES_ZZ; i++) {

			//There are two cases where SPR move won'y change RF distance between current ST and source_trees_array[i],
			//and thus we should ignore them:
			//////////////////////////case 1:
			//if the induced subtree at v, doesn't have any taxa in a source tree, S, then pruning and regrafting
			//v (at any place) will not change the RF distance, i.e.  RF_dist(T,S) = RF_dist(T',S)
			//////////////////////////case 2:
			//if the induced subtree at v contains all the taxa in source_trees_array[i].
			//since the clusters are sorted we can use std:include()

			vector<int> v_cluster = v.find_cluster_int_labels();
			vector<int> current_source_trre_root_cluster = source_trees_array[i]->get_cluster();
			vector<int> intersection;
			set_intersection(v_cluster.begin(), v_cluster.end(), current_source_trre_root_cluster.begin(), current_source_trre_root_cluster.end(), back_inserter(intersection));

			//case 1
			if (intersection.size() == 0 ) {
				//cout << "the induced subtree at v does not include any taxa in current source tree, thus we ignore it.\n";
				continue;
			}

			//case 2
			if (intersection.size() == current_source_trre_root_cluster.size()) {
				//cout << "the induced subtree at v includes all taxa in current source tree, thus we ignore it.\n";
				continue;
			}
			// In other words, if root of current source tree's LCA is v or a descendant of it.
			//Overlooking above case will cause 'double free' error when you try delete S_prime since S_prime is going to be a singletone node
			//which has already been free'ed once in suppress_nodes_with_mapping_in_Rv()


			//find restricted supertrees to current source tree
			Node * temp = build_tree(R->str_subtree());
			///I need the following line cuz:
			//build_tree() ignores extra parenthesis which here represent Node R with one child
			Node * R_for_restricted_supertree = new Node();
			R_for_restricted_supertree->add_child(temp);
			R_for_restricted_supertree->copy_fields_for_supertree(*R);	//prenum is also copied so no adjustTree() needed
			restrict_supertree(*R_for_restricted_supertree, non_shared_taxon_arr[i]);

			//corresponding nodes in restricted R
			v_in_restricted_st = R_for_restricted_supertree->find_by_prenum(v.get_preorder_number());
			Q_in_restricted_st = R_for_restricted_supertree->find_by_prenum(T.get_preorder_number());


			compute_lca_mapping_S_R(*source_trees_array[i], *R_for_restricted_supertree);
			set_cluster_size_in_supertree(R_for_restricted_supertree);

			//Find S' in lemma 12:
			Node* S_prime = build_tree(source_trees_array[i]->str_subtree() + ";"); //make a copy of S (not the best way to mak ea copy, I know, should be fixed)
			S_prime->copy_fields_for_source_tree(*source_trees_array[i]);
			suppress_nodes_with_mapping_in_Rv(*S_prime, *v_in_restricted_st);  //S' in lemma 12 is now  constructed from S

			traverse_S_and_update_alpha_beta_in_Q(*Q, *Q_in_restricted_st, *v_in_restricted_st, *source_trees_array[i], *S_prime, weighted);

			S_prime->delete_tree();
			R_for_restricted_supertree->delete_tree();	//this will delete "restricted_st" as well

			//R clean up at the end since I will need to call find_best_regraft_place() on Q before restoring T
		}//end of for-loop

	}

	//cout << "\nQ is: " << Q->str_subtree()  << endl;
	//cout << "The final alpha beta values in nodes in Q after the algorithm finishes are as follow: " << endl;
	//preorder_traversal(*Q);
	int max_alpha_minus_beta = INT_MIN;
	find_best_regraft_place(*Q, best_regraft_place, max_alpha_minus_beta);  //NOTE Q should be passed

	//restore T and prevent mem-leack
	if (v.get_p() != &T) {
		int which_sibling = 0;
		v.spr(old_sibling, which_sibling);  // Note T has been changed here, and it should be retrieved to its original form
		T.cut_parent(); // after finding best place to regraft and make the corresponding SPR move, R is a node with only one child T, we don't need R anymore.
		//(v.get_p()) -> cut_parent(); // NOTE v is not a descendant of T anymore, they are SEPARATE trees
		R->delete_tree();
	}


	return best_regraft_place;
}


//Given the way R is constructed, lca_mapping of any node in S, is either in T_v OR in T_T. (NOTE this is a recursive function and S plays the role of u in the alg in paper)
//The way I implemented here, avoids some duplicate work (in compare to when I pass R as parameter instead of T and v).
void traverse_S_and_update_alpha_beta_in_Q(Node & Q, Node & Q_in_restricted_st, Node & v_in_restricted_st,
        Node & S, Node & S_prime, bool weighted) { //note S here is u in paper's notation

	//only for internal nodes, i.e. non-root and non-leaf
	if (S.is_leaf()) {  //if leaf
		//do nothing

	} else if (S.get_p() == NULL) { //if root
		list<Node *>::iterator c;
		list<Node *> children = S.get_children();
		for (c = children.begin(); c != children.end(); c++) {
			traverse_S_and_update_alpha_beta_in_Q(Q, Q_in_restricted_st, v_in_restricted_st, **c, S_prime, weighted);
		}
	} else {  //if internal node

		//cout << "++++++++++++++++++++node u in S being considered: " << S.str_subtree() << endl;
		//cout << "its lca_mapping is: " << S.get_lca_mapping()->str_subtree() << endl;
		Node* a;

		//Suppose for each u in I(S), lca_mpping(u)=a
		//There are 4 cases.
		//1- if u satisfies precondition of Lemma 9, i.e. a belongs to V(R_v)
		if (a = v_in_restricted_st.find_by_prenum(S.get_lca_mapping() -> get_preorder_number())) {
			//cout << "lemma 9" <<  endl;
			//do nothing
		}

		//2- if u satisfies precondition of Lemma 11, i.e. a belongs to V(Q_v) AND f_R(S)==0
		else if (a = Q_in_restricted_st.find_by_prenum(S.get_lca_mapping()  -> get_preorder_number() )) {
			if (S.get_cluster_size() != a->get_cluster_size()) { //is f_R(S)==0
				//cout << "lemma 11" <<  endl;
				//do nothing
			} else { //precondition of lemma 10 is true
				//3- if u satisfies precondition of Lemma 10, i.e. a belongs to V(T_v) AND f_R(S)==1
				//cout << "lemma 10" <<  endl;
				if (weighted) {
					int w = S.get_edge_weight();
					//a->increment_beta_in_all_descendants(w);
					(Q.find_by_prenum(a->get_preorder_number()))->increment_beta_in_all_descendants(w);
					//cout << "________w____increment beta on : " << a->str_subtree() << endl;
				} else {
					//a->increment_beta_in_all_descendants(1);
					(Q.find_by_prenum(a->get_preorder_number()))->increment_beta_in_all_descendants(1);
					//cout << "____________increment beta on : " << a->str_subtree() << endl;
					//cout << "node on which increment is called: " << (Q.find_by_prenum(a->get_preorder_number()))->str_subtree() << endl;
				}
			}
		}

		//Lemma 12's precondition: a=T->p() AND |L(R_b)|+|L(R_v)|=|L(S_u)|
		//Note: if S's lca_mapping is not found in T_v or T_T, then it's lca_mapping is definitely parent of T (or v) which is the only child of R
		else {
			//cout << "maybe lemma 12" <<  endl;

			Node* u_in_S_prime = S_prime.find_by_prenum(S.get_preorder_number()); //find u (S in parameters) in S'
			Node* lca_mapping_lemma12;
			vector<int> cluster = u_in_S_prime->find_cluster_int_labels();  //note u's cluster field is not valid in S', and I should find it again
			bool f = false;
			compute_lca_mapping_helper_2(cluster, Q_in_restricted_st.get_p(), lca_mapping_lemma12, f);
			Node* b_in_lemma12 = (Q_in_restricted_st.get_p())->find_by_prenum(lca_mapping_lemma12  -> get_preorder_number() );



			if ( ((b_in_lemma12->find_leaves()).size()) + ((v_in_restricted_st.find_leaves()).size()) == (S.get_cluster_size()) ) { //|L(R_b)|+|L(R_v)|=?|L(S_u)|
				//cout << (b_in_lemma12->find_leaves()).size() << "+" << ((v_in_restricted_st.find_leaves()).size()) << " == " << (S.get_cluster_size()) << endl;
				if (weighted) {
					int w = S.get_edge_weight();
					//b_in_lemma12->increment_alpha_in_all_descendants(w);
					(Q.find_by_prenum(b_in_lemma12->get_preorder_number()))->increment_alpha_in_all_descendants(w);
					//cout << "_____w_______increment alpha on : " << b_in_lemma12->str_subtree() << endl;
				} else {
					//b_in_lemma12->increment_alpha_in_all_descendants(1);
					//cout << "node on which increment is called: " << (Q.find_by_prenum(b_in_lemma12->get_preorder_number()))->str_subtree() << endl;
					(Q.find_by_prenum(b_in_lemma12->get_preorder_number()))->increment_alpha_in_all_descendants(1);
					//cout << "____________increment alpha on b lemma 12: " << b_in_lemma12->str_subtree() << endl;
				}

			}

		}

		list<Node *>::iterator c;
		list<Node *> children = S.get_children();
		for (c = children.begin(); c != children.end(); c++) {
			traverse_S_and_update_alpha_beta_in_Q(Q, Q_in_restricted_st, v_in_restricted_st, **c, S_prime, weighted);
		}


	}
}

//pre-orderly traverses tree S, and suppresses all nodes in it for which lca_mapping is in T_v
void suppress_nodes_with_mapping_in_Rv(Node & S, Node & v) {
	if (v.find_by_prenum(S.get_lca_mapping()  -> get_preorder_number())) {
		//cout << "removed node : " << S.str_subtree() << endl;
		Node* p = S.get_p();
		//cout << "parent before deletion: " << p->str_subtree() << endl;
		p -> delete_child(&S);
		S.delete_tree();  //prevent memory leak!
		//cout << "parent after deletion: " << p->str_subtree() << endl;

		//NOTE the following if statement SHOULD NOT be done, cuz if you contract nodes in S', how can you find corresponding nodes of S in S'??
		/*
		if(p -> get_children().size() == 1) {
		  Node* pp = p->get_p();
		  p -> contract_node();
		}
		*/
	} else {
		list<Node *>::iterator c;
		list<Node *> children = S.get_children();
		for (c = children.begin(); c != children.end(); c++) {
			suppress_nodes_with_mapping_in_Rv(**c, v);
		}
	}


}

//Why do we need all the following 3 parameters:
//S: we need S to make a copy of it and make S' from it (now all prenum, cluster, ... are invalid)
//u: u is the node for which I know I will need to know b_in_lemma12 field (u exists in both S and S')
//R: b_in_lemma12 is the mapping of u in R
//v: cuz S' is constructed from S based on v
void find_b_in_lemma12(Node & u, Node & S, Node & R, Node & v) {
	int prenum_of_u = u.get_preorder_number();
	Node S_prime = Node(S); //make a copy of S

	//TO BE COMPLETED ..

}


//for all u in I(S), find the LCA mapping of u:
//1-Let's cluster(u) be set of leaves in subtree rooted u in tree S,
//2-Find LCA of cluster(u) in tree R, and set its lca_mapping
//NOTE: lca_mapping is the prenum f the LCA node in R
void compute_lca_mapping_S_R(Node & S, Node & R) {

	if (S.is_leaf()) {
		Node* lca_mapping;
		vector<int> cluster = S.get_cluster();
		bool f = false;
		compute_lca_mapping_helper_2(cluster, &R, lca_mapping, f);
		S.set_lca_mapping(lca_mapping);
	} else {
		Node* lca_mapping;
		vector<int> cluster = S.get_cluster();
		bool f = false;
		compute_lca_mapping_helper_2(cluster, &R, lca_mapping, f);
		S.set_lca_mapping(lca_mapping);
		//recursively set children's lca_mapping
		list<Node *>::iterator c;
		list<Node *> children = S.get_children();
		for (c = children.begin(); c != children.end(); c++) {
			compute_lca_mapping_S_R(**c, R);
		}
	}
}

void compute_lca_mapping_helper_1(vector<int>& cluster, Node * R, int& lca) {
	bool lca_must_be_in_lower_levels = false;
	Node* which_child;

	list<Node *>::iterator c;
	list<Node *> children = R->get_children();
	for (c = children.begin(); c != children.end(); c++) {
		vector<int> c_cluster = (*c)-> find_cluster_int_labels();
		if (includes(c_cluster.begin(), c_cluster.end(), cluster.begin(), cluster.end())) {  //if cluster is a subset of c_cluster
			lca_must_be_in_lower_levels = true;
			which_child = *c;
			break;
		}
	}

	if (lca_must_be_in_lower_levels) {
		compute_lca_mapping_helper_1(cluster, which_child, lca);
	} else {
		lca = R -> get_preorder_number();
	}
}

//assuming "cluster" is the labels for which we are lookig for the LCA, the function finds the LCA,
//and assigns its prenum to parameter "lca"
void compute_lca_mapping_helper_2(vector<int>& cluster, Node * R, Node*& lca, bool & found) {
	if (! R->is_leaf()) {
		list<Node *>::iterator c;
		list<Node *> children = R->get_children();
		for (c = children.begin(); c != children.end(); c++) {
			compute_lca_mapping_helper_2(cluster, *c, lca, found);
		}
	}

	if (found) { // to return from nested recursive calls after
		return;
	}

	if (R-> get_lca_hlpr() == cluster.size()) {
		lca = R;
		found = true;
		R->set_lca_hlpr(0); //to avoid for reseting lca_hlpr values to 0 for the next time we compute lca_mapping
		return;
	}

	if (R->is_leaf()) {
		if (cluster.size() == 1) { //leaves can have lca as well
			if (R->get_int_label() == cluster[0]) {
				lca = R;
				found = true;
				return;
			}
		} else {
			//assuming "cluster" is SORTED!!
			if (binary_search(cluster.begin(), cluster.end(), R->get_int_label())) {
				( R->get_p() ) -> increase_lca_hlpr_by(1);
			}
		}
	} else {
		( R->get_p() ) -> increase_lca_hlpr_by(R->get_lca_hlpr());
		R->set_lca_hlpr(0); //to avoid for reseting lca_hlpr values to 0 for the next time we compute lca_mapping
	}
}


void find_f_u_with_regard_to_R(vector<int>& cluster, vector<int>& source_tree_leaf_set, Node * R, Node*& lca, bool & found, int& f_u) {
	if (! R->is_leaf()) {
		list<Node *>::iterator c;
		list<Node *> children = R->get_children();
		for (c = children.begin(); c != children.end(); c++) {
			find_f_u_with_regard_to_R(cluster, source_tree_leaf_set, *c, lca, found, f_u);
		}
	}

	if (found) { // to return from nested recursive calls afterward
		return;
	}

	if (R-> get_lca_hlpr() == cluster.size()) {
		lca = R;
		found = true;
		if (R->get_XXX() == true) {
			f_u = 0;
			R->set_XXX(false);  //reste fo next calls
		} else {
			f_u = 1;
		}
		R->set_lca_hlpr(0); //to avoid for reseting lca_hlpr values to 0 for the next time we compute lca_mapping
		return;
	}

	if (R->is_leaf()) {
		if (cluster.size() == 1) { //leaves can have lca as well
			if (R->get_int_label() == cluster[0]) {
				lca = R;
				found = true;
				return;
			}
		} else {
			//assuming "cluster" is SORTED!!
			if (binary_search(cluster.begin(), cluster.end(), R->get_int_label())) {
				( R->get_p() ) -> increase_lca_hlpr_by(1);
			} else {
				if (binary_search(source_tree_leaf_set.begin(), source_tree_leaf_set.end(), R->get_int_label())) {
					( R->get_p() ) -> set_XXX(true);
				}
			}
		}
	} else {
		( R->get_p() ) -> increase_lca_hlpr_by(R->get_lca_hlpr());
		R->set_lca_hlpr(0); //to avoid for reseting lca_hlpr values to 0 for the next time we compute lca_mapping
		if (R->get_XXX() == true) {
			( R->get_p() ) -> set_XXX(true);
		}
		R->set_XXX(false);  //reste fo next calls
	}
}


//it does two things:
//1- sets int labels for taxa on supertree
//2- saves the mapping for future refrence (labeling source trees, ...)
void find_int_labels_for_leaves_in_supertree(Node * root, int& current_label, unordered_map<string, int>& map) {
	if (root->is_leaf()) {
		map.insert(pair<string, int>(root->get_name(), current_label));
		root->set_int_label(current_label);
		current_label++;
	} else {
		for (auto child : root->get_children()) {
			find_int_labels_for_leaves_in_supertree(child, current_label, map);
		}
	}
}

//set int labels for leaves in one source tree
void set_int_labels_for_leaves_in_source_tree(Node * root, unordered_map<string, int>& map) {
	if (root->is_leaf()) {
		root->set_int_label(map.find(root->get_name())->second);
	} else {
		for (auto child : root->get_children()) {
			set_int_labels_for_leaves_in_source_tree(child, map);
		}
	}
}

//finds clusters associated with each node in each SOURCE tree, and stores it in DS to avoid repeated work
//this method should be called before calling find_lca()
//note a node's cluster is a SORTED vector contining int_labels of the leaves in subtree induced at that node
void set_cluster_and_cluster_size(Node * T) {
	if (T-> is_leaf()) {
		std::vector<int> v {T->get_int_label()};
		T-> set_cluster(v);
		T-> set_cluster_size(1);
	} else {
		T-> set_cluster(T->find_cluster_int_labels());
		T-> set_cluster_size((T->get_cluster()).size());
	}

	list<Node *>::iterator c;
	list<Node *> children = T->get_children();
	for (c = children.begin(); c != children.end(); c++) {
		set_cluster_and_cluster_size(*c);
	}
}

//find and set cluster_size for each node in the supertree
void set_cluster_size_in_supertree(Node * T) {
	if (T-> is_leaf()) {
		T-> set_cluster_size(1);
	} else {
		T-> set_cluster_size((T->find_leaves()).size());
	}

	list<Node *>::iterator c;
	list<Node *> children = T->get_children();
	for (c = children.begin(); c != children.end(); c++) {
		set_cluster_size_in_supertree(*c);
	}
}

//reset invalid fields in ST after each iteration to initial values
void reset_fields_to_initial_values(Node * T) {
	T->set_cluster_size(0);
	//T->set_lca_hlpr(0);
	T->set_alpha(0);
	T->set_beta(0);

	list<Node *>::iterator c;
	list<Node *> children = T->get_children();
	for (c = children.begin(); c != children.end(); c++) {
		reset_fields_to_initial_values(*c);
	}
}




//*******************added for SPR neighborhood*********************
//******************************************************************
//applies all the possible spr's on "tree", and writes them into the file "z_spr_neighbours"
int produce_all_spr_neighbors(Node * myTree, int number_of_taxa) {

	int number_of_neighbors = 0;


	//The following line deletes the contents of z_spr_neighbours
	const char* neighbors_file = "z_spr_neighbours";
	std::ofstream ofile(neighbors_file, ios_base::trunc);

	for (int i = 1; i < number_of_taxa; i++) {

		Node* spr_on = myTree->find_by_prenum(i);
		int which_sibling = 0;
		/*
		cout << "************************************************************" << endl;
		cout << "spr_on's pre-order number = i = " << spr_on->get_preorder_number() << endl;
		cout << "new_sibling's pre-order number = j = " << spr_on->get_preorder_number() << endl;
		cout << "original  tree: " << myTree->str_subtree() << endl;
		cout << "************************************************************" << endl;
		*/

		for (int j = 1; j < number_of_taxa; j++) {

			if (j != i) {

				Node* new_sibling = myTree->find_by_prenum(j);
				bool good_spr = true;
				//cout<< "i = " <<  i << ", j = " << j << endl;



				//bad spr: check whether new_sibling is parent
				Node* parent = spr_on->get_p();
				if (parent->get_preorder_number() == j) {
					//cout << "-----BAD SPR WAS IGNORED----- \n" << endl;
					continue;
				}


				//bad spr: check whether new_sibling is a descendant of spr_on
				if (! spr_on->is_leaf()) {
					vector<Node *> descendants = spr_on->find_descendants();
					vector<Node *>::iterator it;
					for (it = descendants.begin(); it != descendants.end(); ++it) {
						if (new_sibling->get_preorder_number() == (*it)->get_preorder_number()) {
							//cout << "-----BAD SPR WAS IGNORED----- \n" << endl;
							good_spr = false;
							continue;       //goes out of the inner loop which is 4 lines above
						}
					}
				}




				//bad spr: if new sibling is old sibling ignore it cuz this spr-move results to the original tree
				list<Node *>::iterator itr;
				for (itr = (spr_on->parent())->get_children().begin(); itr != (spr_on->parent())->get_children().end(); ++itr) {
					if ((*itr)->get_preorder_number() == j) {
						good_spr = false;
						continue;
					}
				}


				if (good_spr) {
					Node* undo = spr_on->spr(new_sibling, which_sibling);
					adjustTree(myTree);
					number_of_neighbors ++;
					string new_tree = myTree ->str_subtree();
					write_line_to_File(new_tree, neighbors_file);
					//cout<<  "tree after spr: " << new_tree << "\n" <<endl;

					//restore the original tree
					spr_on->spr(undo, which_sibling);
					adjustTree(myTree);
				}

			}

		}

	}



//    cout << "**********************************************************" << endl;
//    cout << "Total Number Of Taxa In Trees: " << number_of_taxa << endl;
//    cout << "Current Node (Supertree) : " << myTree->str_subtree() << endl;
//    cout << "#of spr_neighbors of initial supertree = " << number_of_neighbors << endl;
//    cout << "Trees were written in z_spr_neighbours in newick format." << endl;
//    cout << "**********************************************************" << endl;
	return number_of_neighbors;

}

//writes "s\n" into file "file_name"
void write_line_to_File(string s, const char* file_name) {
	ofstream myFile;
	myFile.open(file_name, std::ios_base::app);
	myFile << s + ";\n";
	myFile.close();
}


void adjustTree(Node * myTree) {
	//myTree->set_depth(0);
	//myTree->fix_depths();
	myTree->preorder_number();
	//myTree->edge_preorder_interval();
}




//returns the total number of nodes in the clade of "node", i.e. (#of taxa)+(#of internal nodes)
//This method is a modification of the set_preorder_number() method in node.h
void total_number_of_nodes(Node * node, int& total_nodes) {
	total_nodes++;
	list<Node *>::iterator c;
	list<Node *> children = node->get_children();
	for (c = children.begin(); c != children.end(); c++) {
		total_number_of_nodes(*c, total_nodes);
	}
}



/////////////////////////////Added for restrict_st()/////////////////////////
//for each "taxon" in "non_shared_taxon_set", traverse the tree and find leaf whose name is "taxon", remove it.
//if parent is left with one child, then contract parent.
void restrict_supertree(Node & supertree, set<string>& non_shared_taxon_set) {
	Node* n;  //this willl be the node whose "name" is "taxon"
	Node* p;  //this will be the parent of "n"

	for (auto taxon : non_shared_taxon_set) {

		n = 0;
		find_node_to_be_removed(supertree, taxon, n);
		//cout << n->get_name() << endl;
		p = n->get_p();
		//cout << "parent before deletion---------" << p->str_subtree() << endl;
		p -> delete_child(n);
		delete n;
		//cout << "parent after delete_child()----" << p->str_subtree() << endl;
		/*
		if (p -> get_children().size() == 1) {
			Node* pp = p->get_p();
			cout << "grand-pa before contrac--------" << pp -> str_subtree() << endl;
			p -> contract_node();
			cout << "grand-pa after contract_node()-" << pp -> str_subtree() << endl;
		}
		*/
	}
	//cout << "final restricted ST: " << supertree.str_subtree() << endl;
}

void find_node_to_be_removed(Node & root, string taxon, Node* & taxan_to_be_removed) {

	if (taxan_to_be_removed) {	//to save some time, BUT NOTE THAT "taxan_to_be_removed" should be initialized to 0 before it's passed.
		return;
	}

	if (root.is_leaf()) {
		if (root.get_name() == taxon) {
			taxan_to_be_removed = &root;
		}
	} else {
		for (auto child : root.get_children()) {
			find_node_to_be_removed(*child, taxon, taxan_to_be_removed);
		}
	}
}




set<string> find_non_common_taxa_set(const string & supertree, const string & source_tree) {
	set<string> supertree_taxon_set;
	set<string> source_tree_taxon_set;
	set<string> non_shared_taxon_set;

	//get taxon set of source tree
	regex reg("[^\\w]+");
	string temp1 = regex_replace(source_tree, reg, " ");  //replacing all parenthesis, comma, and semi-colon with white space
	temp1 = regex_replace(temp1, regex("^ +| +$|( ) +"), "$1"); //removing leading, trailing and extra spaces
	vector<string> source_tree_taxon_vector = split(temp1, ' ');
	for (auto i : source_tree_taxon_vector) {
		source_tree_taxon_set.insert(i);
	}


	//get taxon set of supertree
	string temp2 = regex_replace(supertree, reg, " ");  //replacing all parenthesis, comma, and semi-colon with white space
	temp2 = regex_replace(temp2, regex("^ +| +$|( ) +"), "$1"); //removing leading, trailing and extra spaces
	vector<string> supertree_taxon_vector = split(temp2, ' ');
	for (auto i : supertree_taxon_vector) {
		supertree_taxon_set.insert(i);
	}

	//find set difference
	set_difference( supertree_taxon_set.begin(), supertree_taxon_set.end(), source_tree_taxon_set.begin(), source_tree_taxon_set.end(), inserter(non_shared_taxon_set, non_shared_taxon_set.begin()));

	//for symmetric difference you can use the following
	//set_symmetric_difference( supertree_taxon_set.begin(), supertree_taxon_set.end(), source_tree_taxon_set.begin(), source_tree_taxon_set.end(), inserter(non_shared_taxon_set, non_shared_taxon_set.begin()));

	return non_shared_taxon_set;
}



void split(const string & s, char delim, vector<string> &elems) {
	stringstream ss;
	ss.str(s);
	string item;
	while (getline(ss, item, delim)) {
		elems.push_back(item);
	}
}


vector<string> split(const string & s, char delim) {
	vector<string> elems;
	split(s, delim, elems);
	return elems;
}
