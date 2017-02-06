
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
void find_prenum_of_node_to_be_removed(Node& tree, string taxon, int& prenum_of_node_to_be_removed);
std::string ReplaceAll(std::string str, const std::string& from, const std::string& to);
static unsigned int NUMBER_OF_SOURCE_TREES_ZZ;


////////////////RFS implementation//////////////////////
#include <unordered_map>
void initialize_alpha_beta(Node* Q);
void set_cluster_size(Node* Q);
void find_int_labels_for_leaves_in_supertree(Node* root, int& current_label, unordered_map<string, int>& map);
void set_int_labels_for_leaves_in_source_tree(Node* root, unordered_map<string, int>& map);
void compute_lca_mapping_S_R(Node& S, Node& R);
void compute_lca_mapping_helper_1(vector<int>& cluster, Node* R, int& lca);
void compute_lca_mapping_helper_2(vector<int>& cluster, Node* R, Node*& lca, bool& found);
void set_cluster_in_source_tree(Node* T);
void set_cluster_size_in_supertree(Node* T);
void reset_fields_to_initial_values(Node* T);
Node* apply_SPR_RS_algorithm_to_find_best_regraft_place(Node& T, Node& spr_on, Node& S, bool weighted);
void traverse_S_and_update_alpha_beta_in_Q(Node& T, Node& v, Node& S, Node& S_prime, bool weighted);
void find_b_in_lemma12(Node& S, Node& R); // i'm not using this now
void suppress_nodes_with_mapping_in_Rv(Node& S_prime, Node& v);
void find_best_regraft_place(Node& n, Node*& best_regraft_place, int& max);
void preorder_trversal(Node& n);
int find_best_node_to_prune_and_its_best_regraft_place(Node& T, Node& S, Node* & best_node_to_prune, Node* & best_node_to_regraft, bool weighted);
void find_F_T(Node& S, Node& T, int& best, bool weighted);
void put_internal_nodes_in_vector(Node& n, vector<Node*>& nodes);
void put_all_nodes_in_vector(Node& n, vector<Node*>& nodes);

void print_weighted_tree(Node& n);

template <typename T>
bool IsSubset(std::vector<T> A, std::vector<T> B)
{
  //std::sort(A.begin(), A.end());
  //std::sort(B.begin(), B.end());
  return std::includes(A.begin(), A.end(), B.begin(), B.end());
}


/////////////////////////////     main()     /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
int NUM_SPR_NGHBRS = 0;

int main(int argc, char** argv) {
  srand(time(NULL));

  //running time
  clock_t start_time, finish_time;
  start_time = clock();

  string suptree = "(a,(c,(d,(e,(f,(b,g))))));";
  //string suptree = "(((c,f),(d,e)),(a,(b,g)));";
  //string input_tree = "(c,(f,(e,(d,(g,(b,a))))));";
  string input_tree = "(a,((c,f),(d,(e,(b,g)))));";


  //string suptree = "((Pygoscelis_adeliae,(((Eudyptes_chrysolophus,Eudyptes_chrysocome),Aptenodytes_patagonicus),((Pygoscelis_antarctica,Pygoscelis_papua),((Eudyptes_pachyrhynchus,Megadyptes_antipodes),(Spheniscus_demersus,Eudyptula_minor))))),((Gavia_stellata,Gavia_immer),(((Oceanodroma_melania,(Oceanodroma_tethys,Halocyptena_microsoma)),(((Oceanodroma_castro,Hydrobates_pelagicus),(Oceanodroma_furcata,Oceanodroma_hornbyi)),(Oceanodroma_leucorhoa,Oceanodroma_tristrami))),((((Garrodia_nereis,Pelagodroma_marina),(Fregetta_grallaria,Fregetta_tropica)),Oceanites_oceanicus),(((Pelecanoides_garnotii,((Pelecanoides_georgicus,Pelecanoides_magellani),Pelecanoides_urinatrix)),(((Pterodroma_axillaris,(Pterodroma_nigripennis,Pterodroma_cervicalis)),((((Pterodroma_alba,(((Pterodroma_baraui,Pterodroma_arminjoniana),(Pterodroma_neglecta,Pterodroma_externa)),(Pterodroma_heraldica,(Pterodroma_sandwichensis,Pterodroma_phaeopygia)))),Pterodroma_inexpectata),(Pterodroma_ultima,(Pterodroma_solandri,(Pterodroma_mollis,((Pterodroma_hasitata,((Pterodroma_madeira,Pterodroma_feae),Pterodroma_cahow)),(Pterodroma_magentae,(Pterodroma_incerta,(Pterodroma_lessonii,Pterodroma_macroptera)))))))),(Pterodroma_hypoleuca,((Pterodroma_pycrofti,Pterodroma_longirostris),(Pterodroma_brevipes,(Pterodroma_leucoptera,(Pterodroma_cookii,Pterodroma_defilippiana))))))),((Pagodroma_nivea,((Daption_capense,((Macronectes_halli,Macronectes_giganteus),(Fulmarus_glacialis,Fulmarus_glacialoides))),Thalassoica_antarctica)),((Procellaria_cinerea,((Lugensa_brevirostris,((((Puffinus_gravis,(Puffinus_griseus,(Puffinus_creatopus,Puffinus_carneipes))),(Puffinus_bulleri,Puffinus_pacificus)),Puffinus_tenuirostris),((Calonectris_diomedea,Calonectris_leucomelas),(Puffinus_nativitatis,(((Puffinus_huttoni,Puffinus_gavia),(((Puffinus_opisthomelas,Puffinus_puffinus),Puffinus_auricularis),(Puffinus_assimilis,Puffinus_lherminieri))),(Puffinus_mauretanicus,Puffinus_yelkouan)))))),(Pseudobulweria_aterrima,Pseudobulweria_rostrata))),(((Procellaria_westlandica,(Procellaria_aequinoctialis,Procellaria_parkinsoni)),Bulweria_bulwerii),((Pachyptila_turtur,(Pachyptila_vittata,(Pachyptila_desolata,Pachyptila_salvini))),Halobaena_caerulea)))))),(((Phoebastria_irrorata,((Phoebastria_immutabilis,Phoebastria_nigripes),Phoebastria_albatrus)),((Diomedea_sanfordi,Diomedea_epomophora),(Diomedea_dabbenena,((Diomedea_gibsoni,Diomedea_antipodensis),(Diomedea_amsterdamensis,Diomedea_exulans))))),((Phoebetria_fusca,Phoebetria_palpebrata),((Thalassarche_bassi,Thalassarche_chlororhynchus),((Thalassarche_chrysostoma,(Thalassarche_melanophris,Thalassarche_impavida)),(Thalassarche_bulleri,(Thalassarche_cauta,(Thalassarche_salvini,Thalassarche_eremita))))))))))));";
  //string input_tree = "((Puffinus_lherminieri,Pterodroma_cervicalis),((((Puffinus_nativitatis,((((Procellaria_cinerea,(Lugensa_brevirostris,(((((((((Pygoscelis_adeliae,(Puffinus_auricularis,((Pygoscelis_antarctica,(Pygoscelis_papua,((Eudyptula_minor,Spheniscus_demersus),(Megadyptes_antipodes,Eudyptes_pachyrhynchus)))),(Aptenodytes_patagonicus,(Eudyptes_chrysocome,Eudyptes_chrysolophus))))),(Gavia_immer,Gavia_stellata)),(((Oceanodroma_hornbyi,(Oceanodroma_furcata,(Oceanodroma_castro,Hydrobates_pelagicus))),(Oceanodroma_leucorhoa,Oceanodroma_tristrami)),(Oceanodroma_melania,(Halocyptena_microsoma,Oceanodroma_tethys)))),(Oceanites_oceanicus,((Fregetta_tropica,Fregetta_grallaria),(Pelagodroma_marina,Garrodia_nereis)))),(((((Thalassarche_chrysostoma,(Thalassarche_melanophris,Thalassarche_impavida)),(Thalassarche_bulleri,(Thalassarche_cauta,(Thalassarche_eremita,Thalassarche_salvini)))),(Thalassarche_chlororhynchus,Thalassarche_bassi)),(Phoebetria_fusca,Phoebetria_palpebrata)),(((Diomedea_dabbenena,((Diomedea_exulans,Diomedea_amsterdamensis),(Diomedea_gibsoni,Diomedea_antipodensis))),(Diomedea_epomophora,Diomedea_sanfordi)),(Phoebastria_irrorata,(Phoebastria_albatrus,(Phoebastria_immutabilis,Phoebastria_nigripes)))))),(Pelecanoides_garnotii,(Pelecanoides_urinatrix,(Pelecanoides_georgicus,Pelecanoides_magellani)))),(((((Pterodroma_ultima,(Pterodroma_solandri,(((Pterodroma_magentae,(Pterodroma_incerta,(Pterodroma_lessonii,Pterodroma_macroptera))),(Pterodroma_hasitata,((Pterodroma_feae,Pterodroma_madeira),Pterodroma_cahow))),Pterodroma_mollis))),(Pterodroma_inexpectata,((((Pterodroma_externa,Pterodroma_neglecta),(Pterodroma_baraui,Pterodroma_arminjoniana)),(Pterodroma_heraldica,(Pterodroma_phaeopygia,Pterodroma_sandwichensis))),Pterodroma_alba))),(((Pterodroma_brevipes,(Pterodroma_leucoptera,(Pterodroma_cookii,Pterodroma_defilippiana))),(Pterodroma_pycrofti,Pterodroma_longirostris)),Pterodroma_hypoleuca)),Pterodroma_axillaris),Pterodroma_nigripennis)),(Pagodroma_nivea,(Thalassoica_antarctica,(Daption_capense,((Macronectes_giganteus,Macronectes_halli),(Fulmarus_glacialis,Fulmarus_glacialoides)))))),((Halobaena_caerulea,(Pachyptila_turtur,(Pachyptila_vittata,(Pachyptila_salvini,Pachyptila_desolata)))),(Bulweria_bulwerii,(Procellaria_westlandica,(Procellaria_aequinoctialis,Procellaria_parkinsoni))))))),(Pseudobulweria_aterrima,Pseudobulweria_rostrata)),(Puffinus_tenuirostris,(((Puffinus_griseus,(Puffinus_carneipes,Puffinus_creatopus)),Puffinus_gravis),(Puffinus_pacificus,Puffinus_bulleri)))),(Calonectris_diomedea,Calonectris_leucomelas))),(Puffinus_huttoni,Puffinus_gavia)),((Puffinus_puffinus,Puffinus_opisthomelas),(Puffinus_yelkouan,Puffinus_mauretanicus))),Puffinus_assimilis));";
  //string input_tree = "((((Procellaria_parkinsoni,(Puffinus_nativitatis,Procellaria_aequinoctialis)),((((((Eudyptes_chrysocome,(((Calonectris_diomedea,Pterodroma_externa),Thalassarche_chrysostoma),(Pterodroma_hasitata,Thalassarche_impavida))),(((((Thalassarche_bulleri,Oceanodroma_furcata),Pterodroma_feae),(Halocyptena_microsoma,Puffinus_creatopus)),((Diomedea_amsterdamensis,Puffinus_mauretanicus),Pelecanoides_urinatrix)),Oceanodroma_tristrami)),(Fulmarus_glacialis,Aptenodytes_patagonicus)),((((((Lugensa_brevirostris,Pterodroma_arminjoniana),((Pterodroma_hypoleuca,Pterodroma_incerta),Calonectris_leucomelas)),(((Pterodroma_alba,Puffinus_carneipes),(Pterodroma_phaeopygia,Pterodroma_solandri)),Pterodroma_sandwichensis)),(((Fregetta_tropica,(Pterodroma_leucoptera,(Pterodroma_mollis,Pygoscelis_antarctica))),((Puffinus_opisthomelas,Pseudobulweria_rostrata),Pagodroma_nivea)),(Oceanodroma_castro,Phoebastria_albatrus))),((Pterodroma_axillaris,Hydrobates_pelagicus),(Pelagodroma_marina,(Diomedea_gibsoni,Pygoscelis_adeliae)))),((Garrodia_nereis,Puffinus_tenuirostris),((Gavia_immer,Phoebastria_nigripes),((Pterodroma_cahow,Puffinus_bulleri),Thalassarche_melanophris))))),(((((Spheniscus_demersus,Thalassarche_bassi),((Macronectes_giganteus,(Pterodroma_defilippiana,Pachyptila_desolata)),Diomedea_sanfordi)),((((Pterodroma_heraldica,Diomedea_epomophora),Megadyptes_antipodes),((Bulweria_bulwerii,Pterodroma_brevipes),(Phoebastria_irrorata,Diomedea_exulans))),(Pygoscelis_papua,Phoebetria_palpebrata))),((Pterodroma_magentae,(Fulmarus_glacialoides,Oceanites_oceanicus)),(((Thalassarche_cauta,Eudyptes_pachyrhynchus),Pterodroma_cervicalis),(Thalassarche_salvini,Procellaria_westlandica)))),((Thalassarche_chlororhynchus,(Diomedea_dabbenena,Gavia_stellata)),((Pterodroma_nigripennis,Halobaena_caerulea),Puffinus_assimilis)))),((Puffinus_puffinus,Oceanodroma_leucorhoa),Pelecanoides_garnotii))),(((((((Puffinus_gavia,Daption_capense),(Fregetta_grallaria,Puffinus_auricularis)),(((Puffinus_griseus,(Pterodroma_baraui,Pseudobulweria_aterrima)),Diomedea_antipodensis),(Pachyptila_turtur,Pelecanoides_magellani))),(Phoebastria_immutabilis,((Oceanodroma_tethys,Thalassoica_antarctica),Pachyptila_salvini))),((Puffinus_huttoni,(Pterodroma_ultima,Pterodroma_cookii)),Puffinus_pacificus)),(((Macronectes_halli,(Pterodroma_madeira,Thalassarche_eremita)),(Pachyptila_vittata,(Pelecanoides_georgicus,Pterodroma_lessonii))),Pterodroma_longirostris)),((Oceanodroma_hornbyi,(Eudyptula_minor,Pterodroma_inexpectata)),Eudyptes_chrysolophus))),(((((Puffinus_yelkouan,Oceanodroma_melania),Procellaria_cinerea),(Pterodroma_neglecta,Phoebetria_fusca)),Pterodroma_pycrofti),(Puffinus_gravis,(Pterodroma_macroptera,Puffinus_lherminieri))));";

  Node* T = build_tree(suptree);
  adjustTree(T);

  Node* S = build_tree(input_tree);
  adjustTree(S);

  ///////////setting int labels for leaves, should be done in main()
  unordered_map<string, int> int_label_map;
  int starting_label = 1;
  find_int_labels_for_leaves_in_supertree(T, starting_label, int_label_map);  //finding int labels for all taxa in supertree
  set_int_labels_for_leaves_in_source_tree(S, int_label_map); // set int labels for taxa in source trees
  //for ( auto it = int_label_map.begin(); it != int_label_map.end(); ++it ) cout << it->first << ":" << it->second << endl;

  //since source trees won't change, for each node in each source tree, we can find clusters only ones, and store it in DS
  set_cluster_in_source_tree(S);







  //tttttttttttttttttttttttttttttttesting
  int best_score = 1000;
  for (int i = 0; i < 0; ++i)
  {
    NUM_SPR_NGHBRS = 0;
    start_time = clock();

    cout << "Suppose we want to solve SPR_RS for the following supertree (T) and given sourse tree(S):" << endl;
    cout << "S: " << S->str_subtree() << endl;
    cout << "T: " << T->str_subtree() << "\n" << endl;



    Node* best_node_to_prune;
    Node* best_node_to_regraft;
    bool weighted = false;
    int current_score = find_best_node_to_prune_and_its_best_regraft_place(*T, *S, best_node_to_prune, best_node_to_regraft, weighted);

    if (current_score < best_score) {

      best_score = current_score;


      int which_sibling = 0;
      Node* old_sibling = best_node_to_prune->spr(best_node_to_regraft, which_sibling);

      //perform best possible spr move of the neighbourhood
      /*
      int which_sibling = 0;
      if (! best_node_to_regraft->get_p()) { //if best node to regraft is root
        best_node_to_prune->spr(best_node_to_regraft, which_sibling);
        return *(best_node_to_prune->get_p());
      } else {
        best_node_to_prune->spr(best_node_to_regraft, which_sibling);
        return T;
      }
      */


      T = best_node_to_prune->find_root();
      adjustTree(T);


      finish_time = clock();
      float diff ((float)finish_time - (float)start_time);
      float seconds = diff / CLOCKS_PER_SEC;

      //cout << "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n" << endl;
      //cout << "source tree:" << endl;
      //cout << S->str_subtree() << endl;
      //cout << "intial ST:" << endl;
      //cout << T->str_subtree() << endl;
      //cout << "\n-------------------------  best ST found in SPR neighbourhood:   -------------------------------" << endl;
      cout << "T': " <<  T->str_subtree() << ";" << endl;
      //cout << "--------------------------------------------------------------------------------------------------" << endl;
      cout << "\nNUM_SPR_NGHBRS : " << NUM_SPR_NGHBRS << endl;
      cout << "the running time is: " << seconds << " sec." << endl;
      cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n\n" << endl;

    }
    else {
      break;
    }




  }



  cout << "-------------------------------------------------------------------------------------------" << endl;
  cout << "-------------------------------------------------------------------------------------------" << endl;
  cout << "-------------------------------------------------------------------------------------------" << endl;
  cout << "-------------------------------------------------------------------------------------------" << endl;
  cout << "Now let's see what is the result when WEIGHTED RF is used as optimality criteria: \n" << endl;

  Node* SS = build_tree(input_tree);
  adjustTree(SS);
  Node* TT = build_tree(suptree);
  adjustTree(TT);
  ///////////setting int labels for leaves, should be done in main()
  unordered_map<string, int> int_label_map1;
  int starting_label1 = 1;
  find_int_labels_for_leaves_in_supertree(TT, starting_label1, int_label_map1);  //finding int labels for all taxa in supertree
  set_int_labels_for_leaves_in_source_tree(SS, int_label_map1); // set int labels for taxa in source trees
  //for ( auto it = int_label_map.begin(); it != int_label_map.end(); ++it ) cout << it->first << ":" << it->second << endl;

  //since source trees won't change, for each node in each source tree, we can find clusters only ones, and store it in DS
  set_cluster_in_source_tree(SS);




  int weight = 2;
  int perc = 50;
  SS->reweight_edges_in_source_tree(perc, weight);

  Node* best_node_to_prune1;
  Node* best_node_to_regraft1;
  bool weighted1 = true;
  find_best_node_to_prune_and_its_best_regraft_place(*TT, *SS, best_node_to_prune1, best_node_to_regraft1, weighted1);
  int which_sibling1 = 0;
  Node* old_sibling1 = best_node_to_prune1->spr(best_node_to_regraft1, which_sibling1);



  finish_time = clock();
  float diff2 ((float)finish_time - (float)start_time);
  float  seconds = diff2 / CLOCKS_PER_SEC;

  cout << "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n" << endl;
  cout << "source_tree: " << input_tree << endl;
  cout << "intial   ST: " << suptree << endl;
  cout << "\n-------------------------  best ST found in SPR neighbourhood:   -------------------------------" << endl;
  cout << TT->str_subtree() << ";" << endl;
  cout << "--------------------------------------------------------------------------------------------------" << endl;
  cout << "\nNUM_SPR_NGHBRS : " << NUM_SPR_NGHBRS << endl;
  cout << "the running time is: " << seconds << " sec." << endl;
  cout << "\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n" << endl;

  print_weighted_tree(*SS);


  //don't forget to free memory!!
  T->delete_tree();
  S->delete_tree();

  return 0;
}

//tetsting
void preorder_trversal(Node& n) {
  cout << n.str_subtree() << endl;
  //cout << " prenum is : " << n.get_preorder_number() << endl;
  //cout << "lca mapping node is: " << n.get_lca_mapping().str_subtree() << "\n" << endl;
  cout << n.get_alpha() << "-" << n.get_beta() << "= " << n.get_alpha() - n.get_beta() << "\n" << endl;
  //cout << n.get_edge_weight() << "\n" << endl;

  //just for testing:
  list<Node *>::iterator c;
  list<Node *> children = n.get_children();
  for (c = children.begin(); c != children.end(); c++) {
    preorder_trversal(**c);
  }
}

//tetsting
void print_weighted_tree(Node& n) {
  cout << n.str_subtree() << endl;
  //cout << " prenum is : " << n.get_preorder_number() << endl;
  //cout << "lca mapping node is: " << n.get_lca_mapping().str_subtree() << "\n" << endl;
  //cout << n.get_alpha() << "-" << n.get_beta() << "= " << n.get_alpha() - n.get_beta() << "\n" << endl;
  cout << n.get_edge_weight() << "\n" << endl;

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
int find_best_node_to_prune_and_its_best_regraft_place(Node& T, Node& S, Node* & best_node_to_prune, Node* & best_node_to_regraft, bool weighted) {
  Node best_neighbor;

  int min_F = INT_MAX;
  Node* best_prune;
  Node* best_regraft;

  vector<Node*> internal_nodes;
  put_all_nodes_in_vector(T, internal_nodes);
  vector<Node*>::iterator iter, end;
  //any internal node except root may be pruned, iter plave the role of "v" in the algorithm
  for (iter = internal_nodes.begin(), end = internal_nodes.end(); iter != end;  iter++) {

    if (!(*iter)->get_p()) { //root or its children
      continue;
    }

    //cout << "\n-------------------------------------------------------" << endl;

    Node* best_regraft_place = apply_SPR_RS_algorithm_to_find_best_regraft_place(T, **iter, S, weighted);

    int which_sibling = 0;
    cout << "\n\n------T before: " << T.str_subtree() << endl;
    cout << "-----spr_on(v): " << (*iter)->str_subtree() << endl;
    cout << "--best regraft: " << best_regraft_place->str_subtree() << endl;
    Node* old_sibling = (*iter)->spr(best_regraft_place, which_sibling);
    adjustTree(&T);
    cout << "------T after : " << T.str_subtree() << endl;


    int current_F = 0 ;  //start from -1 instead of 0 since F is ("RF-dist"+1), That's because the edge connecting root to its children is counted twice (if binary)
    find_F_T(S, T, current_F, weighted);
    cout << "best F: " << min_F << ", and curent F: " << current_F << endl;
    if (current_F < min_F) {
      cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>> better SPR move with F: " << current_F << endl;
      //cout << "--best regraft: " << best_regraft_place->str_subtree() <<  endl;

      min_F = current_F;
      best_prune = *iter;
      best_regraft = best_regraft_place;
    }
    (*iter)->spr(old_sibling, which_sibling);
    //reset_alpha_beta(T);
    adjustTree(&T);

  }

  best_node_to_prune = best_prune;
  best_node_to_regraft = best_regraft;

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


//Actually finds RF distance!!
//assuming S has correct values for lca_mapping, finds the number of internal nodes whose corresponding lca does not exist in T
//Note the way I implemented find_best_node_to_prune_and_its_best_regraft_place(), I don't want to change T so that I have to do calculations (lca_mapping, ...) for each v to be pruned
//The node T being passed to this function is T', i.e. best neghibour for given v to be pruned.
//THUS the lca_mappings are NOT valid anymore, and that's why I compute it again here
void find_F_T(Node& S, Node& T_prime, int& num_clusters_not_in_T_prime, bool weighted) {
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
      find_F_T(**c, T_prime, num_clusters_not_in_T_prime, weighted);
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

      cout << ".....................Node u is: " << S.str_subtree() << endl;
      cout << "a, i.e. lca_mapping_in_T_prime: " << lca_mapping_in_T_prime->str_subtree() << endl;

      if ( (cluster.size()) != (lca_mapping_in_T_prime->number_of_leaves()) ) { //f(u)=0 , here S is actually u --->NOTE you CAN'T use get_cluster_size() cuz it will return cluster size of "this" in T not T'
        if (weighted) {
          cout << "weight of this bipartition: " << S.get_edge_weight() << endl;
          num_clusters_not_in_T_prime += S.get_edge_weight();
        } else {
          num_clusters_not_in_T_prime ++;
        }

      }
    }
  }

}


//finds node for which "alpha - beta" is maximized
//ASSUMES T is the root of Q, i.e. only containing VALID regraft places for v under examination
void find_best_regraft_place(Node& T, Node*& best_regraft_place, int& max) {
  NUM_SPR_NGHBRS++;
  int best = T.get_alpha() - T.get_beta();
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
Node* apply_SPR_RS_algorithm_to_find_best_regraft_place(Node & T, Node & v, Node & S, bool weighted) {

  Node* best_regraft_place = 0;

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


  //BUT if v's parent is T, then its a bad SPR and should be handled
  if (v.get_p() == &T) {
    //cout << "v is child of T, and should be handled differently.." << endl;

    Node* Q =  v.get_sibling();

    //cout << "In preprocessing step, the algorithm, initializes alpha, beta, cluster size, lca_mapping, and also constructs R from T as discussed in paper: " << endl;
    //cout << "R: " << T.str_subtree() << "; prenum in R is: " << T.get_preorder_number() << endl;
    //cout << "v: " << v.str_subtree() << "; prenum in R is: " << v.get_preorder_number() << endl;
    //cout << "Q: " << Q->str_subtree() << "; prenum in R is: " << Q->get_preorder_number() << "\n" << endl;

    compute_lca_mapping_S_R(S, T);  //No need for R here, T will do the work
    set_cluster_size_in_supertree(&T);

    //Find S' in lemma 12:
    Node* S_prime = build_tree(S.str_subtree() + ";"); //make a copy of S (not the best way to mak ea copy, I know, should be fixed)
    adjustTree(S_prime);
    S_prime->copy_fields_from_S(S);



    //cout << "-------------------S' before: " << S_prime->str_subtree() << endl;
    suppress_nodes_with_mapping_in_Rv(*S_prime, v);  //S' in lemma 12 is now  constructed from S
    //cout << "-------------------S' after : " << S_prime->str_subtree() << endl;


    traverse_S_and_update_alpha_beta_in_Q(*Q, v, S, *S_prime, weighted);

    //cout << "The final alpha beta values in nodes in Q after the algorithm finishes are as follow: " << endl;
    //preorder_trversal(*Q);

    int max_alpha_minus_beta = INT_MIN;
    find_best_regraft_place(*Q, best_regraft_place, max_alpha_minus_beta);  //NOTE Q should be passed

    S_prime->delete_tree();  //prevent memory leak

  } else { //normal case
    Node* R = new Node();
    R->add_child(&T);
    adjustTree(R); //because .. prenum of nodes in T has been changed (actually +1'ed)

    //cunstructing R which is SPR(v,rt(T))
    int which_sibling = 0;
    Node* old_sibling = v.spr(&T, which_sibling);
    adjustTree(R);//DON'T FORGET THIS!!!
    //cout << "R after SPR(v,T): " << R->str_subtree() << endl;

    //cout << "In preprocessing step, the algorithm, initializes alpha, beta, cluster size, lca_mapping, and also constructs R from T as discussed in paper: " << endl;
    //cout << "R: " << R->str_subtree() << "; prenum in R is: " << R->get_preorder_number() << endl;
    //cout << "v: " << v.str_subtree() << "; prenum in R is: " << v.get_preorder_number() << endl;
    //cout << "Q: " << T.str_subtree() << "; prenum in R is: " << T.get_preorder_number() << "\n" << endl;

    compute_lca_mapping_S_R(S, *R);
    set_cluster_size_in_supertree(R);

    //Find S' in lemma 12:
    Node* S_prime = build_tree(S.str_subtree() + ";"); //make a copy of S (not the best way to mak ea copy, I know, should be fixed)
    adjustTree(S_prime);
    S_prime->copy_fields_from_S(S);


    //cout << "-------------------S' before: " << S_prime->str_subtree() << endl;
    suppress_nodes_with_mapping_in_Rv(*S_prime, v);  //S' in lemma 12 is now  constructed from S
    //cout << "-------------------S' after : " << S_prime->str_subtree() << endl;


    //Note T now is Q in paper's notation
    traverse_S_and_update_alpha_beta_in_Q(T, v, S, *S_prime, weighted);

    //cout << "The final alpha beta values in nodes in Q after the algorithm finishes are as follow: " << endl;
    //preorder_trversal(T);

    int max_alpha_minus_beta = INT_MIN;
    find_best_regraft_place(T, best_regraft_place, max_alpha_minus_beta);  //NOTE here T is actually root of Q in apper notation

    v.spr(old_sibling, which_sibling);  // Note T has been changed here, and it should be retrieved to its original form
    T.cut_parent(); // after finding best place to regraft and make the corresponding SPR move, R is a node with only one child T, we don't need R anymore.
    //(v.get_p()) -> cut_parent(); // NOTE v is not a descendant of T anymore, they are SEPARATE trees
    R->delete_tree();
    S_prime->delete_tree();
  }

  return best_regraft_place;
}


//Given the way R is constructed, lca_mapping of any node in S, is either in T_v OR in T_T. (NOTE this is a recursive function and S plays the role of u in the alg in paper)
//The way I implemented here, avoids some duplicate work (in compare to when I pass R as parameter instead of T and v).
//Note T now is Q in paper's notation
void traverse_S_and_update_alpha_beta_in_Q(Node & T, Node & v, Node & S, Node & S_prime, bool weighted) { //note S here is u in paper's notation

  //only for internal nodes, i.e. non-root and non-leaf
  if (S.is_leaf()) {  //if leaf
    //do nothing

  } else if (S.get_p() == NULL) { //if root
    list<Node *>::iterator c;
    list<Node *> children = S.get_children();
    for (c = children.begin(); c != children.end(); c++) {
      traverse_S_and_update_alpha_beta_in_Q(T, v, **c, S_prime, weighted);
    }
  } else {  //if internal node

    //cout << "++++++++++++++++++++node u in S being considered: " << S.str_subtree() << endl;
    //cout << "its lca_mapping is: " << S.get_lca_mapping()->str_subtree() << endl;
    Node* a;

    //Suppose for each u in I(S), lca_mpping(u)=a
    //There are 4 cases.
    //1- if u satisfies precondition of Lemma 9, i.e. a belongs to V(R_v)
    if (a = v.find_by_prenum(S.get_lca_mapping() -> get_preorder_number())) {
      //cout << "lemma 9" << "\n" <<  endl;
      //do nothing
    }

    //2- if u satisfies precondition of Lemma 11, i.e. a belongs to V(Q_v) -which is T here- AND f_R(S)==0
    else if (a = T.find_by_prenum(S.get_lca_mapping()  -> get_preorder_number() )) {
      if (S.get_cluster_size() != a->get_cluster_size()) { //is f_R(S)==0
        //cout << "lemma 11" << "\n" <<  endl;
        //do nothing
      } else { //precondition of lemma 10 is true
        //3- if u satisfies precondition of Lemma 10, i.e. a belongs to V(T_v) AND f_R(S)==1
        //cout << "lemma 10" << "\n" <<  endl;
        if (weighted) {
          int w = S.get_edge_weight();
          a->increment_beta_in_all_descendants(w);
          //cout << "________w____increment beta on : " << a->str_subtree() << endl;
        } else {
          a->increment_beta_in_all_descendants(1);
          //cout << "____________increment beta on : " << a->str_subtree() << endl;
        }
      }
    }

    //Lemma 12's precondition: a=T->p() AND |L(R_b)|+|L(R_v)|=|L(S_u)|
    //Note: if S's lca_mapping is not found in T_v or T_T, then it's lca_mapping is definitely parent of T (or v) which is the only child of R
    else {
      //cout << "lemma 12" <<  endl;

      Node* u_in_S_prime = S_prime.find_by_prenum(S.get_preorder_number()); //find u (S in parameters) in S'
      Node* lca_mapping_lemma12;
      vector<int> cluster = u_in_S_prime->find_cluster_int_labels();  //note u's cluster field is not valid in S', and I should find it again
      bool f = false;
      compute_lca_mapping_helper_2(cluster, T.get_p(), lca_mapping_lemma12, f);
      Node* b_in_lemma12 = (T.get_p())->find_by_prenum(lca_mapping_lemma12  -> get_preorder_number() );


      if ( ((b_in_lemma12->find_leaves()).size()) + ((v.find_leaves()).size()) == (S.get_cluster_size()) ) { //|L(R_b)|+|L(R_v)|=?|L(S_u)|
        //cout << (b_in_lemma12->find_leaves()).size() << "+" << ((v.find_leaves()).size()) << " == " << (S.get_cluster_size()) << endl;
        if (weighted) {
          int w = S.get_edge_weight();
          b_in_lemma12->increment_alpha_in_all_descendants(w);
          //cout << "_____w_______increment alpha on : " << b_in_lemma12->str_subtree() << endl;
        } else {
          b_in_lemma12->increment_alpha_in_all_descendants(1);
          //cout << "____________increment alpha on b lemma 12: " << b_in_lemma12->str_subtree() << endl;
        }

      }

    }

    list<Node *>::iterator c;
    list<Node *> children = S.get_children();
    for (c = children.begin(); c != children.end(); c++) {
      traverse_S_and_update_alpha_beta_in_Q(T, v, **c, S_prime, weighted);
    }


  }
}

//pre-orderly traverses tree S, and suppresses all nodes in it for which lca_mapping is in T_v
void suppress_nodes_with_mapping_in_Rv(Node & S, Node & v) {
  if (v.find_by_prenum(S.get_lca_mapping()  -> get_preorder_number() )) {
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
void set_cluster_in_source_tree(Node * T) {
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
    set_cluster_in_source_tree(*c);
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



