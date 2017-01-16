
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
void write_line_to_File(string s,const char* file_name);
double calculate_rf_btwn_ST_n_source_tree(int argc, char** argv);
double calculate_total_rf(string source_trees[], set<string> non_shared_taxa[], string supertree, bool is_weighted);



void writeToFile(string s,const char* file_name);
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
void compute_lca_mapping_helper_2(vector<int>& cluster, Node* R, int& lca, bool& found);
void set_cluster_in_source_tree(Node* T);
void set_cluster_size_in_supertree(Node* T);
void reset_fields_to_initial_values(Node* T);
void apply_SPR_RS_algorithm_to_find_alpha_betas(Node& T, Node& spr_on, Node& S);
void find_best_spr_move_and_apply_it(Node& T, Node& spr_on, Node& S);
void traverse_S_and_update_alpha_beta_in_Q(Node& T, Node& v, Node& S);
void find_b_in_lemma12(Node& S, Node& R); // i'm not using this now
void suppress_nodes_with_mapping_in_Rv(Node& S_prime, Node& v);
void find_best_regraft_place(Node& n, Node*& best_regraft_place, int& max);
void preorder_trversal(Node& n);


template <typename T>
bool IsSubset(std::vector<T> A, std::vector<T> B)
{
    //std::sort(A.begin(), A.end());
    //std::sort(B.begin(), B.end());
    return std::includes(A.begin(), A.end(), B.begin(), B.end());
}


/////////////////////////////     main()     /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


int main(int argc, char** argv) {

    //running time
    clock_t start_time,finish_time;
    start_time= clock();

    //two trees with RF distance of 4
    string suptree = "((a,(b,g)),((c,f),(d,e)));";
    //string input_tree = "(c,(f,(e,(d,(g,(b,a))))));";
    string input_tree = "(c,(f,(e,(d,(g,(b,a))))));";

    //string suptree = "((t11,(t1,((t23,(((t3,t12),t4),(((t8,((t20,(((t26,t19),t30),t24)),t6)),t16),t17))),((t28,((t15,t25),((((t9,t18),t10),t5),(((t7,t22),t27),(t21,t2))))),(t13,t29))))),t14);"
    //string input_tree = "(t11,(t14,(t1,((t23,(((t3,t12),t4),(((t8,((t20,(((t26,t19),t30),t24)),t6)),t16),t17))),((t28,((t15,t25),((((t9,t18),t10),t5),(((t7,t22),t27),(t21,t2))))),(t13,t29))))));"

    Node* T= build_tree(suptree);
    adjustTree(T);

    Node* S= build_tree(input_tree);
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
    
    cout << "Suppose we want to solve SPR_RS for the following supertree (T) and given sourse tree(S):" << endl;
    cout << "S: " << S->str_subtree() << endl;
    cout << "T: " << T->str_subtree() << "\n" << endl;

    cout << "Now suppose, we want to prune node v with prenum=7 in T, which is the node: " << endl;    
    Node* v = T->find_by_prenum(6);
    cout << "v: " << v->str_subtree() << "\n" << endl;

    
    cout << "Now we call find_best_spr_move_and_apply_it() which first calls apply_SPR_RS_algorithm_to_find_alpha_betas(),\n" << 
      "and then finds the best place in Q for regrafting, and applies spr_move to T: " << endl;
    cout << "-------------------------------------------------------" << endl;
    find_best_spr_move_and_apply_it(*T, *v, *S);
    
    




    //don't forget to free memory!!
    T->delete_tree();
    S->delete_tree();


    finish_time=clock();
    float diff ((float)finish_time - (float)start_time);
    float seconds= diff / CLOCKS_PER_SEC;

    cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    cout << "the running time is: " << seconds << " sec." << endl; 


    return 0;
}

//tetsting
void preorder_trversal(Node& n) {
    cout << n.str_subtree() << endl;
    //cout << n.get_lca_mapping() << "\n" << endl;
    cout << n.get_alpha() << "-" << n.get_beta() << "= " << n.get_alpha() - n.get_beta() << "\n"<< endl;

    //just for testing:
    list<Node *>::iterator c;
    list<Node *> children= n.get_children();
    for(c = children.begin(); c != children.end(); c++) {
      preorder_trversal(**c);
    }
}



void find_best_spr_move_and_apply_it(Node& T, Node& spr_on, Node& S) {
    apply_SPR_RS_algorithm_to_find_alpha_betas(T, spr_on, S);

    int max_alpha_minus_beta = 0;
    Node* best_regraft_place;
    find_best_regraft_place(T, best_regraft_place, max_alpha_minus_beta);  //NOTE here T is actually root of Q in apper notation
    if(best_regraft_place) {
      int which_sibling=0;
      cout << "*****T before: " << T.str_subtree() << endl;
      cout << "*****spr_on(v): " << spr_on.str_subtree() << endl;
      cout << "*****new_sibling: " << best_regraft_place->str_subtree() << endl;
      spr_on.spr(best_regraft_place, which_sibling);
      cout << "*****T after : " << T.str_subtree() << endl;
    } else {
      //this means there is no better neighbour: LOCAL OPTIMUM!
    }

}

//fins pointer to node for which alpha-beta is maximized
void find_best_regraft_place(Node& T, Node*& best_regraft_place, int& max) {
    int best = T.get_alpha() - T.get_beta();
    if(best > max) { cout << "is better: " << T.str_subtree() << endl;
      max = best;
      best_regraft_place = &T;
    }

    //just for testing:
    list<Node *>::iterator c;
    list<Node *> children= T.get_children();
    for(c = children.begin(); c != children.end(); c++) {
      find_best_regraft_place(**c, best_regraft_place, max);
    }
}


//given trees T, S and node v on T to be pruned, finds, in O(n), the best place to regraft v so that RF-dist(T,S) is minimized,
//i.e. it solves SPR(v) for T (and S as source tree)
//it changes the T to its best SPR neighbor, too.
void apply_SPR_RS_algorithm_to_find_alpha_betas(Node& T, Node& v, Node& S) {

  //1- preprocessing step:
  ////////////////////////////////////////////////////////
    //we make tree R = SPR(v,rt(T)) from T with the following steps:
    //1- create a new node, R, and make rt(T) to be its child
    //2- Make the spr move as v->spr(rt(T),0), i.e. prune p(v) from tree and regraft it between rt(T) and R
    Node* R= new Node();
    R->add_child(&T);
    adjustTree(R); //because .. prenum of nodes in T has been changed (actually +1'ed)
    
    //cunstructing R which is SPR(v,rt(T))
    int which_sibling = 0;
    Node* old_sibling = v.spr(&T,which_sibling);
    adjustTree(R);//DON'T FORGET THIS!!!
    //cout << "R after SPR(v,T): " << R->str_subtree() << endl;
    
    compute_lca_mapping_S_R(S, *R);
    set_cluster_size_in_supertree(R);

    cout << "In preprocessing step, the algorithm, initializes alpha, beta, cluster size, lca_mapping, and also constructs R from T as discussed in paper: " << endl;
    cout << "R: " << R->str_subtree() << "; prenum is: " << R->get_preorder_number() << endl;
    cout << "v: " << v.str_subtree() << "; prenum is: " << v.get_preorder_number() << endl;
    cout << "Q: " << T.str_subtree() << "; prenum is: " << T.get_preorder_number() << "\n" << endl;

  //2- Calculating alpha and beta in Q:
  ///////////////////////////////////////////////////////////
    //The algorithm traverses through S, and updates alpha and beta values in Q (which is the parameter T). 
    //Note I pass T (rather than R) to the following function, since T repreents Q in the paper's notation, and we only care about alpha&beta vaues in Q.
    cout << "**********Main part of Algorithm****************" << endl;
    cout << "The algorithm traverses through S and updates alpha beta i each internal node in Q:" << "\n" << endl;

    traverse_S_and_update_alpha_beta_in_Q(T, v, S);

    cout << "The final alpha beta values in nodes in Q after the algorithm finishes are as follow: " << endl;
    preorder_trversal(T);


    v.spr(old_sibling, which_sibling);
    T.cut_parent(); // after finding best place to regraft and make the corresponding SPR move, R is a node with only one child T, we don't need R anymore. 
    //(v.get_p()) -> cut_parent(); // NOTE v is not a descendant of T anymore, they are SEPARATE trees
    R->delete_tree();
}


//Given the way R is constructed, lca_mapping of any node in S, is either in T_v OR in T_T. (NOTE this is a recursive function and S plays the role of u in the alg in paper)
//The way I implemented here, avoids some duplicate work (in compare to when I pass R as parameter instead of T and v).
void traverse_S_and_update_alpha_beta_in_Q(Node& T, Node& v, Node& S) {   //note S here is u in paper's notation


  if(S.is_leaf()) { 
    //do nothing
  }else {
      cout << "node in S being considered: " << S.str_subtree() << endl;
  cout << "its lca_mapping's prenum is: " << S.get_lca_mapping() << endl;
    Node* a;
    bool lemma10 = false;
    //Suppose for each u in I(S), lca_mpping(u)=a
    //There are 4 cases. 
    //1- if u satisfies precondition of Lemma 9, i.e. a belongs to V(R_v)
    if (a = v.find_by_prenum(S.get_lca_mapping())) { cout << "lemma 9" << "\n" <<  endl;
      //do nothing
    }

    //2- if u satisfies precondition of Lemma 11, i.e. a belongs to V(T_v) AND f_R(S)==0
    else if (a = T.find_by_prenum(S.get_lca_mapping())) {  cout << "lemma 11" << "\n" <<  endl;
      if (S.get_cluster_size() != a->get_cluster_size()) { //is f_R(S)==0
        //do nothing
      }else { //precondition of lemma 10 is true
        lemma10 = true;
      }
    } 

    //3- if u satisfies precondition of Lemma 10, i.e. a belongs to V(T_v) AND f_R(S)==1
    else if (lemma10) { cout << "lemma 10" << "\n" <<  endl;
      a->increment_by1_beta_in_all_descendants();
    } 

    //Lemma 12's precondition: a=T->p() AND |L(R_b)|+|L(R_v)|=|L(S_u)|
    //Note: if S's lca_mapping is not found in T_v or T_T, then it's lca_mapping is definitely parent of T (or v) which is the only child of R
    else { cout << "lemma 12 \n" <<  endl;
      //sooooo the tricky part is to find b here:
      Node S_prime = Node(S); //make a copy of S
      suppress_nodes_with_mapping_in_Rv(S_prime, v);  //S' is now  constructed from S
      Node* u_in_S_prime = S_prime.find_by_prenum(S.get_preorder_number()); //find u (S in parameters) in S'  
      int lca_mapping_lemma12;
      vector<int> cluster = u_in_S_prime->find_cluster_int_labels();  //note u's cluster field is not valid in S', and I should find it again
      bool f = false;
      compute_lca_mapping_helper_2(cluster, T.get_p(), lca_mapping_lemma12, f);
      Node* b_in_lemma12 = (T.get_p())->find_by_prenum(lca_mapping_lemma12);
      
      //S_prime.delete_tree(); // don't forget to free memory

      if ( ((b_in_lemma12->find_leaves()).size()) + ((v.find_leaves()).size()) == (S.get_cluster_size()) ) { //|L(R_b)|+|L(R_v)|=?|L(S_u)|
        b_in_lemma12->increment_by1_alpha_in_all_descendants();
      }

    }



    list<Node *>::iterator c;
    list<Node *> children= S.get_children();
    for(c = children.begin(); c != children.end(); c++) {
      traverse_S_and_update_alpha_beta_in_Q(T, v, **c);
    }
  } 
}

//pre-orderly traverses tree S, and suppresses all nodes in it for which lca_mapping is in T_v
void suppress_nodes_with_mapping_in_Rv(Node& S, Node& v) {
    if(v.find_by_prenum(S.get_lca_mapping())) {
      //cout << "remove node : " << S.str_subtree() << endl;
      Node* p = S.get_p();
      p -> delete_child(&S);
      S.delete_tree();  //prevent memory leak!
      if(p -> get_children().size() == 1) {
        Node* pp = p->get_p();
        p -> contract_node();
      }
    } else {
      list<Node *>::iterator c;
      list<Node *> children= S.get_children();
      for(c = children.begin(); c != children.end(); c++) {
        suppress_nodes_with_mapping_in_Rv(**c, v);
      }
    }


}

//Why do we need all the following 3 parameters:
//S: we need S to make a copy of it and make S' from it (now all prenum, cluster, ... are invalid)
//u: u is the node for which I know I will need to know b_in_lemma12 field (u exists in both S and S')
//R: b_in_lemma12 is the mapping of u in R
//v: cuz S' is constructed from S based on v
void find_b_in_lemma12(Node& u, Node& S, Node& R, Node& v) {
  int prenum_of_u = u.get_preorder_number();
  Node S_prime = Node(S); //make a copy of S

}


//for all u in I(S), find the LCA mapping of u:
//1-Let's cluster(u) be set of leaves in subtree rooted u in tree S,
//2-Find LCA of cluster(u) in tree R, and set its lca_mapping
//NOTE: lca_mapping is the prenum f the LCA node in R
void compute_lca_mapping_S_R(Node& S, Node& R) {
    if(S.is_leaf()) {
      int lca_mapping;
      vector<int> cluster = S.get_cluster();
      bool f = false;
      compute_lca_mapping_helper_2(cluster, &R, lca_mapping, f);
      S.set_lca_mapping(lca_mapping);
    }else {
      int lca_mapping;
      vector<int> cluster = S.get_cluster();
      bool f = false;
      compute_lca_mapping_helper_2(cluster, &R, lca_mapping, f);
      S.set_lca_mapping(lca_mapping);
      //recursively set children's lca_mapping
      list<Node *>::iterator c;
      list<Node *> children= S.get_children();
      for(c = children.begin(); c != children.end(); c++) {
        compute_lca_mapping_S_R(**c, R);
      }
    }
}

void compute_lca_mapping_helper_1(vector<int>& cluster, Node* R, int& lca) {
    bool lca_must_be_in_lower_levels = false;
    Node* which_child;

    list<Node *>::iterator c;
    list<Node *> children= R->get_children();
    for(c = children.begin(); c != children.end(); c++) {
      vector<int> c_cluster = (*c)-> find_cluster_int_labels();
      if (includes(c_cluster.begin(), c_cluster.end(), cluster.begin(), cluster.end())) {  //if cluster is a subset of c_cluster
        lca_must_be_in_lower_levels = true;
        which_child = *c;
        break;
      }
    }

    if(lca_must_be_in_lower_levels) {
      compute_lca_mapping_helper_1(cluster, which_child, lca);
    } else {
      lca = R -> get_preorder_number();
    }
}

//assuming "cluster" is the labels for which we are lookig for the LCA, the function finds the LCA,
//and assigns its prenum to parameter "lca"
void compute_lca_mapping_helper_2(vector<int>& cluster, Node* R, int& lca, bool& found) {
    if(! R->is_leaf()){
      list<Node *>::iterator c;
      list<Node *> children= R->get_children();
      for(c = children.begin(); c != children.end(); c++) {
      compute_lca_mapping_helper_2(cluster, *c, lca, found);
      }
    }

    if(found) { // to return from nested recursive calls after 
      return;
    }

    if(R-> get_lca_hlpr() == cluster.size()) {
      lca = R->get_preorder_number();
      found = true;
      R->set_lca_hlpr(0); //to avoid for reseting lca_hlpr values to 0 for the next time we compute lca_mapping
      return;
    }

    if(R->is_leaf()) {
      if(cluster.size() == 1) { //leaves can have lca as well
        if (R->get_int_label() == cluster[0]) {
          lca = R->get_preorder_number();
          found = true;
          return;
        }
      }else {
        //assuming "cluster" is SORTED!!
        if(binary_search(cluster.begin(), cluster.end(), R->get_int_label())) {
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
void find_int_labels_for_leaves_in_supertree(Node* root, int& current_label, unordered_map<string, int>& map) {
    if (root->is_leaf()) {
        map.insert(pair<string, int>(root->get_name(), current_label));
        root->set_int_label(current_label);
        current_label++;
    }else {
      for (auto child : root->get_children()) {
        find_int_labels_for_leaves_in_supertree(child, current_label, map);
      }
    }  
}

//set int labels for leaves in one source tree
void set_int_labels_for_leaves_in_source_tree(Node* root, unordered_map<string, int>& map) {
    if (root->is_leaf()) {
        root->set_int_label(map.find(root->get_name())->second);
    }else {
      for (auto child : root->get_children()) {
        set_int_labels_for_leaves_in_source_tree(child, map);
      }
    }  
}

//finds clusters associated with each node in each SOURCE tree, and stores it in DS to avoid repeated work
//this method should be called before calling find_lca()
//note a node's cluster is a SORTED vector contining int_labels of the leaves in subtree induced at that node
void set_cluster_in_source_tree(Node* T) {
    if (T-> is_leaf()) {
      std::vector<int> v {T->get_int_label()};
      T-> set_cluster(v);
      T-> set_cluster_size(1);
    }else {
      T-> set_cluster(T->find_cluster_int_labels());
      T-> set_cluster_size((T->get_cluster()).size());
    }

    list<Node *>::iterator c;
    list<Node *> children= T->get_children();
    for(c = children.begin(); c != children.end(); c++) {
      set_cluster_in_source_tree(*c);
    }
}

//find and set cluster_size for each node in the supertree
void set_cluster_size_in_supertree(Node* T) {
    if (T-> is_leaf()) {
      T-> set_cluster_size(1);
    }else {
      T-> set_cluster_size((T->find_leaves()).size());
    }

    list<Node *>::iterator c;
    list<Node *> children= T->get_children();
    for(c = children.begin(); c != children.end(); c++) {
      set_cluster_size_in_supertree(*c);
    }
}

//reset invalid fields in ST after each iteration to initial values
void reset_fields_to_initial_values(Node* T) {
    T->set_cluster_size(0);
    //T->set_lca_hlpr(0);
    T->set_alpha(0);
    T->set_beta(0);

    list<Node *>::iterator c;
    list<Node *> children= T->get_children();
    for(c = children.begin(); c != children.end(); c++) {
      reset_fields_to_initial_values(*c);
    }    
}




//*******************added for SPR neighborhood*********************
//******************************************************************
//applies all the possible spr's on "tree", and writes them into the file "z_spr_neighbours"
int produce_all_spr_neighbors(Node* myTree, int number_of_taxa) {

    int number_of_neighbors= 0;


    //The following line deletes the contents of z_spr_neighbours
    const char* neighbors_file = "z_spr_neighbours";
    std::ofstream ofile(neighbors_file, ios_base::trunc);

    for(int i=1; i<number_of_taxa; i++) {

        Node* spr_on= myTree->find_by_prenum(i);
        int which_sibling=0;
        /*
        cout << "************************************************************" << endl;
        cout << "spr_on's pre-order number = i = " << spr_on->get_preorder_number() << endl;
        cout << "new_sibling's pre-order number = j = " << spr_on->get_preorder_number() << endl;
        cout << "original  tree: " << myTree->str_subtree() << endl;
        cout << "************************************************************" << endl;
        */

        for(int j=1; j<number_of_taxa; j++) {

            if (j != i) {

                Node* new_sibling= myTree->find_by_prenum(j);
                bool good_spr= true;
                //cout<< "i = " <<  i << ", j = " << j << endl;



                //bad spr: check whether new_sibling is parent
                Node* parent= spr_on->get_p();
                if (parent->get_preorder_number() == j) {
                    //cout << "-----BAD SPR WAS IGNORED----- \n" << endl;
                    continue;
                }


                //bad spr: check whether new_sibling is a descendant of spr_on
                if(! spr_on->is_leaf()) {
                    vector<Node *> descendants = spr_on->find_descendants();
                    vector<Node *>::iterator it;
                    for(it = descendants.begin(); it!= descendants.end(); ++it) {
                        if(new_sibling->get_preorder_number() == (*it)->get_preorder_number()) {
                            //cout << "-----BAD SPR WAS IGNORED----- \n" << endl;
                            good_spr= false;
                            continue;       //goes out of the inner loop which is 4 lines above
                        }
                    }
                }


                

                //bad spr: if new sibling is old sibling ignore it cuz this spr-move results to the original tree
                list<Node *>::iterator itr;
                for(itr = (spr_on->parent())->get_children().begin(); itr!= (spr_on->parent())->get_children().end(); ++itr) {
                    if((*itr)->get_preorder_number() == j) {
                        good_spr= false;
                        continue;
                    }
                }


                if(good_spr) {
                  Node* undo= spr_on->spr(new_sibling, which_sibling);
                  adjustTree(myTree);
                  number_of_neighbors ++;
                  string new_tree= myTree ->str_subtree();
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
void write_line_to_File(string s,const char* file_name) {
    ofstream myFile;
    myFile.open(file_name, std::ios_base::app);
    myFile << s + ";\n";
    myFile.close();
}


void adjustTree(Node* myTree) {
    myTree->set_depth(0);
    myTree->fix_depths();
    myTree->preorder_number();
    myTree->edge_preorder_interval();
}




//returns the total number of nodes in the clade of "node", i.e. (#of taxa)+(#of internal nodes)
//This method is a modification of the set_preorder_number() method in node.h
void total_number_of_nodes(Node* node, int& total_nodes) {
    total_nodes++;
    list<Node *>::iterator c;
    list<Node *> children= node->get_children();
    for(c = children.begin(); c != children.end(); c++) {
        total_number_of_nodes(*c, total_nodes);
    }
}



