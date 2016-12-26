
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


/////////////////////////////     main()     /////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


int main(int argc, char** argv) {

    //two trees with RF distance of 4
    string suptree = "((a,(b,g)),((c,f),(d,e)));";
    string input_tree = "(a,(b,(c,(d,(e,(f,g))))));";

    Node* supertree= build_tree(suptree);
    adjustTree(supertree);

    Node* source_tree= build_tree(input_tree);
    adjustTree(source_tree);

    //cout << supertree->str_subtree() << endl;

    for (int i= 0; i < 13; i++){
      cout << i << "- " ;
      Node* n = supertree->find_by_prenum(i);
      cout << n->str_subtree() << endl;      
    }

    //suppose we want to solve SPR(v) on supertree, where v is node with prenum 10, i.e. to find best place to regraft node v
    Node* st_root;




    return 0;
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



