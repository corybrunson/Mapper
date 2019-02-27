// Simple but limited implementation of the Simplex tree data structure using Rcpp + STL
// Original Reference: Boissonnat, Jean-Daniel, and Clement Maria. "The simplex tree: 
// An efficient data structure for general simplicial complexes." Algorithmica 70.3 (2014): 406-427.

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

#include <map>
#include <unordered_map>
#include <queue>
#include <vector>
#include <memory>
#include "utility_rcpp.h"

typedef std::size_t idx_t;
using std::vector;

// Node structure stored by the simplex tree. Contains the following fields:
//  label := unsigned integer representing the id of simplex it represents
//  parent := node pointer to its parent in the trie
//  children := connected simplexes whose labels > the current simplex's label
struct node {
  using node_ptr = s_ptr< node >; 
  idx_t label;
  s_ptr< node > parent;
  std::map< idx_t, node_ptr > children;
  node(idx_t id, node_ptr c_parent) : label(id), parent(c_parent){ 
    children = std::map< idx_t, node_ptr >(); 
  }
  node(const node& other) : label(other.label), parent(other.parent), children(other.children) { }
};
using node_ptr = s_ptr< node >; 

// Forward declarations 
class dfs_iter;
class bfs_iter;

// Simplex tree data structure. 
struct SimplexTree {
  using node_ptr = s_ptr< node >; 
  node_ptr root; // empty face; initialized to id = 0, parent = nullptr
  std::unordered_map< std::string, std::vector< node_ptr > > level_map; // maps strings of the form "<id>-<depth>" to a vector of node pointers
  std::vector<idx_t> n_simplexes;
    
  // Attach friend class iterators  
  friend class dfs_iter;
  friend class bfs_iter;
    
  // Constructor + Destructor 
  SimplexTree();
  ~SimplexTree();
  
  // Utilities
  SEXP as_XPtr();
  void record_new_simplexes(const idx_t k, const idx_t n);// record keeping
  
  // User-facing API 
  vector< idx_t > to_vec(IntegerVector);
  IntegerVector to_ivec(vector< idx_t >);
  int find_vertex(const int v_id);
  vector< idx_t > vertex_available(idx_t n_vertices);
  vector< idx_t > adjacent_vertices(const int v);
  void add_vertices(const idx_t v_i);
  void remove_vertices(vector< idx_t > vertex_ids);
  void remove_edge(vector< idx_t > labels);
  void remove_vertex_cofaces(const int v);
  void insert_simplex(std::vector<idx_t> labels);
  bool find_simplex(const vector< idx_t >& simplex);
  // void remove_simplex(const IntegerVector& simplex);
  void expansion(const idx_t k);
  void print_tree();
  void print_cofaces(int depth);
  IntegerVector get_vertices();

  // Export utilities
  IntegerMatrix as_adjacency_matrix(); // Exports the 1-skeleton as an adjacency matrix 
  List as_adjacency_list(); // Exports the 1-skeleton as an adjacency matrix 
  IntegerMatrix as_edge_list(); // Exports the 1-skeleton as an edgelist 
  
  // Recursive helper functions
  node_ptr insert_child(node_ptr c_parent, node_ptr new_child, idx_t depth);
  node_ptr remove(idx_t* labels, const size_t i, const size_t n_keys, node_ptr c_node, const idx_t depth);
  void remove_child(node_ptr c_parent, idx_t child_label, idx_t depth);
  void add_children(node_ptr c_parent, const std::vector<idx_t>& new_children, idx_t depth);
  void insert(idx_t* labels, const size_t i, const size_t n_keys, node_ptr c_node, const idx_t depth);
  node_ptr find (idx_t label);
  node_ptr find (vector<idx_t> simplex);
  // void as_list_helper(List& res, node_ptr c_node, idx_t depth, idx_t* simplex, const size_t n_keys);
  
  // Utility 
  idx_t get_dfs_height(node_ptr cnode, idx_t c_height);
  void print_level(node_ptr cnode, idx_t level);
  std::vector<idx_t> getLabels(const std::map<idx_t, node_ptr>& level, const idx_t offset = 0);
  idx_t intersection_size(std::vector<idx_t> v1, std::vector<idx_t> v2);
  std::vector<idx_t> intersection(std::vector<idx_t> v1, std::vector<idx_t> v2);
  void expand(std::map<idx_t, node_ptr>& v, const idx_t k, idx_t depth, idx_t* simplex);
  
  // Standrad typedefs for iterating   
  typedef std::ptrdiff_t difference_type;
  typedef std::size_t size_type;
  typedef node_ptr value_type;
};

#include <stack>
class dfs_iter {
private: 
  using node_ptr = s_ptr< node >; 
  size_t size;
  SimplexTree& stree; 
  std::stack< node_ptr > node_stack; 
  node_ptr current, onePastEnd; 
public: 
  dfs_iter(SimplexTree& stree_) : stree(stree_){
    size = std::accumulate(stree.n_simplexes.begin(), stree.n_simplexes.end(), 0);
    node_stack = std::stack< node_ptr >();
    node_stack.emplace(stree.root);
    current = stree.root;
  }  
  // void dfs_r(node_ptr c_node, size_t depth){
  //   if (c_node->children.empty()){ return; } 
  //   else {
  //     std::for_each(c_node->children.begin(), c_node->children.end(), [=](std::pair<idx_t, node_ptr>& node_pair){
  //       dfs_r(node_pair.second, depth + 1);
  //     });
  //     node_stack.emplace(c_node);
  //   }
  // }
  bool operator==(const dfs_iter& other){ return(this->current == other.current); } const
  bool operator!=(const dfs_iter& other){ return(!(this->current == other.current)); } const 
  node_ptr& operator*(){ return(current); } const 
  dfs_iter& operator++(){ // prefix
    // current++;
    if (node_stack.empty()){ return *this; }
    else {
      current = node_stack.top();
      node_stack.pop();
      if (!current->children.empty()){
        for (auto& nv: current->children){ node_stack.emplace(nv.second); }
      }
    }
    return *this; 
  } 
  // dfs_iter operator++ (int) { // postfix operator
  //   dfs_iter n = node_stack.at(offset);
  //   ++offset;
  //   return clone;
  // }

};


