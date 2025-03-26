#ifndef TREE_LD_H
#define TREE_LD_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <deque>
#include <tuple>
#include <random>
#include <cmath>
#include <limits>
#include <utility>
#include <algorithm>
#include <string>
#include <type_traits>
#include <random>
#include <filesystem> // FG 24/02/25 -- in order to make it compatible with WSL
/////////////////////////////
#include "../random/random.h"
#include "XoshiroCpp.hpp"

//gnu gsl library
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>


// Define SEED and S_DIR if not passed as macros
std::string curr_path = std::filesystem::current_path(); // FG 24/02/25 Adding absolute path
const std::string SEED = curr_path + "/../seed"; 
const std::string S_DIR = curr_path + "/../output";

//const std::string SEED = curr_path + "../seed"; 
//const std::string S_DIR = curr_path + "../output";


/*
Implementation of selection mechanism. First attempt 27/02/2025.
Idea: consider death_prob as node feature, not tree feature. 
This idea implies:
* Adapt struct Node
* Adapt tree so that death_prob is no longer a tree variable
* Adapt evolution/mutation function so that single the process of death happens node-wise.
Actually, the decision whether the daughter cells survive is made node wise. 
The point is properly implementing the death prob as a distribution or anything else. 
28/02/25
We try another way, by setting the death_probability as a node feature, 
but having it extracted from a pdf whose parameters (type, notable parameters) are tree features
*/

template <typename T> struct Node {
  T data;
  T parent;
  double death_probability;
  Node () {}
  Node(const T &data, const T &parent, double &death_probability) : data{data}, parent{parent}, death_probability{death_probability} {}
  Node& operator =(const Node &node)
  {
    data = node.data;
    parent = node.parent;
    death_probability = node.death_probability;
    return *this;
  }
};
/* // Old version
template <typename T> struct Node {
  T data;
  T parent;
  Node () {}
  Node(const T &data, const T &parent) : data{data}, parent{parent} {}
  Node& operator =(const Node &node)
  {
    data = node.data;
    parent = node.parent;
    return *this;
  }
};
*/
template <typename T> 
class Tree
{
public:
  Random rnd;
  std::vector<double> freq_mutations;

  Tree(const std::vector<Node<T>> &tree, 
       char death_prob_type,
       double death_prob_p0,
       double death_prob_p1,
       double seq_error,
       double K,
       size_t max_generations,
       size_t n_bases,
       size_t coverage,
       size_t ploidy,
       size_t n_zero):
       rnd(SEED + "/Primes", SEED + "/seed.in"),
       tree{tree},
       death_prob_type{death_prob_type},
       death_prob_p0{death_prob_p0}, 
       death_prob_p1{death_prob_p1},  
       seq_error{seq_error},
       K{K},
       max_generations{max_generations},
       n_bases{n_bases},
       coverage{coverage},
       ploidy{ploidy},
       n_zero{n_zero}
  {

    //passed to start_tree 
    mutation_rate = 0.0;
    min_survived_nodes_number = 0;
    threshold = 0.0;
    ////////////////////////////
    max_node_index = 0; 
    survived_nodes = 0;
    n_mutations = 0;
    count_mutations = 0;
    count_thr_mut = 0;
    all_alive = 0;
    n_mut_prior = 0;
    threshold = 0.0;
    cov_over_p = 0.0;
    ploidy_p = 0;
  }

  //get
  /////////////// // questi sono tutti metodi per accedere alle pvt???
  size_t get_max_node_index() const {return max_node_index;}
  
  std::vector<Node<T>> get_last_layer() const
  {
    return generations.back();
  }

  size_t get_generations_size() const
  {
    return generations.size();
  }

  size_t get_survived_nodes() const
  {
    return survived_nodes;
  }

  size_t get_all_alive() const
  {
    return all_alive;
  }

  char get_death_prob_type() const
  {
    return death_prob_type;
  }

  size_t get_not_mutated_bases(size_t s_zero) const
  {
   size_t n_muts = n_bases;
   if (!freq_mutations.empty())
   {
     n_muts = n_bases - freq_mutations.size() + s_zero;
   } 
   return n_muts;
  }
  ///////////////
  void set_K(double k_)
  {
    K = k_;
  }

  void set_nbases(size_t n_bases)
  {
    n_bases = n_bases;
  }

  void warm_up_random()
  {
    for (size_t i = 0; i < 1000000; ++i)
    {
      rnd.rannyu();
    }
  }

  void set_engine(std::mt19937 &engine)
  {
    engine = engine; 
  }

  void set_variables()
  {
    n_mutations = 0;
    count_mutations = 0;
    all_alive = 0;
    freq_mutations.resize(0);
  }

  void start_tree(long double mut, size_t n_e, size_t gen_i, std::mt19937 &engine)
  {  
    resize_tree_generations();
    max_node_index = 0;
    survived_nodes = 0;
    n_mutations = 0;
    count_mutations = 0;
    count_thr_mut = 0;
    all_alive = 0;
    n_mut_prior = 0;
    ploidy_p = 0;

    mutation_rate = mut;
    min_survived_nodes_number = n_e;
    max_generations = gen_i;


    // Now death_probability is a node feature. We implement this check in main.
    // Alternative/supplementary hypothesis: raise an exception while creating the node (not so easy to build)
    /*if (death_probability > 0.5) 
    {
      std::cout << std::endl;
      std::cout << "Value of death probability not allowed. " << std::endl;
      std::cout << "Tree extinction with probability = 1.   " << std::endl;
      std::cout << "Death probability needs to be <= 0.5.   " << std::endl;
      std::cout << std::endl;
      exit(1);
    }*/
    if (min_survived_nodes_number > 0 && max_generations == 0)
    {
      std::vector<Node<T>> root_layer;
      generations.push_back(root_layer);
      root_node(engine);
      generate_tree_min_size(mutation_rate, min_survived_nodes_number, max_generations, engine);
    }
    else if (min_survived_nodes_number == 0 && max_generations > 0) 
    {
      std::vector<Node<T>> root_layer;
      generations.push_back(root_layer);
      root_node(engine);
      generate_tree_fixed_generations(mutation_rate, min_survived_nodes_number, max_generations, engine);
    }
    else
    {
      std::cout << "Min_survived_nodes_number and max_generations both" << "\n";
      std::cout << "greater or equal to zero are forbidden. Exit." << std::endl;
      exit(1);
    }
  }

  void resize_tree_generations()
  {
    tree.resize(0);
    generations.resize(0); 
    freq_mutations.resize(0);
  }

  
  void root_node(std::mt19937 &engine) 
  {
    //std::cout << "Test" << std::endl;
    double death_prob = generate_death_prob(death_prob_type, death_prob_p0, death_prob_p1, engine);
    //std::cout << "Root node death probability = " << death_prob << std::endl;
    tree.push_back(Node<size_t>(0,0,death_prob));
    generations[0].push_back(Node<size_t>(0,0,death_prob));
  }

  void generate_tree_min_size(long double &mutation_rate, size_t &n_ext, size_t &gen_i, std::mt19937 &engine)
  {
    size_t tree_gen = 0;
    size_t old_max_node_index = 0;
    size_t incremental_index = 0;
    while (survived_nodes < min_survived_nodes_number)
    {
      add_layer(old_max_node_index, incremental_index, tree_gen, engine);
      if (((incremental_index + 1) - old_max_node_index) == 0)
      {
        // Questo passaggio non mi è particolarmente chiaro...
        //std::cout << "Tree is extinct. Tree regeneration... " << std::endl;
        start_tree(mutation_rate, n_ext, gen_i, engine);
        //std::cout << "Tree is extinct. Tree generation aborted." << std::endl;
        //break;
      }

      if (survived_nodes < min_survived_nodes_number)
      {
        survived_nodes = 0;
      }
    }
  }

  void generate_tree_fixed_generations(long double &mutation_rate, size_t &n_ext, size_t &gen_i, std::mt19937 &engine)
  {
    size_t tree_gen = 0;
    size_t old_max_node_index = 0;
    size_t incremental_index = 0;
    while (tree_gen < max_generations)
    {    
      add_layer(old_max_node_index, incremental_index, tree_gen, engine);

      if (((incremental_index + 1) - old_max_node_index) == 0)
      {
        break;
      }

      if (tree_gen < max_generations)
      {
        survived_nodes = 0;
      }
    }  
  }

  void add_layer(size_t &old_max_node_index, size_t &incremental_index, size_t &tree_gen, std::mt19937 &engine)
  {
    for (; old_max_node_index < incremental_index + 1; ++old_max_node_index)
    { 
      childs_nodes(tree[old_max_node_index], tree_gen, engine);
    } 
    
    incremental_index = max_node_index;
    old_max_node_index = (incremental_index + 1) - survived_nodes;
  }
  // Old version
  /*void childs_nodes(const Node<T> &node, size_t &tree_gen)
  { 
    Node temp = node;
    for (size_t i = 0; i < 2; ++i) // binary division. From one parent cell, 0, 1, or 2 daugther cells can result as outcome.
    {
      if ((rnd.rannyu() > death_probability) || death_probability == 0.0) // now the death probability is retrieved from the node
      {
        max_node_index += 1; 
        survived_nodes += 1;
        tree.push_back(Node<size_t>(max_node_index, temp.data)); 
        tree_generations(Node<size_t>(max_node_index, temp.data), tree_gen);
      }
    } 
  }*/

  // The survival decision is made at the node level. 
  // This one should work "as is"
  // What's left now is how to set the death_prob in the node
  void childs_nodes(const Node<T> &node, size_t &tree_gen, std::mt19937 &engine)
  { 
    Node temp = node;
    double new_node_death_prob;
    for (size_t i = 0; i < 2; ++i) // binary division. From one parent cell, 0, 1, or 2 daugther cells can result as outcome.
    {
      if ((rnd.rannyu() > temp.death_probability) || temp.death_probability == 0.0) // now the death probability is retrieved from the node
      {
        max_node_index += 1; 
        survived_nodes += 1;
        new_node_death_prob = generate_death_prob(death_prob_type, death_prob_p0, death_prob_p1, engine);
        //std::cout << "Node " << max_node_index << " has death probability = " << new_node_death_prob << std::endl;
        tree.push_back(Node<size_t>(max_node_index, temp.data, new_node_death_prob)); 
        tree_generations(Node<size_t>(max_node_index, temp.data, new_node_death_prob), tree_gen);
      }
    } 
  }

  // Wrap around different functions using different pdfs to generate death probabilities
  double generate_death_prob (char type, double p0, double p1, std::mt19937 &engine) 
  {
    double death_prob = 0;
    if (type == 'C')
    {
      death_prob = generate_death_prob_constant(p0);
    }
    else if (type == 'U')
    {
      death_prob = generate_death_prob_uniform(p0, p1, engine);
    }
    else
    {
      std::cerr << "Error: Death probability type not recognized. Exit." << std::endl;
      exit(1);
    }
    return death_prob;
  }
  
  // Single value extracted from a uniform pdf
  // Which raises the idea that the death_prob should be taken by the tree in terms of parameters: type, p0, p1
  double generate_death_prob_uniform(double lower_bound, double upper_bound, std::mt19937 &engine) // I'm not sure this is proper sintax
  {
    size_t attempts = max_node_index; 
    std::uniform_real_distribution<double> dis(lower_bound, upper_bound);
    double death_prob = dis(engine);
    return death_prob;
  }

  double generate_death_prob_constant(double p0)
  {
    return p0;
  }


  void fill_tree(const Node<T> &node)
  {
    tree.emplace_back(node);
    if (node.data != 0)
    {
      tree_generations(node);
    }
    else
    {
      generations[node.data].push_back(node);    
    }
  }

  void tree_generations(const Node<T> &node, size_t &tree_gen)
  {
    Node p = node; 
    size_t gen = 0;
    do
    {
      gen += 1;
      p = tree[p.parent];
    }
    while (p.data != 0); 

    if (min_survived_nodes_number > 0)
    {
      if (gen < generations.size())
      {
        generations[gen].push_back(node); 
      }
      else
      {
        std::vector<Node<T>> layer;
        generations.push_back(layer);
        generations[gen].push_back(node);
      }
    }
    else
    {
      if (gen < generations.size())
      {
        generations[gen].push_back(node); 
      }
      else
      {
        std::vector<Node<T>> layer;
        generations.push_back(layer);
        generations[gen].push_back(node);
        tree_gen += 1;
      }
    }
  }
 
  std::string find_cell_in_generation(T cell_id, std::vector<Node<T>> &generation)
  {
    std::string ret_idx = "NONE";
    for (auto it = generation.begin(); it != generation.end(); ++it)
    {
      if ((*it).data == cell_id)
      {
        //std::cout << "CellID = " << (*it).data << "\t" << " Parent = " << (*it).parent << std::endl;
        ret_idx = std::to_string(std::distance(generation.begin(), it));
        return ret_idx;
        break;
      }
    }

    return ret_idx;
  }

  std::tuple<std::string, std::string> find_cell_in_tree(T cell_id)
  {
    std::tuple<std::string, std::string> cell = std::make_tuple("NONE", "NONE");
    for (size_t i = 0; i < generations.size(); ++i)
    {
      if ((cell_id >= (generations[i].front()).data) && (cell_id <= (generations[i].back()).data))
      {
        std::string cell_idx = find_cell_in_generation(cell_id, generations[i]); 
        if (cell_idx != "NONE")
        {
          std::get<0>(cell) = std::to_string(i);
          std::get<1>(cell) = cell_idx;
          return cell;
          break;
        }
        break;
      }
    }
    return cell;
  }

  std::vector<T> descendant_ids(T cell_id)
  {
    std::vector<T> ids;
    std::tuple<std::string, std::string> cell_idxs = find_cell_in_tree(cell_id);
    if (std::get<0>(cell_idxs) != "NONE")
    {
      size_t row = static_cast<size_t>(std::stoi(std::get<0>(cell_idxs)));
      size_t col = static_cast<size_t>(std::stoi(std::get<1>(cell_idxs)));
      std::deque<Node<T>> nodes;
      std::deque<Node<T>> temp;
      nodes.push_back(generations[row][col]);
      for (size_t i = row+1; i < generations.size(); ++i)
      { 
        for (auto j = generations[i].begin(); j != generations[i].end(); ++j)
        {
          for (auto k = nodes.begin(); k != nodes.end(); ++k)
          {
            if ((*j).parent == (*k).data)
            {
              ids.push_back((*j).data);
              temp.push_back(*j);
            }
          }
        } 
        nodes.resize(0);
        nodes = temp;
        temp.resize(0);
        if (nodes.empty())
        {
          return ids;
          break;
        }
      }
    } 

    return ids;
  }

  size_t count_alive_progeny(T cell_id)
  {
    size_t alives = 0; 
    if (cell_id >= ((generations.back()).front()).data)
    {
      alives = 1;
      return alives;
    }
    else
    {
      std::vector<T> descendants = descendant_ids(cell_id);
      for (const auto& data : descendants)
      {
        if (data >= ((generations.back()).front()).data)
        {
          alives += 1;
        }
      }
      return alives;
    }
  }

  void count_all_alive()
  { 

    if (((generations.back()).front()).data != 0)
    {
      std::vector<T> nodes;
      //handles last layer separately
      for (auto it = (generations.back()).begin(); it != (generations.back()).end(); ++it)
      {
        nodes.push_back((*it).data);
        all_alive += 1;
        if (std::find(nodes.begin(), nodes.end(), (*it).parent) == nodes.end())
        {
          nodes.push_back((*it).parent);
        }
      }

      for (auto it = generations.rbegin()+1; it != generations.rend()-1; ++it)
      {
        for (auto iv = (*it).begin(); iv != (*it).end(); ++iv)
        {
          if (std::find(nodes.begin(), nodes.end(), (*iv).data) != nodes.end())
          { 
            nodes.push_back((*iv).parent);
            all_alive += 1;
          }
        }
      }
    }
    else
    {
      all_alive += 0;
    } 
    //std::cout << "All alive cells: " << all_alive << std::endl; 
  }
  // WARNING
  // This function generates the number of mutations happening in a single generation 
  // The number of mutations is a single "int" value, being the extracted number of mutations expected. 
  size_t generate_mutations(std::mt19937 &engine)
  {
    size_t attempts = max_node_index * n_bases; 
    std::binomial_distribution<long long unsigned> bd(attempts, mutation_rate);
    size_t mutations = bd(engine);
    std::cout << mutations << std::endl;
    return mutations;
  }

  size_t random_mutant_cell() // mutates a single base in a single cell --> attempts = max_node_index * n_bases
  { 
    double attempts = static_cast<double>(max_node_index*n_bases); 
    double mutant_idx = 1.0 + rnd.rannyu(0.0, attempts);
    double app = 0.0;
    double app2 = std::modf(mutant_idx, &app); 
    if (app2 > 0.5)
    {
      app += 1.0;
    }
    return (static_cast<size_t>(app) % max_node_index);
  }

  void set_mutations_with_count(std::mt19937 &engine)
  {
    n_mutations = generate_mutations(engine); // this is a single "int"

    //std::cout << "Number of mutations = " << n_mutations << std::endl;
    for (size_t i = 0; i < n_mutations; ++i)
    {
      size_t cell_id = random_mutant_cell();
      //std::cout << "Mutated cell id = " << cell_id << std::endl;
      size_t alive_mutated_progeny = count_alive_progeny(cell_id);
      //std::cout << "Alive progeny mutated cell = " << alive_mutated_progeny << std::endl;
      //std::cout << "Survived nodes = " << survived_nodes << std::endl;
      if (alive_mutated_progeny != 0)
      {
        double frequency = static_cast<double>(alive_mutated_progeny) / static_cast<double>(survived_nodes);
        freq_mutations.push_back(frequency);
        count_mutations += 1;
        //std::cout << "True frequency = " << frequency << std::endl;
      }
    }
  }

  void set_mutations(std::mt19937 &engine)
  {
    n_mutations = generate_mutations(engine); // mutations happening in a single generation (time-step)

    std::cout << "Number of mutations = " << n_mutations << std::endl;
    for (size_t i = 0; i < n_mutations; ++i)
    {
      size_t cell_id = random_mutant_cell(); // it's random_mutant_cell to *actually* perform the mutation, by choosing a random cell&base
      //std::cout << "Mutated cell id = " << cell_id << std::endl;
      // chosen a cell we say it's mutated, all of its progeny will be equally mutated.
      size_t alive_mutated_progeny = count_alive_progeny(cell_id);
      //std::cout << "Alive progeny mutated cell = " << alive_mutated_progeny << std::endl;
      //std::cout << "Survived nodes = " << survived_nodes << std::endl;
      if (alive_mutated_progeny != 0)
      {
        double frequency = static_cast<double>(alive_mutated_progeny) / static_cast<double>(survived_nodes);
        freq_mutations.push_back(frequency);
        //count_mutations += 1;
        //std::cout << "True frequency = " << frequency << std::endl;
      }
      // Test print: --implem FG 21/02/25
      //if (i%10 == 0) {
      //  std::cout << "i =" << i << std::endl;
      //}
    }
  }

  void set_n_mut_prior(std::mt19937 &engine)
  {
    size_t attempts = n_zero * n_bases; 
    std::binomial_distribution<long long unsigned> bd(attempts, mutation_rate);
    n_mut_prior = bd(engine);
    for (size_t i = 0; i < n_mut_prior; ++i)
    {
      freq_mutations.push_back(1.0);
    }
  }

  size_t ploidy_poisson(std::mt19937 &engine)
  {
    size_t lambda = ploidy - 1;
    std::poisson_distribution<long unsigned> poisson(lambda);
    return (1 + poisson(engine));
  }

  size_t sampling_coverage(std::mt19937 &engine, double frequency, bool set_ploidy)
  {
    size_t coverage_p = coverage;
    if (set_ploidy == true)
    {
      coverage_p = coverage / ploidy;
    } 
    std::binomial_distribution<long unsigned> bd(coverage_p, frequency);
    size_t cov = bd(engine); 
    return cov;
  }

  int sampling_C(std::mt19937 &engine, double frequency, size_t coverage)
  { 
    std::binomial_distribution<long unsigned> bd(coverage, frequency);
    size_t cov = bd(engine); 
    return cov;
  }

  size_t sampling_coverage_ploidy(std::mt19937 &engine, double frequency, double pp,  bool set_ploidy)
  {
    double freq_p = frequency/pp;
    std::binomial_distribution<long unsigned> bd(coverage, freq_p);
    size_t cov = bd(engine); 
    return cov;
  }

  size_t set_false_negative(std::mt19937 &engine, size_t r_1)
  {
    std::binomial_distribution<long unsigned> bd(r_1, seq_error);
    size_t n_false_negative = bd(engine);  
    return n_false_negative;
  }

  size_t set_false_positive(std::mt19937 &engine, size_t r_1, bool set_ploidy, size_t pp)
  {
    std::binomial_distribution<long unsigned> bd(coverage - r_1, seq_error);
    size_t n_false_positive = bd(engine);  
    return n_false_positive;
  }

  // Here starts the sequencing part

  size_t error_no_mut_bases(std::mt19937 &engine, const double &e)
  {
    double r = static_cast<double>(coverage)/e;
    double err_bias = 1.0 - binomial_cdf(r, static_cast<int>(coverage), seq_error); 
    std::binomial_distribution<long unsigned> bd(n_bases - freq_mutations.size(), err_bias);
    size_t c = bd(engine);
    return c;
  }

  void sequencing_err_not_mutated_bases(const gsl_rng *r, double pp, 
                                        bool set_ploidy, size_t n_zeros,
                                        std::vector<size_t> &mut_bases_error)
  {    

    double vec_pr[coverage];
    for (size_t i = 0; i < coverage; ++i)
    { 
      double pr = gsl_ran_binomial_pdf(i, seq_error, coverage);
      vec_pr[i] = pr;
    }

    unsigned int vec_values[coverage];
    gsl_ran_multinomial(r, coverage, n_zeros, vec_pr, vec_values);
  
    for (size_t i = 0; i < coverage; ++i)
    {
      mut_bases_error[i] += static_cast<size_t>(vec_values[i]);
  //  }
  //}

  //void parameters_control(std::mt19937 &engine)
  //{
  //  set_mutations(engine);
  //  for (const auto& f : freq_mutations)
  //  {
  //    size_t pp = ploidy_poisson(engine);
  //    size_t r_1 = sampling_coverage_ploidy(engine, f, pp, true);
  //    while (r_1 <= 0)
  //    {
  //      r_1 = sampling_coverage_ploidy(engine, f, pp, true);
  //    }

  //    double f_tot = static_cast<double>(r_1) * 
  //                   static_cast<double>(pp) *
  //                   (1.0/static_cast<double>(coverage));
  //   
  //    
  //    threshold = (1.0 / static_cast<double>(survived_nodes)) + 
  //                (K * (1.0 / std::sqrt(static_cast<double>(survived_nodes))) * 
  //                (1.0/std::sqrt(static_cast<double>(coverage)/static_cast<double>(pp))));
  //    //if (f_tot > 0.0)
  //    if (f_tot >= threshold && f_tot < 1.0)
  //    {
  //      count_mutations += 1;
  //    }
  //  }
  //}

  //void compute_mutations(std::mt19937 &engine)
  //{
  //  set_mutations(engine);
  //  for (const auto& f : freq_mutations)
  //  {
  //      size_t r_1 = sampling_coverage(engine, f, true);

  //      size_t r_2 = 0.0;
  //      if (seq_error > 0.0)
  //      { 
  //        r_2 = set_false_negative(engine, r_1);
  //      }

  //      size_t r_3 = 0.0;
  //      if (seq_error > 0.0)
  //      { 
  //        r_3 = set_false_positive(engine, r_1, true);
  //      }
  //      //size_t r_3 = 0;
  //      //std::cout << "R3 = " << r_3 << std::endl;
  //      
  //      double f_tot = static_cast<double>(r_1 - r_2 + r_3) * 
  //                     static_cast<double>(ploidy) *
  //                     (1.0/static_cast<double>(coverage));

  //      threshold = (1.0 / static_cast<double>(survived_nodes)) + 
  //                  (K * (1.0 / std::sqrt(static_cast<double>(survived_nodes))) * 
  //                  (1.0/std::sqrt(static_cast<double>(coverage)/static_cast<double>(ploidy))));

  //      //if (f_tot > 0.0)
  //      if (f_tot >= threshold && f_tot < 1.0)
  //      {
  //        count_mutations += 1;
  //      }
  //  }

  //  size_t count_err = 0;
  //  if (seq_error > 0.0)
  //  {
  //    count_err = error_no_mut_bases(engine, threshold);
  //  }
  //  count_mutations += count_err; 
  //}


  //void compute_mutation_new(std::mt19937 &engine, const size_t &e)
  //{
  //  count_mutations = 0;
  //  //double thr = (1.0/static_cast<double>(e)) + (K/std::sqrt(static_cast<double>(coverage)/static_cast<double>(ploidy)))*std::sqrt(1.0/static_cast<double>(e));
  //  double thr = ((1.0/static_cast<double>(e))*static_cast<double>(ploidy)) + (K/std::sqrt(static_cast<double>(coverage)/static_cast<double>(ploidy)))*std::sqrt(1.0/static_cast<double>(e));
  //  for (const auto& f : freq_mutations) 
  //  { 
  //    //size_t r_1 = sampling_coverage_ploidy(engine, f, static_cast<double>(ploidy), true);
  //    size_t r_1 = sampling_coverage_ploidy(engine, f, static_cast<double>(ploidy), false);

  //    size_t r_2 = 0.0;
  //    if (seq_error > 0.0)
  //    { 
  //      r_2 = set_false_negative(engine, r_1);
  //    }

  //    size_t r_3 = 0.0;
  //    if (seq_error > 0.0)
  //    { 
  //      r_3 = set_false_positive(engine, r_1, false, ploidy);
  //    }

  //    double f_hat = static_cast<double>(r_1-r_2+r_3) * 
  //                   static_cast<double>(ploidy) *
  //                   (1.0/static_cast<double>(coverage));
  //     
  //    //if (f_hat >= thr)
  //    //if ((f_hat >= thr) && (f_hat < 0.99))
  //    if (f_hat > 0.0)
  //    {
  //      count_mutations += 1;
  //    }
  //  }

  //  size_t count_err = 0;
  //  if (seq_error > 0.0)
  //  {
  //    count_err = error_no_mut_bases(engine, static_cast<double>(e));
  //  }
  //  count_mutations += count_err; 
  //}

  //void compute_frequencies(std::mt19937 &engine, const gsl_rng *r, 
  //                         std::vector<double> &frequencies_sampling, 
  //                         std::vector<size_t> &zeros, 
  //                         std::vector<size_t> &mut_bases_error)
  //{ 
  //  //std::cout << "////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
  //  //std::cout << "////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
  //  //std::cout << "////////////////////////////////////////////////////////////////////////////////////////" << std::endl;

  //  std::vector<double> f_hat_zeros;
  //  for (const auto& f : freq_mutations) 
  //  {

  //    //size_t pp = ploidy_poisson(engine);
  //    size_t pp = ploidy;
  //    //std::cout << "True frequency = " << f << std::endl;
  //    size_t r_1 = sampling_coverage_ploidy(engine, f, static_cast<double>(pp), false);

  //    //std::cout << "Coverage sampling = " << r_1 << std::endl;
  //    size_t r_2 = 0.0;
  //    if (seq_error > 0.0)
  //    { 
  //      r_2 = set_false_negative(engine, r_1);
  //    }

  //    //std::cout << "False negative = " << r_2 << std::endl;
  //    size_t r_3 = 0.0;
  //    if (seq_error > 0.0)
  //    { 
  //      r_3 = set_false_positive(engine, r_1, false, pp);
  //    }

  //    //std::cout << "False positive = " << r_3 << std::endl;
  //    double f_hat = static_cast<double>(r_1-r_2+r_3) * 
  //                   static_cast<double>(pp) *
  //                   (1.0/static_cast<double>(coverage));

  //    //std::cout << "F_HAT = " << f_hat << std::endl;
  //    //if ((f_hat > 0.0) && (f_hat < 0.99))
  //    if (f_hat > 0.0)
  //    {
  //      //frequencies_sampling.push_back(f_hat);
  //      frequencies_sampling.push_back(static_cast<double>(r_1-r_2+r_3));
  //    }

  //    if (std::abs(f_hat) <= 0.0)
  //    {
  //      f_hat_zeros.push_back(f_hat);
  //    }
  //  }
  //  
  //  size_t z = get_not_mutated_bases(f_hat_zeros.size());
  //  zeros.push_back(z);
  //  sequencing_err_not_mutated_bases(r, static_cast<double>(ploidy), false, z, mut_bases_error);
  //}

  //void compute_mutation_new_k(std::mt19937 &engine, const size_t &e, const double &k_)
  //{
  //  count_mutations = 0;
  //  double thr = (1.0/static_cast<double>(e)) + (K/std::sqrt(static_cast<double>(coverage)/static_cast<double>(ploidy)))*std::sqrt(1.0/static_cast<double>(e));
  //  for (const auto& f : freq_mutations) 
  //  { 
  //    size_t r_1 = sampling_C(engine, f, coverage); 

  //    size_t r_2 = 0.0;
  //    if (seq_error > 0.0)
  //    { 
  //      r_2 = set_false_negative(engine, r_1);
  //    }

  //    size_t r_3 = 0.0;
  //    if (seq_error > 0.0)
  //    { 
  //      r_3 = set_false_positive(engine, r_1, false);
  //    }

  //    double f_hat = static_cast<double>(r_1-r_2+r_3) * 
  //                   (1.0/static_cast<double>(coverage));
  //     
  //    //if (f_hat >= thr)
  //    //if ((f_hat >= thr) && (f_hat <= 0.95))
  //    if ((f_hat >= thr) && (f_hat < 1.0))
  //    {
  //      count_mutations += 1;
  //    }
  //  }

  //  size_t count_err = 0;
  //  if (seq_error > 0.0)
  //  {
  //    count_err = error_no_mut_bases(engine, static_cast<double>(e));
  //  }
  //  count_mutations += count_err; 
  //}

  //void compute_mutation_wilson(std::mt19937 &engine, const size_t &e)
  //{
  //  count_mutations = 0;
  //  for (const auto& f : freq_mutations) 
  //  { 
  //    size_t r_1 = sampling_C(engine, f, coverage); 

  //    double f_hat = static_cast<double>(r_1) * 
  //                   static_cast<double>(ploidy) *
  //                   (1.0/static_cast<double>(coverage));

  //    f_hat = (1.0/(1.0 + ((1.96*1.96)/static_cast<double>(coverage))))*(f_hat + ((1.96*1.96)/(2*static_cast<double>(coverage))));
  //    
  //    double thr = (1.0/static_cast<double>(e)) + (K/std::sqrt(static_cast<double>(coverage)))*std::sqrt(1.0/static_cast<double>(e));

  //    if (f_hat >= thr)
  //    //if ((f_hat >= thr) && (f_hat <= 0.95))
  //    //if ((f_hat >= thr) && (f_hat < 1.0))
  //    {
  //      count_mutations += 1;
  //    }
  //  }
  //}

  //void compute_mutation_prior(std::mt19937 &engine, const size_t &e)
  //{
  //  count_mutations = 0;
  //  double thr = (1.0/static_cast<double>(e));
  //  for (const auto& f : freq_mutations) 
  //  { 
  //    size_t r_1 = sampling_C(engine, f, coverage); 

  //    size_t r_2 = 0.0;
  //    if (seq_error > 0.0)
  //    { 
  //      r_2 = set_false_negative(engine, r_1);
  //    }

  //    size_t r_3 = 0.0;
  //    if (seq_error > 0.0)
  //    { 
  //      r_3 = set_false_positive(engine, r_1, true);
  //    }

  //    //double f_hat = (static_cast<double>(r_1)-1.0)*(1.0/static_cast<double>(coverage));
  //    double f_hat = (static_cast<double>(r_1-r_2+r_3)-1.0)*(1.0/static_cast<double>(coverage)); 

  //    if (f_hat >= thr)
  //    //if ((f_hat >= thr) && (f_hat <= 0.95))
  //    //if ((f_hat >= thr) && (f_hat < 1.0))
  //    {
  //      count_mutations += 1;
  //    }
  //  }

  //  size_t count_err = 0;
  //  if (seq_error > 0.0)
  //  {
  //    count_err = error_no_mut_bases(engine, thr);
  //  }
  //  count_mutations += count_err; 
  //}

  //long double mr_estimator()
  //{ 
  //  count_all_alive();
  //  long double mr_est = -std::log(1.0 - (static_cast<long double>(count_mutations)/ 
  //                                 static_cast<long double>(n_bases)))/static_cast<long double>(all_alive);
  //  return mr_est;
  //}

  //long double mr_estimator_correct(const size_t &e) 
  //{
  //  //generations
  //  double G = std::log(static_cast<double>(e))/std::log(2.0*(1.0-death_probability));
  //  //double G = static_cast<double>(generations.size()-1);
  //  //alive cells
  //  double A = (1.0 - death_probability*death_probability)*((std::pow((2.0*(1.0-death_probability)), G - 1.0) - 1.0) /std::log(2.0*(1.0 - death_probability))) + static_cast<double>(e);
  //  long double mr_est = -std::log(1.0 - (static_cast<long double>(count_mutations)/ 
  //                                 static_cast<long double>(n_bases)))/static_cast<long double>(A);
  //  return mr_est;
  //}

  //long double mr_estimator_survived() 
  //{
  //  //generations
  //  double G = static_cast<double>(generations.size()-1);
  //  //alive cells
  //  double A = (1.0 - death_probability*death_probability)*((std::pow((2.0*(1.0-death_probability)), G - 1.0) - 1.0) /std::log(2.0*(1.0 - death_probability))) + static_cast<double>(survived_nodes);
  //  long double mr_est = -std::log(1.0 - (static_cast<long double>(count_mutations)/ 
  //                                 static_cast<long double>(n_bases)))/static_cast<long double>(A);
  //  return mr_est;
  //}

  //long double mr_estimator_new()
  //{ 
  //  long double mr_est = -std::log(1.0 - (static_cast<long double>(count_mutations)/ 
  //                                 static_cast<long double>(n_bases)))/static_cast<long double>(all_alive);
  //  return mr_est;
  //}

  //long double mr_estimator_base()
  //{
  //  long double mr_est = -std::log(1.0 - (static_cast<long double>(count_mutations)/ 
  //                                 static_cast<long double>(n_bases)))/static_cast<long double>(max_node_index);
  //  //std::cout << "mr = " << mr_est << "count mutations = " << count_mutations << std::endl;;
  //  return mr_est;
  //}

  //void mr_estimator_no_ret()
  //{
  //  count_all_alive();
  //  //std::cout << "All alive = " << all_alive << std::endl;
  //  long double mr_est = -std::log(1.0 - (static_cast<long double>(count_mutations)/ 
  //                                 static_cast<long double>(n_bases)))/static_cast<long double>(all_alive);
  //  //std::cout << "mr = " << mr_est << "count mutations = " << count_mutations << std::endl;;
  //}

  void print_tree(std::ostream &out = std::cout) const
  {
    out << std::endl;
    out << "Sequential binary tree" << std::endl;
    out << std::endl;
    for (const auto& layer : generations)
    {
      for (const auto& node : layer)
      {
        out << "(" << node.data << " , " << node.parent << ")" << "\t";
      }
      out << std::endl;
      out << std::endl;
    } 
    out << std::endl;
    out << "Nodes survived: " << survived_nodes << std::endl; 
    out << std::endl;
    out << "Max node index: " << max_node_index << std::endl;
    out << std::endl;
  }

  friend std::ostream& operator<<(std::ostream &os, const Tree &tree)
  {
    tree.print_tree(os);
    return os;
  }

  void print_vector_tree()
  {
    for (const auto& el : tree)
    {
      std::cout << el.data << std::endl;
    }
  }

  void print_node(const Node<T> &node)
  {
    std::cout << node.data << "\t" << node.parent << std::endl; 
  }

  void print_last_layer()
  {
    for (const auto& node : generations.back())
    {
      std::cout << "(" << node.data << " , " << node.parent << ")" << "\t";
    }
    std::cout << std::endl;
  }

  void print_freq_mutations()
  {
    for (const auto& f : freq_mutations)
    {
      std::cout << f << "\t";
    }
    std::cout << std::endl;
  }

  void print_all_alive()
  {
    std::cout << "All alive cells: " << std::endl;
    std::cout << all_alive << std::endl;
  }

  double binomial_cdf(double m, int cov, double p)
  { 
    double n = static_cast<double>(cov) / static_cast<double>(ploidy);
    double cdf = 0.0;
    double b = 0.0;
  
    for (double k = 0.0; k < m + 1; k+=1.0)
    {
      if (k > 0.0)
      {
        b += std::log(n - k + 1.0) - std::log(k);
      }
      double log_pmf_k = b + (k * std::log(p)) + ((n - k) * std::log(1.0 - p));
      cdf += std::exp(log_pmf_k);
    }
    
    return cdf;
  }

private:
  std::vector<std::vector<Node<T>>> generations;
  std::vector<Node<T>> tree;
  long double mutation_rate;
  char death_prob_type;
  double death_prob_p0;
  double death_prob_p1;
  double seq_error;
  double K;
  double threshold;
  double cov_over_p;
  size_t max_node_index;
  size_t min_survived_nodes_number;
  size_t max_generations;
  size_t survived_nodes;
  size_t n_mutations;
  size_t count_mutations;
  size_t count_thr_mut;
  size_t all_alive;
  size_t n_mut_prior;
  size_t n_bases;
  size_t coverage;
  size_t ploidy;
  size_t n_zero;
  size_t ploidy_p;
};

#endif //TREE_LD_H
