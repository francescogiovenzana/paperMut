#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <format>
//#include <libfmt/format.h>
#include <chrono>
#include <vector>
#include <random>
#include "tree.h"
#include "../random/random.h"
#include "XoshiroCpp.hpp"
#include <filesystem>

using namespace std::chrono;

// here set paths
// Define SEED and S_DIR if not passed as macros
// Not needed here bcause already defined in tree.h
//const std::string SEED = "../seed"; 
//const std::string S_DIR = "../output";

// Implementing selection.
// It is now necessary to take as input the preferred choice for death probability distribution. 
// It should be suitable having in input_file on the same row the type of distribution and the parameters.
// Format in input_file
// [C]\t[Value]\t[NULL] --> for single value for all nodes
// [U]\t[Lower_bound]\t[Upper_bound] --> for uniform distribution

// This function verifies wheter the death probability lays within appropriate boundaries. 
// In fact, the death probability should be <= 0.5, in order to not cause tree extinction with prob ==1
int check_death_prob(char death_prob_type, double death_prob_p0, double death_prob_p1)
  {
  int test = 0;
  if (death_prob_type == 'C')
  {
    if (death_prob_p0 > 0.5)
    {
      std::cout << std::endl;
      std::cout << "Value of death probability not allowed. " << std::endl;
      std::cout << "Tree extinction with probability = 1.   " << std::endl;
      std::cout << "Death probability needs to be <= 0.5.   " << std::endl;
      std::cout << std::endl;
      test=1;
    }
  }
  else if (death_prob_type == 'U')
  {
    if (death_prob_p1 > 0.5)
    {
      std::cout << std::endl;
      std::cout << "Value of death probability not allowed. " << std::endl;
      std::cout << "Tree extinction with probability = 1.   " << std::endl;
      std::cout << "Death probability needs to be <= 0.5.   " << std::endl;
      std::cout << std::endl;
      test=1;
    }
  }
  else
  {
    std::cerr << "Error: Death probability type not recognized. Exit." << std::endl;
    test=1;
  }

  return test;
} // end function

int main(int argc, char* argv[])
{
 
  if (argc != 6)
  {
    std::cout  << "Use: " << argv[0] << " input.dat [strain number] [sub-strain number] [depth] [ploidy]" << std::endl;
    return 1;
  }

  std::ifstream data_file(std::string{argv[1]});
  if (!data_file.is_open())
  {
    std::cerr << "File input.dat not opened. Exit" << std::endl;
    return 1;
  }

  // Initialize random number generator
  Random rnd(SEED + "/Primes", SEED + "/seed.in");
  for (size_t i = 0; i < 100000; ++i)
  {
    rnd.rannyu();
  } 

  // Read input parameters
  std::string seq_name;
  size_t n_blocks = 0;
  size_t n_steps = 0; 
  size_t n_bases = 0; //10^8 or 10^9
  size_t coverage = 0;
  size_t ploidy = 1;
  size_t n_zero = 0;

  // Death_probability parameters
  char death_prob_type; // 'C'onst or 'U'nif
  double death_prob_p0 = 0.0; // Const value or lower bound
  double death_prob_p1 = 0.0; // NULL or upper bound
  double seq_error = 0.0;
  double K = 0.0;
  size_t min_survived_nodes_number = 0;
  size_t max_generations = 0;
  long double MUT_rate = 0.0;
  data_file >> seq_name >> n_blocks >> n_steps >> n_bases >> 
               coverage >> ploidy >> n_zero >> death_prob_type >> death_prob_p0;
  if (death_prob_type == 'U') 
  {
    data_file >> death_prob_p1;
   }
  data_file >> seq_error >> K >> min_survived_nodes_number >> max_generations >> MUT_rate;
  data_file.close();

  std::cout << "Simulation parameters:\n";
  std::cout << seq_name << std::endl;
  std::cout << "Number of simulation blocks: " << n_blocks << std::endl;
  std::cout << "Number of simulation steps: " << n_steps << std::endl;
  std::cout << "Number of bases: " << n_bases << std::endl;
  std::cout << "Coverage: " << coverage << std::endl;
  std::cout << "Ploidy: " << ploidy << std::endl;
  std::cout << "N0: " << n_zero << std::endl;
  std::cout << "Sequencing error: " << seq_error << std::endl;
  std::cout << "K: " << K << std::endl;

  if (death_prob_type == 'C')
  {
    std::cout << "Death probability: " << death_prob_p0 << std::endl;
  }
  else if (death_prob_type == 'U')
  {
    std::cout << "Lower bound: " << death_prob_p0 << std::endl;
    std::cout << "Upper bound: " << death_prob_p1 << std::endl;
  }
  else
  {
    std::cerr << "Error: Death probability type not recognized. Exit." << std::endl;
    return 1;
  }
  // end of part added for selection implementation

  std::cout << "Death probability distribution: " << death_prob_type << std::endl;
  std::cout << "Min number of survived cells: " << min_survived_nodes_number << std::endl;
  std::cout << "Max number of generations: " << max_generations << std::endl;
  std::cout << "Mutation rate: " << MUT_rate << std::endl;

  

  // This section verifies wheter the death probability lays within appropriate boundaries. 
  // In fact, the death probability should be <= 0.5, in order to not cause tree extinction with prob ==1
  int test = check_death_prob(death_prob_type, death_prob_p0, death_prob_p1);
  if (test == 1) {
    return 1;
  }

  auto start_timer = high_resolution_clock::now();

  // Initialize Xoshiro RNG
  constexpr std::uint64_t seed_mutations = 12345;
  XoshiroCpp::Xoshiro256PlusPlus rng(seed_mutations);
  auto seed_r = rng(); 
  static std::mt19937 engine(seed_r);

 
  std::vector<Node<size_t>> sequential_tree; 
  Tree<size_t> tree(sequential_tree, death_prob_type, death_prob_p0, death_prob_p1, seq_error, K, 
                    max_generations, n_bases, coverage, ploidy, n_zero);
  
  // Backup version with death_probability at tree level
  /*Tree<size_t> tree(sequential_tree, death_probability, seq_error, K, 
                    max_generations, n_bases, coverage, ploidy, n_zero);
                    */
  tree.warm_up_random();

  //std::cout << "type=" << tree.get_death_prob_type() << std::endl;

  /*
  // Debugging
  std::string curr_path = std::filesystem::current_path();
  std::cout << "Current path =" << curr_path << std::endl;
  */
  
  // Simulation loop

  for (size_t j = 0; j < n_steps; ++j)
  {
    std::cout << "Step = " << j << std::endl;
    // See chat here for some possible improvements
    
    std::string output_dir = S_DIR + "/" + argv[2] + "/" + argv[3] + "/" +
                                 argv[4] + "x/ploidy" + argv[5] + "/n_ext/";

    // Making sure the directory exists. If not, creating it
    if (!std::filesystem::exists(output_dir)) {
      if (std::filesystem::create_directories(output_dir))
      {
        std::cout << "Directories created: " << output_dir << std::endl;
      }
      else
      {
        std::cerr << "Error: Could not create directories: " << output_dir << std::endl;
        return 1;
      }
    }
    std::cout << "Output directory: " << output_dir << std::endl;

    std::string file_freqs = output_dir + "file_frequencies_" + std::to_string(j) + ".dat";
    std::string file_nodes = output_dir + "file_n_nodes_" + std::to_string(j) + ".dat";

    std::ofstream file_out_freqs(file_freqs);
    std::ofstream file_out_nodes(file_nodes);

    if (!file_out_freqs.is_open() || !file_out_nodes.is_open())
      {
        std::cerr << "Error: Could not open output files at step " << j << std::endl;
        return 1;
    }
        //This all reconnects before "start_tree"
    /*

    std::string file_freqs = "file_frequencies_"+std::to_string(j)+".dat";
    std::string file_nodes = "file_n_nodes_"+std::to_string(j)+".dat";

    //std::cout << "Test 1" << std::endl;
    std::ofstream file_out_freqs(S_DIR +std::string{argv[2]}+std::string{"/"}+std::string{argv[3]}+std::string{"/"}
                                           +std::string{argv[4]}+"x/ploidy"+std::string{argv[5]}
                                           +"/n_ext/"+file_freqs, std::ios_base::binary);
    std::ofstream file_out_nodes(S_DIR +std::string{argv[2]}+std::string{"/"}+std::string{argv[3]}+std::string{"/"}
                                           +std::string{argv[4]}+"x/ploidy"+std::string{argv[5]}
                                           +"/n_ext/"+file_nodes, std::ios_base::binary);

    */
    tree.start_tree(MUT_rate, min_survived_nodes_number, max_generations, engine);
    //std::cout << "Tree started" << std::endl;
    tree.set_mutations(engine);

    //std::cout << "Test 1" << std::endl;

    for (const auto& el : tree.freq_mutations)
    {
      //std::cout << "el=" << el << "\t";
      file_out_freqs << el << std::endl;
    }
    
    file_out_nodes << tree.get_survived_nodes() << std::endl;

    file_out_nodes.close();
    file_out_freqs.close();
  }
  
  auto stop_timer = high_resolution_clock::now();
  auto measure_time = duration_cast<minutes>(stop_timer - start_timer);
  std::cout << "Execution time:  " << measure_time.count() << "  min." << std::endl;
  return 0;
}
