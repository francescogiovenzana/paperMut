#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cmath>
#include <vector>
#include <set>
#include <algorithm>
#include <cassert>
#include <utility>

using namespace std::chrono;

void distribution(std::vector<int> &coverages, std::vector<std::pair<int,int>> &counts);

int main (int argc, char* argv[])
{

  if (argc != 3)
  {
    std::cout << "Use: " << argv[0] << " [path_pseudo-pileup.dat] [output_path]" << std::endl; 
    return 1;
  }

  std::string seq_file = std::string{argv[1]};
  std::string output_path = std::string{argv[2]};
  std::string distr_forward_out_file = "coverage_forward_distribution.dat";
  std::string distr_reverse_out_file = "coverage_reverse_distribution.dat";
  std::string mean_coverage_forward_file = "mean_coverage_forward.dat";
  std::string mean_coverage_reverse_file = "mean_coverage_reverse.dat";
  //std::string cov_forward_sorted_file = "coverage_forward_sorted.dat";
  //std::string cov_reverse_sorted_file = "coverage_reverse_sorted.dat";
  //input file
  std::ifstream in_file(seq_file, std::ios_base::binary);
  if (!in_file.is_open())
  {
    std::cerr << "File infile.dat not opened. Exit" << std::endl;
    return 1;
  }

  //output file
  std::ofstream distr_forward_out(output_path+distr_forward_out_file, std::ios_base::binary);
  if (!distr_forward_out.is_open())
  {
    std::cerr << "File distr_out.dat not opened. Exit" << std::endl;
    return 1;
  }

  //output file
  std::ofstream distr_reverse_out(output_path+distr_reverse_out_file, std::ios_base::binary);
  if (!distr_reverse_out.is_open())
  {
    std::cerr << "File distr_out.dat not opened. Exit" << std::endl;
    return 1;
  }

  //output file
  std::ofstream mean_coverage_forward_out(output_path+mean_coverage_forward_file, std::ios_base::binary);
  if (!mean_coverage_forward_out.is_open())
  {
    std::cerr << "File mean_cov_out.dat not opened. Exit" << std::endl;
    return 1;
  }

  //output file
  std::ofstream mean_coverage_reverse_out(output_path+mean_coverage_reverse_file, std::ios_base::binary);
  if (!mean_coverage_reverse_out.is_open())
  {
    std::cerr << "File mean_cov_out.dat not opened. Exit" << std::endl;
    return 1;
  }

  //output file
  //std::ofstream cov_forward_sorted_out(output_path+cov_forward_sorted_file, std::ios_base::binary);
  //if (!cov_forward_sorted_out.is_open())
  //{
  //  std::cerr << "File cov_sorted_out.dat not opened. Exit" << std::endl;
  //  return 1;
  //}

  //output file
  //std::ofstream cov_reverse_sorted_out(output_path+cov_reverse_sorted_file, std::ios_base::binary);
  //if (!cov_reverse_sorted_out.is_open())
  //{
  //  std::cerr << "File cov_sorted_out.dat not opened. Exit" << std::endl;
  //  return 1;
  //}

  auto start_timer = high_resolution_clock::now();

  std::vector<int> aux_empty;

  std::vector<int> coverages_forward;
  std::vector<int> coverages_reverse;
  std::vector<std::pair<int,int>> counts_forward;
  std::vector<std::pair<int,int>> counts_reverse;
  double mean_coverage = 0.0;
  double var_coverage = 0.0;
  double std_coverage = 0.0;
  //reference base
  std::string A_s = "A";
  std::string G_s = "G";
  std::string C_s = "C";
  std::string T_s = "T";

  //ancestor
  std::string line;
  std::string chr;
  std::string base_num;
  std::string ref;
  std::string A_count;
  std::string G_count;
  std::string C_count;
  std::string T_count;
  std::string a_count;
  std::string g_count;
  std::string c_count;
  std::string t_count;
  int sequencing_forward = 0;
  int sequencing_reverse = 0;
  double total_coverage = 0.0;
  double lines = 0.0;

  while (std::getline(in_file, line))
  {
    std::stringstream data_seq(line);
    lines += 1.0;
    data_seq >> chr >> base_num >> ref >> A_count >> G_count >> C_count >> T_count >> 
                                          a_count >> g_count >> c_count >> t_count; 

    sequencing_forward = std::stoi(A_count) + std::stoi(G_count) + std::stoi(C_count) + std::stoi(T_count);

    coverages_forward.push_back(sequencing_forward);
    total_coverage += static_cast<double>(sequencing_forward);
  }

  std::cout << "Input file reading completed!" << std::endl;

  ////////////////FORWARD/////////////////////////////////////////////////////////////////////
  distribution(coverages_forward, counts_forward);
  for (const auto &el : counts_forward)
  {
    distr_forward_out << el.first << "\t" << el.second << std::endl;
  }
  distr_forward_out.close();
  in_file.close();
 
  mean_coverage = total_coverage / lines;
  for (auto& el : coverages_forward)
  {
    var_coverage += (static_cast<double>(el) - mean_coverage)*(static_cast<double>(el) - mean_coverage);
  }
  var_coverage = var_coverage / lines;
  std_coverage = std::sqrt(var_coverage);
  mean_coverage_forward_out << mean_coverage << "\t" << std_coverage << std::endl;
  mean_coverage_forward_out.close();
  ////////////////END FORWARD/////////////////////////////////////////////////////////////////

  ///////////////////////////////
  mean_coverage = 0.0;
  var_coverage = 0.0;
  std_coverage = 0.0;
  lines = 0.0;
  total_coverage = 0.0;
  coverages_forward.clear();
  coverages_forward.shrink_to_fit();
  ///////////////////////////////

  std::cout << "REVERSE" << std::endl;
  std::cout << std::endl;

  in_file.open(seq_file, std::ios_base::binary);
  if (!in_file.is_open())
  {
    std::cerr << "File infile.dat not opened. Exit" << std::endl;
    return 1;
  }

  while (std::getline(in_file, line))
  {
    std::stringstream data_seq(line);
    lines += 1.0;
    data_seq >> chr >> base_num >> ref >> A_count >> G_count >> C_count >> T_count >> 
                                          a_count >> g_count >> c_count >> t_count; 

    sequencing_reverse = std::stoi(a_count) + std::stoi(g_count) + std::stoi(c_count) + std::stoi(t_count);

    coverages_reverse.push_back(sequencing_reverse);
    total_coverage += static_cast<double>(sequencing_reverse);
  }

  std::cout << "Input file reading completed!" << std::endl;

  ////////////////REVERSE/////////////////////////////////////////////////////////////////////
  distribution(coverages_reverse, counts_reverse);
  for (const auto &el : counts_reverse)
  {
    distr_reverse_out << el.first << "\t" << el.second << std::endl;
  }
  distr_reverse_out.close();
  in_file.close();
 
  mean_coverage = total_coverage / lines;
  for (auto& el : coverages_reverse)
  {
    var_coverage += (static_cast<double>(el) - mean_coverage)*(static_cast<double>(el) - mean_coverage);
  }
  var_coverage = var_coverage / lines;
  std_coverage = std::sqrt(var_coverage);
  mean_coverage_reverse_out << mean_coverage << "\t" << std_coverage << std::endl;
  mean_coverage_reverse_out.close();
  ////////////////END REVERSE/////////////////////////////////////////////////////////////////

  auto stop_timer = high_resolution_clock::now();
  auto measure_time = duration_cast<minutes>(stop_timer - start_timer);

  std::cout << "Execution time:  " << measure_time.count() << "  min." << std::endl;

  return 0;
}

void distribution(std::vector<int> &coverages, std::vector<std::pair<int,int>> &counts)
{
  std::sort(coverages.begin(), coverages.end());
  
  std::set<int> temp_set;

  for (auto &el : coverages)
  {
    temp_set.insert(el);
  }

  for (const auto &el : temp_set)
  {
    int mycount = std::count(coverages.begin(), coverages.end(), el);
    auto p_count = std::make_pair(el, mycount);
    counts.push_back(p_count);
  }
}
