#include <iostream>
#include <limits>
#include <iomanip>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cmath>
#include <vector>
#include <tuple>
#include <algorithm>

using namespace std::chrono;

int rounded(double val);
bool same_reeds(std::string A_count, std::string G_count, 
                std::string C_count, std::string T_count,
                std::string a_count, std::string g_count, 
                std::string c_count, std::string t_count, 
                int total_coverage);
///////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char* argv[])
{
  //constexpr double mutation_thr_ancestor = 0.0;
  //p value thr
  constexpr double mutation_thr_ancestor = 0.0;
  double alpha = 0.05;
  //core variables
  int coverage_thr_ancestor = 0;
  int std_dev_cov_ancestor = 0;
  int max_cov_a = 0;
  int min_cov_a = 0;
  int coverage_thr_endpoint = 0;
  int std_dev_cov_endpoint = 0;
  int max_cov_e = 0;
  int min_cov_e = 0;

  if (argc != 9)
  {
    std::cout << "Use: " << argv[0] << " [common_lines_ancestor.dat] [common_lines_endpoint.dat] [output_path] [coverage_modal_ancestor] [coverage_std_dev_ancestor] [coverage_modal_endpoint] [coverage_std_dev_endpoint] [ploidy]" << std::endl; 
    return 1;
  }

  std::ofstream file_L(std::string{argv[3]}+"L_correction_p"+std::string{argv[8]}+".dat", std::ios_base::binary);
  if (!file_L.is_open())
  {
    std::cerr << "File file_L.dat not opened. Exit" << std::endl;
    return 1;
  }

  double ploidy_end = std::stod(argv[8]);
  size_t L1 = 0;
  coverage_thr_ancestor = rounded(std::stod(argv[4]));
  std_dev_cov_ancestor = rounded(std::stod(argv[5]));
  min_cov_a = coverage_thr_ancestor - std_dev_cov_ancestor;
  max_cov_a = coverage_thr_ancestor + std_dev_cov_ancestor;

  coverage_thr_endpoint = rounded(std::stod(argv[6]));
  std_dev_cov_endpoint = rounded(std::stod(argv[7]));
  min_cov_e = coverage_thr_endpoint - std_dev_cov_endpoint;
  max_cov_e = coverage_thr_endpoint + std_dev_cov_endpoint;
  ///////////////////////////////////////////////////////////////////////

  std::cout << "PLOIDY = " << ploidy_end << std::endl;
  std::cout << "Modal coverage ancestor = " << coverage_thr_ancestor << std::endl;
  std::cout << "STD coverage ancestor = " << std_dev_cov_ancestor << std::endl;
  std::cout << "Min coverage ancestor = " << min_cov_a << std::endl;
  std::cout << "Max coverage ancestor = " << max_cov_a << std::endl;
  std::cout << "Modal coverage endpoint = " << coverage_thr_endpoint << std::endl;
  std::cout << "STD coverage endpoint = " << std_dev_cov_endpoint << std::endl;
  std::cout << "Min coverage endpoint = " << min_cov_e << std::endl;
  std::cout << "Max coverage endpoint = " << max_cov_e << std::endl;
  std::cout << "L1 = " << L1 << std::endl;
  std::cout << "alpha = " << alpha << std::endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  auto start_timer = high_resolution_clock::now(); 

  std::string A_s = "A"; 
  std::string G_s = "G";
  std::string C_s = "C";
  std::string T_s = "T";

  //ancestor
  std::string ancestor_line;
  std::string chr_a;
  std::string base_num_a;
  std::string ref_a;
  std::string A_count_a;
  std::string G_count_a;
  std::string C_count_a;
  std::string T_count_a;
  std::string a_count_a;
  std::string g_count_a;
  std::string c_count_a;
  std::string t_count_a;
  //size_t lines_ancestor = 0;
  int ancestor_coverage = 0;
  int ancestor_forward = 0;
  int ancestor_reverse = 0;

  //endpoint
  std::string endpoint_line;
  std::string chr_e;
  std::string base_num_e;
  std::string ref_e;
  std::string A_count_e;
  std::string G_count_e;
  std::string C_count_e;
  std::string T_count_e;
  std::string a_count_e;
  std::string g_count_e;
  std::string c_count_e;
  std::string t_count_e;
  size_t lines_endpoint = 0;
  int endpoint_coverage = 0;
  int endpoint_forward = 0;
  int endpoint_reverse = 0;
  double ancestor_mut_freq = 0.0;

  std::ifstream in_file_ancestor(std::string{argv[1]}, std::ios_base::binary);
  if (!in_file_ancestor.is_open())
  {
    std::cerr << "File ancestor.dat not opened. Exit" << std::endl;
    return 1;
  }

  std::ifstream in_file_endpoint(std::string{argv[2]}, std::ios_base::binary);
  if (!in_file_endpoint.is_open())
  {
    std::cerr << "File endpoint.dat not opened. Exit" << std::endl;
    return 1;
  }
 
  while (std::getline(in_file_endpoint, endpoint_line))
  {
    lines_endpoint += 1;
    //ancestor
    std::stringstream data_endpoint(endpoint_line);
    data_endpoint >> chr_e >> base_num_e >> ref_e >> A_count_e >> G_count_e >> C_count_e >> T_count_e >>
                                                     a_count_e >> g_count_e >> c_count_e >> t_count_e; 

    //ancestor
    std::getline(in_file_ancestor, ancestor_line);
    std::stringstream data_ancestor(ancestor_line);
    data_ancestor >> chr_a >> base_num_a >> ref_a >> A_count_a >> G_count_a >> C_count_a >> T_count_a >>
                                                     a_count_a >> g_count_a >> c_count_a >> t_count_a; 

    if ((ref_a == A_s || ref_a == G_s || ref_a == C_s || ref_a == T_s) && 
        (ref_e == A_s || ref_e == G_s || ref_e == C_s || ref_e == T_s) &&
        ((std::stoi(chr_a) == std::stoi(chr_e)) && (std::stoi(base_num_a) == std::stoi(base_num_e)))) 
    {    
      //endpoint mutation 
      endpoint_forward = std::stoi(A_count_e) + std::stoi(G_count_e) + std::stoi(C_count_e) + std::stoi(T_count_e);
      endpoint_reverse = std::stoi(a_count_e) + std::stoi(g_count_e) + std::stoi(c_count_e) + std::stoi(t_count_e);
      endpoint_coverage = endpoint_forward + endpoint_reverse;
      //ancestor mutation 
      ancestor_forward = std::stoi(A_count_a) + std::stoi(G_count_a) + std::stoi(C_count_a) + std::stoi(T_count_a);
      ancestor_reverse = std::stoi(a_count_a) + std::stoi(g_count_a) + std::stoi(c_count_a) + std::stoi(t_count_a);
      ancestor_coverage = ancestor_forward + ancestor_reverse;
 
      //bool same_reeds_ancestor = same_reeds(A_count_a, G_count_a, C_count_a, T_count_a,
      //                                      a_count_a, g_count_a, c_count_a, t_count_a, ancestor_coverage);
      if (ref_a == A_s)
      {
        int mutated_reads_a = std::max(std::stoi(T_count_a)+std::stoi(t_count_a),
                                       std::stoi(C_count_a)+std::stoi(c_count_a)); 
        mutated_reads_a = std::max(mutated_reads_a, std::stoi(G_count_a)+std::stoi(g_count_a));
        ancestor_mut_freq = static_cast<double>(mutated_reads_a)/static_cast<double>(ancestor_coverage);
      } 

      if (ref_a == G_s)
      {
        int mutated_reads_a = std::max(std::stoi(T_count_a)+std::stoi(t_count_a),
                                       std::stoi(C_count_a)+std::stoi(c_count_a)); 
        mutated_reads_a = std::max(mutated_reads_a, std::stoi(A_count_a)+std::stoi(a_count_a));
        ancestor_mut_freq = static_cast<double>(mutated_reads_a)/static_cast<double>(ancestor_coverage);
      } 

      if (ref_a == C_s)
      {
        int mutated_reads_a = std::max(std::stoi(T_count_a)+std::stoi(t_count_a),
                                       std::stoi(A_count_a)+std::stoi(a_count_a)); 
        mutated_reads_a = std::max(mutated_reads_a, std::stoi(G_count_a)+std::stoi(g_count_a));
        ancestor_mut_freq = static_cast<double>(mutated_reads_a)/static_cast<double>(ancestor_coverage);
      } 

      if (ref_a == T_s)
      {
        int mutated_reads_a = std::max(std::stoi(A_count_a)+std::stoi(a_count_a),
                                       std::stoi(C_count_a)+std::stoi(c_count_a)); 
        mutated_reads_a = std::max(mutated_reads_a, std::stoi(G_count_a)+std::stoi(g_count_a));
        ancestor_mut_freq = static_cast<double>(mutated_reads_a)/static_cast<double>(ancestor_coverage);
      } 
 
      if ((ancestor_coverage >= min_cov_a) && (ancestor_coverage <= max_cov_a) && (ancestor_mut_freq <= mutation_thr_ancestor))
      { 
        if ((endpoint_coverage >= min_cov_e) && (endpoint_coverage <= max_cov_e))
        { 
          L1 += 1;
        }
      }
    }

    if (lines_endpoint % 10000000 == 0)
    {
      std::cout << lines_endpoint << "\t" << chr_e << std::endl; 
    }
  }

  in_file_ancestor.close();
  in_file_endpoint.close();
  file_L << L1 << std::endl;
  file_L.close();

  alpha = alpha/static_cast<double>(L1);

  std::cout << "L1 = " << L1 << std::endl;
  std::cout << "alpha/L1 = " << alpha << std::endl;

  auto stop_timer = high_resolution_clock::now();
  auto measure_time = duration_cast<minutes>(stop_timer - start_timer);

  std::cout << "Execution time:  " << measure_time.count() << "  min." << std::endl;

  return 0;
}

int rounded(double val)
{
  int r_val = 0;
  if ((val > static_cast<double>(std::numeric_limits<int>::min())) && 
      (val < static_cast<double>(std::numeric_limits<int>::max())))
  {
    r_val = static_cast<int>(std::floor(val));
    return r_val;
  }
  else
  {
    exit(1);
  }

  return r_val;
}

bool same_reeds(std::string A_count, std::string G_count, 
                std::string C_count, std::string T_count,
                std::string a_count, std::string g_count, 
                std::string c_count, std::string t_count, 
                int total_coverage)
{
  bool sr = false;

  int ATOT = std::stoi(A_count) + std::stoi(a_count);
  int GTOT = std::stoi(G_count) + std::stoi(g_count);
  int CTOT = std::stoi(C_count) + std::stoi(c_count);
  int TTOT = std::stoi(T_count) + std::stoi(t_count);
  if ((ATOT == total_coverage) || (GTOT == total_coverage) || (CTOT == total_coverage) || (TTOT == total_coverage))
  {
    sr = true;
  }
  return sr;
}

