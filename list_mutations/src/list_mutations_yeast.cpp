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

bool balanced_tot(std::string A_count, std::string G_count, std::string C_count, std::string T_count,
                  std::string a_count, std::string g_count, std::string c_count, std::string t_count);
int rounded(double val);
bool same_reeds(std::string A_count, std::string G_count, 
                std::string C_count, std::string T_count,
                std::string a_count, std::string g_count, 
                std::string c_count, std::string t_count, 
                int total_coverage);
///////////////////////////////////////////////////////////////////////////////////////////////////////
double binomial_probability(int n, int m);
double binomial_p_value(int N, int Nm);

int main (int argc, char* argv[])
{
  constexpr double mutation_thr_ancestor = 0.0;
  //core variables
  int coverage_thr_ancestor = 0;
  int std_dev_cov_ancestor = 0;
  int max_cov_a = 0;
  int min_cov_a = 0;
  int coverage_thr_endpoint = 0;
  int std_dev_cov_endpoint = 0;
  int max_cov = 0;
  int min_cov = 0;

  if (argc != 10)
  {
    std::cout << "Use: " << argv[0] << " [common_lines_control.dat]  [common_lines_evolved.dat]  [output_path]  [coverage_ancestor_modal]  [coverage_ancestor_std_dev]  [coverage_endpoint_modal]  [coverage_endpoint_std_dev]  [copy_number]  [L1]" << std::endl; 
    return 1;
  }

  std::ofstream file_reeds_p2(std::string{argv[3]}+"file_reeds_p"+std::string{argv[8]}+".dat", std::ios_base::binary);
  if (!file_reeds_p2.is_open())
  {
    std::cerr << "File file_reeds.dat not opened. Exit" << std::endl;
    return 1;
  }

  std::ofstream file_sel_p2(std::string{argv[3]}+"file_no_mutations_p"+std::string{argv[8]}+".dat", std::ios_base::binary);
  if (!file_sel_p2.is_open())
  {
    std::cerr << "File file_no_mutations.dat not opened. Exit" << std::endl;
    return 1;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  size_t selected_bases_p2 = 0;
  size_t n_bases_ploidy2_mut = 0;

  auto start_timer = high_resolution_clock::now(); //reference base 
  std::string A_s = "A"; 
  std::string G_s = "G";
  std::string C_s = "C";
  std::string T_s = "T";

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
  int endpoint_coverage = 0;
  int endpoint_forward = 0;
  int endpoint_reverse = 0;
  double endpoint_mut_freq = 0.0;
  std::string endpoint_mutated_in;
  size_t lines_endpoint = 0;

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
  int ancestor_coverage = 0;
  int ancestor_forward = 0;
  int ancestor_reverse = 0;
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

  double ploidy_end = std::stod(argv[8]);
  double L1 = std::stod(argv[9]);
  //double L2 = std::stod(argv[10]);
  double alpha1 = 0.05/L1;
  //double alpha2 = 0.05 / L2;
  coverage_thr_ancestor = rounded(std::stod(argv[4]));
  std_dev_cov_ancestor = rounded(std::stod(argv[5]));
  min_cov_a = coverage_thr_ancestor - std_dev_cov_ancestor;
  max_cov_a = coverage_thr_ancestor + std_dev_cov_ancestor;
  coverage_thr_endpoint = rounded(std::stod(argv[6]));
  std_dev_cov_endpoint = rounded(std::stod(argv[7]));
  min_cov = coverage_thr_endpoint - std_dev_cov_endpoint;
  max_cov = coverage_thr_endpoint + std_dev_cov_endpoint;
  ///////////////////////////////////////////////////////////////////////

  std::cout << "PLOIDY = " << ploidy_end << std::endl;
  std::cout << "Modal coverage ancestor = " << coverage_thr_ancestor << std::endl;
  std::cout << "STD coverage ancestor = " << std_dev_cov_ancestor << std::endl;
  std::cout << "Min coverage ancestor = " << min_cov_a << std::endl;
  std::cout << "Max coverage ancestor = " << max_cov_a << std::endl;
  std::cout << "Modal coverage endpoint = " << coverage_thr_endpoint << std::endl;
  std::cout << "STD coverage endpoint = " << std_dev_cov_endpoint << std::endl;
  std::cout << "Min coverage endpoint = " << min_cov << std::endl;
  std::cout << "Max coverage endpoint = " << max_cov << std::endl;
  std::cout << "L1 = " << L1 << std::endl;
  std::cout << "alpha/L1 = " << alpha1 << std::endl;

  while (std::getline(in_file_endpoint, endpoint_line))
  {
    //endpoint
    int mutated_reads = 0;
    lines_endpoint += 1;
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
      int Nm_endpoint = std::max(endpoint_forward, endpoint_reverse);

      if (ref_e == A_s)
      {

        mutated_reads = std::max(std::stoi(T_count_e)+std::stoi(t_count_e),
                                 std::stoi(C_count_e)+std::stoi(c_count_e)); 
        mutated_reads = std::max(mutated_reads, std::stoi(G_count_e)+std::stoi(g_count_e));
        endpoint_mut_freq = static_cast<double>(mutated_reads)/static_cast<double>(endpoint_coverage);
        if (mutated_reads != 0)
        {
          if (mutated_reads == (std::stoi(G_count_e)+std::stoi(g_count_e))) { endpoint_mutated_in = "G";}
          if (mutated_reads == (std::stoi(C_count_e)+std::stoi(c_count_e))) { endpoint_mutated_in = "C";}
          if (mutated_reads == (std::stoi(T_count_e)+std::stoi(t_count_e))) { endpoint_mutated_in = "T";}
        }
        else
        {
          endpoint_mutated_in = "NN";
        }
      } 

      if (ref_e == G_s)
      {
        mutated_reads = std::max(std::stoi(T_count_e)+std::stoi(t_count_e),
                                 std::stoi(C_count_e)+std::stoi(c_count_e)); 
        mutated_reads = std::max(mutated_reads, std::stoi(A_count_e)+std::stoi(a_count_e));
        endpoint_mut_freq = static_cast<double>(mutated_reads)/static_cast<double>(endpoint_coverage);
        if (mutated_reads != 0)
        {
          if (mutated_reads == (std::stoi(A_count_e)+std::stoi(a_count_e))) { endpoint_mutated_in = "A";}
          if (mutated_reads == (std::stoi(C_count_e)+std::stoi(c_count_e))) { endpoint_mutated_in = "C";}
          if (mutated_reads == (std::stoi(T_count_e)+std::stoi(t_count_e))) { endpoint_mutated_in = "T";}
        }
        else
        {
          endpoint_mutated_in = "NN";
        }
      }

      if (ref_e == C_s)
      {
        mutated_reads = std::max(std::stoi(T_count_e)+std::stoi(t_count_e),
                                 std::stoi(A_count_e)+std::stoi(a_count_e)); 
        mutated_reads = std::max(mutated_reads, std::stoi(G_count_e)+std::stoi(g_count_e));
        endpoint_mut_freq = static_cast<double>(mutated_reads)/static_cast<double>(endpoint_coverage);
        if (mutated_reads != 0)
        {
          if (mutated_reads == (std::stoi(A_count_e)+std::stoi(a_count_e))) { endpoint_mutated_in = "A";}
          if (mutated_reads == (std::stoi(G_count_e)+std::stoi(g_count_e))) { endpoint_mutated_in = "G";}
          if (mutated_reads == (std::stoi(T_count_e)+std::stoi(t_count_e))) { endpoint_mutated_in = "T";}
        }
        else
        {
          endpoint_mutated_in = "NN";
        }
      }

      if (ref_e == T_s)
      {
        mutated_reads = std::max(std::stoi(A_count_e)+std::stoi(a_count_e),
                                 std::stoi(C_count_e)+std::stoi(c_count_e)); 
        mutated_reads = std::max(mutated_reads, std::stoi(G_count_e)+std::stoi(g_count_e));
        endpoint_mut_freq = static_cast<double>(mutated_reads)/static_cast<double>(endpoint_coverage);
        if (mutated_reads != 0)
        {
          if (mutated_reads == (std::stoi(A_count_e)+std::stoi(a_count_e))) { endpoint_mutated_in = "A";}
          if (mutated_reads == (std::stoi(G_count_e)+std::stoi(g_count_e))) { endpoint_mutated_in = "G";}
          if (mutated_reads == (std::stoi(C_count_e)+std::stoi(c_count_e))) { endpoint_mutated_in = "C";}
        }
        else
        {
          endpoint_mutated_in = "NN";
        }
      }
      
      //ancestor mutation 
      ancestor_forward = std::stoi(A_count_a) + std::stoi(G_count_a) + std::stoi(C_count_a) + std::stoi(T_count_a);
      ancestor_reverse = std::stoi(a_count_a) + std::stoi(g_count_a) + std::stoi(c_count_a) + std::stoi(t_count_a);
      ancestor_coverage = ancestor_forward + ancestor_reverse;
      int Nm_ancestor = std::max(ancestor_forward, ancestor_reverse);

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
 
      //bool same_reeds_ancestor = same_reeds(A_count_a, G_count_a, C_count_a, T_count_a,
      //                                      a_count_a, g_count_a, c_count_a, t_count_a, ancestor_coverage);
 
      if ((ancestor_coverage >= min_cov_a) && (ancestor_coverage <= max_cov_a) && (ancestor_mut_freq <= mutation_thr_ancestor))
      { 
        if ((endpoint_coverage >= min_cov) && (endpoint_coverage <= max_cov))
        { 
          double p_value_ancestor = binomial_p_value(ancestor_coverage, Nm_ancestor); 
          if (p_value_ancestor >= alpha1)
          //if (balanced_tot(A_count_a, G_count_a, C_count_a, T_count_a, a_count_a, g_count_a, c_count_a, t_count_a))
          {
            double p_value_endpoint = binomial_p_value(endpoint_coverage, Nm_endpoint);
            if (p_value_endpoint >= alpha1)
            //if (balanced_tot(A_count_e, G_count_e, C_count_e, T_count_e, a_count_e, g_count_e, c_count_e, t_count_e))
            {    
              selected_bases_p2 += 1;

              if (endpoint_mut_freq > 0.0)
              {
                n_bases_ploidy2_mut += 1;
                file_reeds_p2 << endpoint_coverage << "\t" << mutated_reads << "\t" << ploidy_end << std::endl;
              }
            } 
          }
        }
      }
    }

    if (lines_endpoint % 10000000 == 0)
    {
      std::cout << lines_endpoint << "\t" << chr_e << std::endl; 
    }
  }

  
  in_file_endpoint.close();
  in_file_ancestor.close();


  size_t no_mutations_p2 = selected_bases_p2 - n_bases_ploidy2_mut;

  std::cout << "Selected bases ploidy 2: " << selected_bases_p2 << std::endl;

  file_sel_p2 << no_mutations_p2 << std::endl;

  file_reeds_p2.close();
  file_sel_p2.close();

  auto stop_timer = high_resolution_clock::now();
  auto measure_time = duration_cast<minutes>(stop_timer - start_timer);

  std::cout << "Execution time:  " << measure_time.count() << "  min." << std::endl;

  return 0;
}

bool balanced_tot(std::string A_count, std::string G_count, std::string C_count, std::string T_count,
                  std::string a_count, std::string g_count, std::string c_count, std::string t_count)
{
  bool balanced = false;
  double forward = std::stod(A_count) + std::stod(G_count) + std::stod(C_count) + std::stod(T_count);
  double reverse = std::stod(a_count) + std::stod(g_count) + std::stod(c_count) + std::stod(t_count);
  double tot = forward + reverse;
  double x = forward/tot; 
  if ((x > (0.5 - 1.0/std::sqrt(tot))) && (x < (0.5 + 1.0/std::sqrt(tot))))
  {
    balanced = true;
  }
  return balanced;
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

double binomial_probability(int n, int m)
{
  double b = 0.0;
  double b_coeff = 0.0;
  double N = static_cast<double>(n);
  double M = static_cast<double>(m);

  for (double k = 1.0; k <= M; k += 1.0)
  {
    b += std::log(N - k + 1.0) - std::log(k);
  }
 
  double log_pmf_k = b + (N * std::log(0.5));
  b_coeff = std::exp(log_pmf_k);
  //std::cout << "Binomial coeff = " << b_coeff << std::endl; 
  return b_coeff; 
}

// Function to calculate one-tailed p-value for the binomial test (right-tailed)
double binomial_p_value(int N, int Nm) 
{
  double p_value = 0.0;

  // For a right-tailed test, sum the probabilities of getting F or more successes
  for (int m = Nm; m <= N; m++) 
  {
    p_value += binomial_probability(N, m);
  }
  return p_value;
}

