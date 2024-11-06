#include <iostream>
#include <cassert>
#include <iomanip>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std::chrono;

bool select_clonal_hap(double clonal_probability, double endpoint_mut_freq, double error,  double lambda);

int main (int argc, char* argv[])
{

  if (argc != 6)
  {
    std::cout << "Use: " << argv[0] << " [file_reeds_p1.dat] [file_no_mutation_p1.dat] [output_path_mr] [modal_coverage] [generations]" << std::endl; 
    return 1;
  }

  constexpr double seq_error = 0.01;
  //core variables
  double clonal_probability = 1.0;
  double modal_coverage = std::stod(argv[4]);
  double cov_error = 1.0/modal_coverage;
  double generations = std::stod(argv[5]);

  double error {(cov_error > seq_error) ? cov_error : seq_error};

  double lambda_sigma01 = 0.1;
  double lambda_sigma03 = 0.3;
  double lambda_sigma1 = 1.0;
  double lambda_sigma2 = 2.0;
  double lambda_sigma3 = 3.0;
  double lambda_sigma4 = 4.0;
  double lambda_sigma5 = 5.0;
  double lambda_sigma6 = 6.0;

  size_t mutated_bases_sigma01 = 0;
  size_t mutated_bases_sigma03 = 0;
  size_t mutated_bases_sigma1 = 0;
  size_t mutated_bases_sigma2 = 0;
  size_t mutated_bases_sigma3 = 0;
  size_t mutated_bases_sigma4 = 0;
  size_t mutated_bases_sigma5 = 0;
  size_t mutated_bases_sigma6 = 0;
 
  std::ifstream in_file_list_p1(std::string{argv[1]}, std::ios_base::binary);
  if (!in_file_list_p1.is_open())
  {
    std::cerr << "file_reeds_p1.dat not opened. Exit" << std::endl;
    return 1;
  }
 
  std::ifstream in_file_n_bases_p1(std::string{argv[2]}, std::ios_base::binary);
  if (!in_file_n_bases_p1.is_open())
  {
    std::cerr << "file_no_mutations_bases_p1.dat not opened. Exit" << std::endl;
    return 1;
  }

  std::ofstream out_mut_rate(std::string{argv[3]}+"file_mr_clonal.dat", std::ios_base::binary);
  if (!out_mut_rate.is_open())
  {
    std::cerr << "File mr_clonal.dat not opened. Exit." << std::endl;
    return 1;
  }

  std::ofstream out_n_bases(std::string{argv[3]}+"file_n_mutated_bases.dat", std::ios_base::binary);
  if (!out_n_bases.is_open())
  {
    std::cerr << "File n_bases.dat not opened. Exit." << std::endl;
    return 1;
  }

  std::cout << "Number of generations = " << generations << std::endl;
  std::cout << "Error = " << error << std::endl;

  auto start_timer = high_resolution_clock::now();

  size_t count_ll = 0;
  std::string endpoint_line;
  std::string endpoint_coverage;
  std::string endpoint_mutations;
  std::string endpoint_ploidy;
  /////////////////////////////////////////////////////////////////////////////////////////////////////

  //loop over endpoint
  while (std::getline(in_file_list_p1, endpoint_line))
  {

    count_ll += 1;
    std::stringstream data_endpoint(endpoint_line);
    data_endpoint >> endpoint_coverage >> endpoint_mutations >> endpoint_ploidy;
    double coverage = std::stod(endpoint_coverage);
    double mutations = std::stod(endpoint_mutations);
    double f_hat = mutations * (1.0/coverage);

    bool selected_clonal_01 = select_clonal_hap(clonal_probability, f_hat, error, lambda_sigma01);
    bool selected_clonal_03 = select_clonal_hap(clonal_probability, f_hat, error, lambda_sigma03);
    bool selected_clonal_1 = select_clonal_hap(clonal_probability, f_hat, error, lambda_sigma1);
    bool selected_clonal_2 = select_clonal_hap(clonal_probability, f_hat, error, lambda_sigma2);
    bool selected_clonal_3 = select_clonal_hap(clonal_probability, f_hat, error, lambda_sigma3);
    bool selected_clonal_4 = select_clonal_hap(clonal_probability, f_hat, error, lambda_sigma4);
    bool selected_clonal_5 = select_clonal_hap(clonal_probability, f_hat, error, lambda_sigma5);
    bool selected_clonal_6 = select_clonal_hap(clonal_probability, f_hat, error, lambda_sigma6);

    if (selected_clonal_01)
    {
      mutated_bases_sigma01 += 1;
    }

    if (selected_clonal_03)
    {
      mutated_bases_sigma03 += 1;
    }
    
    if (selected_clonal_1)
    {
      mutated_bases_sigma1 += 1;
    }
    
    if (selected_clonal_2)
    {
      mutated_bases_sigma2 += 1;
    }
    
    if (selected_clonal_3)
    {
      mutated_bases_sigma3 += 1;
    }

    if (selected_clonal_4)
    {
      mutated_bases_sigma4 += 1;
    }

    if (selected_clonal_5)
    {
      mutated_bases_sigma5 += 1;
    }

    if (selected_clonal_6)
    {
      mutated_bases_sigma6 += 1;
    }

    if (count_ll % 10000000 == 0)
    {
      std::cout << count_ll << std::endl; 
    }
  }

  std::string n_bases_string_p1;
  std::string endpoint_n_bases_line_p1;
  std::getline(in_file_n_bases_p1, endpoint_n_bases_line_p1);
  std::stringstream data_endpoint_n_bases_p1(endpoint_n_bases_line_p1);
  data_endpoint_n_bases_p1 >> n_bases_string_p1;
  size_t n_bases_p1 = std::stoll(n_bases_string_p1);
  std::cout << "Selected bases ploidy 1 = " << n_bases_p1 << std::endl;
  std::cout << std::endl;

  size_t n_bases_tot = n_bases_p1 + count_ll; 
  std::cout << "Selected bases tot = " << n_bases_tot << std::endl;

  double rate_mut_sigma01 = static_cast<double>(mutated_bases_sigma01)/(static_cast<double>(n_bases_tot) * generations);
  double rate_mut_sigma03 = static_cast<double>(mutated_bases_sigma03)/(static_cast<double>(n_bases_tot) * generations);
  double rate_mut_sigma1 = static_cast<double>(mutated_bases_sigma1)/(static_cast<double>(n_bases_tot) * generations);
  double rate_mut_sigma2 = static_cast<double>(mutated_bases_sigma2)/(static_cast<double>(n_bases_tot) * generations);
  double rate_mut_sigma3 = static_cast<double>(mutated_bases_sigma3)/(static_cast<double>(n_bases_tot) * generations);
  double rate_mut_sigma4 = static_cast<double>(mutated_bases_sigma4)/(static_cast<double>(n_bases_tot) * generations);
  double rate_mut_sigma5 = static_cast<double>(mutated_bases_sigma5)/(static_cast<double>(n_bases_tot) * generations);
  double rate_mut_sigma6 = static_cast<double>(mutated_bases_sigma6)/(static_cast<double>(n_bases_tot) * generations);

  out_mut_rate << std::fixed << std::setprecision(12) << rate_mut_sigma01 << std::endl;
  out_mut_rate << std::fixed << std::setprecision(12) << rate_mut_sigma03 << std::endl;
  out_mut_rate << std::fixed << std::setprecision(12) << rate_mut_sigma1 << std::endl;
  out_mut_rate << std::fixed << std::setprecision(12) << rate_mut_sigma2 << std::endl;
  out_mut_rate << std::fixed << std::setprecision(12) << rate_mut_sigma3 << std::endl;
  out_mut_rate << std::fixed << std::setprecision(12) << rate_mut_sigma4 << std::endl;
  out_mut_rate << std::fixed << std::setprecision(12) << rate_mut_sigma5 << std::endl;
  out_mut_rate << std::fixed << std::setprecision(12) << rate_mut_sigma6 << std::endl;

  out_n_bases << std::fixed << std::setprecision(12) << mutated_bases_sigma01 << std::endl;
  out_n_bases << std::fixed << std::setprecision(12) << mutated_bases_sigma03 << std::endl;
  out_n_bases << std::fixed << std::setprecision(12) << mutated_bases_sigma1 << std::endl;
  out_n_bases << std::fixed << std::setprecision(12) << mutated_bases_sigma2 << std::endl;
  out_n_bases << std::fixed << std::setprecision(12) << mutated_bases_sigma3 << std::endl;
  out_n_bases << std::fixed << std::setprecision(12) << mutated_bases_sigma4 << std::endl;
  out_n_bases << std::fixed << std::setprecision(12) << mutated_bases_sigma5 << std::endl;
  out_n_bases << std::fixed << std::setprecision(12) << mutated_bases_sigma6 << std::endl;

  out_mut_rate.close();
  out_n_bases.close();

  auto stop_timer = high_resolution_clock::now();
  auto measure_time = duration_cast<minutes>(stop_timer - start_timer);

  std::cout << "Execution time:  " << measure_time.count() << "  min." << std::endl;

  return 0;
}


bool select_clonal_hap(double clonal_probability, double endpoint_mut_freq, double error,  double lambda)
{
  bool selection = false;

  double lim_inf = (clonal_probability - (lambda * error));

  if ((endpoint_mut_freq <= clonal_probability) && (endpoint_mut_freq > lim_inf))
  {
    selection = true;
  }

  return selection;
}

