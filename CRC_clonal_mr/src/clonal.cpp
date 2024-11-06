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

bool select_clonal_mut(double clonal_probability, double endpoint_mut_freq, double endpoint_coverage, double lambda1, double lambda2);

int main (int argc, char* argv[])
{

  if (argc != 8)
  {
    std::cout << "Use: " << argv[0] << " [file_reeds_p2.dat] [file_no_mutations_p2.dat] [file_reeds_p3.dat] [file_no_mutation_p3.dat] [output_path_mr] [clone number] [brother number]" << std::endl; 
    return 1;
  }

  double division_number = 0.0;

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

  std::string path_file_reeds_p2 = std::string{argv[1]};
  std::string path_file_selected_bases_p2 = std::string{argv[2]};
  std::string path_file_reeds_p3 = std::string{argv[3]};
  std::string path_file_selected_bases_p3 = std::string{argv[4]};
  std::string path_file_mr = std::string{argv[5]};
  std::string file_mut_rate = "file_mr_clonal.dat";
  std::string file_n_bases = "file_n_mutated_bases.dat";

  //1307
  ////// 
  if (std::strcmp(argv[6],"1307") == 0 && std::strcmp(argv[7],"9") == 0)
  {
    division_number = 221.3449; 
  }
  
  if (std::strcmp(argv[6],"1307") == 0 && std::strcmp(argv[7],"2") == 0)
  {
    division_number = 215.1366; 
  }
  
  if (std::strcmp(argv[6],"1307") == 0 && std::strcmp(argv[7],"8") == 0)
  {
    division_number = 229.6166; 
  }
  
  //1502
  ////// 
  if (std::strcmp(argv[6],"1502") == 0 && std::strcmp(argv[7],"3") == 0)
  {
    division_number = 531.525687; 
  }
  
  if (std::strcmp(argv[6],"1502") == 0 && std::strcmp(argv[7],"8") == 0)
  {
    division_number = 395.0971374; 
  }
  
  if (std::strcmp(argv[6],"1502") == 0 && std::strcmp(argv[7],"9") == 0)
  {
    division_number = 466.04; 
  }
  
  if (std::strcmp(argv[6],"1502") == 0 && std::strcmp(argv[7],"10") == 0)
  {
    division_number = 505.3314122; 
  }

  //0282
  ////// 
  if (std::strcmp(argv[6],"0282") == 0 && std::strcmp(argv[7],"5") == 0)
  {
    division_number = 342.159973974;
  }
  
  if (std::strcmp(argv[6],"0282") == 0 && std::strcmp(argv[7],"7") == 0)
  {
    division_number = 351.519994799;
  }


  std::cout << "b * days: " << division_number << std::endl;
  
  auto start_timer = high_resolution_clock::now();

  std::ifstream in_file_list_p2(path_file_reeds_p2, std::ios_base::binary);
  if (!in_file_list_p2.is_open())
  {
    std::cerr << "file_reeds_p2.dat not opened. Exit" << std::endl;
    return 1;
  }
 
  std::ifstream in_file_n_bases_p2(path_file_selected_bases_p2, std::ios_base::binary);
  if (!in_file_n_bases_p2.is_open())
  {
    std::cerr << "file_no_mutations_bases_p2.dat not opened. Exit" << std::endl;
    return 1;
  }

  std::ifstream in_file_list_p3(path_file_reeds_p3, std::ios_base::binary);
  if (!in_file_list_p3.is_open())
  {
    std::cerr << "file_reeds_p3.dat not opened. Exit" << std::endl;
    return 1;
  }
 
  std::ifstream in_file_n_bases_p3(path_file_selected_bases_p3, std::ios_base::binary);
  if (!in_file_n_bases_p3.is_open())
  {
    std::cerr << "file_no_mutations_bases_p3.dat not opened. Exit" << std::endl;
    return 1;
  }

  std::ofstream out_mut_rate(path_file_mr+file_mut_rate, std::ios_base::binary);
  if (!out_mut_rate.is_open())
  {
    std::cerr << "File mutation_rate.dat not opened. Exit." << std::endl;
    return 1;
  }

  std::ofstream out_n_bases(path_file_mr+file_n_bases, std::ios_base::binary);
  if (!out_n_bases.is_open())
  {
    std::cerr << "File n_mutated_bases.dat not opened. Exit." << std::endl;
    return 1;
  }

  size_t count_ll = 0;
  std::string endpoint_line;
  std::string endpoint_coverage;
  std::string endpoint_mutations;
  std::string endpoint_ploidy;

  while (std::getline(in_file_list_p2, endpoint_line))
  { 
    count_ll += 1;
    std::stringstream data_endpoint(endpoint_line);
    data_endpoint >> endpoint_coverage >> endpoint_mutations >> endpoint_ploidy;
    double coverage = std::stod(endpoint_coverage);
    double mutations = std::stod(endpoint_mutations);
    double ploidy = std::stod(endpoint_ploidy);
    double clonal_probability = 1.0/ploidy; 
    double f_hat = mutations * (1.0/coverage);

    bool selected_clonal_01 = select_clonal_mut(clonal_probability, f_hat, coverage, lambda_sigma01, lambda_sigma01);
    bool selected_clonal_03 = select_clonal_mut(clonal_probability, f_hat, coverage, lambda_sigma03, lambda_sigma03);
    bool selected_clonal_1 = select_clonal_mut(clonal_probability, f_hat, coverage, lambda_sigma1, lambda_sigma1);
    bool selected_clonal_2 = select_clonal_mut(clonal_probability, f_hat, coverage, lambda_sigma2, lambda_sigma2);
    bool selected_clonal_3 = select_clonal_mut(clonal_probability, f_hat, coverage, lambda_sigma3, lambda_sigma3);
    bool selected_clonal_4 = select_clonal_mut(clonal_probability, f_hat, coverage, lambda_sigma4, lambda_sigma4);
    bool selected_clonal_5 = select_clonal_mut(clonal_probability, f_hat, coverage, lambda_sigma5, lambda_sigma5);
    bool selected_clonal_6 = select_clonal_mut(clonal_probability, f_hat, coverage, lambda_sigma6, lambda_sigma6);
    
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

    if ((count_ll % 1'000'000) == 0)
    {
      std::cout << count_ll << std::endl; 
    }
  }

  //Ploidy 3
  //
  while (std::getline(in_file_list_p3, endpoint_line))
  { 
    count_ll += 1;
    std::stringstream data_endpoint(endpoint_line);
    data_endpoint >> endpoint_coverage >> endpoint_mutations >> endpoint_ploidy;
    double coverage = std::stod(endpoint_coverage);
    double mutations = std::stod(endpoint_mutations);
    double ploidy = std::stod(endpoint_ploidy);
    double clonal_probability = 1.0/ploidy; 
    double f_hat = mutations * (1.0/coverage);

    bool selected_clonal_01 = select_clonal_mut(clonal_probability, f_hat, coverage, lambda_sigma01, lambda_sigma01);
    bool selected_clonal_03 = select_clonal_mut(clonal_probability, f_hat, coverage, lambda_sigma03, lambda_sigma03);
    bool selected_clonal_1 = select_clonal_mut(clonal_probability, f_hat, coverage, lambda_sigma1, lambda_sigma1);
    bool selected_clonal_2 = select_clonal_mut(clonal_probability, f_hat, coverage, lambda_sigma2, lambda_sigma2);
    bool selected_clonal_3 = select_clonal_mut(clonal_probability, f_hat, coverage, lambda_sigma3, lambda_sigma3);
    bool selected_clonal_4 = select_clonal_mut(clonal_probability, f_hat, coverage, lambda_sigma4, lambda_sigma4);
    bool selected_clonal_5 = select_clonal_mut(clonal_probability, f_hat, coverage, lambda_sigma5, lambda_sigma5);
    bool selected_clonal_6 = select_clonal_mut(clonal_probability, f_hat, coverage, lambda_sigma6, lambda_sigma6);
    
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

    if ((count_ll % 1'000'000) == 0)
    {
      std::cout << count_ll << std::endl; 
    }
  }

  std::string n_bases_string_p2;
  std::string endpoint_n_bases_line_p2;
  std::getline(in_file_n_bases_p2, endpoint_n_bases_line_p2);
  std::stringstream data_endpoint_n_bases_p2(endpoint_n_bases_line_p2);
  data_endpoint_n_bases_p2 >> n_bases_string_p2;
  size_t n_bases_p2 = std::stoll(n_bases_string_p2);
  std::cout << "Selected bases ploidy 2 = " << n_bases_p2 << std::endl;
  std::cout << std::endl;

  std::string n_bases_string_p3;
  std::string endpoint_n_bases_line_p3;
  std::getline(in_file_n_bases_p3, endpoint_n_bases_line_p3);
  std::stringstream data_endpoint_n_bases_p3(endpoint_n_bases_line_p3);
  data_endpoint_n_bases_p3 >> n_bases_string_p3;
  size_t n_bases_p3 = std::stoll(n_bases_string_p3);
  std::cout << "Selected bases ploidy 3 = " << n_bases_p3 << std::endl;
  std::cout << std::endl;

  size_t n_bases_tot = n_bases_p2 + n_bases_p3 + count_ll; 
  std::cout << "Selected bases tot = " << n_bases_tot << std::endl;

  double rate_mut_sigma01 = static_cast<double>(mutated_bases_sigma01)/(static_cast<double>(n_bases_tot) * division_number);
  double rate_mut_sigma03 = static_cast<double>(mutated_bases_sigma03)/(static_cast<double>(n_bases_tot) * division_number);
  double rate_mut_sigma1 = static_cast<double>(mutated_bases_sigma1)/(static_cast<double>(n_bases_tot) * division_number);
  double rate_mut_sigma2 = static_cast<double>(mutated_bases_sigma2)/(static_cast<double>(n_bases_tot) * division_number);
  double rate_mut_sigma3 = static_cast<double>(mutated_bases_sigma3)/(static_cast<double>(n_bases_tot) * division_number);
  double rate_mut_sigma4 = static_cast<double>(mutated_bases_sigma4)/(static_cast<double>(n_bases_tot) * division_number);
  double rate_mut_sigma5 = static_cast<double>(mutated_bases_sigma5)/(static_cast<double>(n_bases_tot) * division_number);
  double rate_mut_sigma6 = static_cast<double>(mutated_bases_sigma6)/(static_cast<double>(n_bases_tot) * division_number);

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

  auto stop_timer = high_resolution_clock::now();
  auto measure_time = duration_cast<minutes>(stop_timer - start_timer);
  std::cout << "Execution time:  " << measure_time.count() << "  min." << std::endl;
  return 0;
}

 bool select_clonal_mut(double clonal_probability, double endpoint_mut_freq, double endpoint_coverage, double lambda1, double lambda2)
{
  bool selection = false;
  double variance = ((1.0 - clonal_probability)*clonal_probability) / endpoint_coverage;
  double lim_inf = (clonal_probability - (lambda1 * std::sqrt(variance)));
  double lim_sup = (clonal_probability + (lambda2 * std::sqrt(variance)));

  if ((endpoint_mut_freq < lim_sup) && (endpoint_mut_freq > lim_inf))
  {
    selection = true;
  }

  return selection;
}
