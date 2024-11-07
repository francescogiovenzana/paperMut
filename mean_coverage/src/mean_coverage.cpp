#include <iostream>
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

int main (int argc, char* argv[])
{
  
  if (argc != 4)
  {
    std::cout << "Use: " << argv[0] << " [path [ancestor|endpoint]_pseudopileup.dat]  [path [ancestor|endpoint]_mean_coverage.dat]  [max threshold]" << std::endl; 
    return 1;
  }

  int max_thr = std::stoi(argv[3]);

  std::ifstream in_file_(std::string{argv[1]}, std::ios_base::binary);
  if (!in_file_.is_open())
  {
    std::cerr << "File input.dat not opened. Exit" << std::endl;
    return 1;
  }
 
  //output file mean coverage
  std::ofstream file_mean_coverage(std::string{argv[2]}, std::ios_base::binary);
  if (!file_mean_coverage.is_open())
  {
    std::cerr << "File mean_coverage.dat not opened. Exit" << std::endl;
    return 1;
  }

  auto start_timer = high_resolution_clock::now();

  std::cout << "Max coverage used: " << max_thr << std::endl;
  //endpoint
  std::string _line;
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
  int _coverage = 0;
  size_t total_coverage = 0;
  size_t lines_ = 0;

  std::vector<int> single_coverage;
  //loop over endpoint
  while (std::getline(in_file_, _line))
  {
    //endpoint
    lines_ += 1;
    std::stringstream data_(_line);
    data_ >> chr_e >> base_num_e >> ref_e >> A_count_e >> G_count_e >> C_count_e >> T_count_e >>
                                             a_count_e >> g_count_e >> c_count_e >> t_count_e; 


    //endpoint mutation 
    _coverage = std::stoi(A_count_e) + std::stoi(G_count_e) + std::stoi(C_count_e) + std::stoi(T_count_e) +
                std::stoi(a_count_e) + std::stoi(g_count_e) + std::stoi(c_count_e) + std::stoi(t_count_e);

    if (lines_ % 10000000 == 0)
    {
      std::cout << lines_ << "\t" << chr_e << std::endl; 
      std::cout << std::endl;
      std::cout << std::endl;
    }
 
    if (_coverage < max_thr)
    {
      single_coverage.emplace_back(_coverage);
      total_coverage += _coverage;
    }
  }

  double mean_coverage = static_cast<double>(total_coverage) / static_cast<double>(lines_);
  double var_coverage = 0.0;
  double std_coverage = 0.0;

  for (auto& el : single_coverage)
  {
    var_coverage += (static_cast<double>(el) - mean_coverage)*(static_cast<double>(el) - mean_coverage);
  }

  //var
  var_coverage = var_coverage / static_cast<double>(lines_);

  //standard deviation
  std_coverage = std::sqrt(var_coverage);

  std::cout << mean_coverage << "\t" << std_coverage << std::endl;

  file_mean_coverage << mean_coverage << std::endl;
  file_mean_coverage << std_coverage << std::endl;

  std::cout << std::endl;
  file_mean_coverage.close();

  auto stop_timer = high_resolution_clock::now();
  auto measure_time = duration_cast<minutes>(stop_timer - start_timer);

  std::cout << "Execution time:  " << measure_time.count() << "  min." << std::endl;

  return 0;
}

