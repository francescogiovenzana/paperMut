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

int main (int argc, char* argv[])
{

  if (argc != 5)
  {
    std::cout << "Use: " << argv[0] << " [ancestor.dat] [endpoint.dat] [common_lines_ancestor.dat] [common_lines_endpoint.dat]" << std::endl; 
    return 1;
  }

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

  //output file common lines
  std::ofstream common_ancestor_out(std::string{argv[3]}, std::ios_base::binary);
  if (!common_ancestor_out.is_open())
  {
    std::cerr << "File common_lines_ancestor.dat not opened. Exit." << std::endl;
    return 1;
  }
 
  std::ofstream common_endpoint_out(std::string{argv[4]}, std::ios_base::binary);
  if (!common_endpoint_out.is_open())
  {
    std::cerr << "File common_lines_endpoint.dat not opened. Exit." << std::endl;
    return 1;
  }

  auto start_timer = std::chrono::high_resolution_clock::now();
  //reference base
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
  bool ancestor_ends = false;

  size_t common_lines = 0;

  std::getline(in_file_ancestor, ancestor_line);
  std::stringstream data_ancestor(ancestor_line);
  data_ancestor >> chr_a >> base_num_a >> ref_a >> A_count_a >> G_count_a >> C_count_a >> T_count_a >>
                                                   a_count_a >> g_count_a >> c_count_a >> t_count_a; 

  while (std::getline(in_file_endpoint, endpoint_line))
  {
    std::stringstream data_endpoint(endpoint_line);
    lines_endpoint += 1;
    data_endpoint >> chr_e >> base_num_e >> ref_e >> A_count_e >> G_count_e >> C_count_e >> T_count_e >>
                                                     a_count_e >> g_count_e >> c_count_e >> t_count_e; 

    while ((std::stoi(chr_a) < std::stoi(chr_e)) || 
           (std::stoi(chr_a) == std::stoi(chr_e) && 
            std::stoi(base_num_a) < std::stoi(base_num_e)))
    {
      if (std::getline(in_file_ancestor, ancestor_line))
      { 
        std::stringstream data_ancestor(ancestor_line);      
        data_ancestor >> chr_a >> base_num_a >> ref_a >> A_count_a >> G_count_a >> C_count_a >> T_count_a >>
                                                         a_count_a >> g_count_a >> c_count_a >> t_count_a; 
      }
      else
      {
        std::cout << "Ancestor ended before Endpoint..." << std::endl;
        std::cout << "Ancestor last line: " << chr_a << "\t" << base_num_a << std::endl;
        ancestor_ends = true;
        break;
      }
    }

    //check ancestor and endpoint lines to be the same -> only COMMON LINES in ancestor and endpoint.
    if ((ref_a == A_s || ref_a == G_s || ref_a == C_s || ref_a == T_s) &&
        (ref_e == A_s || ref_e == G_s || ref_e == C_s || ref_e == T_s) &&
        (std::stoi(chr_a) == std::stoi(chr_e) && std::stoi(base_num_a) == std::stoi(base_num_e))) 
    {
      common_lines += 1; 

      common_ancestor_out << chr_a << "\t" << base_num_a << "\t" << ref_a << "\t" << 
                             A_count_a << "\t" << G_count_a << "\t" << C_count_a << "\t" << T_count_a << "\t" << 
                             a_count_a << "\t" << g_count_a << "\t" << c_count_a << "\t" << t_count_a << std::endl;

      common_endpoint_out << chr_e << "\t" << base_num_e << "\t" << ref_e << "\t" << 
                             A_count_e << "\t" << G_count_e << "\t" << C_count_e << "\t" << T_count_e << "\t" << 
                             a_count_e << "\t" << g_count_e << "\t" << c_count_e << "\t" << t_count_e << std::endl;
    }

    if (lines_endpoint % 10000000 == 0)
    {
      std::cout << lines_endpoint << "\t" << chr_e << "\n";
    }

    //exit code because ancestor file has reached EOF
    if (ancestor_ends) { break; }
  }
  
  //exit code because endpoint file has reached EOF 
  if (ancestor_ends == false)
  {
    std::cout << "Endpoint ended before ancestor..." << std::endl;
    std::cout << "Endpoint last line: " << chr_e << "\t" << base_num_e << std::endl;
  }

  std::cout << "Ancestor - Endpoint common lines: " << common_lines << std::endl;

  auto stop_timer = std::chrono::high_resolution_clock::now();
  auto measure_time = std::chrono::duration_cast<std::chrono::minutes>(stop_timer - start_timer);

  std::cout << "Execution time:  " << measure_time.count() << "  min." << std::endl;

  return 0;
}

