#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "random.h"

Random::Random(std::string primes_in, std::string seed_in)
{
  int seed[4];
  int p1, p2;
  std::ifstream primes(primes_in);

  if (primes.is_open())
  {
    primes >> p1 >> p2 ;
  } 
  else 
  {
   std::cerr << "PROBLEM: Unable to open " << primes_in << std::endl;
   exit(1);
  }
  primes.close();

  std::ifstream input(seed_in);
  std::string property;
  if (input.is_open())
  {
    while ( !input.eof() )
    {
      input >> property;
      if( property == "RANDOMSEED" )
      {
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        set_random(seed, p1, p2);
      }
    }
    input.close();
  } 
  else 
  {
    std::cerr << "PROBLEM: Unable to open " << seed_in << std::endl;
    exit(1);
  }

}

void Random::save_seed(std::string seed_out)
{
  std::ofstream write_seed;
  write_seed.open(seed_out);
  if (write_seed.is_open())
  {
    write_seed << l1 << " " << l2 << " " << l3 << " " << l4 << std::endl;;
  } 
  else 
  {
    std::cerr << "PROBLEM: Unable to open " << seed_out << std::endl;
  }
  write_seed.close();
  return;
}

double Random::rannyu(double min, double max)
{
  return min+(max-min)*rannyu();
}

long double Random::rannyu(long double min, long double max)
{
  return min+(max-min)*static_cast<long double>(rannyu());
}

size_t Random::rannyu_uint(double min, double max)
{
  double res = min+(max-min)*rannyu();
  double app = 0.0;
  double app2 = std::modf(res, &app); 
  if (app2 > 0.5)
  {
    app += 1.0;
  }

  return static_cast<size_t>(app);
}

double Random::rannyu(void)
{
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random::set_random(int * s, int p1, int p2)
{
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}


