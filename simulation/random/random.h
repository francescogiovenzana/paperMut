#ifndef RANDOM_H
#define RANDOM_H

#include <string>
#include <vector>
#include <functional>

class Random 
{
private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

public:
  // constructor
  Random(std::string primes_in, std::string seed_in);

  // methods
  void set_random(int * , int, int);
  void save_seed(std::string seed_out);
  double rannyu(void);
  double rannyu(double min, double max);
  long double rannyu(long double min, long double max);
  size_t rannyu_uint(double min, double max);
};

#endif // RANDOM_H
