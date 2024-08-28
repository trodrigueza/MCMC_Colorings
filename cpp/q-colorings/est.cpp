#include <iomanip>
#include <random>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
using namespace boost::multiprecision;

typedef number<cpp_dec_float<60>> big_float; // floats with 40 digits of precision
typedef cpp_int big_int; // arbitrarily large integers
typedef pair<int, int> pii; // structure to store graph coordinates

big_float ratios[11] = {0.80559, 0.83633, 0.85889, 0.87616, 0.8897, 0.90054, 0.90948, 0.91697, 0.92335, 0.92874, 0.93349};

int main() {
  cout << "k\\q " << "5 " << "6 " << "7 " << "8 " << "9 " << "10 " << "11 " << "12 " << "13 " << "14 " << "15\n";
  for (int k = 13; k <= 20; k++) {
    int n = k*k;
    int m = 2*k*(k-1);
    cout << k << " ";
    for (int q = 5; q <= 15; q++) {
      big_float est = pow(big_float(q), n);
      est *= pow(ratios[q-5], m);
      big_int rounded = (est + 0.5).convert_to<big_int>();
      cout << "\"" << rounded.convert_to<big_float>() << "\" ";
    }
    cout << "\n";
  }
}
// g++ -std=c++20 -O3 est14..20.cpp -I/opt/homebrew/Cellar/boost/1.85.0/include -L/opt/homebrew/Cellar/boost/1.85.0/lib -lboost_system -o est
