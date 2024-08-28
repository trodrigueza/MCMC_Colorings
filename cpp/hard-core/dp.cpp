/*////////////////////////////////////////////////////////////////////////////
DP algorithm for computing the number of feasible configurations of the Hard
Core model of a kxk lattice.

Approx time:
  k = 16 -> 1sec
  k = 20 -> 7sec
  k = 22 -> 1min
  k = 24 -> 7min
  k = 25 -> 20min

boost/multiprecision is needed -> https://www.boost.org/users/download/
///////////////////////////////////////////////////////////////////////////*/

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
using namespace boost::multiprecision;

typedef number<cpp_dec_float<60>> big_float; // floats with 40 digits of precision
typedef cpp_int big_int; // arbitrarily large integers

vector<int> validMasks;
vector<vector<int>> compatible;
vector<big_int> dp[2];

bool isValid(int mask) {
  return (mask & (mask >> 1)) == 0;
}

void precompute(int k) {
  int fullMask = (1 << k) - 1;
  for (int mask = 0; mask <= fullMask; mask++) {
    if (isValid(mask)) {
      validMasks.push_back(mask);
    }
  }

  compatible.resize(validMasks.size());
  compatible.resize(validMasks.size());
  for (int i = 0; i < validMasks.size(); i++) {
    for (int j = 0; j < validMasks.size(); j++) {
      if ((validMasks[i] & validMasks[j]) == 0) {
        compatible[i].push_back(j);
      }
    }
  }
}

int main() {
  int k;
  cout << "Input k: ";
  cin >> k;

  clock_t start = clock();

  validMasks.clear();
  compatible.clear();

  precompute(k);

  dp[0].assign(validMasks.size(), 0);
  dp[1].assign(validMasks.size(), 0);

  int curr = 0, prev = 1;

  // Base cases
  for (int i = 0; i < validMasks.size(); i++) {
    dp[curr][i] = 1;
  }

  // Transitions
  for (int row = 1; row < k; row++) {
    fill(dp[prev].begin(), dp[prev].end(), 0);
    swap(curr, prev);

    for (int old_idx = 0; old_idx < validMasks.size(); old_idx++) {
      if (dp[prev][old_idx] == 0) continue;

      for (int new_idx : compatible[old_idx]) {
        dp[curr][new_idx] += dp[prev][old_idx];
      }
    }
  }

  // Calculate answer
  big_int total = 0;
  for (int i = 0; i < validMasks.size(); i++) {
    total += dp[curr][i];
  }

  clock_t end = clock();
  double time_taken = double(end - start) / CLOCKS_PER_SEC;

  cout << "Number of feasible Hard Core configurations for a " << k << "x" << k << " lattice: " << setprecision(5) << total.convert_to<big_float>() << "\n";
  // Uncomment for outputing answer without truncating:
  // cout << "Number of feasible Hard Core configurations for a " << k << "x" << k << " lattice: " << total << "\n";
  cout << "\nfinished in ";
  cerr << setprecision(5) << time_taken << " seconds.\n";
}

// g++ -std=c++20 -O3 dp.cpp -I/opt/homebrew/Cellar/boost/1.85.0/include -L/opt/homebrew/Cellar/boost/1.85.0/lib -lboost_system -o dp
