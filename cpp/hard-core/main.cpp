/*////////////////////////////////////////////////////////////////////////////
MCMC algorithm for estimating the number of feasible Hard-Core configurations
of a kxk lattice.

boost/multiprecision is needed -> https://www.boost.org/users/download/
///////////////////////////////////////////////////////////////////////////*/

#include <random>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>

using namespace std;
using namespace boost::multiprecision;

typedef number<cpp_dec_float<60>> big_float; // floats with 40 digits of precision
typedef cpp_int big_int; // arbitrarily large integers
typedef pair<int, int> pii; // structure to store graph coordinates

const int K = 25; // maximum lattice size (KxK)
const int MAX_EDGES = 2 * K * (K - 1); // maximum number of edges in the lattice
const int MAX_COLORS = 2; // maximum number of colors

// Random number generator
unsigned seed = 1234;
mt19937 rng(seed);
// random_device rd;
// mt19937 rng(rd());

// Initiate lattice with all 0s
void create_lattice(int k, int lattice[K][K]) {
  for (int i = 0; i < k; i++)
    for (int j = 0; j < k; j++)
      lattice[i][j] = 0;
}

// Generate lattice edges
void gen_edges(int k, pair<pii, pii> edges[MAX_EDGES], int &edge_count) {
  edge_count = 0;
  for (int i = 0; i < k; ++i)
    for (int j = 0; j < k; ++j) {
      if (i + 1 < k) edges[edge_count++] = make_pair(pii{i, j}, pii{i + 1, j});
      if (j + 1 < k) edges[edge_count++] = make_pair(pii{i, j}, pii{i, j + 1});
    }
}

// Systematic gibbs sampler (k^2 steps)
void gibbs_sampler_step(int k, int lattice[K][K], pii neighbours[K][K][4], int neighbour_count[K][K]) {
  int f = 1;
  for (int i = 0; i < k; i++)
    for (int j = 0; j < k; j++) {
      f = 1;
      for (int n = 0; n < neighbour_count[i][j]; n++) {
        pii ne = neighbours[i][j][n];
        if (lattice[ne.first][ne.second] == 1) {
          f = 0;
          break;
        }
      }

      lattice[i][j] = (f == 1 ? uniform_int_distribution<>(0, 1)(rng) : 0);
    }
}

// Estimate values of the telescopic product
big_float estimate_ratio(int k, int num_simulations, int num_gibbs_steps, const pair<pii, pii> &edge, int lattice[K][K], pii neighbours[K][K][4], int neighbour_count[K][K]) {
  int u_y = edge.first.first, u_x = edge.first.second;
  int v_y = edge.second.first, v_x = edge.second.second;
  int count = 0;
  for (int i = 0; i < num_simulations; i++) {
    for (int j = 0; j < ceil(num_gibbs_steps/(k*k))+1; j++)
      gibbs_sampler_step(k, lattice, neighbours, neighbour_count);
    count += (lattice[u_y][u_x] == 0 || lattice[v_y][v_x] == 0);
  }
  return big_float(count) / num_simulations;
}

int32_t main() {
  int k;
  big_float epsilon;
  cout << "Enter: k eps\n";
  cin >> k >> epsilon;

  int n = k * k; // |V|
  int m; // |E|
  int lattice[K][K];
  create_lattice(k, lattice);
  pair<pii, pii> edges[MAX_EDGES]; // store edges
  gen_edges(k, edges, m); // generate edges
  cout << "n = " << n << ", " << "m = " << m <<".\n";
  // int num_simulations = static_cast<int>(ceil(48 * pow(big_float(4), 3) * pow(big_float(n), 3) / (epsilon * epsilon)));
  // int num_gibbs_steps = max(3, abs(static_cast<int>(ceil(n * ((2 * log(n) + log(1 / epsilon) + log(8)) / log(big_float(q) / 32) + 1)))));
  int num_simulations = static_cast<int>(pow(big_float(n), 3) / (epsilon * epsilon));
  int num_gibbs_steps = max(n, static_cast<int>(n * (log(n) + log(1 / epsilon)) * (1 / abs(log(2.0/64)))));
  cout << "Simulations: " << num_simulations << ", sampler steps: " << num_gibbs_steps << ".\n";

  pii neighbours[K][K][4]; // for storing vertex neighbours
  int neighbour_count[K][K] = {0}; // for storing number of neighbours for each vertex

  big_float Z = pow(big_float(2), n);
  big_float meanRatio = 0;

  for (int i = 0; i < m; i++) {
    const auto &edge = edges[i];
    const auto &x = edge.first;
    const auto &y = edge.second;

    big_float ratio = estimate_ratio(k, num_simulations, num_gibbs_steps, edge, lattice, neighbours, neighbour_count);
    Z *= ratio;
    meanRatio += ratio;

    neighbours[x.first][x.second][neighbour_count[x.first][x.second]++] = y;
    neighbours[y.first][y.second][neighbour_count[y.first][y.second]++] = x;
  }

  meanRatio /= m;
  big_int rounded_Z = (Z + 0.5).convert_to<big_int>();
  big_float r_Z = rounded_Z.convert_to<big_float>();

  cout << "Estimated number of Hard Core model feasible configurations for " << k << "x" << k << " lattice: " << rounded_Z << "\n";
  cerr << "\nfinished in " << clock() * 1.0 / CLOCKS_PER_SEC << " sec\n";
}

// g++ -std=c++20 -O3 main.cpp -I/opt/homebrew/Cellar/boost/1.85.0/include -L/opt/homebrew/Cellar/boost/1.85.0/lib -lboost_system -o array
