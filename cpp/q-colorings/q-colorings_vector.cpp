#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <random>

using namespace std;
using namespace boost::multiprecision;

typedef number<cpp_dec_float<40>> big_float; // floats with 40 digits of precision
typedef cpp_int big_int; // arbitrarely large integers
typedef pair<int, int> pii; // structure to store graph coordinates

// Random number generator
unsigned seed = 1234;
mt19937 rng(seed);

// Initiate lattice with random colors
vector<vector<int>> create_lattice(int k, int q) {
  vector<vector<int>> lattice(k, vector<int>(k));
  for (auto &row : lattice)
    for (auto &cell : row)
      cell = uniform_int_distribution<int>(0, q - 1)(rng);

  return lattice;
}

// Generate lattice edges
vector<pair<pii, pii>> gen_edges(int k) {
  vector<pair<pii, pii>> edges;
  edges.reserve(2 * k * (k - 1));
  for (int i = 0; i < k; i++)
    for (int j = 0; j < k; j++) {
      if (i + 1 < k) edges.emplace_back(pii{i, j}, pii{i + 1, j});
      if (j + 1 < k) edges.emplace_back(pii{i, j}, pii{i, j + 1});
    }

  return edges;
}

// Systematic Gibbs sampler (k^2 steps)
void gibbs_sampler_step(int k, int q, vector<vector<int>> &lattice, const vector<vector<vector<pii>>> &neighbours) {
  int used[q], unused[q], version = 0, count = 0;
  for (int i = 0; i < k; i++)
    for (int j = 0; j < k; j++) {
      version++;
      count = 0;
      for (const auto &ne : neighbours[i][j])
        used[lattice[ne.first][ne.second]] = version;

      for (int color = 0; color < q; color++)
        if (used[color] != version)
          unused[count++] = color;

      lattice[i][j] = (count == 0 ? lattice[i][j] : unused[uniform_int_distribution<int>(0, count-1)(rng)]);
    }
}

// Estimate values of the telescopic product
big_float estimate_ratio(int k, int q, int num_simulations, int num_gibbs_steps, const pair<pii, pii> &edge, vector<vector<int>> &lattice, const vector<vector<vector<pii>>> &neighbours) {
  int u_y = edge.first.first, u_x = edge.first.second;
  int v_y = edge.second.first, v_x = edge.second.second;
  int count = 0;

  for (int i = 0; i < num_simulations; i++) {
    for (int j = 0; j < ceil(num_gibbs_steps/(k*k)); j++)
      gibbs_sampler_step(k, q, lattice, neighbours);
    count += lattice[u_y][u_x] != lattice[v_y][v_x];
  }

  return big_float(count) / num_simulations;
}

int main() {
  int k, q;
  big_float epsilon;
  cin >> k >> q >> epsilon; // read k, q and eps from input
  int n = k * k; // |V|
  auto lattice = create_lattice(k, q);
  auto edges = gen_edges(k);
  vector<vector<vector<pii>>> neighbours(k, vector<vector<pii>>(k));

  int num_simulations = static_cast<int>(pow(big_float(n), 3) / (epsilon * epsilon));
  int num_gibbs_steps = abs(static_cast<int>(n * ((2 * log(n) + log(1 / epsilon) + log(8)) / log(big_float(q) / 32) + 1)));
  cout << "Sims: " << num_simulations << ", steps: " << num_gibbs_steps << "\n";

  big_float Z = pow(big_float(q), n); // Z = Z_0

  for (const auto &edge : edges) { // edge = {x, y}
      const auto &x = edge.first;
      const auto &y = edge.second;
      big_float ratio = estimate_ratio(k, q, num_simulations, num_gibbs_steps, edge, lattice, neighbours);
      Z *= ratio;
      // add edge
      neighbours[x.first][x.second].emplace_back(y);
      neighbours[y.first][y.second].emplace_back(x);
  }

  big_int rounded_Z = (Z + 0.5).convert_to<big_int>();
  cout << "Estimated number of " << q << "-colorings for " << k << "x" << k << " lattice: " << rounded_Z << "\n";

  cerr << "\nfinished in " << clock() * 1.0 / CLOCKS_PER_SEC << " sec\n";
}
// g++ -std=c++20 -O3 q-colorings_vector.cpp -I/opt/homebrew/Cellar/boost/1.85.0/include -L/opt/homebrew/Cellar/boost/1.85.0/lib -lboost_system -o vector
