#include <bits/stdc++.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
using namespace boost::multiprecision;

typedef number<cpp_dec_float<40>> big_float; // floats with 40 digits of precision
typedef cpp_int big_int; // arbitrarely large integers
typedef pair<int, int> pii;

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

inline int random_choice(const int unused[], int size) {
  return unused[uniform_int_distribution<int>(0, size - 1)(rng)];
}

// Initiate lattice with random colors
vector<vector<int>> create_lattice(int k, int q) {
  vector<vector<int>> lattice(k, vector<int>(k));
  for (auto &row : lattice) {
    for (auto &cell : row) {
      cell = uniform_int_distribution<int>(0, q - 1)(rng);
    }
  }
  return lattice;
}

// Generate lattice edges
vector<pair<pii, pii>> gen_edges(int k) {
  vector<pair<pii, pii>> edges;
  edges.reserve(2 * k * (k - 1));
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < k; ++j) {
      if (i + 1 < k) edges.emplace_back(pii{i, j}, pii{i + 1, j});
      if (j + 1 < k) edges.emplace_back(pii{i, j}, pii{i, j + 1});
    }
  }
  return edges;
}

// A systematic gibbs sampler step
void gibbs_sampler_step(int k, int q, vector<vector<int>> &lattice, const vector<vector<vector<pii>>> &neighbours) {
  int used[q], unused[q], version = 0, unused_count = 0;

  for (int i = 0; i < k; i++) {
    for (int j = 0; j < k; j++) {
      version++;
      unused_count = 0;

      for (const auto &ne : neighbours[i][j])
        used[lattice[ne.first][ne.second]] = version;

      for (int color = 0; color < q; color++)
        if (used[color] != version)
          unused[unused_count++] = color;

      lattice[i][j] = (unused_count == 0 ? lattice[i][j] : random_choice(unused, unused_count));
    }
  }
}

// Estimate values of the telescopic product
big_float estimate_ratio(int k, int q, int num_simulations, int num_gibbs_steps, const pair<pii, pii> &edge, vector<vector<int>> &lattice, const vector<vector<vector<pii>>> &neighbours) {
  int u_y = edge.first.first, u_x = edge.first.second;
  int v_y = edge.second.first, v_x = edge.second.second;
  int count = 0;

  for (int i = 0; i < num_simulations; i++) {
    for (int j = 0; j < num_gibbs_steps; j++)
      gibbs_sampler_step(k, q, lattice, neighbours);
    count += lattice[u_y][u_x] != lattice[v_y][v_x];
  }
  big_float ratio = big_float(count) / num_simulations;
  cout << ratio << "\n";
  return ratio;
}

int main() {
  ios::sync_with_stdio(0);
  cin.tie(0);

  int k, q;
  big_float epsilon = 0.5;
  // cin >> k >> q >> epsilon;

  k = 5; q = 7;
  int n = k * k;
  auto lattice = create_lattice(k, q);
  auto edges = gen_edges(k);
  vector<vector<vector<pii>>> neighbours(k, vector<vector<pii>>(k));

  // int num_simulations = static_cast<int>(48 * pow(big_float(4), 3) * pow(big_float(n), 3) / (epsilon * epsilon));
  int num_simulations = static_cast<int>(pow(big_float(n), 3) / (epsilon * epsilon));
  // int num_gibbs_steps = abs(static_cast<int>(n * ((2 * log(n) + log(1 / epsilon) + log(8)) / log(big_float(q) / 32) + 1)));
  int num_gibbs_steps = static_cast<int>(n * (log(n) + log(1 / epsilon)));
  cout << "sims: " << num_simulations << ", steps: " << num_gibbs_steps << "\n";

  big_float Z = pow(big_float(q), n);
  cout << big_int(Z) << "\n";

  for (const auto& edge : edges) {
      const auto &x = edge.first;
      const auto &y = edge.second;
      cout << "(" << x.first << "," << x.second << ")" << " <-> " << "(" << y.first << "," << y.second << ")\n";

      big_float ratio = estimate_ratio(k, q, num_simulations, num_gibbs_steps, edge, lattice, neighbours);
      Z *= ratio;
      cout << big_int(Z) << "\n";

      neighbours[x.first][x.second].push_back(y);
      neighbours[y.first][y.second].push_back(x);
  }

  big_float half(0.5);
  big_int rounded_Z = (Z + half).convert_to<big_int>();
  cout << "Estimated number of " << q << "-colorings for " << k << "x" << k << " lattice: " << rounded_Z << "\n";

  cerr << "\nfinished in " << clock() * 1.0 / CLOCKS_PER_SEC << " sec\n";
}
// g++ -std=c++20 -O3 main.cpp -I/opt/homebrew/Cellar/boost/1.85.0/include -L/opt/homebrew/Cellar/boost/1.85.0/lib -lboost_system -o main
