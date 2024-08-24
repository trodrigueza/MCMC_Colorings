#include <random>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_int.hpp>
using namespace std;
using namespace boost::multiprecision;

typedef number<cpp_dec_float<40>> big_float; // floats with 40 digits of precision
typedef cpp_int big_int; // arbitrarily large integers
typedef pair<int, int> pii; // structure to store graph coordinates

const int K = 20; // maximum lattice size (KxK)
const int MAX_EDGES = 2 * K * (K - 1); // maximum number of edges in the lattice
const int MAX_COLORS = 15; // maximum number of colors

// Random number generator
unsigned seed = 1234;
mt19937 rng(seed);
// random_device rd;
// mt19937 rng(rd());

// Initiate lattice with random colors
void create_lattice(int k, int q, int lattice[K][K]) {
  for (int i = 0; i < k; i++)
    for (int j = 0; j < k; j++)
      lattice[i][j] = uniform_int_distribution<>(0, q - 1)(rng);
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
void gibbs_sampler_step(int k, int q, int lattice[K][K], pii neighbours[K][K][4], int neighbour_count[K][K]) {
  int used[MAX_COLORS] = {-1}, unused[q], cur = 1, count = 0;
  for (int i = 0; i < k; i++)
    for (int j = 0; j < k; j++) {
      count = 0;
      for (int n = 0; n < neighbour_count[i][j]; n++) {
        pii ne = neighbours[i][j][n];
        used[lattice[ne.first][ne.second]] = cur;
      }
      for (int color = 0; color < q; color++) {
        if (used[color] != cur) {
          unused[count++] = color;
        }
      }
      lattice[i][j] = (count == 0 ? lattice[i][j] : unused[uniform_int_distribution<>(0, count - 1)(rng)]);
      cur++;
    }
}

// Estimate values of the telescopic product
big_float estimate_ratio(int k, int q, int num_simulations, int num_gibbs_steps, const pair<pii, pii> &edge, int lattice[K][K], pii neighbours[K][K][4], int neighbour_count[K][K]) {
  int u_y = edge.first.first, u_x = edge.first.second;
  int v_y = edge.second.first, v_x = edge.second.second;
  int count = 0;

  for (int i = 0; i < num_simulations; i++) {
    for (int j = 0; j < ceil(num_gibbs_steps/(k*k)); j++)
      gibbs_sampler_step(k, q, lattice, neighbours, neighbour_count);
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
  big_float epsilon;
  cin >> k >> q >> epsilon;
  int n = k * k; // |V|
  int m; // |E|

  int lattice[K][K];
  create_lattice(k, q, lattice);
  pair<pii, pii> edges[MAX_EDGES]; // store edges
  gen_edges(k, edges, m); // generate edges
  pii neighbours[K][K][4]; // for storing vertex neighbours
  int neighbour_count[K][K] = {0}; // for storing number of neighbours for each vertex

  // int num_simulations = static_cast<int>(48 * pow(big_float(4), 3) * pow(big_float(n), 3) / (epsilon * epsilon));
  int num_simulations = static_cast<int>(pow(big_float(n), 3) / (epsilon * epsilon));
  // int num_gibbs_steps = abs(static_cast<int>(n * ((2 * log(n) + log(1 / epsilon) + log(8)) / log(big_float(q) / 32) + 1)));
  int num_gibbs_steps = static_cast<int>(n * (log(n) + log(1 / epsilon)));
  cout << "Sims: " << num_simulations << ", steps: " << num_gibbs_steps << "\n";

  big_float Z = pow(big_float(q), n);
  cout << big_int(Z) << "\n";

  for (int i = 0; i < m; i++) {
      const auto &edge = edges[i];
      const auto &x = edge.first;
      const auto &y = edge.second;
      cout << "(" << x.first << "," << x.second << ")" << " <-> " << "(" << y.first << "," << y.second << ")\n";

      big_float ratio = estimate_ratio(k, q, num_simulations, num_gibbs_steps, edge, lattice, neighbours, neighbour_count);
      Z *= ratio;
      // meanRatio += ratio;

      cout << big_int(Z) << "\n";

      neighbours[x.first][x.second][neighbour_count[x.first][x.second]++] = y;
      neighbours[y.first][y.second][neighbour_count[y.first][y.second]++] = x;
  }

  big_float half(0.5);
  big_int rounded_Z = (Z + half).convert_to<big_int>();
  clock_t end = clock();
  cout << "Estimated number of " << q << "-colorings for " << k << "x" << k << " lattice: " << rounded_Z << "\n";

  cerr << "\nfinished in " << clock() * 1.0 / CLOCKS_PER_SEC << " sec\n";
}
// g++ -std=c++20 -O3 main.cpp -I/opt/homebrew/Cellar/boost/1.85.0/include -L/opt/homebrew/Cellar/boost/1.85.0/lib -lboost_system -o main
