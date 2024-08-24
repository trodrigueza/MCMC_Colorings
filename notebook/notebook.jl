### A Pluto.jl notebook ###
# v0.19.39

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ b54280b0-f07d-4adc-a238-65dfc475a5cc
using Random,PlutoUI,Printf

# ╔═╡ 21f2ee62-321d-47fe-a2b0-aa3ff2d558f8
html"""
<h1 style="text-align: center;">Tarea 2</h3>
"""

# ╔═╡ 53933472-805f-4ef9-acd6-d3fe44bb95ea
md"#### Integrantes: 
* David Alejandro Alquichire Rincón  -  *dalquichire@unal.edu.co*
* Kevin Felipe Marroquín Olaya  -  *kfmarroquino@unal.edu.co*
* Tomas David Rodríguez Agudelo  -  *trodrigueza@unal.edu.co*
"

# ╔═╡ 7a35759c-a1c7-47d0-8eec-947dcef2eaa9
html"""
<h2 style="text-align: center;">Primer Punto</h3>
"""

# ╔═╡ 7d1dc16a-b983-4778-a646-9dedabae8f31
md"Con base en el siguiente resultado visto en clase:"

# ╔═╡ 84e4ad14-e289-416c-8471-bf70493d1929
md"""
![Pluto (dwarf planet)](https://i.ibb.co/2k8WPcM/Teorema9-1.png)
"""

# ╔═╡ 7b304895-df81-4b98-9a1d-238bd55532b4
md"Realice experimentos que permitan dar valores aproximados del número de $q$-coloraciones de un lattice $k\times k$."

# ╔═╡ 9fc0084a-2770-4821-929a-3de1a93bb92f
md"Considere $3 \leq k \leq 20$ y $2 \leq q \leq 15$."

# ╔═╡ ee3e5859-4947-49c1-8ed1-edc5987a4582
md"""
**(a) Reporte:**
- Valores $\epsilon$ usados.
- Número de simulaciones (muestras) usadas.
- Números de pasos del "Systematic Gibbs Sampler" usados para cada muestra.
**(b)** Use algún paquete para el conteo exacto y compárelo con algunos de los resultados obtenidos en el item anterior.
"""

# ╔═╡ 50bfb720-cf21-438f-ac1f-0a4134729478
html"""
<h3 style="text-align: center;">Solución.</h3>
"""

# ╔═╡ 72d3565e-6997-48a3-a6c1-2de34fcb27f1
md"""
Antes de comenzar a solucionar los numerales, daremos una breve descripción del algoritmo que implementaremos, el cual se muestra en la demostración del Teorema 9.1 (Libro *Finite Markov chains and algorithm applications*):

Sea $G = (V, E)$ el grafo reticular $k\times k$ con $V = \{v_1, v_2, ..., v_{n=k^k}\}$ y $E = \{e_1, e_2, ..., e_{m=2k(k-1)}\}$. Definimos:

$$G_0 = (V, \emptyset)$$
$$G_i = (V, \{e_1, ..., e_i\})$$
para $1\leq i \leq m$. $Z_i$ será el número de q-coloraciones del grafo $G_i$.

Queremos aproximar $Z_{m}$, el cual puede ser re escrito como:

$$Z_{m} = \frac{Z_{m}}{Z_{m-1}} \times \frac{Z_{m-1}}{Z_{m-2}} \times ... \times \frac{Z_2}{Z_1}\times \frac{Z_1}{Z_0}\times Z_0.$$

1. Comenzaremos con $G_0$, el grafo que tiene todos los vértices (en este caso $k \times k$) pero ninguna arista. Esto significa que cualquier asignación de colores es válida, es decir $Z_0 = q^{k^2}$.

Para cada arista que se añade al grafo, estimaremos la proporción $\frac{Z_i}{Z_{i-1}}$. Para esto, podemos notar que las $q$-coloraciones de $G_i$ son aquellas $q$-coloraciones de $G_{i-1}$ en las que los vértices de la arista añadida en $G_i$ tienen un color distinto. Por lo tanto, la proporción de coloraciones de $G_{i-1}$ que también son válidas en $G_i$ se convierte en $\frac{Z_i}{Z_{i-1}}$. Es decir, $\frac{Z_i}{Z_{i-1}}$ es la probabilidad de que dos vértices conectados por la nueva arista tengan diferentes colores. En este sentido:

2. Para cada arista $e_i = \{x_i, y_i\}$, estimamos $\frac{Z_i}{Z_{i-1}}$ realizando múltiples simulaciones en las que, utilizando Gibbs sampler, obtenemos una $q$-coloración del grafo $G_{i-1}$, así $\frac{Z_i}{Z_{i-1}}$ será la proporción de simulaciónes que resulten en una coloración en la que $x_i, y_i$ tienen colores distintos.

3. Construimos la estimación final multiplicando sucesivamente estas razones, partiendo desde $Z_0$ (que como vimos es trivial de calcular).

Realizamos múltiples simulaciones para cada arista añadida, permitiendo que la distribución de $q$-coloraciones converja adecuadamente realizando suficientes pasos para el Gibbs sampler.
"""

# ╔═╡ 17fe42a9-9a29-496d-b0b3-b4f13330c2cd
html"""
<h5>Implementación.</h5>
"""

# ╔═╡ 93025495-2e64-4efb-9bf5-92a3fcf737de
md"""
Para representar el retículo, el vértice superior izquierdo será etiquetado como $(1, 1)$ y el inferior derecho como $(k, k)$, luego una arista entre los vértices $(x_1, y_1)$ y $(x_2, y_2)$ será representada como la pareja $((x_1, y_1), (x_2, y_2))$.
"""

# ╔═╡ 1fc608d9-85bd-48df-b48c-fed092c15031
md"""
|   |   |   | |
|:---: | :---: |:---:| :---: |
| (1,1) | (1,2) | ... | (1, k) |
| (2,1) | (2,2) | ... | (2, k) |
| . | . | . | . |
| . | . | . | . |
| . | . | . | . |
| (k,1) | (k,2) | ... | (k, k) |
"""

# ╔═╡ 0a94da1c-a84b-4803-9c8c-c49af41899b8
md"""
Bibliotecas utilizadas:
"""

# ╔═╡ 6514bc72-96a7-47b6-b95b-560caef27286
md"""
En primer lugar, implementamos las funciones `create_lattice()`, que inicializa un retículo de $k\times k$ con una coloración aleatoria; y `gen_edges()` que retorna una lista con las aristas del grafo reticular. 
"""

# ╔═╡ 6267be34-ae76-4cc4-b13e-cdf1a353c7f9
function create_lattice(k::Int, q::Int)
    lattice = zeros(Int, k, k) 
    for i in 1:k
        for j in 1:k
            lattice[i, j] = rand(0:q-1)
        end
    end
    return lattice
end

# ╔═╡ 5ebc9437-7992-451f-bbf0-1d324b5b13de
function gen_edges(k::Int)
    edges = []
    for i in 1:k
        for j in 1:k
            if i + 1 <= k
                push!(edges, ((i, j), (i + 1, j))) 
            end
            if j + 1 <= k
                push!(edges, ((i, j), (i, j + 1))) 
            end
        end
    end
    return edges
end

# ╔═╡ c9d8d7ff-9a64-4ad7-b78c-e7d93accb936
md"""
A continuación, implementamos el systematic Gibbs sampler mediante la función `gibbs_step()`. Systematic Gibbs sampler es una variante de Gibbs sampler en la que cada vértice $v_i$ se actualiza sistemáticamente en los tiempos $i, n+i, 2n+i, ...$, donde $n = |V|$. En los tiempos específicos de cada vértice, se actualiza el color de dicho vértice de acuerdo a la distribución uniforme sobre el conjunto de colores disponibles (según los vecinos del vértice). Nosotros iteraremos en orden por cada fila del retículo, luego $v_1 = (1,1), v_{k+1} = (2,1), ..., v_{n} = (k, k)$.
"""

# ╔═╡ d89e9000-44fb-484d-9780-a74a13ce8a07
md"""
|   |   |   | |
|:---: | :---: |:---:| :---: |
| v_1 | v_2 | ... | v_k |
| v_{k+1} | v_{k+2} | ... | v_{2k} |
| . | . | . | . |
| . | . | . | . |
| . | . | . | . |
| v_{(k-1)k+1} | v_{(k-1)k+2} | ... | v_{k^2} |
"""

# ╔═╡ 27ada6b2-521d-44b7-8078-aa64c83139b8
function gibbs_step(k::Int, q::Int, lattice::Array{Int,2}, neighbours)
    used = zeros(Int, q)
    for i in 1:k, j in 1:k
		used = zeros(Int, q)
        unused = []
		if !isempty(neighbours[i,j])
        	for ne in neighbours[i,j]
            	used[lattice[ne[1], ne[2]] + 1] = 1
        	end
		end
        for color in 0:q-1
            if used[color + 1] != 1
                push!(unused, color)
            end
        end
        if !isempty(unused)
            lattice[i, j] = rand(unused)
        end
    end
end

# ╔═╡ 71321119-8d3c-44cc-9b41-e921df940eff
md"""
Note que `gibbs_step()` realiza $k^2$ pasos del systematic Gibbs sampler.
"""

# ╔═╡ eda34072-d90e-4eda-9264-47eadc9b5ceb
md"""
Ahora, implementamos la función `estimate_ratio()`, en la que realizamos `num_simulations` número de simulaciones y -mínimo- `num_gibbs_steps` número de pasos del systematic Gibbs sampler para estimar $\frac{Z_i}{Z_{i-1}}$. En cada simulación, aumentamos la cuenta de `count` si en la coloración obtenida los vértices de `edge` ($e_i = (u, v)$), tienen un color distinto.
"""

# ╔═╡ 68e4fb0e-4d6a-44f8-add0-85ad9bfbc769
function estimate_ratio(k::Int, q::Int, num_simulations::Int, num_gibbs_steps::Int, edge, lattice, neighbours)
    u_y, u_x = edge[1]
    v_y, v_x = edge[2]
    count = 0
    for _ in 1:num_simulations
        for _ in 1:ceil(num_gibbs_steps/(k*k))
            gibbs_step(k, q, lattice, neighbours)
        end
        if (lattice[u_y, u_x] != lattice[v_y, v_x])
			count += 1
		end
    end
    ratio = BigFloat(count / num_simulations)
    return ratio
end

# ╔═╡ 77bf627c-92d8-4905-88a8-9cc24d53954e
md"""
Finalmente, establecemos los parámetros:
"""

# ╔═╡ 5f264e21-ca5b-49bf-994e-22e630984e28
@bind k Slider(1:1:20, default=5)

# ╔═╡ 3ae52f7e-8b13-42a9-810a-ebbdd869d38d
@printf("k = %i", k)

# ╔═╡ 668f9bfc-b97e-437f-8525-b9812f705932
n = k*k

# ╔═╡ 202922d1-1f4a-47c4-bad7-31ddd646594e
@bind q Slider(1:1:15, default=5)

# ╔═╡ 77c1f5c0-a44e-4b82-9490-550c320d5fa9
@printf("q = %i", q)

# ╔═╡ f8e4bf99-1a61-46bf-9836-bf79fe1a8a3f
@bind epsilon Slider(0.1:0.1:1, default=5)

# ╔═╡ dc6111e4-9865-47a3-86d3-135de9b098ce
@printf("epsilon = %.1f", epsilon)

# ╔═╡ 2b91afcf-dca9-4826-8a80-a07c9b766a51
md"""
Inicializamos el retículo:
"""

# ╔═╡ 558c29c4-12ab-4139-81d6-ecc2bb35396f
lattice = create_lattice(k, q)

# ╔═╡ bc10486a-93d8-4944-b70c-0c9ddbee3410
md"""
Generamos las aristas del retículo:
"""

# ╔═╡ 3f6b49eb-3cbb-40cd-accf-404665e67609
edges = gen_edges(k)

# ╔═╡ c668d301-1b14-422d-9465-39c134350ae8
md"""
Declaramos el valor inicial de `Z` como $Z_0=q^n$, para esto utilizamos la función `BigFloat()` de Julia que brinda una precisión arbitraria sobre la variable `Z` (dependiendo de $k$ y de $q$ puede tomar valores muy grandes):
"""

# ╔═╡ 6b19e038-7cb8-42cc-83dd-5ae3d2fee0f9
Z = BigFloat(q)^n

# ╔═╡ d6f66d6c-f144-4734-95d1-293ede2f3303
md"""
Finalmente, establecemos el número de simulaciones y el número de pasos para el systematic Gibbs sampler, para esto utilizamos los resultados del teorema mencionado anteriormente. Como observación, según dichos resultados se deberían hacer $\frac{48d^3n^3}{\epsilon^2}$. Sin embargo realizaremos $\frac{n^3}{\epsilon^2}$ (manteniendo el orden de $n^3$ sobre el número de simulaciones) pues de la otra manera para el ejemplo actual de $k=5$ y $q=10$ habían pasado más de diez mil segundos (2.7 horas) y aún no se obtenía respuesta.
"""

# ╔═╡ 94f2ea4a-2425-4e09-bf37-6418cc4f38d9
num_simulations = Int(ceil(n^3/epsilon^2))

# ╔═╡ b7a2edfe-1368-494a-9d8d-3cba500a0697
num_gibbs_steps = Int(ceil(abs(n * ((2 * log(n) + log(1 / epsilon) + log(8)) / log((q) / 32) + 1))))

# ╔═╡ 408c773c-7368-4473-a473-ae87b928e457
begin
Random.seed!(1234)
neighbours = [Vector{Tuple{Int,Int}}() for _ in 1:k, _ in 1:k]
ratios = []
for edge in edges
	u_y, u_x = edge[1]
    v_y, v_x = edge[2]
	ratio = estimate_ratio(k, q, num_simulations, num_gibbs_steps, edge, lattice, neighbours)
    Z *= BigFloat(ratio)
	push!(ratios, BigFloat(ratio))	
	push!(neighbours[u_y, u_x], (v_y, v_x))
    push!(neighbours[v_y, v_x], (u_y, u_x))
end
val_est = round(Z)
@printf("Número estimado de %d-coloraciones para el retículo %dx%d:\n", q, k, k)
@printf("%d\n", val_est)
end

# ╔═╡ d6b4cc79-499b-496c-9cb2-a901c0be78fe
md"""
Para estimar el error del resultado utilizaremos *Mathematica*, que nativamente implementa una función que permite calcular el polinomio cromático $P_G$ de un grafo arbitrario $G$. Un resultado de teoría de grafos dice que $P_G(q)$ dice el número exacto de $q$-coloraciones del grafo $G$. 

**Nota:** En el repositorio de [GitHub](https://github.com/trodrigueza/MCMC_Colorings) en la carpeta */Mathematica*, se encuentra un notebook en el que se realizan los cálculos que se mostrarán en el presente informe, en caso de que un resultado sea obtenido mediante *Mathematica* escribiremos $\star$.
"""

# ╔═╡ 8616ef86-dcd5-4f3d-83e7-f5329b48db4a
md"""
Sea $L_k$ el grafo reticular $k\times k$, entonces

$$P_{L_5}(10) = 151086899096935604867610 \text{ }\star$$
"""

# ╔═╡ cd95cd1f-27a1-4525-893e-378270fc3657
val_real = 151086899096935604867610

# ╔═╡ 183fc6fb-2ef1-4056-a37f-2d606abdab69
@printf("val_est = %d", val_est)

# ╔═╡ 5ed31698-dda6-4050-83a2-dc56acb620bf
md"""
Es decir que la estimación obtenida tiene un error relativo de:
"""

# ╔═╡ 62b3fe4f-6546-4b5b-8187-69573f0e2f45
@printf("%0.2f%%", (abs(val_real-val_est)/val_real)*100)

# ╔═╡ e7b1d579-3300-49cd-8e0d-28f942f99be9
md"""
Confirmando que es una muy buena estimación del número de $10$-coloraciones del grafo reticular $5\times 5$.
"""

# ╔═╡ d4d9a98a-01ac-478c-a787-37c15dc669e6
html"""
<h5>Implementación (C++).</h5>
"""

# ╔═╡ 00f6943d-4bbc-44a5-b020-8ebc4287f2d0
md"""
Como mencionamos anteriormente, si hubiésemos realizado $\frac{48d^3n^3}{\epsilon^2}$ simulaciones, probablemente se habría obtenido la respuesta después de muchas horas. Es por esto que decidimos realizar una implementación optimizada del código anterior en el lenguaje `C++`. Por ejemplo, dejando los mismos parámetros anteriores ($k=5, q = 10, \epsilon = 1$) y realizando $\frac{48d^3n^3}{\epsilon^2}$ simulaciones el programa en `C++` tarda 4463 segundos, lo cual es poco en comparación con las posibles múltiples horas que habría tardado en Julia. Realizando $\frac{n^3}{\epsilon^2}$ a penas se tarda 2.86 segundos versus 19.4 segundos que tardó en Julia.
"""

# ╔═╡ dcd57bc1-dd22-4f6a-82f2-cdc74c5cce49
html"""
<a><img src="https://i.ibb.co/L0wFDHZ/Screenshot-2024-08-24-at-7-36-46-AM.png" alt="L_5 10-colorings" border="0"></a>
"""

# ╔═╡ 0c9552aa-bbf6-40f4-a53e-c55f1b20b526
html"""
<a><img src="https://i.ibb.co/pQ6SWg2/Screenshot-2024-08-24-at-7-42-51-AM.png" alt="L_5 10-colorings" border="0"></a>
"""

# ╔═╡ 3589ae6a-4ec2-4aa0-a87b-08baf508d260
md"""
A continuación mostramos la implementación realizada.

**Nota:** En la carpeta *cpp* del [GitHub](https://github.com/trodrigueza/MCMC_Colorings) también se encuentra la implementación y detalles sobre la misma.
"""

# ╔═╡ 4eec39a2-deb3-49ad-91d4-839a4fe02f33
md"""
```cpp
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
```
"""

# ╔═╡ cc0ef645-3313-4ca8-92a9-7a7f6e3cf53d
html"""
<h5>Reporte de resultados.</h5>
"""

# ╔═╡ a07a0f97-56a7-4577-822a-54734af93229
md"""
Para esta sección utilizaremos la implementación realizada en `C++`.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "0830c8bfd27ddb998b777ca2f373aa88137e65e2"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─21f2ee62-321d-47fe-a2b0-aa3ff2d558f8
# ╟─53933472-805f-4ef9-acd6-d3fe44bb95ea
# ╟─7a35759c-a1c7-47d0-8eec-947dcef2eaa9
# ╟─7d1dc16a-b983-4778-a646-9dedabae8f31
# ╟─84e4ad14-e289-416c-8471-bf70493d1929
# ╟─7b304895-df81-4b98-9a1d-238bd55532b4
# ╟─9fc0084a-2770-4821-929a-3de1a93bb92f
# ╟─ee3e5859-4947-49c1-8ed1-edc5987a4582
# ╟─50bfb720-cf21-438f-ac1f-0a4134729478
# ╟─72d3565e-6997-48a3-a6c1-2de34fcb27f1
# ╟─17fe42a9-9a29-496d-b0b3-b4f13330c2cd
# ╟─93025495-2e64-4efb-9bf5-92a3fcf737de
# ╟─1fc608d9-85bd-48df-b48c-fed092c15031
# ╟─0a94da1c-a84b-4803-9c8c-c49af41899b8
# ╠═b54280b0-f07d-4adc-a238-65dfc475a5cc
# ╟─6514bc72-96a7-47b6-b95b-560caef27286
# ╠═6267be34-ae76-4cc4-b13e-cdf1a353c7f9
# ╠═5ebc9437-7992-451f-bbf0-1d324b5b13de
# ╟─c9d8d7ff-9a64-4ad7-b78c-e7d93accb936
# ╟─d89e9000-44fb-484d-9780-a74a13ce8a07
# ╠═27ada6b2-521d-44b7-8078-aa64c83139b8
# ╟─71321119-8d3c-44cc-9b41-e921df940eff
# ╟─eda34072-d90e-4eda-9264-47eadc9b5ceb
# ╠═68e4fb0e-4d6a-44f8-add0-85ad9bfbc769
# ╟─77bf627c-92d8-4905-88a8-9cc24d53954e
# ╟─3ae52f7e-8b13-42a9-810a-ebbdd869d38d
# ╟─5f264e21-ca5b-49bf-994e-22e630984e28
# ╟─668f9bfc-b97e-437f-8525-b9812f705932
# ╟─77c1f5c0-a44e-4b82-9490-550c320d5fa9
# ╟─202922d1-1f4a-47c4-bad7-31ddd646594e
# ╟─dc6111e4-9865-47a3-86d3-135de9b098ce
# ╟─f8e4bf99-1a61-46bf-9836-bf79fe1a8a3f
# ╟─2b91afcf-dca9-4826-8a80-a07c9b766a51
# ╠═558c29c4-12ab-4139-81d6-ecc2bb35396f
# ╟─bc10486a-93d8-4944-b70c-0c9ddbee3410
# ╠═3f6b49eb-3cbb-40cd-accf-404665e67609
# ╟─c668d301-1b14-422d-9465-39c134350ae8
# ╠═6b19e038-7cb8-42cc-83dd-5ae3d2fee0f9
# ╟─d6f66d6c-f144-4734-95d1-293ede2f3303
# ╠═94f2ea4a-2425-4e09-bf37-6418cc4f38d9
# ╠═b7a2edfe-1368-494a-9d8d-3cba500a0697
# ╠═408c773c-7368-4473-a473-ae87b928e457
# ╟─d6b4cc79-499b-496c-9cb2-a901c0be78fe
# ╟─8616ef86-dcd5-4f3d-83e7-f5329b48db4a
# ╟─cd95cd1f-27a1-4525-893e-378270fc3657
# ╟─183fc6fb-2ef1-4056-a37f-2d606abdab69
# ╟─5ed31698-dda6-4050-83a2-dc56acb620bf
# ╠═62b3fe4f-6546-4b5b-8187-69573f0e2f45
# ╟─e7b1d579-3300-49cd-8e0d-28f942f99be9
# ╟─d4d9a98a-01ac-478c-a787-37c15dc669e6
# ╟─00f6943d-4bbc-44a5-b020-8ebc4287f2d0
# ╟─dcd57bc1-dd22-4f6a-82f2-cdc74c5cce49
# ╟─0c9552aa-bbf6-40f4-a53e-c55f1b20b526
# ╟─3589ae6a-4ec2-4aa0-a87b-08baf508d260
# ╟─4eec39a2-deb3-49ad-91d4-839a4fe02f33
# ╟─cc0ef645-3313-4ca8-92a9-7a7f6e3cf53d
# ╟─a07a0f97-56a7-4577-822a-54734af93229
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
