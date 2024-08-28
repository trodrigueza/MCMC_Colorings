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
using Random,PlutoUI,Printf, Statistics

# ╔═╡ 3a058fac-395e-48da-a045-2f9496714327
using CSV, DataFrames, Plots, StatsPlots

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

# ╔═╡ 1d2dcd3a-2cf8-4d07-8e60-83272ffe418b
html"""
<div style="text-align: center;">
<img src="https://i.ibb.co/2k8WPcM/Teorema9-1.png"
	width="400" class=center>
</div>
"""

# ╔═╡ 7b304895-df81-4b98-9a1d-238bd55532b4
md"Realice experimentos que permitan dar valores aproximados del número de $q$-coloraciones de un lattice $k\times k$."

# ╔═╡ 9fc0084a-2770-4821-929a-3de1a93bb92f
md"Considere $2 \leq k \leq 20$ y $2 \leq q \leq 15$."

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
<h3 style="text-align: center;">Solución:</h3>
"""

# ╔═╡ 609cdee4-0e35-42f8-9b78-44d427d0da12
md"""
Antes de comenzar a solucionar los númerales, una $q$-coloración de un grafo es una asignación de colores a sus vértices, usando exactamente $q$ colores distintos, de tal manera que ningún par de vértices adyacentes (conectados por una arista) comparta el mismo color.
"""

# ╔═╡ 3a69b474-7ea9-494c-a236-3c50413947f3
html"""
<div style="text-align: center;">
<img src="https://i.ibb.co/c16JWFq/3q.png"
	width="100" class=center alt="Ejemplo 3-coloración del retículo 3x3">
</div>
<div style="text-align: center;">
(Ejemplo 3-coloración del grafo reticular 3x3)
</div>
"""

# ╔═╡ 72d3565e-6997-48a3-a6c1-2de34fcb27f1
md"""
Adiciónalmente, daremos una breve descripción del algoritmo que implementaremos, el cual se muestra en la demostración del Teorema 9.1 (Libro *Finite Markov chains and algorithm applications*):

Sea $G = (V, E)$ el grafo reticular $k\times k$ con $V = \{v_1, v_2, ..., v_{n=k^2}\}$ y $E = \{e_1, e_2, ..., e_{m=2k(k-1)}\}$. Definimos:

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
<h4>Implementación:</h4>
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
A continuación, implementamos el systematic Gibbs sampler mediante la función `gibbs_step()`. Systematic Gibbs sampler es una variante de Gibbs sampler en la que el vértice $v_i$ se actualiza sistemáticamente en los tiempos $i, n+i, 2n+i, ...$, donde $n = |V|$. En los tiempos específicos de cada vértice se actualiza el color de dicho vértice de acuerdo a la distribución uniforme sobre el conjunto de colores disponibles (según los vecinos del vértice), a los demás vértices no se les hace cambio alguno. Nosotros iteraremos en orden por cada fila del retículo, luego $v_1 = (1,1), v_{k+1} = (2,1), ..., v_{n} = (k, k)$.
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

# ╔═╡ cf84e4b2-99cb-4e3f-82bc-7fbf12d5bcfb
m = 2*k*(k-1)

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

# ╔═╡ a9ba504b-3ecc-4e80-b063-06242abb23ca
md"""
Algo interesante que podemos notar realizando estos experimentos es que, al menos para $k = 5$ y $q=10$, las razones $\frac{Z_i}{Z_{i-1}}$ en general son muy cercanas independientemente de la arista $e_i$:
"""

# ╔═╡ 6fdef38f-a45a-4cfe-9694-a0ac144f16b9
ratios

# ╔═╡ fdc76661-9c2e-47f0-8716-a3d719e88998
@printf("Desviación estándar de `ratios`: %0.3f", std(ratios))

# ╔═╡ aa182161-2906-4e5d-aa08-3153c19f97e0
md"""
Es decir que, al menos para este caso, podríamos obtener una estimación del resultado simplemente multiplicando $Z_0$ por $[\frac{Z1}{Z0}]^{|E|}$:
"""

# ╔═╡ 57afe92f-88a4-4fdc-83db-5d58cc525e73
approx_z1 = BigInt(round(BigFloat(q)^n * ratios[1]^m))

# ╔═╡ 205399df-28be-4674-a1c2-35bc5db6e765
@printf("Error relativo: %0.3f%%", (abs(approx_z1 - val_real)/val_real)*100)

# ╔═╡ 62e7c451-ba95-496c-8f16-5d7deea3388a
md"""
Utilizando $\frac{Z2}{Z_1}$:
"""

# ╔═╡ 899c7974-3284-4eb8-bf7e-7d5a3028ec04
approx_z2 = BigInt(round(BigFloat(q)^n * ratios[2]^m))

# ╔═╡ 52042204-a585-46b1-9757-694b9f1de30c
@printf("Error relativo: %0.3f%%", (abs(approx_z2 - val_real)/val_real)*100)

# ╔═╡ 0d70f4d5-08d3-43bc-ae1d-eb856ccfb248
md"""
Utilizando el promedio de `ratios`:
"""

# ╔═╡ 1934af33-79b9-430d-8667-4fb9cc99bddf
approx_zmean = BigInt(round(BigFloat(q)^n * mean(ratios)^m))

# ╔═╡ 33cf49f4-8b88-4d6f-a3d0-42a778836013
@printf("Error relativo: %0.3f%%", (abs(approx_zmean - val_real)/val_real)*100)

# ╔═╡ 00aa11ee-8a3d-424e-9ea5-3c8a030c7d0b
md"""
¿Existe algún resultado sobre esto? Por ejemplo para $k=9$ y $q=3$ las aristas horizontales a partir de la segunda fila tienen una razón en general cercana a $0.75$ mientras que para el resto de las aristas a $0.666$.
"""

# ╔═╡ d4d9a98a-01ac-478c-a787-37c15dc669e6
html"""
<h4>Implementación (C++):</h4>
"""

# ╔═╡ 00f6943d-4bbc-44a5-b020-8ebc4287f2d0
md"""
Como mencionamos anteriormente, si hubiésemos realizado $\frac{48d^3n^3}{\epsilon^2}$ simulaciones, probablemente se habría obtenido la respuesta después de muchas horas, tal vez días. Es por esto que decidimos realizar una implementación optimizada del código anterior en el lenguaje `C++`. Por ejemplo, dejando los mismos parámetros anteriores ($k=5, q = 10, \epsilon = 1$) y realizando $\frac{48d^3n^3}{\epsilon^2}$ simulaciones el programa en `C++` tarda 4463 segundos, lo cual es poco en comparación con las posibles múltiples horas que habría tardado en Julia. Realizando $\frac{n^3}{\epsilon^2}$ a penas se tarda 2.86 segundos versus 19.4 segundos que tardó en Julia.
"""

# ╔═╡ 23b80981-5bc0-407d-b9f1-8ebd8c8eda5b
html"""
<div style="text-align: center;">
<img src="https://i.ibb.co/L0wFDHZ/Screenshot-2024-08-24-at-7-36-46-AM.png" alt="L_5 10-colorings" width="400" class=center>
</div>
"""

# ╔═╡ 1d6f351c-4409-4d69-adc8-dce1a961ac82
@printf("Error relativo: %0.3f%%", (abs(val_real-151122468821583104995130)/val_real)*100)

# ╔═╡ 6040e5b2-8428-4808-bba4-586ee7bed8de
html"""
<div style="text-align: center;">
<img src="https://i.ibb.co/pQ6SWg2/Screenshot-2024-08-24-at-7-42-51-AM.png" alt="L_5 10-colorings" width="400" class=center>
</div>
"""

# ╔═╡ 169716b4-0e42-44a2-ac7a-994b8ab517e3
@printf("Error relativo: %0.3f%%", (abs(val_real-151838548342407414490648)/val_real)*100)

# ╔═╡ 3589ae6a-4ec2-4aa0-a87b-08baf508d260
md"""
A continuación mostramos la implementación realizada.

**Nota:** En la carpeta *cpp* del [GitHub](https://github.com/trodrigueza/MCMC_Colorings)e también se encuentra la implementación y detalles sobre la misma.
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
<h3>Reporte de resultados:</h3>
"""

# ╔═╡ a07a0f97-56a7-4577-822a-54734af93229
md"""
Para esta sección utilizamos la implementación realizada en `C++`.

**Convención:**

- |$k$| Dimensión retículo.
- |$q$| Número de colores.
- |$\text{Gibbs}$| Número de pasos del systematic Gibbs sampler.
- |$\text{Est}$| Estimación obtenida.
- |$\text{Res}\star$| Valor real.
- |$\text{R\_err}$| Error relativo (%).
- |$\text{Mean\_r}$| Promedio de las razones obtenidas.
- |$\text{Time}$| Tiempo en segundos que tardó el programa.
"""

# ╔═╡ 6aaa555a-833c-4662-bfd7-249bd1c93a3e
md"""
En primer lugar, consideraremos resultados obtenidos con los siguientes parámetros:

$\text{num\_simulations} = \frac{48d^3n^3}{\epsilon^2}$
$\text{num\_gibbs\_steps} = n\left(\left|\frac{2\log{n}+\log{\epsilon^{-1}}+\log{8}}{\log{\frac{q}{2d^2}}}\right|+1\right)$

**Nota:** Para no extender mucho el documento, solo mostraremos las tablas para el $\epsilon$ más pequeño, el resto de datos se podrá visualizar en las gráficas. Todos los resultados obtenidos se encuentran en la carpeta *results* del [GitHub](https://github.com/trodrigueza/MCMC_Colorings).
"""

# ╔═╡ e47cb888-c040-4ad4-a7f9-2c6c1604b1b5
md"""
$\textbf{k=2}$
"""

# ╔═╡ cd636c43-de57-4eef-a23f-d1eb3eb3f8e0
md"""
$\epsilon = 0.1$

$\text{num\_simulations} = 19660800$
"""

# ╔═╡ 84c216c7-3e28-472b-b50d-7b3cfd307a81
file_path_a = "../results/q-colorings/k=2/eps=0.1.csv"

# ╔═╡ ea51c9cc-ce50-4758-83f2-3240779604f4
dfa = CSV.File(file_path_a, delim=' ', header=true) |> DataFrame

# ╔═╡ 3c70a826-5f3f-45cf-a054-7d947392f694
md"""
$\epsilon = 1$

$\text{num\_simulations} = 196608$
"""

# ╔═╡ 764aa6f1-f3c3-4cbd-a7ba-d0614eee0093
file_path_b = "../results/q-colorings/k=2/eps=1.csv"

# ╔═╡ 65e5e18e-9df3-4dd6-94f4-218fc9b1c694
dfb = CSV.File(file_path_b, delim=' ', header=true) |> DataFrame;

# ╔═╡ c844196d-f3ac-4ce5-b972-f2fd1682fa13
md"""
Pasos del Gibbs para $2\leq q \leq 15$ ($\epsilon=1$), respectivamente:
"""

# ╔═╡ 58473d1b-26e3-40f7-9c37-a0bdf3e07990
dfb.Gibbs

# ╔═╡ 7c777c46-1b5f-4257-9706-f181a08c7e5a
md"""
###### Visualización de los datos:
"""

# ╔═╡ e45ce12c-2a75-49b0-a2e6-b96a8467cd58
begin
pb1 = plot(dfa.q, dfa.Res⭑, label="Respuesta Real", linewidth=0.5, marker=:star, color=:red, markersize=5, size=(550, 350), yscale=:log10, leg=:bottomright, ylimits = (1e0, 1e5), titlefontsize=10, labelfontsize=8, legendfontsize=6)
	
plot!(dfa.q, dfa.Est, label="Estimación eps = 0.1", title="Estimación vs Respuesta Real L_2", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:rect, color=:blue, alpha=0.5, markersize=3.5)
	
plot!(dfb.q, dfb.Est, label="Estimación eps = 1", linewidth=0.7, marker=:diamond, color=:orange, alpha=0.5, markersize=3)

pb2 = plot(dfa.q, dfa.R_err, label="eps = 0.1", color=:blue, xlabel="q", ylabel="Error relativo (%)", marker=:rect, size=(550, 350), title="Comparación errores relativos L_2", titlefontsize=10, labelfontsize=8, legendfontsize=6)

plot!(dfb.q, dfb.R_err, label="eps = 1", marker=:diamond, color=:orange)

pb3 = bar(dfa.q, dfa.Time, label="eps = 0.1", color=:blue, xlabel="q", ylabel="Tiempo (Segundos)", size=(550, 350), bar_width=0.5, title="Comparación tiempos L_2", titlefontsize=10, labelfontsize=8, legendfontsize=6)

bar!(dfb.q, dfb.Time, label="eps = 1", color=:orange, ylabel="Tiempo(segs)", size=(550, 350), bar_width=0.6)
	
layout = @layout [a b; c]
plot(pb1, pb3, pb2, layout=layout, size=(600,400))
end

# ╔═╡ 7b14e5d2-2774-4876-9458-e631a1615447
md"""
Podemos observar como las estimaciones con $\epsilon = 1$ se acercan de muy buena manera a la respuesta teniendo en cuenta que en general las estimaciones con $\epsilon=0.1$ tardan en promedio 175 veces más.
"""

# ╔═╡ aaf90bb9-c2b7-4e5a-980f-c9580f565ccd
md"""
**Nota:** La gráfica "Estimación vs Respuesta" está en escala logarítmica ($\log_{10}$).
"""

# ╔═╡ 4a8d77dc-f16d-4c93-94a4-c2e1f3f6e1ef
md"""

"""

# ╔═╡ ca36d1eb-86bc-473e-9bbe-dac1cfc13fa6
md"""
$\textbf{k=3}$
"""

# ╔═╡ b5f024b3-aa26-479e-b4c0-4535696b8f52
md"""
$\epsilon=0.31$

$\text{num\_simulations} = 22394879$
"""

# ╔═╡ 5184161f-547b-4cfb-94e6-f55b9765fd76
file_path_c = "../results/q-colorings/k=3/eps=0.31.csv"

# ╔═╡ 33b74cac-6dcf-4211-882e-0dc6abd4c435
dfc = CSV.File(file_path_c, delim=' ', header=true) |> DataFrame

# ╔═╡ 2a924600-584c-41e0-b344-43c1e34411a8
md"""
**Nota:** Estos cálculos se hicieron con una versión no optimizada del código, de ahí los elevados tiempos de ejecución.
"""

# ╔═╡ ada44c79-3b92-4be9-8b0f-14c09c7547a0
md"""
$\epsilon=1$

$\text{num\_simulations} = 2239488$
"""

# ╔═╡ c79341d1-13d3-4e71-9297-a8df019a63d9
file_path_d = "../results/q-colorings/k=3/eps=1.csv"

# ╔═╡ f771ff8f-6e39-4fe7-9f9d-a7c9b4191fc5
dfd = CSV.File(file_path_d, delim=' ', header=true) |> DataFrame;

# ╔═╡ 3c858f7b-51b7-43c0-9c34-1c9294973825
md"""
Pasos del Gibbs para $2\leq q \leq 15$ ($\epsilon=1$), respectivamente:
"""

# ╔═╡ 9675b48f-5468-4030-8117-d090d2074d6b
dfd.Gibbs

# ╔═╡ f307b7c6-8350-4d5e-bafe-b893b9de5407
md"""
###### Visualización de los datos:
"""

# ╔═╡ fd195252-d166-44bc-95a8-2280201c1bb8
begin
pa1 = plot(dfc.q, dfc.Res⭑, label="Respuesta Real", linewidth=0.5, marker=:star, color=:red, markersize=5, size=(550, 350), yscale=:log10, leg=:bottomright, ylimits = (1e0, 1e12), titlefontsize=10, labelfontsize=8, legendfontsize=6)
	
plot!(dfc.q, dfc.Est, label="Estimación eps = 0.31", title="Estimación vs Respuesta Real L_3", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:rect, color=:blue, alpha=0.5, markersize=3.5)
	
plot!(dfd.q, dfd.Est, label="Estimación eps = 1", linewidth=0.7, marker=:diamond, color=:orange, alpha=0.5, markersize=3)

pa2 = plot(dfc.q, dfc.R_err, label="eps = 0.31", color=:blue, xlabel="q", ylabel="Error relativo (%)", marker=:rect, size=(550, 350), title="Comparación errores relativos L_3", titlefontsize=10, labelfontsize=8, legendfontsize=6, leg=:top)

plot!(dfd.q, dfd.R_err, label="eps = 1", marker=:diamond, color=:orange)

pa3 = bar(dfc.q, dfc.Time, label="eps = 0.31", xlabel="q", color=:blue, size=(550, 350), bar_width=0.5, title="Comparación tiempos L_3", titlefontsize=10, labelfontsize=8, legendfontsize=6)

bar!(dfd.q, dfd.Time, label="eps = 1", color=:orange, ylabel="Tiempo(segs)", size=(550, 350), bar_width=0.6)
	
plot(pa1, pa3, pa2, layout=layout, size=(600, 400))
end

# ╔═╡ 9c29cb2c-d2c3-4815-ad8e-c53a7d34722d
md"""
**Nota:** Los cálculos para $\epsilon=0.31$ se hicieron con una versión no optimizada del código, de ahí los tiempos de ejecución.
"""

# ╔═╡ aa14e81b-0530-4579-8e97-c62d49046c60
md"""

Ahora, consideraremos resultados obtenidos con los siguientes parámetros:

$\text{num\_simulations} = \frac{n^3}{\epsilon^2}$
$\text{num\_gibbs\_steps} = n\left(\left|\frac{2\log{n}+\log{\epsilon^{-1}}+\log{8}}{\log{\frac{q}{2d^2}}}\right|+1\right)$
"""

# ╔═╡ e14facff-e9d5-49f5-84ae-2cd8626d5986
md"""
El orden de complejidad total sigue siendo de $Cn^5\log{n}$ sin embargo hemos reducido la constante $C.$
"""

# ╔═╡ 4465b037-9586-4439-af42-3755aad96e31
md"""

"""

# ╔═╡ 5ac758b9-5cda-40c5-b829-4603b00aef03
md"""
$\textbf{k=4}$
"""

# ╔═╡ ce6c8449-c44d-4309-b5e1-00b4adf35950
md"""
$\epsilon=0.1$

$\text{num\_simulations} = 409600$
"""

# ╔═╡ 971394c0-4918-4e37-9324-54333cd3816e
file_path_e = "../results/q-colorings/k=4/eps=0.1.csv"

# ╔═╡ 38894819-bfdc-4ee5-a3c8-237f6dbd400e
dfe = CSV.File(file_path_e, delim=' ', header=true) |> DataFrame

# ╔═╡ 03cb0dac-f17d-46e4-9578-8a72df036a5e
md"""
$\epsilon=0.5$

$\text{num\_simulations} = 16384$
"""

# ╔═╡ 4ea0cb09-d6db-4e0c-b844-9ba06ff54048
file_path_f = "../results/q-colorings/k=4/eps=0.5.csv"

# ╔═╡ 5acace06-a72e-4b8e-a019-eb75f1b80fe9
dff = CSV.File(file_path_f, delim=' ', header=true) |> DataFrame;

# ╔═╡ f5d128d4-8f46-4196-803b-0de3b90b7bdc
md"""
Pasos del Gibbs para $2\leq q \leq 15$ ($\epsilon=0.5$), respectivamente:
"""

# ╔═╡ fd746cac-f3a7-4828-8361-b1d4194e0f5d
dff.Gibbs

# ╔═╡ b5395b0e-1c9b-489a-b903-642838f15d97
md"""
$\epsilon=1$

$\text{num\_simulations} = 4096$
"""

# ╔═╡ e56ee173-8e24-4f88-ab3b-6d2dfe3f51d1
file_path_g = "../results/q-colorings/k=4/eps=1.csv"

# ╔═╡ eeb80031-9afc-4aa2-9b86-66199dcf9746
dfg = CSV.File(file_path_g, delim=' ', header=true) |> DataFrame;

# ╔═╡ e28113c4-6e57-486b-b58a-9f380db8485a
md"""
Pasos del Gibbs para $2\leq q \leq 15$ ($\epsilon=1$), respectivamente:
"""

# ╔═╡ 0a81a4a9-7d99-4ea4-b698-f3d907f12b20
dfg.Gibbs

# ╔═╡ 23a04d74-d7be-4b7a-bb11-d9a64a469685
md"""
###### Visualización de los datos:
"""

# ╔═╡ 5d04d261-896b-4fad-8c75-76af8037c9f0
begin
pc1 = plot(dfe.q, dfe.Res⭑, label="Respuesta Real", linewidth=0.5, marker=:star, color=:red, markersize=5, size=(550, 350), yscale=:log10, leg=:bottomright, ylimits = (1e0, 1e20), titlefontsize=10, labelfontsize=8, legendfontsize=6)
	
plot!(dfe.q, dfe.Est, label="Estimación eps = 0.1", title="Estimación vs Respuesta Real L_4", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:rect, color=:blue, alpha=0.5, markersize=3.5)
	
plot!(dff.q, dff.Est, label="Estimación eps = 0.5", linewidth=0.7, marker=:diamond, color=:green, alpha=0.5, markersize=3)

plot!(dfe.q, dfg.Est, label="Estimación eps = 1", title="Estimación vs Respuesta Real L_4", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:circle, color=:orange, alpha=0.5, markersize=3.5)

pc2 = plot(dfe.q, dfe.R_err, label="eps = 0.1", color=:blue, xlabel="q", ylabel="Error relativo (%)", marker=:rect, size=(550, 350), title="Comparación errores relativos L_4", titlefontsize=10, labelfontsize=8, legendfontsize=6)

plot!(dff.q, dff.R_err, label="eps = 0.5", marker=:diamond, color=:green)
	
plot!(dfe.q, dfg.R_err, label="eps = 1", color=:orange, marker=:circle, secondary=true)

pc3 = bar(dfe.q, dfe.Time, label="eps = 0.1", color=:blue, xlabel="q", ylabel="Tiempo (Segundos)", size=(550, 350), bar_width=0.5, title="Comparación tiempos L_4", titlefontsize=10, labelfontsize=8, legendfontsize=6)

bar!(dfe.q, dff.Time, label="eps = 0.5", color=:green, ylabel="Tiempo(seg)", size=(550, 350), bar_width=0.6)
	
bar!(dfe.q, dfg.Time, label="eps = 1", color=:orange, secondary=true, bar_width=0.7)
	
plot(pc1, pc3, pc2, layout=layout, size=(600,400))
end

# ╔═╡ 395d9d7d-93c0-4d35-8bec-a44b030e93a4
md"""
**Observación:** Debido a que redujimos la constante en la complejidad, el error relativo en epsilons grandes aumentó considerablemente.
"""

# ╔═╡ 3888f02c-a7e8-4959-a209-fbed5376ae67
md"""

"""

# ╔═╡ 280be26f-635f-4382-8870-fdf0e09223aa
md"""
$\textbf{k=5}$
"""

# ╔═╡ eaf519b1-2c36-461f-89a8-ad2ba52ab3c4
md"""
$\epsilon=0.1$

$\text{num\_simulations} = 1562500$
"""

# ╔═╡ 55cf5973-1db4-492d-b1a7-e237119c69e0
file_path_h = "../results/q-colorings/k=5/eps=0.1.csv"

# ╔═╡ dd5f08d6-d743-4d41-b48e-77e70c901981
dfh = CSV.File(file_path_h, delim=' ', header=true) |> DataFrame

# ╔═╡ e16f67a2-ef40-4bbf-9b28-23a0ca3c5ff5
md"""
$\epsilon=0.4$

$\text{num\_simulations} = 97656$
"""

# ╔═╡ 816d2f18-60b7-42b9-9e52-c1f2d9151db4
file_path_i = "../results/q-colorings/k=5/eps=0.4.csv"

# ╔═╡ d35643f4-7bbc-4f92-adc2-49f1df506be3
dfi = CSV.File(file_path_i, delim=' ', header=true) |> DataFrame;

# ╔═╡ 40b195f2-31d8-4e6e-a1f7-367df2ede70f
md"""
Pasos del Gibbs para $2\leq q \leq 15$ ($\epsilon=0.4$), respectivamente:
"""

# ╔═╡ cc6826d1-82b9-4e8d-af67-9aa34e6b76ae
dfi.Gibbs

# ╔═╡ 79913af6-060a-4823-843a-d591139df7e4
md"""
$\epsilon=0.8$

$\text{num\_simulations} = 24414$
"""

# ╔═╡ e16423d5-3456-4f53-acbd-4ede3c9b7ca4
file_path_j = "../results/q-colorings/k=5/eps=0.8.csv"

# ╔═╡ 5fb677e9-8353-45fe-a0ac-467b0bb39f37
dfj = CSV.File(file_path_j, delim=' ', header=true) |> DataFrame;

# ╔═╡ c94b11d5-3c2b-4aa4-932e-97c1a5dade17
md"""
Pasos del Gibbs para $2\leq q \leq 15$ ($\epsilon=0.8$), respectivamente:
"""

# ╔═╡ 8d89bf10-0f60-42d2-bcbf-6d020c22d2ca
dfj.Gibbs

# ╔═╡ 844bfa80-bbd8-4ca9-929a-185c4b9cc536
md"""
###### Visualización de los datos:
"""

# ╔═╡ 912ae862-d93b-4f6b-bbff-106e106cb4e5
begin
pd1 = plot(dfh.q, dfh.Res⭑, label="Respuesta Real", linewidth=0.5, marker=:star, color=:red, markersize=5, size=(550, 350), yscale=:log10, leg=:bottomright, ylimits = (1e0, 1e30), titlefontsize=10, labelfontsize=8, legendfontsize=6)
	
plot!(dfh.q, dfh.Est, label="Estimación eps = 0.1", title="Estimación vs Respuesta Real L_5", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:rect, color=:blue, alpha=0.5, markersize=3.5)
	
plot!(dfi.q, dfi.Est, label="Estimación eps = 0.4", linewidth=0.7, marker=:diamond, color=:green, alpha=0.5, markersize=3)

plot!(dfj.q, dfj.Est, label="Estimación eps = 0.8", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:circle, color=:orange, alpha=0.5, markersize=3.5)

pd2 = plot(dfh.q, dfh.R_err, label="eps = 0.1", color=:blue, xlabel="q", ylabel="Error relativo (%)", marker=:rect, size=(550, 350), title="Comparación errores relativos L_5", titlefontsize=10, labelfontsize=8, legendfontsize=6)

plot!(dfi.q, dfi.R_err, label="eps = 0.4", marker=:diamond, color=:green)
	
plot!(dfj.q, dfj.R_err, label="eps = 0.8", color=:orange, marker=:circle, secondary=true)

pd3 = bar(dfh.q, dfh.Time, label="eps = 0.1", color=:blue, xlabel="q", ylabel="Tiempo (Segundos)", size=(550, 350), bar_width=0.5, title="Comparación tiempos L_5", titlefontsize=10, labelfontsize=8, legendfontsize=6)

bar!(dfi.q, dfi.Time, label="eps = 0.4", color=:green, ylabel="Tiempo(segs)", size=(550, 350), bar_width=0.6)
	
bar!(dfj.q, dfj.Time, label="eps = 0.8", color=:orange, secondary=true, bar_width=0.7)
	
plot(pd1, pd3, pd2, layout=layout, size=(600,400))
end

# ╔═╡ 78c00185-8dee-434b-b6a3-7573e14a6cf4
md"""
**Observación:** Aunque el teorema pide que $q>2d=8$, hasta el momento hemos conseguido muy buenos resultados para $q \leq 8$.
"""

# ╔═╡ 0d364111-35f0-4a3f-9c7a-1fcc426102bb
md"""

"""

# ╔═╡ d3eff505-597f-478d-8e01-5c2f6f365f19
md"""
$\textbf{k=6}$
"""

# ╔═╡ 2c6110b5-cc07-4bf8-9515-7d7585e7b25a
md"""
$\epsilon=0.2$

$\text{num\_simulations} = 1166400$
"""

# ╔═╡ 2703165f-1c4a-4ae0-8723-4fa00b54560b
file_path_k = "../results/q-colorings/k=6/eps=0.2.csv"

# ╔═╡ 0f282d75-f954-4611-ab7e-d62d7d588e10
dfk = CSV.File(file_path_k, delim=' ', header=true) |> DataFrame

# ╔═╡ 5b9def3d-79e3-4c11-8f9f-27ed9f9a2e82
md"""
$\epsilon=0.5$

$\text{num\_simulations} = 186624$
"""

# ╔═╡ b0165c9e-2828-418d-8265-298803cd721f
file_path_l = "../results/q-colorings/k=6/eps=0.5.csv"

# ╔═╡ 462fc69d-e11e-4f54-8182-ad42ff6321ef
dfl = CSV.File(file_path_l, delim=' ', header=true) |> DataFrame;

# ╔═╡ 74de5a7c-aa7c-44ee-b752-82f232f5c00a
md"""
Pasos del Gibbs para $2\leq q \leq 15$ ($\epsilon=0.5$), respectivamente:
"""

# ╔═╡ b5fc3fb7-a5ce-4107-9cf2-4b43999b6c9a
dfl.Gibbs

# ╔═╡ 21b1c011-8039-4347-9c80-101ad2fdb07a
md"""
$\epsilon=0.8$

$\text{num\_simulations} = 72900$
"""

# ╔═╡ 1dab9389-aec3-49d6-bcc2-57bbe367f84b
file_path_m = "../results/q-colorings/k=6/eps=0.8.csv"

# ╔═╡ a90fbfa5-49f4-4308-b977-bb271d1c76d9
dfm = CSV.File(file_path_m, delim=' ', header=true) |> DataFrame;

# ╔═╡ 29e3b0ff-5a4a-418c-95c0-6d9861abbb31
md"""
Pasos del Gibbs para $2\leq q \leq 15$ ($\epsilon=0.8$), respectivamente:
"""

# ╔═╡ 49146722-3a28-4610-a625-e7b06d013e87
dfm.Gibbs

# ╔═╡ d691c2ba-6759-4b53-901d-165ff14fad01
md"""
###### Visualización de los datos:
"""

# ╔═╡ 8a394b92-8c1b-4653-b2c4-1722a257ec60
begin
pe1 = plot(dfk.q, dfk.Res⭑, label="Respuesta Real", linewidth=0.5, marker=:star, color=:red, markersize=5, size=(550, 350), yscale=:log10, leg=:bottomright, ylimits = (1e0, 1e42), titlefontsize=10, labelfontsize=8, legendfontsize=6)
	
plot!(dfk.q, dfk.Est, label="Estimación eps = 0.2", title="Estimación vs Respuesta Real L_6", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:rect, color=:blue, alpha=0.5, markersize=3.5)
	
plot!(dfl.q, dfl.Est, label="Estimación eps = 0.5", linewidth=0.7, marker=:diamond, color=:green, alpha=0.5, markersize=3)

plot!(dfm.q, dfm.Est, label="Estimación eps = 0.8", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:circle, color=:orange, alpha=0.5, markersize=3.5)

pe2 = plot(dfk.q, dfk.R_err, label="eps = 0.2", color=:blue, xlabel="q", ylabel="Error relativo (%)", marker=:rect, size=(550, 350), title="Comparación errores relativos L_6", titlefontsize=10, labelfontsize=8, legendfontsize=6)

plot!(dfl.q, dfl.R_err, label="eps = 0.5", marker=:diamond, color=:green)
	
plot!(dfm.q, dfm.R_err, label="eps = 0.8", color=:orange, marker=:circle, secondary=true)

pe3 = bar(dfk.q, dfk.Time, label="eps = 0.2", color=:blue, xlabel="q", ylabel="Tiempo (Segundos)", size=(550, 350), bar_width=0.5, title="Comparación tiempos L_6", titlefontsize=10, labelfontsize=8, legendfontsize=6)

bar!(dfl.q, dfl.Time, label="eps = 0.5", color=:green, ylabel="Tiempo(segs)", size=(550, 350), bar_width=0.6)
	
bar!(dfm.q, dfm.Time, label="eps = 0.8", color=:orange, secondary=true, bar_width=0.7)
	
plot(pe1, pe3, pe2, layout=layout, size=(600,350))
end

# ╔═╡ 188a3b33-3bb7-486e-a232-06bd50bf778d
md"""
$\textbf{k=7}$
"""

# ╔═╡ 5c2dabb8-2d77-4afa-bd0f-233f636d863e
md"""
$\epsilon=0.4$

$\text{num\_simulations} = 735306$
"""

# ╔═╡ 4cab9347-40b6-421b-8ee3-c5306d60fffe
file_path_n = "../results/q-colorings/k=7/eps=0.4.csv"

# ╔═╡ dee1520c-80eb-4611-9cb0-ad93cb32f9aa
dfn = CSV.File(file_path_n, delim=' ', header=true) |> DataFrame

# ╔═╡ 794b9a20-951f-4ea1-8c74-018cefb47ae4
md"""
$\epsilon=0.6$

$\text{num\_simulations} = 326802$
"""

# ╔═╡ 79ae96ec-ede3-4ca8-ba1c-2139e1d4b2c5
file_path_o = "../results/q-colorings/k=7/eps=0.6.csv"

# ╔═╡ 01121bc2-60c6-4cba-9c68-e99def46736c
dfo = CSV.File(file_path_o, delim=' ', header=true) |> DataFrame;

# ╔═╡ ea157415-146c-4755-90b4-7060cdff3302
md"""
Pasos del Gibbs para $2\leq q \leq 15$ ($\epsilon=0.6$), respectivamente:
"""

# ╔═╡ 599bcd25-e939-4e3f-9d03-cdf27badeef9
dfo.Gibbs

# ╔═╡ 4804a758-7548-42a6-9c21-5e82668d82d8
md"""
$\epsilon=0.8$

$\text{num\_simulations} = 183826$
"""

# ╔═╡ 3d85d992-dec5-4984-ab61-2fc71fe592c9
file_path_p = "../results/q-colorings/k=7/eps=0.8.csv"

# ╔═╡ 7c248a7c-22ae-4752-9ae9-6d66b42dc1da
dfp = CSV.File(file_path_p, delim=' ', header=true) |> DataFrame;

# ╔═╡ b699f6a1-f065-4a2c-b20d-bdf571102e98
md"""
Pasos del Gibbs para $2\leq q \leq 15$ ($\epsilon=0.8$), respectivamente:
"""

# ╔═╡ 0e91e96c-5ffd-4972-a1dd-c96524da29a4
dfp.Gibbs

# ╔═╡ b867bac9-3984-49c4-9afa-de729db21f7b
md"""
###### Visualización de los datos:
"""

# ╔═╡ ff80f7b7-d5e6-473a-8626-190d9f3eedca
begin
pf1 = plot(dfn.q, dfn.Res⭑, label="Respuesta Real", linewidth=0.5, marker=:star, color=:red, markersize=5, size=(550, 350), yscale=:log10, leg=:bottomright, ylimits = (1e0, 1e60), titlefontsize=10, labelfontsize=8, legendfontsize=6)
	
plot!(dfn.q, dfn.Est, label="Estimación eps = 0.4", title="Estimación vs Respuesta Real L_7", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:rect, color=:blue, alpha=0.5, markersize=3.5)
	
plot!(dfo.q, dfo.Est, label="Estimación eps = 0.6", linewidth=0.7, marker=:diamond, color=:green, alpha=0.5, markersize=3)

plot!(dfp.q, dfp.Est, label="Estimación eps = 0.8", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:circle, color=:orange, alpha=0.5, markersize=3.5)

pf2 = plot(dfn.q, dfn.R_err, label="eps = 0.4", color=:blue, xlabel="q", ylabel="Error relativo (%)", marker=:rect, size=(550, 350), title="Comparación errores relativos L_7", titlefontsize=10, labelfontsize=8, legendfontsize=6)

plot!(dfo.q, dfo.R_err, label="eps = 0.6", marker=:diamond, color=:green)
	
plot!(dfp.q, dfp.R_err, label="eps = 0.8", color=:orange, marker=:circle, secondary=true)

pf3 = bar(dfn.q, dfn.Time, label="eps = 0.4", color=:blue, xlabel="q", ylabel="Tiempo (Segundos)", size=(550, 350), bar_width=0.5, title="Comparación tiempos L_7", titlefontsize=10, labelfontsize=8, legendfontsize=6)

bar!(dfo.q, dfo.Time, label="eps = 0.6", color=:green, ylabel="Tiempo(segs)", size=(550, 350), bar_width=0.6)
	
bar!(dfp.q, dfp.Time, label="eps = 0.8", color=:orange, secondary=true, bar_width=0.7)
	
plot(pf1, pf3, pf2, layout=layout, size=(600,400))
end

# ╔═╡ 7515f748-d9ad-4444-afba-43dcabcc165c
md"""

"""

# ╔═╡ 4b6ea471-822d-4333-af3b-600a5cdac961
md"""

En este punto tenemos que decir que para $k\geq8$ *Mathematica* ya no permite calcular el polinomio cromático debido al límite computacional sobre el plan de la prueba gratuita. Intentamos algunos recursos como calcular directamente el polinomio cromático con la librería `NetworkX` de `Python` o calcular el [polinomio de Tutte](https://arminda.whitman.edu/_flysystem/fedora/2021-10/The_Tutte_polynomial_and_applications.pdf) con un [programa](https://github.com/thorehusfeldt/tutte_bhkk) que encontramos en lenguaje `C`, sin embargo al parecer $L_{k\geq 8}$ ya es un grafo lo suficientemente grande como para considerar estas opciones viables. Es por esto que para los siguientes resultados compararemos las estimaciones obtenidas con las respuestas que se muestran en [esta página](https://oeis.org/wiki/Colorings_of_grid_graphs), que contiene resultados hasta $k=13$ y $q=10$.

Ahora, consideraremos resultados obtenidos con los siguientes parámetros:

$\text{num\_simulations} = \frac{n^3}{\epsilon^2}$
$\text{num\_gibbs\_steps} = n\left(\frac{\log{n} + \log{\epsilon^{-1}}}{\log{\frac{q}{4d^2}}}\right)$
"""

# ╔═╡ aee25458-9214-4b0f-ab00-ae6060efb2a5
md"""
De nuevo, el orden de complejidad total sigue siendo de $Cn^5\log{n}$ para alguna constante $C$.
"""

# ╔═╡ 6fad3bd0-76e9-48c3-87fc-4139d083e65d
md"""
$\textbf{k=8}$
"""

# ╔═╡ 75eb8374-d24b-4ab3-9412-6a9255b7090e
md"""
$\epsilon=0.4$

$\text{num\_simulations} = 1638400$
"""

# ╔═╡ 2d21efe6-773d-4255-b04b-be20dd609942
file_path_q = "../results/q-colorings/k=8/eps=0.4.csv"

# ╔═╡ 02d28a01-4ade-47c5-9381-0d9806f5399d
dfq = CSV.File(file_path_q, delim=' ', header=true) |> DataFrame

# ╔═╡ a3152eb3-9b54-4d1e-ae39-4c51a922af59
md"""
$\epsilon=0.6$

$\text{num\_simulations} = 728177$
"""

# ╔═╡ 54770db9-45ed-46d9-ab07-c78aff4d3b23
file_path_r = "../results/q-colorings/k=8/eps=0.6.csv"

# ╔═╡ 2db53eb3-b9c4-4fc1-8b87-2189b3c835a3
dfr = CSV.File(file_path_r, delim=' ', header=true) |> DataFrame;

# ╔═╡ 8f4391d4-8f1e-4e30-9229-196a30b4f0d7
md"""
Pasos del Gibbs para $2\leq q \leq 10$ ($\epsilon=0.6$), respectivamente:
"""

# ╔═╡ 391c4da7-82b8-4e92-81db-9e56bd32f8a5
dfr.Gibbs

# ╔═╡ d59a90ca-f666-4cfd-b8c5-ea81e56d22d3
md"""
$\epsilon=0.8$

$\text{num\_simulations} = 409600$
"""

# ╔═╡ 319cc440-fa4c-4d09-bfa3-aec77d6a79bf
file_path_s = "../results/q-colorings/k=8/eps=0.8.csv"

# ╔═╡ e0de149f-0e42-41fc-9173-53f85eae9acd
dfs = CSV.File(file_path_s, delim=' ', header=true) |> DataFrame;

# ╔═╡ b815249b-b1ee-42fb-a7db-144ae967eab8
md"""
Pasos del Gibbs para $2\leq q \leq 10$ ($\epsilon=0.8$), respectivamente:
"""

# ╔═╡ 125190f1-0731-4757-95af-0a64b891b294
dfs.Gibbs

# ╔═╡ 9a6ce4d5-0095-4a42-a18c-059360750f53
md"""
###### Visualización de los datos:
"""

# ╔═╡ d254202d-162c-43ff-96c6-c6774e843f5c
begin
pg1 = plot(dfq.q, dfq.Res⭑, label="Respuesta Real", linewidth=0.5, marker=:star, color=:red, markersize=5, size=(550, 350), yscale=:log10, leg=:bottomright, ylimits = (1e0, 1e60), titlefontsize=10, labelfontsize=8, legendfontsize=6)
	
plot!(dfq.q, dfq.Est, label="Estimación eps = 0.4", title="Estimación vs Respuesta Real L_8", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:rect, color=:blue, alpha=0.5, markersize=3.5)
	
plot!(dfr.q, dfr.Est, label="Estimación eps = 0.6", linewidth=0.7, marker=:diamond, color=:green, alpha=0.5, markersize=3)

plot!(dfs.q, dfs.Est, label="Estimación eps = 0.8", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:circle, color=:orange, alpha=0.5, markersize=3.5)

pg2 = plot(dfq.q, dfq.R_err, label="eps = 0.4", color=:blue, xlabel="q", ylabel="Error relativo (%)", marker=:rect, size=(550, 350), title="Comparación errores relativos L_8", titlefontsize=10, labelfontsize=8, legendfontsize=6)

plot!(dfr.q, dfr.R_err, label="eps = 0.6", marker=:diamond, color=:green)
	
plot!(dfs.q, dfs.R_err, label="eps = 0.8", color=:orange, marker=:circle, secondary=true)

pg3 = bar(dfq.q, dfq.Time, label="eps = 0.4", color=:blue, xlabel="q", ylabel="Tiempo (Segundos)", size=(550, 350), bar_width=0.5, title="Comparación tiempos L_8", titlefontsize=10, labelfontsize=8, legendfontsize=6)

bar!(dfr.q, dfr.Time, label="eps = 0.6", color=:green, ylabel="Tiempo(segs)", size=(550, 350), bar_width=0.6)
	
bar!(dfs.q, dfs.Time, label="eps = 0.8", color=:orange, secondary=true, bar_width=0.7)
	
plot(pg1, pg3, pg2, layout=layout, size=(600,400))
end

# ╔═╡ aad8f02d-3843-4187-b007-381e535b4c42
md"""

"""

# ╔═╡ 305a180b-e65b-44ed-8754-00e3b32a139a
md"""
$\textbf{k=9}$
"""

# ╔═╡ 8737cb3e-f390-4dbd-ab53-1c1671ea36e4
md"""
$\epsilon=0.5$

$\text{num\_simulations} = 2125764$
"""

# ╔═╡ 745aa527-46cb-4777-b996-0e4b54633afd
file_path_t = "../results/q-colorings/k=9/eps=0.5.csv"

# ╔═╡ b2ed57d7-a199-48a1-9951-1d0cdef6d076
dft = CSV.File(file_path_t, delim=' ', header=true) |> DataFrame

# ╔═╡ a24dd754-63b8-4fcc-8a76-38e085c8dc96
md"""
$\epsilon=0.7$

$\text{num\_simulations} = 1084573$
"""

# ╔═╡ f243e139-380e-49cb-955d-7a580c8a7b33
file_path_u = "../results/q-colorings/k=9/eps=0.7.csv"

# ╔═╡ 370fa37f-368c-49e1-9b3f-65f76993b09c
dfu = CSV.File(file_path_u, delim=' ', header=true) |> DataFrame;

# ╔═╡ d099d0df-ae74-4ceb-af53-305e6029e5d2
md"""
Pasos del Gibbs para $2\leq q \leq 10$ ($\epsilon=0.7$), respectivamente:
"""

# ╔═╡ d7f72084-8fcb-4d1f-994b-b0e31ecaccff
dfu.Gibbs

# ╔═╡ 3747cc76-bbe2-48a9-bc38-69529af96100
md"""
$$\epsilon=0.9$$

$\text{num\_simulations} = 656099$
"""

# ╔═╡ 93e9755a-5919-4d01-8f69-e2f363088c2d
file_path_v = "../results/q-colorings/k=9/eps=0.9.csv"

# ╔═╡ 994c886e-e03a-47f5-b04b-647f446567b3
dfv = CSV.File(file_path_v, delim=' ', header=true) |> DataFrame;

# ╔═╡ 4cfec33d-3abd-42bf-927a-0059c6767e14
md"""
Pasos del Gibbs para $2\leq q \leq 10$ ($\epsilon=0.9$), respectivamente:
"""

# ╔═╡ ba5eba06-c713-4eee-88d6-be9a7a222043
md"""
###### Visualización de los datos:
"""

# ╔═╡ 01d77131-86c4-4ffc-bb05-0f3550b09409
begin
pt1 = plot(dft.q, dft.Res⭑, label="Respuesta Real", linewidth=0.5, marker=:star, color=:red, markersize=5, size=(550, 350), yscale=:log10, leg=:bottomright, ylimits = (1e0, 1e80), title="Estimación vs Respuesta Real L_9", titlefontsize=10, labelfontsize=8, legendfontsize=6);
	
plot!(dft.q, dft.Est, label="Estimación eps = 0.5", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:rect, color=:blue, alpha=0.5, markersize=3.5);
	
plot!(dfu.q, dfu.Est, label="Estimación eps = 0.7", linewidth=0.7, marker=:diamond, color=:green, alpha=0.5, markersize=3);

plot!(dfv.q, dfv.Est, label="Estimación eps = 0.9", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:circle, color=:orange, alpha=0.5, markersize=3.5);


pt2 = plot(dft.q, dft.R_err, label="eps = 0.5", color=:blue, xlabel="q", ylabel="Error relativo (%)", marker=:rect, size=(550, 350), title="Comparación errores relativos L_9", titlefontsize=10, labelfontsize=8, legendfontsize=6)

plot!(dfu.q, dfu.R_err, label="eps = 0.7", marker=:diamond, color=:green)
	
plot!(dfv.q, dfv.R_err, label="eps = 0.9", color=:orange, marker=:circle, secondary=true)


pt3 = bar(dft.q, dft.Time, label="eps = 0.5", color=:blue, xlabel="q", ylabel="Tiempo(seg)", size=(550, 350), bar_width=0.5, title="Comparación tiempos L_9", titlefontsize=10, labelfontsize=8, legendfontsize=6)

bar!(dfu.q, dfu.Time, label="eps = 0.7", color=:green, ylabel="Tiempo(seg)", size=(550, 350), bar_width=0.55)
	
bar!(dfv.q, dfv.Time, label="eps = 0.9", color=:orange, secondary=true, bar_width=0.6)

plot(pt1, pt3, pt2, layout=layout, size=(600, 400))
end

# ╔═╡ 645f3a10-0af1-4573-b5ec-d698ff5d12a2
md"""

"""

# ╔═╡ eb41a25c-7e74-4144-aaaf-6ba05e3097f0
md"""
$\textbf{k=10}$
"""

# ╔═╡ 5daa5df3-8aea-4f1a-ba38-753a7f7b6f9b
md"""
$$\epsilon=1$$

$\text{num\_simulations} = 1000000$
"""

# ╔═╡ 4f7530c9-4744-4c67-91f5-e77948053e3b
file_path_x = "../results/q-colorings/k=10/eps=1.csv"

# ╔═╡ 200eaef4-f500-4d78-9725-fe6183563943
dfx = CSV.File(file_path_x, delim=' ', header=true) |> DataFrame

# ╔═╡ 2337aaa1-820c-469e-85d8-0de47c26223b
md"""
Para este ejemplo, fijaremos el número de simulaciones independiente del epsilon (solo afectará el número de pasos del Gibbs sampler).

$$\epsilon=1$$

$\text{num\_simulations} = 100000$
"""

# ╔═╡ e8ff12e9-d637-4f81-ae2d-cd4120d3b127
file_path_y = "../results/q-colorings/k=10/sims=1e5.csv"

# ╔═╡ 9e17f7b5-95b6-4945-99e7-965bc65ce9b6
dfy = CSV.File(file_path_y, delim=' ', header=true) |> DataFrame;

# ╔═╡ cbee50c5-d31b-4538-9211-ab2d622a7373
md"""
Pasos del Gibbs para $2\leq q \leq 10$ ($\epsilon=0.1$), respectivamente:
"""

# ╔═╡ da129dbf-bca4-4ae8-b20a-a4232e8538e2
dfy.Gibbs

# ╔═╡ feab2f26-13b2-4aca-a2ad-8adc7c5eafcc
md"""
Para este ejemplo, fijaremos el número de simulaciones independiente del epsilon (solo afectará el número de pasos del Gibbs sampler).

$$\epsilon=1$$

$\text{num\_simulations} = 50000$
"""

# ╔═╡ a1257066-754e-4a38-a711-ec5af82207c0
file_path_z = "../results/q-colorings/k=10/sims=5e4.csv"

# ╔═╡ 045d40c6-deab-45c8-8995-8407d4832276
dfz = CSV.File(file_path_z, delim=' ', header=true) |> DataFrame;

# ╔═╡ 315a6845-7b3c-49c5-add7-42c88b063ab7
md"""
Pasos del Gibbs para $2\leq q \leq 10$ ($\epsilon=0.1$), respectivamente:
"""

# ╔═╡ 687f00b7-5fae-402a-9475-2e74c4a28256
dfz.Gibbs

# ╔═╡ 26241cd6-619c-4e8e-8182-55b0d4360238
md"""
###### Visualización de los datos:
"""

# ╔═╡ 2a712fb7-673a-4c9a-9aff-f5f882a4761c
begin
pu1 = plot(dfx.q, dfx.Res⭑, label="Respuesta Real", linewidth=0.5, marker=:star, color=:red, markersize=5, size=(550, 350), yscale=:log10, leg=:bottomright, ylimits = (1e0, 1e100), title="Estimación vs Respuesta Real L_10", titlefontsize=10, labelfontsize=8, legendfontsize=6);
	
plot!(dfx.q, dfx.Est, label="Estimación eps = 1", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:rect, color=:blue, alpha=0.5, markersize=3.5);
	
plot!(dfy.q, dfy.Est, label="Estimación num_sims = 1e5", linewidth=0.7, marker=:diamond, color=:green, alpha=0.5, markersize=3);

plot!(dfz.q, dfz.Est, label="Estimación num_sims = 5e4", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:circle, color=:orange, alpha=0.5, markersize=3.5);


pu2 = plot(dfx.q, dfx.R_err, label="eps = 1", color=:blue, xlabel="q", ylabel="Error relativo (%)", marker=:rect, size=(550, 350), title="Comparación errores relativos L_10", titlefontsize=10, labelfontsize=8, legendfontsize=6)

plot!(dfy.q, dfy.R_err, label="num_sims = 1e5", marker=:diamond, color=:green)
	
plot!(dfz.q, dfz.R_err, label="num_sims = 5e4", color=:orange, marker=:circle, secondary=true)


pu3 = bar(dfx.q, dfx.Time, label="eps = 1", color=:blue, xlabel="q", ylabel="Tiempo(seg)", size=(550, 350), bar_width=0.5, title="Comparación tiempos L_10", titlefontsize=10, labelfontsize=8, legendfontsize=6)

bar!(dfy.q, dfy.Time, label="num_sims = 1e5", color=:green, ylabel="Tiempo(seg)", size=(550, 350), bar_width=0.55)
	
bar!(dfz.q, dfz.Time, label="num_sims = 5e4", color=:orange, secondary=true, bar_width=0.6)

plot(pu1, pu3, pu2, layout=layout, size=(600, 400))
end

# ╔═╡ 0e6a2398-b99f-487d-ac8c-ef28e38b7739
md"""

"""

# ╔═╡ b3772b0e-3022-46ed-8ea7-a4fa64c01a43
md"""
$\textbf{k=11}$
"""

# ╔═╡ 086cd178-b2c0-4dd8-b4a8-5f22ad0718d0
md"""
$$\epsilon=1$$

$\text{num\_simulations} = 1771561$
"""

# ╔═╡ 1fbc185b-e7af-4a49-b2d5-a582bc3f658e
file_path_aa = "../results/q-colorings/k=11/eps=1.csv"

# ╔═╡ 42976fef-fb3c-4cb4-94a9-8cac5db8fce9
dfaa = CSV.File(file_path_aa, delim=' ', header=true) |> DataFrame

# ╔═╡ 9c06873d-60f9-4850-b789-ec85030bec12
md"""
Para este ejemplo, fijaremos el número de simulaciones independiente del epsilon (solo afectará el número de pasos del Gibbs sampler).

$$\epsilon=1$$

$\text{num\_simulations} = 100000$
"""

# ╔═╡ 061737cd-0d64-4e74-800d-47cddb3c6368
file_path_ab = "../results/q-colorings/k=11/sims=1e5.csv"

# ╔═╡ 7460ae37-8214-4a6e-bbec-a7016baa0320
dfab = CSV.File(file_path_ab, delim=' ', header=true) |> DataFrame;

# ╔═╡ 0b19a92d-d4c0-4735-a953-7e227ce51fd7
md"""
Pasos del Gibbs para $2\leq q \leq 10$ ($\epsilon=0.1$), respectivamente:
"""

# ╔═╡ 99f0881a-636c-4ada-8a9a-679719b7bbe7
dfab.Gibbs

# ╔═╡ d6ebfb27-e7d7-43ca-8ac9-864e5e43373e
md"""
Para este ejemplo, fijaremos el número de simulaciones independiente del epsilon (solo afectará el número de pasos del Gibbs sampler).

$$\epsilon=1$$

$\text{num\_simulations} = 50000$
"""

# ╔═╡ 4670de0d-ed22-4e12-b238-3d7176d74b08
file_path_ac = "../results/q-colorings/k=11/sims=5e4.csv"

# ╔═╡ c8382cfb-7046-471a-9a79-2b7022697deb
dfac = CSV.File(file_path_ac, delim=' ', header=true) |> DataFrame;

# ╔═╡ 55352d1c-1307-41ad-810b-da857f340f39
md"""
Pasos del Gibbs para $2\leq q \leq 10$ ($\epsilon=0.1$), respectivamente:
"""

# ╔═╡ 58919029-abbc-4c50-8d54-8bcddb71a848
dfac.Gibbs

# ╔═╡ 8447a741-d1ec-4117-89b7-6999e748381b
md"""
###### Visualización de los datos:
"""

# ╔═╡ 5d0280b4-da6d-40e9-b5cb-b26dab1b273c
begin
pv1 = plot(dfaa.q, dfaa.Res⭑, label="Respuesta Real", linewidth=0.5, marker=:star, color=:red, markersize=5, size=(550, 350), yscale=:log10, leg=:bottomright, ylimits = (1e0, 1e115), title="Estimación vs Respuesta Real L_11", titlefontsize=10, labelfontsize=8, legendfontsize=6);
	
plot!(dfaa.q, dfaa.Est, label="Estimación eps = 1", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:rect, color=:blue, alpha=0.5, markersize=3.5);
	
plot!(dfab.q, dfab.Est, label="Estimación num_sims = 1e5", linewidth=0.7, marker=:diamond, color=:green, alpha=0.5, markersize=3);

plot!(dfac.q, dfac.Est, label="Estimación num_sims = 5e4", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:circle, color=:orange, alpha=0.5, markersize=3.5);


pv2 = plot(dfaa.q, dfaa.R_err, label="eps = 1", color=:blue, xlabel="q", ylabel="Error relativo (%)", marker=:rect, size=(550, 350), title="Comparación errores relativos L_11", titlefontsize=10, labelfontsize=8, legendfontsize=6)

plot!(dfab.q, dfab.R_err, label="num_sims = 1e5", marker=:diamond, color=:green)
	
plot!(dfac.q, dfac.R_err, label="num_sims = 5e4", color=:orange, marker=:circle, secondary=true)


pv3 = bar(dfaa.q, dfaa.Time, label="eps = 1", color=:blue, xlabel="q", ylabel="Tiempo(seg)", size=(550, 350), bar_width=0.5, title="Comparación tiempos L_11", titlefontsize=10, labelfontsize=8, legendfontsize=6)

bar!(dfab.q, dfab.Time, label="num_sims = 1e5", color=:green, ylabel="Tiempo(seg)", size=(550, 350), bar_width=0.55)
	
bar!(dfac.q, dfac.Time, label="num_sims = 5e4", color=:orange, secondary=true, bar_width=0.6)

plot(pv1, pv3, pv2, layout=layout, size=(600, 400))
end

# ╔═╡ 3cce1304-c75e-484c-bee7-20ab999709d0
md"""

"""

# ╔═╡ 3eda40e9-31f2-48e5-9be4-88a960a0c65f
md"""
$\textbf{k=12}$
"""

# ╔═╡ 41e33e79-7c01-4b65-9e79-718ddb1645e2
md"""
Para este ejemplo, fijaremos el número de simulaciones independiente del epsilon (solo afectará el número de pasos del Gibbs sampler).

$$\epsilon=1$$

$\text{num\_simulations} = 1000000$
"""

# ╔═╡ fa5dabf7-8f0b-4f4d-afe7-9635d6cc1769
file_path_ad = "../results/q-colorings/k=12/sims=1e6.csv"

# ╔═╡ e21273e1-e6e9-45c2-b99e-1dda86c386c9
dfad = CSV.File(file_path_ad, delim=' ', header=true) |> DataFrame

# ╔═╡ 55eb1b17-0bb6-4774-92cf-03e2c67aba0c
md"""
Para este ejemplo, fijaremos el número de simulaciones independiente del epsilon (solo afectará el número de pasos del Gibbs sampler).

$$\epsilon=1$$

$\text{num\_simulations} = 100000$
"""

# ╔═╡ 67ba42cd-04cc-49d4-a490-d56662d80e03
file_path_ae = "../results/q-colorings/k=12/sims=1e5.csv"

# ╔═╡ 616768f1-e52f-4f26-beea-8f085b812786
dfae = CSV.File(file_path_ae, delim=' ', header=true) |> DataFrame;

# ╔═╡ d27853ed-346f-427e-9ed9-8683620a6114
md"""
Pasos del Gibbs para $2\leq q \leq 10$ ($\epsilon=1$), respectivamente:
"""

# ╔═╡ b479a379-e7de-46a4-a8b4-3e0ce1464424
dfae.Gibbs

# ╔═╡ 32d541bf-0b7a-4eb8-b17f-597e04c1aaad
md"""
Para este ejemplo, fijaremos el número de simulaciones independiente del epsilon (solo afectará el número de pasos del Gibbs sampler).

$$\epsilon=1$$

$\text{num\_simulations} = 50000$
"""

# ╔═╡ 225815b8-9b2e-4a7f-814f-89f6adf2e8f6
file_path_af = "../results/q-colorings/k=12/sims=5e4.csv"

# ╔═╡ a034a1a7-3af2-4f5d-9897-3e77a67f7568
dfaf = CSV.File(file_path_af, delim=' ', header=true) |> DataFrame;

# ╔═╡ 91d659a4-5286-4476-86f5-f805268c88fe
md"""
Pasos del Gibbs para $2\leq q \leq 10$ ($\epsilon=1$), respectivamente:
"""

# ╔═╡ 8457bb11-126e-401a-ae53-560381ce2ef4
dfaf.Gibbs

# ╔═╡ abe3a1f2-2ab4-456b-8dfb-f809b3f3ef05
md"""
###### Visualización de los datos:
"""

# ╔═╡ e2612b27-cc9b-4de7-bce0-1a5224d45bb4
begin
px1 = plot(dfad.q, dfad.Res⭑, label="Respuesta Real", linewidth=0.5, marker=:star, color=:red, markersize=5, size=(550, 350), yscale=:log10, leg=:bottomright, ylimits = (1e0, 1e140), title="Estimación vs Respuesta Real L_12", titlefontsize=10, labelfontsize=8, legendfontsize=6);
	
plot!(dfad.q, dfad.Est, label="Estimación num_sims = 1e6", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:rect, color=:blue, alpha=0.5, markersize=3.5);
	
plot!(dfae.q, dfae.Est, label="Estimación num_sims = 1e5", linewidth=0.7, marker=:diamond, color=:green, alpha=0.5, markersize=3);

plot!(dfaf.q, dfaf.Est, label="Estimación num_sims = 5e4", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:circle, color=:orange, alpha=0.5, markersize=3.5);


px2 = plot(dfad.q, dfad.R_err, label="num_sims = 1e6", color=:blue, xlabel="q", ylabel="Error relativo (%)", marker=:rect, size=(550, 350), title="Comparación errores relativos L_12", titlefontsize=10, labelfontsize=8, legendfontsize=6)

plot!(dfae.q, dfae.R_err, label="num_sims = 1e5", marker=:diamond, color=:green)
	
plot!(dfaf.q, dfaf.R_err, label="num_sims = 5e4", color=:orange, marker=:circle, secondary=true)


px3 = bar(dfad.q, dfad.Time, label="num_sims = 1e6", color=:blue, xlabel="q", ylabel="Tiempo(seg)", size=(550, 350), bar_width=0.5, title="Comparación tiempos L_12", titlefontsize=10, labelfontsize=8, legendfontsize=6)

bar!(dfae.q, dfae.Time, label="num_sims = 1e5", color=:green, ylabel="Tiempo(seg)", size=(550, 350), bar_width=0.55)
	
bar!(dfaf.q, dfaf.Time, label="num_sims = 5e4", color=:orange, secondary=true, bar_width=0.6)

plot(px1, px3, px2, layout=layout, size=(600, 400))
end

# ╔═╡ fc529c7f-24d1-4fd9-a72d-1a757b6d18da
md"""
$\textbf{k=13}$
"""

# ╔═╡ 123f5380-ca96-40e8-bf03-9386e7f165fa
md"""
Para este ejemplo, fijaremos el número de simulaciones independiente del epsilon (solo afectará el número de pasos del Gibbs sampler).

$$\epsilon=1$$

$\text{num\_simulations} = 1000000$
"""

# ╔═╡ a40a4249-bf07-4642-ab41-407315cc4e45
file_path_ag = "../results/q-colorings/k=13/sims=1e6.csv"

# ╔═╡ dbcc9b76-b38d-4b86-a356-7fc85daa04ec
dfag = CSV.File(file_path_ag, delim=' ', header=true) |> DataFrame

# ╔═╡ 6fee85b2-61e0-4f09-9c77-88981cd0be7a
md"""
Para este ejemplo, fijaremos el número de simulaciones independiente del epsilon (solo afectará el número de pasos del Gibbs sampler).

$$\epsilon=1$$

$\text{num\_simulations} = 100000$
"""

# ╔═╡ f4c13240-4562-4ca3-bffd-d9010b83deea
file_path_ah = "../results/q-colorings/k=13/sims=1e5.csv"

# ╔═╡ 6a38ddfa-b2ac-45fb-af00-baf09220f5bc
dfah = CSV.File(file_path_ah, delim=' ', header=true) |> DataFrame;

# ╔═╡ 3675225d-299a-4181-a3e9-d86e64e470b2
md"""
Pasos del Gibbs para $2\leq q \leq 10$ ($\epsilon=1$), respectivamente:
"""

# ╔═╡ 78f3a968-71ac-484c-9347-f0d47cc5de5a
dfah.Gibbs

# ╔═╡ 96e155ef-1159-4474-9737-798b57a62591
md"""
Para este ejemplo, fijaremos el número de simulaciones independiente del epsilon (solo afectará el número de pasos del Gibbs sampler).

$$\epsilon=1$$

$\text{num\_simulations} = 50000$
"""

# ╔═╡ 951a9d02-3321-46d9-97c3-f90eda11bb55
file_path_ai = "../results/q-colorings/k=13/sims=5e4.csv"

# ╔═╡ 757ef5e3-8ea3-4a44-84ea-c8035e619cb8
dfai = CSV.File(file_path_ai, delim=' ', header=true) |> DataFrame;

# ╔═╡ 546283f0-7bde-40dc-a040-0d468b50bb84
md"""
Pasos del Gibbs para $2\leq q \leq 10$ ($\epsilon=1$), respectivamente:
"""

# ╔═╡ c259b3e5-3786-41e4-a801-585e09dd0944
dfai.Gibbs

# ╔═╡ 53580e32-14df-4392-b3ce-5aaa6ccf7e0a
md"""
###### Visualización de los datos:
"""

# ╔═╡ c3bb759b-ef5e-4569-be23-a664a775e714
begin
py1 = plot(dfag.q, dfag.Res⭑, label="Respuesta Real", linewidth=0.5, marker=:star, color=:red, markersize=5, size=(550, 350), yscale=:log10, leg=:bottomright, ylimits = (1e0, 1e162), title="Estimación vs Respuesta Real L_13", titlefontsize=10, labelfontsize=8, legendfontsize=6);
	
plot!(dfag.q, dfag.Est, label="Estimación num_sims = 1e6", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:rect, color=:blue, alpha=0.5, markersize=3.5);
	
plot!(dfah.q, dfah.Est, label="Estimación num_sims = 1e5", linewidth=0.7, marker=:diamond, color=:green, alpha=0.5, markersize=3);

plot!(dfai.q, dfai.Est, label="Estimación num_sims = 5e4", xlabel="q", ylabel="Valor", linewidth=0.7, marker=:circle, color=:orange, alpha=0.5, markersize=3.5);


py2 = plot(dfag.q, dfag.R_err, label="num_sims = 1e6", color=:blue, xlabel="q", ylabel="Error relativo (%)", marker=:rect, size=(550, 350), title="Comparación errores relativos L_13", titlefontsize=10, labelfontsize=8, legendfontsize=6)

plot!(dfah.q, dfah.R_err, label="num_sims = 1e5", marker=:diamond, color=:green)
	
plot!(dfai.q, dfai.R_err, label="num_sims = 5e4", color=:orange, marker=:circle, secondary=true)


py3 = bar(dfag.q, dfag.Time, label="num_sims = 1e6", color=:blue, xlabel="q", ylabel="Tiempo(seg)", size=(550, 350), bar_width=0.5, title="Comparación tiempos L_13", titlefontsize=10, labelfontsize=8, legendfontsize=6)

bar!(dfah.q, dfah.Time, label="num_sims = 1e5", color=:green, ylabel="Tiempo(seg)", size=(550, 350), bar_width=0.55)
	
bar!(dfai.q, dfai.Time, label="num_sims = 5e4", color=:orange, secondary=true, bar_width=0.6)

plot(py1, py3, py2, layout=layout, size=(600, 400))
end

# ╔═╡ 1c18d0c0-d515-420c-8f10-4be6d5abfaeb
md"""
**Nota:** No consideramos $k>13$ pues no encontramos resultados exactos sobre estos valores para realizar comparaciones, y como mencionamos anteriormente calcular dichos polinomios cromáticos ya no es una opción muy viable. Además, los tiempos de ejecución de los programas ya estaban siendo altos incluso para $\epsilon\approx 1$ (Aunque tiempos aún bajos realizando un número de simulaciónes fijo de hasta $1e5$, obteniendo errores relativos en su mayoría por debajo de $4\%$ -varios por debajo del $2\%$-).
"""

# ╔═╡ b28683c5-f7a5-497e-9558-b15f27d69c8a
md"""
---
"""

# ╔═╡ c7bba1dc-6b77-4250-bfd6-8051914d7cbb
html"""
<h4>Promedio de los ratios:</h4>
"""

# ╔═╡ be3186c4-c5c9-45f5-9050-4856b70bb2e8
begin
plot(dfe.q, dfa.Mean_r, label="k = 2", color=:blue, ylabel="Promedio ratio", marker=:circle, size=(600, 400), alpha = 0.5, leg=:right, titlefontsize=10, labelfontsize=8, legendfontsize=7)

plot!(dff.q, dfd.Mean_r, label="k = 3", marker=:utriangle, color=:green, alpha = 0.5)
	
plot!(dfe.q, dfe.Mean_r, label="k = 4", color=:orange, marker=:rect, secondary=true, alpha = 0.5)

plot!(dfe.q, dfh.Mean_r, label="k = 5", color=:violet, marker=:pentagon, secondary=true, alpha = 0.5)

plot!(dfe.q, dfk.Mean_r, label="k = 6", color=:red, marker=:hexagon, secondary=true, alpha = 0.5)

plot!(dfe.q, dfn.Mean_r, label="k = 7", color=:pink, marker=:heptagon, secondary=true, alpha = 0.5)

plot!(dfe.q, vcat(dfq.Mean_r, [0.90949, 0.91696, 0.92336, 0.92875, 0.93349]), label="k = 8", color=:brown, marker=:octagon, secondary=true, alpha = 0.3)

plot!(dfe.q, vcat(dft.Mean_r, [0.90948, 0.91697, 0.92335, 0.92874, 0.93349]), label="k = 9", color=:cyan, marker=:star, secondary=true, alpha = 0.3)

plot!(dfe.q, vcat(dfx.Mean_r, [0.90952, 0.91702, 0.9233, 0.92876, 0.93348]), label="k = 10", color=:blue3, marker=:star6, secondary=true, alpha = 0.3)

plot!(dfe.q, vcat(dfaa.Mean_r, [0.90948, 0.91695, 0.92333, 0.92877, 0.93352]), label="k = 11", color=:red3, marker=:xcross, secondary=true, alpha = 0.3)

plot!(dfe.q, vcat(dfad.Mean_r, [0.90949, 0.91697, 0.92332, 0.92877, 0.9335,]), label="k = 12", color=:orange3, marker=:star8, secondary=true, alpha = 0.3)

plot!(dfe.q, vcat(dfag.Mean_r, [0.90948, 0.91698, 0.92329, 0.92871, 0.93351]), label="k = 13", color=:orange3, marker=:diamond, secondary=true, alpha = 0.3)

title!("Comparación promedios de los ratios")
xlabel!("q")
end

# ╔═╡ 21580099-c127-46c7-917f-84edbf446b80
md"""
Según estos datos podríamos decir que el promedio de los ratios es independiente de $k$ (al menos para $5\leq q \leq 15$ y $2\leq k \leq 13$), ¿Se mantiene esto para cualquier $k > 13$ y cualquier $q>15$?

Teniendo esto en cuenta, a priori podemos realizar las siguientes estimaciones para el número de $q$-coloraciones ($5\leq q \leq 15$) de un lattice $k\times k$ ($k>13$) con los valores calculados del promedio de los ratios como:

$q^{n} * (r_{q})^{m}$

donde $r_{q}$ es algún promedio de los ratios para el color $q$ en algún $k$ calculado anteriormente. Para la siguiente muestra tomamos los promedios anteriormente calculados para $k=9$.
"""

# ╔═╡ d51c4084-8255-420c-84d2-c0a908f8e7ed
file_path_est = "../results/q-colorings/k=14..20/est.csv"

# ╔═╡ 55a0d112-5c54-4bbf-931c-0fa930c0196c
begin
dfest1 = CSV.File(file_path_est, delim=' ', header=true, types=String) |> DataFrame;
dfest1[!, names(dfest1)[1]] = parse.(Int, dfest1[!, names(dfest1)[1]])
for col in names(dfest1)[2:12]
    dfest1[!, col] = BigFloat.(dfest1[!, col])
end
end

# ╔═╡ d4bcc27a-d2f5-4e58-b6e4-29ea7b0d4eeb
dfest1

# ╔═╡ 76acfab3-3e3d-4413-9b57-7ff763b32084
md"""
Algunas gráficas de los datos obtenidos:
"""

# ╔═╡ 0d07bf65-32b8-4569-9539-9f3124ed2322
begin
lim = BigFloat("1e400")
qs1 = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
layout2 = @layout [a b; c d]
paa1 = plot(qs1, collect(dfest1[1, :][2:12]), label="Estimación", linewidth=0.5, marker=:circle, color=:orange, markersize=5, size=(550, 350), yscale=:log10, leg=:bottomright, ylimits = (1e80, 1e200), title="Estimación q-coloraciones L_13", titlefontsize=10, labelfontsize=8, legendfontsize=7, xlabel ="q", ylabel="Valor");

paa2 = plot(qs1, collect(dfest1[4, :][2:12]), label="Estimación", linewidth=0.5, marker=:rect, color=:green, markersize=5, size=(550, 350), yscale=:log10, leg=:bottomright, ylimits = (1e100, 1e300), title="Estimación q-coloraciones L_16", titlefontsize=10, labelfontsize=8, legendfontsize=7, xlabel ="q", ylabel="Valor");

paa3 = plot(qs1, (collect(dfest1[3, :][2:12])), label="Estimación", linewidth=0.5, marker=:diamond, color=:blue, markersize=5, size=(550, 350), yscale=:log10, leg=:bottomright, title="Estimación q-coloraciones L_15", titlefontsize=10, labelfontsize=8, legendfontsize=7, xlabel ="q", ylabel="Valor (log10)");

paa4 = plot(qs1, log10.(collect(dfest1[8, :][2:12])), label="Estimación", linewidth=0.5, marker=:hexagon, color=:pink, markersize=5, size=(550, 350), yscale=:log10, leg=:bottomright, title="Estimación q-coloraciones L_20", titlefontsize=10, labelfontsize=8, legendfontsize=7, xlabel ="q", ylabel="Valor (log10)");

plot(paa1, paa3, paa2, paa4, layout=layout2, size = (800, 600))
end

# ╔═╡ d3f70173-3ae7-4d46-b383-b703b15a0483
md"""
### Conclusiones:
"""

# ╔═╡ 9009028e-4bd2-46c9-a47c-c795070ee18f
md"""
- El método aproximado demostró ser efectivo y eficiente, comparado con métodos exactos que pueden ser computacionalmente inviables para grafos grandes.

- Los errores relativos fueron generalmente bajos, lo que sugiere que el método de aproximación, a pesar de su naturaleza estocástica, proporciona estimaciones razonablemente precisas.

- La implementación logró resultados útiles en tiempos razonables incluso para configuraciones de grafo más grandes, destacando la utilidad práctica del algoritmo en el problema de contar las $q$-coloraciones de un grafo reticular.

- El algoritmo proporciona un control razonable sobre el error, lo que permite confiar en la precisión de las estimaciones para la toma de decisiones.

- La variación en el número de simulaciones y los pasos del sampler muestra que el algoritmo es sensible a estos parámetros, lo que subraya la importancia de una configuración cuidadosa para maximizar la eficiencia y la efectividad.

- Existe espacio para optimizar aún más el algoritmo ajustando dinámicamente los parámetros en respuesta a las características del grafo y los resultados intermedios de las simulaciones, incluso se podría hacer uso de computación paralela para implementar varios hilos al momento de realizar las estimaciones de los ratios (esta probablemente sería una optimización clave y como anotación, un compañero nuestro ajeno a la clase lo está intentando hacer).
"""

# ╔═╡ 0f699a12-b72a-4074-8771-95f650b68a25
html"""
<h2 style="text-align: center;">Segundo Punto</h3>
"""

# ╔═╡ 8c8e1e48-3656-4bdb-9759-c88a55e08cda
md"""
Obtenga valores aproximados para el número de configuraciones factibles en el modelo "Hard Core" para lattices $k\times k$ con $2\leq k \leq 20$.

Reporte de manera similar a lo hecho en el item **(a)** del ejercicio anterior.
"""

# ╔═╡ d22a515c-fc00-428d-9f36-846295efc210
html"""
<h3 style="text-align: center;">Solución:</h3>
"""

# ╔═╡ 0c0f8f82-179c-433b-b570-41100a74bd37
md"""
Antes de solucionar el punto, en el modelo Hard Core sobre un grafo $G = (V, E)$ se asigna de manera aleatoria el valor de $0$ o $1$ a cada uno de los vértices, de tal manera que si dos vértices son adyacentes no pueden tomar el valor de $1$ simultáneamente. Las asignaciones de ceros y unos son llamadas configuraciones y aquellas que cumplan la anterior condición son configuraciones factibles.
"""

# ╔═╡ 50519776-3063-49f3-a7fb-a3949ac1ee99
html"""
<div style="text-align: center;">
<img src="https://i.ibb.co/MD82nSD/2f.jpg"
	width="300" class=center alt="Ejemplo 3-coloración del retículo 3x3">
</div>
<div style="text-align: center;">
(Configuraciones factibles del grafo reticular 2x2)
</div>
"""

# ╔═╡ 458760be-de62-44f2-af1f-a642d831e2e1
md"""
Buscando en internet no encontramos resultados sobre las respuestas exactas a este problema, es por esto que implementamos una algoritmo utilizando [programación dinámica](https://es.wikipedia.org/wiki/Programación_dinámica) logrando encontrar respuestas hasta $k=25$ en un tiempo relativamente moderado. Al final del notebook hablamos sobre este acercamiento.
"""

# ╔═╡ a645728d-3627-4f5f-b01e-321d1db730cb
file_path_reshc = "../results/hard-core/dpAnswers/out.csv"

# ╔═╡ 52387ede-862c-452a-9afa-01eb9cbc6943
dfreshc = CSV.File(file_path_reshc, delim=' ', header=true) |> DataFrame

# ╔═╡ b4401e0b-7df8-44b5-a992-aa703f9a4145
md"""
### Solución mediante programación dinámica:
"""

# ╔═╡ 53a52203-2882-4c6b-81fe-6d4e4a698efe
md"""
Cada fila del grafo reticular $k\times k$ puede ser representada como un número binario de longitud $k$, donde cada bit indica si el valor del vértice es $1$ o $0$ (esto se conoce como una máscara de bits).

**Estados:**

$$\text{dp}[i][\textit{mask}]$$ 

representará el número de configuraciones factibles para las primeras $i$ filas, donde la fila $i$ tiene la configuración representada por $\textit{mask}$.

**Casos base:**

$\text{dp}[0][\textit{mask}] = 1$ 

para toda máscara válida $\textit{mask}$ que no tenga bits adyacentes prendidos, esto se puede verificar comprobando que `mask&(mask>>1) == 0`, donde `&` representa la operación bit a bit *and* y *>>* representa la operación de corrimiento de bits hacia la derecha. El resto de estados tendrá el valor inicial de $0$.

**Transiciones:**

Para cada posible máscara de la fila actual $i$, `cur_mask`, y cada posible máscara de la fila anterior $i-1$, `old_mask`: Si `cur_mask` y `old_mask` son válidas, y adicionalmente son compatibles (no hay unos adyacentes verticalmente), entonces 

$$\text{dp}[i][\textit{cur\_mask}] = \text{dp}[i][\textit{cur\_mask}] + \text{dp}[i-1][\textit{old\_mask}].$$

**Respuesta final:**

$$\text{Res} = \sum_{\textit{mask}}dp[k-1][\textit{mask}]$$
"""

# ╔═╡ 1724968f-2a3f-48f2-9be0-b829a8e20f41
md"""
Veamos el ejemplo más sencillo para $k=2$:

**Casos base:**

- Posibles máscaras fila $0$: $\{\color{green}00 \color{black},\color{green} 01\color{black}, \color{green}10\color{black}, \color{red}{11} \color{black}\}$

$\text{dp}[0][00] = \text{dp}[0][01] = \text{dp}[0][10] = 1$

"""

# ╔═╡ 04c2cf76-bd83-400e-b86e-ec44589249a0
md"""
**Transiciones:**
- Posibles máscaras fila $1$: $\{\color{green}00 \color{black},\color{green} 01\color{black}, \color{green}10\color{black}, \color{red}{11} \color{black}\}$

$\text{dp}[1][\color{orange}00\color{black}] = \text{dp}[0][00] + \text{dp}[0][01] + \text{dp}[0][10] = 3$

$\begin{matrix}
0 & 0\\
\color{orange}0 & \color{orange}0 
\end{matrix}$

$\begin{matrix}
0 & 1\\
\color{orange}0 & \color{orange}0 
\end{matrix}$

$\begin{matrix}
1 & 0\\
\color{orange}0 & \color{orange}0 
\end{matrix}$

$\text{dp}[1][\color{orange}01\color{black}] = \text{dp}[0][00] + \text{dp}[0][10] = 2$

$\begin{matrix}
0 & 0\\
\color{orange}0 & \color{orange}1
\end{matrix}$

$\begin{matrix}
1 & 0\\
\color{orange}0 & \color{orange}1
\end{matrix}$

$\text{dp}[1][\color{orange}10\color{black}] = \text{dp}[0][00] + \text{dp}[0][01] = 2$

$\begin{matrix}
0 & 0\\
\color{orange}1 & \color{orange}0
\end{matrix}$

$\begin{matrix}
0 & 1\\
\color{orange}1 & \color{orange}0
\end{matrix}$
"""

# ╔═╡ f63d0a79-667f-4f2c-9d49-fc0acfc0250a
md"""
**Respuesta:**

$\text{Res} = \sum_{\textit{mask}} \text{dp}[1][\textit{mask}] = \text{dp}[1][00] + \text{dp}[1][01] + \text{dp}[1][10] + \text{dp}[1][11] = 7$
"""

# ╔═╡ f5a4937d-d4e1-4ec0-9eb9-f22052c48ed1
md"""
**Complejidad:**
- Estados: Cada fila se representa con una `bitmask` de $k$ bits. Hay $2^k$ posibles configuraciones (máscaras) para cada fila.

- Transiciones: Para cada fila se consideran todas las posibles configuraciones y se verifica su compatibilídad con todas las configuraciones de la fila anterior. Así, para cada una de las $k$ filas la transición involucra comparar $2^k\times 2^k = 4^k$ pares de máscaras.

Es decir que la complejidad es de $O(k\times 4^k)$, permitiéndonos calcular en un tiempo razonable máximo hasta $k=20$.

Podemos optimizar de la siguiente manera:

- Hay $2^k$ posibles configuraciones para cada fila. Pero hay $l<2^k$ posibles configuraciones factibles para cada fila. Podemos precalcularlas y guardarlas ($O(2^k)$).

- Para cada máscara válida, podemoslos determinar las máscaras que son compatibles con ella (verticalmente). Esto también puede ser precalculado, guardando en una lista las máscaras compatibles para cada máscara válida ($O(l^2)$).

- Para las transiciones, utilizamos la lista precalculada.

Es decir que la complejidad bajaría de $O(k\times4^k)$ a $O(2^k + kl^2)$, -no es mucho pero es trabajo honesto-.

Observación: $l = F(k+2)$ donde $F(i)$ es el $i$-ésimo número de Fibonacci. Teniendo en cuenta que $F(k+2) \approx \frac{\phi^{k+2}}{\sqrt{5}}$, entonces podemos calcular en tiempo razonable (menos de $20$ minutos) hasta $k=25$.
"""

# ╔═╡ b81d1973-cf1c-41d5-9e60-d4ac0f2ecf21
md"""
#### Implementación:
"""

# ╔═╡ 59ded85c-0ef8-418d-abc9-28b0ddfd8afd
md"""
Por facilidad realizamos la implementación en C++, esta se encuentra en el [GitHub](https://github.com/trodrigueza/MCMC_Colorings) en *cpp/hard-core/dp.cpp*.
"""

# ╔═╡ e1c064db-c1f5-4f83-8b00-cd2a48021068
md"""
```cpp
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

  validMasks.clear();
  compatible.clear();

  precompute(k);

  dp[0].assign(validMasks.size(), 0);
  dp[1].assign(validMasks.size(), 0);

  int curr = 0, prev = 1;

  // base cases
  for (int i = 0; i < validMasks.size(); i++) {
    dp[curr][i] = 1;
  }

  // transitions
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

  // calculate answer
  big_int total = 0;
  for (int i = 0; i < validMasks.size(); i++) {
    total += dp[curr][i];
  }

  cout << "Number of feasible Hard Core configurations for a " << k << "x" << k << " lattice: " << setprecision(5) << total.convert_to<big_float>() << "\n";
}
```
"""

# ╔═╡ Cell order:
# ╟─21f2ee62-321d-47fe-a2b0-aa3ff2d558f8
# ╟─53933472-805f-4ef9-acd6-d3fe44bb95ea
# ╟─7a35759c-a1c7-47d0-8eec-947dcef2eaa9
# ╟─7d1dc16a-b983-4778-a646-9dedabae8f31
# ╟─1d2dcd3a-2cf8-4d07-8e60-83272ffe418b
# ╟─7b304895-df81-4b98-9a1d-238bd55532b4
# ╟─9fc0084a-2770-4821-929a-3de1a93bb92f
# ╟─ee3e5859-4947-49c1-8ed1-edc5987a4582
# ╟─50bfb720-cf21-438f-ac1f-0a4134729478
# ╟─609cdee4-0e35-42f8-9b78-44d427d0da12
# ╟─3a69b474-7ea9-494c-a236-3c50413947f3
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
# ╟─cf84e4b2-99cb-4e3f-82bc-7fbf12d5bcfb
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
# ╟─a9ba504b-3ecc-4e80-b063-06242abb23ca
# ╠═6fdef38f-a45a-4cfe-9694-a0ac144f16b9
# ╠═fdc76661-9c2e-47f0-8716-a3d719e88998
# ╟─aa182161-2906-4e5d-aa08-3153c19f97e0
# ╠═57afe92f-88a4-4fdc-83db-5d58cc525e73
# ╟─205399df-28be-4674-a1c2-35bc5db6e765
# ╟─62e7c451-ba95-496c-8f16-5d7deea3388a
# ╠═899c7974-3284-4eb8-bf7e-7d5a3028ec04
# ╠═52042204-a585-46b1-9757-694b9f1de30c
# ╟─0d70f4d5-08d3-43bc-ae1d-eb856ccfb248
# ╠═1934af33-79b9-430d-8667-4fb9cc99bddf
# ╟─33cf49f4-8b88-4d6f-a3d0-42a778836013
# ╟─00aa11ee-8a3d-424e-9ea5-3c8a030c7d0b
# ╟─d4d9a98a-01ac-478c-a787-37c15dc669e6
# ╟─00f6943d-4bbc-44a5-b020-8ebc4287f2d0
# ╟─23b80981-5bc0-407d-b9f1-8ebd8c8eda5b
# ╟─1d6f351c-4409-4d69-adc8-dce1a961ac82
# ╟─6040e5b2-8428-4808-bba4-586ee7bed8de
# ╟─169716b4-0e42-44a2-ac7a-994b8ab517e3
# ╟─3589ae6a-4ec2-4aa0-a87b-08baf508d260
# ╟─4eec39a2-deb3-49ad-91d4-839a4fe02f33
# ╟─cc0ef645-3313-4ca8-92a9-7a7f6e3cf53d
# ╠═3a058fac-395e-48da-a045-2f9496714327
# ╟─a07a0f97-56a7-4577-822a-54734af93229
# ╟─6aaa555a-833c-4662-bfd7-249bd1c93a3e
# ╟─e47cb888-c040-4ad4-a7f9-2c6c1604b1b5
# ╟─cd636c43-de57-4eef-a23f-d1eb3eb3f8e0
# ╟─84c216c7-3e28-472b-b50d-7b3cfd307a81
# ╟─ea51c9cc-ce50-4758-83f2-3240779604f4
# ╟─3c70a826-5f3f-45cf-a054-7d947392f694
# ╟─764aa6f1-f3c3-4cbd-a7ba-d0614eee0093
# ╠═65e5e18e-9df3-4dd6-94f4-218fc9b1c694
# ╟─c844196d-f3ac-4ce5-b972-f2fd1682fa13
# ╟─58473d1b-26e3-40f7-9c37-a0bdf3e07990
# ╟─7c777c46-1b5f-4257-9706-f181a08c7e5a
# ╟─e45ce12c-2a75-49b0-a2e6-b96a8467cd58
# ╟─7b14e5d2-2774-4876-9458-e631a1615447
# ╟─aaf90bb9-c2b7-4e5a-980f-c9580f565ccd
# ╟─4a8d77dc-f16d-4c93-94a4-c2e1f3f6e1ef
# ╟─ca36d1eb-86bc-473e-9bbe-dac1cfc13fa6
# ╟─b5f024b3-aa26-479e-b4c0-4535696b8f52
# ╟─5184161f-547b-4cfb-94e6-f55b9765fd76
# ╟─33b74cac-6dcf-4211-882e-0dc6abd4c435
# ╟─2a924600-584c-41e0-b344-43c1e34411a8
# ╟─ada44c79-3b92-4be9-8b0f-14c09c7547a0
# ╟─c79341d1-13d3-4e71-9297-a8df019a63d9
# ╠═f771ff8f-6e39-4fe7-9f9d-a7c9b4191fc5
# ╟─3c858f7b-51b7-43c0-9c34-1c9294973825
# ╟─9675b48f-5468-4030-8117-d090d2074d6b
# ╟─f307b7c6-8350-4d5e-bafe-b893b9de5407
# ╟─fd195252-d166-44bc-95a8-2280201c1bb8
# ╟─9c29cb2c-d2c3-4815-ad8e-c53a7d34722d
# ╟─aa14e81b-0530-4579-8e97-c62d49046c60
# ╟─e14facff-e9d5-49f5-84ae-2cd8626d5986
# ╟─4465b037-9586-4439-af42-3755aad96e31
# ╟─5ac758b9-5cda-40c5-b829-4603b00aef03
# ╟─ce6c8449-c44d-4309-b5e1-00b4adf35950
# ╟─971394c0-4918-4e37-9324-54333cd3816e
# ╟─38894819-bfdc-4ee5-a3c8-237f6dbd400e
# ╟─03cb0dac-f17d-46e4-9578-8a72df036a5e
# ╟─4ea0cb09-d6db-4e0c-b844-9ba06ff54048
# ╠═5acace06-a72e-4b8e-a019-eb75f1b80fe9
# ╟─f5d128d4-8f46-4196-803b-0de3b90b7bdc
# ╟─fd746cac-f3a7-4828-8361-b1d4194e0f5d
# ╟─b5395b0e-1c9b-489a-b903-642838f15d97
# ╟─e56ee173-8e24-4f88-ab3b-6d2dfe3f51d1
# ╠═eeb80031-9afc-4aa2-9b86-66199dcf9746
# ╟─e28113c4-6e57-486b-b58a-9f380db8485a
# ╟─0a81a4a9-7d99-4ea4-b698-f3d907f12b20
# ╟─23a04d74-d7be-4b7a-bb11-d9a64a469685
# ╟─5d04d261-896b-4fad-8c75-76af8037c9f0
# ╟─395d9d7d-93c0-4d35-8bec-a44b030e93a4
# ╟─3888f02c-a7e8-4959-a209-fbed5376ae67
# ╟─280be26f-635f-4382-8870-fdf0e09223aa
# ╟─eaf519b1-2c36-461f-89a8-ad2ba52ab3c4
# ╟─55cf5973-1db4-492d-b1a7-e237119c69e0
# ╟─dd5f08d6-d743-4d41-b48e-77e70c901981
# ╟─e16f67a2-ef40-4bbf-9b28-23a0ca3c5ff5
# ╟─816d2f18-60b7-42b9-9e52-c1f2d9151db4
# ╠═d35643f4-7bbc-4f92-adc2-49f1df506be3
# ╟─40b195f2-31d8-4e6e-a1f7-367df2ede70f
# ╟─cc6826d1-82b9-4e8d-af67-9aa34e6b76ae
# ╟─79913af6-060a-4823-843a-d591139df7e4
# ╟─e16423d5-3456-4f53-acbd-4ede3c9b7ca4
# ╠═5fb677e9-8353-45fe-a0ac-467b0bb39f37
# ╟─c94b11d5-3c2b-4aa4-932e-97c1a5dade17
# ╟─8d89bf10-0f60-42d2-bcbf-6d020c22d2ca
# ╟─844bfa80-bbd8-4ca9-929a-185c4b9cc536
# ╟─912ae862-d93b-4f6b-bbff-106e106cb4e5
# ╟─78c00185-8dee-434b-b6a3-7573e14a6cf4
# ╟─0d364111-35f0-4a3f-9c7a-1fcc426102bb
# ╟─d3eff505-597f-478d-8e01-5c2f6f365f19
# ╟─2c6110b5-cc07-4bf8-9515-7d7585e7b25a
# ╟─2703165f-1c4a-4ae0-8723-4fa00b54560b
# ╟─0f282d75-f954-4611-ab7e-d62d7d588e10
# ╟─5b9def3d-79e3-4c11-8f9f-27ed9f9a2e82
# ╟─b0165c9e-2828-418d-8265-298803cd721f
# ╠═462fc69d-e11e-4f54-8182-ad42ff6321ef
# ╟─74de5a7c-aa7c-44ee-b752-82f232f5c00a
# ╟─b5fc3fb7-a5ce-4107-9cf2-4b43999b6c9a
# ╟─21b1c011-8039-4347-9c80-101ad2fdb07a
# ╟─1dab9389-aec3-49d6-bcc2-57bbe367f84b
# ╠═a90fbfa5-49f4-4308-b977-bb271d1c76d9
# ╟─29e3b0ff-5a4a-418c-95c0-6d9861abbb31
# ╟─49146722-3a28-4610-a625-e7b06d013e87
# ╟─d691c2ba-6759-4b53-901d-165ff14fad01
# ╟─8a394b92-8c1b-4653-b2c4-1722a257ec60
# ╟─188a3b33-3bb7-486e-a232-06bd50bf778d
# ╟─5c2dabb8-2d77-4afa-bd0f-233f636d863e
# ╟─4cab9347-40b6-421b-8ee3-c5306d60fffe
# ╟─dee1520c-80eb-4611-9cb0-ad93cb32f9aa
# ╟─794b9a20-951f-4ea1-8c74-018cefb47ae4
# ╟─79ae96ec-ede3-4ca8-ba1c-2139e1d4b2c5
# ╠═01121bc2-60c6-4cba-9c68-e99def46736c
# ╟─ea157415-146c-4755-90b4-7060cdff3302
# ╟─599bcd25-e939-4e3f-9d03-cdf27badeef9
# ╟─4804a758-7548-42a6-9c21-5e82668d82d8
# ╟─3d85d992-dec5-4984-ab61-2fc71fe592c9
# ╠═7c248a7c-22ae-4752-9ae9-6d66b42dc1da
# ╟─b699f6a1-f065-4a2c-b20d-bdf571102e98
# ╟─0e91e96c-5ffd-4972-a1dd-c96524da29a4
# ╟─b867bac9-3984-49c4-9afa-de729db21f7b
# ╟─ff80f7b7-d5e6-473a-8626-190d9f3eedca
# ╟─7515f748-d9ad-4444-afba-43dcabcc165c
# ╟─4b6ea471-822d-4333-af3b-600a5cdac961
# ╟─aee25458-9214-4b0f-ab00-ae6060efb2a5
# ╟─6fad3bd0-76e9-48c3-87fc-4139d083e65d
# ╟─75eb8374-d24b-4ab3-9412-6a9255b7090e
# ╟─2d21efe6-773d-4255-b04b-be20dd609942
# ╟─02d28a01-4ade-47c5-9381-0d9806f5399d
# ╟─a3152eb3-9b54-4d1e-ae39-4c51a922af59
# ╟─54770db9-45ed-46d9-ab07-c78aff4d3b23
# ╠═2db53eb3-b9c4-4fc1-8b87-2189b3c835a3
# ╟─8f4391d4-8f1e-4e30-9229-196a30b4f0d7
# ╟─391c4da7-82b8-4e92-81db-9e56bd32f8a5
# ╟─d59a90ca-f666-4cfd-b8c5-ea81e56d22d3
# ╟─319cc440-fa4c-4d09-bfa3-aec77d6a79bf
# ╠═e0de149f-0e42-41fc-9173-53f85eae9acd
# ╟─b815249b-b1ee-42fb-a7db-144ae967eab8
# ╟─125190f1-0731-4757-95af-0a64b891b294
# ╟─9a6ce4d5-0095-4a42-a18c-059360750f53
# ╟─d254202d-162c-43ff-96c6-c6774e843f5c
# ╟─aad8f02d-3843-4187-b007-381e535b4c42
# ╟─305a180b-e65b-44ed-8754-00e3b32a139a
# ╟─8737cb3e-f390-4dbd-ab53-1c1671ea36e4
# ╟─745aa527-46cb-4777-b996-0e4b54633afd
# ╟─b2ed57d7-a199-48a1-9951-1d0cdef6d076
# ╟─a24dd754-63b8-4fcc-8a76-38e085c8dc96
# ╟─f243e139-380e-49cb-955d-7a580c8a7b33
# ╠═370fa37f-368c-49e1-9b3f-65f76993b09c
# ╟─d099d0df-ae74-4ceb-af53-305e6029e5d2
# ╟─d7f72084-8fcb-4d1f-994b-b0e31ecaccff
# ╟─3747cc76-bbe2-48a9-bc38-69529af96100
# ╟─93e9755a-5919-4d01-8f69-e2f363088c2d
# ╠═994c886e-e03a-47f5-b04b-647f446567b3
# ╟─4cfec33d-3abd-42bf-927a-0059c6767e14
# ╟─ba5eba06-c713-4eee-88d6-be9a7a222043
# ╟─01d77131-86c4-4ffc-bb05-0f3550b09409
# ╟─645f3a10-0af1-4573-b5ec-d698ff5d12a2
# ╟─eb41a25c-7e74-4144-aaaf-6ba05e3097f0
# ╟─5daa5df3-8aea-4f1a-ba38-753a7f7b6f9b
# ╟─4f7530c9-4744-4c67-91f5-e77948053e3b
# ╟─200eaef4-f500-4d78-9725-fe6183563943
# ╟─2337aaa1-820c-469e-85d8-0de47c26223b
# ╟─e8ff12e9-d637-4f81-ae2d-cd4120d3b127
# ╠═9e17f7b5-95b6-4945-99e7-965bc65ce9b6
# ╟─cbee50c5-d31b-4538-9211-ab2d622a7373
# ╟─da129dbf-bca4-4ae8-b20a-a4232e8538e2
# ╟─feab2f26-13b2-4aca-a2ad-8adc7c5eafcc
# ╟─a1257066-754e-4a38-a711-ec5af82207c0
# ╠═045d40c6-deab-45c8-8995-8407d4832276
# ╟─315a6845-7b3c-49c5-add7-42c88b063ab7
# ╟─687f00b7-5fae-402a-9475-2e74c4a28256
# ╟─26241cd6-619c-4e8e-8182-55b0d4360238
# ╟─2a712fb7-673a-4c9a-9aff-f5f882a4761c
# ╟─0e6a2398-b99f-487d-ac8c-ef28e38b7739
# ╟─b3772b0e-3022-46ed-8ea7-a4fa64c01a43
# ╟─086cd178-b2c0-4dd8-b4a8-5f22ad0718d0
# ╟─1fbc185b-e7af-4a49-b2d5-a582bc3f658e
# ╟─42976fef-fb3c-4cb4-94a9-8cac5db8fce9
# ╟─9c06873d-60f9-4850-b789-ec85030bec12
# ╟─061737cd-0d64-4e74-800d-47cddb3c6368
# ╠═7460ae37-8214-4a6e-bbec-a7016baa0320
# ╟─0b19a92d-d4c0-4735-a953-7e227ce51fd7
# ╟─99f0881a-636c-4ada-8a9a-679719b7bbe7
# ╟─d6ebfb27-e7d7-43ca-8ac9-864e5e43373e
# ╟─4670de0d-ed22-4e12-b238-3d7176d74b08
# ╠═c8382cfb-7046-471a-9a79-2b7022697deb
# ╟─55352d1c-1307-41ad-810b-da857f340f39
# ╟─58919029-abbc-4c50-8d54-8bcddb71a848
# ╟─8447a741-d1ec-4117-89b7-6999e748381b
# ╟─5d0280b4-da6d-40e9-b5cb-b26dab1b273c
# ╟─3cce1304-c75e-484c-bee7-20ab999709d0
# ╟─3eda40e9-31f2-48e5-9be4-88a960a0c65f
# ╟─41e33e79-7c01-4b65-9e79-718ddb1645e2
# ╟─fa5dabf7-8f0b-4f4d-afe7-9635d6cc1769
# ╟─e21273e1-e6e9-45c2-b99e-1dda86c386c9
# ╟─55eb1b17-0bb6-4774-92cf-03e2c67aba0c
# ╟─67ba42cd-04cc-49d4-a490-d56662d80e03
# ╠═616768f1-e52f-4f26-beea-8f085b812786
# ╟─d27853ed-346f-427e-9ed9-8683620a6114
# ╟─b479a379-e7de-46a4-a8b4-3e0ce1464424
# ╟─32d541bf-0b7a-4eb8-b17f-597e04c1aaad
# ╟─225815b8-9b2e-4a7f-814f-89f6adf2e8f6
# ╠═a034a1a7-3af2-4f5d-9897-3e77a67f7568
# ╟─91d659a4-5286-4476-86f5-f805268c88fe
# ╟─8457bb11-126e-401a-ae53-560381ce2ef4
# ╟─abe3a1f2-2ab4-456b-8dfb-f809b3f3ef05
# ╟─e2612b27-cc9b-4de7-bce0-1a5224d45bb4
# ╟─fc529c7f-24d1-4fd9-a72d-1a757b6d18da
# ╟─123f5380-ca96-40e8-bf03-9386e7f165fa
# ╟─a40a4249-bf07-4642-ab41-407315cc4e45
# ╟─dbcc9b76-b38d-4b86-a356-7fc85daa04ec
# ╟─6fee85b2-61e0-4f09-9c77-88981cd0be7a
# ╟─f4c13240-4562-4ca3-bffd-d9010b83deea
# ╠═6a38ddfa-b2ac-45fb-af00-baf09220f5bc
# ╟─3675225d-299a-4181-a3e9-d86e64e470b2
# ╟─78f3a968-71ac-484c-9347-f0d47cc5de5a
# ╟─96e155ef-1159-4474-9737-798b57a62591
# ╟─951a9d02-3321-46d9-97c3-f90eda11bb55
# ╠═757ef5e3-8ea3-4a44-84ea-c8035e619cb8
# ╟─546283f0-7bde-40dc-a040-0d468b50bb84
# ╟─c259b3e5-3786-41e4-a801-585e09dd0944
# ╟─53580e32-14df-4392-b3ce-5aaa6ccf7e0a
# ╟─c3bb759b-ef5e-4569-be23-a664a775e714
# ╟─1c18d0c0-d515-420c-8f10-4be6d5abfaeb
# ╟─b28683c5-f7a5-497e-9558-b15f27d69c8a
# ╟─c7bba1dc-6b77-4250-bfd6-8051914d7cbb
# ╟─be3186c4-c5c9-45f5-9050-4856b70bb2e8
# ╟─21580099-c127-46c7-917f-84edbf446b80
# ╟─d51c4084-8255-420c-84d2-c0a908f8e7ed
# ╟─55a0d112-5c54-4bbf-931c-0fa930c0196c
# ╟─d4bcc27a-d2f5-4e58-b6e4-29ea7b0d4eeb
# ╟─76acfab3-3e3d-4413-9b57-7ff763b32084
# ╟─0d07bf65-32b8-4569-9539-9f3124ed2322
# ╟─d3f70173-3ae7-4d46-b383-b703b15a0483
# ╟─9009028e-4bd2-46c9-a47c-c795070ee18f
# ╟─0f699a12-b72a-4074-8771-95f650b68a25
# ╟─8c8e1e48-3656-4bdb-9759-c88a55e08cda
# ╟─d22a515c-fc00-428d-9f36-846295efc210
# ╟─0c0f8f82-179c-433b-b570-41100a74bd37
# ╟─50519776-3063-49f3-a7fb-a3949ac1ee99
# ╟─458760be-de62-44f2-af1f-a642d831e2e1
# ╟─a645728d-3627-4f5f-b01e-321d1db730cb
# ╟─52387ede-862c-452a-9afa-01eb9cbc6943
# ╟─b4401e0b-7df8-44b5-a992-aa703f9a4145
# ╟─53a52203-2882-4c6b-81fe-6d4e4a698efe
# ╟─1724968f-2a3f-48f2-9be0-b829a8e20f41
# ╟─04c2cf76-bd83-400e-b86e-ec44589249a0
# ╟─f63d0a79-667f-4f2c-9d49-fc0acfc0250a
# ╟─f5a4937d-d4e1-4ec0-9eb9-f22052c48ed1
# ╟─b81d1973-cf1c-41d5-9e60-d4ac0f2ecf21
# ╟─59ded85c-0ef8-418d-abc9-28b0ddfd8afd
# ╟─e1c064db-c1f5-4f83-8b00-cd2a48021068
