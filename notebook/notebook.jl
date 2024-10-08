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

1. Comenzamos con $G_0$, el grafo que tiene todos los vértices pero ninguna arista. Esto significa que cualquier asignación de colores es válida, es decir $Z_0 = q^{k^2}$.

2. Para cada arista $e_i = \{x_i, y_i\}$ que se añade al grafo, estimaremos la proporción $\frac{Z_i}{Z_{i-1}}$. Para esto, podemos notar que las $q$-coloraciones de $G_i$ son aquellas $q$-coloraciones de $G_{i-1}$ en las que los vértices de la arista añadida en $G_i$ tienen un color distinto. Es decir, $\frac{Z_i}{Z_{i-1}}$ es la probabilidad de que dos vértices conectados por la nueva arista tengan diferentes colores. Esto es:

$\frac{Z_i}{Z_{i-1}} = P_{G_{i-1}}(X(x_i) \neq X(y_i))$

   donde $X$ es una $q$-coloración aleatoria elegida uniformemente entre las $q$-coloraciones de $G_{i-1}$.

3. Para estimar $\frac{Z_i}{Z_{i-1}}$, realizamos $N$ simulaciones utilizando el algoritmo de Gibbs sampler:

   a) Partimos de una $q$-coloración inicial aleatoria.

   b) Realizamos $T$ pasos de Gibbs sampler:
      - Elejimos un vértice $v\in V$ uniformemente al azar (Esto lo adaptaremos para el Systematic Gibbs Sampler que será descrito más abajo).
      - Determinamos los colores disponibles para $v$.
      - Elegimos uniformemente un color entre los disponibles y lo asignamos a $v$.

   c) Después de los $T$ pasos, contamos si $x_i$ y $y_i$ tienen colores diferentes.
   
   d) Estimamos $\frac{Z_i}{Z_{i-1}}$ como la proporción de las $N$ simulaciones donde $x_i$ y $y_i$ tienen colores diferentes.

4. Construimos la estimación final multiplicando sucesivamente estas razones, partiendo desde $Z_0$ (que como vimos es trivial de calcular).

$${Z}_m = Z_0 \prod_{i=1}^m {\frac{Z_i}{Z_{i-1}}}$$
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
@bind q Slider(1:1:15, default=10)

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
Finalmente, establecemos el número de simulaciones y el número de pasos para el systematic Gibbs sampler, para esto utilizamos los resultados del teorema mencionado anteriormente. Como observación, según dichos resultados se deberían hacer $\frac{48d^3n^3}{\epsilon^2}$ simulaciones, sin embargo realizaremos $\frac{n^3}{\epsilon^2}$ (manteniendo el orden de $n^3$ sobre el número de simulaciones) pues de la otra manera para el ejemplo actual de $k=5$ y $q=10$ habían pasado más de diez mil segundos (2.7 horas) y aún no se obtenía respuesta.
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

**Nota:** En el repositorio de [GitHub](https://github.com/trodrigueza/MCMC_Colorings) en la carpeta */Mathematica*, se encuentra un notebook en el que se realizan los cálculos que se mostrarán en el presente informe, en caso de que un resultado sea obtenido mediante *Mathematica* (o de una fuente exterior) escribiremos $\star$.
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
Es decir que la estimación obtenida tiene un error relativo $\left(\frac{|x-\hat{x}|}{x}\cdot 100\right)$ de:
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
¿Existe algún resultado sobre esto?
"""

# ╔═╡ d4d9a98a-01ac-478c-a787-37c15dc669e6
md"""
---
#### Implementación (C++):
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

**Nota:** En la carpeta *cpp/q-colorings* del [GitHub](https://github.com/trodrigueza/MCMC_Colorings) también se encuentra la implementación y detalles sobre la misma.
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
md"""
---
### Reporte de resultados:
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

- El algoritmo proporciona un control razonable sobre el error, lo que permite confiar en la precisión de las estimaciones para la toma de decisiones.

- La variación en el número de simulaciones y los pasos del sampler muestra que el algoritmo es sensible a estos parámetros, lo que subraya la importancia de una configuración cuidadosa para maximizar la eficiencia y la efectividad.

- Existe espacio para optimizar aún más el algoritmo ajustando dinámicamente los parámetros en respuesta a las características del grafo y los resultados intermedios de las simulaciones, incluso se podría hacer uso de computación paralela para implementar varios hilos al momento de realizar las estimaciones de los ratios (esta probablemente sería una optimización clave y como anotación, un compañero nuestro ajeno a la clase lo intentará hacer).
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
Buscando en internet no encontramos resultados sobre las respuestas exactas a este problema, es por esto que implementamos una algoritmo utilizando [programación dinámica](https://es.wikipedia.org/wiki/Programación_dinámica) (dp) logrando encontrar respuestas hasta $k=25$ en un tiempo relativamente moderado. Al final del notebook explicamos este enfoque. Las respuestas obtenidas utilizando fueron:
"""

# ╔═╡ a645728d-3627-4f5f-b01e-321d1db730cb
file_path_reshc = "../results/hard-core/dpAnswers/out.csv"

# ╔═╡ 52387ede-862c-452a-9afa-01eb9cbc6943
dfreshc = CSV.File(file_path_reshc, delim=' ', header=true) |> DataFrame

# ╔═╡ 5e67e81c-638f-449f-af07-2487021240e6
md"""
Para estimar el número de configuraciones factibles del modelo Hard Core podemos adaptar el algoritmo anterior de conteo de $q$-coloraciones: El modelo Hard-Core puede ser visto como un problema de $2$-coloración, donde el $0$ representa un vértice sin ocupar. La restricción de que dos vértices adyacentes no tengan ambos el valor de $1$ es análogo a la restricción del problema de $q$-coloraciones en donde dos vértices adyacentes no pueden tener el mismo color.

**Descripción formal** (Análoga a la del anterior ejercicio):
"""

# ╔═╡ 86bbaa49-8a57-4e4f-abfd-454eb57c0463
md"""
Sea $G = (V, E)$ el grafo reticular $k\times k$ con $V = \{v_1, v_2, ..., v_{n=k^2}\}$ y $E = \{e_1, e_2, ..., e_{m=2k(k-1)}\}$. Definimos:

$$G_0 = (V, \emptyset)$$
$$G_i = (V, \{e_1, ..., e_i\})$$
para $1\leq i \leq m$. $Z_i$ será el número de configuraciones factibles del modelo hard-core en el grafo $G_i$.

Queremos aproximar $Z_{m}$, el cual puede ser reescrito como:

$$Z_{m} = \frac{Z_{m}}{Z_{m-1}} \times \frac{Z_{m-1}}{Z_{m-2}} \times ... \times \frac{Z_2}{Z_1}\times \frac{Z_1}{Z_0}\times Z_0.$$

1. Comenzamos con $G_0$, el grafo que tiene todos los vértices pero ninguna arista. En este caso, cualquier asignación de $0$s y $1$s es válida, es decir $Z_0 = 2^{k^2}$.

2. Para cada arista $e_i = \{x_i, y_i\}$ que se añade al grafo, estimamos la proporción $\frac{Z_i}{Z_{i-1}}$. Las configuraciones factibles de $G_i$ son aquellas configuraciones factibles de $G_{i-1}$ en las que los vértices de la arista añadida en $G_i$ no tienen ambos el valor 1. Por lo tanto:

   $$\frac{Z_i}{Z_{i-1}} = P_{G_{i-1}}(X(x_i) = 0 \text{ o } X(y_i) = 0)$$

   donde $X$ es una configuración aleatoria elegida uniformemente entre las configuraciones factibles de $G_{i-1}$.

3. Para estimar $\frac{Z_i}{Z_{i-1}}$, realizamos múltiples simulaciones utilizando el algoritmo MCMC descrito en el Ejemplo $7.2$ del libro *Finite Markov chains and algorithm applications* (adaptado a **Systematic Gibbs Sampler**):

   a) Partimos de una configuración factible inicial (por ejemplo, todos $0$s).

   b) En cada paso del MCMC:
      - Elegimos un vértice $v \in V$ uniformemente al azar (Esto cambia en la implementación con Systematic Gibbs Sampler).
      - Lanzamos una moneda justa.
      - Si sale cara y todos los vecinos de $v$ tienen valor $0$, asignamos $X_{n+1}(v) = 1$; de lo contrario, $X_{n+1}(v) = 0$.
      - Para todos los demás vértices $w \neq v$, mantenemos $X_{n+1}(w) = X_n(w)$.

   c) Después de un número suficiente de pasos para permitir que la cadena se mezcle, contamos la proporción de configuraciones en las que $x_i = 0$ o $y_i = 0$.

4. Construimos la estimación final multiplicando sucesivamente estas razones, partiendo desde $Z_0 = 2^{k^2}$:

   $$Z_m = Z_0 \prod_{i=1}^m {\frac{Z_i}{Z_{i-1}}}$$

"""

# ╔═╡ dc6d37bf-d8ce-4db2-98be-f69a2ca4338b
html"""
<h4>Implementación:</h4>
"""

# ╔═╡ c6313a94-4403-415a-ab1e-b5b649a284ed
md"""
En realidad no cambian muchas cosas.

**Parámetros:**
"""

# ╔═╡ 4900b650-4b0c-41ac-968d-f450a56a7422
@bind k1 Slider(1:1:20, default=6)

# ╔═╡ 33830c92-e0af-4457-ace6-d13410589517
@printf("k = %i", k1)

# ╔═╡ bf5db19f-4a01-489a-b080-ffb0a88e411b
@bind eps1 Slider(0.1:0.1:1, default=0.7)

# ╔═╡ 0c28834e-68a5-4293-b58b-6395635a3bc4
@printf("epsilon = %0.1f", eps1)

# ╔═╡ d3e90ebf-f8fa-4d45-b958-3bf365f8bc1f
num_simulations_hc = Int(ceil((k1*k1)^3/eps1^2))

# ╔═╡ 6c1075bf-abb6-4f30-83fd-0a1b2fc1aea7
num_gibbs_steps_hc = Int(ceil(abs(k1*k1 * ((2 * log(k1*k1) + log(1 / eps1) + log(8)) / log((2) / 32) + 1))))

# ╔═╡ ef1c4b5d-9b79-4a1b-93ce-6487e84870df
md"""
Primero generamos un lattice con la configuración de todos los vértices con valor $0$ (factible).
"""

# ╔═╡ 9eef9473-3848-4588-b863-ee9b795e87eb
lattice_hc = [0 for i in 1:k1, j in 1:k1]

# ╔═╡ a7b406ff-3878-4024-b5ff-34067b89bbb2
md"""
Generamos las aristas:
"""

# ╔═╡ 4112e338-ddd2-47e5-8bb3-3128bd168d41
edges_hc = gen_edges(k1)

# ╔═╡ 2c971739-0ed9-4ea5-8393-594cfeaa7236
md"""
Implementamos el Systematic Gibbs Sampler ($k^2$ pasos):
"""

# ╔═╡ c8d140d0-d0f9-4e45-806f-9e54a0b60ca2
function gibbs_step_hc(k1::Int, lattice::Array{Int,2}, neighbours)
    for i in 1:k1, j in 1:k1
		f = 1
		if !isempty(neighbours[i,j])
        	for ne in neighbours[i,j]
            	if lattice[ne[1], ne[2]] == 1
					f = 0
					break
				end
        	end
		end

        if f == 1 && rand(0:1) == 1
            lattice[i, j] = 1
		else
			lattice[i, j] = 0
		end
    end
end

# ╔═╡ dd932974-95d2-432b-81b3-92bdc2ac5a0c
md"""
Función para estimar los ratios:
"""

# ╔═╡ 1a79db3a-c3dd-4940-892a-31f6117616f1
function estimate_ratio_hc(k1::Int, num_simulations::Int, num_gibbs_steps::Int, edge, lattice, neighbours)
    uy, ux = edge[1]
    vy, vx = edge[2]
    count = 0
    for _ in 1:num_simulations
        for _ in 1:ceil(num_gibbs_steps/(k1*k1))+1
            gibbs_step_hc(k1, lattice, neighbours)
        end
        if (lattice[uy, ux] == 0 || lattice[vy, vx] == 0)
			count += 1
		end
    end

    return BigFloat(count) / num_simulations
end

# ╔═╡ 52e7f7f5-59e8-4820-8bbc-ab1063181f13
md"""
Estimamos la respuesta:
"""

# ╔═╡ b1f31edd-a8bc-432b-89b8-b9c4287c7f68
begin
Random.seed!(1234)
Z_hc = BigFloat(2)^(k1*k1)
neighbours_hc = [Vector{Tuple{Int,Int}}() for _ in 1:k1, _ in 1:k1]
ratios_hc = []
for edge in edges_hc
	uy, ux = edge[1]
    vy, vx = edge[2]
	ratio_hc = estimate_ratio_hc(k1, num_simulations_hc, num_gibbs_steps_hc, edge, lattice_hc, neighbours_hc)
    Z_hc *= BigFloat(ratio_hc)
	push!(ratios_hc, BigFloat(ratio_hc))	
	push!(neighbours_hc[uy, ux], (vy, vx))
    push!(neighbours_hc[vy, vx], (uy, ux))
end
val_est_hc = round(Z_hc)
@printf("Número estimado de configuraciones factibles para el retículo %dx%d:\n", k1, k1)
@printf("%d\n", val_est_hc)
end

# ╔═╡ 19063be6-7607-4189-8677-3e75017a3efa
md"""
El número exacto de configuraciones factibles para el grafo reticular $6\times 6$ es $5598861$. Es decir que la estimación obtenida tiene un error relativo del:
"""

# ╔═╡ 408f4220-16b6-4ed4-be86-fff30f50e4ff
@printf("%0.3f%%", 100*abs(val_est_hc-5598861)/5598861)

# ╔═╡ e877328e-1f29-4e88-b720-a0ea9ac29b74
md"""
Luego fue una muy buena estimación.
"""

# ╔═╡ a604f419-a73c-41b1-8bc2-891c0ed3e582
html"""
<h3>Reporte de resultados:</h3>
"""

# ╔═╡ 57e8c338-544d-4bd9-81ef-e3440a6d3aad
md"""
Presentaremos los resultados obtenidos de la implementación en C++ que se encuentra en la carpeta `cpp/hard-core` en el [GitHub](https://github.com/trodrigueza/MCMC_Colorings).
"""

# ╔═╡ 21930e01-c3ae-4c7f-9b84-4b0cba242538
md"""
**Convención:**

- |$k$| Dimensión retículo.
- |$\text{Sims}$| Número de simulaciones.
- |$\text{Gibbs}$| Número de pasos del systematic Gibbs sampler.
- |$\text{Est}$| Estimación obtenida.
- |$\text{Res}$| Valor real.
- |$\text{R\_err}$| Error relativo (%).
- |$\text{Mean\_r}$| Promedio de las razones obtenidas.
- |$\text{Time}$| Tiempo en segundos que tardó el programa.
"""

# ╔═╡ 00c5adbc-9e30-4638-adfd-3cea76c0fdde
md"""
Para los siguientes resultados se utilizó:

$\text{num\_sims} = \frac{n^3}{\epsilon^2}$

$\text{gibbs\_steps} = n \left(\frac{\log{n} + \log{\epsilon^{-1}}}{\log{\frac{2}{4d^2}}}\right)$
"""

# ╔═╡ 8b9e5883-8986-45c1-8d04-d565de087a1a
md"""
$k=2...9$
"""

# ╔═╡ 0f2a76c8-61e9-474a-b4ab-1ad9463d67d1
md"""
$\epsilon = 0.2$
"""

# ╔═╡ 73b6ae2c-edb4-4d6a-a28d-d353a94d325e
file_path_hca = "../results/hard-core/estimations/k=2..9/eps=0.2.csv"

# ╔═╡ 28e4c64d-154e-4e84-96cd-80f48fcbec26
dfhca = CSV.File(file_path_hca, delim=' ', header=true) |> DataFrame

# ╔═╡ 1cb3d3a4-fc56-41e5-a39d-e20a02de0288
md"""
$\epsilon = 0.7$
"""

# ╔═╡ 188f8098-ff39-4885-b6c3-ce8bad162c15
file_path_hcb = "../results/hard-core/estimations/k=2..9/eps=0.7.csv"

# ╔═╡ a03ffcc8-7e85-4208-a56b-c0fbed144487
dfhcb = CSV.File(file_path_hcb, delim=' ', header=true) |> DataFrame

# ╔═╡ 8ef9401b-c354-4ac9-b6b9-a386ed3c5612
md"""
$\epsilon = 1$
"""

# ╔═╡ 44dd1186-6afe-4d89-b73b-843170c44a27
file_path_hcc = "../results/hard-core/estimations/k=2..9/eps=1.csv"

# ╔═╡ 16ee1499-8048-4603-ada0-f8bac2b1fd0c
dfhcc = CSV.File(file_path_hcc, delim=' ', header=true) |> DataFrame

# ╔═╡ cf2bd14f-7e0e-4eff-af92-6cde28108a58
begin
phca1 = plot(dfreshc.k[1:8], dfreshc.Res[1:8], label="Respuesta Real", linewidth=0.5, marker=:star, color=:red, markersize=3, size=(550, 350), yscale=:log10, leg=:bottomright, ylimits = (1e0, 1e18), title="Estimación vs Respuesta Real", titlefontsize=10, labelfontsize=8, legendfontsize=6);
	
plot!(dfhca.k[1:8], dfhca.Est[1:8], label="Estimación eps = 0.2", xlabel="k", ylabel="Valor", linewidth=0.7, marker=:rect, color=:blue, alpha=0.5, markersize=4);

plot!(dfhcb.k[1:8], dfhcb.Est[1:8], label="Estimación eps = 0.7", xlabel="k", ylabel="Valor", linewidth=0.7, marker=:hexagon, color=:green, alpha=0.5, markersize=4);

plot!(dfhcc.k[1:8], dfhcc.Est[1:8], label="Estimación eps = 1", xlabel="k", ylabel="Valor", linewidth=0.7, marker=:diamond, color=:orange, alpha=0.5, markersize=5);

phca2 = plot(dfhca.k[1:8], dfhca.R_err[1:8], label="eps = 0.2", color=:blue, xlabel="k", ylabel="Error relativo (%)", marker=:rect, size=(550, 350), title="Comparación errores relativos", titlefontsize=10, labelfontsize=8, legendfontsize=6, ylimits=(0, 8))

plot!(dfhcb.k[1:8], dfhcb.R_err[1:8], label="eps = 0.7", xlabel="k", linewidth=0.7, marker=:hexagon, color=:green, alpha=0.5, markersize=4);
	
plot!(dfhcc.k[1:8], dfhcc.R_err[1:8], label="eps = 1", color=:orange, marker=:diamond, secondary=true)

phca3 = bar(dfhca.k[1:8], dfhca.Time[1:8], label="eps = 0.2", color=:blue, xlabel="q", ylabel="Tiempo(seg)", size=(550, 350), bar_width=0.5, title="Comparación tiempos", titlefontsize=10, labelfontsize=8, legendfontsize=6, ylimits = (0, 500), leg=:topleft)

bar!(dfhcb.k[1:8], dfhcb.Time[1:8], label="eps = 0.7", color=:green, ylabel="Tiempo(seg)", size=(550, 350), bar_width=0.55)
	
bar!(dfhcc.k[1:8], dfhcc.Time[1:8], label="eps = 1", color=:orange, secondary=true, bar_width=0.6)
	
bar!(dfreshc.k[1:8], dfreshc.Time[1:8], label="dp", color=:red, secondary=true, bar_width=0.6)
	
plot(phca1, phca3, phca2, layout=layout, size=(700, 500))
end

# ╔═╡ e8999de2-9d93-437e-806a-994daa9928db
md"""
Para los siguientes resultados fijamos el número de simulación a $1e6$.
"""

# ╔═╡ 86372603-e49f-45a2-a798-8956dbb87051
md"""
$k=10...25$
"""

# ╔═╡ 3b39aadd-3ea2-4bed-848e-d661b2e6a141
md"""
$\epsilon=0.7$
"""

# ╔═╡ 98833e3a-7041-46d6-9553-cffd8fec5a0a
md"""
$\text{num\_sims} = 1000000$
"""

# ╔═╡ 0e5e9587-02a9-4034-a142-1291c6396792
file_path_hcd = "../results/hard-core/estimations/k=10..25/eps=0.7.csv"

# ╔═╡ 604ce5c3-ac1b-4236-ae33-a2d14c4a35f6
dfhcd = CSV.File(file_path_hcd, delim=' ', header=true) |> DataFrame

# ╔═╡ d06b5dfc-1713-44a3-bb7c-b2195b8582c4
md"""
$\epsilon=1$
"""

# ╔═╡ 73f27435-bde2-4a39-82d5-d45717a59e21
md"""
$\text{num\_sims} = 1000000$
"""

# ╔═╡ 9bf8014a-9de0-407d-a591-978a7daed43f
file_path_hce = "../results/hard-core/estimations/k=10..25/eps=1.csv"

# ╔═╡ bd10a9a0-5a58-48ab-8cd2-eb9a900d0873
dfhce = CSV.File(file_path_hce, delim=' ', header=true) |> DataFrame

# ╔═╡ 2211c46e-c4ba-4b02-96a8-2af0624f87b8
begin
phcb1 = plot(dfreshc.k[9:24], dfreshc.Res[9:24], label="Respuesta Real", linewidth=0.5, marker=:star, color=:red, markersize=3, size=(550, 350), yscale=:log10, leg=:bottomright, ylimits = (1e15, 1e115), title="Estimación vs Respuesta Real", titlefontsize=10, labelfontsize=8, legendfontsize=6);
	
plot!(dfhcd.k[1:16], dfhcd.Est[1:16], label="Estimación eps = 0.7", xlabel="k", ylabel="Valor", linewidth=0.7, marker=:hexagon, color=:green, alpha=0.5, markersize=4);

plot!(dfhce.k[1:16], dfhce.Est[1:16], label="Estimación eps = 1", xlabel="k", ylabel="Valor", linewidth=0.7, marker=:diamond, color=:orange, alpha=0.5, markersize=5);

phcb2 = plot(dfhcd.k[1:16], dfhcd.R_err[1:16], label="eps = 0.7", xlabel="k", linewidth=0.7, marker=:hexagon, color=:green, alpha=0.5, markersize=4, ylabel="Error relativo (%)", size=(550, 350), title="Comparación errores relativos", titlefontsize=10, labelfontsize=8, legendfontsize=6, ylimits=(0, 2.3));
	
plot!(dfhce.k[1:16], dfhce.R_err[1:16], label="eps = 1", color=:orange, marker=:diamond, secondary=true)

phcb3 = bar(dfhcd.k[1:16], dfhcd.Time[1:16], xlabel="k", ylabel="Tiempo(seg)", label="eps = 0.7", color=:green, size=(550, 350), bar_width=0.5, title="Comparación tiempos", titlefontsize=10, labelfontsize=8, legendfontsize=6)
	
bar!(dfhce.k[1:16], dfhce.Time[1:16], label="eps = 1", color=:orange, secondary=true, bar_width=0.6)
	
bar!(dfreshc.k[9:24], dfreshc.Time[9:24], label="dp", color=:red, secondary=true, bar_width=0.6)
	
plot(phcb1, phcb3, phcb2, layout=layout, size=(600, 400))
end

# ╔═╡ 97a7625c-6d23-4cf3-b2b1-408680751942
md"""
#### Promedio ratios:
"""

# ╔═╡ 1093ca76-cd23-4688-bd8c-80bc9d18b045
md"""
Considerando $5\leq k \leq 25$:
"""

# ╔═╡ d7ca2f97-0c8c-463c-b898-bf61317f7284
begin
plhcr1 = plot(dfhca.k[4:8], dfhca.Mean_r[4:8], label="eps = 0.2", color=:blue, ylabel="Promedios", marker=:circle, size=(600, 400), alpha = 0.5, leg=:right, titlefontsize=10, labelfontsize=8, legendfontsize=7)

plot!(dfhcd.k, dfhcd.Mean_r, label="eps = 0.7", marker=:utriangle, color=:green, alpha = 0.5)

title!("Comparación promedios de los ratios")
xlabel!("k")

data = vcat(dfhca.Mean_r[4:8], dfhcd.Mean_r)
plhcr2 = boxplot(["Datos"], data, title="Promedio ratios - Box Plot", ylabel="Promedios", legend=false, size=(400,300), titlefontsize=10, labelfontsize=8, legendfontsize=7, color = :orange, alpha=0.5)

layout3 = @layout [a; b]
plot(plhcr1, plhcr2, layout = layout3, size = (600, 400))
end

# ╔═╡ 0521270d-6318-41ca-b01b-e0773544acc0
md"""
Podemos observar que los valores son muy consistentes, sugiriendo una cierta independecia al valor de $k$, como lo notamos en el ejercicio de las $q$-coloraciones.
"""

# ╔═╡ 36df72a6-2b28-4ecb-b11d-3d297500a8b2
md"""
### Conclusiones:
"""

# ╔═╡ 3fefdf33-231f-4d58-a6e3-5c171487b00e
md"""
- Aunque la solución con programación dinámica (dp) parece ser más rápida, eventualmente para $k$s grandes esto dejará de ser así pues la complejidad de la dp es exponencial, abajo hablamos de esto.

- Fijando el número de simulaciones -poniéndole el tope de $1e6$- obtenemos resultados excelentes para $\epsilon\approx 1$, con seguridad podríamos intentar bajar este tope según los parámetros para obtener respuestas confiables en tiempos mucho más rápidos.

- Es un excelente método para estimar las respuestas al problema de contar las configuraciones factibles del modelo Hard Core de un grafo reticular.
"""

# ╔═╡ b4401e0b-7df8-44b5-a992-aa703f9a4145
md"""
---
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

En cuanto al número de estados, cada fila se representa con una `bitmask` de $k$ bits. Hay $2^k$ posibles configuraciones (máscaras) para cada fila. En cuanto a las transiciones, en cada fila se consideran todas las posibles configuraciones y se verifica su compatibilídad con todas las configuraciones de la fila anterior. Así, para cada una de las $k$ filas la transición involucra comparar $2^k\times 2^k = 4^k$ pares de máscaras.

Es decir que la complejidad es de $O(k\times 4^k)$, permitiéndonos calcular en un tiempo razonable máximo hasta $k=20$.

Podemos optimizar de la siguiente manera:

- Hay $2^k$ posibles configuraciones para cada fila. Pero hay $l<2^k$ posibles configuraciones factibles para cada fila. Podemos precalcularlas y guardarlas ($O(2^k)$).

- Para cada máscara válida, podemos determinar las máscaras que son compatibles con ella (verticalmente). Esto también puede ser precalculado, guardando en una lista las máscaras compatibles para cada máscara válida ($O(l^2)$).

- Para las transiciones, utilizamos la lista precalculada.

Es decir que la complejidad bajaría de $O(k\times4^k)$ a $O(2^k + kl^2)$, -que no es mucho pero es trabajo honesto-.

Observación: $l$ resulta ser $F(k+2)$ donde $F(i)$ es el $i$-ésimo número de Fibonacci. Teniendo en cuenta que $F(k+2) \approx \frac{\phi^{k+2}}{\sqrt{5}}$, entonces podemos calcular en tiempo razonable (menos de $20$ minutos) hasta $k=25$.
"""

# ╔═╡ b81d1973-cf1c-41d5-9e60-d4ac0f2ecf21
md"""
#### Implementación:
"""

# ╔═╡ 59ded85c-0ef8-418d-abc9-28b0ddfd8afd
md"""
Por facilidad realizamos la implementación en C++, esta se encuentra en el [GitHub](https://github.com/trodrigueza/MCMC_Colorings) en *cpp/hard-core/dp.cpp*.
"""

# ╔═╡ 8d9c27ec-5331-4240-9942-3c21c4773207
md"""
## Bibliografía:
- Häggström O. Finite Markov Chains and Algorithmic Applications. Cambridge University Press; 2002.
- Colorings of grid graphs - OeisWiki. (n.d.). https://oeis.org/wiki/Colorings_of_grid_graphs
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"

[compat]
CSV = "~0.10.14"
DataFrames = "~1.6.1"
Plots = "~1.40.5"
PlutoUI = "~0.7.59"
StatsPlots = "~0.15.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.4"
manifest_format = "2.0"
project_hash = "d432fc5ea84416649ea4cc1c0577c2c3e273cc8f"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "6c834533dc1fabd820c1db03c839bf97e45a3fab"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.14"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "71acdbf594aab5bbb2cec89b208c41b4c411e49f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.24.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "9ebb045901e9bbf58767a9f34ff89831ed711aae"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.7"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "66c4c81f259586e8f002eacebc177e1fb06363b0"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.11"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "0e0a1264b0942f1f3abb2b30891f2a590cc652ac"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.110"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "9f00e42f8d99fdde64d40c8ea5d14269a2e2c1aa"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.21"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fd0002c0b5362d7eb952450ad5eb742443340d6e"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.12.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "532f9126ad901533af1d4f5c198867227a7bb077"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+1"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "629693584cef594c3f6f99e76e7a7ad17e60e8d5"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a8863b69c2a0859f2c2c87ebdc4c6712e88bdf0d"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.7+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "7c4195be1649ae622304031ed46a2f4df989f1eb"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.24"

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

[[deps.InlineStrings]]
git-tree-sha1 = "45521d31238e87ee9f9732561bfee12d4eebd52d"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.2"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14eb2b542e748570b56446f4c50fbfb2306ebc45"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "39d64b09147620f5ffbf6b2d3255be3c901bec63"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.8"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "7d703202e65efa1369de1279c162b915e245eed1"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.9"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e16271d212accd09d52ee0ae98956b8a05c4b626"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "17.0.6+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

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

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MultivariateStats]]
deps = ["Arpack", "Distributions", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "816620e3aac93e5b5359e4fdaf23ca4525b00ddf"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "91a67b4d73842da90b526011fa85c5c4c9343fe0"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.18"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "1a27764e945a152f7ca7efa04de513d473e9542e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.1"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a028ee3cb5641cccc4c24e90c36b0a4f7707bdf5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9dd97171646850ee607593965ce1f55063d8d3f9"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "082f0c4b70c202c37784ce4bfbc33c9f437685bf"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.5"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

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

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "66b20dd35966a748321d3b2537c4584cf40387c7"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PtrArrays]]
git-tree-sha1 = "f011fbb92c4d401059b2212c05c0601b70f8b759"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "e237232771fdafbae3db5c31275303e056afaa9f"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.10.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e60724fd3beea548353984dc61c943ecddb0e29a"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.3+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "ff11acffdb082493657550959d4feb4b6149e73a"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.5"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "cef0472124fab0695b58ca35a77c6fb942fdab8a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.1"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "NaNMath", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "3b1dcbf62e469a67f6733ae493401e53d92ff543"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.15.7"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "e84b3a11b9bece70d14cce63406bbc79ed3464d2"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.2"

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

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "fcdae142c1cfc7d89de2d11e08721d0f2f86c98a"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.6"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "1165b0443d0eca63ac1e32b8c0eb69ed2f4f8127"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "936081b536ae4aa65415d869287d43ef3cb576b2"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.53.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
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
# ╟─62b3fe4f-6546-4b5b-8187-69573f0e2f45
# ╟─e7b1d579-3300-49cd-8e0d-28f942f99be9
# ╟─a9ba504b-3ecc-4e80-b063-06242abb23ca
# ╠═6fdef38f-a45a-4cfe-9694-a0ac144f16b9
# ╟─fdc76661-9c2e-47f0-8716-a3d719e88998
# ╟─aa182161-2906-4e5d-aa08-3153c19f97e0
# ╠═57afe92f-88a4-4fdc-83db-5d58cc525e73
# ╟─205399df-28be-4674-a1c2-35bc5db6e765
# ╟─62e7c451-ba95-496c-8f16-5d7deea3388a
# ╠═899c7974-3284-4eb8-bf7e-7d5a3028ec04
# ╟─52042204-a585-46b1-9757-694b9f1de30c
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
# ╟─5e67e81c-638f-449f-af07-2487021240e6
# ╟─86bbaa49-8a57-4e4f-abfd-454eb57c0463
# ╟─dc6d37bf-d8ce-4db2-98be-f69a2ca4338b
# ╟─c6313a94-4403-415a-ab1e-b5b649a284ed
# ╟─33830c92-e0af-4457-ace6-d13410589517
# ╟─4900b650-4b0c-41ac-968d-f450a56a7422
# ╟─0c28834e-68a5-4293-b58b-6395635a3bc4
# ╟─bf5db19f-4a01-489a-b080-ffb0a88e411b
# ╠═d3e90ebf-f8fa-4d45-b958-3bf365f8bc1f
# ╠═6c1075bf-abb6-4f30-83fd-0a1b2fc1aea7
# ╟─ef1c4b5d-9b79-4a1b-93ce-6487e84870df
# ╠═9eef9473-3848-4588-b863-ee9b795e87eb
# ╟─a7b406ff-3878-4024-b5ff-34067b89bbb2
# ╠═4112e338-ddd2-47e5-8bb3-3128bd168d41
# ╟─2c971739-0ed9-4ea5-8393-594cfeaa7236
# ╠═c8d140d0-d0f9-4e45-806f-9e54a0b60ca2
# ╟─dd932974-95d2-432b-81b3-92bdc2ac5a0c
# ╠═1a79db3a-c3dd-4940-892a-31f6117616f1
# ╟─52e7f7f5-59e8-4820-8bbc-ab1063181f13
# ╠═b1f31edd-a8bc-432b-89b8-b9c4287c7f68
# ╟─19063be6-7607-4189-8677-3e75017a3efa
# ╠═408f4220-16b6-4ed4-be86-fff30f50e4ff
# ╟─e877328e-1f29-4e88-b720-a0ea9ac29b74
# ╟─a604f419-a73c-41b1-8bc2-891c0ed3e582
# ╟─57e8c338-544d-4bd9-81ef-e3440a6d3aad
# ╟─21930e01-c3ae-4c7f-9b84-4b0cba242538
# ╟─00c5adbc-9e30-4638-adfd-3cea76c0fdde
# ╟─8b9e5883-8986-45c1-8d04-d565de087a1a
# ╟─0f2a76c8-61e9-474a-b4ab-1ad9463d67d1
# ╟─73b6ae2c-edb4-4d6a-a28d-d353a94d325e
# ╟─28e4c64d-154e-4e84-96cd-80f48fcbec26
# ╟─1cb3d3a4-fc56-41e5-a39d-e20a02de0288
# ╟─188f8098-ff39-4885-b6c3-ce8bad162c15
# ╟─a03ffcc8-7e85-4208-a56b-c0fbed144487
# ╟─8ef9401b-c354-4ac9-b6b9-a386ed3c5612
# ╟─44dd1186-6afe-4d89-b73b-843170c44a27
# ╟─16ee1499-8048-4603-ada0-f8bac2b1fd0c
# ╟─cf2bd14f-7e0e-4eff-af92-6cde28108a58
# ╟─e8999de2-9d93-437e-806a-994daa9928db
# ╟─86372603-e49f-45a2-a798-8956dbb87051
# ╟─3b39aadd-3ea2-4bed-848e-d661b2e6a141
# ╟─98833e3a-7041-46d6-9553-cffd8fec5a0a
# ╟─0e5e9587-02a9-4034-a142-1291c6396792
# ╟─604ce5c3-ac1b-4236-ae33-a2d14c4a35f6
# ╟─d06b5dfc-1713-44a3-bb7c-b2195b8582c4
# ╟─73f27435-bde2-4a39-82d5-d45717a59e21
# ╟─9bf8014a-9de0-407d-a591-978a7daed43f
# ╟─bd10a9a0-5a58-48ab-8cd2-eb9a900d0873
# ╟─2211c46e-c4ba-4b02-96a8-2af0624f87b8
# ╟─97a7625c-6d23-4cf3-b2b1-408680751942
# ╟─1093ca76-cd23-4688-bd8c-80bc9d18b045
# ╟─d7ca2f97-0c8c-463c-b898-bf61317f7284
# ╟─0521270d-6318-41ca-b01b-e0773544acc0
# ╟─36df72a6-2b28-4ecb-b11d-3d297500a8b2
# ╟─3fefdf33-231f-4d58-a6e3-5c171487b00e
# ╟─b4401e0b-7df8-44b5-a992-aa703f9a4145
# ╟─53a52203-2882-4c6b-81fe-6d4e4a698efe
# ╟─1724968f-2a3f-48f2-9be0-b829a8e20f41
# ╟─04c2cf76-bd83-400e-b86e-ec44589249a0
# ╟─f63d0a79-667f-4f2c-9d49-fc0acfc0250a
# ╟─f5a4937d-d4e1-4ec0-9eb9-f22052c48ed1
# ╟─b81d1973-cf1c-41d5-9e60-d4ac0f2ecf21
# ╟─59ded85c-0ef8-418d-abc9-28b0ddfd8afd
# ╟─8d9c27ec-5331-4240-9942-3c21c4773207
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
