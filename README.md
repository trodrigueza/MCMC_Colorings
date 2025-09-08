For information about the implementations, see the `notebook` folder. There you will find the [Pluto](https://plutojl.org) notebook (also available as a PDF).

Considerations:

- The Julia implementations are included in the notebook.
- The images in the `images` folder were generated using Julia (and are also in the notebook).
- [Course website](https://sites.google.com/unal.edu.co/fohernandezr/docencia/materias/cadenas-de-markov-y-aplicaciones).

For q-colorings:

- In the `cpp/q-colorings` folder there are two implementations in C++: one using arrays (`main.cpp`) and one using vectors (`vector.cpp`). The file `est.cpp` is a small program for the estimates made using the average of the ratios ($k=14..20$).

- In the `results/q-colorings` folder you will find all results in CSV format.

- In the `mathematica` folder there are exact calculations up to $k=7$. The other results were taken from [here](https://oeis.org/wiki/Colorings_of_grid_graphs).

- For the estimates, we mostly used the array-based implementation, which is slightly faster.

For hard-core:

- In the `cpp/hard-core` folder you will find the C++ implementation (`main.cpp`).

- Additionally, in `cpp/hard-core` there is an implementation using dynamic programming to calculate the exact values (`dp.cpp`).

- In the `results/hard-core` folder you will find all results in CSV format.
