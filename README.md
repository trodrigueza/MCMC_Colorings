Para información sobre la implementación ver la carpeta `notebook` allí se encuentra el notebook de [Pluto](https://plutojl.org) (También un pdf del mismo).

Consideraciones:

- Las implementaciones en Julia se encuentran en el notebook.
- [Sitio web de la clase](https://sites.google.com/unal.edu.co/fohernandezr/docencia/materias/cadenas-de-markov-y-aplicaciones).

Para las q-coloraciones:

- En la carpeta `cpp/q-colorings` se encuentran dos implementaciones realizadas en C++, una utilizando arrays (`main.cpp`) y una utilizando vectores (`vector.cpp`). `est.cpp` es un pequeño programa para las estimaciones realizadas con el promedio de los ratios ($k=14..20$).

- En la carpeta `results/q-colorings` se encuentran todos los resultados obtenidos en formato csv.

- Para el cálculo de resultados se utilizó la implementación realizada con arrays, que es un poco más rápida.

Para hard-core:

- En la carpeta `cpp/hard-core` se encuentra la implementación realizada en C++ (`main.cpp`).

- Adicionalmente, en `cpp/hard-core` se encuentra la implementación utilizando programación dinámica para calcular los valores reales (`dp.cpp`).

- En la carpeta `results/hard-core` se encuentran todos los resultados obtenidos en formato csv.
