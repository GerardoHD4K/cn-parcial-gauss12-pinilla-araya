# Informe técnico — El pivote: eliminación gaussiana con mínimo costo y máxima exactitud

---

## 1. Identificación

| Campo | Detalle |
| :--- | :--- |
| **Asignatura** | Computación Numérica (INF-321) |
| **Carrera** | Ingeniería Civil Informática, Universidad Católica del Maule |
| **Docentes** | Juan Manuel Cabezas Contardo · Pedro Pesante |
| **Estudiantes** | Gerardo Pinilla · Matías Araya |
| **Fecha** | Talca, noviembre de 2026 |
| **Repositorio** | https://github.com/GerardoHD4K/cn-parcial-gauss12-pinilla-araya |

---

## 2. Planteamiento

Tenemos que resolver un sistema de 12 ecuaciones con 12 incógnitas, **Ax = b**, usando nuestra propia eliminación gaussiana en doble precisión (`float64`, que guarda 53 bits de mantisa y tiene un error de redondeo de unos 2.22×10⁻¹⁶). La matriz A no tiene ningún cero (144 valores con 9 cifras cada uno) y b trae 20 cifras. Los cargamos completos desde `datos/`, sin cortar ninguna cifra.

Como no hay ceros, no podemos ahorrar operaciones saltándonos entradas. El ahorro solo puede venir de no repetir trabajo. Además, la matriz tiene tres problemas puestos a propósito para que se pierdan dígitos:

1. **Las filas tienen tamaños muy distintos.** El mayor valor de cada fila va desde 8.05×10⁻⁵ (fila 5) hasta 9.67×10⁵ (fila 6), es decir, una diferencia de unas 10¹⁰ veces.

   | Tamaño de la fila | Filas |
   | :--- | :--- |
   | Del orden de 10⁻⁵ | F1, F5, F9 |
   | Del orden de 10⁻³ a 10⁻² | F3, F7, F11 |
   | Del orden de 10³ | F4, F10 |
   | Del orden de 10⁴ a 10⁵ | F2, F6, F8, F12 |

2. **El primer pivote es diminuto.** a₁₁ = 1.3407×10⁻⁸ es mucho más chico que el resto de su fila y 13 órdenes de magnitud menor que a₆₁ = 3.69×10⁵. No es cero, así que se puede usar, pero se pierden muchos dígitos.
3. **La fila 7 es casi una combinación de las filas 3 y 5.** Se diferencia de ella solo en el orden de 10⁻⁴. La matriz no es singular (det(A) ≈ 1.365×10⁹), pero queda muy mal condicionada: κ₂(A) ≈ 3.085×10¹³.

El problema es lograr el menor costo posible sin perder dígitos por culpa de estas tres trampas.

---

## 3. Estrategia de eliminación

Hicimos una sola función, `eliminacion_gaussiana(A, b, pivoteo)`, que se puede correr de tres formas: sin pivoteo, con pivoteo parcial, y con **pivoteo parcial escalado**, que es la que elegimos como principal. Guarda además L, U y las filas que se intercambiaron.

**Cómo funciona el pivoteo escalado.** Antes de empezar se calcula, para cada fila, su valor más grande en tamaño (sᵢ). En cada paso k se mira la columna k y se elige como pivote la fila con mayor |uᵢₖ| / sᵢ. Así se compara cada número contra el tamaño de su propia fila y no contra el de toda la matriz. Si la fila elegida no es la que toca, se intercambian.

**Por qué no pasamos del mínimo de operaciones:**
- **Un multiplicador por fila.** Calculamos mᵢₖ = uᵢₖ / uₖₖ una sola vez (1 división) y lo reutilizamos en toda la fila y en b.
- **No tocamos lo que ya es cero.** En el paso k el barrido empieza en la columna k+1 y dejamos `U[i][k] = 0.0` directamente. Eso evita 2×66 = 132 operaciones inútiles y deja exactamente 0 en los 66 lugares bajo la diagonal.
- **Una sola factorización.** Nunca calculamos la inversa. El refinamiento reutiliza L y U.

El orden de pivotes que resultó fue: F8, F2, F1, F6, F9, F11, F12, F5, F4, F10, F3, F7. Se ve que F1 se descarta en el paso 1 por su a₁₁ chiquito, pero entra en el paso 3, donde sí es un buen pivote para su fila (|u₃₃|/s₁ = 0.951). Y la fila casi dependiente (F7) queda al final. La traza de los 11 pasos está en `resultados/traza_eliminacion.txt` y confirma 66 de 66 ceros bajo la diagonal.

---

## 4. Conteo de operaciones

**Convención:** 1 flop es una operación elemental (+, −, ×, ÷), contando la actualización de b y la sustitución regresiva. El conteo sale de un contador dentro del código (`src/conteo.py`), no de una fórmula escrita a mano. El costo de buscar pivotes (comparaciones, divisiones por escala e intercambios) se cuenta aparte.

En el paso k hay r = 12 − k filas por debajo del pivote, y se hacen r divisiones, r² multiplicaciones y r² restas en U, más r multiplicaciones y r restas en b.

| Paso | Fila pivote | Mul/Div | Sum/Res | Total |
| :---: | :---: | :---: | :---: | :---: |
| 1 | F8 | 143 | 132 | **275** |
| 2 | F2 | 120 | 110 | **230** |
| 3 | F1 | 99 | 90 | **189** |
| 4 | F6 | 80 | 72 | **152** |
| 5 | F9 | 63 | 56 | **119** |
| 6 | F11 | 48 | 42 | **90** |
| 7 | F12 | 35 | 30 | **65** |
| 8 | F5 | 24 | 20 | **44** |
| 9 | F4 | 15 | 12 | **27** |
| 10 | F10 | 8 | 6 | **14** |
| 11 | F3 | 3 | 2 | **5** |
| **Eliminación** | — | **638** | **572** | **1210** |
| **Sustitución regresiva** | — | **78** | **66** | **144** |
| **TOTAL** | — | **716** | **638** | **1354** |

**Costo extra del pivoteo escalado** (aparte del total): 132 comparaciones para las escalas iniciales, 66 comparaciones de pivote, 77 divisiones por escala y 8 intercambios de fila.

**Comparación con las referencias:**

| Estrategia | Flops | Diferencia vs 1486 |
| :--- | :---: | :---: |
| Referencia del enunciado (eliminación densa) | 1486 | 0 % |
| **Nuestro código (sin pivoteo, con parcial o con escalado)** | **1354** | **−8.88 %** |
| Mismo conteo sumando las 77 divisiones de escala | 1431 | −3.70 % |
| Multiplicadores recalculados por columna | ≈ 1750 | +17.77 % |
| Usando la inversa | ≥ 7000 | muy por encima |

Nuestro total da 132 menos que la referencia, justo las operaciones que no hacemos sobre la columna que ya quedó en cero. Está dentro del ±10 % pedido. El total es el mismo con y sin pivoteo, porque intercambiar filas no agrega operaciones. El reporte completo se genera en `resultados/conteo_flops.md`.

---

## 5. Análisis de mantisa y pivoteo

Corrimos el mismo código sin pivoteo y con pivoteo parcial escalado. Agregamos el pivoteo parcial común (sin escalar) como control, para ver si escalar de verdad ayuda.

| Medida | Sin pivoteo | Parcial sin escalar | **Parcial escalado** |
| :--- | :---: | :---: | :---: |
| Multiplicador máximo | 2.753×10¹³ | 0.620 | 2.466×10⁹ |
| Multiplicador máximo "escalado" (m·sₖ/sᵢ) | 6.918×10³ | 8.225 | **0.962** |
| Factor de crecimiento (valor más grande durante la eliminación ÷ el de A) | 2.641×10³ | 1.000 | **1.783** |
| Norma infinito del residuo | 1.772×10⁻⁷ | 5.373×10⁻¹⁰ | **2.291×10⁻¹⁰** |
| Error relativo final | 2.080×10⁻¹⁰ | 5.520×10⁻¹³ | **2.357×10⁻¹³** |
| Dígitos correctos (mín. a máx.) | 7 a 10 | 9 a 12 | **9 a 13** |
| Bits de mantisa perdidos, de 53 (mín. a máx., promedio) | 19.2 a 29.3 (22.8) | 10.6 a 20.8 (14.3) | **9.4 a 19.6 (13.0)** |

**Por qué el pivote a₁₁ = 1.3407×10⁻⁸ es inaceptable.** Aunque no vale cero, es muy chico comparado con su propia fila (|a₁₁|/s₁ = 1.446×10⁻⁴). Si lo usamos para eliminar la fila 6 (a₆₁ = 3.69×10⁵), el multiplicador sale m₆₁ = 2.753×10¹³. Al restar m₆₁·F1 a F6, por ejemplo en la columna 9, la cuenta es a₆₉ − m₆₁·a₁₉ ≈ −8.652×10⁵ − 2.553×10⁹. El segundo número es unas 2951 veces más grande que el primero. Para sumarlos hay que alinear los exponentes, y los bits bajos del número chico (unos log₂(2951) ≈ 11.5 bits) se salen de la mantisa y se pierden para siempre, ya en el primer paso. Esto pasa también en las otras filas grandes, y la matriz crece hasta 2.554×10⁹. Al final se pierden entre 19 y 29 bits de 53, más de la mitad de la mantisa en x₃ y x₆, y el error queda unas 882 veces peor que con pivoteo escalado.

**Por qué el multiplicador 2.466×10⁹ del pivoteo escalado no hace daño.** Aparece en el paso 3, cuando F1 (tamaño 10⁻⁵) elimina a F6 (tamaño 10⁵). Se ve enorme solo porque las dos filas tienen tamaños muy distintos. Al multiplicarlo por los valores de F1 (10⁻⁴), el producto es del orden de 10⁵, igual que los valores de F6. Por eso ya no hay números que se coman a otros: el factor de crecimiento queda en 1.783 y el multiplicador escalado en 0.236 (menor que 1). No se pierden bits por alineamiento.

**Por qué conviene escalar y no usar el pivoteo parcial común.** El parcial común mira solo el tamaño absoluto y elige primero las filas grandes (F6, F8, F12, F2…), sin darse cuenta de que algunos valores son chicos dentro de su propia fila. Eso hace que el crecimiento respecto de cada fila llegue a 14.77, contra 3.18 con escalado. El error y el residuo también salen más de dos veces peores (5.52×10⁻¹³ contra 2.36×10⁻¹³ y 5.37×10⁻¹⁰ contra 2.29×10⁻¹⁰).

---

## 6. Condicionamiento

El número de condición κ(A) indica cuánto puede amplificarse un error pequeño en los datos al resolver el sistema. Estos son los valores que calculó `src/verificacion.py`:

| Medida | Valor | Qué significa |
| :--- | :---: | :--- |
| det(A) | 1.364551×10⁹ | No es singular |
| κ₂(A), matriz original | 3.084744×10¹³ | Muy alto |
| κ∞(A), matriz original | 5.367827×10¹³ | Igual de alto |
| κ₂(D⁻¹A), filas escaladas | 8.132348×10³ | Mucho más bajo |
| κ∞(D⁻¹A), filas escaladas | 1.490302×10⁴ | Mucho más bajo |

(D⁻¹A es la matriz que resulta de dividir cada fila por su valor más grande.)

**Cuántos dígitos se pueden recuperar.** La cota clásica dice que el error relativo es a lo más κ₂(A)·ε ≈ 3.085×10¹³ × 2.22×10⁻¹⁶ ≈ 6.85×10⁻³. Con eso no se podrían asegurar ni 2 dígitos. Sin embargo, en la práctica obtuvimos entre 9 y 13 dígitos correctos (error relativo cercano a 2.4×10⁻¹³).

**Por qué la cota resulta tan pesimista aquí.** La mayor parte de ese κ alto viene de que las filas tienen tamaños muy distintos, no de que el sistema sea realmente difícil. Si dividimos cada ecuación por el valor más grande de su fila, la solución x es exactamente la misma, pero el número de condición baja a κ₂(D⁻¹A) ≈ 8.13×10³. Como la computadora redondea de forma relativa a cada número, el error real se parece más a esta segunda cota: 8.132×10³ × 2.22×10⁻¹⁶ ≈ 1.81×10⁻¹². Esa cota predice entre 10 y 12 dígitos recuperables, que es lo que vimos.

**La casi dependencia de la fila 7.** Aun con las filas escaladas, κ sigue siendo de miles y no cercano a 1. La causa es que F7 es casi una combinación de F3 y F5. Ajustando por mínimos cuadrados (`np.linalg.lstsq`) obtuvimos F7 ≈ 4.177250·F3 + 583.246740·F5, con una diferencia de solo 1.942×10⁻⁴ (el 0.228 % del tamaño de F7). En los valores singulares de la matriz escalada se nota como un salto fuerte entre σ₁₁ = 0.138 y σ₁₂ = 4.831×10⁻⁴. También se ve en el paso 11 de la traza, cuando F3 elimina a F7 con m = 4.2132:

u₁₂,₁₂ = −6.6232×10⁻³ − (4.2132)(−1.6093×10⁻³) = 1.5692×10⁻⁴

Aquí se restan dos números casi iguales y se pierden cerca de dos dígitos en el último pivote.

---

## 7. Resultados

Para no afirmar más dígitos de los que se pueden defender, usamos la corrección de la primera pasada de refinamiento (|Δx|) como estimación del error de la solución directa. Afirmamos un dígito solo si ese error es menor que media unidad en ese dígito. Después de una pasada de refinamiento, la segunda corrección es menor que 10⁻¹⁶ en términos relativos, así que podemos defender **15 dígitos** en las 12 incógnitas. El detalle está en `resultados/solucion.csv`.

| Incógnita | Valor directo | Dígitos correctos | Error relativo estimado | Bits perdidos | Valor refinado (15 dígitos) |
| :---: | :---: | :---: | :---: | :---: | :---: |
| x₁ | `0.738491562700` | **12** | 7.26×10⁻¹⁴ | 9.35 | `0.738491562700000` |
| x₂ | `-1.20593741680` | **12** | 5.56×10⁻¹³ | 12.29 | `-1.20593741680000` |
| x₃ | `0.004387195200` | **10** | 8.30×10⁻¹¹ | 19.51 | `0.00438719520000000` |
| x₄ | `-0.091726453800` | **11** | 2.55×10⁻¹² | 14.48 | `-0.0917264538000000` |
| x₅ | `2.846109372500` | **13** | 1.21×10⁻¹³ | 10.09 | `2.84610937250000` |
| x₆ | `-0.00628154090` | **9** | 8.54×10⁻¹¹ | 19.55 | `-0.00628154090000000` |
| x₇ | `1.473820591600` | **13** | 3.10×10⁻¹³ | 11.45 | `1.47382059160000` |
| x₈ | `-0.038492716400` | **11** | 1.63×10⁻¹² | 13.84 | `-0.0384927164000000` |
| x₉ | `0.915264803700` | **12** | 1.87×10⁻¹³ | 10.71 | `0.915264803700000` |
| x₁₀ | `-2.607139485200` | **13** | 1.40×10⁻¹³ | 10.30 | `-2.60713948520000` |
| x₁₁ | `0.051873629400` | **11** | 1.25×10⁻¹² | 13.45 | `0.0518736294000000` |
| x₁₂ | `-1.084937162500` | **13** | 3.13×10⁻¹³ | 11.46 | `-1.08493716250000` |

En la corrida directa, x₃ y x₆ conservan menos dígitos (10 y 9) porque son números chicos (cerca de 10⁻³). El error absoluto es parecido en todas las incógnitas (unos 10⁻¹³), pero al dividirlo por un valor pequeño el error relativo se hace mayor. Aun así, las 12 incógnitas tienen al menos 9 dígitos correctos de forma directa, y llegan a 15 con el refinamiento.

---

## 8. Verificación

Todo está en `src/verificacion.py` y se guarda en `resultados/residuo.txt`.

**1. Residuo en precisión extendida.** Calculamos r = b − Ax con `Decimal` de 50 dígitos y los datos exactos. Con pivoteo escalado, ‖r‖∞ = 2.290876×10⁻¹⁰ (y ‖r‖₂ = 2.369583×10⁻¹⁰), unas 773 veces menor que sin pivoteo (1.771752×10⁻⁷).

**2. Refinamiento iterativo.** Calculamos el residuo en precisión extendida, resolvemos con los mismos L y U (276 flops por pasada) y corregimos x.

| Pasada | Variante | ‖r‖∞ antes | ‖Δx‖∞ | ‖r‖∞ después |
| :---: | :--- | :---: | :---: | :---: |
| 1 | Escalado, x en `float64` | 2.291×10⁻¹⁰ | 6.709×10⁻¹³ | 2.603×10⁻¹⁰ |
| 2 | Escalado, x en `float64` | 2.603×10⁻¹⁰ | 1.712×10⁻¹⁶ | 2.603×10⁻¹⁰ |
| 1 | Escalado, x acumulado en `Decimal` | 2.291×10⁻¹⁰ | 6.709×10⁻¹³ | **1.366×10⁻²²** |
| 2 | Escalado, x acumulado en `Decimal` | 1.366×10⁻²² | 1.757×10⁻²⁵ | **3.928×10⁻³⁵** |
| 1 | Sin pivoteo, x en `float64` | 1.772×10⁻⁷ | 5.921×10⁻¹⁰ | 2.603×10⁻¹⁰ |

Nos llamó la atención que en `float64` el residuo no baja de 2.603×10⁻¹⁰. La explicación está en la pasada 2: la corrección es de 6×10⁻¹⁷ en términos relativos, menor que el error de la máquina, así que x ya llegó al límite de lo que cabe en 53 bits. Ese x tiene un redondeo de unos 5×10⁻¹⁷, y al multiplicarlo por los valores de la fila 6 (cerca de 10⁶) queda un residuo inevitable del orden de 10⁻¹⁰. Para comprobarlo, repetimos guardando x en `Decimal` (con los mismos L y U), y ahí el residuo bajó a 1.366×10⁻²² y luego a 3.928×10⁻³⁵. Es decir, el refinamiento sí mejora la solución; lo que no cabe es mostrar la mejora en un `float64`.

**3. Control con biblioteca.** Solo como comparación, `numpy.linalg.solve` da ‖r‖∞ = 7.756×10⁻¹¹. La mayor diferencia relativa entre nuestra solución directa y la de NumPy es 7.307×10⁻¹², y con nuestra solución refinada es 9.267×10⁻¹¹. Esa diferencia más grande viene del error propio de NumPy en x₆, no del nuestro.

---

## 9. Conclusiones

**Qué elegimos y por qué.** Elegimos eliminación gaussiana con pivoteo parcial escalado más una pasada de refinamiento iterativo sobre los mismos L y U. El costo de las operaciones queda en 1354 flops, el mínimo sin trabajo repetido y un 8.88 % bajo la referencia de 1486. Se hace mas en el pivoteo (198 comparaciones, 77 divisiones y 8 intercambios) pero a cambio evitamos el problema del pivote a₁₁: el error baja unas 882 veces respecto de no pivotear y 2.3 veces respecto del pivoteo parcial común. Con eso obtenemos de 9 a 13 dígitos correctos en la corrida directa y 15 después de 276 flops más de refinamiento.

**Qué haríamos distinto con n = 1000.** El costo crece como 2n³/3, o sea unos 6.7×10⁸ flops, casi medio millón de veces más que con n = 12. Cambiaríamos varias cosas:
- **No usar ciclos de Python puro.** Tardarían minutos. Usaríamos NumPy con operaciones por bloques o librerías en C/Fortran (BLAS/LAPACK).
- **Escalar las filas una sola vez al principio.** En vez de dividir por sᵢ en cada paso, dividiríamos A y b una vez y usaríamos pivoteo parcial normal.
- **Guardar L y U en una sola matriz** (unos 8 MB en `float64`) y no guardar una copia de la matriz en cada paso, que serían casi 8 GB.
- **Calcular el residuo en precisión extendida con algo más rápido que `Decimal`**, porque 10⁶ productos en `Decimal` serían muy lentos.
- **Si la matriz tuviera muchos ceros**, usar formatos para matrices dispersas o métodos iterativos como Gauss-Seidel o gradiente conjugado.

---

## 10. Reproducibilidad

- **Python:** 3.14.6 (sirve cualquier versión 3.10 o superior).
- **Librerías:** `numpy==2.5.2` (más los módulos estándar `csv`, `decimal`, `math`, `pathlib`).
- **Cómo correrlo**, desde la raíz del repositorio:

```bash
pip install -r requirements.txt
python src/main.py
```
ss
Con eso se regeneran, sin pasos manuales, los cuatro archivos de `resultados/`: `solucion.csv`, `traza_eliminacion.txt`, `conteo_flops.md` y `residuo.txt`.

### Uso de asistentes de IA

Declaracion:
- Gemini-neon 4: para creacion de la matriz para el csv, organizar de mejor manera la explicación del informe y ayuda para guiarnos en un error en el main. 


### Referencias

1. Burden, R. L., & Faires, J. D. (2011). *Análisis Numérico* (9.ª ed.). Cengage Learning.
2. Chapra, S. C., & Canale, R. P. (2015). *Métodos Numéricos para Ingenieros* (7.ª ed.). McGraw-Hill.
3. Golub, G. H., & Van Loan, C. F. (2013). *Matrix Computations* (4th ed.). Johns Hopkins University Press.
4. IEEE Computer Society. (2019). *IEEE Standard for Floating-Point Arithmetic (IEEE Std 754-2019)*.
