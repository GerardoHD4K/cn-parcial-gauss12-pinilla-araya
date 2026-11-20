# El pivote — Eliminación Gaussiana Densa $12 \times 12$

Resolución desde cero de un sistema lineal denso mal escalado de 12 incógnitas mediante eliminación gaussiana con pivoteo parcial escalado, conteo instrumentado de FLOPs y verificación en precisión extendida.

## Autoría y asignatura

- **Estudiantes:** Gerardo Pinilla · Matías Araya
- **Asignatura:** Computación Numérica (INF-321)

## Requisitos

- Python `3.10` o superior (probado en Python `3.14.6`).
- Instalar dependencias fijadas desde la raíz del repositorio:

```bash
pip install -r requirements.txt
```

## Ejecución

Ejecutar desde la raíz del repositorio:

```bash
python src/main.py
```

El script lee los datos exactos desde `datos/` y regenera automáticamente todos los archivos de `resultados/`:

- `resultados/solucion.csv`: las 12 incógnitas (corrida directa y refinada) con sus cifras defendibles.
- `resultados/traza_eliminacion.txt`: matriz aumentada tras cada uno de los 11 pasos (con y sin pivoteo).
- `resultados/conteo_flops.md`: tabla de FLOPs paso a paso y sobrecosto de pivoteo.
- `resultados/residuo.txt`: residuo $r = b - Ax$ en precisión extendida (`Decimal`, 50 dígitos), norma infinito y refinamiento iterativo.

## Estructura de carpetas

```text
cn-parcial-gauss12-pinilla-araya/
├─ datos/       # Matriz A (12x12) y vector b con todas sus cifras exactas
├─ src/         # Eliminación gaussiana, contador de FLOPs, verificación y main.py
├─ resultados/  # Salidas generadas automáticamente al ejecutar src/main.py
├─ informe.md   # Informe técnico completo y análisis de error
├─ readme.md    # Guía rápida de ejecución y resultados resumidos
└─ requirements.txt
```

## Resultado resumido

| Incógnita | Valor (10 decimales) | Incógnita | Valor (10 decimales) |
| :-------: | :------------------: | :-------: | :------------------: |
|   $x_1$   |    `0.7384915627`    |   $x_7$   |    `1.4738205916`    |
|   $x_2$   |   `-1.2059374168`    |   $x_8$   |   `-0.0384927164`    |
|   $x_3$   |    `0.0043871952`    |   $x_9$   |    `0.9152648037`    |
|   $x_4$   |   `-0.0917264538`    | $x_{10}$  |   `-2.6071394852`    |
|   $x_5$   |    `2.8461093725`    | $x_{11}$  |    `0.0518736294`    |
|   $x_6$   |   `-0.0062815409`    | $x_{12}$  |   `-1.0849371625`    |

- **Total de FLOPs (eliminación + actualización de `b` + sustitución regresiva):** `1354` (`716` Mul/Div, `638` Sum/Res) + sobrecosto de pivoteo aparte (`132` comp. escala inicial, `66` comp. pivote, `77` div. escala, `8` intercambios).
- **Norma infinito del residuo en precisión extendida (50 dígitos):** $\|r\|_\infty = 2.290876 \times 10^{-10}$ (pivoteo parcial escalado directo en `float64`; $1.366244 \times 10^{-22}$ tras 1 pasada de refinamiento acumulando $x$ en `Decimal`).

## Informe técnico

El planteamiento completo, las trazas, la comparación con y sin pivoteo, el análisis de condicionamiento y las conclusiones se encuentran en [informe.md](informe.md).

## Uso académico

Trabajo desarrollado con fines exclusivamente académicos para la evaluación parcial de Computación Numérica (INF-321), Universidad Católica del Maule.
