# Conteo de FLOPs — Eliminación Gaussiana (n = 12)

**Convención:** 1 FLOP = 1 operación elemental (+, -, *, /). Multiplicadores reutilizados por fila; sobrecosto de pivoteo aparte.

| Paso | Pivote | Mul/Div | Sum/Res | Total FLOPs | Comp. | Div. escala | Swaps |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 1 | F8 | 143 | 132 | 275 | 11 | 12 | 1 |
| 2 | F2 | 120 | 110 | 230 | 10 | 11 | 0 |
| 3 | F1 | 99 | 90 | 189 | 9 | 10 | 1 |
| 4 | F6 | 80 | 72 | 152 | 8 | 9 | 1 |
| 5 | F9 | 63 | 56 | 119 | 7 | 8 | 1 |
| 6 | F11 | 48 | 42 | 90 | 6 | 7 | 1 |
| 7 | F12 | 35 | 30 | 65 | 5 | 6 | 1 |
| 8 | F5 | 24 | 20 | 44 | 4 | 5 | 1 |
| 9 | F4 | 15 | 12 | 27 | 3 | 4 | 1 |
| 10 | F10 | 8 | 6 | 14 | 2 | 3 | 0 |
| 11 | F3 | 3 | 2 | 5 | 1 | 2 | 0 |
| **Eliminación** | — | 638 | 572 | 1210 | 66 | 77 | 8 |
| **Sust. regresiva** | — | 78 | 66 | 144 | 0 | 0 | 0 |
| **TOTAL** | — | **716** | **638** | **1354** | **66** | **77** | **8** |

- **Escalas iniciales:** 132 comparaciones (0 FLOPs).
- **Cota 2n³/3:** 1152 FLOPs | **Referencia enunciado:** 1486 FLOPs.
- **Desviación vs referencia:** -8.88 % (dentro de ±10 %).
