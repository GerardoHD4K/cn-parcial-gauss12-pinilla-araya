import sys
sys.dont_write_bytecode = True

import csv
import platform
from decimal import Decimal
from pathlib import Path
import numpy as np

# Rutas del proyecto (funciona igual en Windows y Linux desde cualquier carpeta)
DIR_SRC = Path(__file__).resolve().parent
if str(DIR_SRC) not in sys.path:
    sys.path.insert(0, str(DIR_SRC))

from conteo import guardar_conteo_md, obtener_totales
from gauss import eliminacion_gaussiana, guardar_traza_txt
from verificacion import (
    analizar_condicionamiento,
    calcular_residuo_extendido,
    control_numpy,
    estimar_digitos_y_bits,
    guardar_residuo_txt,
    guardar_solucion_csv,
    refinamiento_iterativo,
)

RAIZ = DIR_SRC.parent
DIR_DATOS = RAIZ / "datos"
DIR_RESULTADOS = RAIZ / "resultados"


def cargar_sistema(ruta_a, ruta_b):
    # Leemos A y b como texto primero para no perder ninguna cifra en Decimal
    A_str = []
    with Path(ruta_a).open("r", encoding="utf-8", newline="") as f:
        for fila in csv.reader(f):
            celdas = [c.strip() for c in fila if c.strip()]
            if celdas:
                A_str.append(celdas)

    b_str = []
    with Path(ruta_b).open("r", encoding="utf-8", newline="") as f:
        for fila in csv.reader(f):
            celdas = [c.strip() for c in fila if c.strip()]
            if celdas:
                b_str.append(celdas[0])

    if len(A_str) != 12 or any(len(f) != 12 for f in A_str) or len(b_str) != 12:
        raise ValueError("Se espera una matriz 12x12 y un vector b de largo 12 en datos/.")

    A_dec = [[Decimal(v) for v in fila] for fila in A_str]
    b_dec = [Decimal(v) for v in b_str]
    A_f64 = [[float(v) for v in fila] for fila in A_str]
    b_f64 = [float(v) for v in b_str]

    return A_f64, b_f64, A_dec, b_dec


def main():
    # Asegurar carpeta de salida
    DIR_RESULTADOS.mkdir(parents=True, exist_ok=True)

    # 1. Cargar matriz A y vector b sin truncar
    A_f64, b_f64, A_dec, b_dec = cargar_sistema(
        DIR_DATOS / "matriz_a.csv",
        DIR_DATOS / "vector_b.csv",
    )
    A_np = np.array(A_f64, dtype=np.float64)
    b_np = np.array(b_f64, dtype=np.float64)

    # 2. Ejecutar las corridas sobre el mismo código de eliminación gaussiana
    res_sin = eliminacion_gaussiana(A_f64, b_f64, pivoteo="ninguno")
    res_par = eliminacion_gaussiana(A_f64, b_f64, pivoteo="parcial")
    res_esc = eliminacion_gaussiana(A_f64, b_f64, pivoteo="escalado")

    # 3. Residuos en precisión extendida (Decimal 50 dígitos)
    r_sin = calcular_residuo_extendido(A_dec, b_dec, res_sin["x"])
    r_par = calcular_residuo_extendido(A_dec, b_dec, res_par["x"])
    r_esc = calcular_residuo_extendido(A_dec, b_dec, res_esc["x"])

    # 4. Refinamiento iterativo reutilizando L y U
    ref_esc = refinamiento_iterativo(A_dec, b_dec, res_esc, iteraciones=2)
    ref_sin = refinamiento_iterativo(A_dec, b_dec, res_sin, iteraciones=2)

    # 5. Dígitos defendibles y bits perdidos a partir de las correcciones dx
    metricas_esc_dir = estimar_digitos_y_bits(res_esc["x"], ref_esc["pasos_f64"][0]["dx"])
    metricas_esc_ref = estimar_digitos_y_bits(ref_esc["x_final_f64"], ref_esc["pasos_f64"][1]["dx"])
    metricas_sin_dir = estimar_digitos_y_bits(res_sin["x"], ref_sin["pasos_f64"][0]["dx"])

    # 6. Condicionamiento y control con numpy.linalg.solve
    cond = analizar_condicionamiento(A_np)
    ctrl_np = control_numpy(A_np, b_np, A_dec, b_dec, res_esc["x"], ref_esc["x_final_f64"])

    # 7. Guardar los 4 archivos de resultados/
    guardar_traza_txt(A_f64, b_f64, res_esc, res_sin, DIR_RESULTADOS / "traza_eliminacion.txt")
    guardar_conteo_md(
        {
            "ninguno": res_sin["contador"],
            "parcial": res_par["contador"],
            "escalado": res_esc["contador"],
        },
        DIR_RESULTADOS / "conteo_flops.md",
    )
    guardar_solucion_csv(metricas_esc_dir, metricas_esc_ref, DIR_RESULTADOS / "solucion.csv")
    guardar_residuo_txt(
        r_sin, r_par, r_esc, ref_esc, ref_sin, cond, ctrl_np,
        DIR_RESULTADOS / "residuo.txt",
    )

    # 8. Resumen en consola
    tot = obtener_totales(res_esc["contador"])
    digs_sin = [m["digitos"] for m in metricas_sin_dir]
    digs_esc = [m["digitos"] for m in metricas_esc_dir]

    print("=" * 76)
    print("EL PIVOTE — Eliminación Gaussiana 12x12 (INF-321, UCM)")
    print(f"Python {platform.python_version()} | NumPy {np.__version__}")
    print("=" * 76)
    print(f"Ceros bajo diagonal  : {res_esc['ceros_bajo_diagonal']}/{res_esc['total_bajo_diagonal']} exactos")
    print(f"FLOPs aritméticos    : {tot['total_flops']} (Mul/Div={tot['total_mul_div']}, Sum/Res={tot['total_sum_res']})")
    print(
        f"Sobrecosto pivoteo   : {tot['comp_escala_inicial']} comp. escala + "
        f"{tot['comp_pivote']} comp. pivote + {tot['div_pivote']} div. escala + {tot['swaps']} swaps"
    )
    print("-" * 76)
    print(
        f"Sin pivoteo      -> max|m|={res_sin['max_mult']:.3e}, ρ={res_sin['factor_crecimiento']:.3e}, "
        f"||r||_inf={r_sin['norma_inf']:.3e}, digs={min(digs_sin)}..{max(digs_sin)}"
    )
    print(
        f"Pivoteo escalado -> max|m|={res_esc['max_mult']:.3e} (esc={res_esc['max_mult_escalado']:.3f}), "
        f"ρ={res_esc['factor_crecimiento']:.3f}, ||r||_inf={r_esc['norma_inf']:.3e}, digs={min(digs_esc)}..{max(digs_esc)}"
    )
    print(
        f"Tras 1 refinam.  -> ||dx1||_inf={ref_esc['pasos_f64'][0]['norma_dx']:.3e}, "
        f"||dx2||_inf={ref_esc['pasos_f64'][1]['norma_dx']:.3e}, digs=15 en todas"
    )
    print("-" * 76)
    for i, (m_dir, m_ref) in enumerate(zip(metricas_esc_dir, metricas_esc_ref), start=1):
        print(
            f"  x_{i:<2d} = {m_dir['formateado']:>18s} ({m_dir['digitos']:2d} d.s.)  |  "
            f"refinado: {m_ref['formateado']:>20s} ({m_ref['digitos']:2d} d.s.)"
        )
    print("-" * 76)
    print(f"Resultados generados en: {DIR_RESULTADOS}")


if __name__ == "__main__":
    main()
