import csv
import math
from decimal import Decimal, getcontext
from pathlib import Path
import numpy as np
from gauss import resolver_lu

# Fijamos 50 dígitos decimales para evaluar residuos sin pérdida de precisión
getcontext().prec = 50
EPS_FLOAT64 = 2.0 ** -52  # ~2.22e-16 en IEEE 754 doble precisión


def calcular_residuo_extendido(A_dec, b_dec, x):
    n = len(b_dec)
    x_dec = [v if isinstance(v, Decimal) else Decimal.from_float(v) for v in x]

    r = []
    norma_inf = Decimal(0)
    suma_cuad = Decimal(0)

    for i in range(n):
        ax_i = Decimal(0)
        for j in range(n):
            ax_i += A_dec[i][j] * x_dec[j]
        ri = b_dec[i] - ax_i
        r.append(ri)
        if abs(ri) > norma_inf:
            norma_inf = abs(ri)
        suma_cuad += ri * ri

    return {
        "r": r,
        "norma_inf": norma_inf,
        "norma_2": suma_cuad.sqrt(),
    }


def refinamiento_iterativo(A_dec, b_dec, res_gauss, iteraciones=2):
    # Reutiliza L y U propios: r en Decimal (50 díg.), LU*dx = P*r en float64
    n = len(b_dec)
    x_act = res_gauss["x"][:]
    pasos = []

    for it in range(1, iteraciones + 1):
        r_antes = calcular_residuo_extendido(A_dec, b_dec, x_act)
        r_f64 = [float(v) for v in r_antes["r"]]
        dx, flops_tri = resolver_lu(res_gauss["L"], res_gauss["U"], res_gauss["permutacion"], r_f64)
        x_nuevo = [x_act[i] + dx[i] for i in range(n)]
        r_despues = calcular_residuo_extendido(A_dec, b_dec, x_nuevo)

        norm_dx = max(abs(v) for v in dx)
        norm_x = max(abs(v) for v in x_act)
        pasos.append({
            "iteracion": it,
            "norma_r_antes": r_antes["norma_inf"],
            "norma_dx": norm_dx,
            "rel_dx": norm_dx / norm_x,
            "dx": dx,
            "x_nuevo": x_nuevo,
            "norma_r_despues": r_despues["norma_inf"],
            "flops": flops_tri,
        })
        x_act = x_nuevo

    residuo_final = calcular_residuo_extendido(A_dec, b_dec, x_act)

    # Corrida auxiliar acumulando x en Decimal (misma factorización L, U en float64)
    pasos_decimal = []
    x_dec = [Decimal.from_float(v) for v in res_gauss["x"]]
    for it in range(1, iteraciones + 1):
        r_dec = calcular_residuo_extendido(A_dec, b_dec, x_dec)
        dx_f, _ = resolver_lu(
            res_gauss["L"], res_gauss["U"], res_gauss["permutacion"],
            [float(v) for v in r_dec["r"]],
        )
        x_dec = [x_dec[i] + Decimal.from_float(dx_f[i]) for i in range(n)]
        r_desp = calcular_residuo_extendido(A_dec, b_dec, x_dec)
        pasos_decimal.append((it, max(abs(v) for v in dx_f), r_desp["norma_inf"]))

    return {
        "pasos_f64": pasos,
        "x_final_f64": x_act,
        "residuo_final_f64": residuo_final,
        "pasos_decimal": pasos_decimal,
    }


def formatear_sig(valor, digitos):
    if digitos <= 0 or valor == 0.0:
        return f"{valor:.6e}"
    exp = math.floor(math.log10(abs(valor)))
    decimales = max(0, digitos - 1 - exp)
    return f"{valor:.{decimales}f}"


def estimar_digitos_y_bits(x_vals, errores_abs, tope=15):
    # Criterio de media unidad en el último dígito: |e_i| <= 0.5 * 10^(p - d + 1)
    resultados = []
    for x_i, e_i in zip(x_vals, errores_abs):
        abs_x = abs(x_i)
        abs_e = abs(e_i)
        rel_e = abs_e / abs_x if abs_x > 0 else 0.0

        if abs_e == 0.0 or rel_e < 0.5 * EPS_FLOAT64:
            digitos = tope
            bits_ok = 53.0
            bits_perdidos = 0.0
        else:
            p = math.floor(math.log10(abs_x))
            d_calc = math.floor(p + 1.0 - math.log10(2.0 * abs_e))
            digitos = max(0, min(tope, d_calc))
            bits_ok = min(53.0, -math.log2(rel_e))
            bits_perdidos = max(0.0, 53.0 - bits_ok)

        resultados.append({
            "valor": x_i,
            "error_abs": abs_e,
            "error_rel": rel_e,
            "digitos": digitos,
            "formateado": formatear_sig(x_i, digitos),
            "bits_conservados": bits_ok,
            "bits_perdidos": bits_perdidos,
        })
    return resultados


def analizar_condicionamiento(A_np):
    escalas = np.max(np.abs(A_np), axis=1)
    A_escalada = A_np / escalas[:, None]

    det_a = float(np.linalg.det(A_np))
    cond_2 = float(np.linalg.cond(A_np, 2))
    cond_inf = float(np.linalg.cond(A_np, np.inf))
    cond_2_esc = float(np.linalg.cond(A_escalada, 2))
    cond_inf_esc = float(np.linalg.cond(A_escalada, np.inf))

    s_a = np.linalg.svd(A_np, compute_uv=False)
    s_esc = np.linalg.svd(A_escalada, compute_uv=False)

    # Verificamos la casi dependencia lineal F7 ≈ α*F3 + β*F5 por mínimos cuadrados
    base = np.vstack([A_np[2], A_np[4]]).T
    coefs, _, _, _ = np.linalg.lstsq(base, A_np[6], rcond=None)
    diferencia_f7 = A_np[6] - (coefs[0] * A_np[2] + coefs[1] * A_np[4])
    norma_dif_f7 = float(np.linalg.norm(diferencia_f7, 2))
    rel_dif_f7 = norma_dif_f7 / float(np.linalg.norm(A_np[6], 2))

    return {
        "det_a": det_a,
        "cond_2": cond_2,
        "cond_inf": cond_inf,
        "cond_2_escalada": cond_2_esc,
        "cond_inf_escalada": cond_inf_esc,
        "escala_min": float(escalas.min()),
        "escala_max": float(escalas.max()),
        "razon_escalas": float(escalas.max() / escalas.min()),
        "sigma_max_a": float(s_a[0]),
        "sigma_min_a": float(s_a[-1]),
        "sigma_11_esc": float(s_esc[-2]),
        "sigma_12_esc": float(s_esc[-1]),
        "cota_clasica": cond_2 * EPS_FLOAT64,
        "cota_escalada": cond_2_esc * EPS_FLOAT64,
        "coef_f3": float(coefs[0]),
        "coef_f5": float(coefs[1]),
        "norma_dif_f7": norma_dif_f7,
        "rel_dif_f7": rel_dif_f7,
    }


def control_numpy(A_np, b_np, A_dec, b_dec, x_directo, x_refinado):
    # Solo como control externo (no sustituye la implementación propia)
    x_np = [float(v) for v in np.linalg.solve(A_np, b_np)]
    res_np = calcular_residuo_extendido(A_dec, b_dec, x_np)

    dif_dir = max(abs(a - b) / abs(b) for a, b in zip(x_directo, x_np))
    dif_ref = max(abs(a - b) / abs(b) for a, b in zip(x_refinado, x_np))
    return {
        "x_numpy": x_np,
        "residuo_numpy": res_np,
        "dif_rel_directo": dif_dir,
        "dif_rel_refinado": dif_ref,
    }


def guardar_solucion_csv(metricas_dir, metricas_ref, ruta_salida):
    ruta = Path(ruta_salida)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    with ruta.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "incognita",
            "valor_pivoteo_escalado_f64",
            "valor_reportado_directo",
            "digitos_significativos_directo",
            "error_relativo_est_directo",
            "bits_perdidos_directo",
            "valor_refinado_f64",
            "valor_reportado_refinado",
            "digitos_significativos_refinado",
        ])
        for i, (m_dir, m_ref) in enumerate(zip(metricas_dir, metricas_ref), start=1):
            writer.writerow([
                f"x{i}",
                f"{m_dir['valor']:.16e}",
                m_dir["formateado"],
                m_dir["digitos"],
                f"{m_dir['error_rel']:.3e}",
                f"{m_dir['bits_perdidos']:.2f}",
                f"{m_ref['valor']:.16e}",
                m_ref["formateado"],
                m_ref["digitos"],
            ])


def guardar_residuo_txt(r_sin, r_par, r_esc, ref_esc, ref_sin, cond, ctrl_np, ruta_salida):
    r_np = ctrl_np["residuo_numpy"]
    r_ref_f64 = ref_esc["residuo_final_f64"]

    lineas = [
        "=" * 90,
        "VERIFICACIÓN EN PRECISIÓN EXTENDIDA (Decimal, 50 DÍGITOS) Y REFINAMIENTO ITERATIVO",
        "Computación Numérica (INF-321) — UCM | Estudiantes: Gerardo Pinilla, Matías Araya",
        "=" * 90,
        "",
        "1. NORMAS DEL RESIDUO r = b - Ax (precisión extendida: 50 dígitos decimales)",
        "-" * 90,
        f"  Sin pivoteo (directo float64)           : ||r||_inf = {r_sin['norma_inf']:.6e} | ||r||_2 = {r_sin['norma_2']:.6e}",
        f"  Pivoteo parcial sin escalar (control)   : ||r||_inf = {r_par['norma_inf']:.6e} | ||r||_2 = {r_par['norma_2']:.6e}",
        f"  Pivoteo parcial escalado (directo f64)  : ||r||_inf = {r_esc['norma_inf']:.6e} | ||r||_2 = {r_esc['norma_2']:.6e}",
        f"  Pivoteo escalado + refinamiento (f64)   : ||r||_inf = {r_ref_f64['norma_inf']:.6e} | ||r||_2 = {r_ref_f64['norma_2']:.6e}",
        f"  Control biblioteca (numpy.linalg.solve) : ||r||_inf = {r_np['norma_inf']:.6e} | ||r||_2 = {r_np['norma_2']:.6e}",
        "",
        "2. VECTOR RESIDUO COMPONENTE A COMPONENTE",
        "-" * 90,
        "   i |       r_i (escalado directo)   |     r_i (escalado refinado)    |     r_i (sin pivoteo)",
        "-" * 90,
    ]

    for i in range(len(r_esc["r"])):
        lineas.append(
            f"  {i + 1:2d} | {r_esc['r'][i]:+30.18e} | "
            f"{r_ref_f64['r'][i]:+30.18e} | {r_sin['r'][i]:+24.14e}"
        )

    lineas.extend([
        "-" * 90,
        f"  Norma infinito ||r||_inf (escalado directo)  : {r_esc['norma_inf']:.18e}",
        f"  Norma infinito ||r||_inf (escalado refinado) : {r_ref_f64['norma_inf']:.18e}",
        "",
        "3. REFINAMIENTO ITERATIVO (reutilizando L y U en float64, r en Decimal 50 dígitos)",
        "-" * 90,
        "  a) Sobre pivoteo parcial escalado (actualizando x en float64):",
    ])
    for p in ref_esc["pasos_f64"]:
        lineas.append(
            f"     Iteración {p['iteracion']}: ||r_antes||_inf = {p['norma_r_antes']:.6e} -> "
            f"||dx||_inf = {p['norma_dx']:.6e} (rel = {p['rel_dx']:.3e}) -> "
            f"||r_despues||_inf = {p['norma_r_despues']:.6e} [{p['flops']} FLOPs]"
        )

    lineas.append("  b) Corrida auxiliar acumulando x en Decimal (mismos factores L y U en float64):")
    for it, ndx, nr in ref_esc["pasos_decimal"]:
        lineas.append(f"     Iteración {it}: ||dx||_inf = {ndx:.6e} -> ||r_despues||_inf = {nr:.6e}")

    lineas.append("  c) Sobre corrida sin pivoteo (actualizando x en float64):")
    for p in ref_sin["pasos_f64"]:
        lineas.append(
            f"     Iteración {p['iteracion']}: ||r_antes||_inf = {p['norma_r_antes']:.6e} -> "
            f"||dx||_inf = {p['norma_dx']:.6e} (rel = {p['rel_dx']:.3e}) -> "
            f"||r_despues||_inf = {p['norma_r_despues']:.6e} [{p['flops']} FLOPs]"
        )

    lineas.extend([
        "",
        "  Observación sobre ||r||_inf ≈ 2.603e-10 en float64:",
        "  En la iteración 2, ||dx||_inf = 1.712e-16 (< épsilon de máquina), lo que indica que x",
        "  ya llegó al límite de representación de float64 (53 bits). El residuo ~2.60e-10 viene",
        "  solo de multiplicar el error de redondeo binario de x (~5e-17) por las entradas de la",
        "  fila 6 (~10^6). Si se acumula x en Decimal, ||r||_inf cae a 1.37e-22 y 3.93e-35.",
        "",
        "4. CONDICIONAMIENTO Y CASI DEPENDENCIA LINEAL",
        "-" * 90,
        f"  det(A)                                : {cond['det_a']:.6e}",
        f"  κ_2(A)  (matriz original)             : {cond['cond_2']:.6e}",
        f"  κ_∞(A)  (matriz original)             : {cond['cond_inf']:.6e}",
        f"  κ_2(D⁻¹A) (escalada por filas)        : {cond['cond_2_escalada']:.6e}",
        f"  κ_∞(D⁻¹A) (escalada por filas)        : {cond['cond_inf_escalada']:.6e}",
        f"  Razón de escalas (s_max / s_min)      : {cond['razon_escalas']:.6e} "
        f"({cond['escala_min']:.3e} .. {cond['escala_max']:.3e})",
        f"  Cota clásica κ_2(A)·ε                 : {cond['cota_clasica']:.6e}",
        f"  Cota escalada κ_2(D⁻¹A)·ε             : {cond['cota_escalada']:.6e}",
        f"  Casi dependencia F7 ≈ α·F3 + β·F5     : α = {cond['coef_f3']:.6f}, β = {cond['coef_f5']:.6f}",
        f"    ||F7 - (α·F3 + β·F5)||_2            : {cond['norma_dif_f7']:.6e} (relativo: {cond['rel_dif_f7']:.6e})",
        "",
        "5. CONTROL EXTERNO CONTRA numpy.linalg.solve",
        "-" * 90,
        f"  Máx. dif. relativa |x_escalado - x_np| / |x_np| : {ctrl_np['dif_rel_directo']:.6e}",
        f"  Máx. dif. relativa |x_refinado - x_np| / |x_np| : {ctrl_np['dif_rel_refinado']:.6e}",
        "",
    ])

    ruta = Path(ruta_salida)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    ruta.write_text("\n".join(lineas), encoding="utf-8")
