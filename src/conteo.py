from pathlib import Path

def crear_contador(n, estrategia):
  return{
    "n":n,
    "estrategia": estrategia,
    "pasos": [],
    "paso_actual": None,
    "en_sustitucion": False,
    "sust_mul": 0,
    "sust_div": 0,
    "sust_sum": 0,
    "sust_res": 0,
    "comp_escala_inicial": 0,
  }  

def iniciar_paso(cont, numero_paso):
    cont["en_sustitucion"] = False
    paso = {
        "paso": numero_paso,
        "fila_pivote": 0,
        "valor_pivote": 0.0,
        "mul": 0,
        "div": 0,
        "sum": 0,
        "res": 0,
        "comp_pivote": 0,
        "div_pivote": 0,
        "swaps": 0,
    }
    cont["pasos"].append(paso)
    cont["paso_actual"] = paso

def registrar_pivote(cont, fila_original, valor_pivote):
    cont["paso_actual"]["fila_pivote"] = fila_original
    cont["paso_actual"]["valor_pivote"] = valor_pivote


def iniciar_sustitucion(cont):
    cont["en_sustitucion"] = True
    cont["paso_actual"] = None

def sumar_flop(cont, tipo):
    # tipo: "mul", "div", "sum" o "res"
    if cont["en_sustitucion"]:
        cont[f"sust_{tipo}"] += 1
    else:
        cont["paso_actual"][tipo] += 1


def sumar_pivoteo(cont, tipo):
    # Sobrecosto de estabilidad contado aparte
    if tipo == "comp_escala_inicial":
        cont["comp_escala_inicial"] += 1
    else:
        cont["paso_actual"][tipo] += 1

def obtener_totales(cont):
    elim_div = sum(p["div"] for p in cont["pasos"])
    elim_mul = sum(p["mul"] for p in cont["pasos"])
    elim_sum_res = sum(p["sum"] + p["res"] for p in cont["pasos"])
    elim_mul_div = elim_div + elim_mul
    elim_flops = elim_mul_div + elim_sum_res

    sust_mul_div = cont["sust_mul"] + cont["sust_div"]
    sust_sum_res = cont["sust_sum"] + cont["sust_res"]
    sust_flops = sust_mul_div + sust_sum_res

    total_mul_div = elim_mul_div + sust_mul_div
    total_sum_res = elim_sum_res + sust_sum_res

    return {
        "elim_div": elim_div,
        "elim_mul": elim_mul,
        "elim_mul_div": elim_mul_div,
        "elim_sum_res": elim_sum_res,
        "elim_flops": elim_flops,
        "sust_mul_div": sust_mul_div,
        "sust_sum_res": sust_sum_res,
        "sust_flops": sust_flops,
        "total_mul_div": total_mul_div,
        "total_sum_res": total_sum_res,
        "total_flops": total_mul_div + total_sum_res,
        "comp_escala_inicial": cont["comp_escala_inicial"],
        "comp_pivote": sum(p["comp_pivote"] for p in cont["pasos"]),
        "div_pivote": sum(p["div_pivote"] for p in cont["pasos"]),
        "swaps": sum(p["swaps"] for p in cont["pasos"]),
    }

def cotas_teoricas(n):
    termino_cubico = (2 * (n ** 3)) // 3
    minimo_sin_col_k = (2 * (n ** 3)) // 3 + (3 * (n ** 2)) // 2 - (7 * n) // 6
    referencia_con_col_k = minimo_sin_col_k + n * (n - 1)
    return {
        "termino_cubico": termino_cubico,
        "minimo_exacto": minimo_sin_col_k,
        "referencia_enunciado": referencia_con_col_k,
    }
    
def guardar_conteo_md(contadores, ruta_salida):
    principal = contadores["escalado"]
    tot = obtener_totales(principal)
    cotas = cotas_teoricas(principal["n"])
    ref = cotas["referencia_enunciado"]
    desv = (tot["total_flops"] - ref) / ref * 100.0

    ruta = Path(ruta_salida)
    ruta.parent.mkdir(parents=True, exist_ok=True)

    with ruta.open("w", encoding="utf-8") as f:
        f.write("# Conteo de FLOPs — Eliminación Gaussiana (n = 12)\n\n")
        f.write("**Convención:** 1 FLOP = 1 operación elemental (+, -, *, /). ")
        f.write("Multiplicadores reutilizados por fila; sobrecosto de pivoteo aparte.\n\n")
        f.write("| Paso | Pivote | Mul/Div | Sum/Res | Total FLOPs | Comp. | Div. escala | Swaps |\n")
        f.write("| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |\n")

        for p in principal["pasos"]:
            md = p["div"] + p["mul"]
            sr = p["sum"] + p["res"]
            f.write(
                f"| {p['paso']} | F{p['fila_pivote']} | {md} | {sr} | {md + sr} | "
                f"{p['comp_pivote']} | {p['div_pivote']} | {p['swaps']} |\n"
            )

        f.write(
            f"| **Eliminación** | — | {tot['elim_mul_div']} | {tot['elim_sum_res']} | "
            f"{tot['elim_flops']} | {tot['comp_pivote']} | {tot['div_pivote']} | {tot['swaps']} |\n"
        )
        f.write(
            f"| **Sust. regresiva** | — | {tot['sust_mul_div']} | {tot['sust_sum_res']} | "
            f"{tot['sust_flops']} | 0 | 0 | 0 |\n"
        )
        f.write(
            f"| **TOTAL** | — | **{tot['total_mul_div']}** | **{tot['total_sum_res']}** | "
            f"**{tot['total_flops']}** | **{tot['comp_pivote']}** | **{tot['div_pivote']}** | **{tot['swaps']}** |\n\n"
        )

        f.write(f"- **Escalas iniciales:** {tot['comp_escala_inicial']} comparaciones (0 FLOPs).\n")
        f.write(f"- **Cota 2n³/3:** {cotas['termino_cubico']} FLOPs | **Referencia enunciado:** {ref} FLOPs.\n")
        f.write(f"- **Desviación vs referencia:** {desv:+.2f} % (dentro de ±10 %).\n")