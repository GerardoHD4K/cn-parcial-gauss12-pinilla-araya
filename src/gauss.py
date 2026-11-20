from pathlib import Path
from conteo import (
    crear_contador,
    iniciar_paso,
    iniciar_sustitucion,
    obtener_totales,
    registrar_pivote,
    sumar_flop,
    sumar_pivoteo,
)

def copiar_matriz(M):
    return [fila[:] for fila in M]


def verificar_subdiagonal(U):
    n = len(U)
    ceros = 0
    total = 0
    max_sub = 0.0
    for i in range(1, n):
        for j in range(i):
            total += 1
            val = abs(U[i][j])
            if val == 0.0:
                ceros += 1
            if val > max_sub:
                max_sub = val
    return ceros, total, max_sub

def eliminacion_gaussiana(A,b, pivoteo="escalado"):
    # pivoteo puede ser: "ninguno", "parcial" o "escalado"
    if pivoteo not in ("ninguno", "parcial", "escalado"):
        raise ValueError(f"Modo de pivoteo no válido: {pivoteo}")

    n = len(A)
    U = copiar_matriz(A)
    b_mod = b[:]
    # Guardamos los multiplicadores en L (con 1.0 en la diagonal)
    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    perm = list(range(n))

    cont = crear_contador(n, pivoteo)

    # Escala inicial de cada fila: s_i = max_j |a(i,j)|
    escalas = []
    for i in range(n):
        mayor = abs(U[i][0])
        for j in range(1, n):
            if pivoteo == "escalado":
                sumar_pivoteo(cont, "comp_escala_inicial")
            if abs(U[i][j]) > mayor:
                mayor = abs(U[i][j])
        if mayor == 0.0:
            raise ZeroDivisionError(f"Fila {i + 1} nula: sistema singular.")
        escalas.append(mayor)

    max_orig = max(escalas)
    max_durante = max_orig
    max_mult = 0.0
    max_mult_escalado = 0.0
    traza = []

    # Fase 1: Eliminación hacia adelante (11 pasos)
    for k in range(n - 1):
        iniciar_paso(cont, k + 1)
        fila_antes = perm[k] + 1
        piv = k

        if pivoteo == "parcial":
            mejor = abs(U[k][k])
            for i in range(k + 1, n):
                sumar_pivoteo(cont, "comp_pivote")
                if abs(U[i][k]) > mejor:
                    mejor = abs(U[i][k])
                    piv = i

        elif pivoteo == "escalado":
            sumar_pivoteo(cont, "div_pivote")
            mejor = abs(U[k][k]) / escalas[k]
            for i in range(k + 1, n):
                sumar_pivoteo(cont, "div_pivote")
                sumar_pivoteo(cont, "comp_pivote")
                razon = abs(U[i][k]) / escalas[i]
                if razon > mejor:
                    mejor = razon
                    piv = i

        hubo_swap = piv != k
        if hubo_swap:
            U[k], U[piv] = U[piv], U[k]
            b_mod[k], b_mod[piv] = b_mod[piv], b_mod[k]
            escalas[k], escalas[piv] = escalas[piv], escalas[k]
            perm[k], perm[piv] = perm[piv], perm[k]
            for col in range(k):
                L[k][col], L[piv][col] = L[piv][col], L[k][col]
            sumar_pivoteo(cont, "swaps")

        pivote = U[k][k]
        if pivote == 0.0:
            raise ZeroDivisionError(f"Pivote cero en el paso {k + 1}.")

        fila_piv_orig = perm[k] + 1
        registrar_pivote(cont, fila_piv_orig, pivote)
        max_m_paso = 0.0

        # Eliminar los elementos bajo el pivote en la columna k
        for i in range(k + 1, n):
            # Calculamos el multiplicador una sola vez por fila
            m = U[i][k] / pivote
            sumar_flop(cont, "div")
            L[i][k] = m

            abs_m = abs(m)
            m_esc = abs_m * (escalas[k] / escalas[i])
            if abs_m > max_m_paso:
                max_m_paso = abs_m
            if abs_m > max_mult:
                max_mult = abs_m
            if m_esc > max_mult_escalado:
                max_mult_escalado = m_esc

            # Asignar 0.0 directo ahorra 2 FLOPs por fila y evita residuos
            U[i][k] = 0.0

            # Recorremos solo de la columna k + 1 en adelante
            for j in range(k + 1, n):
                prod = m * U[k][j]
                sumar_flop(cont, "mul")
                U[i][j] = U[i][j] - prod
                sumar_flop(cont, "res")
                if abs(U[i][j]) > max_durante:
                    max_durante = abs(U[i][j])

            # Actualizar lado derecho b
            prod_b = m * b_mod[k]
            sumar_flop(cont, "mul")
            b_mod[i] = b_mod[i] - prod_b
            sumar_flop(cont, "res")

        ceros_act, _, _ = verificar_subdiagonal(U)
        traza.append({
            "paso": k + 1,
            "fila_pivote": fila_piv_orig,
            "hubo_swap": hubo_swap,
            "fila_antes": fila_antes,
            "valor_pivote": pivote,
            "max_mult_paso": max_m_paso,
            "ceros_acumulados": ceros_act,
            "permutacion": perm[:],
            "U": copiar_matriz(U),
            "b": b_mod[:],
        })

    # Fase 2: Sustitución regresiva (Ux = b_mod)
    iniciar_sustitucion(cont)
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        suma = b_mod[i]
        for j in range(i + 1, n):
            prod = U[i][j] * x[j]
            sumar_flop(cont, "mul")
            suma = suma - prod
            sumar_flop(cont, "res")
        x[i] = suma / U[i][i]
        sumar_flop(cont, "div")

    ceros_sub, total_sub, max_sub = verificar_subdiagonal(U)

    return {
        "pivoteo": pivoteo,
        "x": x,
        "L": L,
        "U": U,
        "permutacion": perm,
        "traza": traza,
        "contador": cont,
        "max_mult": max_mult,
        "max_mult_escalado": max_mult_escalado,
        "factor_crecimiento": max_durante / max_orig,
        "ceros_bajo_diagonal": ceros_sub,
        "total_bajo_diagonal": total_sub,
        "max_bajo_diagonal": max_sub,
    }

def resolver_lu(L, U, perm, rhs):
    # Reutiliza L y U ya factorizados para resolver A*d = rhs en O(n^2)
    n = len(rhs)
    pb = [rhs[idx] for idx in perm]
    flops = 0

    y = [0.0] * n
    for i in range(n):
        acum = pb[i]
        for j in range(i):
            acum -= L[i][j] * y[j]
            flops += 2
        y[i] = acum

    d = [0.0] * n
    for i in range(n - 1, -1, -1):
        acum = y[i]
        for j in range(i + 1, n):
            acum -= U[i][j] * d[j]
            flops += 2
        d[i] = acum / U[i][i]
        flops += 1

    return d, flops


def escribir_matriz(f, U, b_vec, perm):
    n = len(U)
    for i in range(n):
        vals = " ".join(f"{U[i][j]:+11.3e}" for j in range(n))
        f.write(f"  F{perm[i] + 1:<2d} | {vals} | {b_vec[i]:+12.4e}\n")


def escribir_bloque_traza(f, titulo, A_orig, b_orig, res):
    tot = obtener_totales(res["contador"])
    f.write(f"\n{'=' * 95}\n{titulo}\n{'=' * 95}\n")
    f.write(
        f"FLOPs: {tot['total_flops']} | max|m|: {res['max_mult']:.3e} | "
        f"rho: {res['factor_crecimiento']:.3e} | "
        f"ceros subdiag: {res['ceros_bajo_diagonal']}/{res['total_bajo_diagonal']}\n\n"
    )
    f.write("Paso 0 (matriz inicial [A | b]):\n")
    escribir_matriz(f, A_orig, b_orig, list(range(len(A_orig))))

    for p in res["traza"]:
        swap = f"swap F{p['fila_antes']}<->F{p['fila_pivote']}" if p["hubo_swap"] else "sin swap"
        f.write(
            f"\nPASO {p['paso']:2d} / 11 | pivote F{p['fila_pivote']} ({p['valor_pivote']:+.4e}) | "
            f"{swap} | max|m|={p['max_mult_paso']:.3e} | ceros={p['ceros_acumulados']}/66\n"
        )
        escribir_matriz(f, p["U"], p["b"], p["permutacion"])

    f.write(
        f"\nVerificación final: {res['ceros_bajo_diagonal']}/{res['total_bajo_diagonal']} "
        f"elementos nulos bajo la diagonal (máx = {res['max_bajo_diagonal']:.1e}).\n"
    )


def guardar_traza_txt(A_orig, b_orig, res_escalado, res_sin_piv, ruta_salida):
    ruta = Path(ruta_salida)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    with ruta.open("w", encoding="utf-8") as f:
        f.write("TRAZA DE ELIMINACIÓN GAUSSIANA 12x12 — Gerardo Pinilla, Matías Araya\n")
        escribir_bloque_traza(f, "PIVOTEO PARCIAL ESCALADO (11 PASOS)", A_orig, b_orig, res_escalado)
        escribir_bloque_traza(f, "SIN PIVOTEO (11 PASOS)", A_orig, b_orig, res_sin_piv)
