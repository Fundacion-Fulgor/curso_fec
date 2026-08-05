import cocotb
from cocotb.clock import Clock
from cocotb.triggers import RisingEdge, Timer

from galois import GFPoly, GaloisField

F = GaloisField(3, 0b011)          # GF(2^3), P(x) = x^3 + x + 1
G = GFPoly(F, [1, 0, 1, 1])        # g(x) = x^3 + x + 1
A = F.alpha
K = F.power(A, 3)                  # la H del RTL da s = a^3 * r(a)

LATENCIA = 3                       # reg de entrada + los 2 ciclos del decoder


# el bit j del entero pasa a ser el coeficiente de x^j
def poly(valor, ancho):
    return GFPoly(F, [(valor >> i) & 1 for i in reversed(range(ancho))])

# operacion inversa: el coeficiente de x^j vuelve a ser el bit j
def bits(p, ancho):
    return sum(p.coeff(i) << i for i in range(ancho))

# encoder de referencia: la paridad es el resto de dividir por g(x)
def codificar(u):
    return u | (bits((poly(u, 4) * GFPoly(F, [1, 1])) % G, 3) << 4)

# decoder de referencia: si el sindrome no es nulo, su logaritmo es la posicion
def decodificar(r):
    s = F.mul(K, poly(r, 7).eval(A))
    if s == 0:
        return r & 0xF, 0
    c = r ^ (1 << F.log(F.div(s, K)))
    return c & 0xF, 1

# separa los 3 bits de paridad de los 4 de informacion
def fmt(v):
    return f"{v >> 4:03b}_{v & 0xF:04b}"


@cocotb.test()
async def test_hamming_7_4(dut):
    # creamos un clock de periodo 10ns
    cocotb.start_soon(Clock(dut.i_clk, 10, "ns").start())
    dut.i_rst.value = 1
    dut.i_info.value = 0
    dut.i_error_vector.value = 0
    await RisingEdge(dut.i_clk)
    await RisingEdge(dut.i_clk)
    dut.i_rst.value = 0

    print(f"\n{F}")
    print(f"g(x) = {G}   a^3 = {K}   ->  s = a^3 * r(a)")

    # todos los mensajes contra todos los patrones de error posibles
    casos = [(u, e) for u in range(16) for e in range(128)]

    print(f"\nBarrido: {len(casos)} casos (16 mensajes x 128 vectores de error)")
    print(f"  {'u':>4} | {'e':>10} | {'v transmitido':>13} | {'r recibido':>10} | "
          f"{'o_info':>6} | {'is_corr':>7}")
    print(f"  {'-'*4}-+-{'-'*10}-+-{'-'*13}-+-{'-'*10}-+-{'-'*6}-+-{'-'*7}")

    ok = {}
    total = {}

    for k in range(len(casos) + LATENCIA - 1):
        if k < len(casos):
            u, e = casos[k]
            dut.i_info.value = u
            dut.i_error_vector.value = e
        # esperamos un ciclo de clock
        await RisingEdge(dut.i_clk)
        await Timer(1, "ns")

        # el encoder es combinacional: su salida ya corresponde al caso k
        if k < len(casos):
            u, e = casos[k]
            assert int(dut.message.value) == codificar(u), f"u={u:04b}: encoder"
            assert int(dut.received.value) == codificar(u) ^ e, f"u={u:04b} e={e:07b}: canal"

        # la salida del decoder llega LATENCIA-1 iteraciones despues
        if k >= LATENCIA - 1:
            u, e = casos[k - LATENCIA + 1]
            r = codificar(u) ^ e
            info_esp, flag_esp = decodificar(r)
            info = int(dut.o_info.value)
            flag = int(dut.o_is_corrected.value)
            assert info == info_esp, f"u={u:04b} e={e:07b}: o_info {info:04b} vs {info_esp:04b}"
            assert flag == flag_esp, f"u={u:04b} e={e:07b}: o_is_corrected {flag}"

            peso = bin(e).count("1")
            total[peso] = total.get(peso, 0) + 1
            ok[peso] = ok.get(peso, 0) + (info == u)
            if u == 0b1011 and peso <= 1:
                print(f"  {u:04b} | {fmt(e):>10} | {fmt(codificar(u)):>13} | {fmt(r):>10} | "
                      f"{info:04b}   | {flag:>7}")

    print("\nRecuperacion del mensaje segun el peso del vector de error")
    print(f"  {'peso':>4} | {'casos':>5} | {'recuperados':>11} | {'%':>5}")
    print(f"  {'-'*4}-+-{'-'*5}-+-{'-'*11}-+-{'-'*5}")
    for peso in sorted(total):
        print(f"  {peso:>4} | {total[peso]:>5} | {ok[peso]:>11} | "
              f"{100*ok[peso]/total[peso]:>5.1f}")

    # sin errores y con un solo error el mensaje siempre se recupera
    assert ok[0] == total[0] and ok[1] == total[1]
    # con dos errores nunca, porque el decoder corrige un tercer bit
    assert ok[2] == 0
    print(f"\n{sum(total.values())} casos verificados contra el modelo")