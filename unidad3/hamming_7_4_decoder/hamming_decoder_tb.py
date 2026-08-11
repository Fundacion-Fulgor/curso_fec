import cocotb
from cocotb.clock import Clock
from cocotb.triggers import RisingEdge, Timer

from galois import GFPoly, GaloisField

F = GaloisField(3, 0b011)          # GF(2^3), P(x) = x^3 + x + 1
G = GFPoly(F, [1, 0, 1, 1])        # g(x) = x^3 + x + 1
A = F.alpha
K = F.power(A, 3)                  # la H del RTL da s = a^3 * r(a)


# el bit j del entero pasa a ser el coeficiente de x^j
def poly(valor, ancho):
    return GFPoly(F, [(valor >> i) & 1 for i in reversed(range(ancho))])

# chequeo de paridad: evaluamos el mensaje recibido en alpha
def sindrome(r):
    return F.mul(K, poly(r, 7).eval(A))

# decoder de referencia: sindrome nulo => la palabra es valida
def modelo(r):
    s = sindrome(r)
    if s == 0:
        return s, r, r & 0xF, 0
    # si el sindrome no es nulo, el logaritmo del sindrome es la posicion del bit erroneo
    j = F.log(F.div(s, K))
    c = r ^ (1 << j)
    return s, c, c & 0xF, 1

# separa los 3 bits de paridad de los 4 de informacion
def fmt(v):
    return f"{v >> 4:03b}_{v & 0xF:04b}"


@cocotb.test()
async def test_hamming_decoder(dut):
    # creamos un clock de periodo 10ns
    cocotb.start_soon(Clock(dut.i_clk, 10, "ns").start())
    dut.i_rst.value = 1
    dut.i_message.value = 0
    await RisingEdge(dut.i_clk)
    await RisingEdge(dut.i_clk)
    dut.i_rst.value = 0

    print(f"\n{F}")
    print(f"g(x) = {G}   a^3 = {K}   ->  s = a^3 * r(a)")

    print("\nBarrido de las 128 palabras de 7 bits")
    print(f"  {'r recibido':>10} | {'s':>3} | {'zero':>4} | {'corregido':>10} | "
          f"{'o_info':>6} | {'is_corr':>7}")
    print(f"  {'-'*10}-+-{'-'*3}-+-{'-'*4}-+-{'-'*10}-+-{'-'*6}-+-{'-'*7}")

    for k in range(129):
        if k < 128:
            dut.i_message.value = k
        # esperamos un ciclo de clock
        await RisingEdge(dut.i_clk)
        await Timer(1, "ns")

        if k < 128:
            s_esp, _, _, _ = modelo(k)
            s = int(dut.u_syndrome_calculator.o_syndrome.value)
            z = int(dut.u_syndrome_calculator.o_zero_syndrome.value)
            assert s == s_esp, f"r={k:07b}: sindrome {s}, modelo {s_esp}"
            assert z == (s_esp == 0), f"r={k:07b}: o_zero_syndrome {z}"

        if k >= 1:
            r = k - 1 # nos desplazamos un ciclo por la latencia del modulo
            s_esp, corr_esp, info_esp, flag_esp = modelo(r)
            info = int(dut.o_info.value)
            flag = int(dut.o_is_corrected.value)
            assert info == info_esp, f"r={r:07b}: o_info {info:04b}, modelo {info_esp:04b}"
            assert flag == flag_esp, f"r={r:07b}: o_is_corrected {flag}"

            print(f"  {fmt(r):>10} | {s_esp:>3} | {int(s_esp == 0):>4} | {fmt(corr_esp):>10} | "
                  f"{info:04b}   | {flag:>7}")