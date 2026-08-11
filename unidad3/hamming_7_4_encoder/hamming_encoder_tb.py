import cocotb
from cocotb.triggers import Timer

from galois import GFPoly, GaloisField

F = GaloisField(3, 0b011)          # GF(2^3), P(x) = x^3 + x + 1
G = GFPoly(F, [1, 0, 1, 1])        # g(x) = x^3 + x + 1


# el bit j del entero pasa a ser el coeficiente de x^j
def poly(valor, ancho):
    return GFPoly(F, [(valor >> i) & 1 for i in reversed(range(ancho))])

# operacion inversa: el coeficiente de x^j vuelve a ser el bit j
def bits(p, ancho):
    return sum(p.coeff(i) << i for i in range(ancho))

# encoder de referencia: la paridad es el resto de dividir por g(x)
def modelo(u):
    return u | (bits((poly(u, 4) * GFPoly(F, [1, 1])) % G, 3) << 4)

# separa los 3 bits de paridad de los 4 de informacion
def fmt(v):
    return f"{v >> 4:03b}_{v & 0xF:04b}"


@cocotb.test()
async def test_hamming_encoder(dut):
    # arrancamos el tiempo antes de imprimir, asi el aviso del VCD no parte la tabla
    await Timer(1, "ns")

    print(f"\n{F}")
    print(f"  {'i':>2} | {'a^i':>3} | {'bin':>3} | polinomio")
    print(f"  {'-'*2}-+-{'-'*3}-+-{'-'*3}-+-{'-'*15}")
    for i in range(F.n):
        a = F.exp(i)
        print(f"  {i:>2} | {a:>3} | {a:03b} | {F.poly_str(a)}")
    print(f"\ng(x) = {G}")

    print("\nBarrido de los 16 mensajes")
    print(f"  {'u':>4} | {'o_message':>9} | {'modelo':>9} | {'v(x)':>35}")
    print(f"  {'-'*4}-+-{'-'*9}-+-{'-'*9}-+-{'-'*35}")

    for u in range(16):
        # el encoder es combinacional: aplicamos y esperamos que propague
        dut.i_info.value = u
        await Timer(2, "ns")
        v = int(dut.o_message.value)
        v_esp = modelo(u)

        assert v == v_esp, f"u={u:04b}: DUT={v:07b}, modelo={v_esp:07b}"
        print(f"  {u:04b} | {fmt(v):>9} | {fmt(v_esp):>9} | {str(poly(v, 7)):>35}")