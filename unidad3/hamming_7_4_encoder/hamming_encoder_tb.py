import cocotb
from cocotb.triggers import Timer

from galois import GFPoly, GaloisField

F = GaloisField(3, 0b011)          # GF(2^3), P(x) = x^3 + x + 1
G = GFPoly(F, [1, 0, 1, 1])        # g(x) = x^3 + x + 1
A = F.alpha


# el bit j del entero pasa a ser el coeficiente de x^j
def poly(valor, ancho):
    return GFPoly(F, [(valor >> i) & 1 for i in reversed(range(ancho))])

# operacion inversa: el coeficiente de x^j vuelve a ser el bit j
def bits(p, ancho):
    return sum(p.coeff(i) << i for i in range(ancho))

# encoder de referencia: la paridad es el resto de dividir por g(x)
def codificar(u):
    return u | (bits((poly(u, 4) * GFPoly(F, [1, 1])) % G, 3) << 4)

# chequeo de paridad: toda palabra del codigo se anula en alpha
def sindrome(v):
    return poly(v, 7).eval(A)

# separa los 3 bits de paridad de los 4 de informacion
def fmt(v):
    return f"{v >> 4:03b}_{v & 0xF:04b}"


# el encoder es combinacional: aplicamos y esperamos que propague
async def aplicar(dut, u):
    dut.i_info.value = u
    await Timer(2, "ns")
    return int(dut.o_message.value)


@cocotb.test()
async def test_hamming_encoder(dut):
    await Timer(1, "ns")

    print(f"\n{F}")
    print(f"  {'i':>2} | {'a^i':>3} | {'bin':>3} | polinomio")
    print(f"  {'-'*2}-+-{'-'*3}-+-{'-'*3}-+-{'-'*15}")
    for i in range(F.n):
        a = F.exp(i)
        print(f"  {i:>2} | {a:>3} | {a:03b} | {F.poly_str(a)}")

    # alpha es raiz de g y genera el grupo multiplicativo
    assert G.eval(A) == 0 and F.power(A, F.n) == 1
    print(f"\ng(x) = {G}   g(a) = {G.eval(A)}   a^7 = {F.power(A, F.n)}")

    print("\nBarrido de las 16 palabras de informacion")
    print(f"  {'u':>4} | {'v = o_message':>13} | {'v(x)':>35} | {'v(a)':>4} | {'v mod g':>7}")
    print(f"  {'-'*4}-+-{'-'*13}-+-{'-'*35}-+-{'-'*4}-+-{'-'*7}")

    codigo = {}

    for u in range(16):
        v = await aplicar(dut, u)
        vp = poly(v, 7)
        assert v == codificar(u), f"u={u:04b}: DUT={v:07b} modelo={codificar(u):07b}"
        assert v & 0xF == u, f"u={u:04b}: el codigo no es sistematico"
        assert sindrome(v) == 0, f"u={u:04b}: v(a) = {sindrome(v)}"
        assert (vp % G).is_zero, f"u={u:04b}: g(x) no divide a v(x)"
        codigo[u] = v
        print(f"  {u:04b} | {fmt(v):>13} | {str(vp):>35} | {sindrome(v):>4} | {str(vp % G):>7}")

    # el codigo es un subespacio lineal de dimension 4
    assert len(set(codigo.values())) == 16
    for a in range(16):
        for b in range(16):
            assert codigo[a ^ b] == codigo[a] ^ codigo[b], f"linealidad rota en {a},{b}"
    print("\n16 palabras distintas y enc(a^b) = enc(a)^enc(b)")

    # la distancia minima acota lo que despues puede corregir el decoder
    dmin = min(bin(codigo[a] ^ codigo[b]).count("1")
               for a in range(16) for b in range(a + 1, 16))
    assert dmin == 3
    print(f"d_min = {dmin}  =>  detecta 2 errores, corrige 1")

    print("\nBarrido sobre los elementos no nulos del campo")
    print(f"  {'j':>2} | {'a^j':>3} | {'bin':>3} | {'v con error en j':>16} | {'v(a)':>4}")
    print(f"  {'-'*2}-+-{'-'*3}-+-{'-'*3}-+-{'-'*16}-+-{'-'*4}")

    v = codigo[0b1011]
    for j in range(7):
        r = v ^ (1 << j)
        # dar vuelta el bit j deja un sindrome igual a alpha^j
        assert sindrome(r) == F.power(A, j) and F.log(sindrome(r)) == j
        print(f"  {j:>2} | {F.power(A, j):>3} | {F.power(A, j):03b} | {fmt(r):>16} | "
              f"{sindrome(r):>4}")

    # cada posicion de error da un sindrome distinto: por eso se puede corregir
    assert len({F.power(A, j) for j in range(7)}) == 7
    for u in range(16):
        for j in range(7):
            assert F.log(sindrome(codigo[u] ^ (1 << j))) == j
    print("\n112 casos de error simple: el logaritmo del sindrome es la posicion")