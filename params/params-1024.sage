N = 2**129
L = 2*N**2

a = 2**128
a = 2**128
while a % 24 != 19:
    a += 1

params1024 = []

from tqdm import tqdm
for _ in tqdm(range(10000000)):
    a += 24
    if not is_pseudoprime(a):
        continue
    a = a
    d = L - a^2
    assert d > 0
    #f = 4*(a^2 - d)*a = 4*(2*a^2 - 2*N^2)*a = 8(a - N)(a + N)a
    f1 = int(a - N)
    f2 = int((a + N)/3)

    if all([is_pseudoprime(abs(f_fac)) for f_fac in [a, f1, f2]]):
        assert all([kronecker_symbol(-d, abs(part)) == 1 for part in [a, f1, f2]])
        
        params1024.append((d, a, f1, f2))
        K.<sq_d> = NumberField(x^2 + d)
        w = a + sq_d
        omega = w^4
        _, f = list(omega)

        print("Params:")
        print(f"d = {d}")
        print(f"d = {factor(d, limit = 10000)}")
        print(f"omega = {omega}")
        print(f"N(omega) = {factor(int(omega.norm()), limit=100)}")
        assert f == 8*3*a*f1*f2
        print(f"f = {f} = 8*3*{a}*{f1}*{f2}")
        print(f"size d: {round(log(d, 2), 5)}")
        print(f"size a: {round(log(a, 2), 5)}")
        print(f"size f1: {round(log(abs(f1), 2), 5)}")
        print(f"size f2: {round(log(abs(f2), 2), 5)}")
        print("factoriszations:")
        print(f"d = {factor(d, limit = 10000)}")
        print(f"a-1 = {factor(a-1)}")
        print(f"f1-1 = {factor(f1-1)}")
        print(f"f2-1 = {factor(f2-1)}")




