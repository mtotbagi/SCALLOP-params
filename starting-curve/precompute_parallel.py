from sage.all import *
from utilities import BiDLP
from ec import *
from xonly import xPoint
from subprocess import run
from ast import literal_eval
import time
from multiprocessing import Pool
import os
import sys
from random import shuffle
####################################
#                                  #
#     cost.py from DeuringFTP      #
#                                  #
####################################



def smallest_div(q, l):
    mi = -1
    for k in l.keys():
        if(q.divides(l[k])):
            if(mi == -1): mi = k
            elif(cost(k) < cost(mi)): mi = k
    return mi

def smooth_part(n, B):
    r"""
    Compute the maximal odd divisor of n which is B-smooth.
    """
    n = ZZ(n)
    n = n.prime_to_m_part(2)
    return n.factor(limit=B)

def which_pol_from_extdeg(extdeg):
    if(extdeg % 4 == 0): return extdeg
    elif (extdeg % 2 == 1): return 2*extdeg
    else: return extdeg//2

R = ZZ['x']
def cyclpol_from_extdeg(extdeg):
    return R.cyclotomic_polynomial(which_pol_from_extdeg(extdeg))

def cost(k):
    if k == 1: return 0
    k = Integer(k)
    poldeg = euler_phi(which_pol_from_extdeg(k))
    mul = poldeg * k * float(sqrt(k)) / 2.58
    return mul

def choose_torsion(p, lowbound, max_ext = 100, max_tors = 50000):
    facToExt = dict()
    ords = {i: cyclpol_from_extdeg(i)(p) for i in range(1, max_ext)}
    le = ZZ.one()
    while lcm(le for le in facToExt.keys()) < lowbound*2**20 or le < max_tors:
        le += 1
        if p.divides(le):
            continue
        if le % 2 == 0:
            continue
        if not is_pseudoprime(radical(le)):
            continue
        k = smallest_div(le, ords)
        if k == -1: continue
        facToExt[le] = k

    
    #print("Candidates for torsions (+extdegs), sorted by cost")
    #print(tups)
    #print()

    extToFac = dict()
    for k in range(1, 100):
        extToFac[k] = []
        for le in facToExt.keys():
            if facToExt[le] == k: extToFac[k].append(le)
        if extToFac[k] == []: del extToFac[k]

    tups = sorted(extToFac.items(), key=lambda item: cost(item[0])/lcm(item[1]))
    T = ZZ.one()
    it = 0
    for k, torsions in tups:
        T = lcm(T, lcm(torsions))
        it += 1
        if(T >= lowbound):
            break
    
    assert T >= lowbound
    print("Extdegs with their torsions")
    extToFac = dict(tups[:it])
    for extdeg in extToFac.keys():
        print(f"{extdeg}")
        print(f"{extToFac[extdeg]}")

    def shave_excess(T, extToFac, B):
        tups = sorted(extToFac.items(), key=lambda item: -cost(item[0])/lcm(item[1]))
        for k, l in tups:
            if T//lcm(l) > B:
                print(f"Removing torsions {l}, from extdeg {k}")
                del extToFac[k]
                return shave_excess(T//lcm(l), extToFac, B)
        return T, extToFac
    
    #We picked each extension by cost, and it is possible that with the last one we have extra torsion
    #which we can shave off
    #It is possible that this will not be "optimal", but should be close enough
    print("Shaving off excess torsions...")
    T, extToFac = shave_excess(T, extToFac, lowbound)
    
    return T, extToFac

####################################
#                                  #
#    Action matrices on basis      #
#                                  #
####################################

def endo_i(P):
    E0 = P.curve()
    F = E0.base_field()
    x, y = P.xy()
    return E0(-x, F(sqrtm1)*y)

def endo_j(P):
    E0 = P.curve()
    pi = E0.base_field().frobenius_endomorphism()
    x,y = P.xy()
    return E0(pi(x), pi(y))

def EvalEndomorphism(alpha, P, ord):
    d = lcm(c.denominator() for c in alpha)
    E = P.curve()
    assert P*ord == 0
    assert d == 1 or d == 2
    if gcd(d, ord) == 2:
        alpha = d*alpha
        Fbig, _ = E.base_field().extension(4,'A').objgen()
        Ebig = E.base_extend(Fbig)
        P = Ebig(P).division_points(2)[0]
    iP = endo_i(P)
    jP = endo_j(P)
    kP = endo_i(jP)
    coeffs = [coeff % ord for coeff in alpha]
    return E(sum(c*Q for c, Q in zip(coeffs, [P, iP, jP, kP])))

def ActionMatrix(alpha, basis, ord):
    P, Q = basis
    alphaP = EvalEndomorphism(alpha, P, ord)
    alphaQ = EvalEndomorphism(alpha, Q, ord)
    #if(alpha == B(1)):
    #    print(P, alphaP)
    a, c = BiDLP(alphaP, P, Q, ord)
    b, d = BiDLP(alphaQ, P, Q, ord)
    return [[a, b],[c,d]]

def BasisAct(base, ord):
    #t0 = time.time_ns()
    i,j,k = B.gens()
    M_0 = ActionMatrix(B(1), base, ord)
    M_1 = ActionMatrix(i, base, ord)
    M_2 = ActionMatrix((i+j)/2, base, ord)
    M_3 = ActionMatrix((1+k)/2, base, ord)
    M_theta = ActionMatrix(j + (1+k)/2, base, ord)
    #t1 = time.time_ns()
    #t = float((t1-t0) / (10 ** 9))
    #print(f"Basis action time for {ord}: {t}")
    return (M_0, M_1, M_2, M_3, M_theta)

start = time.time()
proof.all(False)

try: bitsize = int(sys.argv[1])
except: bitsize = 512
print(bitsize)

match bitsize:
    case 512:
        p = Integer(167286323857221689112346016933207258999493176647479781908348180838625562682489086996433613891517156716513168813389511283523303877305870043625319556074092350622078682027926127121250390130545025689333161800187261717666012808114506310888652381411304126074312458239)
        lowbound = 2**(0.75*RR(log(p,2)))
    case 1024:
        p = Integer(15922486913903522274738876334491173327471544233065982291422062886364771917539984633619947369979657884209891416505022976305949903782081999018601968658399212422459232014682162496366184150793214330561748936802196212819401893052051149876176865225396801900263676014808365383967944863381855436796098378947797615849584416454584458385591707302057385868132351)
        lowbound = 2**(0.6*RR(log(p,2)))
    case 2048:
        p = Integer(596189944777372638172818300422290385024394468583357632698904620718401201446423451780730293885911387310424468466907078162034493044851014467790119468529446035454449745387614406711198662080037561127345507482289795632805587649490909091103733371042789817522706777642319645749180265150818764030351009256231058220919495563789230174870643587414401676910226850987845426748837131252505799030848442553318641506513527604875273309202232111785774057997494359022678091161441638118008960288159517233147661046896608323640602694819387549541102995842204231072991679669787422227224111599861798724874762075043130922161141224785272355703617441872269175696231180105996744064900599544211922383800484570456691509756280863371782883761347598535069684820049133567)
        lowbound = 2**(0.6*RR(log(p,2)))

F = GF((p,2), name='z2', modulus=var('x')**2 + 1)
sqrtm1 = F.gens()[0]

filename = f"Precomputed/75-primes-{bitsize}.py"
try:
    open(filename, "r")
    print('precomputed file already exists!')
    quit()
except:
    print("Precompute starting:")

assert is_pseudoprime(p)


print("Choosing torsion....")

T, extToFac = choose_torsion(p, lowbound)
facToExt = dict()
for extdeg in extToFac.keys():
    for D in extToFac[extdeg]:
        facToExt[D] = extdeg
T = ZZ(T)

f = (p+1).valuation(2)

sec_param = 128

f2 = min(f, 128)
if f2 < sec_param:
    f3 = Integer(p + 1).valuation(3)
    D_chall = 2**f2*3**f3
else:
    f3 = 0
    D_chall = 2**f2
D_com = T.prime_to_m_part(3**f3)

#assert D_com > p
assert D_chall >= 2**sec_param

assert p % 4 == 3
B = QuaternionAlgebra(-1, -p)
i, j, k = B.gens()
O0 = B.quaternion_order([1, i, (i+j)/2, (1+k)/2])

E0 = EllipticCurve(F, [1,0])
E0.set_order((p+1)**2, num_checks=0)

facToBasis = {}
facToAction = {}

B_2 = TorsionBasis(E0, 2**f, xOnly = 1)
P2, Q2 = B_2
P2.set_order(2**f)
Q2.set_order(2**f)
facToAction[2**f] = BasisAct(B_2, 2**f)

B_chall = TorsionBasis(E0, D_chall, xOnly = 1)
Pc, Qc = B_chall
Pc.set_order(D_chall)
Qc.set_order(D_chall)
facToAction[D_chall] = BasisAct(B_chall, D_chall)

degToField = {1 : F}
degToModulus = {1 : str(var('x')**2 + 1)}
degToi = {1 : sqrtm1}
for extdeg in extToFac.keys():
    print(f"{extdeg}")
    print(f"{lcm(extToFac[extdeg])}")

#Testing that every torsion exists, and we find big enough torsion
TT = 1
for k in extToFac.keys():
    assert (lcm(extToFac[k])).divides(cyclpol_from_extdeg(k)(p))
    TT = lcm(T, lcm(extToFac[k]))
assert TT > lowbound

def parallel1(extdeg_D):
    extdeg, D = extdeg_D
    s = time.time()
    print(f"Finding torsion {D} basis in extension {extdeg}")
    i=1
    while True:
        cp = run(["./c_torsion/torsion_basis", f"{bitsize}", f"{i}", f"{extdeg}", f"{D}"], text=True, capture_output=True)
        if(cp.returncode == 0): 
            print(f"Torsion {D} basis in extension {extdeg} found in {time.time()-s:.1f}s")
            return cp.stdout, extdeg, D
        if("stack overflows" not in cp.stderr):
            print(f"Error in the C code, couldn't find torsion {D} in extension {extdeg}")
            print(" Error message ", cp.stderr)
            print(" Stdout: ", cp.stdout)
            return
        print(f"Memory overflow in PARI, trying again with bigger memory (extdeg={extdeg})")
        i += 1

def nonparallel(result, extdeg, D):
    print(f"Doing nonparallel stuff for extension {extdeg}")
    s = time.time()
    lines = [l for l in result.splitlines()]
    if extdeg in degToField.keys():
        Fbig = degToField[extdeg]
    else:
        Fbig = GF((p,2*extdeg), name='z' + str(2*extdeg), modulus=literal_eval(lines[0]))
        degToField[extdeg] = Fbig
        degToModulus[extdeg] = str(Fbig.modulus())
    i = Fbig(literal_eval(lines[1]))
    xP = Fbig(literal_eval(lines[2]))
    yP = Fbig(literal_eval(lines[3]))
    xQ = Fbig(literal_eval(lines[4]))
    yQ = Fbig(literal_eval(lines[5]))

    try:
        Fbig.register_coercion(F.hom([i], codomain=Fbig, check=False))
        degToi[extdeg] = i
    except: pass

    Ebig = E0.base_extend(Fbig)
    order = p**extdeg - (-1)**extdeg
    Ebig.set_order(order**2, num_checks=0)
    Pbig = Ebig(xP, yP)
    Qbig = Ebig(xQ, yQ)
    for l, ee in factor(D):
        for e in range(1, ee+1):
            if(l**e not in extToFac[extdeg]): continue
        le = l**e
        P = Pbig * (D/le)
        Q = Qbig * (D/le)
        PmQ = P-Q
        xP, xQ, xPmQ = xPoint(P.xy()[0], Ebig), xPoint(Q.xy()[0], Ebig), xPoint(PmQ.xy()[0], Ebig)
        Basis = P, Q

        facToBasis[D] = [xP, xQ, xPmQ]
        facToAction[le] = BasisAct(Basis, le)
    print(f"Finished everything in extension {extdeg} in {time.time()-s:.1f}s")
    return i
    

extdeg_D_pairs = []
for D in extToFac[1]:
    result = parallel1((1, D))
    res, extdeg, D = result
    nonparallel(res, extdeg, D)

for extdeg in extToFac.keys():
    if(extdeg == 1): continue
    extdeg_D_pairs.append((extdeg,lcm(extToFac[extdeg])))

print(extdeg_D_pairs)
#the extdeg_D_pairs get processed in chunks. Each chunk is processed by one core
#They get shuffled so each chunk takes roughly the same time
shuffle(extdeg_D_pairs)
#Pool param is the number of cores to use
with Pool(os.cpu_count()-1) as pool:
    for result in pool.imap(parallel1, extdeg_D_pairs, 6):
        if(result is None): continue
        res, extdeg, D = result
        degToi[extdeg] = nonparallel(res, extdeg, D)


printable_facToBasis = {}
for key in facToBasis.keys():
    printable_facToBasis[key] = [[str(c) for c in P.X] for P in facToBasis[key]]

with open(filename, 'w') as f:
    f.write(f'{p}\n')
    f.write(f'{(p+1).valuation(2)}\n') #f
    f.write(f'{D_com}\n')
    f.write(f'{D_chall}\n')
    f.write(f'{T}\n')
    f.write(f'{[[str(c) for c in P.xy()] for P in B_2]}\n')
    f.write(f'{[[str(c) for c in P.xy()] for P in B_chall]}\n')
    f.write(f'{facToExt}\n')
    f.write(f'{degToModulus}\n')
    f.write(f'{printable_facToBasis}\n')
    f.write(f'{facToAction}\n')

    #This is needed for fast coercion between F_p^2 and F_p^2k, the sage coercion finding takes ages
    #Use it like this:
    #i = degToi[extdeg]
    #try:
    #    Fbig.register_coercion(F.hom([i], codomain=Fbig, check=False))
    #except: pass
    f.write(f'{degToi}\n')

end = time.time()
print(end - start)