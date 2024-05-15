from sage.all import *
from xonly import isMontgomery, xPoint, sqrtDeterministic
import time
#############################
#                           #
#        Basis-stuff        #
#                           #
#############################


def is_independent(R, S, D):
    facD = factor(D)
    Drad = radical(D)
    Rsmalls = []
    Rsmall = R*(D/Drad)
    for (l, e) in facD:
        Rsmalli = Rsmall*(Drad/l)
        assert Rsmalli
        Rsmalls.append(Rsmalli)

    Ssmall = S*(D/Drad)
    for index, (l, e) in enumerate(facD):
        Ssmalli = Ssmall*(Drad/l)
        if (Rsmalls[index].weil_pairing(Ssmalli, l) == 1):
            return False
    return True

def has_order(P, D):
    Drad = radical(D)
    Dfac = factor(D)
    Psmall = P*(D/Drad)
    for (l, e) in Dfac:
        if not Psmall*(Drad/l):
            return False
    return not Psmall*Drad

def frob(P):
    p = P.x().parent().characteristic()
    xQ = P.x() ** (p ** 2)
    yQ = P.y() ** (p ** 2)
    return P.curve()(xQ, yQ)

def mult_with_pol(P, pol):
    if isinstance(pol, int) or isinstance(pol, Integer):
        return pol * P
    frobs = [P]
    for j in range(pol.degree()):
        frobs.append(-frob(frobs[-1]))
    return sum([a*Q for a, Q in zip(pol.list(), frobs)])

def mult_w_frob(P, N):
    R = ZZ['x']
    D = Integer(sqrt(P.curve().order()) // N)
    p = P.x().parent().characteristic()

    degs = []
    rem = 0

    extdeg = P.x().parent().degree() // 2
    if(extdeg % 4 == 0):
        rem = R.cyclotomic_polynomial(extdeg)(p) // D
        degs = divisors(extdeg)
        degs.pop(len(degs)-1)
    elif (extdeg % 2 == 1):
        rem = R.cyclotomic_polynomial(2*extdeg)(p) // D
        degs = [2*i for i in divisors(extdeg)]
        degs.pop(len(degs)-1)
    else:
        rem = R.cyclotomic_polynomial(extdeg//2)(p) // D
        degs = [i for i in divisors(extdeg)]
        degs.pop(len(degs)-2)

    pol = prod([R.cyclotomic_polynomial(i) for i in degs])
    res = mult_with_pol(P, pol)
    res = rem * res
    return res

    
def Torsion_Point(E, D):
    F = E.base_field()
    k = F.degree() // 2
    cof = Integer(sqrt(E.order())/D)
    i = F.gens()[0]
    x = 1 + ZZ.random_element(1000) * i
    while True:
        x += i
        try:
            P = E.lift_x(x)
        except:
            continue
        T = mult_w_frob(P, cof)
        #assert E.is_on_curve(T.x(), T.y())
        #assert has_order(T, D)
        #assert T == cof * P
        if has_order(T, D):
            T.set_order(D)
            return T

def CompleteBasis(R, D, x = 1, E0 = None):
    E = R.curve()
    F = E.base_field()
    i = F.gens()[0]

    if(E0 is not None):
        a = E0.base_field().gens()[0]
        xx, yy = R.xy()
        v = -xx
        w = (E.base_field())(a) * yy
        S = E([v, w])
        if(is_independent(R, S, D) and has_order(S, D)): 
            return S
        p = E0.base_field().characteristic()
        v =(E.base_field())(a) * yy
        S = E(xx**(p*p), yy**(p*p))
        if(is_independent(R, S, D) and has_order(S, D)): 
            return S

    for _ in range(1000):
        S = Torsion_Point(E, D)
        if is_independent(R, S, D):
            RmS = point_difference(R, S)
            if RmS != (R - S).xy()[0]:
                S = -S
                assert RmS == (R - S).xy()[0]
            S.set_order(D)
            return S
    assert False, "Something went wrong in Complete Basis..."

def TorsionBasis(E, D, xOnly = 0, E0 = None):
    start = time.time()
    P = Torsion_Point(E, D)
    Q = CompleteBasis(P, D, E0 = E0)
    end = time.time()
    print(f"Basis in torsion {D} found in: {end - start}")
    if xOnly == 1:
        return P, Q
    PmQ = P-Q
    xP, xQ, xPmQ = xPoint(P.xy()[0], E), xPoint(Q.xy()[0], E), xPoint(PmQ.xy()[0], E)
    if xOnly == 2:
        return P, Q, xP, xQ, xPmQ
    return xP, xQ, xPmQ


def point_difference(P, Q):
	#follows code from SQIsign NIST version 1.0
	#input: affine points P,Q
	#		affine curve constant ProjA = [A, [1,0]]
	#output: x-coordinate of P-Q = PmQ

	#check if all inputs are affine
    E = P.curve()
    assert isMontgomery(E)
    A = E.a2()
    Px, Qx = P.xy()[0], Q.xy()[0]
    PmQZ = Px - Qx
    t2 = Px*Qx
    t3 = t2 - 1
    t0 = PmQZ*t3
    PmQZ = PmQZ**2
    t0 = t0**2
    t1 = t2 + 1
    t3 = Px + Qx
    t1 = t1*t3
    t2 = t2*A
    t2 = 2*t2
    t1 = t1 + t2
    t2 = t1**2
    t0 = t2-t0
    assert t0.is_square()
    t0 = sqrtDeterministic(t0)
    PmQX = t0 + t1
    return PmQX/PmQZ

#############################
#                           #
#     Hashing to chal       #
#                           #
#############################

from hashlib import sha256

def hashToPoint(D, msg, E, seeds=None, small_ns=None, small_s=None):
    H = sha256()
    H.update(bytes(msg, "utf-8"))
    if small_ns and small_s:
        P, Q, seeds = TorsionBasis(E, D, small_ns=small_ns, small_s=small_s)
        s = int.from_bytes(H.digest())%D
        return P + s*Q, seeds
    else:
        if seeds:
            P, Q = TorsionBasis(E, D, seeds=seeds, small_ns=small_ns, small_s=small_s)
        else:
            P, Q = TorsionBasis(E, D, xOnly = 1)
        s = int.from_bytes(H.digest())%D
        return P + s*Q