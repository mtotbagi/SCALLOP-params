from sage.rings.finite_rings.integer_mod import square_root_mod_prime #<- use this, otherwise takes forever

def legendre(a, p):
    return pow(a, (p - 1) // 2, p)

def tonelli(n, p):
    assert legendre(n, p) == 1, "not a square (mod p)"
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    if s == 1:
        return pow(n, (p + 1) // 4, p)
    for z in range(2, p):
        if p - 1 == legendre(z, p):
            break
    c = pow(z, q, p)
    r = pow(n, (q + 1) // 2, p)
    t = pow(n, q, p)
    m = s
    t2 = 0
    while (t - 1) % p != 0:
        t2 = (t * t) % p
        for i in range(1, m):
            if (t2 - 1) % p == 0:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        r = (r * b) % p
        c = (b * b) % p
        t = (t * c) % p
        m = i
    return r

def custom_cornacchia(d, m): #<- sage's built in gets stuck on factoring something
    # Find all solutions to x^2 + dy^2 = m
    assert is_pseudoprime(m)
    
    Z_M = Integers(m)
    
    if pow(int(-d), int((m-1)/2), m) != 1:
        return []
    #rs = [ZZ(r) for r in (X^2 + d).roots(multiplicities=False)] #This sometimes get stuck...
    #print("Sage built in:")
    #print(rs)
    proof.all(False)
    r = tonelli(-int(d), int(m))
    rs = [r, m-r]
    #print("This shit")
    #print(rs)
    bound = 2**(RR(log(m, 2))/2)
    for r in rs:
        n = m
        while r > bound:
            n, r = r, n%r
        s = sqrt((m - r^2)/d)
        if s in ZZ:
            return [r, s]
    return []

def Cornacchia(m):
    #m_prime = prod([l**e for l, e in factor(m, limit=1000) if l < 100])
    if not ((m % 4 == 1) and is_pseudoprime(m)):
        return None, None, False
    print("Potential solution??")
    sol = custom_cornacchia(1, m)
    if len(sol) == 0:
        return None, None, False
    return sol[0], sol[1], True
    

if __name__=="__main__":
    proof.all(False)

    ### 512 ####
    f = Integer(24)
    e = 256
    N = 60196370408045148140207657188808862146370697220999958612942173757586470395546138019678238953106463774301234671040320507948845836764922766902389652812026810195778680612171536245745535
    d = 113481268849242306765395827905964600750727332653427994103871703137003013203711
    conductor = 343661038852806362288976515503125130369
    fullG = 3 * 5^2 * 7^2 * 11^2 * 13 * 17^2 * 19 * 23 * 29 * 31 * 37 * 41 * 43 * 47 * 53 * 61 * 67 * 71 * 79 * 101 * 109 * 113 * 131 * 137 * 149 * 157 * 163 * 181 * 193 * 199 * 223 * 227 * 229 * 239 * 251 * 257 * 269 * 271 * 283 * 313 * 331 * 337 * 349 * 373 * 419 * 431 * 443 * 449 * 457 * 463 * 479 * 491 * 503 * 523 * 557 * 563 * 571 * 587 * 593 * 601 * 613 * 643 * 647 * 653 * 683 * 739 * 743 * 751 * 769 * 773 * 797 * 811 * 839 * 853 * 857 * 863 * 877 * 881 * 883 * 911 * 937 * 941

    ### 1024 ###
    #f = Integer(817)
    #f = Integer(83)
    #e = 518
    #N = 22711775683133590811692237376342429238176828137256920271016701278225953126746258369874072200439826127375102732028076826684141799893454573890942633774296721382453081077340536135990593000334249
    #d = 810544624661213367964996895060809843164181099028487384400155560066253019789927
    #conductor = 340282366920938463463374607431776310331*340282366920938463463374607431760112581*340282366920938463463374607431770911081
    #fullG = 3 * 5^3 * 7^2 * 11^3 * 13^2 * 17 * 19 * 23 * 29 * 31 * 37 * 41 * 43 * 47 * 53 * 61 * 71 * 73 * 79 * 89 * 97 * 137 * 139 * 151 * 163 * 167 * 181 * 193 * 199 * 223 * 239 * 241 * 257 * 281 * 311 * 317 * 331 * 349 * 353 * 367 * 373 * 397 * 401 * 409 * 419 * 421 * 433 * 457 * 461 * 463 * 487 * 499 * 509 * 541 * 547 * 569 * 571 * 577 * 587 * 593 * 617 * 619 * 631 * 641 * 659 * 691 * 719 * 727 * 739 * 743 * 751 * 757 * 761 * 773 * 787 * 797 * 827 * 829 * 853 * 857 * 863 * 881 * 941 * 953 * 967 * 971

    ### 2048 ###
    #f = Integer(371)
    #e = 1554
    #d = 810544624661213367964996895054442262638792765109240414998632320479971130861487
    #conductor = 340282366920938463463374607422403793911*340282366920938463463374607434889683971*115792089237316195423570985002314815519172857115397030761886816432492353021871*115792089237316195423570985007413289386450559155591857383943430492828974316323*340282366920938463463374607441132629001
    #fullG = 5^4 * 7^4 * 11^2 * 13^2 * 17^2 * 19^2 * 23 * 29^2 * 31^2 * 37^2 * 41^2 * 43^2 * 47^2 * 53^3 * 59^2 * 61^2 * 67^2 * 71^2 * 73^2 * 79^2 * 83^2 * 89 * 97 * 101 * 103 * 107 * 109 * 113 * 127 * 131 * 137 * 139 * 149 * 151 * 157 * 163 * 167 * 173 * 179 * 181 * 191 * 193 * 197 * 199 * 211 * 223 * 227 * 229 * 233 * 239 * 241 * 251 * 257 * 263 * 269 * 271 * 277 * 281 * 283 * 293 * 307 * 311 * 313 * 317 * 331 * 337 * 347 * 349 * 353 * 359 * 367 * 373 * 379 * 383 * 389 * 397 * 401 * 409 * 419 * 421 * 431 * 433 * 439 * 443 * 449 * 457 * 461 * 463 * 467 * 479 * 487 * 491 * 499 * 503 * 509 * 521 * 523 * 541 * 547 * 557 * 563 * 569 * 571 * 601 * 613 * 617 * 743 * 853 * 937 * 967 * 1031 * 1171 * 1181 * 1201 * 1321 * 1741 * 1801 * 1873 * 1901 * 2129 * 2297 * 2311 * 2833 * 3121 * 3169 * 3673 * 4019 * 4561 * 4591 * 4831 * 4889 * 5153 * 5591 * 5623 * 5669 * 9649 * 10651 * 10711 * 15607 * 16561 * 21313

    p = Integer(2)**e*f*N - 1

    print(p)
    assert is_pseudoprime(p)
   
    B = QuaternionAlgebra(-1, -p)
    i, j, k = B.gens()
    O0 = B.quaternion_order([1, i, (i+j)/2, (1+k)/2])
    
    Z_p = Integers(p)
    Bound = 1000

    omega = None

    print(RR(log(d*(conductor)**2, p)))

    
    found = False
    for div in range(1, 1000000):
        if (fullG) % div != 0:
            continue
        path = (fullG)/div
        n = d*(conductor*path)**2
        print(f"Should be more than 2: {RR(log(n, p))}")
        #print(RR(log(n, p)))
        #rts = [ZZ(r) for r in Z_p(n).sqrt(all=True)]
        if pow(int(n), int((p-1)/2), p) != 1:
            print("Dont think this should happen")
            continue
        r = square_root_mod_prime(Integers(p)(n))
        if not r:
            print("Dont think this should happen")
            continue
        rts = [ZZ(r), ZZ(p-r)]
        for x in rts:
            M = ZZ((n - x**2)/p)
            if M > 0:
                print(RR(log(abs(M), 2)))
            y, z, sol = Cornacchia(M)
            if sol:
                print("Cornacchia gave a sol")
                omega = x*i + y*j + z*k
                if omega/2 not in O0 and omega/3 not in O0:
                    print("We got it!!!")
                    found = True
                    break
        if found:
            break

    assert found
    
    I = O0*path + O0*omega


    print(I)
    print(RR(log(I.norm(), p)))
    assert omega/path in I.right_order()
    print(omega/path)
    print(omega)
    print(omega/2 in I.left_order())
    print(omega/2 in O0)
    print(omega/(path*2) in I.right_order())
    