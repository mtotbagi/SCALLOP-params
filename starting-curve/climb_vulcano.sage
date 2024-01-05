def Cornacchia(QF, m):
    m_prime = prod([l**e for l, e in factor(m, limit=1000) if l < 1000])
    if not is_pseudoprime(m/m_prime):
        return None, None, False
    print("Potential solution??")
    sol = QF.solve_integer(m)
    if not sol:
        return None, None, False
    return sol[0], sol[1], True


if __name__=="__main__":
    f = Integer(817)
    N = Integer(22711775683133590811692237376342429238176828137256920271016701278225953126746258369874072200439826127375102732028076826684141799893454573890942633774296721382453081077340536135990593000334249)
    p = Integer(2)**518*f*N - 1

    d = 810544624661213367964996895060809843164181099028487384400155560066253019789927
    conductor = 945648148713467501094696962403454233607167517444788408618448739064461646654337307344171931120486640510821138175364584

    B = QuaternionAlgebra(-1, -p)
    i, j, k = B.gens()
    O0 = B.quaternion_order([1, i, (i+j)/2, (1+k)/2])
    
    Z_p = Integers(p)
    Bound = 1000

    qf = BinaryQF([1, 0, 1])

    omega = None

    # Chosen torsion from precompute
    fullG = 3 * 5^3 * 7^2 * 11^3 * 13^2 * 17 * 19 * 23 * 29 * 31 * 37 * 41 * 43 * 47 * 53 * 61 * 71 * 73 * 79 * 89 * 97 * 137 * 139 * 151 * 163 * 167 * 181 * 193 * 199 * 223 * 239 * 241 * 257 * 281 * 311 * 317 * 331 * 349 * 353 * 367 * 373 * 397 * 401 * 409 * 419 * 421 * 433 * 457 * 461 * 463 * 487 * 499 * 509 * 541 * 547 * 569 * 571 * 577 * 587 * 593 * 617 * 619 * 631 * 641 * 659 * 691 * 719 * 727 * 739 * 743 * 751 * 757 * 761 * 773 * 787 * 797 * 827 * 829 * 853 * 857 * 863 * 881 * 941 * 953 * 967 * 971

    for div in range(1, 100000):
        if (fullG) % div != 0:
            continue
        path = (fullG)/div
        n = d*(conductor*path)**2
        print(RR(log(n, p)))
        rts = [ZZ(r) for r in Z_p(n).sqrt(all=True)]
        for x in rts:
            M = ZZ((n - x**2)/p)
            print(M)
            y, z, sol = Cornacchia(qf, M)
            if sol:
                print("We got it!!!")
                omega = x*i + y*j + z*k
                break
        if omega:
            break
    
    I = O0*path + O0*omega


    print(I)
    print(RR(log(I.norm(), p)))
    assert omega/path in I.right_order()
    print(omega/path)
    