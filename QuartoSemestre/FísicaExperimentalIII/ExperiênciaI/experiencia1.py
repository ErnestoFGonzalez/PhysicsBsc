ms = [182.36*(10**-3), 232.36*(10**-3), 282.37*(10**-3)]
Ls = [55*(10**-2), 60*(10**-2), 70*(10**-2)]
densidade_linear = 3.6*(10**-4)
ns = [1,2,3,4]

for m in ms:
    for L in Ls:
        for n in ns:
            f_n = (n/(2*L)) * ((m*9.81)/(densidade_linear))**0.5
            print("n={}\tm={}\tL={}\tf={}\t\tT={} ms".format(n,m,L,f_n,(f_n**-1)*(10**3)))
