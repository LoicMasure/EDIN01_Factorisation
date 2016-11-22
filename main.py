# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 14:54:02 2016
@author: loicm
Problems for now: the higher is N, the smaller the factor base has to be.
Otherwise, the computation of (x, y) such that x^2=y^2 % N is wrong...
"""
import math
import numpy as np
import naiv
import os
from time import sleep, perf_counter as pc
t0 = pc()
N = 392742364277          # The number to factorize
L = 1000  # The number of relations we want to have to compute x and y
L_more = 10  # The number or additional rows compare to L
list_p = naiv.generate_primes(L)

print("B= {}".format(list_p[-1]))


# We need to generate a list of numbers whose the fators are B-smooth.
m = 2
j = m
k = 0
list_r2 = []
list_all = []
while len(list_all) < L+L_more:
    r = math.floor(math.sqrt(k * N)) + j
    r2 = r**2 % N
    # We compute naïvely the prime factors up to B
    row, b = naiv.decompose_opt(r2, list_p)
#    verif = naiv.recompose(row, list_p, N)
#    print("Vérification: r={},  r2={}, verif={}".format(r, r2, verif))
    if b and r2 not in list_r2:
        test = False
        row %= 2
        for i in range(len(list_all)):
            if (row == list_all[i][2]).all():
                test = True
                break
        if not test:
            list_r2.append(r2)  # For unicity
            list_all.append((r, r2, row))
            # if len(list_all) % 20 == 0:
                # print(len(list_all))
#            print("Ajout de {}, ligne={}".format(r, row))
    # Increment in diagonals, is arbitrary.
    if k == m:
        m += 1
        j = m
        k = 0
    else:
        j -= 1
        k += 1
print("Fin de la boucle.")
list_all = sorted(list_all, key=lambda tutuple: tutuple[0])
list_r, list_r2, list_factor = zip(*list_all)
# print(list_r, list_r2)
for r, r2 in zip(list_r, list_r2):
    assert(r2 == r**2 % N)
M = np.array(list_factor, dtype=np.int)

header = "{} {}".format(L + L_more, L)
np.savetxt("essai", M, fmt="%d")
with open("essai", 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(header.rstrip('\r\n') + '\n' + content)
# We call GaussBin to compute the solutions to the binary linear system.
os.system("GaussBin.exe essai sortie")

# We load the solution in the program.
X = np.loadtxt("sortie", dtype=np.int, skiprows=1)
for x in range(X.shape[0]):
    # We get the x such that x^2=y^2
    x1 = naiv.recompose(X[x, :], list_r, N)

    # We get the y such that x^2=y^2
    row = np.zeros(L)
    for i, c in enumerate(X[x]):
        if c == 1:
            plus, b = naiv.decompose_opt(list_r2[i], list_p)
            if b:
                row += plus
            else:
                raise Exception("Bad decomposition: {} in {} gives {}".format(list_r[i],list_p, row))
    assert((row % 2 == 0).all())
    row //= 2
    row = row.astype(np.int)
    x2 = naiv.recompose(row, list_p, N)
    print(x)
    print("Solution possible: x = {} y = {}, x2 = {}, y2 = {}".format(x1,
                                                                      x2,
                                                                      x1**2%N,
                                                                      x2**2%N))
    assert((x1 ** 2) % N == (x2 ** 2) % N)

    d = math.gcd(max(x1, x2) - min(x1, x2), N)
    if d != N and d != 1:
        p = d
        q = int(N/d)
        break

print("Final solution: N = {}, p = {}, q = {}".format(N, p, q))
assert (N == p*q)
print("Execution time: {}".format(pc()-t0))
