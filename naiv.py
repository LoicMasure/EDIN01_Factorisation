# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 16:52:21 2016
@author: loicm
"""
import numpy as np
import itertools


def check_prime(n, primes):
    for p in primes:
        if not n % p:
            return False
    return True


def prime_sieve():
    primes = set()
    for n in itertools.count(2):
        if check_prime(n, primes):
            primes.add(n)
            yield n


def generate_primes(L):
    """
    Generates a list of the L first prime numbers .
    """
    # We need to compute the Bound of the factor set.
    i = 0
    list_p = []
    for p in prime_sieve():
        i += 1
        list_p.append(p)
        if i >= L:
            break
    return list_p


def decompose_opt(r, list_p):
    """
    Decompose naïvely r as a product of factors in list_p and returns an array
    with the powers corresponding to the factors.
    """
#    print("Decompose pour {}".format(r))
    res = np.zeros(len(list_p), dtype=np.int)
    i = 0
    while i < len(list_p) and r != 1:
        if r % list_p[i] == 0:
            res[i] += 1
            r //= list_p[i]
        else:
            i += 1
    return res, (r < list_p[-1])

def recompose(x, list_p, N):
    """
    Reconstruct the integer decomposed with the prime factor basis given in
    argument.
    """
    res = 1
    for i in zip(x, list_p):
        plus = 1
        for j in range(i[0]):
            plus *= i[1]
            plus %= N
        res *= plus
        res %= N
    return int(res % N)
