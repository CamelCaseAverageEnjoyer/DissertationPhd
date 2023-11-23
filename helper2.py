import numpy as np
import matplotlib.pyplot as plt

def if_prime_number(a: int):
    k = 0
    for i in range(2, a // 2 + 1):
        if a % i == 0:
            k = k + 1
    return k <= 0

def get_divisors(a: int):
    k = 0
    for i in range(2, a // 2 + 1):
        if a % i == 0:
            k = k + 1
    return k


x = [i+1 for i in range(100)]
# y = [int(if_prime_number(i)) for i in x]
y = [get_divisors(i) for i in x]
plt.plot(x, y)
'''for i in range(len(x)):
    if y[i] > 0:
        plt.scatter(x[i], y[i], c='navy')'''
plt.xlabel("Года")
plt.show()
