#!/usr/bin/env python
# coding: utf-8




import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt


def American_option(K, S, T, r, sigma, step, Type):
    K = float(K)
    S = float(S)
    T = float(T)
    r = float(r)
    sigma = float(sigma)
    step = int(step)
    dt = T / step  # delta_t
    DF = np.e ** (-r * dt)  # discount factor
    up = np.e ** (sigma * np.sqrt(dt))  # u
    down = 1 / up  # d
    p = (np.e ** (r * dt) - down) / (up - down)  # probability
    stock = np.zeros((step + 1, step + 1))  # record asset price
    option = np.zeros((step + 1, step + 1))  # record option price
    stock[0][0] = S  # initialize asset price

    # according to related position, calculate stock price
    for i in range(1, step + 1):
        for j in range(i):
            stock[j][i] = stock[j][i - 1] * up
            stock[j + 1][i] = stock[j][i - 1] * down
    # Put
    if Type == 'Put':
        for i in range(step + 1):  # calculate option price at maturity
            option[i][step] = max(K - stock[i][step], 0)
        inverse = step

        # backforward calculate option price at each step
        while (inverse > 0):
            for i in range(inverse):
                price = DF * (p * option[i][inverse] + (1 - p) * option[i + 1][inverse])
                option[i][inverse - 1] = max((K - stock[i][inverse - 1]), price)
            inverse -= 1
    # American Call is the same as European
    else:
        for i in range(step + 1):  # calculate option price at maturity
            option[i][step] = max(stock[i][step] - K, 0)
        inverse = step

        # backforward calculate option price at each step
        while (inverse > 0):
            for i in range(inverse):
                price = DF * (p * option[i][inverse] + (1 - p) * option[i + 1][inverse])
                option[i][inverse - 1] = price
            inverse -= 1
    return option[0][0]