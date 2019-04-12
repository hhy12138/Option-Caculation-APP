#!/usr/bin/env python
# coding: utf-8

# In[16]:


import numpy as np
import math
from scipy.stats import norm
import random


# In[17]:


# Calculate the Geometric Asian Option
# S is asset price at time zero, 
# sigma is volatility
# r is risk free rate
# T is time from now to the maturity data, the unit is day
# K is strike price
# n is number of times for geometric average 
# option_type is the option type, call or put
def Geometric_Asian_Option(S, sigma, r, T, K, n, option_type):
    sigma_hat = sigma * math.sqrt((n + 1) * (2 * n + 1)/ (6 * n** 2))
    mu_hat = (r - 0.5 * sigma** 2)*(n + 1)/(2*n) + 0.5 * sigma_hat** 2
    d1_hat = (math.log(S / K) + (mu_hat + 0.5 * sigma_hat** 2) * T)/(sigma_hat * math.sqrt(T))
    d2_hat = d1_hat - sigma_hat * math.sqrt(T)
    if option_type == 'Call':
        return math.exp(-r * T) * (S * math.exp(mu_hat * T) * norm.cdf(d1_hat) - K * norm.cdf(d2_hat))
    elif option_type == 'Put':
        return math.exp(-r * T) * (K * norm.cdf(-d2_hat) - S * math.exp(mu_hat * T) * norm.cdf(-d1_hat))


# In[18]:


# Geometric_Asian_option Test
if __name__=='__main__':
    r = 0.05
    T = 3
    S = 100
    case1 = {'sigma' : 0.3, 'K' : 100, 'n' : 50, 'Type' : 'Put'}
    case2 = {'sigma' : 0.3, 'K' : 100, 'n' : 100, 'Type' : 'Put'}
    case3 = {'sigma' : 0.4, 'K' : 100, 'n' : 50, 'Type' : 'Put'}
    case4 = {'sigma' : 0.3, 'K' : 100, 'n' : 50, 'Type' : 'Call'}
    case5 = {'sigma' : 0.3, 'K' : 100, 'n' : 100, 'Type' : 'Call'}
    case6 = {'sigma' : 0.4, 'K' : 100, 'n' : 50, 'Type' : 'Call'}
    check_list = [case1, case2, case3, case4, case5, case6]
    for item in check_list:
        print(Geometric_Asian_Option(S, item['sigma'], r, T, item['K'], item['n'], item['Type']))


# In[19]:


# Calculate the Arithmetic Asian Option
# S is asset price at time zero, 
# sigma is volatility
# r is risk free rate
# T is time from now to the maturity data, the unit is day
# K is strike price
# n is number of times for geometric average 
# option_type is the option type, call or put
# path is number of paths in monte carlo simulation
# ctr_method is to specify the control_varite method
def Arithmetic_Asian_Option(S, sigma, r, T, K, n, option_type, path, ctr_method = None):
    geo = Geometric_Asian_Option(S, sigma, r, T, K, n, option_type)
    Dt = T / n
    drift = math.exp((r - 0.5 * sigma** 2) * Dt)
    random.seed(1)
    spath = []
    for i in range(path):
        sample = []
        growthFactor = drift * math.exp(sigma * math.sqrt(Dt) * random.gauss(0,1))
        sample.append(S * growthFactor)
        for j in range(1, n):
            growthFactor = drift * math.exp(sigma * math.sqrt(Dt) * random.gauss(0,1))
            sample.append(sample[j-1] * growthFactor)
        spath.append(sample)
    
    Spath = np.array(spath)
    
    #Arithmetic Mean
    arithMean = np.mean(Spath, axis = 1)
    if(option_type == 'Call'):
        arithPayOff = math.exp(-r * T) * np.maximum(arithMean - K , 0)
    elif(option_type == 'Put'):
        arithPayOff = math.exp(-r * T) * np.maximum(K - arithMean , 0)
    #Geomrteic Mean
    geoMean = np.exp((1/n) * np.sum(np.log(Spath), axis = 1))
    if(option_type == 'Call'):
        geoPayOff = math.exp(-r * T) * np.maximum(geoMean - K , 0)
    elif(option_type == 'Put'):
        geoPayOff = math.exp(-r * T) * np.maximum(K - geoMean , 0)
        
    #Standard Monte Carlo for Asian Geometric
    Geo_mean = np.mean(geoPayOff)
    Geo_std = np.std(geoPayOff)
    confmc_geo = [Geo_mean - 1.96 * Geo_std/math.sqrt(path), Geo_mean + 1.96 * Geo_std/math.sqrt(path)]
    
    #Standard Monte Carlo
    Pmean = np.mean(arithPayOff)
    Pstd = np.std(arithPayOff)
    confmc = [Pmean - 1.96 * Pstd/math.sqrt(path), Pmean + 1.96 * Pstd/math.sqrt(path)]
    
    # control variate
    covXY = np.mean(arithPayOff * geoPayOff) - np.mean(arithPayOff) * np.mean(geoPayOff)
    theta = covXY / (np.var(geoPayOff))
    
    #control varite version
    Z = arithPayOff + theta * (geo - geoPayOff)
    Zmean = np.mean(Z)
    Zstd = np.std(Z)
    confcv = [Zmean - 1.96 * Zstd/math.sqrt(path), Zmean + 1.96 * Zstd/math.sqrt(path)]
    
    if(ctr_method == None):
        return Pmean, confmc, Geo_mean, confmc_geo
    return Zmean, confcv, Geo_mean, confmc_geo


# In[20]:


#Test arithmetic with control variate
if __name__=='__main__':
    c1 = {'sigma' : 0.3, 'K' : 100, 'n' : 50, 'Type' : 'Put', 'ctr_method': 'control variate'}
    c2 = {'sigma' : 0.3, 'K' : 100, 'n' : 100, 'Type' : 'Put', 'ctr_method': 'control variate'}
    c3 = {'sigma' : 0.4, 'K' : 100, 'n' : 50, 'Type' : 'Put', 'ctr_method': 'control variate'}
    c4 = {'sigma' : 0.3, 'K' : 100, 'n' : 50, 'Type' : 'Call', 'ctr_method': 'control variate'}
    c5 = {'sigma' : 0.3, 'K' : 100, 'n' : 100, 'Type' : 'Call', 'ctr_method': 'control variate'}
    c6 = {'sigma' : 0.4, 'K' : 100, 'n' : 50, 'Type' : 'Call', 'ctr_method': 'control variate'}
    clist = [c1,c2,c3,c4,c5,c6]
    for item in clist:
        value1, interval1, value2, interval2 = Arithmetic_Asian_Option(S, item['sigma'], r, T, item['K'], item['n'], item['Type'], 100000, item['ctr_method'])
        print('Atirthmetic Asian with control variate {}'.format(value1))
        print(interval1)
        print('Geometric Asian with Monte Carlo {}'.format(value2))
        print(interval2)
        print('')


# In[21]:


#Test arithmetic with no control variate
    for item in check_list:
        value1, interval1, value2, interval2 = Arithmetic_Asian_Option(S, item['sigma'], r, T, item['K'], item['n'], item['Type'],100000)
        print('Atirthmetic Asian with Monte Carlo {}'.format(value1))
        print(interval1)
        print('Geometric Asian with Monte Carlo {}'.format(value2))
        print(interval2)
        print('')

