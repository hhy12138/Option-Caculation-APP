#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np
import math
from scipy.stats import norm
import random


# In[8]:


# Calculate the Geometric basket options
# S1, s2 are asset prices at time zero, 
# sigma is volatility
# r is risk free rate
# T is time from now to the maturity data, the unit is day
# K is strike price
# rho is correlation 
# option_type is the option type, call or put
def Geo_Basket_option(S1, S2, sigma1, sigma2, r, T, K, rho, option_type):
    sigma_B_g = math.sqrt(sigma1** 2 + 2 * sigma1 * sigma2 * rho + sigma2** 2) / 2
    mu_B_g = r - 0.5 * (sigma1** 2 + sigma2** 2) / 2 + 0.5 * sigma_B_g** 2
    B_g_0 = (S1 * S2)** (1/2)
    d1_hat = (math.log(B_g_0 / K) + (mu_B_g + 0.5 * sigma_B_g** 2)* T) / (sigma_B_g * math.sqrt(T))
    d2_hat = d1_hat - sigma_B_g * math.sqrt(T)
    if(option_type == 'Call'):
        return math.exp(-r * T) * (B_g_0 * math.exp(mu_B_g * T) * norm.cdf(d1_hat) - K * norm.cdf(d2_hat))
    elif(option_type == 'Put'):
        return math.exp(-r * T) * (K * norm.cdf(-d2_hat) - B_g_0 * math.exp(mu_B_g * T) * norm.cdf(-d1_hat))


# In[9]:


# Test Geometric basket options
if __name__=='__main__':
    c1 = {'S1' : 100, 'S2' : 100, 'K' : 100, 'sigma1' : 0.3, 'sigma2' : 0.3, 'rho': 0.5, 'option' : 'Put'}
    c2 = {'S1' : 100, 'S2' : 100, 'K' : 100, 'sigma1' : 0.3, 'sigma2' : 0.3, 'rho': 0.9, 'option' : 'Put'}
    c3 = {'S1' : 100, 'S2' : 100, 'K' : 100, 'sigma1' : 0.3, 'sigma2' : 0.3, 'rho': 0.5, 'option' : 'Call'}
    c5 = {'S1' : 100, 'S2' : 100, 'K' : 100, 'sigma1' : 0.3, 'sigma2' : 0.3, 'rho': 0.9, 'option' : 'Call'}
    bucket = [c1, c2, c3, c5]
    for item in bucket:
        print(Geo_Basket_option(item['S1'], item['S2'], item['sigma1'], item['sigma2'], 0.05, 3, item['K'], item['rho'], item['option']))


# In[10]:


# Calculate the Arithmetic basket options
# S1, s2 are asset prices at time zero, 
# sigma is volatility
# r is risk free rate
# T is time from now to the maturity data, the unit is day
# K is strike price
# rho is correlation 
# option_type is the option type, call or put
# path is number of simulations
# ctr_method is control variate method
def Ari_Basket_option(S1, S2, sigma1, sigma2, r, T, K, rho, option_type, path, ctr_method = None):
    geo_basket = Geo_Basket_option(S1, S2, sigma1, sigma2, r, T, K, rho, option_type)
    drift1 = math.exp((r - 0.5 * sigma1** 2) * T)
    drift2 = math.exp((r - 0.5 * sigma2** 2) * T)
    bucket = []
    random.seed(1)
    for i in range(path):
        sample = []
        Z1 = random.gauss(0,1)
        Y = random.gauss(0,1)
        Z2 = rho * Z1 + math.sqrt(1 - rho** 2) * Y
        growthFactor1 = drift1 * math.exp(sigma1 * math.sqrt(T) * Z1)
        growthFactor2 = drift2 * math.exp(sigma2 * math.sqrt(T) * Z2)
        sample.append(S1 * growthFactor1)
        sample.append(S2 * growthFactor2)
        bucket.append(sample)
    
    Bucket = np.array(bucket)
    
    #Arithmetic Mean
    arithMean = np.mean(Bucket, axis = 1)
    if(option_type == 'Call'):
        arithPayOff = math.exp(-r * T) * np.maximum(arithMean - K , 0)
    elif(option_type == 'Put'):
        arithPayOff = math.exp(-r * T) * np.maximum(K - arithMean , 0)
    
    #Geomrteic Mean
    geoMean = np.exp(0.5 * np.sum(np.log(Bucket), axis = 1))
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
    Z = arithPayOff + theta * (geo_basket - geoPayOff)
    Zmean = np.mean(Z)
    Zstd = np.std(Z)
    confcv = [Zmean - 1.96 * Zstd/math.sqrt(path), Zmean + 1.96 * Zstd/math.sqrt(path)]
    
    if(ctr_method == None):
        return Pmean, confmc, Geo_mean, confmc_geo
    return Zmean, confcv, Geo_mean, confmc_geo


# In[11]:

if __name__=='__main__':
    for item in bucket:
        value1, interval1, value2, interval2 = Ari_Basket_option(item['S1'], item['S2'], item['sigma1'], item['sigma2'], 0.05, 3, item['K'], item['rho'], item['option'],100000)
        print('Atirthmetic Asian with Monte Carlo {}'.format(value1))
        print(interval1)
        print('Geometric Asian with Monte Carlo {}'.format(value2))
        print(interval2)
        print('')


    # In[12]:


    # Test Geometric basket options
    cc1 = {'S1' : 100, 'S2' : 100, 'K' : 100, 'sigma1' : 0.3, 'sigma2' : 0.3, 'rho': 0.5, 'option' : 'Put', 'method':'ctr'}
    cc2 = {'S1' : 100, 'S2' : 100, 'K' : 100, 'sigma1' : 0.3, 'sigma2' : 0.3, 'rho': 0.9, 'option' : 'Put', 'method':'ctr'}
    cc3 = {'S1' : 100, 'S2' : 100, 'K' : 100, 'sigma1' : 0.3, 'sigma2' : 0.3, 'rho': 0.5, 'option' : 'Call', 'method':'ctr'}
    cc4 = {'S1' : 100, 'S2' : 100, 'K' : 100, 'sigma1' : 0.3, 'sigma2' : 0.3, 'rho': 0.9, 'option' : 'Call', 'method':'ctr'}
    check = [cc1,cc2,cc3,cc4]
    for item in check:
        value1, interval1, value2, interval2 = Ari_Basket_option(item['S1'], item['S2'], item['sigma1'], item['sigma2'], 0.05, 3, item['K'], item['rho'], item['option'],100000, item['method'])
        print('Atirthmetic Asian with control variate {}'.format(value1))
        print(interval1)
        print('Geometric Asian with Monte Carlo {}'.format(value2))
        print(interval2)
        print('')


    # In[ ]:




