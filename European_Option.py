#!/usr/bin/env python
# coding: utf-8
from scipy.stats import norm
import math
def d1(stock,strike,T,sigma,r):
    return (math.log(stock/strike,math.e)+(r+0.5*math.pow(sigma,2))*(T))/(sigma*math.sqrt(T))

def European_option(stock,strike,T,sigma,r,type):
    D1 = d1(stock, strike, T, sigma, r)
    D2 = D1 - sigma * math.sqrt(T)
    price = stock * norm.cdf(D1) - strike * math.exp(-r * (T)) * norm.cdf(D2)
    if type=='Put':
        price = price+strike*math.exp(-r*(T))-stock
    return price

def Implied_volatility(stock,strike,T,r,true_price,type):
    threshold=1e-8
    nmax=5000
    sigma_hat=math.sqrt(2*abs((math.log(stock/strike)+(r)*T)/(T)))
    sigma=sigma_hat
    sigmadiff=1
    n=1
    judge=0
    while sigmadiff > threshold and n<nmax:
        price = European_option(stock,strike,T,sigma,r,type)
        D1=d1(stock,strike,T,sigma,r)
        temp=-(math.log(stock/strike)+(r)*T)/(sigma*sigma*math.sqrt(T))+0.5*math.sqrt(T)
        vega=stock*math.sqrt(T)
        vega=vega*(1/math.sqrt(2*math.pi))*math.exp(-0.5*(D1*D1))*(-D1)*temp
        increment=(price-true_price)/(vega+10.0)
        sigma=sigma-increment
        n=n+1
        sigmadiff=abs(increment)
        if sigmadiff<=threshold:
            judge=1
    if judge==0:
        return "nan"
    return sigma

if __name__=="__main__":
    print(European_option(100,120,0.5,0.2,0.01,'Put'))
    print(Implied_volatility(70.0,71.0,40/365,0.2,1.2,'Call'))







