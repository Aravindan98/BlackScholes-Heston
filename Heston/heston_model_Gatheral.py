# Implements the Heston Model using the discretization from Gatheral
import numpy as np
import math

# t: Date at which to compute process
# S: Price of asset at time t
# V: Volatility at time t
# K: Strike price
# T: Maturity
# r: risk free rate
# k : Speed of mean reversion
# vBar : Long term variance mean
# eta : Volatility of volatility
# rho: Relation between Brownian processes
def HestonProcess(t,S0,V0,K,T,r,k,vBar,eta,rho,n,numSamples):
    samples = np.zeros(numSamples)
    for i in range(0,len(samples)):
        samples[i] = simulate(t,r,V0,T,k,vBar,eta,rho,n)
    samples = S0 * np.exp(samples)
    samples = samples-K
    samples = np.clip(samples,0,None)
    price = np.mean(samples) * math.exp(-r*(T-t))
    return price

def simulate(t,r,V0,T,k,vBar,eta,rho,n):
    x = 0
    v = V0
    deltaT = (T-t)/n;
    for i in range(1,n):
        x,v = step(r,x,v,k,vBar,eta,rho,deltaT)
    return x

def step(r,x,v,k,vBar,eta,rho,deltaT):
    cov = np.array([[1,rho],[rho,1]])
    mu = np.array([0,0])
    z,w = np.random.multivariate_normal(mu,cov).T
    nv = (math.sqrt(v) + 0.5*eta*math.sqrt(deltaT)*z)**2 - k*(v - vBar)*deltaT - 0.25*eta*eta*deltaT
    if (nv < 0):
        nv = 0
    nx = x - 0.5*v*deltaT + math.sqrt(v*deltaT)*w
    return (nx,nv)
