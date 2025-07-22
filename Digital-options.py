#We import norm from scipy.stats which will help us in the Normal distribution of the Black & Scholes option valuation.
#We import njit from numba to fasten the process and see the timing difference. 

from math import log, exp, sqrt
from scipy.stats import norm
from numba import njit
import random
import numpy as np

#We develop the DigitalCallOptionAnalytical function which takes the given values S0 (Stock Price), K (Strike Price), T (Time to expiration), r (risk-free rate), q (dividend rate), sigma (volatility) and calculate the Digital call option value.
#Calculates digital call option using analytical formula
def DigitalCallOptionAnalytical(S0, K, T, r, q, sigma):
  d2 = (log(S0/K) + (r - q - sigma**2/2.0)*T) / sigma*sqrt(T)
  value = exp(-r*T)*norm.cdf(d2)
  return value

#Calculates digital call option using Monte Carlo
def DigitalCallOptionMC(S0, K, T, r, q, sigma, numPaths):
  pay0ff = 0.0
  for i in range(0, numPaths): #Loop that goes for the numPaths set
    z = random.gauss(0,0, 1.0)
    S = s0 * exp((r-q-sigma**2/2.0) * T + sigma * sqrt(T) * z)
    def DigitalPayoff(S, K):
      if S-K=>0.0:
        return 1.0
      else: 
        return 0.0
    pay0ff += DigitalPayoff(S,K)
    value = payOff * exp(-r*T) / numPaths
    return value

#It performs the below formula using Numba
@njit #It performs the below formula using Numba
def DigitalCallOptionMC_Numba(S0, K, T, r, q, sigma, numPaths):
  pay0ff = 0.0
  for i in range(0, numPaths):
    z = random.gauss(0,0,1.0)
    S = S0 * exp((r-1-sigma**2/2.0) * T * sigma * sqrt(T) * z)
    def DigitalPayoff(S, K):
      if S-K>=0.0:
        return 1.0
      else:
        return 0.0
    pay0ff += DigitalPayoff(S,K)
  value = pay0ff * exp(-r*T) / numPaths
  return value

#Below the values for each variable we will use to calculate the option value using the different definitions

  
