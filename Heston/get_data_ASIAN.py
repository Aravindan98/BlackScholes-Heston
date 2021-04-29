# Gets the data
from heston_model_Euler import HestonProcess
import math
import pandas as pd
import numpy as np
from datetime import datetime
from scipy.optimize import differential_evolution

df = 0

def read_data():
    global df
    # Date of retrieval
    date_format = "%d/%m/%Y"
    ref = datetime.strptime('20/5/2020', date_format)

    # Read the data. Columns are:
    # SPOT PRICE, MATURITY, STRIKE PRICE, BID, ASK, IV
    df = pd.read_csv('./data/ASIANdata.csv')

    # Store Maturity as number of remaining days instead of date
    for i in range(len(df["MATURITY"])):
        df["MATURITY"][i] = (datetime.strptime(df["MATURITY"][i],date_format) - ref).days

def marketPrices(df):
    a = [df["PREMIUM"][i] for i in range(len(df))]
    return np.array(a);

def evaluateParams(x):
    print("Evaluating...",flush=True)
    [v0,vBar,k,eta,rho] = x
    f = 0.0
    for i in range(len(df)):
        t = 0
        r = 0.0575
        S0 = df["SPOT PRICE"][i]
        V0 = 0
        K = df["STRIKE PRICE"][i]
        T = df["MATURITY"][i]
        n = T
        numSamples = 1000
        actualPrice = df["PREMIUM"][i]
        hestonPrice = HestonProcess(t,S0,V0,K,T,r,k,vBar,eta,rho,n,numSamples)
        # print("K T A P")
        # print(K)
        # print(T)
        print(actualPrice)
        print(hestonPrice)
        print("---")
        f += (actualPrice - hestonPrice)**2
    f /= len(df)
    f = math.sqrt(f)
    return f

if __name__ == '__main__':
    read_data()
    # bounds = ([0,0,0,0,-1],[1,1,np.inf,5,1])
    # bounds = [(0,1),(0,1),(0,1),(0,5),(-1,1)]
    # result = differential_evolution(evaluateParams,bounds,popsize=10,workers=-1,maxiter=20,disp=True,)
    # print(result)
    # x = [0.73887546,  0.00410743,  0.60635024,  0.16099308, -0.5377229]
    x = [0.01126412,  0.00838991,  0.00614512,  0.06720342, -0.75]
    # x = [0.01,0.005,0.02,0.05,-0.75]
    print(evaluateParams(x))
