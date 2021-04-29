# Gets the data
import QuantLib as ql
import math
import pandas as pd
import numpy as np
from datetime import datetime
from scipy.optimize import differential_evolution

df1 = 0
df2 = 0
df3 = 0

def read_data():
    global df1, df2, df3
    # Read the data. Columns are:
    # SPOT PRICE, MATURITY, STRIKE PRICE, BID, ASK, IV
    df1 = pd.read_csv('./data/NIFTY1.csv')
    df2 = pd.read_csv('./data/NIFTY2.csv')
    df3 = pd.read_csv('./data/NIFTY3.csv')

def evaluateParams(x):
    # print("Evaluating...", flush=True)
    f = 0
    [v0,theta,kappa,sigma,rho] = x
    for i in range(len(df1)):
        strike_price  = float(df1["STRIKE PRICE"][i])
        spot_price = float(df1["SPOT PRICE"][i])
        maturity_date = ql.Date(29, 4, 2021)
        payoff = ql.PlainVanillaPayoff(ql.Option.Call, strike_price)
        dividend_rate =  0.0575
        option_type = ql.Option.Call
        risk_free_rate = 0.0575
        day_count = ql.Actual365Fixed()
        calendar = ql.UnitedStates()
        calculation_date = ql.Date(13, 4, 2021)
        ql.Settings.instance().evaluationDate = calculation_date
        # construct the European Option
        payoff = ql.PlainVanillaPayoff(option_type, strike_price)
        exercise = ql.EuropeanExercise(maturity_date)
        european_option = ql.VanillaOption(payoff, exercise)
        # construct the Heston process
        spot_handle = ql.QuoteHandle(
            ql.SimpleQuote(spot_price)
            )
        flat_ts = ql.YieldTermStructureHandle(
            ql.FlatForward(calculation_date, risk_free_rate, day_count)
            )
        dividend_yield = ql.YieldTermStructureHandle(
            ql.FlatForward(calculation_date, dividend_rate, day_count)
            )
        heston_process = ql.HestonProcess(flat_ts,
            dividend_yield,
            spot_handle,
            v0,
            kappa,
            theta,
            sigma,
            rho)
        engine = ql.AnalyticHestonEngine(ql.HestonModel(heston_process),0.01, 100000)
        european_option.setPricingEngine(engine)
        h_price = european_option.NPV()
        f += ((df1["BID"][i]+df1["ASK"][i])*0.5 - h_price)**2
    for i in range(len(df2)):
        strike_price  = float(df2["STRIKE PRICE"][i])
        spot_price = float(df2["SPOT PRICE"][i])
        maturity_date = ql.Date(27, 5, 2021)
        payoff = ql.PlainVanillaPayoff(ql.Option.Call, strike_price)
        dividend_rate =  0.0575
        option_type = ql.Option.Call
        risk_free_rate = 0.0575
        day_count = ql.Actual365Fixed()
        calendar = ql.UnitedStates()
        calculation_date = ql.Date(13, 4, 2021)
        ql.Settings.instance().evaluationDate = calculation_date
        # construct the European Option
        payoff = ql.PlainVanillaPayoff(option_type, strike_price)
        exercise = ql.EuropeanExercise(maturity_date)
        european_option = ql.VanillaOption(payoff, exercise)
        # construct the Heston process
        spot_handle = ql.QuoteHandle(
            ql.SimpleQuote(spot_price)
            )
        flat_ts = ql.YieldTermStructureHandle(
            ql.FlatForward(calculation_date, risk_free_rate, day_count)
            )
        dividend_yield = ql.YieldTermStructureHandle(
            ql.FlatForward(calculation_date, dividend_rate, day_count)
            )
        heston_process = ql.HestonProcess(flat_ts,
            dividend_yield,
            spot_handle,
            v0,
            kappa,
            theta,
            sigma,
            rho)
        engine = ql.AnalyticHestonEngine(ql.HestonModel(heston_process),0.01, 100000)
        european_option.setPricingEngine(engine)
        h_price = european_option.NPV()
        f += ((df2["BID"][i]+df2["ASK"][i])*0.5 - h_price)**2
    for i in range(len(df3)):
        strike_price  = float(df3["STRIKE PRICE"][i])
        spot_price = float(df3["SPOT PRICE"][i])
        maturity_date = ql.Date(22, 4, 2021)
        payoff = ql.PlainVanillaPayoff(ql.Option.Call, strike_price)
        dividend_rate =  0.0575
        option_type = ql.Option.Call
        risk_free_rate = 0.0575
        day_count = ql.Actual365Fixed()
        calendar = ql.UnitedStates()
        calculation_date = ql.Date(13, 4, 2021)
        ql.Settings.instance().evaluationDate = calculation_date
        # construct the European Option
        payoff = ql.PlainVanillaPayoff(option_type, strike_price)
        exercise = ql.EuropeanExercise(maturity_date)
        european_option = ql.VanillaOption(payoff, exercise)
        # construct the Heston process
        spot_handle = ql.QuoteHandle(
            ql.SimpleQuote(spot_price)
            )
        flat_ts = ql.YieldTermStructureHandle(
            ql.FlatForward(calculation_date, risk_free_rate, day_count)
            )
        dividend_yield = ql.YieldTermStructureHandle(
            ql.FlatForward(calculation_date, dividend_rate, day_count)
            )
        heston_process = ql.HestonProcess(flat_ts,
            dividend_yield,
            spot_handle,
            v0,
            kappa,
            theta,
            sigma,
            rho)
        engine = ql.AnalyticHestonEngine(ql.HestonModel(heston_process),0.01, 100000)
        european_option.setPricingEngine(engine)
        h_price = european_option.NPV()
        f += ((df3["BID"][i]+df3["ASK"][i])*0.5 - h_price)**2
    f /= len(df1) + len(df2) + len(df3)
    f = math.sqrt(f)
    return f

if __name__ == '__main__':
    read_data()
    bounds = [(0.01,1),(0.01,1),(0.01,1),(0,1),(-1,1)]
    result = differential_evolution(evaluateParams,bounds,popsize=300,workers=-1,maxiter=500,disp=True,)
    print(result)
    # x=[0.15942742, 0.85085073, 0.09458799, 0.62937492, 0.95141181]
    # print(evaluateParams(x))

# NIFTY:::
# DE:
# [0.04624947,  0.41215684,  0.60967619,  0.81901855, -0.99993523]
# 

