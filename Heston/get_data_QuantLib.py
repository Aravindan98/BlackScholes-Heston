# Gets the data
import QuantLib as ql
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


def evaluateParams(x):
    f = 0
    [v0,theta,kappa,sigma,rho] = x
    for i in range(len(df)):
        strike_price  = float(df["STRIKE PRICE"][i])
        spot_price = float(df["SPOT PRICE"][i])
        maturity_date = ql.Date(25, 6, 2021)
        payoff = ql.PlainVanillaPayoff(ql.Option.Call, strike_price)

        dividend_rate =  0.0575
        option_type = ql.Option.Call

        risk_free_rate = 0.0575
        day_count = ql.Actual365Fixed()
        calendar = ql.UnitedStates()

        calculation_date = ql.Date(20, 5, 2021)
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
        # print(df["PREMIUM"][i])
        # print(h_price)
        # print("---")
        f += (df["PREMIUM"][i] - h_price)**2;
    f /= len(df)
    f = math.sqrt(f)
    return f

if __name__ == '__main__':
    read_data()
    bounds = [(0.01,1),(0.01,1),(0.01,1),(0,1),(-1,1)]
    result = differential_evolution(evaluateParams,bounds,popsize=300,workers=-1,maxiter=500,disp=True,)
    print(result)
    # x=[0.15942742, 0.85085073, 0.09458799, 0.62937492, 0.95141181]
    # print(evaluateParams(x))

# ASIAN:::
# DE:
# [0.13034587, 0.83903986, 0.57130949, 0.98225277, 0.9937405]
#
