# Gets the data
import pandas as pd
from datetime import datetime

# Date of retrieval
date_format = "%d/%m/%Y"
ref = datetime.strptime('13/4/2021', date_format)

# Read the data. Columns are:
# SPOT PRICE, MATURITY, STRIKE PRICE, BID, ASK, IV
df = pd.read_csv('../data/NIFTYDATA.csv')

# Store Maturity as number of remaining days instead of date
for i in range(len(df["MATURITY"])):
    df["MATURITY"][i] = (datetime.strptime(df["MATURITY"][i],date_format) - ref).days

print(df)

