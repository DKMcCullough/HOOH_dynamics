#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 16:05:08 2023

location: '/Users/dkm/Documents/Talmy_research/Zinser_lab/Projects/Monocultures/src/Monocultures.py'


@author: dkm
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import odeint


xl = pd.ExcelFile('../data/Monocultures.xlsx')

strains = xl.sheet_names  # see all sheet names
print(strains)


for s in strains: 
    df = xl.parse(strains[s])  # read a specific sheet to DataFrame
    print(df)

print('Done')
#P9313 = pd.read_csv("../data/Mono-P_mit9313.csv")
#Pmed4 = pd.read_csv("../data/Mono-P_med4.csv")
#P9312 = pd.read_csv("../data/Mono-P_mit9312.csv")
#P9215 = pd.read_csv("../data/Mono-P_mit9215.csv")
