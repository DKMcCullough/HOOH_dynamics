import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy 

df_all = pd.read_csv("../data/BCC_1-31-dataset.csv",header=1)

f,ax = plt.subplots(1,2,figsize=[8,5])

ax[0].plot(np.mean((df_all.rep1+df_all.rep3),axis=0),np.std((df_all.rep1+df_all.rep3),axis=0))
ax[0].plot(np.mean((df_all.rep1+df_all.rep3),axis=0),np.std((df_all.rep1+df_all.rep3),axis=0))

plt.show()


