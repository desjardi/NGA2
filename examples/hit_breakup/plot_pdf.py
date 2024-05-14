import numpy as np
import matplotlib as plt
import seaborn as sns
import pandas as pd

diam=[np.genfromtxt('diameter.elvira'),np.genfromtxt('diameter.plicnet')]

sns.set_theme()
# Histogram-style pdf
#sns.displot(diam,bins=40,log_scale=True,stat='density',multiple='dodge',common_norm=False)

# Log-log fitted pdf
sns.displot(diam,kind='kde',log_scale=[True,True],common_norm=False)#,fill=True)

# Log-lin fitted pdf
#sns.displot(diam,kind='kde',log_scale=True,common_norm=False,fill=True)

plt.pyplot.show()
