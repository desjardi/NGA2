import numpy as np
import matplotlib as plt
import seaborn as sns

diam=[np.genfromtxt('diameter.elvira'),np.genfromtxt('diameter.plicnet')]

sns.set_theme()
#sns.displot(diam,bins=40,log_scale=True,stat="density",multiple="dodge",common_norm=False)
sns.displot(diam,kind='kde',log_scale=True,common_norm=False,fill=True)
plt.pyplot.show()
