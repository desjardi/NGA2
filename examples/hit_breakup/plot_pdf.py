import numpy as np
import matplotlib as plt
import seaborn as sns
import pandas as pd

# Import our diameter data
elvira=pd.DataFrame({'Normalized diameter':np.genfromtxt('elvira/diameter.192.elvira/diameter_1.00000E+01')})
elvira['Method']='ELVIRA'
elvira['Weight']=0.2
lvira=pd.DataFrame({'Normalized diameter':np.genfromtxt('lvira/diameter.192.lvira/diameter_1.00000E+01')})
lvira['Method']='LVIRA'
lvira['Weight']=0.115
plicnet=pd.DataFrame({'Normalized diameter':np.genfromtxt('plicnet/diameter.192.plicnet/diameter_1.00000E+01')})
plicnet['Method']='PLICnet'
plicnet['Weight']=0.1

data=pd.concat([elvira,lvira,plicnet],ignore_index=True)

# Set plotting style
sns.set_theme()

# Histogram-style pdf
#g=sns.displot(data,x='Normalized diameter',hue='Method',weights='Weight',bins=32,log_scale=True,stat='density',multiple='layer',common_norm=False)#,kde=True)#element='poly'
g=sns.displot(data,bins=26,x='Normalized diameter',hue='Method',weights='Weight',log_scale=True,multiple='dodge')#,kde=True)
g.set(yscale='log')
g.set(ylabel='Scaled density')
g.set(xbound=[1e-4,1])
g.set(ybound=[1e-1,1e3])

# Log-log fitted pdf
#sns.displot(diam,kind='kde',log_scale=[True,True],common_norm=False)#,fill=True)

# Log-lin fitted pdf
#sns.displot(diam,kind='kde',log_scale=True,common_norm=False,fill=True)

plt.pyplot.show()
