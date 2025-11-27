

import pandas as pd
import numpy as np
import seaborn as sns
import math as m
import matplotlib.pyplot as plt

# Import dataframe

df = pd.read_csv('Kil_pyMN_dfout.csv')

# Calculete log 10 dpdt 

df['log10 dpdt'] = np.log10(df['dpdt'].values)


# Create plot
plt.figure(figsize=(15, 10))
sns.violinplot(data = df, x = 'Sample', y = 'log10 dpdt', inner = 'quart', gap = 0.5, linewidth = 0.8)

# set plot area

plt.savefig('Ferris_0.5_violin.pdf')
plt.savefig('Ferris_0.5_violin.png')
plt.close()


# Create KDE of all distributions

sns.kdeplot(data= df, x="log10 dpdt", bw_method = 'silverman', fill = True)

plt.savefig('Ferris_0.5_kde.pdf')
plt.savefig('Ferris_0.5_kde.png')
plt.close()
