## process sim output

from pathlib import Path
import re, pickle
import numpy as np
import pandas as pd
# from scipy.stats import chi2, kstest
import matplotlib.pyplot as plt


# save_path = './sim/260712/'
save_path = './sim/260713/'
files = [p.name for p in Path(save_path).iterdir() if p.is_file()]
files = [save_path + file for file in files if re.search('^save.*pkl$',file)]

# file = files[1]
ns = []; ks = []; ks_bartlett = []
for file in files:
    with open(file, 'rb') as f:
        save_data = pickle.load(f)
    ns += save_data['ns']
    ks += save_data['ks']
    ks_bartlett += save_data['ks_bartlett']
# ns, ks, ks_bartlett += [save_data[k] for k in ['ns','ks', 'ks_bartlett']]
# rate = np.polyfit(np.log(list(ns)),np.log(ks),deg=1)[0]
# rate_bartlett = np.polyfit(np.log(list(ns)),np.log(ks_bartlett),deg=1)[0]
long = pd.DataFrame({'n':ns,'ks':ks,'ks_bartlett':ks_bartlett})
distances = long.groupby('n',as_index=False)[['ks','ks_bartlett']].median()
rate = np.polyfit(list(np.log(distances['n'])), list(np.log(distances['ks'])), deg=1)[0]
rate_bartlett = np.polyfit(list(np.log(distances['n'])), list(np.log(distances['ks_bartlett'])), deg=1)[0]
print(rate,rate_bartlett)
distances.plot(x='n',loglog=True)
plt.show()






