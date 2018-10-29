import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

res = [line.strip().split('\t') for line in open("sizes.txt").readlines()]

a = []

for line in res:
    sizes = [int(i) for i in line[1].split(',')]
    strands = [j for j in line[2].split(',')]
    for si, st in zip(sizes, strands):
        a.append((line[0], si, st))

df = pd.DataFrame(a, columns=['name', 'length', 'strand1'])
df.loc[df.strand1 == "True", 'strand'] = '-'
df.loc[df.strand1 == "False", 'strand'] = '+'
df.replace(dict(d6843='Patient1',
                d5945='Patient2'),
           regex=True,
           inplace=True
           )

plt.close("all")
sns.swarmplot(x="name", y="length", data=df, hue="strand")
plt.legend(loc='upper left', title="strand", frameon=False)
plt.show()
