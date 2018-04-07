#%%
import numpy as np
import re

#%%
B = list(); b = list()
with open('Btest', mode='r') as f:
    lines = f.readlines()

def getLines(i):
    return list(map(lambda s: float(s), re.sub(r'\s\s+', ' ', lines[i]).strip().split()))

for i in range(1, len(lines)):
    B.append((getLines(i)[:-1]))
    b.append(getLines(i)[-1])

B = np.matrix(B)
b = np.matrix(b)

for i in range(B.shape[0]):
    for j in range(B.shape[1]):
        B[i,j] = -B[i,j] / B[i,i] if i!=j else B[i,j]
    b[0,i] = b[0,i] / B[i,i]
    B[i,i] = 0.
#%%
print(B)
#%%
print(b)
#%%
x = np.matrix(np.zeros(b.shape))
#%%
x = (B * x.T + b.T).T
print(x)
