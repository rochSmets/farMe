#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
sys.path.append("/home/smets/codes/far/PHARE/pyphare/")
import pyphare
import matplotlib.pyplot as plt
from pyphare.pharesee.run import Run
import numpy as np
get_ipython().run_line_magic('matplotlib', 'widget')


# In[2]:


ls /home/smets/sherpa/weibel/test/run


# In[4]:


run = Run("/home/smets/sherpa/weibel/test/run")


# In[23]:


time = 1


# In[24]:


B = run.GetB(time)


# In[25]:


type(B)


# In[26]:


for ilvl, lvl in B.levels(time).items():
    print(ilvl, lvl)


# In[27]:


for ilvl, lvl in B.levels(time).items():
    for patch in lvl.patches:
        print(patch.id)


# In[28]:


list_of_ranks = []
n_rank = 12  # could be retrieved from the hierarchy...
ranked = {rank: [] for rank in list(range(n_rank))}

for ilvl, lvl in B.levels(time).items():
    for patch in lvl.patches:
        mpi_rank, patch_id = patch.id.strip('p').split('#')

        ranked[int(mpi_rank)] += [ patch.box.nCells() ]

        if int(mpi_rank) not in list_of_ranks:
            list_of_ranks.append(int(mpi_rank))


# In[29]:


num_of_cells = [sum(ranked[i]) for i in list(range(n_rank))]

print(len(num_of_cells))
print(num_of_cells)

fig, ax = plt.subplots()

ax.plot(num_of_cells)


# In[74]:


pops = ["background", "toleft", "toright"]

for pop in pops:
    parts = run.GetParticles(time, pop)
    print(pop)
    num_of_part = 0
    for ilvl, lvl in parts.levels(time).items():
        for patch in lvl.patches:
            mpi_rank, patch_id = patch.id.strip('p').split('#')
            # print(mpi_rank, patch_id)
            zob = patch.patch_datas
            name = pop+"_domain"
            # print(name)
            # print(dir(zob[name]))
            pipi = zob[name].dataset.iCells
            popo = zob[name].dataset.weights
            # num_of_part += pipi.shape[0]
            # num_of_part += popo.shape
            print(len(popo))
    print(num_of_part)


# In[ ]:





# In[ ]:




