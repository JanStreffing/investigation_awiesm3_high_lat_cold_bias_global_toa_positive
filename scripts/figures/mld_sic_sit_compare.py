#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import os
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pyfesom2 as pf


# In[48]:


config = {
    'src_path_off_on_cavityon_fesom_core3'             : '/work/ab0246/a270234/runtime/fesom-2.7/test001_core3_momixoff_sppon/outdata/fesom/',
    'src_path_off_off_cavityon_fesom_core3'            : '/work/ab0246/a270234/runtime/fesom-2.7/test002_core3_restart/outdata/fesom/',
    'src_path_off_off_cavityon_esm_core3'   : '/work/bb1469/a270270/runtime/awiesm3-v3.4/AWIESM700_CMIP7_SPINUP_TCO95_CORE3/outdata/fesom/',
    'src_path_off_on_cavityon_esm_hr'     : '/work/bb1469/a270089/runtime/awiesm3-v3.4.1/AWI-ESM3-VEG-HR-CMIP7-Spinup_cont2/outdata/fesom/',
    'mesh_diag_path'     : '/work/ab0246/a270092/input/fesom2/core3/',
    'mesh_path'          : '/work/ab0246/a270092/input/fesom2/core3/',
    'mesh_hr_path'          : '/work/ab0246/a270234/mesh/a270092/dars2/',
    'mesh_hr_diag_path'     : '/work/ab0246/a270234/mesh/a270092/dars2/',
}


# In[9]:


mesh_core3 = pf.load_mesh(config['mesh_diag_path'], usepickle=False, usejoblib=False)
mesh_core3


# In[49]:


mesh_hr = pf.load_mesh(config['mesh_hr_diag_path'], usepickle=False, usejoblib=False)
mesh_hr


# In[12]:


mld_off_on = pf.get_data(config['src_path_off_on_cavityon_fesom_core3'], 'MLD3', range(1972,1982), mesh_core3, how=None, compute=False)
mldw_off_on = mld_off_on[8::12,:].mean(dim='time')
mlds_off_on = mld_off_on[3::12,:].mean(dim='time')


# In[13]:


mld_on_on = pf.get_data(config['src_path_off_off_cavityon_fesom_core3'], 'MLD3', range(1972,1982), mesh_core3, how=None, compute=False)
mldw_on_on = mld_on_on[8::12,:].mean(dim='time')
mlds_on_on = mld_on_on[3::12,:].mean(dim='time')


# In[46]:


mld_on_off = pf.get_data(config['src_path_off_off_cavityon_esm_core3'], 'MLD3', range(1839,1849), mesh_core3, how=None, compute=False)
mldw_on_off = mld_on_off[8::12,:].mean(dim='time')
mlds_on_off = mld_on_off[3::12,:].mean(dim='time')


# In[50]:


mld_off_off = pf.get_data(config['src_path_off_on_cavityon_esm_hr'], 'MLD3', range(1621,1631), mesh_hr, how=None, compute=False)
mldw_off_off = mld_off_off[8::12,:].mean(dim='time')
mlds_off_off = mld_off_off[3::12,:].mean(dim='time')


# In[63]:


pf.plot(mesh_core3, 
        [mldw_off_on.values, mldw_on_on.values, mldw_on_off.values] ,
        mapproj='sp',
        box=[-180, 180, -90, -60],
        rowscol=[1,3],
        figsize=(14,21),
        titles=["off_on_fesom_core3_cavity", "off_off_fesom_core3_cavity", "off_off_esm_core3_cavity"],
        levels=(-2000,-100,20),)
#plt.savefig('MLD2_winter.png')


# In[64]:


pf.plot(mesh_hr, 
        [mldw_off_off.values] ,
        mapproj='sp',
        box=[-180, 180, -90, -60],
        rowscol=[1,3],
        figsize=(14,21),
        titles=["off_on_esm_hr_cavity"],
        levels=(-2000,-100,20),)
#plt.savefig('MLD2_winter.png')


# In[30]:


sic_off_on = pf.get_data(config['src_path_off_on_cavityon_fesom_core3'], 'a_ice', range(1972,1982), mesh_core3, how=None, compute=False)
sicw_off_on = sic_off_on[8::12,:].mean(dim='time')
sics_off_on = sic_off_on[3::12,:].mean(dim='time')


# In[31]:


sic_on_on = pf.get_data(config['src_path_off_off_cavityon_fesom_core3'], 'a_ice', range(1972,1982), mesh_core3, how=None, compute=False)
sicw_on_on = sic_on_on[8::12,:].mean(dim='time')
sics_on_on = sic_on_on[3::12,:].mean(dim='time')


# In[32]:


sic_on_off = pf.get_data(config['src_path_off_off_cavityoff_esm_core3'], 'a_ice', range(1840,1850), mesh_core3, how=None, compute=False)
sicw_on_off = sic_on_off[8::12,:].mean(dim='time')
sics_on_off = sic_on_off[3::12,:].mean(dim='time')


# In[56]:


sic_off_off = pf.get_data(config['src_path_off_on_cavityon_esm_hr'], 'a_ice', range(1622,1632), mesh_hr, how=None, compute=False)
sicw_off_off = sic_off_off[8::12,:].mean(dim='time')
sics_off_off = sic_off_off[3::12,:].mean(dim='time')


# In[36]:


pf.plot(mesh_core3, 
        [sicw_off_on.values, sicw_on_on.values, sicw_on_off.values] ,
        mapproj='sp',
        box=[-180, 180, -90, -50],
        rowscol=[1,3],
        figsize=(14,21),
        titles=["off_on_fesom_core3_cavity", "off_off_fesom_core3_cavity", "off_off_esm_core3"],
        levels=(0.05,1,20),)
#plt.savefig('SIC2_winter.png')


# In[58]:


pf.plot(mesh_hr, 
        [sicw_off_off.values] ,
        mapproj='sp',
        box=[-180, 180, -90, -60],
        rowscol=[1,3],
        figsize=(14,21),
        titles=["off_on_esm_hr_cavity"],
        levels=(0.05,1,20),)
#plt.savefig('MLD2_winter.png')


# In[39]:


sit_off_on = pf.get_data(config['src_path_off_on_cavityon_fesom_core3'], 'm_ice', range(1972,1982), mesh_core3, how=None, compute=False)
sitw_off_on = sit_off_on[8::12,:].mean(dim='time')
sits_off_on = sit_off_on[3::12,:].mean(dim='time')


# In[40]:


sit_on_on = pf.get_data(config['src_path_off_off_cavityon_fesom_core3'], 'm_ice', range(1972,1982), mesh_core3, how=None, compute=False)
sitw_on_on = sit_on_on[8::12,:].mean(dim='time')
sits_on_on = sit_on_on[3::12,:].mean(dim='time')


# In[41]:


sit_on_off = pf.get_data(config['src_path_off_off_cavityoff_esm_core3'], 'm_ice', range(1840,1850), mesh_core3, how=None, compute=False)
sitw_on_off = sit_on_off[8::12,:].mean(dim='time')
sits_on_off = sit_on_off[3::12,:].mean(dim='time')


# In[59]:


sit_off_off = pf.get_data(config['src_path_off_on_cavityon_esm_hr'], 'm_ice', range(1622,1632), mesh_hr, how=None, compute=False)
sitw_off_off = sit_off_off[8::12,:].mean(dim='time')
sits_off_off = sit_off_off[3::12,:].mean(dim='time')


# In[43]:


pf.plot(mesh_core3, 
        [sitw_off_on.values, sitw_on_on.values, sitw_on_off.values] ,
        mapproj='sp',
        box=[-180, 180, -90, -50],
        rowscol=[1,3],
        figsize=(14,21),
        titles=["off_on_fesom_core3_cavity", "off_off_fesom_core3_cavity", "off_off_esm_core3"],
        levels=(0.1,2,20),)
#plt.savefig('SIC2_summer.png')


# In[60]:


pf.plot(mesh_hr, 
        [ sitw_off_off.values] ,
        mapproj='sp',
        box=[-180, 180, -90, -50],
        rowscol=[1,3],
        figsize=(14,21),
        titles=["off_on_esm_hr_cavity"],
        levels=(0.1,2,20),)

