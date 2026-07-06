#!/usr/bin/env python
# coding: utf-8

# In[1]:


import tripyview.sub_notebookheader as nb_setup
# which_matplotlib = 'inline' (default), 
#                    'notebook'(jupyter notebook), 
#                    'widget'(jupyterlab)
nb_setup.init_notebook(which_matplotlib="inline")
# centralized autoimport of: 
# import os
# import warnings
# import time as clock
# import numpy as np
# import xarray as xr
# import shapefile as shp
# import tripyview as tpv
# client, use_existing_client = None, "tcp://0.0.0.0:0000"


# In[2]:


config = {
    'src_path_off_on'    : '/work/ab0246/a270234/runtime/fesom-2.7/test001_core3_momixoff_sppon/outdata/fesom/',
    'src_path_on_on'     : '/work/ab0246/a270234/runtime/fesom-2.7/test001_core3_momixon_sppon/outdata/fesom/',
    'src_path_on_off'    : '/work/ab0246/a270234/runtime/fesom-2.7/test001_core3_momixon_sppoff/outdata/fesom/',
    'src_path_off_off'   : '/work/ab0246/a270234/runtime/fesom-2.7/test002_core3_restart/outdata/fesom/',
    'src_path_new_on_off': '/work/ab0246/a270234/runtime/fesom-2.7/test001_new_core3_momixon_sppoff/outdata/fesom/',
    'mesh_diag_path'     : '/work/ab0246/a270092/input/fesom2/core3/',
    'mesh_path'          : '/work/ab0246/a270092/input/fesom2/core3/',
}


# In[3]:


# parameters
#___Dask Client Parameters____________________________________________________________
do_papermill      = False
do_parallel       = True
parallel_nprc     = 64   # number of total processes
parallel_nthread  = 2    # number of threads per worker --> number of workers then is parallel_nprc/parallel_nthread
parallel_tmem     = 256  # max. available RAM

#___Mesh Path & Save Path_____________________________________________________________
mesh_path         = config['mesh_path']
save_path         = None #'~/figures/test_papermill/'
save_fname        = None # filename from papermill come in through save_fname
tripyrun_name     = None # papermill workflow name of notebook 
tripyrun_analysis = None # papermill diagnostic driver
tripyrun_spath_nb = None # papermill path to processed notebooks
tripyrun_spath_fig= None # papermill path to processed figures

#___Data Path & Input Names___________________________________________________________
input_paths   = list()
input_paths.append(config['src_path_off_on'])
input_paths.append(config['src_path_on_on'])
input_paths.append(config['src_path_on_off'])
input_paths.append(config['src_path_off_off'])
# input_paths.append(config['src_path_new_on_off'])

input_names   = list()
input_names.append('off_on')
input_names.append('on_on')
input_names.append('on_off')
input_names.append('off_off')
# input_names.append('new_on_off')

# n_cycl: which spinupcycle should be plottet if do_allcycl all spinupcycles from [1...n_cycle] are plottet, if None path is directly used
n_cycl         = None
do_allcycl     = False
vname          = 'temp'
year           = [1972, 1981]
mon            = None
day            = None
record         = None 
box            = None
depth          = None
do_edgevec_r2g = True
do_datavec_r2g = False  # set to False if u,v data are already in geo-coordinates
# do_bolus       = True  # add bolus velocity to vec+u+v

#___Define Reference Data, Year, Mon ...______________________________________________
# do anomaly plots in case ref_path is not None
ref_path  = None # '/albedo/work/projects/p_fesom/pscholz/project_TRR181/trr181_tke_ctrl_ck0.1/5/'
ref_name  = None # 'TKE'
ref_year  = None #[1979, 2019]
ref_mon   = None
ref_day   = None
ref_record= None

#___Define Transects via [ [lon], [lat], transect-name ]______________________________
input_transect = list()
input_transect.append([[-80,-70,-60,-50,-40,-35],
                       [-80,-79,-77,-75,-70,-65], 
                       'West Weddell Sea', False])


#___Define Climatology________________________________________________________________
which_clim= 'phc3'
clim_path = '/albedo/work/projects/p_fesom/FROM-OLLIE/FESOM2/hydrography/phc3.0/phc3.0_annual.nc'
# clim_path = '/pool/data/AWICM/FESOM2/INITIAL/phc3.0/phc3.0_annual.nc'

#___Define Colormap Parameters________________________________________________________
# papermill doesnt like multi variable alignment in a single line
cstr      = 'blue2red'
cnum      = 15
cref      = None
crange    = None
cmin      = None
cmax      = None
cfac      = None
climit    = None
chist     = True
ctresh    = 0.995

ref_cstr  = 'wbgyr'
ref_cnum  = 15
ref_cref  = None
ref_crange= None
ref_cmin  = None
ref_cmax  = None
ref_cfac  = None
ref_climit= None
ref_chist = True
ref_ctresh= 0.995

#___Define Plot Parameters____________________________________________________________
ncol              = 2      # number of pannel columns in figure
nrow              = None
box               = [-180, 180, -90, 90]
do_plt            = 'tcf'  # plot pcolor (tpc) or contourf (tcf)
plt_contb         = True   # background contour line (thin)
plt_contf         = False  # contour line of main colorbar steps 
plt_contr         = False  # contour line of reference value 
plt_contl         = False  # label contourline of main colorbar steps 
do_rescale        = None   # rescale data: None, 'log10', 'slog10', np.array(...)
do_lsm            ='fesom' # 'fesom', 'bluemarble', 'etopo', 'stock'
do_mesh, mesh_opt = False, dict({'color':'k', 'linewidth':0.10})
do_enum           = False  # do enumeration of panels
do_reffig         = True   # plot reference fig when doing anomalies 
do_clim           = False   # plot climatolgy values when doing absoluts
ax_title          = None
cb_label          = None
save_dpi          = 300
save_fmt          = ['png']
do_edgevec_r2g    = False


# In[4]:


# start parallel dask client if do_parallel=True
client = tpv.shortcut_setup_daskclient(client, 
                                       use_existing_client, 
                                       do_parallel, 
                                       parallel_nprc, 
                                       parallel_tmem,
                                       threads_per_worker=parallel_nthread)


# In[5]:


#___LOAD FESOM2 MESH___________________________________________________________________________________
mesh=tpv.load_mesh_fesom2(mesh_path, do_rot='None', do_info=True)

#______________________________________________________________________________________________________
# create input_path spinupcycle structure
input_paths, input_names, _ , _ = tpv.shortcut_setup_pathwithspinupcycles(input_paths, input_names, ref_path, ref_name, n_cycl, do_allcycl)
        
#______________________________________________________________________________________________________        
cinfo=tpv.set_cinfo(cstr, cnum, crange, cmin, cmax, cref, cfac, climit, chist, ctresh)
ref_cinfo=None
if (ref_path is not None): 
    if ref_year   is None: ref_year   = year
    if ref_mon    is None: ref_mon    = mon
    if ref_record is None: ref_record = record    
    ref_cinfo=tpv.set_cinfo(ref_cstr, ref_cnum, ref_crange, ref_cmin, ref_cmax, ref_cref, ref_cfac, ref_climit, ref_chist, ref_ctresh)    
    if do_rescale not in ['log10', 'slog10']: cinfo['cref'], ref_cinfo['cref']=0.0, 0.0
    
#______________________________________________________________________________________________________    
# concatenate ref_path and input_path together if is not None,  concatenate list = list1+list2
input_paths, input_names = tpv.shortcut_setup_concatinputrefpath(input_paths, input_names, ref_path, ref_name)

# #______________________________________________________________________________________________________
# # define index regions --> reading shape files
# box = tpv.shortcut_setup_boxregion(box_region)

#______________________________________________________________________________________________________
# set predefined chunks size here! The optimized worker memory dependent chunk size is computed internally. 
# see def compute_optimal_chunks(path, client=None, varname=None, opti_dim='h', opti_chunkfrac=0.10):
# The here presetted values are used when tpv.load_data_fesom2( ..., opti_dim=None', ...), otherwise the 
# chunks are choosen to be not larger than 10% of the worker memory tpv.load_data_fesom2( ..., 
# opti_dim='hori', opti_chunkfrax=0.1, ...). Optimized can be the horizontal, vertical or time 
# dimension opti_dim: 'h', 'v', 'hv', 'vh', 't', 'off', None
chunks = dict({
               'elem' : 'auto', 'nod2' : 'auto', 'edg_n': 'auto',
               'nz1'  : 'auto', 'nz'   : 'auto', 
               'time' : 'auto', 
               }) 


# In[6]:


#______________________________________________________________________________________________________    
# load information about edges 
ts = clock.time()
datapath = input_paths[0]
fname    = 'fesom.mesh.diag.nc'
# check for directory with diagnostic file
if   os.path.isfile( os.path.join(datapath, fname) ): 
    dname = datapath
elif os.path.isfile( os.path.join( os.path.join(os.path.dirname(os.path.normpath(datapath)),'1/'), fname) ): 
    dname = os.path.join(os.path.dirname(os.path.normpath(datapath)),'1/')
elif os.path.isfile( os.path.join(mesh.path,fname) ): 
    dname = mesh.path
else:
    raise ValueError('could not find directory with...mesh.diag.nc file')    

#___________________________________________________________________________
# set specfic type when loading --> #convert to specific precision
# drop unnecessary variables:  Based on the documentation (but a bit hidden), 
# the "data_vars=" argument only works with Python 3.9.
from functools import partial
var_keep = ['edges', 'edge_tri', 'edge_cross_dxdy', 'nod_in_elem2D', 'edge_nodes', 'edge_face_links' ]
def _preprocess(x):
    for var in list(x.keys()):
        if var not in var_keep: 
            x = x.drop_vars(var)
            continue
            
        if x[var].dtype=='float64': x[var] = x[var].astype('float32')
        # there were some changings in the fesom.mesh.diag.nc variable naming from fesom 2.5-->2.6
        # for some reason he load now edge_nodes and edge_face_links as float32, while they cant
        # be used as indices by default and need to be converted by hand into int32
        if   var in ['edge_nodes'     ]: x = x.rename({'edge_nodes'     :'edges'})#.astype('int32')                
        elif var in ['edge_face_links']: x = x.rename({'edge_face_links':'edge_tri'})#.astype('int32')                
    return x
partial_func = partial(_preprocess)

#___________________________________________________________________________
# load diag file --> apply drop variables by preprocessor function
# decode_cf=False is need so that the original fillvalue of -999 is used instead of NaN
# or np.iinfo(np.int32).max()
mdiag = xr.open_mfdataset(os.path.join(dname,fname), parallel=False, 
                             chunks=dict({'edg_n':'auto'}), engine='netcdf4', 
                             preprocess=partial_func,  decode_cf=False)
mdiag = mdiag.drop_vars(list(mdiag.coords)).load()

# node indices of edge points [2 x n2ded]
edge       = mdiag['edges'].values-1
# element indices of triangles that are left and right of edg: [2 x n2ded]
edge_tri   = mdiag['edge_tri'].values-1
# dx & dy of edge midpoints towards element centroid of left and right triangle
edge_dxdy  = mdiag['edge_cross_dxdy'].values[:]

# Be sure that the edge_cross_dxdy variable is in the same rotational frame as your velocities. By default
# edge_cross_dxdy is in rotated coordinates. So if you velocities are also in rotated coordinates things are fine.
# If your velocities should be in geo coordinates than edge_cross_dxdy needs to be rotated as well into geo
# coordinates (do_edgevec_r2g=True)
# if (do_edgevec_r2g):
#     edm_x = mesh.n_x[edge].sum(axis=0)/2.0
#     edm_y = mesh.n_y[edge].sum(axis=0)/2.0
#     edge_dxdy[0,:], edge_dxdy[1,:] = tpv.vec_r2g(mesh.abg, edm_x, edm_y, edge_dxdy[0,:], edge_dxdy[1,:], gridis='geo', do_info=False )
#     edge_dxdy[2,:], edge_dxdy[3,:] = tpv.vec_r2g(mesh.abg, edm_x, edm_y, edge_dxdy[2,:], edge_dxdy[3,:], gridis='geo', do_info=False )
#     del(edm_x, edm_y)
    
# [L]eft  triangle: dx, dy
edge_dxdy_l= np.array([ edge_dxdy[0,:], edge_dxdy[1,:]])
# [R]ight triangle: dx, dy
edge_dxdy_r= np.array([ edge_dxdy[2,:], edge_dxdy[3,:]])
edge_dxdy_r[:, edge_tri[1,:]<0]=0. # if boundarie edge --> right triangle doesnot exist
del(edge_dxdy)

# only needed for plotting transects when scalar data are on elements 
nodeinelem = mdiag['nod_in_elem2D'].values[:,:]-1
if np.sum(nodeinelem==-1)> nodeinelem.shape[0]*nodeinelem.shape[1]-nodeinelem.shape[1]: 
    print(' --> Warning: nod_in_elem2D variable in fesom.mesh.diag.nc might be corrupted, its mostly -1 (old bug, creater actual mesh.diag file), in moment use fall back method')
    nodeinelem = None
del(mdiag)
if client is not None: client.run(gc.collect)
gc.collect()    
print(' --> Load: fesom.mesh.diag.nc --> elasped time: {} min.'.format( (clock.time()-ts)/60  ))

#______________________________________________________________________________________________________    
# analyse transects computes all neccesary arrays
ts = clock.time()
transects = tpv.do_analyse_transects(input_transect, mesh, edge, edge_tri, edge_dxdy_l, edge_dxdy_r, do_rot=do_edgevec_r2g, do_info=False)
for transect in transects:
    print(transect['Name'])
    fig, ax = tpv.plot_transect_position(mesh, transect, edge=edge, do_grid=False)

print(' --> Analyse transects --> elasped time: {} min.'.format( (clock.time()-ts)/60  ))


# In[7]:


t0 = clock.time()
# clean up garabage on workers before the party starts!
if client is not None: client.run(gc.collect)

data_list = list()
#___LOAD FESOM2 DATA___________________________________________________________________________________
for ii, (datapath, descript) in enumerate(zip(input_paths, input_names)):
    print(ii, datapath, descript)
    ts = clock.time()
    #__________________________________________________________________________________________________
    yeari, moni, dayi, recordi = year, mon, day, record
    if (ii==0) and (ref_path != None): yeari, moni, dayi, recordi = ref_year, ref_mon, ref_day, ref_record

    #__________________________________________________________________________________________________
    # input parameter shortcut
    input_dict = dict({'year':yeari, 'mon':moni, 'descript':descript, 'do_rot':False, 
                       'do_info':False, 'do_tarithm':None, 'do_zarithm':None, 'do_ie2n':False, 'do_nan':False, 
                       'chunks':chunks, 'do_load':False, 'do_persist':False, 'do_parallel':do_parallel,
                       'client':client, 'opti_dim':'v' , 'opti_chunkfrac':0.06,
                       'do_info':False,})
    
    # load data
    if vname=='KvN2/N2':
        data  = tpv.load_data_fesom2(mesh, datapath, vname='KvN2' , **input_dict)
        data2 = tpv.load_data_fesom2(mesh, datapath, vname='N2'   , **input_dict)
        data['KvN2'].data = data['KvN2'].data / data2['N2'].data
        data  = data.rename(dict({'KvN2':'KvN2/N2'}))
        data['KvN2/N2'].attrs['units'], data['KvN2/N2'].attrs['description'], data['KvN2/N2'].attrs['long_name'] = '$m^2/s$', '(Kv*N)/N2', '$\\overline{{Kv \\cdot N^2}} / \\overline{{N^2}}$'
        del(data2) 
    else:         
        data  = tpv.load_data_fesom2(mesh, datapath, vname=vname, **input_dict)

    # check if data where loaded
    if data is None: raise ValueError(f'data == None, data could not be readed, your path:{datapath} might be wrong!!!')
    print(' --> elasped time to load {:s} data: {:3.2f} min.'.format(vname, (clock.time()-ts)/60  ))        
    print(' --> data uses {:3.2f} Gb:'.format(data.nbytes/(1024**3)))
    print('')

    #__________________________________________________________________________________________________
    # compute section on reference data if given 
    ts = clock.time()
    if (ii==0) and (ref_path != None):
        csect_ref = tpv.calc_transect_scalar(mesh, data, transects , 
                                             nodeinelem=nodeinelem , 
                                             do_tarithm='mean'     , 
                                             client=client)
        print(' --> elasped time to comp. ref. transect.: {:3.2f} min.'.format( (clock.time()-ts)/60  ))
        for ii, data_ii in enumerate(csect_ref):
            print(' --> csect[{:s}] uses {:3.2f} Mb:'.format(csect_ref[ii][list(csect_ref[ii].keys())[0]].attrs['transect_name'], csect_ref[ii].nbytes/(1024**2)))
        print('')    
        if do_reffig: data_list.append(csect_ref) 
        del(data)    
        continue
    
    # compute section on data
    data_dict = data['temp'][8::12,:,:].to_dataset(name='temp')
    csect = tpv.calc_transect_scalar(mesh, data_dict, transects , 
                                     nodeinelem=nodeinelem , 
                                     do_tarithm='mean'     , 
                                     client=client)
    print(' --> elasped time to comp. data transect.: {:3.2f} min.'.format( (clock.time()-ts)/60  ))   
    for ii, data_ii in enumerate(csect):
        print(' --> csect[{:s}] uses {:3.2f} Mb:'.format(csect[ii][list(csect[ii].keys())[0]].attrs['transect_name'], csect[ii].nbytes/(1024**2)))
    print('')    
    del(data)
    
    #__________________________________________________________________________________________________    
    # compute anomaly/absolute 
    if (ref_path != None): data_list.append(tpv.do_transect_anomaly(csect, csect_ref))
    else                 : data_list.append(csect)
if (ref_path != None): del(csect_ref)        
    
        
#___APPEND ABS CLIMATOLOGY_____________________________________________________________________________    
if (vname in ['temp', 'salt', 'pdens'] or 'sigma' in vname) and do_clim and  (ref_path is None): 
    # load climatology data
    ts = clock.time()
    clim_vname= vname
    if   vname=='temp' and  which_clim.lower()=='woa18': clim_vname = 't00an1'
    elif vname=='salt' and  which_clim.lower()=='woa18': clim_vname = 's00an1'
    clim          = tpv.load_climatology(mesh, clim_path, clim_vname)
    print(' --> elasped time to load clim: {:3.2f} min.'.format( (clock.time()-ts)/60  ))        
    print(' --> clim uses {:3.2f} Gb:'.format(clim.nbytes/(1024**3)))
    print('')

    # compute transect on climatology data
    ts = clock.time()
    clim_transect = tpv.calc_transect_scalar(mesh, clim, transects, nodeinelem=nodeinelem)
    print(' --> elasped time to comp. clim transect.: {:3.2f} sec.'.format( (clock.time()-ts)/60  ))   
    for ii, data_ii in enumerate(clim_transect):
        print(' --> clim_transect[{:s}] uses {:3.2f} Mb:'.format(clim_transect[ii][list(clim_transect[ii].keys())[0]].attrs['transect_name'], clim_transect[ii].nbytes/(1024**2)))
    print('')    
    data_list.append(clim_transect)    
print(' --> total elasped time to process data: {:3.2f} min.'.format( (clock.time()-t0)/60  ))  


# In[8]:


#___PLOT TRANSECT______________________________________________________________________________________
ts = clock.time()
ntrs, ndat = len(transects), len(data_list)
if   ncol != None: 
    ncol0  = np.min([ncol,ndat])    
    nrow0  = np.ceil(ndat/ncol0).astype('int')
elif nrow != None: 
    nrow0  = np.min([nrow,ndat])    
    ncol0  = np.ceil(ndat/nrow0).astype('int')

for trs_idx in range(ntrs):
    svname = list(data_list[0][trs_idx].data_vars)[0]
    slabel = data_list[0][trs_idx][svname].attrs['str_lsave']
    stname = data_list[0][trs_idx][svname].attrs['transect_name'].replace(' ','_').lower()
    #__________________________________________________________________________________________________
    # do save filename path
    spath  = save_path
    sfpath = None
    if spath!=None: 
        sfpath=list()
        for sfmt in save_fmt: sfpath.append( os.path.join(spath,'{:s}_{:s}_{:s}_{:s}.{:s}'.format(svname, 'transect', stname ,slabel, sfmt)) )
    if save_fname!=None: sfpath = [save_fname] # --> needed for diagrun papermille functionality
   
    #__________________________________________________________________________________________________
    # do colorbar either single cbar or ref_cbar + anom_cbar
    if (ref_path != None) and do_reffig: cb_plt, cb_plt_single, cinfo0 = [1]+[2]*(nrow0*ncol0-1), False, [ref_cinfo.copy(), cinfo.copy()]
    else: cb_plt, cb_plt_single, cinfo0 = True, True, cinfo.copy()
    
    #__________________________________________________________________________________________________    
    hfig, hax, hcb = tpv.plot_vslice(mesh, data_list, nrow=nrow0, ncol=ncol0, box_idx=trs_idx, 
                                     cinfo=cinfo0, do_rescale=do_rescale, 
                                     do_plt=do_plt, plt_contb=plt_contb, plt_contf=plt_contf, plt_contr=plt_contr, plt_contl=plt_contl, do_enum=do_enum, 
                                     ax_opt=dict({'fig_sizefac':2.0, 'cb_plt':cb_plt, 'cb_plt_single':cb_plt_single, 'cb_pos':'vertical', 'cb_h':'auto',}), # 'fs_label':14, 'fs_ticks':14, 'ax_dt':1.0}),
                                     cbl_opt=dict(), cb_label=cb_label, cbtl_opt=dict(),
                                     do_save=sfpath, save_dpi=save_dpi )

print(' --> elasped time to plot data: {:3.2f} min.'.format( (clock.time()-ts)/60  ))      


# In[9]:


if do_papermill and do_parallel and client is not None: client.shutdown()

