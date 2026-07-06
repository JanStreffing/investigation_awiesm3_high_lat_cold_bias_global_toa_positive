#!/usr/bin/env python
# Adapted from a270234's temp_fesom_core3_WWS.ipynb (tripyview) to the coupled
# Tuning_test_06 runs: austral-winter (Sep) potential-temperature vertical section
# along the West Weddell Sea transect, Baseline vs the ocean-mixing/coupling runs.
import os, gc, time as clock, numpy as np, xarray as xr
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import tripyview as tpv
from functools import partial

BASE='/work/bb1469/a270270/runtime/awiesm3-v3.4'
MESHP='/work/ab0246/a270092/input/fesom2/core3/'
OUTP='/work/ab0246/a270092/software/release_evaluation_tool2/output/Tuning_test_06_overview/plots'
os.makedirs(OUTP, exist_ok=True)
client=None
RUNS=[('Tuning_test_06_Baseline','Baseline'),
      ('Tuning_test_06O_1hcpl_mospp','06O 1hcpl+mospp'),
      ('Tuning_test_06T_1hcpl_mospp_kpplow','06T +kpplow'),
      ('Tuning_test_06V_1hcpl_mospp_kpplow_entstpc3_1','06V aggressive')]
input_paths=[f'{BASE}/{r}/outdata/fesom/' for r,_ in RUNS]
input_names=[lab for _,lab in RUNS]
year=[1370,1379]; vname='temp'
input_transect=[[[-80,-70,-60,-50,-40,-35],[-80,-79,-77,-75,-70,-65],'West Weddell Sea',False]]

print("load mesh ...")
mesh=tpv.load_mesh_fesom2(MESHP, do_rot='None', do_info=False)
cinfo=tpv.set_cinfo('blue2red',15,None,None,None,None,None,None,True,0.995)

# --- load edge info from fesom.mesh.diag.nc (needed for transects) ---
dname=MESHP; fname='fesom.mesh.diag.nc'
var_keep=['edges','edge_tri','edge_cross_dxdy','nod_in_elem2D','edge_nodes','edge_face_links']
def _pre(x):
    for var in list(x.keys()):
        if var not in var_keep: x=x.drop_vars(var); continue
        if x[var].dtype=='float64': x[var]=x[var].astype('float32')
        if var=='edge_nodes': x=x.rename({'edge_nodes':'edges'})
        elif var=='edge_face_links': x=x.rename({'edge_face_links':'edge_tri'})
    return x
mdiag=xr.open_mfdataset(os.path.join(dname,fname),parallel=False,chunks=dict({'edg_n':'auto'}),
                        engine='netcdf4',preprocess=partial(_pre),decode_cf=False)
mdiag=mdiag.drop_vars(list(mdiag.coords)).load()
edge=mdiag['edges'].values-1
edge_tri=mdiag['edge_tri'].values-1
edge_dxdy=mdiag['edge_cross_dxdy'].values[:]
edge_dxdy_l=np.array([edge_dxdy[0,:],edge_dxdy[1,:]])
edge_dxdy_r=np.array([edge_dxdy[2,:],edge_dxdy[3,:]])
edge_dxdy_r[:,edge_tri[1,:]<0]=0.
nodeinelem=mdiag['nod_in_elem2D'].values[:,:]-1
if np.sum(nodeinelem==-1)>nodeinelem.shape[0]*nodeinelem.shape[1]-nodeinelem.shape[1]: nodeinelem=None
del mdiag; gc.collect()

transects=tpv.do_analyse_transects(input_transect,mesh,edge,edge_tri,edge_dxdy_l,edge_dxdy_r,do_rot=False,do_info=False)

data_list=[]
for datapath,descript in zip(input_paths,input_names):
    print("load",descript)
    idict=dict(year=year,mon=None,descript=descript,do_rot=False,do_info=False,do_tarithm=None,
               do_zarithm=None,do_ie2n=False,do_nan=False,do_load=False,do_persist=False,
               do_parallel=False,client=client,opti_dim='v',opti_chunkfrac=0.06)
    data=tpv.load_data_fesom2(mesh,datapath,vname=vname,**idict)
    data_dict=data['temp'][8::12,:,:].to_dataset(name='temp')   # September
    csect=tpv.calc_transect_scalar(mesh,data_dict,transects,nodeinelem=nodeinelem,do_tarithm='mean',client=client)
    data_list.append(csect); del data; gc.collect()

ntrs=len(transects); ndat=len(data_list); ncol=min(2,ndat); nrow=int(np.ceil(ndat/ncol))
for ti in range(ntrs):
    svname=list(data_list[0][ti].data_vars)[0]
    stname=data_list[0][ti][svname].attrs['transect_name'].replace(' ','_').lower()
    sfpath=[os.path.join(OUTP,f'temp_transect_{stname}_tuning.png')]
    tpv.plot_vslice(mesh,data_list,nrow=nrow,ncol=ncol,box_idx=ti,cinfo=cinfo.copy(),
                    do_plt='tcf',plt_contb=True,
                    ax_opt=dict(fig_sizefac=2.0,cb_plt=True,cb_plt_single=True,cb_pos='vertical',cb_h='auto'),
                    do_save=sfpath,save_dpi=160)
    print("saved",sfpath[0])
print("DONE")
