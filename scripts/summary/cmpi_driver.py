import sys, os, json
REPO='/work/ab0246/a270092/software/release_evaluation_tool2/bg_routines/cmpitool'
sys.path.insert(0, REPO)
from cmpitool import cmpitool, cmpisetup
R='/work/ab0246/a270092/software/release_evaluation_tool2/output/Tuning_test_06_overview'
variable, region, climate_model, siconc, tas, clt, pr, rlut, uas, vas, ua, zg, zos, mlotst, thetao, so = cmpisetup()
mv=[siconc, tas, clt, pr, rlut, uas, vas, ua, zg, zos, mlotst, thetao, so]
names=['Baseline','06A_albpnd028','06D_HRlike','06H_combo','06O_1hcpl_mospp','06T_kpplow','06V_entstpc3']
models=[climate_model(name=n, variables=mv) for n in names]
fixed_limits={'siconc':60.0,'tas':5.0,'clt':30.0,'pr':5.0,'rlut':20.0,'uas':3.0,'vas':3.0,'ua':5.0,'zg':100.0,'zos':0.3,'mlotst':100.0,'thetao':3.0,'so':1.0}
res=cmpitool(R+'/cmpi_input_shared', models, verbose=True, biasmaps=True, biasmap_limits=fixed_limits,
             obs_path=os.path.join(REPO,'obs'), eval_path=os.path.join(REPO,'eval','ERA5'),
             out_path=R+'/cmpi', use_for_eval=False)
try:
    with open(R+'/data/cmpi_result.json','w') as f: json.dump({k:(v if isinstance(v,(int,float,str)) else str(v)) for k,v in (res.items() if hasattr(res,'items') else {})}, f, indent=2)
except Exception as e:
    print("json dump issue:", e)
print("CMPI DONE", res)
