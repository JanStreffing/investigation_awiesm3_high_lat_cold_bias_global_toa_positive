import glob, os
OUT='/work/ab0246/a270092/software/release_evaluation_tool2/output'
R='/work/bb1469/a270092/eval'
# run -> macro
M={'Tuning_test_06_Baseline':'CMPIBASE','Tuning_test_06A_fesomA_albpnd028':'CMPIA',
   'Tuning_test_06D_HRlike':'CMPID','Tuning_test_06H_fesomH_combo_g_rvice018':'CMPIH',
   'Tuning_test_06O_1hcpl_mospp':'CMPIO','Tuning_test_06T_1hcpl_mospp_kpplow':'CMPIT',
   'Tuning_test_06V_1hcpl_mospp_kpplow_entstpc3_1':'CMPIV'}
score={}
for run in M:
    f=f"{OUT}/{run}/cmpi/frac/{run}_fraction.csv"
    if os.path.exists(f):
        for line in open(f):
            if line.startswith('CMPI'): score[run]=float(line.split()[-1])
lines=[]
for run,mac in M.items():
    if run in score:
        lines.append(f"\\newcommand{{\\{mac}}}{{{score[run]:.3f}}}")
    else:
        lines.append(f"\\newcommand{{\\{mac}}}{{\\textit{{run.}}}}")
# CMPI heatmap grid for runs that have cmpi.png
labels={'Tuning_test_06_Baseline':'Baseline','Tuning_test_06A_fesomA_albpnd028':'06A',
        'Tuning_test_06D_HRlike':'06D HRlike','Tuning_test_06H_fesomH_combo_g_rvice018':'06H',
        'Tuning_test_06O_1hcpl_mospp':'06O','Tuning_test_06T_1hcpl_mospp_kpplow':'06T',
        'Tuning_test_06V_1hcpl_mospp_kpplow_entstpc3_1':'06V'}
have=[r for r in M if os.path.exists(f"{OUT}/{r}/cmpi.png")]
order=['Tuning_test_06_Baseline','Tuning_test_06A_fesomA_albpnd028','Tuning_test_06D_HRlike',
       'Tuning_test_06H_fesomH_combo_g_rvice018','Tuning_test_06O_1hcpl_mospp',
       'Tuning_test_06T_1hcpl_mospp_kpplow','Tuning_test_06V_1hcpl_mospp_kpplow_entstpc3_1']
have=[r for r in order if r in have]
fig=[]
if have:
    fig.append("\\begin{figure}[H]\\centering")
    fig.append("\\begin{tabular}{@{}c@{\\hskip2pt}c@{}}")
    for i in range(0,len(have),2):
        row=have[i:i+2]
        cells=[]
        for r in row:
            cells.append(f"\\includegraphics[width=0.47\\textwidth]{{{r}/cmpi.png}}")
        cap=[f"\\footnotesize {labels[r]} (CMPI {score.get(r,'?'):.3f})" if r in score else f"\\footnotesize {labels[r]}" for r in row]
        if len(row)==1:
            fig.append(cells[0]+" \\\\"); fig.append(cap[0]+" \\\\")
        else:
            fig.append(cells[0]+" & "+cells[1]+" \\\\")
            fig.append(cap[0]+" & "+cap[1]+" \\\\")
    fig.append("\\end{tabular}")
    fig.append("\\caption{\\texttt{part4} CMPI heat maps (per variable/region/season, normalised to CMIP6).}")
    fig.append("\\end{figure}")
else:
    fig.append("\\textit{(CMPI heat maps pending.)}")
lines.append("\\newcommand{\\CMPIFIG}{%\n"+"\n".join(fig)+"\n}")
open(R+'/report_values.tex','w').write("\n".join(lines)+"\n")
print("report_values.tex written; scores:", {labels[k]:round(v,3) for k,v in score.items()})
