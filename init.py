# pyrefly: ignore [missing-import]
import signac
# pyrefly: ignore [missing-import]
import numpy as np
import os
from files.python_files import names

project = signac.init_project()

local_eleLam_ljLam_to_initLam = names.eleLam_ljLam_to_initLam
#lambda_ELE = sorted({ele for ele, lj in local_eleLam_ljLam_to_initLam})
#lambda_LJ  = sorted({lj for ele, lj in local_eleLam_ljLam_to_initLam})
# metal = ['La','Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Er','Tm','Lu'] # to test
metal = ['Nd']
polypeptide = ['LBT3-', 'LBT5-']
replicate = [0] # , 1, 2]
lambda_LJ   = [1.0] # [0.0, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 1.0] 
lambda_ELE  = [1.0] # [0.0, 0.075, 0.150, 0.225, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

total_statepoints = list()
legend = open('legend.txt','w')
legend.write('job \t sp \n')
print('job \t sp')

for i in range(len(metal)):
    for j in range(len(replicate)):
        for k in range(len(lambda_ELE)):
            if lambda_ELE[k] == lambda_ELE[-1]:
                for m in range(len(polypeptide)):
                    for l in range(len(lambda_LJ)):
                        statepoint = {
                            "metal": metal[i],
                            "replicate": replicate[j],
                            "lambda_LJ": lambda_LJ[l],
                            "lambda_ELE": lambda_ELE[k],
                            "polypeptide": polypeptide[m]
                        }
                        total_statepoints.append(statepoint)
            else:
                for m in range(len(polypeptide)):
                    statepoint = {
                        "metal": metal[i],
                        "replicate": replicate[j],
                        "lambda_LJ": lambda_LJ[0],
                        "lambda_ELE": lambda_ELE[k],
                        "polypeptide": polypeptide[m]
                    }
                total_statepoints.append(statepoint)

## for i in range(len(metal)):
##     for j in range(len(replicate)):
##         for k in range(len(lambda_LJ)):
##             if lambda_LJ[k] == lambda_LJ[-1]:
##                 for l in range(len(lambda_ELE)):
##                     statepoint = {
##                         "metal": metal[i],
##                         "replicate": replicate[j],
##                         "lambda_LJ": lambda_LJ[k],
##                         "lambda_ELE": lambda_ELE[l]
##                     }
##                     total_statepoints.append(statepoint)
##             else:
##                 statepoint = {
##                     "metal": metal[i],
##                     "replicate": replicate[j],
##                     "lambda_LJ": lambda_LJ[k],
##                     "lambda_ELE": lambda_ELE[0]
##                 }
##                 total_statepoints.append(statepoint)
            
        


for sp in total_statepoints:
    job=project.open_job(
        statepoint=sp,
    ).init()
    legend.write(f' {job} \t\t {sp}\n')
    print(f'{job} \t\t {sp}')
 
 
legend.close()

    
    
