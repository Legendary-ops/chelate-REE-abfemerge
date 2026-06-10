# my goal is to run the worklow manually calling the functions directly.
#from alchemtest.gmx import load_ABFE
#from alchemlyb.workflows import ABFE

# my imports
import pdb
from alchemlyb.parsing.gmx import extract_u_nk
from alchemlyb.parsing.gmx import extract_dHdl
from alchemlyb.preprocessing.subsampling import decorrelate_u_nk
from alchemlyb.preprocessing.subsampling import decorrelate_dhdl
import pandas as pd
import os
from pathlib import Path
from alchemlyb.estimators import MBAR, TI
from alchemlyb.postprocessors.units import get_unit_converter
#from alchemlyb.visualisation.mbar_matrix import plot_mbar_overlap_matrix
#from alchemlyb.visualisation.dF_state import plot_dF_state
from alchemlyb.convergence import forward_backward_convergence
from alchemlyb.visualisation import plot_convergence, plot_dF_state, plot_mbar_overlap_matrix, plot_ti_dhdl

#dir = os.path.dirname(load_ABFE()['data']['complex'][0]) # just a dir name
dir = '.'

UNIT_INPUT = 'kJ/mol'
SIMULATION_TEMPERATURE = 298.0

file_pattern = '**/PRO_CANON_*.xvg'

GENERAL_FILE_PREFIX = 'alchemlyb_HFE'

TI_SUFFIX = 'TI'; MBAR_SUFFIX = 'MBAR'

CONVERGNCE_ROOT = 'convergence'; OVERLAP_ROOT = 'overlap'
STATE_ROOT = 'state'; DH_DL_ROOT = 'dhdl'

list_of_xvg_files = list(map(str, Path(dir).glob(file_pattern)))

u_nk_list = []; decorrelated_u_nk_list = [] 
dHdl_list = []; decorrelated_dhdl_list = [] 

for file in list_of_xvg_files:
    u_nk_list.append(extract_u_nk(file, T=SIMULATION_TEMPERATURE))
    dHdl_list.append(extract_dHdl(file, T=SIMULATION_TEMPERATURE))
    
column_names = u_nk_list[0].columns.values.tolist()

list_to_ensure_files_allign_with_columnLambas = sorted(
    range(len(list_of_xvg_files)),
    key = lambda x: column_names.index(u_nk_list[x].reset_index("time").index.values[0])
)

# apply the sorting so that the files appear in the order of increasing lambda state.
u_nk_list = [u_nk_list[i] for i in list_to_ensure_files_allign_with_columnLambas]
dHdl_list = [dHdl_list[i] for i in list_to_ensure_files_allign_with_columnLambas]

for index_of_unk, u_nk in enumerate(u_nk_list):
    decorrelated_u_nk_list.append(decorrelate_u_nk(u_nk))
    
for index_of_unk, dHdl in enumerate(dHdl_list):
    decorrelated_dhdl_list.append(decorrelate_dhdl(dHdl))
    
cat_decorrelated_u_nk_list = pd.concat(decorrelated_u_nk_list)
cat_decorrelated_dhdl_list = pd.concat(decorrelated_dhdl_list)

mbar_result = MBAR().fit(cat_decorrelated_u_nk_list)
dHdl_result = TI().fit(cat_decorrelated_dhdl_list)

## PERFORM UNFATHOMABLE SORTING OF EACH SIMULATION AND ITS LAMBDA NEIGHBORS
## lines 504-611 of abfe.py; from alchemlyb/workflows/abfe.py
## black box of not understanding the code starts here.

num_mbar_states = len(mbar_result.states_)
num_dHdl_states = len(dHdl_result.states_)

# Verify coupling state from any .mdp file in the directory
import glob
mdp_files = glob.glob('*.mdp')
couple_lambda0 = None
couple_lambda1 = None

if mdp_files:
    mdp_file = mdp_files[0]
    print(f"Reading MDP file: {mdp_file}")
    with open(mdp_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(';') or line.startswith('#'):
                continue
            if 'couple-lambda0' in line:
                parts = line.split('=')
                if len(parts) > 1:
                    couple_lambda0 = parts[1].split(';')[0].strip().lower()
            elif 'couple-lambda1' in line:
                parts = line.split('=')
                if len(parts) > 1:
                    couple_lambda1 = parts[1].split(';')[0].strip().lower()
    print(f"Found couple-lambda0: {couple_lambda0}")
    print(f"Found couple-lambda1: {couple_lambda1}")
else:
    print("Warning: No .mdp files found in the current directory. Cannot automatically verify the coupling state order.")
    print("Defaulting to sign flip (assuming state 0 is coupled and state 18 is decoupled).")

# Determine sign flip
sign_flip = True  # Default assuming state 18 is decoupled
if couple_lambda0 is not None or couple_lambda1 is not None:
    if couple_lambda0 == 'none':
        sign_flip = False
    elif couple_lambda1 == 'none':
        sign_flip = True

# Conversion factors: kT to kJ/mol
conversion_factor = SIMULATION_TEMPERATURE * 8.314 / 1000.0
sign_mult = -1.0 if sign_flip else 1.0
total_factor = conversion_factor * sign_mult

print(f"Conversion factor: {conversion_factor:.6f} kJ/mol per kT")
print(f"Sign flip applied: {sign_flip} (multiplier: {sign_mult})")

data_dict = {"name": [], "state": []}

for i in range(num_mbar_states-1):
    data_dict["name"].append(str(i) + " -- " + str(i + 1))
    data_dict["state"].append("States")

stages =  u_nk_list[0].reset_index("time").index.names

for stage in stages:
    data_dict["name"].append(stage.split("-")[0])
    data_dict["state"].append("Stages")
    
data_dict["name"].append("TOTAL")
data_dict["state"].append("Stages")

col_names = ['MBAR','MBAR_Error']
mbar_delta_f_ = mbar_result.delta_f_
mbar_d_delta_f_ = mbar_result.d_delta_f_

data_dict['MBAR'] = []
data_dict['MBAR_Error'] = []

for index in range(1, num_mbar_states):
    data_dict['MBAR'].append(mbar_delta_f_.iloc[index - 1, index] * total_factor)
    data_dict['MBAR_Error'].append(mbar_d_delta_f_.iloc[index - 1, index] * conversion_factor)
    
for index, stage in enumerate(stages):
    if len(stages) == 1:
        start = 0
        end = num_mbar_states - 1
    else:
        lambda_min = min([state[index] for state in mbar_result.states_])
        lambda_max = max([state[index] for state in mbar_result.states_])
        if lambda_min == lambda_max:
            start = 0
            end = 0
        else:
            states = [state[index] for state in mbar_result.states_]    
            start = list(reversed(states)).index(lambda_min)
            start = num_mbar_states - start - 1
            end = states.index(lambda_max)
            
    result = mbar_delta_f_.iloc[start, end] * total_factor
    error = mbar_d_delta_f_.iloc[start, end] * conversion_factor
    data_dict['MBAR'].append(result)
    data_dict['MBAR_Error'].append(error)

# Appending TOTAL (the last row)
total_result = mbar_delta_f_.iloc[0, num_mbar_states - 1] * total_factor
total_error = mbar_d_delta_f_.iloc[0, num_mbar_states - 1] * conversion_factor
data_dict['MBAR'].append(total_result)
data_dict['MBAR_Error'].append(total_error)

summary = pd.DataFrame.from_dict(data_dict)
summary = summary.set_index(["state", "name"])
summary.index.names = [None, None]
summary = summary.to_string()

summary_txt = open(f'{GENERAL_FILE_PREFIX}.txt','w')
summary_txt.write(f'Free energy results in kJ/mol (converted from kT at T = {SIMULATION_TEMPERATURE} K).\n')
if not mdp_files:
    summary_txt.write(f'WARNING: No .mdp files found in the current directory. Assuming state 0 is coupled and state 18 is decoupled (sign flip applied).\n')
summary_txt.write(f'Sign flip applied: {sign_flip} (Target: Solvation/Hydration free energy).\n\n')
summary_txt.write(f'{summary}\n\n\n______________________________________________________________\n\n\n')

######______________________________________________________________######

data_dict_ti = {"name": [], "state": []}

for i in range(num_dHdl_states-1):
    data_dict_ti["name"].append(str(i) + " -- " + str(i + 1))
    data_dict_ti["state"].append("States")

for stage in stages:
    data_dict_ti["name"].append(stage.split("-")[0])
    data_dict_ti["state"].append("Stages")
    
data_dict_ti["name"].append("TOTAL")
data_dict_ti["state"].append("Stages")

col_names = ['TI','TI_Error']
dHdl_delta_f_ = dHdl_result.delta_f_
dHdl_d_delta_f_ = dHdl_result.d_delta_f_

data_dict_ti['TI'] = []
data_dict_ti['TI_Error'] = []

for index in range(1, num_dHdl_states):
    data_dict_ti['TI'].append(dHdl_delta_f_.iloc[index - 1, index] * total_factor)
    data_dict_ti['TI_Error'].append(dHdl_d_delta_f_.iloc[index - 1, index] * conversion_factor)
    
for index, stage in enumerate(stages):
    if len(stages) == 1:
        start = 0
        end = num_dHdl_states - 1
    else:
        lambda_min = min([state[index] for state in dHdl_result.states_])
        lambda_max = max([state[index] for state in dHdl_result.states_])
        if lambda_min == lambda_max:
            start = 0
            end = 0
        else:
            states = [state[index] for state in dHdl_result.states_]    
            start = list(reversed(states)).index(lambda_min)
            start = num_dHdl_states - start - 1
            end = states.index(lambda_max)
            
    result = dHdl_delta_f_.iloc[start, end] * total_factor
    error = dHdl_d_delta_f_.iloc[start, end] * conversion_factor
    data_dict_ti['TI'].append(result)
    data_dict_ti['TI_Error'].append(error)

# Appending TOTAL (the last row)
total_result_ti = dHdl_delta_f_.iloc[0, num_dHdl_states - 1] * total_factor
total_error_ti = dHdl_d_delta_f_.iloc[0, num_dHdl_states - 1] * conversion_factor
data_dict_ti['TI'].append(total_result_ti)
data_dict_ti['TI_Error'].append(total_error_ti)

summary_ti = pd.DataFrame.from_dict(data_dict_ti)
summary_ti = summary_ti.set_index(["state", "name"])
summary_ti.index.names = [None, None]
summary_ti = summary_ti.to_string()

summary_txt.write(f'\n\n\n')
summary_txt.write(f'______________________________________________________________')
summary_txt.write(f'\n\n\n')
summary_txt.write(f'{summary_ti}')
summary_txt.close()

#summary_txt = open(f'{GENERAL_FILE_PREFIX}.txt','w')

## print(f"len(data_dict['name']): {len(data_dict['name'])}")
## print(f"len(data_dict['state']): {len(data_dict['state'])}")
## 
## print(f"len(data_dict['TI']): {len(data_dict['TI'])}")
## print(f"len(data_dict['TI_Error']): {len(data_dict['TI_Error'])}")

## black box of not understanding the code ends here.

ax = plot_mbar_overlap_matrix(mbar_result.overlap_matrix)
ax.figure.savefig(f'{GENERAL_FILE_PREFIX}_{OVERLAP_ROOT}_{MBAR_SUFFIX}.png')
fig = plot_dF_state([mbar_result,dHdl_result],units=UNIT_INPUT)
fig.savefig(f'{GENERAL_FILE_PREFIX}_{STATE_ROOT}_{TI_SUFFIX}with{MBAR_SUFFIX}.png')

del ax, fig

convergence = forward_backward_convergence(u_nk_list, 'MBAR', num=min(10, len(u_nk_list[0])))

unit_converted_convergence = get_unit_converter(UNIT_INPUT)(convergence)
unit_converted_convergence["data_fraction"] = convergence["data_fraction"]

ax = plot_convergence(convergence, units=UNIT_INPUT)
ax.figure.savefig(f'{GENERAL_FILE_PREFIX}_{CONVERGNCE_ROOT}_{MBAR_SUFFIX}.png')

del ax

ax = plot_ti_dhdl(dHdl_result, units=UNIT_INPUT)
ax.figure.savefig(f'{GENERAL_FILE_PREFIX}_{DH_DL_ROOT}_{TI_SUFFIX}.png')


## states_mbar = mbar_result.states_; mbar_HFE = mbar_result.delta_f_; mbar_fluctuation = mbar_result.d_delta_f_
## states_dHdl = dHdl_result.states_; dHdl_HFE = dHdl_result.delta_f_; dHdl_fluctuation = dHdl_result.d_delta_f_
## 
## 
## print(f'using MBAR')
## for index, value in enumerate(states_mbar): #range(0, num_states_mbar):
##     print(f'init_lambda_state: {index} states_mbar: {states_mbar.iloc[]} HFE: {mbar_HFE} fluctuation: {mbar_fluctuation}')
##     
##     
## print(f'using dHdl (TI)')
## for index, value in enumerate(states_dHdl): #range(0, num_states_dHdl):
##     print(f'init_lambda_state: {index} HFE: {dHdl_HFE} fluctuation: {dHdl_fluctuation}')
