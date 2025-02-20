#!/usr/intel/bin/python3
import UsrIntel.R1
import subprocess
import argparse
import os
import pandas as pd
from tabulate import tabulate
parser = argparse.ArgumentParser()
parser.add_argument('-p','--path',default = os.getcwd())
parser.add_argument('-c','--coho_rerun',action = 'store_const',const = True)
parser.add_argument('-dcm',action= 'store_const',const = True)
args = parser.parse_args()

def run_cmd(cmd):
    try:
        return subprocess.check_output(cmd,universal_newlines = True,shell = True).splitlines()
    except subprocess.CalledProcessError as err:
        print(err.output.rstrip())
        return None	
          

def get_original_run_directory(path):
    if 'coho_rerun' in path and os.path.exists(path):
        for directory in os.listdir(path):
            if directory.startswith('bcs_') and directory != 'bcs_policy_sc':
                return directory
    print('Original run directory does not exists!!')		
    return None	

def get_trace_details(core_te_dbg_file,coho_rt_file,cores):
    cmd = f"/nfs/site/disks/BDC_core_validation.reg.01/USERS/smadugun/infra_works/ajith_bare_bones_rtl_vs_coho.py -l {core_te_dbg_file} -c {coho_rt_file} -n {cores}"
    trace_details = dict()
    result = run_cmd(cmd)
    #print(result)
    if cores == 2:    
        max_rtl_cycles = max_coho_cycles = 0    
        for line in result:
           
            line_split = line.split(',')
            if 'c0' in line:
                trace_details[f'{line_split[0]}_{line_split[1]}'] = {'coho':line_split[2],'rtl':line_split[3]}
            elif 'c1' in line:
                trace_details[f'{line_split[0]}_{line_split[1]}'] = {'coho':line_split[2],'rtl':line_split[3]}
    #print(trace_details)
    return trace_details
def find_stat(pattern,file):
    try:
        value = round(float(subprocess.check_output(["zgrep",f'{pattern}',file], stderr=subprocess.DEVNULL, universal_newlines=True).split()[1]),2)
    except Exception as e:
        value = 0
    return value

def get_apcl(trace):

    if '4x16B' in trace:
        LPC = 4
    elif '1x16B' in trace:
        LPC = 1
    elif '2x32B' in trace:
        LPC = 2
    else:
        LPC = 1
    return LPC	

def display(results_list,display_option):
    if 'dcm' in display_option:
        print(tabulate(results_list,
                tablefmt = 'grid',
                headers = ['Trace','COHO_RDS','COHO_WTS','RTL_RDS','RTL_WTS']))
    else:
        print(tabulate(results,
                tablefmt='grid',
                headers=['Trace','C0 COHO_RDS','C1 COHO_RDS','C0 RTL_RDS','C1 RTL_RDS','C0 COHO_WTS','C1 COHO_WTS','C0 RTL_WTS','C1 RTL_WTS']))    
    
if 'coho_rerun' in args.path:
    original_run_dir = os.path.join(os.path.join(args.path,get_original_run_directory(args.path)),'runs','exp')
    runs_path_for_coho_rerun = os.path.join(args.path,'runs')    

runs_path = os.path.join(args.path,'runs','exp') if 'coho_rerun' not in args.path else original_run_dir   
cores = 2 if 'DCM' in args.path else 1
results = []
for trace in os.listdir(runs_path):
    trace_path = os.path.join(runs_path,trace)
    if os.path.isdir(trace_path) and not trace.startswith('##'):

        indigo_stats_file = os.path.join(trace_path,'indigo_files','indigo.stats.stats.gz')
        core_te_dbg_file = os.path.join(trace_path,'core_te_dbg.elog.gz')
        coho_skl_file = os.path.join(trace_path,trace+'.coho_skl.gz') if 'coho_rerun' not in args.path else os.path.join(runs_path_for_coho_rerun,trace,trace+'.coho_skl.gz')	
        coho_rt_file = os.path.join(trace_path,trace+'.coho_skl.rt.gz') if 'coho_rerun' not in args.path else os.path.join(runs_path_for_coho_rerun,trace,trace+'.coho_skl.rt.gz')	
        rtl_reads_data = data_read = data_read_pref = code_read = code_read_pref = 0
        rtl_writes = rtl_stores_16 = rtl_stores_32 = rtl_stores_64 = rtl_other_writes = 0
        coho_reads_data = 0 
        coho_writes = coho_stores_16 = coho_stores_32 = coho_stores_64 = coho_other_writes = 0	
        if os.path.isfile(indigo_stats_file) and os.path.isfile(core_te_dbg_file) and os.path.isfile(coho_skl_file) and os.path.isfile(coho_rt_file):
            trace_details = get_trace_details(core_te_dbg_file,coho_rt_file,cores)
            core0_coho_cycles = int(trace_details[f'c0_{trace}']['coho'])
            core0_rtl_cycles = int(trace_details[f'c0_{trace}']['rtl'])
            core1_coho_cycles = int(trace_details[f'c1_{trace}']['coho'])
            core1_rtl_cycles = int(trace_details[f'c1_{trace}']['rtl'])
            print(trace,core0_coho_cycles,core1_coho_cycles,core0_rtl_cycles,core1_rtl_cycles)
            
            data_read_rtl = find_stat('S*.p0.iil1.idi_0.data_read',indigo_stats_file)
            data_read_pref_rtl = find_stat('S*.p0.iil1.idi_0.data_read_pref',indigo_stats_file)
            data_read_pte_rtl = find_stat('S*.p0.iil1.idi_0.data_read_pte',indigo_stats_file)
            code_read_rtl = find_stat('S*.p0.iil1.idi_0.code_read',indigo_stats_file)
            code_read_pref_rtl = find_stat('S*.p0.iil1.idi_0.code_read_pref',indigo_stats_file)
            rtl_reads_data = (data_read_rtl + data_read_pref_rtl + data_read_pte_rtl + code_read_rtl + code_read_pref_rtl)*64

        
            core0_rtl_reads_data = rtl_reads_data/core0_rtl_cycles
            core1_rtl_reads_data = rtl_reads_data/core1_rtl_cycles

            rtl_writes_data = find_stat('S*.p0.iil1.idi_0.go_writepull',indigo_stats_file)
            core0_rtl_writes_data = rtl_writes_data/core0_rtl_cycles
            core1_rtl_writes_data = rtl_writes_data/core1_rtl_cycles

            data_read_coho = find_stat('S*.p0.iil1.idi_0.data_read',coho_skl_file)
            data_read_pref_coho = find_stat('S*.p0.iil1.idi_0.data_read_pref',coho_skl_file)
            data_read_pte_coho = find_stat('S*.p0.iil1.idi_0.data_read_pte',coho_skl_file)
            code_read_coho = find_stat('S*.p0.iil1.idi_0.code_read',coho_skl_file)
            code_read_pref_coho = find_stat('S*.p0.iil1.idi_0.code_read_pref',coho_skl_file)
            coho_reads_data = (data_read_coho + data_read_pref_coho + data_read_pte_coho + code_read_coho + code_read_pref_coho)*64

            	    
        
            core0_coho_reads_data = coho_reads_data/core0_coho_cycles
            core1_coho_reads_data = coho_reads_data/core1_coho_cycles

            coho_writes_data = find_stat('S*.p0.iil1.idi_0.go_writepull',coho_skl_file)	
            core0_coho_writes_data = coho_writes_data/core0_coho_cycles
            core1_coho_writes_data = coho_writes_data/core1_coho_cycles            

            if args.dcm:
                coho_reads = coho_reads_data/max(core0_coho_cycles,core1_coho_cycles)
                rtl_reads = rtl_reads_data/max(core0_rtl_cycles,core1_rtl_cycles)
                coho_writes = coho_writes_data/max(core0_coho_cycles,core1_coho_cycles)
                rtl_writes = rtl_writes_data/max(core0_rtl_cycles,core1_rtl_cycles)

                results.append([trace,coho_reads,coho_writes,rtl_reads,rtl_writes])
            else:
                results.append([trace,core0_coho_reads_data,core1_coho_reads_data,core0_rtl_reads_data,core1_rtl_reads_data,
                                  core0_coho_writes_data,core1_coho_writes_data,core0_rtl_writes_data,core1_rtl_writes_data])
if args.dcm:
    display(results,display_option='dcm')
else:
    display(results,display_option='single_core')
  
    






