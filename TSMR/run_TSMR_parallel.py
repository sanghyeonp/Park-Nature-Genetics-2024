from multiprocessing import Pool, Manager
import subprocess
import os
import pandas as pd
import code
# code.interact(local=dict(globals(), **locals()))


def run_bash(bash_cmd):
    """
    Run bash command.
    Return a list containing standard output, line by line.  
    """
    popen = subprocess.Popen(bash_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, _ = popen.communicate()
    return str(stdout, 'utf-8').strip().split('\n')


def mr(exposure_name, outcome_name,
        harmonized_rdsfile, savedir_individual):
    print("MR running for EXPOSURE ({}) and OUTCOME ({})...".format(exposure_name, outcome_name))
        
    r_script = "func_TSMR.R"
    cmd = "Rscript {} --exposure_name '{}' --outcome_name '{}' --harmonized_rdsfile {} --savedir_individual {}".format(
        r_script, exposure_name, outcome_name, harmonized_rdsfile, savedir_individual)
    
    _ = run_bash(cmd)


def mr_parallel(exposure_name_list, outcome_name_list, harmonized_rdsfile_list, savedir_individual_list, fnc, cores=1):
    pool = Pool(processes=cores)
    inputs = [(exposure_name, outcome_name_list[idx], harmonized_rdsfile_list[idx], savedir_individual_list[idx]) for idx, exposure_name in enumerate(exposure_name_list)]
    pool.starmap(fnc, inputs)
    pool.close()
    pool.join()


if __name__ == "__main__":
    metadata_file = "../01_metadata_prep/TSMR_metadata_for_analysis.tsv"
    savedir_individual = "./tsmr_result_individual"
    harmonized_rdsfile_dir = "./harmonized_rds"
    
    df = pd.read_csv(metadata_file, sep="\t", index_col=False)
    
    outcome_name_list = df['Phenotype from PRS-PheWAS'].tolist()
    exposure_name_list = ['Metabolic syndrome'] * len(df)
    harmonized_rdsfile_list = [os.path.join(harmonized_rdsfile_dir, s+".harmonized.rds") for s in df['File'].tolist()]
    savedir_individual_list = [savedir_individual] * len(df)
    
    # code.interact(local=dict(globals(), **locals()))
    
    mr_parallel(exposure_name_list, outcome_name_list, harmonized_rdsfile_list, savedir_individual_list, fnc=mr, cores=12)
