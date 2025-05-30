import os, re, sys
import subprocess
import argparse
import yaml
import stat

parser = argparse.ArgumentParser()

parser.add_argument('--sgn', dest='sgn', action='store_true', help='Signal', default=False)
parser.add_argument('--relwidth', dest='relwidth', type=str, help='Relative Reso Width', default='0')
parser.add_argument('--year', dest='year', type=str, help='Year', default='2017')
parser.add_argument('--isData', dest='isData', action='store_true', help='Data', default=False)

args = parser.parse_args()

listOfDataSets = []
listOutputDir = []
yaml_file_dict = {}

ERAS = {
    '2018': ['A', 'B', 'C', 'D'], '2017': ['B', 'C', 'D', 'E', 'F'],
    '2016': ['F', 'G', 'H'],
    '2016preVFP': ['B', 'C', 'D', 'E', 'F'], #'B"
    '2022': ['C', 'D'], '2022EE': ['E', 'F', 'G']
} # ,

data_sample = ''

MASS_VALUES = ['500', '750', '1000', '1250', '1500', '1750', '2000', '2500', '3000', '4000'] #
relwidth_values = ['4', '10','20', '50'] #'4', 

with open('file_list.yaml') as yaml_f:
    try:
        yaml_file_dict = yaml.safe_load(yaml_f)
    except yaml.YAMLError as exc:
        print(exc)

if args.sgn:
    for width in relwidth_values:
        for mass in MASS_VALUES:
            listOfDataSets.append(yaml_file_dict['sgn'][args.year][0].replace("WIDTH", width).replace("MASS", mass))
            listOutputDir.append(yaml_file_dict['sgn'][args.year][1].replace("WIDTH", width).replace("MASS", mass))
elif args.isData: 
    for data_sample in yaml_file_dict['data'][args.year].keys():
        for era in ERAS[args.year]:
            if 'SingleMuon' in data_sample:
                if args.year == '2016preVFP' and era == 'B':
                    era = era + "2"
                listOfDataSets.append(yaml_file_dict['data'][args.year][data_sample][0].replace("ERA", era))
            else:
                if era == 'E' and args.year == '2022EE':
                    if data_sample == 'DoubleMuon':
                        listOfDataSets.append(yaml_file_dict['data']['2022']['Muon'][0]+'_'+era) # era E is processed with ReReco GT and stored in 2022
                    else:
                        listOfDataSets.append(yaml_file_dict['data']['2022'][data_sample][0]+'_'+era) # era E is processed with ReReco GT and stored in 2022
                else:
                    listOfDataSets.append(yaml_file_dict['data'][args.year][data_sample][0]+'_'+era)
            listOutputDir.append(yaml_file_dict['data'][args.year][data_sample][1]+'_'+era)
else:
    for sample in yaml_file_dict['bkg'][args.year].keys():
        listOfDataSets.append(yaml_file_dict['bkg'][args.year][sample][0])
        listOutputDir.append(yaml_file_dict['bkg'][args.year][sample][1])  

print(listOfDataSets)
print(listOutputDir)


# --- EXECUTABLES FILES
# current_dir = os.getcwd()
# executables_dir = os.path.join(current_dir, 'executable_files', args.year)
# if not os.path.exists(executables_dir):
#     os.makedirs(executables_dir)
# os.chmod(executables_dir, stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)  # equivalent to chmod 777
# subprocess.call(['chmod', '-R', 'a+w+x+r', executables_dir])

# with open('executable_jobs_hotvr_tmp.sh', 'r') as exe_f_tmp:
#     exe_file = exe_f_tmp.read()
# exe_file=exe_file.replace('YEAR', args.year)

# executable_name = 'executable_jobs_' + args.year + '_hotvr.sh'

# with open(executables_dir + '/' + executable_name, 'w') as exe_f:
#     exe_f.write(exe_file)
# os.chmod(executables_dir + '/' + executable_name, 0o777)   
# ----


# ---- CONDOR FILES
with open('condor_submission', 'r') as condor_f:
    condor_sub_file = condor_f.read()
    # executable_name = 'executable_jobs_' + args.year + '_hotvr.sh'
    # condor_sub_file = condor_sub_file.replace('EXE', executables_dir + '/' + executable_name)

command_str, files_str = '', ''
command_dict = dict()

CONDOR_SUBMIT_FILES_DIR = '{}/condor_submit_files'.format(os.getcwd())
if not os.path.isdir(CONDOR_SUBMIT_FILES_DIR):
    print('Making dirs: \n{}'.format(CONDOR_SUBMIT_FILES_DIR))
    os.makedirs(CONDOR_SUBMIT_FILES_DIR)

# --- I/O DIRECTORY PATH
NFS_PATH = os.environ.get('NFS', '/data/dust/user/gmilella')
CONDOR_REPO_LOGFILE = os.path.join(NFS_PATH, 'ttX_ntuplizer/log_files/')
NFS_OUTPUT_REPO = os.path.join(NFS_PATH, 'ttX_ntuplizer/')

PNFS_REPO = os.environ.get('PNFS', '/pnfs/desy.de/cms/tier2/store/user/gmilella')

DAS_COMMAND = f'dasgoclient -query="file dataset=DATASET instance=prod/phys03"'
# ---

output_dir = "bkg_{}_hotvr".format(args.year)
if args.isData: output_dir = "data_{}_hotvr".format(args.year)
elif args.sgn: output_dir = "sgn_{}_central_hotvr".format(args.year)

#command_list = condor_sub_file
#print(command_list)

batch_number = 1 
batch_count = 0

files_batches = {}
files_batch_number = 0
total_files, total_files_per_process = 0, 0 # upper limit on 50k files --> it will result in more than 5k jobs
total_files_in_condor_job = 0

files_for_lxplus = ''

processes_list, processes_dirs_list = [], []

for iprocess, process_dir in enumerate(listOfDataSets):
    processes_list.append(listOutputDir[iprocess])
    processes_dirs_list.append(process_dir)

for process, process_dir in zip(processes_list, processes_dirs_list):
    files_to_process = []

    if 'SingleMuon' in process:
        das_query = DAS_COMMAND.replace("DATASET", process_dir)
        print(f"Querying DAS for {process_dir}: {das_query}")

        try:
            result = subprocess.run(das_query, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            das_files = result.stdout.strip().split("\n")
            das_files = [f"root://cms-xrd-global.cern.ch//{file}" for file in das_files if file]
            print(f"Retrieved {len(das_files)} files from DAS")
            files_to_process.extend(das_files)
        except subprocess.CalledProcessError as e:
            print(f"Error running dasgoclient for {process_dir}: {e.stderr}")
    else:
        if os.path.exists(f"{PNFS_REPO}/{process_dir}"):
            for root, subdirs, files in os.walk(f"{PNFS_REPO}/{process_dir}"):
                if files:
                    files_to_process.extend([os.path.join(root, file) for file in files])
                    print(f"Dataset {process_dir}, number of local files: {len(files)}")

        # Merge DAS and local files
        # all_files = set(das_files + local_files)
        # cms-xrd-global.cern.ch

    if files_to_process:
        if process not in files_batches.keys(): 
            files_batches[process] = {}

        files_batch_number += 1
        command_str = ""
        total_files_per_process = 0
        files_str = ""

        # for root, subdirs, files in os.walk(f'{PNFS_REPO}/{process_dir}'):
        for ifile, file_path in enumerate(files_to_process):
            if ".root" in file_path:
                total_files += 1
                total_files_per_process += 1

                # if args.isData:
                #     if args.year == '2022' or args.year == '2022EE':
                #         # print(file_path)
                #         if 'EGamma' in file_path and '_G' in file_path:
                #             # print(file_path)
                #             if '241217_' not in file_path: continue
                #         elif '241212_' not in file_path: continue

                # if total_files_per_process > 120:
                #     break
                files_str += f" -i {file_path} "
                if ifile % 3 == 0: #processising more file in one single condor job
                    files_batch_number += 1
                    files_str = f" -i {file_path} "
                files_batches[process][files_batch_number] = files_str

                # files_for_lxplus += '{}/{} \n'.format(root.split("/pnfs/desy.de/cms/tier2", 1)[1], file)

        for batch_number in files_batches[process]:
            command_str += f'\narguments = " --year {args.year}'
            command_str += files_batches[process][batch_number]
            command_str += f' {NFS_OUTPUT_REPO}{output_dir}/{process} '

            if args.sgn:
                command_str += ' --isSignal "'
            elif args.isData:
                command_str += ' --isData "'
            else: command_str += ' "'
    
            command_str += f"\nOutput = {CONDOR_REPO_LOGFILE}{output_dir}/{process}/log_{process}$(Cluster).$(Process).out"
            command_str += f"\nError = {CONDOR_REPO_LOGFILE}{output_dir}/{process}/log_{process}$(Cluster).$(Process).err"
            command_str += f"\nLog = {CONDOR_REPO_LOGFILE}{output_dir}/{process}/log_{process}$(Cluster).$(Process).log"
            # command_str += '\n+MyProject = "cms"\nqueue\n'
            command_str += '\nqueue\n'

        base_log_path = f"{CONDOR_REPO_LOGFILE}{output_dir}/{process}"
        os.makedirs(base_log_path, exist_ok=True)
        os.makedirs(os.path.join(NFS_OUTPUT_REPO, f"{output_dir}/{process}"), exist_ok=True)

        if args.sgn: 
            condor_out_name = 'condor_submission_sgn_new_{}_{}.sub'.format(process, args.year)
        elif args.isData: 
            condor_out_name = 'condor_submission_data_new_{}.sub'.format(process)
        else: 
            condor_out_name = 'condor_submission_bkg_new_{}.sub'.format(process)

        with open('{}/{}'.format(CONDOR_SUBMIT_FILES_DIR, condor_out_name), 'w+') as condor_f_new: 
            condor_f_new.write(condor_sub_file)
            condor_f_new.write(command_str)  

    else:
        print(f"Files not found in {process_dir}")
        continue

print('Files in total: {}'.format(total_files))

with open('files_for_lxplus.txt', 'w') as file:
    file.write(files_for_lxplus)
