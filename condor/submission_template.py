import os
import subprocess
import sys
import yaml
from collections import OrderedDict

listOfDataSets = []
listOutputDir = []
yaml_file_dict = {}

with open('file_list.yaml') as yaml_f:
	try:
	    	yaml_file_dict = yaml.safe_load(yaml_f)
		#yam_file_list = yaml_file.split("\n")
	except yaml.YAMLError as exc:
		print(exc)

for key in yaml_file_dict.keys():
	listOfDataSets.append(yaml_file_dict[key][0])
	listOutputDir.append(yaml_file_dict[key][1])

#print(listOfDataSets)
#print(listOutputDir)

dasGoCommand = 'dasgoclient -query="file dataset=file dataset={DATASET} {INSTANCE}"'

queueStr = ""


with open('condor_submission', 'r') as condor_f:
	condor_sub_file = condor_f.read()
	condor_sub_file = condor_sub_file.replace('EXE','executable_jobs.sh')
	#condor_sub_file = condor_sub_file.replace('QUEUE',queueStr)

with open('condor_submission_new', 'w+') as condor_f_new: 
	condor_f_new.write(condor_sub_file)


	for index, datasetName in enumerate(listOfDataSets):
		if 'USER' in datasetName:
			listOfFiles = subprocess.check_output(dasGoCommand.format(DATASET=datasetName,INSTANCE="instance=prod/phys03"), shell = True).split("\n")[:-1]
		else:	
			#dasGoCommand += "
			listOfFiles = subprocess.check_output(dasGoCommand.format(DATASET=datasetName,INSTANCE=" "), shell = True).split("\n")[:-1]
		
		
		for i, inFile in enumerate(listOfFiles):
		
			if os.path.exists('/pnfs/desy.de/cms/tier2/'+inFile):
				#print("capocchione al sugo ", inFile)
				condor_f_new.write('arguments = "-i /pnfs/desy.de/cms/tier2/'+inFile+' /nfs/dust/cms/user/gmilella/'+listOutputDir[index]+'"\n')
				
			else:
				#print("salsizz arraganat ", inFile)
				condor_f_new.write('arguments = "-i root://cms-xrd-global.cern.ch/'+inFile+' /nfs/dust/cms/user/gmilella/'+listOutputDir[index]+'"\n')
				
			condor_f_new.write('Output       = /nfs/dust/cms/user/gmilella/'+listOutputDir[index]+'/log/log.$(Process).out\nError       = /nfs/dust/cms/user/gmilella/'+listOutputDir[index]+'/log/log.$(Process).err\nLog       = /nfs/dust/cms/user/gmilella/'+listOutputDir[index]+'/log/log.$(Process).log\nqueue\n')

			
			if not os.path.isdir("/nfs/dust/cms/user/gmilella/"+listOutputDir[index]):
				os.makedirs("/nfs/dust/cms/user/gmilella/"+listOutputDir[index])
				os.makedirs("/nfs/dust/cms/user/gmilella/"+listOutputDir[index]+"/log")
			else: continue
			
				
			

		

     
