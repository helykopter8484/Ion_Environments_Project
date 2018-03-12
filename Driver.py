import os
import sys
import time
import glob
import pypdb
import shutil
import subprocess
import Subprocess
import Calculations
from PDB import *

pairDict = {
	'BR' : 'Bromium',
	'CA' : 'Calcium',
	'CD' : 'Cadmium',
	'CO' : 'Cobalt',
	'CU' : 'Copper',
	'CL' : 'Chlorine',
	'FE' : 'Iron',
	'K' : 'Potassium',
	'NA' : 'Sodium',
	'HG' : 'Mercury',
	'MG' : 'Magnesium',
	'MN' : 'Manganese',
	'ZN' : 'Zinc',
	'NI' : 'Nickel'
}


# Sequence Identity Parameter: Can be changed.
# Available alternatives: all, 100, 95, 90, 70 - as strings.
seqidCutoff = 'all'

# Setting Home Directory for the Project
home = os.path.dirname(os.path.realpath(__file__)) + "/"
os.chdir(home)

# # System dependent
seqIDPath = '/home/emerald/Documents/Ion_Environments_2018/Sequence_Identities/'
pp_file = '/home/emerald/Java_PDC/potential_1417.out'


Head_List = ['Quad', 'Quad', 'Typ', 'V1', 'V2', 'V3', 'V4', 'SumD', 'Chains', 'L1', 'L2', 'L3', 'L4', 
'L5', 'L6', 'SumL', 'AvgL', 'DevL', 'DevT', 'Volume', 'TF1', 'TF2', 'TF3', 'TF4', 'SStruct']

New_List = ['PDBID', 'QuadA', 'QuadB', 'Typ', 'V1', 'V2', 'V3', 'V4', 'SumD', 'Chains', 'L1', 'L2', 'L3', 'L4', 
'L5', 'L6', 'SumL', 'AvgL', 'DevL', 'DevT', 'Volume', 'TF1', 'TF2', 'TF3', 'TF4', 'SStruct']


for key, value in pairDict.items():
	abbr = key
	name = value
	seqIDFile = seqIDPath + name + '/' + abbr + '_' + seqidCutoff + '.txt'
	Fh = open(seqIDFile, 'r')
	PDBPool = [i.strip() for i in Fh.readlines()]
	Fh.close()

	PDBPool = PDBPool[0:100]

	# System Independent Folders: created in the loop.
	downloaded_pdb = home + name + '/' + '01_Downloaded_PDB/'
	required_chains = home + name + '/' + '02_Required_Chains/'
	filter_output = home + name + '/' + '03_Filter/'
	good_chains = home + name + '/' + '04_Good_Chains/'
	multi_ion_chains = home + name + '/' + '05_Multiple_Ion_Chains/'
	potential_profile_output = home + name + '/' + '06_Potential_Profile/'
	simplex_output = home + name + '/' + '07_Simplex/'
	final_data = home + name + '/' + '08_Final_Data/'
	FreeSASA_Input = home + name + '/' + '09_FreeSASA_Input/'
	FreeSASA_Output = home + name + '/' + '010_FreeSASA_Output/'
	AAAB14_dataFile = final_data + name + "_AAAB14.csv"
	AAAB14_SimCount = final_data + name + "_AAAB14_SimCount.csv"
	AAAB14_Final = final_data + name + "_AAAB14_Final.csv"

	os.mkdir(home + name + '/')
	os.mkdir(downloaded_pdb)
	os.mkdir(required_chains)
	os.mkdir(filter_output)
	os.mkdir(good_chains)
	os.mkdir(multi_ion_chains)
	os.mkdir(potential_profile_output)
	os.mkdir(simplex_output)
	os.mkdir(final_data)
	os.mkdir(FreeSASA_Input)
	os.mkdir(FreeSASA_Output)


	Information = []
	chainGroups = []
	for item in PDBPool:
		item = item.lower()
		ClassPDB = PDB(item, False, downloaded_pdb, 'pdb', True, abbr)
		ClassPDB.GetPDB()
		entry = 'pdb' + item + '.ent'
		RawPDBPath = downloaded_pdb + item
		chains = ClassPDB.GetChainNames(abbr)
		reqChains, chainGroup = ClassPDB.GetUniqueChains(downloaded_pdb, item, chains)
		chainGroups.append(chainGroup)
		info = item, chains, reqChains, chainGroup
		Information.append(info)
		for chain in reqChains:
			chain = chain.split('_')
			pdb, chain = chain[0], chain[1]
			outfname = pdb + "_" + chain + '.pdb'
			outPDB = required_chains + outfname
			ClassPDB = PDB(pdb, False, downloaded_pdb, 'pdb', True, abbr)
			ClassPDB.ChainExtractor(pdb, downloaded_pdb, chain, outPDB)

	ReqChains_Folder = glob.glob(required_chains + '/*.pdb')
	for i in ReqChains_Folder:
		outname = i.split('/')
		chain = outname[-1][5]
		outname = outname[-1][0:6] + '.filter'
		outpath = filter_output + outname
		ClassFilter = Subprocess.ExecuteFilter(i, chain, outpath, pp_file)
		ClassFilter.Filter()

	FilteredChains_Folder = glob.glob(filter_output + '/*.filter')
	for item in FilteredChains_Folder:
		with open(item) as fobj:
			text = fobj.read().replace('\n', '')
			if "BAD:" in text:
				continue
			else:
				item = item.split('/')
				item = item[-1]
				item = item.split('.')
				shutil.copy(required_chains + item[0] + '.pdb', good_chains)

	GoodChains_Folder = glob.glob(good_chains + '/*.pdb')
	for item in GoodChains_Folder:
		e = item.split('/')
		e = e[-1]
		c = e[5]
		PDB.MultiHETATM(e, c, item, multi_ion_chains, abbr)

	MultiHET_Folder = glob.glob(multi_ion_chains + '/*.pdb')
	for multiHet in MultiHET_Folder:
		outname = multiHet.split('/')
		chain = outname[-1][5]
		profile_Outname = outname[-1].replace('.pdb', '.profile')
		simplex_Outname = outname[-1].replace('.pdb', '.simplex')
		profile_Outpath = potential_profile_output + profile_Outname
		simplex_Outpath = simplex_output + simplex_Outname
		ClassPP = Subprocess.ExecutePotentialProfile(multiHet, chain, profile_Outpath, pp_file)
		ClassPP.PotentialProfile()
		ClassSim = Subprocess.ExecuteSimplex(multiHet, chain, simplex_Outpath, pp_file)
		ClassSim.Simplex()

	MultiIonsIn = glob.glob (multi_ion_chains + '/*.pdb')
	Calculations.CreateSASAInput (MultiIonsIn, FreeSASA_Input, abbr)
	FreeSASA_In = glob.glob(FreeSASA_Input + '/*.pdb')
	Calculations.CreateSASA (FreeSASA_In, FreeSASA_Output)

	Simplex_Folder = glob.glob(simplex_output + '/*.simplex')
	AAAB14 = Calculations.AAAB14(Simplex_Folder, abbr, Head_List, New_List, AAAB14_dataFile, AAAB14_SimCount)
	AAAB14 = Calculations.AvgTF(AAAB14)
	AAAB14 = Calculations.Height(AAAB14)
	AAAB14 = Calculations.Equaliterality(AAAB14)
	AAAB14 = Calculations.Omega(AAAB14)
	AAAB14 = Calculations.Triangle1(AAAB14)
	AAAB14 = Calculations.Triangle2(AAAB14)
	AAAB14 = Calculations.Triangle3(AAAB14)
	AAAB14 = Calculations.Triangle4(AAAB14)
	AAAB14 = Calculations.SumSA(AAAB14)
	AAAB14 = Calculations.AEI(AAAB14)
	AAAB14 = Calculations.EIS(AAAB14)
	AAAB14 = Calculations.FLP(AAAB14)
	AAAB14 = Calculations.Triplet(AAAB14)
	AAAB14 = Calculations.Header_Data (AAAB14, downloaded_pdb, name)
	AAAB14 = Calculations.InsertSASA(AAAB14, FreeSASA_Output, name, abbr)
	Calculations.List_to_CSV(AAAB14, AAAB14_Final)
