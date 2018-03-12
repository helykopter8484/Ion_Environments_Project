import os
import amino_dict
import copy
from Bio import PDB
from Bio.PDB import PDBList, PDBParser
from collections import Counter
from itertools import groupby
from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.PDBIO import Select


class PDB(object):

	def __init__(self, pdbID, isobsolete, pdir, file_format, is_overwrite, ion_abbr):
		self.pdbID = pdbID
		self.isobsolete = isobsolete
		self.pdir = pdir
		self.file_format = file_format
		self.is_overwrite = is_overwrite
		self.ion_abbr = ion_abbr


	def GetPDB(self):
		""" Download the pdb file from the server.
			Format will be pdb1crn.ent"""
		pdbl = PDBList()
		pdbl.retrieve_pdb_file(self.pdbID, self.isobsolete, self.pdir,
			self.file_format, self.is_overwrite)


	def PDBParser(self):
		""" Returns a structure object for a given pdb id."""
		parser = PDBParser()
		structure = parser.get_structure(self.pdbID, 'pdb' + self.pdbID + '.ent')
		return structure


	def GetChainNames(self, ion_abbr):
		""" Returns names of chains that contain the required HETATM as a list."""
		Fh = open(self.pdir + 'pdb' + self.pdbID + '.ent', 'r')
		lines = [i.strip() for i in Fh.readlines()]
		chain = []
		for item in lines:
			if item.startswith('HETATM') and ion_abbr in item[12:15]:
				chain.append(item[21])
		chain = set(chain)
		return sorted(list(chain))


	def GetUniqueChains (self, pdir, pdbID, chains_to_check):
		""" Returns a List Unique Chains based on the C-alpha atom information.
			Structure based, not sequence based """
		e = 'pdb' + self.pdbID + '.ent'
		BioParser = PDBParser(PERMISSIVE=True, QUIET = True)
		BioStructure = BioParser.get_structure (self.pdbID, pdir + 'pdb' + self.pdbID + '.ent')
		BioModel = BioStructure[0]
		Chain_AtomSeq = []
		listMatches = []
		for item in chains_to_check:
			pdbid_chain =  e[3:7] + '_' + item
			BioChain = BioModel[item]
			residues = []
			for residue in BioChain:
				for atom in residue:
					if atom.name == 'CA':
						aa1 = amino_dict.replace_all (residue.resname, amino_dict.one_letter)
						residues.append(aa1)
			req_res = [x for x in residues if x in amino_dict.amino]
			atom = "".join(req_res)
			Chain_AtomSeq.append((pdbid_chain, atom))

		Chain_Dict = {}
		for k,v in Chain_AtomSeq:
			Chain_Dict.setdefault(k, v)
		# print (Chain_Dict)

		allChains = [i for i in Chain_Dict.values()]
		set_allChains = list(set(allChains))
		# print (set_allChains)
		groups = {}
		for k, v in Chain_Dict.items():
			groups.setdefault(v, []).append(k)
		matches = {
			k: v 
			for k, v in groups.items()
			}
		list_of_matches = [i for i in matches.values()]
		# print (list_of_matches)
		listMatches.append(list_of_matches)
		req_matches = [i[0] for i in matches.values()]
		return sorted(req_matches), sorted(list_of_matches)


	def DescribePDB (self, pdbID):
		""" Presents a detailed description of a given PDB file. """
		entry = 'pdb' + self.pdbID + '.ent'
		BioParser = PDBParser(PERMISSIVE=True, QUIET = True)
		BioStructure = BioParser.get_structure (entry[3:7], entry)
		for model in BioStructure.get_models():
			print("model", model, "has {} chains".format(len(model)))
			for chain in model:
				print(" - chain ", chain, "has {} residues".format(len(chain)))
				for residue in chain:
					print ("- residue", residue.get_resname(), "has {} atoms".format(len(residue)))
					for atom in residue:
						x,y,z = atom.get_coord()
						print("- atom:", atom.get_name(), "x: {} y:{} z:{}".format(x,y,z))


	def ChainExtractor(self, pdbID, pdir, reqChain_id, outLoc):
		""" Extracts a specific chain from a given PDB file. """
		io = PDBIO()
		BioParser = PDBParser(PERMISSIVE=True, QUIET = True)
		BioStructure = BioParser.get_structure(pdbID, self.pdir + 'pdb' + pdbID + '.ent')
		BioStructure = BioStructure[0]
		for chain in BioStructure.get_chains():
			if chain.id == reqChain_id:
				io.set_structure(chain)
				io.save(outLoc)

	
	def MultiHETATM (pdbID, chainID, pdir, outLoc, abbr):
		""" Creates unique PDB files for the specified HETATM based on the
			number of HETATM entries in the specific chain.
			Specific to the Ion Environments project. """
		datafile = open (pdir, 'r')
		ID = [ p.strip() for p in datafile.readlines() ]
		req_chain_id = 'A' # Converts all non A chains to A for AAAX.
		ATOM = []
		NEW_ATOM = []
		HETATM = []
		for item in ID:
			if item.startswith ('ATOM'):
				ATOM.append (item)
			elif item.startswith('HETATM') and abbr in item[12:15]:
				HETATM.append(item)
		for item in ATOM:
			if item[21] != req_chain_id:
				new_item = item[0:21] + req_chain_id + item[22:]
				NEW_ATOM.append(new_item)
			else:
				NEW_ATOM.append(item)
		NEW_HET = []
		for i in HETATM:
			new_i = ('ATOM   '+ i[7:11] + '  CA  GLY ' + 'B' + '   1      ' + i[32:])
			NEW_HET.append (new_i)
		for index, val in enumerate (NEW_HET):
			index = index + 1
			newfname = pdbID[0:6] + '_' + str(index) + '.pdb'
			current = copy.deepcopy (NEW_ATOM)
			current.append ('TER')
			current.append (val)
			current.append ('TER')
			current.append('END')
			with open (outLoc + newfname, 'w') as FH:
				for i in current:
					FH.write (i + '\n')
				FH.close