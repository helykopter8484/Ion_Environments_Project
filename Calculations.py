
import csv
import glob
import numpy
import math
import pypdb
import amino_dict
import subprocess
import urllib.request
from math import *
from decimal import Decimal as dc
import xml.etree.ElementTree as ET
from Bio.PDB import PDBParser, PDBIO, parse_pdb_header


def AAAB_to_File(Outloc, Fname, DataList):
	""" Writes List to file	"""
	Fw = open(loc + fname, 'w')
	for item in List:
		item = ', '.join(str(e) for e in item)
		Fw.write(item + '\n')
	Fw.close()


def List_to_CSV(DataList, OutLocData):
	""" Writes List to CSV File """
	myfile = open(OutLocData,'w')
	wr = csv.writer(myfile, delimiter=';')
	wr.writerows(line for line in DataList)
	myfile.close()


def GetInd(Lst, item):
	""" Gets the index of a specific string
		from a list of strings. """
	return Lst.index(item)


def TriangleArea(A, B, C):
	""" Calculates area of a triangle
		using lengths of sides as input."""
	A, B, C = dc(A), dc(B), dc(C)
	S = (A+B+C) / 2
	SA = S - A
	SB = S - B
	SC = S - C
	a = ( S * (SA * SB * SC) )
	area = numpy.sqrt(a)
	Area = round(area, 4)
	return Area


def FindEC(chn, key, dictionary):
	""" Scans a Nested Dictionary and returns
		the EC Number associated with the input
		file. """
	try:
		for k, v in dictionary.items():
			if isinstance(v, dict):
				if chn in v['chain']:
					ec = v['ec_number']
	except KeyError:
		ec = 'NA'
	return ec


def CreateSASAInput (mult_loc, SASAIn_loc, ion_abbrev):
	""" Creates PDB file suitable for FreeSASA input
		Specific for Ion Environments project. """
	for inpath in mult_loc:
		entry = inpath.split('/')
		entry = entry[-1]
		dotloc = entry.find('.')
		identifier = entry[0:dotloc]
		outpath = SASAIn_loc + identifier + '.pdb'
		a = identifier.split('_')
		pdb, chain, num  = a[0], a[1], a[2]
		srt = ' ' + ion_abbrev + ' A'
		datafile = open(inpath, 'r')
		newItem = []
		ID = [ p.strip() for p in datafile.readlines() ]
		for item in ID:
			if "GLY B" in item:
				inter = item
			else:
				newItem.append(item)
		newItem = newItem[0:-3]
		y = inter.replace('ATOM  ', 'HETATM')
		z = y.replace('GLY B', srt)
		zz = z.replace(z[12:15], ion_abbrev + ' ')
		newItem.append(zz)
		newItem.append('TER')
		newItem.append('END')

		with open (outpath, 'w') as FH:
			for i in newItem:
				FH.write (i + '\n')
			FH.close()


def CreateSASA (sasa_in, sasa_out):
	""" Subprocess call to FreeSASA program. """
	for inpath in sasa_in:
		entry = inpath.split('/')
		entry = entry[-1]
		dotloc = entry.find('.')
		identifier = entry[0:dotloc]
		outpath = sasa_out + identifier + '.sasa'
		arglist = ["freesasa", "--format=pdb", "--hetatm", inpath, "-o", outpath]
		process = subprocess.Popen(arglist, stdout=subprocess.PIPE)
		stdout, stderr = process.communicate()


#################################################################################


def AAAB14(listFilePaths, ion_abbr, HeadList, NewList, OutLocDataFile, OutLocCountFile):
	""" Returns a list of lists where the vertex is AAAB abd the edge lengths
		are less than or equal to 14 A.
		Specific for Ion Environments project. """
	AAAB14 = []
	Counts = []
	for item in listFilePaths:
		entry = item.split('/')
		entry = entry[-1]
		pdbid = entry[0:4]
		chain = entry[5]
		dotloc = entry.index('.')
		pdbid_het = entry[0:dotloc]
		L1 = (HeadList.index('L1')) 
		L2 = (HeadList.index('L2')) 
		L3 = (HeadList.index('L3')) 
		L4 = (HeadList.index('L4')) 
		L5 = (HeadList.index('L5')) 
		L6 = (HeadList.index('L6')) 

		Fh = open (item, "r")
		lines = [i.strip () for i in Fh.readlines()]
		lineCount = 0
		AAABCount = 0
		AAAB14Count = 0
		for line in lines:
			lineCount = lineCount + 1
			if "AAAB" in line:
				AAABCount = AAABCount + 1
				line = line.split()
				if (dc(line[L1]) <= 14.00) and (dc(line[L2]) <= 14.00) and (dc(line[L3]) <= 14.00) and (dc(line[L4]) <= 14.00) and (dc(line[L5]) <= 14.00) and (dc(line[L6]) <= 14.00):
					line.insert(0,pdbid_het)
					AAAB14Count = AAAB14Count + 1
					AAAB14.append(line)
		Fh.close()
		Cnt_List = [pdbid_het, lineCount, AAABCount, AAAB14Count]
		Counts.append(Cnt_List)
	AAAB14.insert(0, NewList)
	Counts.insert(0, ["PDBID", "Line_Count", "AAABCount", "AAAB14Count"])
	List_to_CSV(AAAB14, OutLocDataFile)
	List_to_CSV(Counts, OutLocCountFile)
	return AAAB14


def AvgTF(AAAB14):
	""" Calculates the Average of Temperature Factors.
		Specific for Ion Environments project. """
	Data = []
	Head = AAAB14[0]
	Tail = AAAB14[1:]
	# First get Index and then insert into header
	TF1 = GetInd(Head, 'TF1')
	TF2 = GetInd(Head, 'TF2')
	TF3 = GetInd(Head, 'TF3')
	TF4 = GetInd(Head, 'TF4')
	SStruct = GetInd(Head, 'SStruct')
	Head.insert(SStruct, 'AvgTF')
	Data.append(Head)
	for item in Tail:
		A, B, C, D = dc(item[TF1]), dc(item[TF2]), dc(item[TF3]), dc(item[TF4])
		Avg = (A + B + C + D) / dc(4.0)
		Avg = round(Avg, 4)
		item.insert(SStruct, str(Avg))
		Data.append(item)
	return Data


def Height(AAAB14):
	""" Calculates the Height of Tetrahedra.
		Specific for Ion Environments project. """
	Data = []
	Head = AAAB14[0]
	Tail = AAAB14[1:]
	L1 = GetInd(Head, 'L1')
	L2 = GetInd(Head, 'L2')
	L4 = GetInd(Head, 'L4')
	Vl = GetInd(Head, 'Volume')
	SStruct = GetInd(Head, 'SStruct')
	Head.insert(SStruct, 'Height')
	Data.append(Head)
	for item in Tail:
		A = dc(item[L1])
		B = dc(item[L2])
		C = dc(item[L4])
		Vol = dc(item[Vl])
		S = (A + B + C) / dc(2)
		SA = S - A
		SB = S - B
		SC = S - C
		a = (S*(SA*SB*SC))
		area = numpy.sqrt(a)
		ht = (3*Vol)/area
		Height = round(ht,4)
		item.insert(SStruct,str(Height))
		Data.append(item)
	return Data


def Equaliterality(AAAB14):
	""" Calculates the Equaliterality of the Tetrahedra.
		Specific for Ion Environments project. """
	Data = []
	Head = AAAB14[0]
	Tail = AAAB14[1:]
	L1 = GetInd(Head, 'L1')
	L2 = GetInd(Head, 'L2')
	L4 = GetInd(Head, 'L4')
	SStruct = GetInd(Head, 'SStruct')
	Head.insert(SStruct, 'Equilaterality')
	Data.append(Head)
	for item in Tail:
		A = dc(item[L1])
		B = dc(item[L2])
		C = dc(item[L4])
		# print A, B, C
		Mn = (A + B + C) / dc(3)
		SqMn = Mn ** 2
		AA = (A-B) ** 2
		BB = (A-C) ** 2
		CC = (B-C) ** 2
		Num = AA + BB + CC
		Den = dc(3) * SqMn
		ans = Num / Den
		Equa = round(ans,4)
		item.insert(SStruct, str(Equa))
		# print (item)
		Data.append(item)
	return Data


def Omega(AAAB14):
	""" Calculates the Solid Angle Omega for the Vertex that represents
		the HETATM.
		Specific for Ion Environments project. """
	Data = []
	Head = AAAB14[0]
	Tail = AAAB14[1:]
	L1 = GetInd(Head, 'L1')
	L2 = GetInd(Head, 'L2')
	L3 = GetInd(Head, 'L3')
	L4 = GetInd(Head, 'L4')
	L5 = GetInd(Head, 'L5')
	L6 = GetInd(Head, 'L6')
	SStruct = GetInd(Head, 'SStruct')
	Head.insert(SStruct, 'Omega')
	Data.append(Head)
	for item in Tail:
		l1 = dc(item[L1])
		l2 = dc(item[L2])
		l3 = dc(item[L3])
		l4 = dc(item[L4])
		l5 = dc(item[L5])
		l6 = dc(item[L6])
		# print (l1, l2, l3, l4, l5, l6)
		try:
			A = ((l5 ** 2) + (l6 ** 2) - (l4 ** 2))
			a = (2 * l5 * l6)
			A = acos(A/a)
			B = ((l5 ** 2) + (l3 ** 2) - (l1 ** 2))
			b = (2 * l5 * l3)
			B = acos(B/b)
			C = ((l6 ** 2) + (l3 ** 2) - (l2 ** 2))
			c = (2 * l3 * l6)
			C = acos(C/c)
			S = (A + B + C) / 2.0
			val = tan(S / 2) * tan((S-A) / 2) * tan((S-B) / 2) * tan((S-C) / 2)
			if val <= 0:
				val = abs(val)
		except ValueError:
			val = 0
		omega = sqrt(val)
		omega = 4 * atan(omega)
		Omega = degrees(omega)
		Omega = round(omega,4)
		item.insert(SStruct, str(Omega))
		Data.append(item)
	return Data


def Triangle1(AAAB14):
	""" Calculates the Area of the Triangle.
		Specific for Ion Environments project."""
	Data = []
	Head = AAAB14[0]
	Tail = AAAB14[1:]
	L1 = GetInd(Head, 'L1')
	L3 = GetInd(Head, 'L3')
	L5 = GetInd(Head, 'L5')
	SStruct = GetInd(Head, 'SStruct')
	Head.insert(SStruct, 'Triangle1')
	Data.append(Head)
	for item in Tail:
		A = dc(item[L1])
		B = dc(item[L3])
		C = dc(item[L5])
		Area = TriangleArea(A, B, C)
		# print (Area)
		item.insert(SStruct, str(Area))
		Data.append(item)
	return Data


def Triangle2(AAAB14):
	""" Calculates the Area of the Triangle.
		Specific for Ion Environments project."""
	Data = []
	Head = AAAB14[0]
	Tail = AAAB14[1:]
	L1 = GetInd(Head, 'L4')
	L3 = GetInd(Head, 'L5')
	L5 = GetInd(Head, 'L6')
	SStruct = GetInd(Head, 'SStruct')
	Head.insert(SStruct, 'Triangle2')
	Data.append(Head)
	for item in Tail:
		A = dc(item[L1])
		B = dc(item[L3])
		C = dc(item[L5])
		Area = TriangleArea(A, B, C)
		# print (Area)
		item.insert(SStruct, str(Area))
		Data.append(item)
	return Data


def Triangle3(AAAB14):
	""" Calculates the Area of the Triangle.
		Specific for Ion Environments project."""
	Data = []
	Head = AAAB14[0]
	Tail = AAAB14[1:]
	L1 = GetInd(Head, 'L2')
	L3 = GetInd(Head, 'L3')
	L5 = GetInd(Head, 'L6')
	SStruct = GetInd(Head, 'SStruct')
	Head.insert(SStruct, 'Triangle3')
	Data.append(Head)
	for item in Tail:
		A = dc(item[L1])
		B = dc(item[L3])
		C = dc(item[L5])
		Area = TriangleArea(A, B, C)
		# print (Area)
		item.insert(SStruct, str(Area))
		Data.append(item)
	return Data


def Triangle4(AAAB14):
	""" Calculates the Area of the Triangle.
		Specific for Ion Environments project."""
	Data = []
	Head = AAAB14[0]
	Tail = AAAB14[1:]
	L1 = GetInd(Head, 'L1')
	L3 = GetInd(Head, 'L2')
	L5 = GetInd(Head, 'L4')
	SStruct = GetInd(Head, 'SStruct')
	Head.insert(SStruct, 'Triangle4')
	Data.append(Head)
	for item in Tail:
		A = dc(item[L1])
		B = dc(item[L3])
		C = dc(item[L5])
		Area = TriangleArea(A, B, C)
		# print (Area)
		item.insert(SStruct, str(Area))
		Data.append(item)
	return Data


def SumSA (AAAB14):
	""" Calculates the Surface Area of the Tetrahedra.
		Specific for the Ion Environments Project. """
	Data = []
	Head = AAAB14[0]
	Tail = AAAB14[1:]
	ta1 = GetInd(Head, 'Triangle1')
	ta2 = GetInd(Head, 'Triangle2')
	ta3 = GetInd(Head, 'Triangle3')
	ta4 = GetInd(Head, 'Triangle4')
	SStruct = GetInd(Head, 'SStruct')
	Head.insert(SStruct, 'Simplex_SurfArea')
	Data.append(Head)
	for item in Tail:
		A, B, C, D = dc(item[ta1]), dc(item[ta2]), dc(item[ta3]), dc(item[ta4])
		Sum = (A + B + C + D)
		Sum = round(Sum,4)
		item.insert(SStruct, str(Sum))
		Data.append(item)
	return Data


def FLP (AAAB14):
	""" Returns the Reduced Amino Acid Representation for the
		amino acid vertices participating in the Delaunay
		Tessellation. 
		Specific for Ion Environments project."""
	Data = []
	Head = AAAB14[0]
	Tail = AAAB14[1:]
	pos = GetInd(Head, 'QuadB')
	Head.insert(pos+1, 'FLP_Group')
	Head.insert(pos+2, 'Sort_FLP')
	Data.append(Head)
	for item in Tail:
		point = item[pos]
		point = point[0:3]
		Vais = amino_dict.replace_all(point, amino_dict.Red_Vaisman)
		Vais = Vais.upper()
		Vais_sort = ''.join(sorted(Vais))
		item.insert(pos+1, Vais)
		item.insert(pos+2, Vais_sort)
		Data.append(item)
	return Data


def EIS (AAAB14):
	""" Returns the Reduced Amino Acid Representation for the
		amino acid vertices participating in the Delaunay
		Tessellation. 
		Specific for Ion Environments project."""
	Data = []
	Head = AAAB14[0]
	Tail = AAAB14[1:]
	pos = GetInd(Head, 'QuadB')
	Head.insert(pos+1, 'EIS_Group')
	Head.insert(pos+2, 'Sort_EIS')
	Data.append(Head)
	for item in Tail:
		point = item[pos]
		point = point[0:3]
		Vais = amino_dict.replace_all(point, amino_dict.Red_Li)
		Vais = Vais.upper()
		Vais_sort = ''.join(sorted(Vais))
		item.insert(pos+1, Vais)
		item.insert(pos+2, Vais_sort)
		Data.append(item)
	return Data


def AEI (AAAB14):
	""" Returns the Reduced Amino Acid Representation for the
		amino acid vertices participating in the Delaunay
		Tessellation. 
		Specific for Ion Environments project."""
	Data = []
	Head = AAAB14[0]
	Tail = AAAB14[1:]
	pos = GetInd(Head, 'QuadB')
	Head.insert(pos+1, 'AEI_Group')
	Head.insert(pos+2, 'Sort_AEI')
	Data.append(Head)
	for item in Tail:
		point = item[pos]
		point = point[0:3]
		Vais = amino_dict.replace_all(point, amino_dict.Red_Wang)
		Vais = Vais.upper()
		Vais_sort = ''.join(sorted(Vais))
		item.insert(pos+1, Vais)
		item.insert(pos+2, Vais_sort)
		Data.append(item)
	return Data


def Triplet(AAAB14):
	""" Returns the Amino Acid triplet for a given Tetrahedra.
		Specific for the Ion Environments project. """
	Data = []
	Head = AAAB14[0]
	Tail = AAAB14[1:]
	pos = GetInd(Head, 'QuadB')
	Head.insert(pos + 1, 'Triplet')
	Data.append(Head)
	for item in Tail:
		point = item[pos]
		point = point[0:3]
		item.insert(pos + 1, point)
		Data.append(item)
	return Data

def GO_Data(pdbid, req_chain):
	""" Returns the Gene Ontology Information for a specific 
		pdbid and chain."""
	GO_url = 'http://www.rcsb.org/pdb/rest/goTerms?structureId=' + pdbid.upper() + '.' + req_chain.upper()
	local_filename, headers = urllib.request.urlretrieve(GO_url)
	tree = ET.parse(local_filename)
	root = tree.getroot()
	for child in root:
		if child.attrib['structureId'] == pdbid.upper() and child.attrib['chainId'] == req_chain:
			goID = child.attrib['id']
			for c in child:
				goName = c.attrib['name']
	return goID, goName


def Pfam_Data (pdbid, req_chain):
	""" Returns the Pfam Information for a specific 
		pdbid and chain."""
	Pfam_url = 'http://www.rcsb.org/pdb/rest/hmmer?structureId=' + pdbid
	local_filename, headers = urllib.request.urlretrieve(Pfam_url)
	tree = ET.parse(local_filename)
	root = tree.getroot()
	for child in root:
		if child.attrib['structureId'] == pdbid.upper() and child.attrib['chainId'] == req_chain:
			name1 = str(child.attrib['pfamName'])
			acces = str(child.attrib['pfamAcc'])
			start = str(child.attrib['pdbResNumStart'])
			end = str(child.attrib['pdbResNumEnd'])
			desc = str(child.attrib['pfamDesc'])
			evalue = str(child.attrib['eValue'])
	return (name1, acces, start, end, desc, evalue)


def Header_Data (AAAB14, pdbloc, ion_name):
	""" Returns crucial Header Data for a specific pdbid.
		Somewhat non specific to Ion Environments project. """
	Data = []
	Head = AAAB14[0]
	Tail = AAAB14[1:]
	Head.insert(1, 'StructMethod')
	Head.insert(1, 'Resolution')
	Head.insert(len(Head), 'ECNum')
	Head.insert(len(Head), 'Type')
	Head.insert(len(Head), 'Ion')
	Data.append(Head)
	for item in Tail:
		fname = item[0]
		fname = fname.split('_')
		pdbid = fname[0]
		pdb = 'pdb' + pdbid + '.ent'
		chain = fname[1]
		handle = open(pdbloc + pdb,'r')
		header_dict = parse_pdb_header(handle)
		handle.close()
		name = header_dict['name']
		head = header_dict['head']
		method = header_dict['structure_method']
		reso = header_dict['resolution']
		Comp = header_dict['compound']
		ec = FindEC(chain.lower(), 'ec_number', Comp)
		if ec.startswith('1.'):
			txt = 'Oxidoreductase'
		elif ec.startswith('2.'):
			txt = 'Transferase'
		elif ec.startswith('3.'):
			txt = 'Hydrolase'
		elif ec.startswith('4.'):
			txt = 'Lyase'
		elif ec.startswith('5.'):
			txt = 'Isomerase'
		elif ec.startswith('6.'):
			txt = 'Ligase'
		else:
			txt = 'Non_Enzyme'
		item.insert(1, method)
		item.insert(1, reso)
		item.insert(len(item), ec)
		item.insert(len(item), txt)
		item.insert(len(item), ion_name)
		Data.append(item)
	return Data


def SASAParser(fPath, a, b, c, d, ion_abbrev):
	""" Parses the FreeSASA output file and returns the SASA and
		ionic radius of C-alpha and HETATM. """
	Fh = open (fPath, "r")
	lines = [i.strip().split() for i in Fh.readlines()]
	a1 = ''
	a2 = ''
	b1 = ''
	b2 = ''
	c1 = ''
	c2 = ''
	d1 = ''
	d2 = ''
	for line in lines:
		if line[0] == "ATOM" and line[2] == "CA" and line[5] == a:
			a1, a2 = line[-1], line[-2]
		if line[0] == "ATOM" and line[2] == "CA" and line[5] == b:
			b1, b2 = line[-1], line[-2]
		if line[0] == "ATOM" and line[2] == "CA" and line[5] == c:
			c1, c2 = line[-1], line[-2]
		if line[0] == "HETATM" and line[2] == ion_abbrev:
			d1, d2 = line[-1], line[-2]
	return a1,  a2, b1, b2, c1, c2, d1, d2


def InsertSASA(AAAB14, sasaFileLoc, ion_name, ion_abbrev):
	""" Adds the required information obtained from SASAParser to the
		data structure of the ion environments project. """
	Data = []
	Head = AAAB14[0]
	Tail = AAAB14[1:]
	v1 = GetInd(Head, 'V1')
	v2 = GetInd(Head, 'V2')
	v3 = GetInd(Head, 'V3')
	v4 = GetInd(Head, 'V4')
	SStruct = GetInd(Head, 'SStruct')
	Head.insert(SStruct, 'SASA_V4')
	Head.insert(SStruct, 'SASA_V3')
	Head.insert(SStruct, 'SASA_V2')
	Head.insert(SStruct, 'SASA_V1')
	Data.append(Head)
	for item in Tail:
		fileID, V1, V2, V3, V4 = item[0], item[v1], item[v2], item[v3], item[v4]
		inFilePath = sasaFileLoc + fileID + '.sasa'
		a1, a2, b1, b2, c1, c2, d1, d2 = SASAParser(inFilePath, V1, V2, V3, V4, ion_abbrev)
		item.insert(SStruct, d1)
		item.insert(SStruct, c1)
		item.insert(SStruct, b1)
		item.insert(SStruct, a1)
		Data.append(item)
	return Data
