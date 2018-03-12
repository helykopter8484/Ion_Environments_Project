import os
import sys
import subprocess
import shutil


class ExecuteFilter:
	""" Subprocess call to the Filter program written in Java.
		Specific to the Ion Environments project. """
	def __init__ (self, inPath, chain, outPath, ppPath = None):
		self.inPath = inPath
		self.chain = chain
		self.outPath = outPath

	def Filter(self):
		arglist = ["java", "Filter", "-f", self.inPath, "-c", self.chain, "-o", self.outPath, "-v"]
		process = subprocess.Popen(arglist, stdout=subprocess.PIPE)
		stdout, stderr = process.communicate()


class ExecutePotentialProfile:
	""" Subprocess call to the Potential Profile program written in Java.
		Specific to the Ion Environments project. """
	def __init__ (self, inPath, chain, outPath, ppPath = None):
		self.inPath = inPath
		self.chain = chain
		self.outPath = outPath
		self.ppPath = ppPath

	def PotentialProfile(self):
		arglist = ["java", "PotentialProfile", "-p", self.ppPath, "-f", self.inPath, "-c", self.chain, "-o", self.outPath]
		process = subprocess.Popen(arglist, stdout=subprocess.PIPE)
		stdout, stderr = process.communicate()


class ExecuteSimplex:
	""" Subprocess call to the Simplex program written in Java.
		Specific to the Ion Environments project. """
	def __init__ (self, inPath, chain, outPath, ppPath = None):
		self.inPath = inPath
		self.chain = chain
		self.outPath = outPath

	def Simplex(self):
		arglist = ["java", "Simplexn", "-f", self.inPath, "-c", "@", "-o", self.outPath]
		process = subprocess.Popen(arglist, stdout=subprocess.PIPE)
		stdout, stderr = process.communicate()
