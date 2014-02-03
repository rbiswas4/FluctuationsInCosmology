#!/usr/bin/env python
#Module for IO functionality that is relatively general. The functions in 
#this module are tested versions of ioutilst.py in sntools.
#R. Biswas,  Thu Jan 30 19:42:19 CST 2014

def tokenizeline (line 
	, delimstrings = "" 
	, ignorestrings = ["#"]):
	"""splits the string line into two substrings before and after the 
	first instance of a string in the list ignorestrings, and returns
	a tuple of a list of tokens obtained by tokenizing the first 
	substring on the delimiter delimstrings and the second substring.
	args:
		line: mandatory, string 
			string to tokenize
		delimstrings: optional, defaults to ""
			string of characters (other than whitespace) to 
			be used as a delimiter for tokenizing the line. 
			Example: for a line of TSV, "\t"
		ignorestrings: optional, defaults to ["#"]
			list of strings, occurances of any of which in a 
			line indicates the remainder of the line is a 
			comment which should not be tokenized
	returns:
		tuple:
			list of token strings 
			string of comments
	example usage:
		(lst , comment ) = tokenizeline(myline, delimstrings = ".")
	status: 
		Tested, seems to work correctly,
		testio.py
		#Section: Test tokenization:
		R. Biswas, July 17, 2012
		Rewritten to work for multiple ignorestrings in list to fix bug,
		R. Biswas, Sep 15, 2012
	Notes:
		TODO: allow multiple delimiter strings. 
	"""	
	tokens=[]
	comments = ''

	tmp = line.strip()
	if tmp: 
		minlengthforst = -1
		actualignorestring = None
		lengthofline = len(tmp)

			#Find the ignore string that occurs first

		for st in ignorestrings:
			linelist = tmp.split(st)
			lengthforst = len(linelist[0])
			if lengthforst  < lengthofline:

			#These strings are on the line
				if lengthforst < minlengthforst or -1 == minlengthforst:
					actualignorestring = st
					minlengthforst = lengthforst 

		tokstring = ""

		if actualignorestring:	
			linelist = tmp.split(actualignorestring)
			if len(linelist[1])>1:
				comments = actualignorestring + actualignorestring.join(linelist[1:])
			tokstring = linelist[0]
		else:
			tokstring = tmp
		if delimstrings== "":
			tokens = tokstring.split()
		else:
			#print "delimstring " , delimstrings
			tokens = map(lambda x: x.strip(), tokstring.split(delimstrings))
	ret = ( tokens , comments)
	return ret

def builddict(fname,
	ignorestrings=['#'],
	dictdelim='=',
	startblock = None, 
	endblock =None):

	"""builddict (fname) reads in the file with filename
	fname, and builds a dictionary of keys vs values from
	it
	args: 
		fname:	mandatory, string 
			filename from which the dictionary is to be built
		ignorestring: optional, string, defaults to ["#"]
			list of strings, after which the remaining part of 
			the line should be ignored. 
		dictdelim: optional, string, defaults to '='
			delimiter used to separate keys, values
			in building the dictionary
		startblock = optional, string, defaults to None
			Can do a replace within only the starting and ending
			blocks but both must be provided. These blocks can 
			start with a comment string 
		endblock = string, optional, defaults to None
			Can do a replace within only the starting and ending
			blocks but both must be provided. These blocks can 
			start with a comment string 
	returns:
		dictionary of keys and values (in strings)
	example usage :
		builddict ( fname)  
	status: 
		Seems to work correctly, tested on CAMB params.ini,
		R. Biswas, July 08, 2012 
		That was in configdict. Rewritten to use ioutilst, not tested 
		yet,
		R. Biswas, Aug 09, 2012
	"""
	f = open(fname, "r")
	line = f.readline()
	i = 0
	
	#print ignorestrings
	paramdict={}
	readin = False
        while line != '': 
		if startblock: 
			if (readin ==False):
				if line.find(startblock) !=-1:
					readin = True
		else:
			readin  =True
		if readin == False:
			line = f.readline()
			continue	
	#while line != '':
		tmp = tokenizeline(line, ignorestrings = ignorestrings , 
			delimstrings = dictdelim)
		
		#print line , tmp
		tmp = tmp[0]
		if len(tmp) >1:
				key = tmp[0].strip()
				#print key, tmp
				val = tmp[1].strip()
				paramdict[str(key)] = str(val)  
		line=f.readline()
		if endblock and line.find(endblock) !=-1:
			readin = False
			#print "FOUND ENDBLOCK"
			continue
	
	f.close()
	return paramdict

if __name__ == "__main__":

	myline = "KJAHS KH AKJHS jjhJH. JH HJ   JHH JH #tests "

	(lst , comment ) = tokenizeline(myline, delimstrings = ".")
	print lst 
	print comment

	print "######################################################\n"
	print "######################################################\n"
	print "######################################################\n"
	print "######################################################\n"

	
	print "Test build dict"
	haccdict = builddict (fname  = "example_data/indat.params",
		dictdelim = " ")
		
	print haccdict
