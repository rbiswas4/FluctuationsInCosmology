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
			print "delimstring " , delimstrings
			tokens = map(lambda x: x.strip(), tokstring.split(delimstrings))
	ret = ( tokens , comments)
	return ret

if __name__ == "__main__":

	myline = "KJAHS KH AKJHS jjhJH. JH HJ   JHH JH #tests "

	(lst , comment ) = tokenizeline(myline, delimstrings = ".")
	print lst 
	print comment
