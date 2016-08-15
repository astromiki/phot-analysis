#! /usr/bin/env python
# -*- coding: utf-8 -*-

#  commentary
## test module

import sys, json

def arguments():
	import argparse
	global filein, fileout
	parser = argparse.ArgumentParser(description="Program converts astrnomical JSON data to a more comfortable format."
																"\n(f.e. Gaia alerts photometry data is being published in this format.)"
																"\n\nRequirements: Python 2.7+"
																"\n\nVersion: 2016.07.02, (PM)", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-f", "--file", type=str, help="JSON file to be processed", required=True)
	parser.add_argument("-o", "--out", type=str, help="output filename (default: <input_file>.data)")
	args = parser.parse_args()
	filein = args.file
	fileout = args.out if args.out > 0 else str(filein).split(".")[0]+".data"

arguments()
with open(filein,'r') as data_file: 
	data = json.load(data_file)
records = len(data['mjd'])
with open(fileout,'w') as data_file:
	data_file.write('{:^12}  {:^6} {:^9}   {:^8}   {:^8}\n'.format("MJD","FILTER","MAG","MAGERR","CALIB_ERR"))
	for i in range (0, records):
		data_file.write('{:f}    {:s}    {:f}   {:f}    {:f}\n'.format(data['mjd'][i], data['filter'][i], data['mag'][i], data['magerr'][i], data['caliberr'][i]))
print filein+" > "+fileout
exit()
