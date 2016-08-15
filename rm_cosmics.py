#! /usr/bin/env python

#  commentary
## test module

import sys
import pyfits
import cosmics, f2n
import os
import argparse

global default_keyword, default_iter
default_keyword, default_iter_no = ("GAIN", 3)	# default keyword for FITS header readtime field & number of iterations for removing cosmics

def arguments():
	global test, filename, keyword, listf, iter_no
	parser = argparse.ArgumentParser(description=" Program removes cosmic rays from FITS image."
																"\n Requirements: Python 2.7+, f2n, cosmics.py."
																"\n\n Version: 2016.08.12, (PM)", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-f", "--file", type=str, help="FITS file (list of files) to be processed", required=True)
	parser.add_argument("-i", "--iter", type=int, help="max. number of iterations (default: "+str(default_iter_no)+")")
	parser.add_argument("-l", "--list", help="processing list of files", action="store_true")
	parser.add_argument("-t", "--test", help="output files are 'test.fits' & 'mask.fits'", action="store_true")
	parser.add_argument("-k", "--key", type=str, help="specify own FITS header keyword for readtime (default: 'GAIN')")
	parser.add_argument("-c", "--clean", help="deletes 'mask.fits' & 'test.fits'", action="store_true")
	args = parser.parse_args()
	[os.system("rm test.fits mask.fits 2> /dev/null"), exit()] if args.clean else None
	iter_no = int(args.iter) if args.iter > 0 else default_iter_no
	test = 1 if args.test else 0
	listf = 1 if args.list else 0
	keyword = args.key if args.key else default_keyword
	filename = str(args.file)
##	print keyword, filename, test, listf, iter_no

def files(fname):
	global file_in, file_out, mask_out
	file_in = str(fname)
	print file_in
	file_out, mask_out = ('test.fits', 'mask.fits') if test > 0 else (file_in.replace('.fits', 'c.fits'), file_in.replace('.fits', '-mask.fits'))

def get_info(fname):
	global mode, sat_level, gain_level, ron_level
	mode = pyfits.getval(fname, keyword)
	sat_level, gain_level, ron_level = (45000, 6.7, 0.9) if (mode == 16) else (31000, 9.0, 1.8)
	print "\nFile:", str(fname), "\nSaturation:", sat_level, "\nGain:", gain_level, "\nRON:", ron_level

def remove_cosmics(fname):
	(array, header) = cosmics.fromfits(fname, verbose = False)
	c = cosmics.cosmicsimage(array, pssl = 0.0, gain=gain_level, readnoise=ron_level,
	sigclip = 6.0, sigfrac = 0.3, objlim = 5.5,
	satlevel = sat_level, verbose = False)
	c.run(maxiter = iter_no, verbose = False)
	cosmics.tofits(file_out, c.cleanarray, header)
	cosmics.tofits(mask_out, c.mask, header)

arguments()
if listf == 1:
	with open(filename) as data:
		for line in data:
			line = line.strip()
			files(line)
			get_info(line)
			remove_cosmics(line)
else:
	files(filename)
	get_info(filename)
	remove_cosmics(filename)

exit()
