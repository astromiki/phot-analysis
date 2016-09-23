#! /usr/bin/env python
# -*- coding: utf-8 -*-

#  commentary
## test module

import sys, os, time, io

start_time = time.time()
default_sigma, default_patt, files = 3.5, 'CALIB', []
program_path = str(os.path.realpath(__file__))
config_path = program_path.split(".")[0]+".conf"
config_status = os.path.isfile(config_path)
version = '{:s}'.format(time.ctime(os.path.getmtime(program_path)))
dependencies_list = ['gvim', 'dao2ds9', 'ds92dao', 'mfwhm.bash', 'sdb_xy', 'ap2lst', 'join_list_reg_xy.sh', 'rm_cosmics.py']

class bcolors:
	HEADER = '\033[95m'
	OKBLUE = '\033[94m$ '
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'
	BOLD = '\033[1m'
	UNDERLINE = '\033[4m'
	def disable(self):
		self.HEADER, self.OKBLUE, self.OKGREEN, self.WARNING = '', '', '', '', 
		self.FAIL, self.ENDC, self.BOLD, self.UNDERLINE = '', '', '', ''

def clean(all=False):
	os.system("rm *.is.* %1 *.psf *.ap* *.nei *.err *.fwhm *.par *.opt *.als *.tar.gz *.coo *.lst *.*_stars .dir* *.reg *.srt .t *-sub.fits 2> /dev/null")
	os.system("rm *.fit* 2> /dev/null") if all == True else None
	print bcolors.OKBLUE + "Cleaned." + bcolors.ENDC if all == False else bcolors.OKBLUE + "Cleaned (including FITS)." + bcolors.ENDC
	exit()

def popup(text):
	import Tkinter as tk
	import tkMessageBox
	print bcolors.HEADER + "> " + text + bcolors.ENDC
	root = tk.Tk()
	root.withdraw()
	tkMessageBox.showinfo("Information", text.strip("."))

def args():
	import argparse
	import textwrap
	global object, sigma, resume, tree, patt, edit, delete
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=textwrap.dedent('''\
		          'Make list of stars' for Linux (mk_lst.py)	
		-----------------------------------------------------------------
		Program prepares lists of stars for further photometric analysis:
		a list of all stars, list for PSF fitting & list constant stars.

		Requirements: Python 2.7+ (with modules: argparse, io, matplotlib,
		              numpy, os, scipy, sys, re, Tkinter, tkMessageBox),
		              gvim, SAODS9, DAOPHOT2 (+ALLSTAR), dao2ds9, ds92dao,
		              join_list_reg_xy.sh, ap2lst, rm_cosmics.py (all PM),	
		              mfwhm.bash (ZK: kolaczkowski(at)astro.uni.wroc.pl), 
		              sdb_xy (PB: brus(at)astro.uni.wroc.pl)
		'''),
		epilog="Version: "+version)
	parser.add_argument("-o", "--object", type=str, metavar='OBJ', help="object name (f.e.: M31)", required=True)
	parser.add_argument("-s", "--sigma", type=float, metavar='SIG', help="sigma threshold for DAOPHOT (f.e.: 3.5)")
	parser.add_argument("-r", "--resume", help="resume working with list of stars", action="store_true")
	parser.add_argument("-t", "--tree", help="look for images in many directories", action="store_true")
	parser.add_argument("-p", "--patt", type=str, metavar='...', help="pattern for searching directories ")
	parser.add_argument("-e", "--edit", help="open configuration file for editing (with 'gvim')", action="store_true")
	parser.add_argument("-d", "--delete", help="delete configuration file", action="store_true")
	parser.add_argument("-c", "--clean", help="deletes all created files", action="store_true")
	parser.add_argument("-cc", "--cleanall", help="deletes all created files (including FITS images)", action="store_true")
	args = parser.parse_args()
	object = str(args.object)
	edit = True if args.edit else False
	delete = True if args.delete else False
	clean() if args.clean else None
	clean(True) if args.cleanall else None
	sigma = float(args.sigma) if args.sigma else default_sigma
	resume = True if args.resume else False
	tree = True if args.tree else False
	patt = str(args.patt) if args.patt else default_patt

def dependencies():
	import subprocess
	missing, missing_c = [], 0
	for prog in dependencies_list:
		try:
			subprocess.check_output('which ' + prog ,shell=True)
		except: 
			missing_c += 1
			missing.append(prog)
	print "$ Dependencies status       :",
	print bcolors.OKGREEN + "OK!" + bcolors.ENDC if missing_c == 0 else bcolors.FAIL + "MISSING! (" + str(missing_c) + ")" + bcolors.ENDC
	if missing_c > 0:
		print "# You must install          :",
		for item in missing:
			print "'" + item + "'",
		exit()

def configuration():
	if delete == True:
		os.system("rm "+config_path+" 2> /dev/null")
		print bcolors.OKBLUE + "Configuration file deleted." + bcolors.ENDC
		exit()
	else:
		None
	if config_status == True:
		print "$ Configuration file status : " +  bcolors.OKGREEN + "OK!" + bcolors.ENDC
		[popup("Now edit configuration file. Remember to save changes."), os.system("gvim -f "+config_path), popup("Now relaunch program."), exit()] if edit == True else None
		par = []
		with open(config_path) as f:
			for line in f:
				line = line.split(":")[1].strip("\n")
				par.append(line.strip())
		f.close()
		global path_calib, path_phot, path_objects, offset, daophot, allstar, default_sigma, default_patt, keyword
		path_calib, path_phot, path_objects, offset = par[0], par[1], par[2], int(par[3])
		daophot, allstar, default_sigma, default_patt, keyword = par[4], par[5], float(par[6]), par[7], par[8]
##		print path_calib, path_phot, path_objects, offset, daophot, allstar, default_sigma, default_patt, exit()
	else:
		print "$ Configuration file status: " +  bcolors.FAIL + "NONE!" + bcolors.ENDC
		popup("It looks like you are running program for the first time. Please fill in the configuration file.")
		os.system("touch "+config_path)
		with open(config_path,"w") as f:
			f.write("Path to calib-bialkow                                      : \n")
			f.write("Path to phot-bialkow                                       : \n")
			f.write("Path to objects                                            : \n")
			f.write("Offset from image edge (in pixels)                         : 2        \n")
			f.write("DAOPHOT running command                                    : sdaophot \n")
			f.write("ALLSTAR running command                                    : sallstar \n")
			f.write("Default sigma threshold for DAOPHOT                        : 3.5      \n")
			f.write("Default pattern for searching images in directory tree     : CALIB    \n")
			f.write("FITS keyword for readtime                                  : GAIN     \n")
		os.system("gvim -f "+config_path)
		popup("Now relaunch program.")
		exit()

def get_info(file):
	global gain, xsize, ysize
	import pyfits
	hdulist = pyfits.open(file)
##	hdulist.verify('fix')
	gain = hdulist[0].header[keyword]
	xsize = hdulist[0].header["NAXIS1"]
	ysize = hdulist[0].header["NAXIS2"]
	hdulist.close()

def ask(prompt, reminder='# Please try again!'):
	while True:
		ok = raw_input("> " + prompt + " [Y/n] ",)
		if ok in ('y', 'ye', 'yes', 'Y', 'Yes', 'Ye', ''):
			return False
		if ok in ('n', 'no', 'nop', 'nope', 'N', 'No'):
			print bcolors.FAIL + "$ Program stopped." + bcolors.ENDC
			clean()
		print bcolors.WARNING + reminder + bcolors.ENDC

def choose_psf_stars(reminder='# Please try again!'):
	global psf_stars
	while True:
		psf_stars = raw_input("> Number of stars for PSF: ",)
		if psf_stars.isdigit() == True:
			return False
		print bcolors.WARNING + reminder + bcolors.ENDC

def run_ds9(image,region=None):
	os.system("ds9 " + image + " -zoom to fit -cmap Heat -linear -scale mode zscale") if region == None else os.system("ds9 " + image + " -regions " + region + " -zoom to fit -cmap Heat -linear -scale mode zscale")

def dao2ds9(infile,outfile,back=False,iters=False):
	if back == False:
		os.system("dao2ds9 " + infile + " > .t && rm .t")
		os.system("mv " + FILE_REG + " " + outfile)
	if back == True:
		os.system("ds92dao " + outfile + " .temp > .t && rm .t")
		os.system("head -3 " + FILE_APS + " > .t")
		os.system("sdb_xy " + FILE_APS + " .temp -r 3.5 -h 3 3 -m 2 >> .t") if iters == False else os.system("sdb_xy " + FILE_ALS + " .temp -r 3.5 -h 3 3 -m 2 >> .t")
		os.system("mv .t " + FILE_LST + " && rm .temp")
	

def ap2lst(infile,outfile):
	os.system("ap2lst " + infile + " > " + outfile)	

def run_daophot(action):
	f = open('daophot.par', 'w')
	print ">",
	if action == 'find':
		print "Finding stars in reference image...",
		command = "op\n\nTH=" + str(sigma) + "\n\nat " + FILE_FITS + "\nfi\n1,1\n" + FILE_COO + "\ny\nex"
	if action == 'phot':
		print "Performing aperture photometry...",
		command = "at " + FILE_FITS + "\nph\n\nph\n\n\nex"
	if action == 'pick-psf':
		print "Picking stars for PSF...",
		command = "at " + FILE_FITS + "\npi\n" + FILE_APS + "\n" + psf_stars + ",20\n" + FILE_LST + "\nex"
	if action == 'psf':
		print "Performing PSF fitting...",
		command = "at " + FILE_FITS + "\npsf\n" + FILE_APS + "\n" + FILE_LST + "\n" + FILE_PSF + "\nex"
	if action == 'sort':
		print "Sorting stars by magnitude...",
		command = "so\n+4\n" + FILE_APX + "\n" + FILE_APS + "\ny\nex"
	if action == 'sort2':
		print "Sorting stars by magnitude...",
		command = "so\n+4\n" + FILE_ALS + "\n\ny\nex"
	f.write(command)
	f.close()
	sys.stdout.flush()
	os.system(daophot + ' < daophot.par > .t && rm .t')
	print bcolors.OKGREEN + "DONE!" + bcolors.ENDC
	sys.stdout.flush()
	drop('daophot.par')

def run_allstar(action):
	f = open('allstar.par', 'w')
	print ">",
	if action == 'pphot':
		print "Performing profile photometry...",
		command = "at\n" + FILE_FITS + "\n" + FILE_PSF + "\n" + FILE_APS + "\n" + FILE_ALS + "\n" + FILE_SUB + "\nex"
	f.write(command)
	f.close()
	os.system(allstar + ' < allstar.par > .t && rm .t')
	print bcolors.OKGREEN + "DONE!" + bcolors.ENDC
	drop('allstar.par')

def compare_lists():
	os.system("mv " + FILE_SRT + " " + FILE_ALS)
	os.system("head -3 " + FILE_LST + " > .temp")
	os.system("sdb_xy " + FILE_ALS + " " + FILE_LST + " -r 3.5 -h 3 3 -m 2 >> .temp")
	os.system("mv .temp " + FILE_LST)

def drop(what):
	os.system("rm " + what + " 2> /dev/null")

def input_data(rtime,path,const=False):
	os.system("tar xf "+path+"input/daophot_input.tar.gz && rm *.par")
	os.system("mv photo-Bialkow.opt photo.opt")
	os.system("mv allstar-Bialkow.opt allstar.opt")
	if const == True:
		os.system("mv daophot-Bialkow-"+str(rtime)+"_const.opt daophot.opt")
		print "$ New input data gathered."
	else:
		os.system("mv daophot-Bialkow-"+str(rtime)+"_var.opt daophot.opt")
		print "$ Input data gathered."
	os.system("rm *16*.opt *2*.opt 2> /dev/null")

def list_files(obj,indir=True):
	if indir == True:
		os.system("ls "+obj+"-*.fits > .t")
	else:
		None
	NOF = 0
	for item in open('.t', 'r'):
		files.append(item.strip('\n'))
		NOF += 1
	print '$ Files found: '+str(NOF)
	drop('.t')

def set_filenames(obj):
	global FILE_FITS, FILE_AP, FILE_APX, FILE_APS, FILE_SRT, FILE_LST, FILE_COO, FILE_SUB, FILE_CONSTSTARS
	global FILE_PSF, FILE_NEI, FILE_ALS, FILE_REG, FILE_ALL_REG, FILE_PSF_REG, FILE_ALLSTARS, FILE_PSFSTARS
	FILE_FITS = obj + ".fits"
	FILE_AP = obj + ".ap"
	FILE_APX = obj + ".apx"
	FILE_APS = obj + ".aps"
	FILE_SRT = obj + ".srt"
	FILE_LST = obj + ".lst"
	FILE_COO = obj + ".coo"
	FILE_SUB = obj + "-sub.fits"
	FILE_PSF = obj + ".psf"
	FILE_NEI = obj + ".nei"
	FILE_ALS = obj + ".als"
	FILE_REG = obj + ".reg"
	FILE_ALL_REG = obj + "-all.reg"
	FILE_PSF_REG = obj + "-psf.reg"
	FILE_ALLSTARS = obj + ".all_stars"
	FILE_PSFSTARS = obj + ".psf_stars"
	FILE_CONSTSTARS = obj + ".const_stars"

def fwhm(file):
	os.system("mfwhm.bash "+file+" >> .t")

def ref_img():
	results, res = [], {}
	res['name'], res['stars'], res['fwhm'] = [], [], []
	for item in open('.t', 'r'):
		results.append(item.split())
	for i in range(len(results)):
		res['name'].append(results[i][0])
		res['stars'].append(int(results[i][4]))
		res['fwhm'].append(float(results[i][7]))
	best_s = max(res['stars'])
	best_f = res['fwhm'][res['stars'].index(best_s)]
	best_n = res['name'][res['stars'].index(best_s)]
	print "$ Best file is '" + bcolors.BOLD + best_n + bcolors.ENDC + "' with " + bcolors.BOLD + str(best_s) + bcolors.ENDC + " stars found and FWHM = " + bcolors.BOLD + str(best_f) + bcolors.ENDC + "."
	os.system("cp " + best_n + " " + FILE_FITS)
	print bcolors.OKGREEN + bcolors.BOLD + "$ Reference image '" + FILE_FITS + "' has been created!" + bcolors.ENDC
	drop('.t')

def clean_edge(file):
	x_offset_l = str(offset)
	x_offset_r = str(xsize - offset)
	y_offset_d = str(offset)
	y_offset_u = str(ysize - offset)
	print "> Removing stars too close to edges (" + str(offset) + " pixels margin)...",
	os.system("head -3 " + file + " > .temp")
	os.system("cat "+ file + " | awk -v a=" + x_offset_l + " -v b=" + x_offset_r + " -v c=" + y_offset_d + " -v d=" + y_offset_u + " 'NR > 3 && $2 > a && $2 < b && $3 > c && $3 < d {print $0}' >> .temp")
	os.system("mv .temp " + file)
	print bcolors.OKGREEN + "DONE!" + bcolors.ENDC

def check_var_psf():
	global psf_status
	psf_status = True
	with open(FILE_PSF, 'r') as f:
		line = f.readline()
		col = len(line.split())
	if col == 0:
		print bcolors.WARNING + "# Variable PSF is impossible! Trying constant PSF." + bcolors.ENDC
		input_data(gain,path_phot,True)
		os.system("rm i.err "+ FILE_PSF)
		psf_status = False
	else:
		psf_status = True

def join_lists(lst,reg):
	os.system("join_list_reg_xy.sh " + lst + " " + reg + " .temp > .t && rm .t")
	os.system("mv .temp " + lst)

def prepare_for_next_iter(iter_c):
	os.system("cp " + FILE_ALS + " " + FILE_APS)
	os.system("mv " + FILE_SUB + " iter-" + str(iter_c) + "-" + FILE_SUB)
	os.system("mv " + FILE_ALS + " iter-" + str(iter_c) + "-" + FILE_ALS)
	os.system("rm " + FILE_PSF + " " + FILE_NEI + " i.err *.reg")

def prepare_output():
	print "> Preparing your lists of stars (& DS9 regions)...",
	dao2ds9(FILE_ALS, FILE_ALL_REG)
	os.system("mv " + FILE_ALS + " " + FILE_ALLSTARS)
	dao2ds9(FILE_LST, FILE_PSF_REG)
	os.system("mv " + FILE_LST + " " + FILE_PSFSTARS)
	if const_stars_switch == True:
		with open(FILE_CONSTSTARS, 'w') as f:
			for i in constant_stars:
				f.write(i + "\n")
		f.close()
	os.system("rm *ap* *.par *.nei *.coo *.opt *.fwhm *.err %1 2> /dev/null")
	print bcolors.OKGREEN + "DONE!" + bcolors.ENDC

def ask_const_stars():
	global constant_stars, const_stars_switch
	constant_stars = []
	while True:
		ok = raw_input("> Do you wish to enter constant stars? [Y/n] ",)
		if ok in ('y', 'ye', 'yes', 'Y', 'Yes', 'Ye', ''):
			const_stars_switch = True
			s_no = raw_input("> Please enter numbers: ",)
			s_no = s_no.split(" ")
			constant_stars = [s for s in s_no if s.isdigit()]
			break
		if ok in ('n', 'no', 'nop', 'nope', 'N', 'No'):
			const_stars_switch = False
			print "$ No constant stars added."
			break
		print bcolors.WARNING + "# Please try again!" + bcolors.ENDC

def copy_output():
	print "> Copying your output to " + path_objects + "...",
	os.system("tar czvf " + object +".tar.gz " + FILE_FITS + " " + FILE_SUB + " " + FILE_ALLSTARS + " " + FILE_ALL_REG + " " + FILE_PSFSTARS + " " + FILE_PSF_REG + " > .t && rm .t") if const_stars_switch == False else os.system("tar czvf " + object +".tar.gz " + FILE_FITS + " " + FILE_SUB + " " + FILE_ALLSTARS + " " + FILE_ALL_REG + " " + FILE_PSFSTARS + " " + FILE_PSF_REG + " " + FILE_CONSTSTARS + " > .t && rm .t")
	os.system("cp " + object + ".tar.gz " + path_objects)
	print bcolors.OKGREEN + "DONE!" + bcolors.ENDC

def prepare_for_resume(name):
	print "> Resuming your previous work with lists of stars for " + bcolors.BOLD + "'" + name + "'" + bcolors.ENDC + "...",
	if os.path.isfile(path_objects + name + '.tar.gz') == True:
		os.system("cp " + path_objects + name + ".tar.gz .")
		os.system("tar xf " + path_objects+name + ".tar.gz")
		os.system("mv " + FILE_ALLSTARS + " " + FILE_APS)
		os.system("mv " + FILE_PSFSTARS + " " + FILE_LST)
		print bcolors.OKGREEN + "DONE!" + bcolors.ENDC
	else:
		print bcolors.FAIL + "FAILED!" + bcolors.ENDC
		print bcolors.WARNING + "# No suitable data exists!" + bcolors.ENDC
		print "$ Correct & relaunch."
		exit()

def remove_cosmics(image):
	print "> Removing cosmic rays...",
	sys.stdout.flush()
	os.system("rm_cosmics.py -f " + image + " -i 4 -t > .t && rm .t mask.fits")
	os.system("mv test.fits " + image)
	print bcolors.OKGREEN + "DONE!" + bcolors.ENDC

args()
dependencies()
configuration()

print "\n==================================\n" + bcolors.BOLD + "       STEP I: PREPARATION        " + bcolors.ENDC + "\n==================================\n"
set_filenames(object)
if resume == False:
	list_files(object)
	print "> Measuring FWHMs...",
	sys.stdout.flush()
	[fwhm(file) for file in files]
	print bcolors.OKGREEN + "DONE!" + bcolors.ENDC
	ref_img()
	popup("Please inspect reference image.")
	run_ds9(FILE_FITS)
	ask("Do you accept?")
	get_info(FILE_FITS)
	input_data(gain,path_phot)
	remove_cosmics(FILE_FITS)
	run_daophot('find')
	clean_edge(FILE_COO)
	stars_count = sum(1 for line in open(FILE_COO))
	print bcolors.OKGREEN + "$ " + bcolors.BOLD + str(stars_count) + bcolors.ENDC + bcolors.OKGREEN + " stars found." + bcolors.ENDC
	run_daophot('phot')
	ap2lst(FILE_AP, FILE_APX)
	run_daophot('sort')
	choose_psf_stars()
	run_daophot('pick-psf')
	popup("Please check & correct your PSF stars list...")
	dao2ds9(FILE_LST, FILE_PSF_REG)
	run_ds9(FILE_FITS, FILE_PSF_REG)
	dao2ds9(FILE_LST, FILE_PSF_REG, back=True)
	print bcolors.BOLD + "$ Preeliminary lists of stars have been created!" + bcolors.ENDC
else:
	prepare_for_resume(object)
	get_info(FILE_FITS)
	input_data(gain,path_phot)
	stars_count = sum(1 for line in open(FILE_APS))
	stars_count_psf = sum(1 for line in open(FILE_LST))
	print bcolors.OKGREEN + "$ " + bcolors.BOLD + str(stars_count) + bcolors.ENDC + bcolors.OKGREEN + " stars found (" + bcolors.BOLD + str(stars_count_psf) + bcolors.ENDC + bcolors.OKGREEN + " stars for PSF)." + bcolors.ENDC
	
print "\n==================================\n" + bcolors.BOLD + "       STEP II: ITERATIONS        " + bcolors.ENDC + "\n==================================\n"
i = 1
while True:
	print "--- " + bcolors.UNDERLINE + "Iteration no. " + str(i) + bcolors.ENDC + ":\n"
	run_daophot('psf')
	check_var_psf()
	run_daophot('psf') if psf_status == False else None
	run_allstar('pphot')
	run_daophot('sort2')
	compare_lists()
	popup("Please inspect subtracted image...")
	run_ds9(FILE_SUB)
	while True:
		ok = raw_input("> Are you happy? [Y/n] ",)
		if ok in ('y', 'ye', 'yes', 'Y', 'Yes', 'Ye', ''):
			print bcolors.BOLD + bcolors.OKGREEN + "$ Very well then!" + bcolors.ENDC
			go = True
			break
		if ok in ('n', 'no', 'nop', 'nope', 'N', 'No'):
			go = False
			break
		print bcolors.WARNING + "# Please try again!" + bcolors.ENDC
	if go == True:
		break
	else:
		popup("Please check & correct your ALL stars list (on raw image)...")
		dao2ds9(FILE_ALS, FILE_ALL_REG)
		run_ds9(FILE_FITS, FILE_ALL_REG)
		popup("Please check & correct your ALL stars list (on sub image)...")
		run_ds9(FILE_SUB, FILE_ALL_REG)
		join_lists(FILE_ALS,FILE_ALL_REG)
		popup("Please check & correct your PSF stars list...")
		dao2ds9(FILE_LST, FILE_PSF_REG)
		run_ds9(FILE_FITS, FILE_PSF_REG)
		dao2ds9(FILE_LST, FILE_PSF_REG, back=True, iters=True)
		prepare_for_next_iter(i)
		i += 1
		print ""
		continue
	
print "\n==================================\n" + bcolors.BOLD + "         STEP III: ENDING         " + bcolors.ENDC + "\n==================================\n"
ask_const_stars()
prepare_output()
copy_output()

print "\n==================================\n" + bcolors.BOLD + "         STEP IV: SUMMARY         " + bcolors.ENDC + "\n==================================\n"
print "$ Object name: '" + bcolors.BOLD + object + bcolors.ENDC + "'"
print "$ All stars: " + str(sum(1 for line in open(FILE_ALLSTARS)))
print "$ Stars for PSF fitting: " + str(sum(1 for line in open(FILE_PSFSTARS)))
print "$ Constant stars entered:",
if const_stars_switch == True:
	print str(sum(1 for line in open(FILE_CONSTSTARS))) 
	print "$ Files created: " + FILE_FITS + ", " + FILE_ALLSTARS + ", " + FILE_PSFSTARS + " & " + FILE_CONSTSTARS
else:
	print str(0)
	print "$ Files created: " + FILE_FITS + ", " + FILE_ALLSTARS + ", " + FILE_PSFSTARS
print "$ Number of iterations needed: " + str(i)
print "$ Time elapsed: " + bcolors.BOLD + '{:.2f}'.format(time.time() - start_time) + bcolors.ENDC + " s."

exit()
