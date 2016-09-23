#! /usr/bin/env python
# -*- coding: utf-8 -*-

#  commentary
## test module

import sys, os, time, io
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib import rc

start_time = time.time()
global program_path, config_path, version
program_path = str(os.path.realpath(__file__))
config_path = program_path.split(".")[0]+".conf"
version = '{:s}'.format(time.ctime(os.path.getmtime(program_path)))
config_status = os.path.isfile(config_path)
## print program_path, "\n", config_path, config_status, "\n", version
par, files, d, real_filters = [], [], {}, []
file_av_mag, file_init_ext, file_final_ext, file_par = 'av_mag.dat', 'ext_', 'final_ext_', 'fit_'

def arguments():
	import argparse
	global filename, edit, object_name, update, delete
	parser = argparse.ArgumentParser(description="               Atmospheric Extinction for Linux (atm_ext.py)                 \n"
																"               ---------------------------------------------               \n\n"
																"Program improves the light curves to the influence of atmospheric extinction,\n"
																"plots important data, performs essential statistics and is fully interactive.\n\n"
																"Requirements: Python 2.7+ (with modules: argparse, io, matplotlib, numpy, os,\n"
																"              scipy, sys, re, Tkinter, tkMessageBox), 'gvim'                 \n"
																"\n\nVersion: "+version, formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-f", "--file", type=str, help="list of light curves' files to process", required=True)
	parser.add_argument("-o", "--obj", type=str, help="object name (f.e.: M31)", required=True)
	parser.add_argument("-u", "--update", help="update light curves for influence of extinction", action="store_true")
	parser.add_argument("-e", "--edit", help="open configuration file for editing (with 'gvim')", action="store_true")
	parser.add_argument("-d", "--delete", help="delete configuration file", action="store_true")
	parser.add_argument("-c", "--clean", help="deletes all created files", action="store_true")
	args = parser.parse_args()
	filename = str(args.file)
	object_name = str(args.obj)
	edit = 1 if args.edit else 0
	update = 1 if args.update else 0
	delete = True if args.delete else False
	[os.system("rm extinction*.dat *ext*.png av_mag.dat filter*.lst *fit* diff*obs *noext* 2> /dev/null"), exit("Cleaned.")] if args.clean else None

if len(sys.argv) < 2:
		import os
		os.system(program_path+" -h")
		exit()

def popup(text):
	import Tkinter as tk
	import tkMessageBox
	root = tk.Tk()
	root.withdraw()
	tkMessageBox.showinfo("Information", text)

def configuration():
	[os.system("rm "+config_path+" 2> /dev/null"), exit("Configuration file deleted.")] if delete == True else None
##	print config_status, edit
	if config_status == True:
		print "Configuration file status  : OK!"
		[popup("Now edit configuration file. Remember to save changes."), os.system("gvim -f "+config_path), popup("Now relaunch program."), exit("Now relaunch program.")] if edit == 1 else None
		with open(config_path) as f:
			for line in f:
				line = line.split(":")[1].strip("\n")
				par.append(line.strip())
		f.close()
		global col_jd, col_air, col_phot, col_phot, col_err, col_oth, iter, filters
		col_jd, col_air, col_phot, col_err, col_oth, iter, filters = int(par[0]) - 1, int(par[1]) - 1, int(par[2]) - 1, int(par[3]) - 1, par[4].split(" "), int(par[5]), par[6].split(" ")
		col_oth = [ int(x) - 1 for x in col_oth ]
##		print col_jd, col_air, col_phot, col_err, col_oth, iter, filters, exit()
	else:
		"Configuration file status: NONE!"
		popup("It looks like you are running program for the first time. Please fill in the configuration file.")
		os.system("touch "+config_path)
		with open(config_path,"w") as f:
			f.write("Number  of column  with julian date             : \n")
			f.write("Number  of column  with air mass                : \n")
			f.write("Number  of column  with photometry              : \n")
			f.write("Number  of column  with photometry error        : \n")
			f.write("Numbers of columns with other photometry        : \n")
			f.write("Number of iterations for statistics             : \n")
			f.write("Filter system (i.e. B V R I Haw Han)            : \n")
		os.system("gvim -f "+config_path)
		popup("Now relaunch program.")
		exit("Now relaunch program.")

def display_input_parameters(show=True):
	if show != False:
		print "Filter system              :", 
		for p in filters: print p,
		print "\nIterations for statistics  :", iter, "\nObject name                :", object_name

def files_by_band(file):
	with open(file) as data:
		for line in data:
			files.append(line.strip("\n"))
		for filter in filters:
			matching = [s for s in files if filter in s]
			newfile = "filter_"+filter+".lst"
			if len(matching) != 0:
				real_filters.append(filter)
				with io.FileIO(newfile, "w") as f:
					for item in matching:
						f.write(item+"\n") 
				f.close()
	data.close()

def load_data():
	import re
	global stars
	stars = []
	for filter in filters:
		with open("filter_"+filter+".lst", "r") as f:
			for line in f:
				data_jd, data_phot, data_err, data_air = [], [], [], []
				star = int(re.sub("\D", "", line.strip("\n")))
##				print filter, star
##				time.sleep(1)
				stars.append(star)
				d['jd_'+str(star)+'_'+filter], d['air_'+str(star)+'_'+filter], d['phot_'+str(star)+'_'+filter],	d['err_'+str(star)+'_'+filter] = [], [], [], []
				with open(line.strip("\n"), "r") as g:
					for line in g:
						if line.split()[col_phot] != '-9.9999':
							data_jd.append(float(line.split()[col_jd])), data_air.append(float(line.split()[col_air]))
							data_phot.append(float(line.split()[col_phot])), data_err.append(float(line.split()[col_err]))
				d['jd_'+str(star)+'_'+filter], d['air_'+str(star)+'_'+filter], d['phot_'+str(star)+'_'+filter],	d['err_'+str(star)+'_'+filter] = data_jd, data_air, data_phot, data_err
	stars = list(set(stars))
	stars.sort()

def av_mag(stats=False):
	for filter in filters:
		for star in stars:
				d['av_mag_'+str(star)+'_'+filter] = np.average(d['phot_'+str(star)+'_'+filter])
				d['av_err_'+str(star)+'_'+filter] = np.average(d['err_'+str(star)+'_'+filter])
	if stats == False:
		with io.FileIO(file_av_mag, "w") as f:
			f.write('  {:s}  {:4s}   {:6s}   {:6s}  {:6s}\n'.format('#','BAND','AV_MAG','AV_ERR','POINTS'))
			for star in stars:
				for filter in filters:
					f.write('{:3d}   {:2s}   {:.4f}   {:.4f}   {:3d}\n'.format(star, filter, d['av_mag_'+str(star)+'_'+filter], d['av_err_'+str(star)+'_'+filter], len(d['phot_'+str(star)+'_'+filter])))

def av_mag_offset():
	for filter in filters:
		for star in stars:
			d['diff_'+str(star)+'_'+filter] = []
			for item in d['phot_'+str(star)+'_'+filter]:
				diff = item - d['av_mag_'+str(star)+'_'+filter]
				d['diff_'+str(star)+'_'+filter].append(diff)

def plotting(x,y,filter,min_x,max_x,min_y,max_y):
	plt.clf()
	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	rc('text', usetex=True)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig = plt.figure()
	fig.suptitle('Atmospheric extinction for '+object_name, fontsize=14)
	ax = fig.add_subplot(111)
	fig.subplots_adjust(top=0.90)
	ax.set_title('Band '+filter, fontsize=12, fontweight='bold')
	ax.set_xlabel(r'\textit{air mass}',fontsize=12)
	ax.set_ylabel(r'\textit{m - \overline{m}}',fontsize=12)
##	ax.grid(True)
	fit = np.polyfit(x,y,1)
	global slope, intercept, r_value, p_value, std_err
	slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	std_err = np.std(y)
	d['stats_'+filter] = []
	d['stats_'+filter].append(slope)
	d['stats_'+filter].append(intercept)
	d['stats_'+filter].append(std_err)
	d['stats_'+filter].append(r_value)
	fit_fn = np.poly1d(fit)
	x_fit = np.linspace(min_x, max_x, 100) 
	y_fit = fit_fn(x_fit)
	plt.xlim(min_x, max_x)
	plt.ylim(min_y, max_y)
	plt.gca().invert_yaxis()
	ax.plot(x, y, 'ko', ms='2.0', label='data')
	ax.plot(x_fit, y_fit, '--', lw=.5, label='fit')
	ax.plot(x_fit, y_fit - 2 * std_err, 'r', lw=.5, label=r'$2\sigma$')
	ax.plot(x_fit, y_fit + 2 * std_err, 'r', lw=.5, label='')
	plt.legend(fontsize=10, loc='lower left', frameon=False)
	sign = "+" if intercept > 0 else	""
	ax.text(0.98 * max_x, 0.6 * min_y, r'$y= {:.3f} x {:1s} {:.3f}$'.format(slope, sign, intercept), fontsize=12)
	ax.text(0.98 * max_x, 0.5 * min_y, r'$R= {:.3f}$'.format(r_value), fontsize=12)
	ax.text(0.98 * max_x, 0.4 * min_y, r'$2\sigma= {:.3f}$'.format(2 * std_err), fontsize=12)
	plt.savefig(file_init_ext+filter+'.png', bbox_inches='tight')

def plotting_final(x,y,filter,min_x,max_x,min_y,max_y):
	plt.clf()
	rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	rc('text', usetex=True)
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	fig = plt.figure()
	fig.suptitle('Atmospheric extinction for '+object_name, fontsize=14)
	ax = fig.add_subplot(111)
	fig.subplots_adjust(top=0.90)
	ax.set_title('Band '+filter+' (filtered)', fontsize=12, fontweight='bold')
	ax.set_xlabel(r'\textit{air mass}',fontsize=12)
	ax.set_ylabel(r'\textit{m - \overline{m}}',fontsize=12)
	plt.xlim(min_x, max_x)
	plt.ylim(min_y, max_y)
	plt.gca().invert_yaxis()
	ax.plot(x, y, 'ko', ms='2.0', label='data')
	x_fit = np.linspace(min_x, max_x, 100)
	slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
	y_fit = slope * x_fit + intercept
	ax.plot(x_fit, y_fit, '--', lw=.5, label='fit')
	ax.plot(x_fit, y_fit - 2 * d['stats_'+filter][2], 'r', lw=.5, label=r'$2\sigma$')
	ax.plot(x_fit, y_fit + 2 * d['stats_'+filter][2], 'r', lw=.5, label='')
	plt.legend(fontsize=10, loc='lower left', frameon=False)
	sign = "+" if d['stats_'+filter][1] > 0 else	""
	ax.text(0.98 * max_x, 0.6 * min_y, r'$y= {:.3f} x {:1s} {:.3f}$'.format(slope, sign, intercept), fontsize=12)
	ax.text(0.98 * max_x, 0.5 * min_y, r'$R= {:.3f}$'.format(r_value), fontsize=12)
	ax.text(0.98 * max_x, 0.4 * min_y, r'$2\sigma= {:.3f}$'.format(2 * d['stats_'+filter][2]), fontsize=12)
	plt.savefig(file_final_ext+filter+'.png', bbox_inches='tight')
	d['stats_'+filter][0], d['stats_'+filter][1], d['stats_'+filter][3] = slope, intercept, r_value

def plot_ext(stats=False):
	xr, yr = [], []
	for filter in filters:
		for star in stars:
			for item in d['air_'+str(star)+'_'+filter]:
				xr.append(item)
			for item in d['diff_'+str(star)+'_'+filter]:
				yr.append(item) 
	for filter in filters:
		x_data, y_data = [], []
		for star in stars:
			for item in d['air_'+str(star)+'_'+filter]:
				x_data.append(item)
			for item in d['diff_'+str(star)+'_'+filter]:
				y_data.append(item)
		plotting(x_data, y_data, filter, min(xr)-0.01, max(xr)+0.01, min(yr)-0.01, max(yr)+0.01) if stats == False else plotting_final(x_data, y_data, filter, min(xr)-0.01, max(xr)+0.01, min(yr)-0.01, max(yr)+0.01)

def statistics():
	i = 0
	for filter in filters:
		slope, intercept, std_err, r_value = d['stats_'+filter][0], d['stats_'+filter][1], d['stats_'+filter][2], d['stats_'+filter][3]
		for star in stars:
			d['resid_'+str(star)+'_'+filter] = []
			for item in d['phot_'+str(star)+'_'+filter]:
				resid = d['diff_'+str(star)+'_'+filter][d['phot_'+str(star)+'_'+filter].index(item)] - (slope * d['air_'+str(star)+'_'+filter][d['phot_'+str(star)+'_'+filter].index(item)] + intercept)
				d['resid_'+str(star)+'_'+filter].append(resid)
		x, y = [], []
		for star in stars:
			for item in d['air_'+str(star)+'_'+filter]:
				x.append(item)
			for item in d['resid_'+str(star)+'_'+filter]:
				y.append(item)
		x, y = np.array(x), np.array(y)
		slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
		std_err = np.std(y)
		d['stats_'+filter][0], d['stats_'+filter][1], d['stats_'+filter][2], d['stats_'+filter][3] = slope, intercept, std_err, r_value
		for item in d['resid_'+str(star)+'_'+filter]:
			if abs(item) > (2 * std_err):
				d['air_'+str(star)+'_'+filter].pop(d['resid_'+str(star)+'_'+filter].index(item))
				d['phot_'+str(star)+'_'+filter].pop(d['resid_'+str(star)+'_'+filter].index(item))
				d['err_'+str(star)+'_'+filter].pop(d['resid_'+str(star)+'_'+filter].index(item))
				d['diff_'+str(star)+'_'+filter].pop(d['resid_'+str(star)+'_'+filter].index(item))
				d['resid_'+str(star)+'_'+filter].pop(d['resid_'+str(star)+'_'+filter].index(item))
				i += 1
	av_mag(True), av_mag_offset()
	print i, "data points removed!"

def save_ext():
	for filter in filters:
		newfile = file_par+filter+'.par'
		with io.FileIO(newfile, "w") as f:
			for item in d['stats_'+filter]:
				f.write('{:.4f}  '.format(item))

def load_ext():
	global bands
	bands = []
	for filter in filters:
		loadfile = file_par+filter+'.par'
		if os.path.isfile(loadfile) == True:
			slope, intercept, std_err, r_value = np.loadtxt(loadfile)		
			bands.append(filter)
			d['stats_'+filter] = []
			d['stats_'+filter].append(slope)
			d['stats_'+filter].append(intercept)
			d['stats_'+filter].append(std_err)
			d['stats_'+filter].append(r_value)
		else:
			bands.append('No '+filter+'!')
 	for item in bands: print item,
	print ""
	
def check_col_no(line):
	with open(line) as g:
		line = g.readline()
		line = line.split()
		return len(line)

def update_lc():	
	global data_points
	data_points = {}
	with open(filename) as f:
		line = f.readline()
		line = line.strip("\n")
	col_no = check_col_no(line)
	for filter in filters:
		slope = d['stats_'+filter][0]
		with open('filter_'+filter+'.lst', 'r') as data:
			for line in data:
				line = line.strip("\n")
				out = line.split(".")[0]+"_noext.obs"
				with open(line, 'r') as f:
					data_points['data_'+filter] = np.loadtxt(line, dtype=str)
					for i in range(len(data_points['data_'+filter])):
						air = float(data_points['data_'+filter][i][col_air])
						print line
						phot = float(data_points['data_'+filter][i][col_phot])
						if phot != -9.9999:
							phot_noext = phot - slope * (air - 1.0)
							data_points['data_'+filter][i][col_phot] = '{:.4f}'.format(phot_noext)
						for j in range(len(col_oth)):
							oth = float(data_points['data_'+filter][i][col_oth[j]])
							if oth != -9.9999:
								oth_noext = oth - slope * (air - 1.0)
								data_points['data_'+filter][i][col_oth[j]] = '{:.4f}'.format(oth_noext)
				with io.FileIO(out, 'w') as g:
					for item in data_points['data_'+filter]:
						for element in item:
							g.write('{:s}  '.format(element))
						g.write("\n")

arguments()
configuration()
display_input_parameters()
files_by_band(filename)
filters = real_filters
if update == 0:
	print '{:40s}'.format('> Loading data...'),
	load_data()
	print '{:5s}'.format('DONE!')
	print '{:40s}'.format('> Average magnitudes...'),
	av_mag()
	av_mag_offset()
	print '{:5s}' ' ({:10s})'.format('DONE!', file_av_mag)
	print '{:40s}'.format('> Plotting initial extinction...'),
	plot_ext(False)
	print '{:5s}' ' ({:9s})'.format('DONE!', file_init_ext+'*.png')
	for i in range(iter):
		print '{:24s}({:1d}){:13s}'.format('> Performing statistics ', i + 1, '...'),
		statistics()
	print '{:40s}'.format('> Plotting final extinction...'),
	plot_ext(True)
	print '{:5s}' ' ({:9s})'.format('DONE!', file_final_ext+'*.png')
	print '{:40s}'.format('> Saving fitting data...'),
	save_ext()
	print '{:5s}' ' ({:9s})'.format('DONE!', file_par+'*.par')
	print '{:13s} {:.2f} {:8s}'.format('$ ALL DONE in', time.time() - start_time, 'seconds!')
else:
	print '{:28s}'.format('> Loading fitting data...'),
	load_ext()
	print '{:28s}'.format('> Updating light curves...'),
	update_lc()
	print '{:5s}' ' ({:10s})'.format('DONE!', '*noext.obs')
	print '{:13s} {:.2f} {:8s}'.format('$ ALL DONE in', time.time() - start_time, 'seconds!')
exit()
