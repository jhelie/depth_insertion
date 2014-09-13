#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.1"
parser = argparse.ArgumentParser(prog = 'depth_insertion', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
**********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/depth_insertion
DOI: 
**********************************************

[ DESCRIPTION ]
 
This script calculate the average distance between each particle of a peptide
and the selected leaflet.


[ REQUIREMENTS ]

The following python modules are needed :
 - MDAnalysis
 - matplotlib
 - np
 - sp
 - networkX


[ NOTES ]

1. The distance is calculated with respect to the z axis, not the bilayer normal. So the
   more your system deforms the noiser and the less meaningful the results get.

 
[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: structure file [.gro] (required)
-x			: trajectory file [.xtc]
-o			: name of output folder
-b			: beginning time (ns) (the bilayer must exist by then!)
-e			: ending time (ns)	
-t 		1	: process every t-frames
--leaflet		: reference leaflet ('upper' or 'lower')

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#data options
parser.add_argument('-f', nargs=1, dest='grofilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-x', nargs=1, dest='xtcfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-o', nargs=1, dest='output_folder', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-b', nargs=1, dest='t_start', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-e', nargs=1, dest='t_end', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-t', nargs=1, dest='frames_dt', default=[1], type=int, help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
#data options
args.grofilename = args.grofilename[0]
args.xtcfilename = args.xtcfilename[0]
args.output_folder = args.output_folder[0]
args.t_start = args.t_start[0]
args.t_end = args.t_end[0]
args.frames_dt = args.frames_dt[0]

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================
#generic science modules
try:
	import math
except:
	print "Error: you need to install the maths module."
	sys.exit(1)
try:
	import numpy as np
except:
	print "Error: you need to install the numpy module."
	sys.exit(1)
try:
	import scipy as sp
	import scipy.stats
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	mcolorconv = mcolors.ColorConverter()
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)

#MDAnalysis module
try:
	import MDAnalysis
	from MDAnalysis import *
	import MDAnalysis.analysis
	import MDAnalysis.analysis.leaflet
	import MDAnalysis.analysis.distances
	#set MDAnalysis to use periodic boundary conditions
	MDAnalysis.core.flags['use_periodic_selections'] = True
	MDAnalysis.core.flags['use_KDTree_routines'] = False
except:
	print "Error: you need to install the MDAnalysis module first. See http://mdanalysis.googlecode.com"
	sys.exit(1)

#=========================================================================================
# sanity check
#=========================================================================================
if not os.path.isfile(args.grofilename):
	print "Error: file " + str(args.grofilename) + " not found."
	sys.exit(1)
if args.t_end != -1 and args.t_end < args.t_start:
	print "Error: the starting time (" + str(args.t_start) + "ns) for analysis is later than the ending time (" + str(args.t_end) + "ns)."
	sys.exit(1)

if args.xtcfilename == "no":
	if '-t' in sys.argv:
		print "Error: -t option specified but no xtc file specified."
		sys.exit(1)
	elif '-b' in sys.argv:
		print "Error: -b option specified but no xtc file specified."
		sys.exit(1)
	elif '-e' in sys.argv:
		print "Error: -e option specified but no xtc file specified."
		sys.exit(1)
elif not os.path.isfile(args.xtcfilename):
	print "Error: file " + str(args.xtcfilename) + " not found."
	sys.exit(1)

#=========================================================================================
# create folders and log file
#=========================================================================================
if args.output_folder == "no":
	if args.xtcfilename == "no":
		args.output_folder = "depth_insertion_" + args.grofilename[:-4]
	else:
		args.output_folder = "depth_insertion_" + args.xtcfilename[:-4]
if os.path.isdir(args.output_folder):
	print "Error: folder " + str(args.output_folder) + " already exists, choose a different output name via -o."
	sys.exit(1)
else:
	#create folder
	os.mkdir(args.output_folder)
		
	#create log
	filename_log = os.getcwd() + '/' + str(args.output_folder) + '/depth_insertion.log'
	output_log = open(filename_log, 'w')		
	output_log.write("[depth_insertion v" + str(version_nb) + "]\n")
	output_log.write("\nThis folder and its content were created using the following command:\n\n")
	tmp_log="python depth_insertion.py"
	for c in sys.argv[1:]:
		tmp_log += " " + c
	output_log.write(tmp_log + "\n")
	output_log.close()
	

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

def load_MDA_universe():												#DONE
	
	global U
	global all_atoms
	global nb_atoms
	global nb_frames_xtc
	global frames_to_process
	global nb_frames_to_process
	global nb_frames_processed
	global f_start
	global f_end
	global residues_list
	f_start = 0
	nb_frames_processed = 0
		
	#load universe
	#-------------
	if args.xtcfilename == "no":
		print "\nLoading file..."
		U = Universe(args.grofilename)
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = 1
		frames_to_process = [0]
		nb_frames_to_process = 1
	else:
		print "\nLoading trajectory..."
		U = Universe(args.grofilename, args.xtcfilename)
		U_timestep = U.trajectory.dt
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = U.trajectory.numframes

		U.trajectory.rewind()
		#sanity check
		if U.trajectory[nb_frames_xtc-1].time/float(1000) < args.t_start:
			print "Error: the trajectory duration (" + str(U.trajectory.time/float(1000)) + "ns) is shorted than the starting stime specified (" + str(args.t_start) + "ns)."
			sys.exit(1)
		if U.trajectory.numframes < args.frames_dt:
			print "Warning: the trajectory contains fewer frames (" + str(nb_frames_xtc) + ") than the frame step specified (" + str(args.frames_dt) + ")."

		#create list of index of frames to process
		if args.t_end != -1:
			f_end = int((args.t_end*1000 - U.trajectory[0].time) / float(U_timestep))
			if f_end < 0:
				print "Error: the starting time specified is before the beginning of the xtc."
				sys.exit(1)
		else:
			f_end = nb_frames_xtc - 1
		if args.t_start != -1:
			f_start = int((args.t_start*1000 - U.trajectory[0].time) / float(U_timestep))
			if f_start > f_end:
				print "Error: the starting time specified is after the end of the xtc."
				sys.exit(1)
		if (f_end - f_start)%args.frames_dt == 0:
			tmp_offset = 0
		else:
			tmp_offset = 1
		frames_to_process = map(lambda f:f_start + args.frames_dt*f, range(0,(f_end - f_start)//args.frames_dt+tmp_offset))
		nb_frames_to_process = len(frames_to_process)

		
	return
def identify_proteins():												#DONE
	print "\nIdentifying proteins..."
	
	#declare variables
	global proteins_nb
	global proteins_sele
	global proteins_sele_string
	global proteins_sele_string_VMD
	global proteins_boundaries
	global proteins_nb_atoms
	global nb_atom_per_protein
	proteins_nb = 0
	proteins_sele = {}
	proteins_sele_string = {}
	proteins_sele_string_VMD = {}
	proteins_boundaries = {}
	
	#check for protein presence
	if U.selectAtoms("protein").numberOfAtoms() == 0:
		print "Error: no protein detected."
		sys.exit(1)
	
	#declare local variables
	proteins_ca_nb = {}
	proteins_ca_nmax = 0
	proteins_ca_group = {}
	proteins_boundaries = {}

	#retrieve 1st atom info
	proteins_sele["all"] = U.selectAtoms("protein")
	proteins_nb_atoms = proteins_sele["all"].numberOfAtoms()
	prec_resnum = proteins_sele["all"][0].resnum
	prec_segid = proteins_sele["all"][0].segid
	prec_atnum = proteins_sele["all"][0].number+1
	prev_atnum = proteins_sele["all"][0].number+1						#atom corresponding to the beginning of the current protein
	#browse following atoms
	for a in proteins_sele["all"][1:]:
		delta_res = a.resnum-prec_resnum
		delta_atm = a.number+1-prec_atnum
		if delta_res < 0 or a.segid != prec_segid or delta_atm > 1:
			proteins_boundaries[proteins_nb] = [prev_atnum,prec_atnum]
			proteins_nb += 1
			prev_atnum = a.number + 1
		prec_resnum = a.resnum
		prec_atnum = a.number + 1
		prec_segid = a.segid		
	#add last protein section
	if prev_atnum < proteins_sele["all"][proteins_nb_atoms-1].number:
		proteins_boundaries[proteins_nb] = [prev_atnum,proteins_sele["all"][proteins_nb_atoms-1].number+1]
		proteins_nb += 1
	
	#display results
	print " -protein found:", proteins_nb
	print " -protein boundaries (atom numbers): see protein.sele file"
	#create protein selections and save into a txt file
	filename_sele=os.getcwd() + '/' + str(args.output_folder) + '/proteins.sele'
	output_stat = open(filename_sele, 'w')	
	output_stat.write("#This file was generated by the script depth_nisertion v" + str(version_nb) +"\n")
	output_stat.write("#The lines below correspond to MDAnalysis section string, e.g. U.selectAtoms(LINE)\n")
	output_stat.write("\n")	
	for p_index in range(0, proteins_nb):
		progress='\r -creating proteins selections: ' + str(p_index+1) + '/' + str(proteins_nb) + '        '
		sys.stdout.flush()
		sys.stdout.write(progress)
		proteins_sele_string[p_index] = "bynum " + str(proteins_boundaries[p_index][0]) + ":" + str(proteins_boundaries[p_index][1])
		proteins_sele_string_VMD[p_index] = "serial " + str(proteins_boundaries[p_index][0]) + " to " + str(proteins_boundaries[p_index][1])
		proteins_sele[p_index] = U.selectAtoms(proteins_sele_string[p_index])
		output_stat.write(proteins_sele_string[p_index] + "\n")
	output_stat.close()

	nb_atom_per_protein = proteins_sele[0].numberOfAtoms()
	print ""

	#create selection for each particle
	global part_sele
	part_sele = {}
	part_type = []
	for p in range(0,nb_atom_per_protein):
		part_sele[p] = U.selectAtoms("bynum " + str(p+1))
		part_type.append(part_sele[p].resnames()[0])
	
	#create position index of each residue type
	global type_pos
	global res_colour
	type_pos = {}
	res_list = {}
	res_list["basic"] = ['ARG','LYS']
	res_list["polar"] = ['SER','THR','ASN','GLN','HIS']		
	res_list["hydrophobic"] = ['VAL','ILE','LEU','MET','PHE','PRO','CYS','TYR','TRP']
	res_list["backbone"] = ['ALA','GLY']	
	res_colour = {}
	res_colour["basic"] = '#253494'						#blue
	res_colour["polar"] = '#a1d99b'						#greenish
	res_colour["hydrophobic"] = '#993404'				#orange brown
	res_colour["backbone"] = '#969696'					#light grey
	for t in ["basic","polar","hydrophobic","backbone"]:
		type_pos[t] = np.asarray({p: part_type[p] in res_list[t] for p in range(0,nb_atom_per_protein)}.values())

	return
def identify_leaflets():												#DONE
	print "\nIdentifying leaflets..."
	
	#declare variables
	global leaflet_sele
	global leaflet_sele_atoms
	global leaflet_peptide
	
	leaflet_sele = {}
	leaflet_sele_atoms = {}
	leaflet_peptide = "tbd"
	for l in ["lower","upper","both"]:
		leaflet_sele[l] = {}
		leaflet_sele_atoms[l] = {}
	
	#use LeafletFinder:
	leaflet_sele_string = "name PO4 or name PO3 or name B1A"
	ts = U.trajectory[frames_to_process[0]]
	cutoff_value = MDAnalysis.analysis.leaflet.optimize_cutoff(U, leaflet_sele_string)
	L = MDAnalysis.analysis.leaflet.LeafletFinder(U, leaflet_sele_string, cutoff_value[0])
	if np.shape(L.groups())[0]<2:
		print "Error: imposssible to identify 2 leaflets."
		sys.exit(1)
	if L.group(0).centerOfGeometry()[2] > L.group(1).centerOfGeometry()[2]:
		leaflet_sele["upper"] = L.group(0)
		leaflet_sele["lower"] = L.group(1)
	else:
		leaflet_sele["upper"] = L.group(1)
		leaflet_sele["lower"] = L.group(0)
	leaflet_sele["both"] = leaflet_sele["lower"] + leaflet_sele["upper"]
	if np.shape(L.groups())[0] == 2:
		print " -found 2 leaflets: ", leaflet_sele["upper"].numberOfResidues(), '(upper) and ', leaflet_sele["lower"].numberOfResidues(), '(lower) lipids'
	else:
		other_lipids=0
		for g in range(2, np.shape(L.groups())[0]):
			other_lipids += L.group(g).numberOfResidues()
		print " -found " + str(np.shape(L.groups())[0]) + " groups: " + str(leaflet_sele["upper"].numberOfResidues()) + "(upper), " + str(leaflet_sele["lower"].numberOfResidues()) + "(lower) and " + str(other_lipids) + " (others) lipids respectively"
				
	return

#=========================================================================================
# data structures
#=========================================================================================

def struct_particles():
	global z_ref
	global z_part
	z_ref = 0
	z_part = np.zeros(nb_atom_per_protein)

	return
	
#=========================================================================================
# core functions
#=========================================================================================

def calculate_depth():													#DONE
	
	global z_part
	global leaflet_peptide
	global nb_frames_processed
	
	if leaflet_peptide == "tbd":	
		#check whether peptide has reached interfacial state and if so on which leaflet
		dist_min_lower = np.min(MDAnalysis.analysis.distances.distance_array(proteins_sele[0].coordinates(), leaflet_sele["lower"].coordinates(), U.trajectory.ts.dimensions), axis = 1)
		dist_min_upper = np.min(MDAnalysis.analysis.distances.distance_array(proteins_sele[0].coordinates(), leaflet_sele["upper"].coordinates(), U.trajectory.ts.dimensions), axis = 1)
		dist = dist_min_upper - dist_min_lower
		if np.size(dist[dist>0]) == np.size(dist):
			leaflet_peptide = "lower"
		elif np.size(dist[dist>0]) == 0:
			leaflet_peptide = "upper"
			
	elif leaflet_peptide == "upper":
		nb_frames_processed += 1
		z_ref = leaflet_sele["upper"].centerOfGeometry()[2]
		for p in range(0,nb_atom_per_protein):
			tmp = part_sele[p].centerOfGeometry()[2]-z_ref
			z_part[p] += tmp
	else:
		nb_frames_processed += 1
		z_ref = leaflet_sele["lower"].centerOfGeometry()[2]
		for p in range(0,nb_atom_per_protein):
			tmp = z_ref - part_sele[p].centerOfGeometry()[2]
			z_part[p] += tmp
		
	return
def calculate_stats():													#DONE
	
	for p in range(0,nb_atom_per_protein):	
		z_part[p] /= float(nb_frames_processed)
	
	return

#=========================================================================================
# outputs
#=========================================================================================

def depth_graph():
			
	#filenames
	filename_png = os.getcwd() + '/' + str(args.output_folder) + '/depth_insertion.png'
	filename_svg = os.getcwd() + '/' + str(args.output_folder) + '/depth_insertion.svg'

	#create figure
	fig = plt.figure(figsize=(8, 6.2))

	#plot data
	ax = fig.add_subplot(111)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.set_xlim(0, nb_atom_per_protein)
	for t in ["basic","polar","hydrophobic","backbone"]:
		tmp_z = np.zeros(nb_atom_per_protein)
		tmp_z[type_pos[t]] = z_part[type_pos[t]]
		plt.bar(np.arange(0,nb_atom_per_protein), tmp_z, color = res_colour[t], label = t)
	fontP.set_size("small")
	ax.legend(prop=fontP)
	plt.hlines(0, 0, nb_atom_per_protein,)
	plt.xlabel('sequence')
	plt.ylabel('z distance to leaflet')
	
	#save figure
	fig.savefig(filename_png)
	fig.savefig(filename_svg)
	plt.close()
				
	return
def depth_write():
			
	#filenames
	filename_xvg = os.getcwd() + '/' + str(args.output_folder) + '/depth_insertion.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [particles interfacial depth of insertion - written by depth_insertion v" + str(version_nb) + "]\n")
	output_xvg.write("# side of bilayer:\n")
	output_xvg.write("#  -> leaflet = " + str(leaflet_peptide) + "\n")
	output_xvg.write("# nb of frames which contributed to this profile:\n")
	output_xvg.write("# -> weight = " + str(nb_frames_processed) + "\n")
	
	#xvg metadata
	output_xvg.write("@ title \"Particles depth of insertion\"\n")
	output_xvg.write("@ xaxis label \"particles (Nter to Cter)\"\n")
	output_xvg.write("@ yaxis label \"z distance to bilayer (Angstrom)\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 1\n")
	output_xvg.write("@ s0 legend \"avg  z distance\"\n")
	
	#data
	for p in range(0,nb_atom_per_protein):
		results = str(p) + "	" + str(z_part[p])
		output_xvg.write(results + "\n")	
	output_xvg.close()

	return

##########################################################################################
# ALGORITHM
##########################################################################################

#=========================================================================================
#process inputs
#=========================================================================================

#data loading
load_MDA_universe()
identify_proteins()
identify_leaflets()

#create data structures
struct_particles()

#=========================================================================================
# generate data
#=========================================================================================
print "\nCalculating depth of insertion..."

#case: structure only
#--------------------
if args.xtcfilename=="no":
	calculate_density(U.trajectory.ts.dimensions)

#case: browse xtc frames
#-----------------------
else:
	for f_index in range(0,nb_frames_to_process):
		ts = U.trajectory[frames_to_process[f_index]]
		progress = '\r -processing frame ' + str(f_index+1) + '/' + str(nb_frames_to_process) + ' (every ' + str(args.frames_dt) + ' frame(s) from frame ' + str(f_start) + ' to frame ' + str(f_end) + ' out of ' + str(nb_frames_xtc) + ')      '  
		sys.stdout.flush()
		sys.stdout.write(progress)							
		calculate_depth()
	print ""

#=========================================================================================
# process data
#=========================================================================================
calculate_stats()

#=========================================================================================
# produce outputs
#=========================================================================================

print "\nWriting outputs..."
depth_graph()
depth_write()

#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check output in ./" + args.output_folder + "/"
print ""
sys.exit(0)
