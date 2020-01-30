### Daniel Vogler
###
### evaluate contact stress data
import numpy as np

import matplotlib
import matplotlib.pyplot as pl
import matplotlib.mlab as ml
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.collections as mc

from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D

import math
import sys
import glob, os
import re
import csv

##################

### experiment number 1, 2, 3
experiment_number = "1"
### maximum detected stress magnitudes of 0.0, 0.5, 2.5 and 10.0 MPa
stress_intervals = [0.0, 0.5, 2.5, 10]
### folder with data of multiple experiments
file_path = "/home/davogler/projects/contact_stress_evaluation/data/"
### string to only plot specific load stages
file_locations = ['']
###
file_path = str(file_path + "vd" + experiment_number + "/")

print('\nEvaluating:\n')
print(" file_path: %s \n experiments: %s \n stress intervals: %s\n sub strings: %s\n" %(file_path, experiment_number, stress_intervals, file_locations) )


##################

legendLocation = "lower right"
radius = 0.061
searchName = "ethz_vd"
saveString = './' #file_path

##################

print("\n\nProcessing aperture field data\n\n")

# define measurement locations on sample
fileType = ['X', 'Y', 'Z']
# convert from meters to mm
conversionScale = 1

stringFile = [[] for _ in range(len(file_locations))]
loadStageMPa = []


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)', text) ]


def findFilesInFolder( variableName, file_path ):
	# find all respective files in folder
	counter = -1
	occurrenceCounter = [0]*len(file_locations)
	pathFiles = [[] for _ in range(len(file_locations))]

	# find all files with given string
	os.chdir(file_path)
	for dl in file_locations:
		counter = counter + 1
		variableName = "ethz_vd"
		searchString = str(variableName+"*"+dl+"*Z.txt")
		print("\t Searching partitions for %s \n" %dl)
		#
		for file in glob.glob(searchString):
			occurrenceCounter[counter] = occurrenceCounter[counter] + 1
			file = re.sub('Z.txt', "", file)
			pathFiles[counter].append(file)
			print("\t\t Found partition file %s " %file)
		print("\n\t\t --> %d partition files found for %s \n" %(len(pathFiles[counter]), dl) )

	return pathFiles

# find all files in folder
pathFiles = findFilesInFolder(searchName, file_path)
pathFiles[0].sort(key=natural_keys)


def loadArray(file_path, filename, fileType ):

	# loop through aperture table files and load individually
	for j in range(len(fileType)):
		# find all files with given string
		fileToLoad = str(file_path+filename+fileType[j]+".txt")
		print("\t Opening file %s \n" %fileToLoad)
		# loading aperture files
		values = np.genfromtxt(fileToLoad, delimiter=" ")
		# sorting in X, Y and Z coordinates
		if( j == 0 ):
			arrayX = values
		elif( j == 1):
			arrayY = values
		elif( j == 2):
			arrayZ = values

	return arrayX, arrayY, arrayZ


def computeContactArea(file_path,fileType,pathFiles,stress_intervals):

	contactAreaPercent = [[] for _ in range(len(stress_intervals))]

	# loop through film stress intensity levels
	for k in range(len(pathFiles[0])):

		print("\n Computing contact area for:" )
		print(" path: %s \n file: %s\n\n" %(file_path, pathFiles[0][k]) )

		# extract load in mpa
		loadStage = re.sub('^(.*)(?=_)',"", pathFiles[0][k])
		loadStage = re.sub("\_", "", loadStage)
		loadStage = re.sub("\.", "", loadStage)
		loadStage = float(loadStage)
		loadStage = loadStage/3.14159/radius**2/1000.0
		loadStageMPa.append(loadStage)

		# load individual grid arrays
		csFieldGridX, csFieldGridY, csFieldGridZ = loadArray(file_path,pathFiles[0][k], fileType)

		# size of loaded array
		gridSize = np.shape(csFieldGridZ)

		# initialize
		fractureInContact = [[] for _ in range(len(stress_intervals))]
		fractureTotalArea = [[] for _ in range(len(stress_intervals))]

		for l in range(len(stress_intervals)):
			# compute contact area
			for i in range(gridSize[0]):
				for j in range(gridSize[1]):
					if (csFieldGridX[i][j]**2 + csFieldGridY[i][j]**2)**(0.5) < radius and csFieldGridZ[i][j] > stress_intervals[l]: 
						fractureInContact[l].append(csFieldGridZ[i][j])
					if (csFieldGridX[i][j]**2 + csFieldGridY[i][j]**2)**(0.5) < radius: 
						fractureTotalArea[l].append(csFieldGridZ[i][j])


			contactAreaPercent[l].append(float(len(fractureInContact[l]))/float(len(fractureTotalArea[l]))*100.0 )

	plotFigure(loadStageMPa, contactAreaPercent, file_path, saveString, file_locations)

	return


def plotFigure(loadStageMPa, contactAreaPercent, file_path, saveString, file_locations):

	print("\n\n\nPlotting Figures \n")

	fontSize = 25
	markerSize = 10
	lineStyle = ""

	pl.figure(num=None, figsize=(10, 8), dpi=80, facecolor='w', edgecolor='k')
	# settings
	font = {'size'   : fontSize}
	matplotlib.rc('font', **font)
	# colormap
	colormap = cm.gnuplot
	normalize = mcolors.Normalize( vmin=0.0, vmax=len(stress_intervals) )

	for l in range(len(stress_intervals)):
		label = "$\sigma_{Min}$ = "+str(stress_intervals[l])+" MPa"
		pl.plot( loadStageMPa, contactAreaPercent[l], c=colormap(normalize(l)), marker='o', label=label, markersize=markerSize, linestyle=lineStyle)

	pl.ylabel('contact area [%]', fontsize = fontSize)
	pl.xlabel('vertical stress $\sigma_Z$ [MPa]', fontsize = fontSize)
	pl.legend(loc=legendLocation, numpoints = 1)
	pl.grid(b=True, which='major', color='lightgrey', linestyle='-')
	pl.tick_params(axis='both', which='major', labelsize=fontSize)
	pl.savefig( str(file_path+"contactAreaPercent.png"), bbox_inches = 'tight' )

	fileToSave = str(file_path+"contactAreaPercent_loadStageMPa.txt")
	with open(fileToSave, 'w') as csvfile:
		writer = csv.writer(csvfile,delimiter=' ')
		writer.writerow(['{:1.2f}'.format(x) for x in loadStageMPa])

	fileToSave = str(file_path+"contactAreaPercent_contactAreaPercent.txt")
	with open(fileToSave, 'w') as csvfile:
		writer = csv.writer(csvfile,delimiter=' ')
		writer.writerow(['{:1.2f}'.format(x) for x in contactAreaPercent[l] ])

	return


computeContactArea(file_path,fileType,pathFiles,stress_intervals)


pl.show()
exit()

