import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import csv
import math
import glob
import os


def main():

    # Define the number of data columns for each file (hardcoded now, switch to make more generic)
    nRun = 10 # number of runs
    nCol = 12 # umber of columns; (or parameters)

    # Get Filenames For Model Outputs
    mypath = os.path.dirname(os.path.realpath(__file__)) # sets the directory containing the python script as the directory containing the model Outputs
    extension = 'csv' # sets the extension to .csv
    os.chdir(mypath) # sets the directory to the directory defined by mypath
    filename = [i for i in glob.glob('*turingBatch.{}'.format(extension))] # creats a list containing all .csv files in directory

    # Initialize data storage arrays
    dataArray_double = np.full((nRun,nCol,),np.nan) # data array for the signal output (-999)
    dataArray_fileID = np.full((nRun),np.nan) # create a separate array for the strings

    # Generate Data Arrays
    dataArray_fileID = np.genfromtxt(filename[0], delimiter=",",skip_header=1,usecols=0, dtype=str) # temporary data variable
    dataArray_double = np.genfromtxt(filename[0], delimiter=",",skip_header=1,usecols=np.arange(1,nCol)) # temporary data variable
    headerNames = ["fileID","runType","nPhotons","FOV","refreal","refIm","a","b","Dmin","Dmax","K","Junge"]

    for i in range(nRun):
        fileOut = 'lidarMCinput'+str(i+1)+'.csv'
        with open(fileOut,'w') as csvfile:
            writer = csv.writer(csvfile, delimiter=',')
            for j in range(nCol):
                if j==0:
                    writer.writerow([headerNames[j],dataArray_fileID[i]]) # [headernames[j]] this wraps each string instead of printing char
                else:
                    writer.writerow([headerNames[j],dataArray_double[i,j-1]]) # [headernames[j]] this wraps each string instead of printing char
if __name__ == "__main__":
    main()
