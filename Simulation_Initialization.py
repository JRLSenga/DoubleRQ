import matplotlib as mpl
mpl.use('TkAgg')
import glob
from E_UnifyingCode import Main

def Initialize():
    FilePath = "/Simulation Parameters/*.csv" # change file path here
    FileCounter = 1
    for file in glob.glob(FilePath):
        print file
        Main(file, FileCounter)
        FileCounter += 1

Initialize()