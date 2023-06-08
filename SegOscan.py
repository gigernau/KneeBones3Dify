#KneeBones3Dify main

import sys
import os
import argparse

import SimpleITK as sitk
import pyvista as pv

import cupy as cp
import numpy as np

import cucim.skimage as cusk
from cucim.skimage.morphology import (cube,ball,square)
import cupyx.scipy.ndimage as cuSci
import scipy.io as sciio
from cucim.skimage import filters,measure,feature,segmentation
from cucim.skimage.exposure import rescale_intensity


#my functions
from MRIcropCoordsRev import MRIcropCoordsRev
from VolumeViewer import VolumeViewer
from MyReadDICOM import MyReadDICOM,swap
from EdgeBySlice import EdgeBySlice
from CHImage import CHImage
import Gui as g

import timeit
import tkinter as tk
from tkinter import filedialog,PhotoImage,Canvas
from PIL import Image, ImageTk
from subprocess import run,check_output


#Imports the correct "morphology" package based on the nvcc NVIDIA compiler version
#- executes the "nvcc --version" command and captures the output
output = check_output(["nvcc", "--version"])
#- decodes the output into a string
output_str = output.decode("utf-8")
#- looks for the string "release" in the decoded output string
release_line = [line for line in output_str.split("\n") if "release" in line][0]
#- extracts the sixth token (version), with the second and third version characters
version = release_line.split(" ")[5][1:3]
#- imports the correct package
if int(version) <= 11:
    from cupyx.scipy.ndimage import morphology as morphology
else:
    from cupyx.scipy.ndimage import _morphology as morphology 

import ctypes
import os
#Handling of C++ code
dir_path = os.path.dirname(os.path.realpath(__file__))
if sys.platform == "linux" or sys.platform == "linux2":
    # linux
    handle = ctypes.CDLL(dir_path + "/smoothPatch.so")
elif sys.platform == "win32":
    # Windows...
    handle = ctypes.CDLL(dir_path + "\\smoothPatch.dll")

    
handle.smoothPatch.argtypes = [np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),ctypes.c_int,
                                np.ctypeslib.ndpointer(dtype=np.int64, ndim=1, flags='C_CONTIGUOUS'),ctypes.c_int]
handle.smoothPatch.restype = ctypes.c_void_p




#Functions
def smoothPatch(K,c,K1,c1):
    return handle.smoothPatch(K,c,K1,c1) 


#Define global vars
if sys.platform == "linux" or sys.platform == "linux2": #Linux
    dataset = "./Datasets/"

elif sys.platform == "win32": #Windows
    dataset =".\\Datasets\\"

elif sys.platform == "darwin": # OS X
    dataset = "./Datasets/"

SogliaCrop = 600
CHadd = 6
FinalClosing = 8
Protrus = 3
Edges = 1


first = True
fin = True

def main():

    
    global dataset, SogliaCrop, CHadd, FinalClosing, Protrus, Edges, flags,fin

    

    while(fin):
        sG()
        global Tstart
        Tstart = timeit.default_timer()
        
        try:
            for f in flags:
                f()
        except:
            print("                        ################")
            print("                        #### ERROR! ####")
            print("                        ################\n\n")
            g.GuiError()
            continue
        fin = g.GuiFin()


def sG():
    global dataset, SogliaCrop, CHadd, FinalClosing, Protrus, Edges,flags,first
    datasetG, SogliaCropG, CHaddG, FinalClosingG, ProtrusG, EdgesG = g.Gui(dataset, SogliaCrop, CHadd, FinalClosing, Protrus, Edges)
    
    ##############################################################################
    ###  START MAIN and PARAMETERS
    ###############################################################################
    print("\n\n###############################")
    print("### Input Data")
    print("###############################\n")
    print("Dataset:", datasetG)
    print("SogliaCrop:", SogliaCropG)
    print("CHadd:", CHaddG)
    print("FinalClosing:", FinalClosingG)
    print("Protrus:", ProtrusG)
    print("Edges:", EdgesG)

    msg = f"""\n\n###############################"
    ### Input Data"
    ###############################\n"
    Dataset: {datasetG}
    SogliaCrop: {SogliaCropG}
    CHadd: {CHaddG}
    FinalClosing: {FinalClosingG}
    Protrus: {ProtrusG}
    Edges: {EdgesG}"""

    if(first == False):
        print("\n\n###############################")
        if(datasetG != dataset):
            flags = np.asarray([s0,s1,s2,s3,s4,s5,s6,s7,s8])
            print("### Restart from state 0")
        elif(SogliaCropG != SogliaCrop):
            flags = np.asarray([s1,s2,s3,s4,s5,s6,s7,s8])
            print("### Restart from state 1")
        elif(CHaddG != CHadd):
            flags = np.asarray([s2,s3,s4,s5,s6,s7,s8])
            print("### Restart from state 2")
        elif(FinalClosingG != FinalClosing):
            flags = np.asarray([s3,s4,s5,s6,s7,s8])
            print("### Restart from state 3")
        elif(ProtrusG != Protrus):
            flags = np.asarray([s4,s5,s6,s7,s8])
            print("### Restart from state 4")
        elif(EdgesG != Edges):
            flags = np.asarray([s5,s6,s7,s8])
            print("### Restart from state 5")            
        print("###############################\n")
                
    
    else:
        first = False
        flags = np.asarray([s0,s1,s2,s3,s4,s5,s6,s7,s8])


    #Structuring elements:
    global SEbordiScuri,SEseparaOssa,SS3s
    SEbordiScuri = cusk.morphology.cube(2,dtype=cp.bool_) #to enhance dark edges around the bones
    SEseparaOssa = cusk.morphology.ball(7,dtype=cp.bool_) #to separate bones from the other structures (e.g., ligaments and fat)
    SS3s = cusk.morphology.square(3,dtype=cp.bool_)       #for a finer segmentation (needed for patella)
    global SEpp1 #for Post-Processing 1 (imclose) to refine border regions
    global SEpp2 #for Post-Processing 2 (imopen)  to eliminate final protrusions (calcifications) 
    global SEpp3 #for Post-Processing 3 (imdilate) to dilate border regions

    dataset, SogliaCrop, CHadd, FinalClosing, Protrus, Edges = datasetG, SogliaCropG, CHaddG, FinalClosingG, ProtrusG, EdgesG

    ##########################################################################


def s0():
    ##############################################################################
    ###  1) Read DICOM
    ##############################################################################
    t0 = timeit.default_timer()
    global V, StrelRotula, w, h, d, dataset,spacing
    V, StrelRotula, w, h, d,spacing = MyReadDICOM(dataset)
    t1 = timeit.default_timer()
    print("\n\n###############################")
    print("### SEGMENTATION TIMES")
    print("###############################")
    print(f"\nTime for reading MRI file: {t1 - t0}s")

def s1():
    ##############################################################################
    ###  2) Pre-processing
    ##############################################################################   
    global spacing
    ######################################################################
    # 2.1) Crops dark areas around the knee
    Vcrop,x1,x2,y1,y2,z1,z2 = MRIcropCoordsRev(V,SogliaCrop,w,h,d)
    Vcrop = cp.rot90(Vcrop,axes=(2,0))
    spacing = np.asarray(spacing)
    spacing[0],spacing[2] = swap(spacing[0],spacing[2])

    ######################################################################
    # 2.2) Enhances dark edges around the bones
    structCube0 = cp.asarray(SEbordiScuri)
    erode = morphology.grey_erosion(Vcrop, footprint=structCube0)
    Img3D = cp.asarray(cp.double(erode.get()))
    # Normalize in [0,1]
    Img3D = rescale_intensity(Img3D, out_range=(0, 1))


    ##############################################################################
    ###  3) Segmentation and extraction of the 3 bones (--> CCa2, CCb2, CCc2)
    ##############################################################################

    ######################################################################
    # 3.1) Segmentation with Otsu thresholding  

    t0 = timeit.default_timer() 

    T = filters.threshold_otsu(Img3D)
    global OutSeg
    OutSeg = Img3D>T

    print(f"Time for segmentation: {timeit.default_timer() - t0}s")

    ######################################################################
    # 3.2) Extraction of a "rough" version of tibia and femur (--> CCa, CCb)

    ##########################################################
    # 3.2.1) "Strong" erosion, to separate the bones from the remaining structures (e.g., ligaments and fat)
    t0 = timeit.default_timer() 

    OutSeg_1 = morphology.binary_erosion(OutSeg, structure=SEseparaOssa)    

    print(f"Time for erosion: {timeit.default_timer() - t0}s")

    ##########################################################
    # 3.2.2) Detects edges and eliminates them from the previous segmentation
    t0 = timeit.default_timer() 

    global edge3D
    edge3D = EdgeBySlice(Img3D)

    print(f"Time for computing edges: {timeit.default_timer() - t0}s")

    OutSeg_2 = OutSeg_1*np.logical_not(edge3D)

    ##########################################################
    #3.2.3) Computes the connected components (CCs) in OutSeg_2
    L0,nr_objects0 = cusk.measure.label(OutSeg_2,return_num=True, connectivity=1) 

    ##########################################################
    # 3.2.4) Selects the 3 CCs with larger volume
    t0 = timeit.default_timer() 

    stats0 = (measure.regionprops(L0))

    print(f"Time for computing the 3 widest CCs: {timeit.default_timer() - t0}s")

    stats0.sort(key=lambda x: x.area, reverse=True)

    cc012 = cp.zeros(L0.shape, dtype=cp.bool_)
    for i in range(0,3):
      cc012[L0.get() == stats0[i].label] = True
    Vol13 = cp.asarray(cc012)

    ##########################################################
    # 3.2.5) Selects the 2 CCs with higher Extent (--> CCa, CCb)
    L1,nr_objects1 = cusk.measure.label(Vol13,return_num=True, connectivity=1)
    stats1 = (measure.regionprops(L1))

    stats1.sort(key=lambda x: x.extent, reverse=True)

    CCa= cp.zeros(L1.shape, dtype=cp.bool_)
    CCa[L1.get() == stats1[0].label] = True
    CCb= cp.zeros(L1.shape, dtype=cp.bool_)
    CCb[L1.get() == stats1[1].label] = True
    CCa = cp.asarray(CCa)
    CCb = cp.asarray(CCb)

    ##########################################################
    # 3.2.6) Computes the convex hull (CH) of tibia and femur (needed in step 3.3 for patella)

    BBa = cp.zeros(6)
    for i in range(6):
      BBa[i] = round(stats1[0].bbox[i])

    t0 = timeit.default_timer()
    CIa = CHImage(CCa, BBa)
    print(f"Time for computing convex image: {timeit.default_timer() - t0}s")

    global Mtfa
    Mtfa = cp.zeros(Vol13.shape,dtype="bool")
    Mtfa[BBa[0]:BBa[3], BBa[1]:BBa[4], BBa[2]:BBa[5]]=CIa

    BBb = cp.zeros(6)
    for i in range(6):
      BBb[i] = round(stats1[1].bbox[i])

    t0 = timeit.default_timer()
    CIb = CHImage(CCb, BBb)
    print(f"Time for computing convex image CIb: {timeit.default_timer() - t0}s")
    global Mtfb
    Mtfb = cp.zeros(Vol13.shape,dtype="bool")
    Mtfb[BBb[0]:BBb[3], BBb[1]:BBb[4], BBb[2]:BBb[5]]=CIb


    ######################################################################
    # 3.3) Extracts a "rough" version of the patella (--> CCcs, in the subvolume)

    ##########################################################
    # 3.3.1) Constructs the volume regions where to look for patella (excluding tibia and femur)

    BW1CCab = cp.logical_xor(OutSeg, Mtfa | Mtfb)

    ##########################################################
    # 3.3.2) Extracts the subvolume SubVol for the patella (top-left related to CCa and CCb)  
    if (BBa[1]<BBb[1]):
      BBtibia = BBb
      BBfemore = BBa
    else:
      BBtibia = BBa
      BBfemore = BBb;

    SubVol = Img3D[0:BBfemore[0],0:BBtibia[1],:]

    # Patella: operations in the SubVol (all the variables ending with 's')
    BW1s = BW1CCab[0:BBfemore[0],0:BBtibia[1],:] 

    ##########################################################
    #3.3.3) Enforces the separation of the CCs (as done for tibia and femur)   
    t0 = timeit.default_timer()
    edge3Ds = EdgeBySlice(SubVol); #Calcolo degli edge nel SubVol
    print(f"Time for computing edges in SubVol: {timeit.default_timer() - t0}s")
      
    BW2s = BW1s*cp.logical_not(edge3Ds) #Segmentazione con separazione fra CC
    BW2s = cp.asarray(BW2s, dtype="float64")
    BW2es = morphology.binary_erosion(BW2s, structure=StrelRotula)
       
    ##########################################################
    # 3.3.4) Selects the CC with highest volume 
    L0s,nr_objects_cc0 = cusk.measure.label(BW2es,return_num=True, connectivity=1)
    stats0s = (measure.regionprops(L0s))

    stats0s.sort(key=lambda x: x.area, reverse=True)
    CCcs = cp.zeros(L0s.shape, dtype=cp.bool_)
    CCcs[L0s.get() == stats0s[0].label] = True
       
    #Patella: end of operations in the SubVol

    ##########################################################
    # 3.3.5) Computes the patella CH (needed in Step 3.4)
    Ls,nr_objects_cc0 = cusk.measure.label(CCcs,return_num=True, connectivity=1)
    Propsc = (measure.regionprops(Ls))
    BBc = cp.zeros(6)
    for i in range(6):
      BBc[i] = round(Propsc[0].bbox[i])

    t0 = timeit.default_timer()
    CIc = CHImage(CCcs, BBc)
    print(f"Time for computing convex image CIc: {timeit.default_timer() - t0}s")
    global Mtfc
    Mtfc = cp.zeros(Vol13.shape,dtype="bool")
    Mtfc[BBc[0]:BBc[3], BBc[1]:BBc[4], BBc[2]:BBc[5]]=CIc


def s2():
    ######################################################################
    # 3.4) Extracts a refined version of the 3 bones (--> CCa2, CCb2, CCc2)

    ##########################################################
    # 3.4.1) Dilates the 3 CHs
    sphMT = cusk.morphology.ball(CHadd)
    Mtfa2 = cusk.morphology.binary_dilation(Mtfa,footprint=sphMT)
    Mtfb2 = cusk.morphology.binary_dilation(Mtfb,footprint=sphMT)
    Mtfc2 = cusk.morphology.binary_dilation(Mtfc,footprint=sphMT)

    ##########################################################
    # 3.4.2) Re-computes OutSeg_4 (eroded version of OutSeg without edges, using SS3) 
    OutSeg_3 = cp.asarray(OutSeg*cp.logical_not(edge3D),dtype="float64")
    OutSeg_4 = OutSeg_3

    for slice in range(0,OutSeg_3.shape[2]):
      OutSeg_4[:,:,slice] = morphology.grey_erosion(OutSeg_3[:,:,slice], footprint=SS3s)

    OutSeg_4 = cp.asarray(OutSeg_4,dtype="bool")
    ##########################################################
    #3.4.3) Extracts tibia and femur, eliminating further CCs introduced in OutSeg_4
    BonesNotCloseabCH = OutSeg_4 & (Mtfa2 | Mtfb2)

    L02,nr_objects_cc02 = cusk.measure.label(BonesNotCloseabCH,return_num=True, connectivity=1)
    stats02 = (measure.regionprops(L02))
    stats02.sort(key=lambda x: x.area, reverse=True)

    global CCa2, CCb2
    CCa2 = cp.zeros(L02.shape, dtype=cp.bool_)
    CCa2[L02.get() == stats02[0].label] = True
    CCb2 = cp.zeros(L02.shape, dtype=cp.bool_)
    CCb2[L02.get() == stats02[1].label] = True

    ##########################################################
    #3.4.4) Extracts the patella, eliminating further CCs introduced in OutSeg_4
    BonesNotClosecCH = OutSeg_4 & Mtfc2; 
    L03,nr_objects_cc03 = cusk.measure.label(BonesNotClosecCH,return_num=True, connectivity=1)
    stats03 = (measure.regionprops(L03))

    stats03.sort(key=lambda x: x.area, reverse=True)

    global CCc2
    CCc2 = cp.zeros(L03.shape, dtype=cp.bool_)
    CCc2[L03.get() == stats03[0].label] = True

    CCa2 = np.pad(CCa2,15, 'constant', constant_values=0)
    CCb2 = np.pad(CCb2,15, 'constant', constant_values=0)
    CCc2 = np.pad(CCc2,15, 'constant', constant_values=0)


    ##############################################################################
    ###  4) Post-processing delle 3 ossa CCa2, CCb2, CCc2 (--> BonesClose)
    ##############################################################################
    print('\nPost-processing:**************************************************')



def s3():
    ######################################################################
    # 4.1) Closing of the 3 CCs separately to refine border regions
    SEpp1 = cusk.morphology.ball(FinalClosing,dtype=cp.bool_)
    t0 = timeit.default_timer()
    global CCa2close,CCb2close,CCc2close
    CCa2close = cusk.morphology.binary_closing(CCa2,footprint=SEpp1)
    CCb2close = cusk.morphology.binary_closing(CCb2,footprint=SEpp1)
    CCc2close = cusk.morphology.binary_closing(CCc2,footprint=SEpp1)
    print(f"Time for Post-Processing 1 (imclose): {timeit.default_timer() - t0}s")



def s4():
    ######################################################################
    # 4.2) Opening of the 3 CCs separately to eliminate protrusions
    SEpp2 = cusk.morphology.ball(Protrus,dtype=cp.bool_) 
    t0 = timeit.default_timer()
    global Sopena
    Sopena = cusk.morphology.binary_opening(CCa2close,footprint=SEpp2)
    global Sopenb
    Sopenb = cusk.morphology.binary_opening(CCb2close,footprint=SEpp2)
    global Sopenc
    #Sopenc = cusk.morphology.binary_opening(CCc2close,footprint=SEpp2)
    Sopenc = CCc2close
    print(f"Time for Post-Processing 2 (imopen): {timeit.default_timer() - t0}s")


def s5():

    #####################################################################
    # 4.3) Dilation to fill-in edges
    SEpp3 = cusk.morphology.ball(Edges,dtype=cp.bool_)
    t0 = timeit.default_timer()
    SopenDila = cusk.morphology.binary_dilation(Sopena,footprint=SEpp3)
    SopenDilb = cusk.morphology.binary_dilation(Sopenb,footprint=SEpp3)
    SopenDilc = cusk.morphology.binary_dilation(Sopenc,footprint=SEpp3)
    global BonesClose
    BonesClose = SopenDila | SopenDilb | SopenDilc
    print(f"Time for Post-Processing 3 (imdilate): {timeit.default_timer() - t0}s")


def s6():

    ##############################################################################
    ###  5) Volume Viewer with PyVista
    ##############################################################################

    BonesClose2 = cp.zeros([BonesClose.shape[0] +2, BonesClose.shape[1],BonesClose.shape[2]],dtype="bool")
    BonesClose2[1:BonesClose2.shape[0]-1,:,:]=BonesClose

    ##############################################################################
    ###  6) Computes the isosurface
    ##############################################################################
    global data1,spacing
    data1 = np.invert(np.asarray(BonesClose2.get()*1))

    data1 = pv.wrap(data1)
    data1.spacing=spacing    

    #Triangulation
    t0 = timeit.default_timer()
    data1 = data1.contour(isosurfaces = 3)
    print(f"\nTime for computing contour (triangulation): {timeit.default_timer() - t0}s")


def s7():

    ##############################################################################
    ###  7) Smoothpatch
    ##############################################################################

    vertices = data1.points.flatten()
    faces = data1.faces.flatten()

    t0 = timeit.default_timer()
    data1Smoothed = smoothPatch(vertices,len(vertices),faces,len(faces))
    print(f"Time for computing smootPatch: {timeit.default_timer() - t0}s")

    data1.points = vertices.reshape(-1,3)
    data1.faces = faces.reshape(-1,4)


def s8():
    ##############################################################################
    ###  8) STL View & Save
    ##############################################################################
    global Tstart,dataset,data1
    #TIME
    print(f"\nTime for total execution: {timeit.default_timer() - Tstart}s")


    data1.flip_normals()
    data1.compute_normals(cell_normals=True, point_normals=False, split_vertices=False, 
                          flip_normals=True, consistent_normals=True, auto_orient_normals=True, 
                          non_manifold_traversal=True,feature_angle=30.0, inplace=False, progress_bar=False)

    
    #STL save...
    if( str(dataset)[-1] == '/' or str(dataset)[-1] == '\\'):
        nm = str(f"{dataset[:-1]}.stl")
    else:
        nm = str(f"{dataset}.stl")
    data1.save(nm)
    
    pl = pv.Plotter()
    pl.add_mesh(data1,color="white")

    pl.show()
    pl.close()


flags = np.asarray([s0,s1,s2,s3,s4,s5,s6,s7,s8])


main()
