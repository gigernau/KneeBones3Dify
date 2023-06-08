#MyReadDICOM: reads a DICOM dataset and returns:
#V           = 3D image volume of the MRI in the sagittal view, with the patella on the top-left side
#w, h, d     = dimensions of the (transformed) volume V
#spacing     = slice spacing of the original DICOM
#StrelRotula = structuring element adopted for patella segmentation

import sys
import SimpleITK as sitk
import cupy as cp
import cucim.skimage as cusk
from cucim.skimage.morphology import (cube,ball,square)


def swap(a,b):
    return b,a
  
  
def MyReadDICOM(dataset):

  print("\n###############################")
  print("### DICOM INFO")
  print("###############################")
  print("\nReading DICOM directory:", dataset)
  reader = sitk.ImageSeriesReader()

  dicom_names = reader.GetGDCMSeriesFileNames(dataset)
  reader.SetFileNames(dicom_names)

  Vorig = reader.Execute()
  w = Vorig.GetWidth()
  h = Vorig.GetHeight()
  d = Vorig.GetDepth()
  spacing = Vorig.GetSpacing()
  Ori= Vorig.GetDirection()
  Vorig = sitk.GetArrayFromImage(Vorig)

  print("DICOM dimensions:",Vorig.shape)

  Ori = tuple([int(round(x,2)) if isinstance(x, float) else x for x in Ori])

  if (Ori[0]): #Coronal or Axial
    if (Ori[4]): #Axial
      print("MRI type: Axial")
      Vorig = cp.rot90(Vorig, k=2,axes=(2,0))
      w,h,d = h,d,w
      VS = cp.flipud(cp.rot90(cp.transpose(Vorig,[0,2,1])))
      V=VS 
      

    else: #%Coronal
      print("MRI type: Coronal")
      Vorig = cp.rot90(cp.rot90(Vorig,axes=(2,1)),axes=(0,1))
      w,h,d = d,h,w
      VS=cp.fliplr(cp.transpose(Vorig,[0,2,1]))
      V=VS 
    
    StrelRotula = cusk.morphology.cube(3,dtype=cp.bool_) 

  else: #Sagittal
          print("MRI type: Sagittal")
          StrelRotula = cusk.morphology.cube(4,dtype=cp.bool_)
          V=Vorig
          
  V = cp.asarray(V)

  return V, StrelRotula, w, h, d, spacing
