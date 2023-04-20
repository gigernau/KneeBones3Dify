import cupy
from VolumeViewer import VolumeViewer
#Crop di un volume 3D per eliminare i bordi con valori inferiori a SogliaCrop
def MRIcropCoordsRev(V, SogliaCrop,w,h,d):

	V0 = cupy.zeros(V.shape)
	cupy.copyto(V0,V)
	V0[V0 < SogliaCrop] = 0


	##############################################################################
	###  Ricerca indici per regione di soli 0
	###############################################################################
	# Crop di un volume 3D per eliminare i bordi tutti nulli
	x1= 0
	x2= 0
	y1= 0
	y2= 0
	z1= 0
	z2 = 0

	for i in range(d):
	  if (cupy.argwhere((V0[i,:,:]) > 0).size > 0):
	    x1 = i
	    break


	for i in range(d-1,0,-1):
	  if (cupy.argwhere((V0[i,:,:]) > 0).size > 0):
	    x2 = i
	    break


	for i in range(h):
	  if (cupy.argwhere((V0[:,i,:]) > 0).size > 0):
	     y1 = i
	     break


	for i in range(h-1,0,-1):
	  if (cupy.argwhere((V0[:,i,:]) > 0).size > 0):
	    y2 = i
	    break


	for i in range(w):
	  if (cupy.argwhere((V0[:,:,i]) > 0).size > 0):
	    z1 = i
	    break


	for i in range(w-1,0,-1):
	  if (cupy.argwhere((V0[:,:,i]) > 0).size > 0):
	    z2 = i
	    break

	
	Vcropped = V[x1:x2+1,y1:y2+1,z1:z2+1]

	return Vcropped,x1,x2,y1,y2,z1,z2
