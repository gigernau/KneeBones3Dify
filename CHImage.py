import pyvista as pv
import numpy as np
import cupy as cp
import cucim.skimage as cusk
import scipy.spatial as scisp


def _check_coords_in_hull(stats, coords, hull0):

    M = stats.bbox[3] - stats.bbox[0]
    N = stats.bbox[4] - stats.bbox[1]
    P = stats.bbox[5] - stats.bbox[2]

    firstRow   = stats.bbox[0] 
    firstCol   = stats.bbox[1]
    firstPlane = stats.bbox[2]

    bx = cp.meshgrid(cp.linspace(firstRow, M + firstRow -1, M),cp.linspace(firstCol, N + firstCol - 1, N),cp.linspace(firstPlane, P + firstPlane -1, P),indexing='ij' )
    bx = cp.asarray(bx)
    bx = cp.reshape(bx, [3,bx.shape[3]*bx.shape[2]*bx.shape[1]])

    #fig = plt.figure()
    #ax = fig.add_subplot(projection='3d')
    #ax.scatter(coords[:,0],coords[:,1],coords[:,2] )
    #plt.show()
    
    dt = scisp.Delaunay(coords)

    #fig = plt.figure()
    #ax = plt.axes(projection='3d')
    #plot_tri_2(ax, coords, dt)
    #plt.show()

    #Get indices of internal points (non -1)
    idx = cp.asarray(scisp.Delaunay.find_simplex(dt,cp.transpose(bx).get()))

    #non -1 indices are internal points
    # idx[idx != -1]=1
    # idx[idx == -1]=0
    idx = cp.where(idx == -1, 0, 1)

    return idx


def CHImage (CCa, BBa):
    
    LCCa,nr_objects_cca = cusk.measure.label(CCa,return_num=True, connectivity=1)
    Propsa = (cusk.measure.regionprops(LCCa))
    CIa = cp.asarray(Propsa[0].image)
    
    #fig, ax = plt.subplots()
    ## create an IndexTracker and make sure it lives during the whole
    ## lifetime of the figure by assigning it to a variable
    #tracker = IndexTracker(ax, CIa.get())
    #
    #fig.canvas.mpl_connect('scroll_event', tracker.on_scroll)
    #plt.show()

    offset=cp.asarray([BBa[0].get(),BBa[1].get(),BBa[2].get()])
    nze = cp.asarray(cp.nonzero(CIa))
    coords = cp.rot90(nze)
    coords = cp.add(coords, offset)

    hull0 = scisp.ConvexHull(coords.get())
    coords = hull0.points[hull0.vertices]
    coords_in_hull = _check_coords_in_hull(Propsa[0],coords,hull0)
    
    CIa = np.reshape(coords_in_hull, CIa.shape)
    return cp.asarray(CIa, dtype="bool")
