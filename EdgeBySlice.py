import cupy as cp
import numpy as np
from cucim.skimage import feature,exposure
from cupyx.scipy.ndimage import gaussian_gradient_magnitude,convolve,convolve1d,gaussian_filter




def  smoothGradient(I, sigma):

    filterExtent = np.ceil(4*sigma)
    x = np.asarray(range(int(-filterExtent),int(filterExtent+1)))

    
    c = 1/(np.sqrt(2*np.pi)*sigma)
    gaussKernel = c * np.exp(-(np.square(x)/(2*np.square(sigma))))
    
    
    gaussKernel = gaussKernel/np.sum(gaussKernel)

    
    derivGaussKernel = cp.gradient(gaussKernel)

    # gaussKernel = cp.asarray(cp.reshape(gaussKernel,[1,gaussKernel.shape[0]]))
    # derivGaussKernel = cp.asarray(cp.reshape(derivGaussKernel,[1,derivGaussKernel.shape[0]]))

    gaussKernel = cp.asarray(gaussKernel)
    derivGaussKernel = cp.asarray(derivGaussKernel)

    negVals = derivGaussKernel < 0
    posVals = derivGaussKernel > 0
    derivGaussKernel[posVals] = derivGaussKernel[posVals]/np.sum(derivGaussKernel[posVals])
    derivGaussKernel[negVals] = derivGaussKernel[negVals]/np.abs(np.sum(derivGaussKernel[negVals]))
    
    
    GX = convolve1d(I, gaussKernel,mode="nearest", axis=1)
    GX = convolve1d(GX, derivGaussKernel, mode="nearest",axis=0)


    GY = convolve1d(I, gaussKernel, mode="nearest",axis=0)
    GY  = convolve1d(GY, derivGaussKernel,mode="nearest",axis=1)

    return GX, GY

def OurCanny(I):
    PercentOfPixelsNotEdges = .7
    ThresholdRatio = .4
    sigma = np.sqrt(2)

    #magGrad = gaussian_gradient_magnitude(I, sigma)
    dx,dy = smoothGradient(I, sigma)
    
    magGrad = np.hypot(dx, dy)
    magmax = cp.max(magGrad[:])
    if magmax > 0:
        magGrad = magGrad / magmax

    m = magGrad.shape[0]
    n = magGrad.shape[1]
    
    counts , loc = exposure.histogram(magGrad, nbins=64)
    #print(counts)
    ind = cp.asarray(cp.where(cp.cumsum(counts) > PercentOfPixelsNotEdges * m * n))

    highTresh = ind[0,0]/64
    lowThresh = ThresholdRatio*highTresh

    edge = feature.canny(I, sigma, mode='nearest', low_threshold=lowThresh, high_threshold=highTresh)
    return edge


#Calcola gli edge del volume I slice per slice
def EdgeBySlice(I):
    #print(f"shape per edge {I.shape}")
    edge3D=cp.zeros(I.shape,dtype="bool")

    #edge3D[:,:,50] = OurCanny(I[:,:,50])
    for slice in range(0,I.shape[2]):
        edge3D[:,:,slice] = OurCanny(I[:,:,slice])
        #edge3D[:,:,slice] = feature.canny(I[:,:,slice], sigma = np.sqrt(2),mode='nearest')

    
    return edge3D
