import mineazimutalAverage
import numpy
import matplotlib.pyplot as pyplot
from correlate2d import correlate2d
def PSD2(image, image2=None, binsize=1./16.):
    """
    Two-dimensional Power Spectral Density.
    NAN values are treated as zero.

    azbins - Number of azimuthal (angular) bins to include.  Default is 1, or
        all 360 degrees.  If azbins>1, the data will be split into [azbins]
        equally sized pie pieces.  Azbins can also be a numpy array.  See
        AG_image_tools.azimuthalAverageBins for details
    """
    
    # prevent modification of input image (i.e., the next two lines of active code)
    image = image.copy()
    Nx=len(image)
    Ny=len(image[0])
    
    # remove NANs (but not inf's)
    image[image!=image] = 0
    F1 = numpy.fft.fft2(image)
    
    if image2 is None:
       image2 = image
       F1bis = F1
    else:
       image2 = image2.copy()
       image2[image2!=image2] = 0
       F1bis = numpy.fft.fft2(image2)

    #Deconvolve the window function

#    pxwind=numpy.zeros(len(image)*len(image[0])).reshape((len(image), len(image[0])))
#    for i in range(0,len(pxwind)):
#         sinci = numpy.sinc((i-Nx/2.)/Nx)
#         for j in range(0,len(pxwind[0])):
#            sincj = numpy.sinc((j-Ny/2.)/Ny)            
#            pxwind[i,j] =  sinci*sincj              

#    pyplot.imshow( numpy.log10( pxwind**2 ))
#    pyplot.colorbar()
#    pyplot.savefig( 'mask.pdf')
#    pyplot.close()


    F2 = numpy.fft.fftshift( (F1) )
    F2bis = numpy.fft.fftshift( (F1bis) )

    psd2 = numpy.abs( (F2*F2bis) ) 

    # normalization is approximately (numpy.abs(image).sum()*numpy.abs(image2).sum())

    
    return psd2

def pspec(psd2, binsize=3.0):
    """
    Create a Power Spectrum (radial profile of a PSD) from a Power Spectral Density image
    """
    (freq,zz,err) = mineazimutalAverage.azimuthalAverage(psd2, binsize=binsize)
    freq = freq.astype('float') 
    return_vals = list(zz)

    return freq,return_vals,err

