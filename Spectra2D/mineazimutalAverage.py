import numpy as np

def azimuthalAverage(image, center=None, binsize=3.):
   """
   Calculate the azimuthally averaged radial profile.
   image - The 2D image
   center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fractional pixels).
   
   """

   # Calculate the indices from the image
   y, x = np.indices(image.shape)
   if center is None:
      center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
   r = np.hypot(x - center[0], y - center[1])

  # the 'bins' as initially defined are lower/upper bounds for each bin
  # so that values will be in [lower,upper] 
   nbins = int(np.round(r.max() / binsize)+1)
   maxbin = nbins * binsize
   bins = np.logspace(np.log10(1.),np.log10(maxbin),nbins+1)
#   bins = np.linspace(1,maxbin,nbins+1)
  # but we're probably more interested in the bin centers than their left or right sides...
   bin_centers = (bins[1:]+bins[:-1])/2.0
  # how many per bin (i.e., histogram)?
  # there are never any in bin 0, because the lowest index returned by digitize is 1
   nr = np.histogram(r,bins)[0]
  
   radial_prof = np.histogram(r, bins, weights=image)[0] / nr
   radial_prof[radial_prof!=radial_prof] = 0
  
   indx=np.where(radial_prof==0)
   radialnew=np.delete(radial_prof, indx)
   binnew=np.delete(bin_centers, indx)
   nrnew=np.delete(nr, indx)

  
   return binnew, radialnew, np.sqrt(2./nrnew)*radialnew
  
