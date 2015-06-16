import numpy as np


def gfield2d(N,k,Pk,V):
 # Returns a N times N Gaussian field with specified spectrum (k,Pk)
 # k and Pk are physical ; quantities in 1/length
 # V is the physical volume of the box, in length^2

 # Generates real and imaginary parts of the Fourier transform of the field and transform back
 # Uses periodic boundary conditions

 # Not optimized at all.
 # JCarron 16.12.13

 # Forces even number of points for convenience :
 if (N % 2)!=0:
    N += 1
 N2 = N/2
 L = np.sqrt(V)
 # values of k for a 1d box, k_i = 2pi*i*/L :
 value=np.arange(N2-1)+1
 kid = np.concatenate([1.*np.arange(N2+1),-value[::-1]])*2.*np.pi/L
 # variance of the re and im parts of modes
 Spk2 = np.sqrt(Pk/2.*V) 
 #;---------------------
 #; generates Fourier modes 
 #;------------------------
 #; First we build a map of the magnitude of the wavevector :
 normmap = np.zeros((N,N)).reshape((N,N)) # The norm of k at each point 
 for i in range(0, N):
     for j in range(0, N):
         normmap[i,j] = np.sqrt(kid[j]*kid[j] + kid[i]*kid[i])
 # Now we generate Gaussin variables with zero mean and unit variances :
 re = np.random.standard_normal([N,N]) # real part
 im = np.random.standard_normal([N,N]) # imaginary part
 re[0,0] = 0. # we force zero power at zero, i.e. zero mean field
 im[0,0] = 0. # idem
 im[N2,0] = 0.
 re[N2,0]  *= np.sqrt(2.) # this one is real, so the variance if VP and not VP/2
 im[0,N2] = 0.  
 re[0,N2]  *= np.sqrt(2.) 
 im[N2,N2] = 0.
 re[N2,N2] *= np.sqrt(2.)

 # We have to enforce the reality conditions, f(-k) = f^*(k)
 # I think the lines below are correct but I am not even sure :
 re[N2+1:N-1,0] =  np.flipud(re[1:N2-1,0])
 im[N2+1:N-1,0] = -np.flipud(im[1:N2-1,0])
 re[0,N2+1:N-1] =  np.flipud(re[0,1:N2-1])
 im[0,N2+1:N-1] = -np.flipud(im[0,1:N2-1])
 
 for i in range(1,N2-1):
     re[N-i,N2+1:N-1] =  np.flipud(re[i,1:N2-1])
     im[N-i,N2+1:N-1] = -np.flipud(im[i,1:N2-1])
     re[i,N2+1:N-1]   =  np.flipud(re[N-i,1:N2-1])
     im[i,N2+1:N-1]   = -np.flipud(im[N-i,1:N2-1])
  
 re[N2,N2+1:N-1] =  np.flipud(re[N2,1:N2-1]) 
 im[N2,N2+1:N-1] = -np.flipud(im[N2,1:N2-1])
 re[N2+1:N-1,N2] =  np.flipud(re[1:N2-1,N2])
 im[N2+1:N-1,N2] = -np.flipud(im[1:N2-1,N2])
 #print re[:,0]
 # We now build the map of variance to mulitply the Gaussian variables with.
 # If the magnitude k does not fall into one bin provided in input then the power is   set to zero

 facmap = np.zeros((N,N)).reshape((N,N))
 ii = (normmap>=np.min(k)) & (normmap<np.max(k))
 if (len(normmap[ii])>0): 
    facmap[ii] = np.exp(np.interp(np.log(normmap[ii]), np.log(k), np.log(Spk2)))
 re *= facmap
 im *= facmap
 # Perform FFT back, the output should be real 
 tab=re + im*complex(0,1)
 #print np.fft.fft2(tab)
 print np.fft.fft2(tab).imag
 #print np.abs(np.fft.fft2(tab))
 return np.fft.fft2(tab).real/V

