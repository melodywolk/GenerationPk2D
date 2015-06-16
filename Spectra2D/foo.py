import pyfits
import numpy as np
import pylab as py
import matplotlib.pyplot as plt
import convertlog
import proj
import projold
from PAstar2 import predictPAstar

zmin=0.6
zmax=0.8
Nz=50

L = 7.4652025*(np.pi/180.)
Volume1 = L**2
bias=1.27
Nboxsize=128
xibar1 = 0.14412392
Ngal1 = 9.6835435*Nboxsize**2

k, Pk,err1=np.loadtxt(open('/home/melody/Desktop/Pk/Outputs/PkW10608.txt'),unpack='True')
kl, Pkl,err1l=np.loadtxt(open('/home/melody/Desktop/Pk/Outputs/PklogW10608.txt'),unpack='True')

k1, Pk1 = np.loadtxt(open('gW10608_Nres5000'),unpack='True')
kA, PkA = np.loadtxt(open('asW10608_Nres5000'),unpack='True')
k1 = (2*np.pi/L)*k1
kA = (2*np.pi/L)*kA
Pk1=Volume1*Pk1
PkA=Volume1*PkA
zdis1=np.loadtxt(open('zdisW10210.txt'),unpack='True')

dndz1, z1 = np.histogram(zdis1, range=(zmin,zmax), bins = Nz)
Nofz1 = np.zeros((len(dndz1),2)).reshape((len(dndz1),2))
Nofz1[:,0] = np.array(z1[1:,])
Nofz1[:,1] = np.array(dndz1)*(Ngal1/np.sum(dndz1))

lA, ClA = proj.projection(Nofz1, (zmax+zmin)/2., 1.27, Nboxsize, Volume1, 200, 0, 1)
l1, Cl1, k2Dlog1, Pk2Dlog1 = predictPAstar(lA,ClA, k1, 50, xibar1,Nboxsize, Ngal1, Volume1)

plt.loglog(k1, Pk1, 'b')
plt.loglog(l1, Cl1, 'b--')
plt.loglog(kA, PkA, 'r')
plt.loglog(k2Dlog1, Pk2Dlog1, 'r--')
plt.loglog(k, Pk)
plt.loglog(kl, Pkl)
plt.show()
