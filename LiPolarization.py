#3D model of Li on bcc lattice, with s orbitals only
#python3 style print
from __future__ import print_function
#import TB model class
from pythtb import *
import matplotlib.pyplot as plt
#define function to set up model for a given parameter set
def set_model(t, del_t, Delta):
#define lattice vectors
    lat=[[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]
#define coordinates of orbitals
    orb=[[0.0, 0.0, 0.0]]
#make 3D model
    my_model=tb_model(3, 3, lat, orb)
#set on-site energy
    my_model.set_onsite([Es])
#set hoppings along 4 unique bonds
#note that neighboring cell must be specified in lattice coordinates
#the corresponding Cartesian coordinates are given for reference
    my_model.set_hop(t, 0, 0, [1, 0, 0]) #[-0.5, 0.5, 0.5] cartesian
    my_model.set_hop(t, 0, 0, [0, 1, 0]) #[0.5, -0.5, 0.5] cartesian
    my_model.set_hop(t, 0, 0, [0, 0, 1]) #[0.5, 0.5, -0.5] cartesian
    my_model.set_hop(t, 0, 0, [1, 1, 1]) #[0.5, 0.5, 0.5] cartesian
    return my_model
#set model parameters
#lattice parameter implicitly set to a=1
#site energy
Es=4.5
#hopping parameter
t=-1.4
#bond strength alternation
del_t=-0.3
#site energy alternation
Delta=0.4
my_model=set_model(t, del_t, Delta)
#print TB model
my_model.display()
#generate k-point path and labels specified in reciprocal lattice coordinates
k_P=[0.25, 0.25, 0.25] #[0.5, 0.5, 0.5] cartesian
k_Gamma=[0.0, 0.0, 0.0] #[0.0, 0.0, 0.0] cartesian
k_H=[-0.5, 0.5, 0.5] #[1.0, 0.0, 0.0] cartesian
path=[k_P, k_Gamma, k_H]
label=(r'$P$', r'$\Gamma$', r'$H$')
#set up and solve the model on a discretized k mesh
nk=101 #100 equal intervals around the unit circle
#(k_vec, k_dist, k_node)=my_model.k_path(path, 101)
(k_vec, k_dist, k_node)=my_model.k_path(path, nk, report=False)
print('--------------------')
print('Starting calculation')
print('--------------------')
print('Calculating bands...')
#solve for eigenvalue of Hamiltonian on the set of k-points from above
#evals=my_model.solve_all(k_vec)
(evals, evec)=my_model.solve_all(k_vec, eig_vectors=True)
#pick band=0 from evec[band, kpoint, orbital]
evec=evec[0] #now just evec[kpoint, orbital]
#compute Berry phase of lowest band
prod=1.+0.j
#<evec_0 | evec_1>...<evec_98 | evec_99>
for i in range(1, nk-1):
#a *= b means a = a*b    
    prod *= np.vdot(evec[i-1], evec[i])
#compute the phase factors needed for last inner product
orb=np.array([0.5]) #relative coordinates of orbitals
phase=np.exp((-2.j)*np.pi*orb) #construct phase factors
evec_last=phase*evec[0] #evec[100] constructed from evec[0]
prod *= np.vdot(evec[-2], evec_last) #include <evec_99 | evec_last>
print("Berry phase is %f" %(-np.angle(prod)))
#plotting of band structure
print('Plotting band structure...')
#make a figure object
fig, ax=plt.subplots(figsize=(4., 3.))
#specify horizontal axis details
ax.set_xlim([0, k_node[-1]])
ax.set_xticks(k_node)
ax.set_xticklabels(label)
for n in range(len(k_node)):
    ax.axvline(x=k_node[n], linewidth=0.5, color='k')
#plot bands
ax.plot(k_dist, evals[0], color='k')
#put title
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy")
#make a png figure of a plot
fig.tight_layout()
fig.savefig("LiPolarization")
print('Done.\n')