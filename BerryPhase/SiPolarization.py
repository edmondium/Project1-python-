#import TB model class
from pythtb import *
import matplotlib.pyplot as plt
#read output from Wannier90
silicon=w90(r"example_a", r"silicon")
#get TB model without hopping terms above 0.01 eV
my_model=silicon.model(min_hopping_norm=0.01)
my_model.ignore_position_operator_offdiagonal()
#solve model on a path and plot it
path=[[0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [0.5, -0.5, 0.0], [0.375, -0.375, 0.0], [0.0, 0.0, 0.0]]
#labels of the nodes
k_label=(r'$L$', r'$\Gamma$', r'$X$', r'$K$', r'$\Gamma$')
#call function k_path to construct the actual path
(k_vec, k_dist, k_node)=my_model.k_path(path, 101, report=True)
#evals=my_model.solve_all(k_vec)
(evals, evec)=my_model.solve_all(k_vec, eig_vectors=True)
#pick band=0 from evec[band, kpoint, orbital]
evec=evec[0] #now just evec[kpoint, orbital]
#compute Berry phase of lowest band
prod=1.+0.j
#<evec_0 | evec_1>...<evec_98 | evec_99>
for i in range(1, 101-1):
#a *= b means a = a*b
    prod *= np.vdot(evec[i-1], evec[i])
#compute the phase factors needed for last inner product
orb=np.array([0.5]) #relative coordinates of orbitals
phase=np.exp((-2.j)*np.pi*orb) #construct phase factors
evec_last=phase*evec[0] #evec[100] constructed from evec[0]
prod *= np.vdot(evec[-2], evec_last) #include <evec_99 | evec_last>
print("Berry phase is %f" %(-np.angle(prod)))
#Berry phase via the wf_array method
evec_array=wf_array(my_model, [101]) #set array dimension
#evec_array.solve_on_grid([0.]) #fill with eigensolutions
berry_phase=evec_array.berry_phase([0]) #Berry phase of bottom band
print("Berry phase is %f" %berry_phase)
#plotting of band structure
fig, ax=plt.subplots()
for i in range(evals.shape[0]):
    ax.plot(k_dist, evals[i], "k-")
for n in range(len(k_node)):
    ax.axvline(x=k_node[n], linewidth=0.5, color='k')
ax.set_xlabel("Path in k-space")
ax.set_ylabel("Band energy (eV)")
ax.set_xlim(k_dist[0], k_dist[-1])
ax.set_xticks(k_node)
ax.set_xticklabels(k_label)
fig.tight_layout()
fig.savefig("SiPolarization")