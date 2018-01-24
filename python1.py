from qutip import *
from pylab import *
from scipy import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import time

# Count program run time
tstart = int(time.time())

# Define parameters
Nd =20            #number of basis states of phonon 1
omega1 = 1          #frequency of boson mode 1
omegaq=1            #frequency of qubit
g1=0.02              #coupling strength mode 1

# Define indintity for mode 1 and qubit
idaq = qeye(2);
idav = qeye(Nd);

# Define vibration operators
b1 = tensor(destroy(Nd),idaq);
# Define down operators for qubit
sm= tensor(idav,sigmam());
# Define spin operators.
sz=sm.dag()*sm-sm*sm.dag();

# Collapse operators decay for the mode 1, qubit and dephasing for the qubit. if none, leave the Collapse operators in mesolve empty with
kappa1=0.01
gammad=0.01
gammaf=0.001
C1 =sqrt(kappa1)*b1;
Dd1 = sqrt(gammad)*sm;
Ddf = sqrt(gammaf)*sz;

# Time list
t0 = 0.0
t1 = 100
Nt = 1000
times = linspace(t0,t1,Nt+1)        # 60 is very long time for the system

# Initial state
psi0 = tensor(basis(Nd,1),basis(2,1));

# Hamiltonian
H =0.5*omegaq*sz+omega1*b1.dag()*b1+g1*(b1.dag()+b1)*(sm+sm.dag());

# master equation slover
data = mesolve(H, psi0, times,[C1,Dd1,Ddf], [sm.dag()*sm])
#obtain the expection values
popexc=data.expect[0]

# save them in the .txt
data = open("data1.txt","w")     # Open(or create) a .txt document named by medata, as "write" style
for j in range(1000) :
    # Translate the data to "string" type
    tt   = '%f' %times[j]         # tt is short for t_txt, which represents txt data for times
    p1t  = '%f' %popexc[j]            # p1t is short for p1_txt
    data.write(tt+'	'+p1t+'	\n')    # write in the txt
data.close()

#plot the figure
fig, ax = plt.subplots(figsize=(8,6))
ax.plot(times, popexc)
ax.set_xlabel('Time')
ax.set_ylabel('P_ex')

ax.set_title('Excitation probabilty of qubit');
#show()

###########################################wigner fun#########################################################################

# obtain steady state
rho = steadystate(H,[C1,Dd1,Ddf])
#Calculating the time-dependent wignerfuction
xvec = linspace(-3,3,100)

ro = ptrace(rho, 0) #partial trace of density matrix
W_coherent = wigner(ro, xvec, xvec)
# another data save method----with python
file_data_store('wig.dat', W_coherent)

fig, axe =plt.subplots(figsize=(6,6))
axe.contourf(xvec, xvec, W_coherent, 100)
show()
###########################################wigner fun#########################################################################


# run time
tend = int(time.time())
tspent = tend - tstart
print('Total run time:',tspent,'\n')


