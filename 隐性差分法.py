import numpy as np

X = 50 #Strike price
S0 = 48 #spot price
N = N = input ('Please input the steps：')  
N = int (N)   #step
I = N#STEP
T = 0.5
r = 0.03
sigma = 0.3
Smax = 3*X
dt = T/N
dS = Smax/I
VGrid = np.zeros((I+1,N+1))

for i in range (N+1):
 
  VGrid[0][i] = (Smax-X)*np.exp(-r*(T-i*dt))# at S=0;
  VGrid[I][i] = 0   # at Smax, unlikely to have positive payoff  
for j in range (I+1):
  VGrid[j][N] = max(X-j*dS,0)     # at S=0;
  
  
  #tridiag
def tridiag(a, b, c):
    return np.diag(a, 0) + np.diag(b, 1) + np.diag(c, -1)   #1up -1 d   0 diag

aj = np.zeros((1,I-1))
for i in range (1,I):
  aj[0][i-1] = 1 + sigma**2*i**2*dt + r*dt
aj = aj.flatten()

bj = np.zeros((1,I-2))
for i in range (1,I-1):
  bj[0][i-1] = -0.5*(sigma**2*i**2 - r*i)*dt
bj = bj.flatten()

cj = np.zeros((1,I-2))
for i in range (1,I-1):
  cj[0][i-1] = -0.5*(sigma**2*i**2 + r*i)*dt
cj = cj.flatten()

M = tridiag(aj, bj, cj).reshape(I-1,I-1)
M1 = np.matrix(M).I

def Rsh(a):
	D = VGrid[:,a][::-1]
	return D[:][1:I].reshape(I-1,1)
	
def Fn(a):
    E = VGrid[:,a][::-1]
    E = E.flatten()
    R = np.zeros((1,I-1))
    R[0][0] = -0.5*(sigma**2 - r)*dt*E[1]
    R[0][-1] = (-0.5*(sigma**2*(I-1)**2 + r*(I-1))*dt)*E[-1]
    return  R.reshape(I-1,1)

A=Rsh(I)
print (A)
def Rshq(n):
 if n==I:
  return A
 else:
   return M1*(Rshq(n+1)-Fn(n+1))
F = int (I/2)
print (Rshq(1)[F][0])

