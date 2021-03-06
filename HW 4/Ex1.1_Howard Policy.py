import numpy as np
from timeit import default_timer as timer


# Paremeters
b=0.988   #beta
th=.679 #theta
d=.013  #delta
kmax = 5.24 #kappa
sig = 2 # consumption elasticity 
h=1

# Steady state capital
ks=(1/(1-th)*((1/b)+d-1))**(-1/th) 

# k grid 
k=np.linspace(0.01,2*ks,10) #close to steady state
x,y=np.meshgrid(k,k)

# Value function
V = np.empty(shape=[10, 45])
V[:,0]=np.zeros((10))  # initial guess of 0 

# Defining return matrix
def return_mat(k1,k2):
    return k1**(1-th) + (1-d)*k1 - k2
N = return_mat(x,y)

# Utility function
def utility(c1,c2):
    for i in range(40):
        for j in range (10):
                if N[i,j]>=0:
                    return np.log(c1**(1-th)*h**(th) + (1-d)*c1 - c2) - kmax*((h**(1+(1/sig)))/(1+(1/sig)))
            
M = utility(x,y)
M[np.isnan(M)] = -10


## part f: Howard policy 
 
# Starting time
start= timer()

# Value function iterations
X = np.empty(shape=[10, 40])
P= np.empty(shape=[10, 45])
iter=0 # starting iteration
vals=1 
err=0.01
while vals>err:
   for n in range(0,40):
       V_old=np.copy(V) # storing old value function
       if iter<100:   # we begin at 100 iterations
           for i in range(10):
               for j in range(10):
                   X[i,j]=M[i,j]+(b*V[:,n][j])
           for i in range(0,10):
               V[:,n+1][i]= np.amax(X[:,i]) 
           vals= abs(np.amax(V-V_old))    
           iter +=1
   else:
       for i in range(10):
               for j in range(10):
                   X[i,j]=M[i,j]+(b*V[:,n][j])
       for i in range(0,10):
               V[:,n+1][i]= np.amax(X[:,i]) 
               P[:,n][i]=np.argmax(X[:,i]) 
       vals= abs(np.amax(V-V_old))    
       iter +=1
        
          
           
     

#end time
end= timer()
iter_time= end - start
print('Howard policy iterations time', iter_time)
print('Howard policy iteration',iter) 


# part g: Different policy iterations

# 5 steps

# Starting time
start= timer()

# Value function iterations
X = np.empty(shape=[10, 40])
P= np.empty(shape=[10, 45])
iter=0 # starting iteration
vals=1 
err=0.01
while vals>err:
   for n in range(0,40):
       V_old=np.copy(V) # storing old value function
       for i in range(10):
           for j in range(10):
               X[i,j]=M[i,j]+(b*V[:,n][j])
       for i in range(0,10):
           V[:,n+1][i]= np.amax(X[:,i])
           if iter%5==0:
               P[:,n][i]=np.argmax(X[:,i]) 
       vals= abs(np.amax(V-V_old))    
       iter=iter+1
        
          
           
     

#end time
end= timer()
iter_time= end - start
print('Howard policy iterations time: 5 steps', iter_time)
print('Howard policy iteration: 5 steps ',iter) 

# 10 steps

# Starting time
start= timer()

# Value function iterations
X = np.empty(shape=[10, 40])
P= np.empty(shape=[10, 45])
iter=0 # starting iteration
vals=1 
err=0.01
while vals>err:
   for n in range(0,40):
       V_old=np.copy(V) # storing old value function
       for i in range(10):
           for j in range(10):
               X[i,j]=M[i,j]+(b*V[:,n][j])
       for i in range(0,10):
           V[:,n+1][i]= np.amax(X[:,i])
           if iter%10==0:
               P[:,n][i]=np.argmax(X[:,i]) 
       vals= abs(np.amax(V-V_old))    
       iter=iter+1
        
          
           
     

#end time
end= timer()
iter_time= end - start
print('Howard policy iterations time: 10 steps', iter_time)
print('Howard policy iteration: 10 steps ',iter)

# 20 steps

# Starting time
start= timer()

# Value function iterations
X = np.empty(shape=[10, 40])
P= np.empty(shape=[10, 45])
iter=0 # starting iteration
vals=1 
err=0.01
while vals>err:
   for n in range(0,40):
       V_old=np.copy(V) # storing old value function
       for i in range(10):
           for j in range(10):
               X[i,j]=M[i,j]+(b*V[:,n][j])
       for i in range(0,10):
           V[:,n+1][i]= np.amax(X[:,i])
           if iter%20==0:
               P[:,n][i]=np.argmax(X[:,i]) 
       vals= abs(np.amax(V-V_old))    
       iter=iter+1
        
          
           
     

#end time
end= timer()
iter_time= end - start
print('Howard policy iterations time: 20 steps', iter_time)
print('Howard policy iteration: 20 steps ',iter)

# 50 steps

# Starting time
start= timer()

# Value function iterations
X = np.empty(shape=[10, 40])
P= np.empty(shape=[10, 45])
iter=0 # starting iteration
vals=1 
err=0.01
while vals>err:
   for n in range(0,40):
       V_old=np.copy(V) # storing old value function
       for i in range(10):
           for j in range(10):
               X[i,j]=M[i,j]+(b*V[:,n][j])
       for i in range(0,10):
           V[:,n+1][i]= np.amax(X[:,i])
           if iter%50==0:
               P[:,n][i]=np.argmax(X[:,i]) 
       vals= abs(np.amax(V-V_old))    
       iter=iter+1
        
          
           
     

#end time
end= timer()
iter_time= end - start
print('Howard policy iterations time: 50 steps', iter_time)
print('Howard policy iteration: 50 steps ',iter)
