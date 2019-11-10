import numpy as np
import matplotlib.pyplot as plt

R = 2 #Resistance
L = 1 #Inductance
C = 0.1 #Capacitance
Vm = 5 #Max.Source Voltage
T_f = 20 #Final Time
i_0 = 0 #Initial Current
i_dot_0 = 100 #Initial i_dot 
w = 1 #Omega

#Define a function which returns the second derivative of current
def i_double_dot(R,L,C,Vm,w,t,i,i_dot):
    return Vm*w*np.cos(w*t)/L +(-1)*((R*i_dot/L)+(i/(L*C)))

#Implement the Euler's method

def CurrentVec (T_f,nt):
    tvec = np.linspace(0,T_f,nt+1)      #Vector of the x_values of the graph
    ivec = np.zeros(tvec.size)          #Vector of the current at corresponding time
    i_dot_vec = np.zeros(tvec.size)
    ivec[0] = i_0
    i_dot_vec[0] = i_dot_0
    for j in range(1,tvec.size):        #For each time step looping to find the curent flow
        delta_t = tvec[j] - tvec[j-1]
        Current_double_dot = i_double_dot(R,L,C,Vm,w,tvec[j-1],ivec[j-1],i_dot_vec[j-1])
        i_dot_vec[j] = Current_double_dot*delta_t + i_dot_vec[j-1]
        ivec[j] = i_dot_vec[j]*delta_t + ivec[j-1]
    return [tvec,ivec]

def ExperimentalCurrentVec (T_f,nt):
    tvec = np.linspace(0,T_f,nt+1)      #Vector of the x_values of the graph
    ivec = np.zeros(tvec.size)          #Vector of the current at corresponding time
    i_dot_vec = np.zeros(tvec.size)
    ivec[0] = i_0
    i_dot_vec[0] = i_dot_0
    for j in range(1,tvec.size):        #For each time step looping to find the curent flow
        delta_t = tvec[j] - tvec[j-1] 
        Current_double_dot = i_double_dot(R,L,C,Vm,w,tvec[j-1],ivec[j-1],i_dot_vec[j-1]) + np.random.randint(-3, 5, size=2)[0]
        i_dot_vec[j] = Current_double_dot*delta_t + i_dot_vec[j-1]
        ivec[j] = i_dot_vec[j]*delta_t + ivec[j-1]
    return [tvec,ivec]

#####################################
#           Main plot               #
#####################################

iextact_nt1000 = CurrentVec(T_f,1000)
iEuler_nt500 = ExperimentalCurrentVec(T_f,900)
iEuler_nt1000 = CurrentVec(T_f,1000)

fig, ax = plt.subplots()
fig.suptitle('Current - RLC Circuit', fontsize=20)
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Current (A)', fontsize=12)

#ax.plot(iexact[0],iexact[1],'k--',label='Exact Solution')
#ax.plot(iEuler_nt50[0],iEuler_nt50[1],'r',
#         label='Numerical Solution with 50 time steps')
#ax.plot(iEuler_nt100[0],iEuler_nt100[1],'g',
#         label='Numerical Solution with 100 time steps')
ax.plot(iEuler_nt500[0],iEuler_nt500[1],'b',
         label='Experimental current')
ax.plot(iEuler_nt1000[0],iEuler_nt1000[1],'m',
         label='Numerical Solution with 1000 time steps')  
legend = ax.legend(loc='upper center', shadow=True, fontsize='medium')
plt.show()

#####################################
#           Error bars              #
#####################################

Err_nt1000 = abs(iextact_nt1000[1] - iEuler_nt1000[1])
iextact_nt500 = iextact_nt1000[1][0:901]
Err_nt500 = abs(iextact_nt500 - iEuler_nt500[1])

fig, ax = plt.subplots()
fig.suptitle('Error - RLC Circuit', fontsize=20)
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Error (%)', fontsize=12)
plt.ylim(0,10)

ax.plot(iEuler_nt500[0],Err_nt500,'b',label='Error plot with 500 time steps')
ax.plot(iEuler_nt1000[0],Err_nt1000,'m',label='Error plot with 1000 time steps')

legend = ax.legend(loc='upper center', shadow=True, fontsize='medium')

dtVec = [iEuler_nt500[0][1]-iEuler_nt500[0][0],
         iEuler_nt1000[0][1]-iEuler_nt1000[0][0]]

ErrVec=[sum(Err_nt500)*(dtVec[0]),
        sum(Err_nt1000)*(dtVec[1])]


plt.plot(dtVec,ErrVec,'k-x') #Plot the Error with dt 
plt.show()