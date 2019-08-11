import scipy.optimize as optimize
from scipy.optimize import fsolve
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
from pylab import rcParams


# include the elastic energy implictly
def f(P):
    global E1 , E2 , T1 , a1 , b , Cph , T0 , a , gamma
    return ( 0.5 * a * P**2.0 + 0.25 * b * P**4.0 + gamma * P * Pd - P * E1 )


def equation1(p):
    global E1 , E2 , T1 , a1 , b , Cph , T0 , a , gamma , EE
    PE = (p)
    return( a * PE + b *  PE**3.0 + gamma * Pd - EE  )

# Define the figure size
rcParams['figure.figsize'] = 13 , 13
rcParams['axes.linewidth'] = 1.8

#%  %  %  strength of the coupling between defect and normal polarization
gamma = 3.0
#%  %  %  parameters of P**2.0
a1 = 1.0 
#%  %  %  curie temperature
T0 = 1.0
#%  %  %  strength of the irreversible
#alpha = 1.0

#%  %  %  initial temperature
#T1 = 0.8
#%  %  %  parameters of P**4.0
b = 1.0 / 3.0
#%  %  %  parameters of P^6
c = 1.0 / 3.0
#%  %  %  heat capacity
Cph = 15.0
#######defect dipole
#Pd = 0.5/3.0
Pd = 0.3 / 6.0
#%  %  %  critical field
Ecp = 6.0 * b**2.0 *  np.sqrt( 3.0 * abs(b) / 10.0 / c ) / 25.0 / c
#%  %  %  initial temperature
Ti = np.linspace(0.1,2.0,50)
#Ti = np.array([ 0.6 , 0.8 , 1.0] )
#Ti = np.array([ 1.1 , 1.2 ])
## electric field
#E = np.linspace( -10.0 * Ecp , 10.0 * Ecp , 42)
E = np.array( [ 0.2 , 0.0] )
# initialize the polarization
P1 = np.zeros( ( len(E) , len(Ti) ) )
P2 = np.zeros( ( len(E) , len(Ti) ) )
EE = np.zeros( ( len(E) , len(Ti) ) )
TT1 = np.zeros( ( len(E) , len(Ti) ) )
print np.shape(Ti)
print Ti
# find the minimum value
i = 0
### save the data



for E1 in E:  
  j = 0
  # %  %  %  parameters of P**2.0
  for T1 in Ti :
    a = a1 * (T1 - T0) 
    P1[i,j] = optimize.fmin( f , [  1.0 ] )
    P2[i,j] = optimize.fmin( f , [  -1.0 ] )
    EE[i,j] = E1
    TT1[i,j] = T1
    #dataset = np.array( [ P1[i,j] , P2[i,j] , EE[i,j] , TT[i,j] ] )
    #print dataset    
    #np.savetxt(f,a)
    j = j + 1 
  
  i = i + 1

# initialize the size of the data  
j = 0  
Tend = np.zeros(len(Ti))
Wloss = np.zeros(len(Ti))
Stotal = np.zeros(len(Ti))
Sconf = np.zeros(len(Ti))
Svib = np.zeros(len(Ti))
jend = len(Ti)
Eend = E1 * np.ones(len(Ti))


for T1 in Ti :
  if T1 < T0 :
    # strength of the irreversible process
    #alpha = 2.0 / ( T0 - T1 )
    alpha = 1.0
    #alpha = 1.0 / ( T0 - T1 )
    # prefactor of P2
    a = a1 * (T1 - T0) 
    PF = - np.sqrt(-a/(3*b))          # negative polarization 
    PD = np.sqrt(-a/(3*b))            # positve polarization 
    EF = a * PF + b * PF**3.0 + gamma * Pd   # positive field 
    EB = EF
    ED = a * PD + b * PD**3.0 + gamma * Pd   # negative field
    EE = ED
    PE = fsolve ( equation1, -1 ) 
    PB = -PE   
    #print 'PE=',PE,'EE=',EE,T1
    # judge when P1 and P2 have opposite sign, choose P1[1] as start,P1[2] as end        
    if P1[1,j] * P2[1,j] < 0 : # S-shape has multi-solution
      PN = P1[1,j]
      PC = P1[0,j]
      print 'case1',T1
      PMC = ( a * PC + b * PC**3.0 + 2 * b * PB**3.0  ) / ( a + 3 * b * PB**2.0 )
      Wloss1_CB =  alpha * ( 0.5 * a * PC**2.0 + 0.25 * b * PC**4.0  - 0.5 * a * PMC**2.0 - 1.5 * b * PB**2.0 * PMC**2.0 + 2 * b * PB**3.0 * PMC - 0.75 * b * PB**4.0 + gamma * Pd * ( PC - PMC ) ) # the starting point to point B
      print 'Wloss1_CB=',Wloss1_CB, 'T1=',T1
      if ( E[0] < EB ) : # the field at the end point is lower than that at point B
        PM2 = ( a * PN + b * PN**3.0 + 2 * b * PB**3.0  ) / ( a + 3 * b * PB**2.0 )
        Wloss1_NB =   alpha * ( 0.5 * a * PN**2.0 + 0.25 * b * PN**4.0  - 0.5 * a * PM2**2.0 - 1.5 * b * PB**2.0 * PM2**2.0 + 2 * b * PB**3.0 * PM2 - 0.75 * b * PB**4.0 + gamma * Pd * ( PN - PM2 ) ) # the end point to point B
        Wloss1_CN =  Wloss1_CB - Wloss1_NB # the starting point to the end point
        Wloss[j] = Wloss1_CN
        Stotal[j] = Wloss[j] / TT1[0,j]
        Sconf[j] = - 0.5 * a1 * PN**2.0 + 0.5 * a1 * PC **2.0
        Svib[j] = Stotal[j] - Sconf[j]
        Tend[j] = np.exp( Svib[j] / Cph ) * TT1[0,j] # irreversible 
        #Tend[j] = np.exp( ( 0.5 * a1 * P1[1,j]**2.0 - 0.5 * a1 * P1[0,j] **2.0 ) / Cph ) * TT1[0,j] # test:reversible process
        print 'case1.1','Wloss1_CB=',Wloss1_CB,'Wloss1_NB=',Wloss1_NB,'Wloss1_CN=',Wloss1_CN,'Tend[j]=',Tend[j],'P1[1,j]=',P1[1,j],'P1[0,j]=',P1[0,j]
      else : ## the field at the end point is lower than that at point B
        PM2 = ( a * PN + b * PN**3.0 + 2 * b * PB**3.0  ) / ( a + 3 * b * PB**2.0 )
        Wloss1_NB =   alpha * ( 0.5 * a * PN**2.0 + 0.25 * b * PN**4.0  - 0.5 * a * PM2**2.0 - 1.5 * b * PB**2.0 * PM2**2.0 + 2 * b * PB**3.0 * PM2 - 0.75 * b * PB**4.0 + gamma * Pd * ( PN - PM2 ) ) # the end point to point B
        Wloss[j] = Wloss1_NB
        Stotal[j] = Wloss[j] / TT1[0,j]
        Sconf[j] = - 0.5 * a1 * PN**2.0 + 0.5 * a1 * PC **2.0
        Svib[j] = Stotal[j] - Sconf[j]
        Tend[j] = np.exp( Svib[j] / Cph ) * TT1[0,j] # irreversible 
        #Tend[j] = np.exp( ( 0.5 * a1 * P1[1,j]**2.0 - 0.5 * a1 * P1[0,j] **2.0 ) / Cph ) * TT1[0,j] # test:reversible process
    elif P1[1,j] * P2[1,j] > 0 : # the S-shape function has one solution
      PN =  P2[1,j]
      print 'case2'
      if ( E[0] < EE ) : # the field at end point is smaller than EE
        print 'case2.1'
        Pstart = P2[0,j]
        # change of the configurational entropy
        dS_conf = - 0.5 * a1 * PN**2.0 + 0.5 * a1 * Pstart **2.0     
        Sconf[j] = dS_conf
        Svib[j] = Stotal[j] - Sconf[j]
        Tend[j] = np.exp( Svib[j] / Cph ) * TT1[0,j]
      elif EE < E[0] and  E[0] < EF : # the field at end point is between EE and EF
        print 'case2.2'
        Pstart = P2[0,j]
        # change of the configurational entropy
        dS_conf = - 0.5 * a1 * PN**2.0 + 0.5 * a1 * Pstart **2.0     
        Sconf[j] = dS_conf
        # change the Pstart
        PM1 = ( a * Pstart + b * Pstart**3.0 + 2 * b * PE**3.0  ) / ( a + 3 * b * PE**2.0 )
        Wloss1_NB =   alpha * ( 0.5 * a * Pstart**2.0 + 0.25 * b * Pstart**4.0  - 0.5 * a * PM1**2.0 - 1.5 * b * PE**2.0 * PM1**2.0 + 2 * b * PE**3.0 * PM1 - 0.75 * b * PE**4.0 + gamma * Pd * ( Pstart - PM1 ) )
        Wloss[j] = Wloss1_NB
        Stotal[j] = Wloss[j] / TT1[0,j]
        Svib[j] = Stotal[j] - Sconf[j]
        Tend[j] = np.exp( Svib[j] / Cph ) * TT1[0,j] # irreversible 
        print 'W_nb', Wloss1_NB
      elif EF < E[0] : # the field at end point is bigger than EF
        print 'case2.3'
        Pstart = P1[0,j]
        print 'Pstart=',Pstart,'T1=',T1
        # change of the configurational entropy
        dS_conf = - 0.5 * a1 * PN**2.0 + 0.5 * a1 * Pstart **2.0     
        Sconf[j] = dS_conf
        PM = ( a * PD + b * PD**3.0 + 2 * b * PB**3.0  ) / ( a + 3 * b * PB**2.0 )
        Wloss1 = alpha * ( 0.5 * a * PD**2.0 + 0.25 * b * PD**4.0  - 0.5 * a * PM**2.0 - 1.5 * b * PB**2.0 * PM**2.0 + 2 * b * PB**3.0 * PM - 0.75 * b * PB**4.0 + gamma * Pd * ( PD - PM ) )
        Wloss2 = (PE - PD) * EE
        print 'Wloss1=',Wloss1,'Wloss2=',Wloss2
        #PM2 = ( a * PN + b * PN**3.0 + 2 * b * PE**3.0  ) / ( a + 3 * b * PE**2.0 )
        print ED
        Wloss[j] = Wloss1 + Wloss2
        Stotal[j] = Wloss[j] / TT1[0,j]
        Svib[j] = Stotal[j] - Sconf[j]
        Tend[j] = np.exp( Svib[j] / Cph ) * TT1[0,j] # irreversible 
  else :
    print 'case3',T1
    Sconf[j] = -0.5 * a1 * P2[1,j]**2.0 + 0.5 * a1 * P2[0,j]**2.0
    Svib[j] = Stotal[j] - Sconf[j]
    Tend[j] = np.exp( Svib[j] / Cph ) * TT1[0,j] # irreversible 
  #print EB, ED, EE, EF , E[1]
  if j< len(Ti)-2 :
    j = j + 1




### save the data
np.savetxt('d'+str(E[0])+'.txt', np.transpose( np.array([Ti,Tend,Eend]) ), header = 'Ts,Tend,Eend', fmt='%f  %f  %f'  )
np.savetxt('S'+str(E[1])+'.txt', np.transpose( np.array([Wloss,Stotal,Sconf]) ), header = 'Wloss,Stotal,Sconf', fmt='%f  %f  %f'  )


# plot the figure
# Plot 12 figures
fig, ( (ax1, ax2) ) = plt.subplots(2, 1, sharex='col', sharey='row'  )


line1 = ax1.plot( Ti[0:(jend-2)] , Tend[0:(jend-2)]-Ti[0:(jend-2)], color='red' , marker='d'  , markersize = 16 , markevery = 2,  linewidth = 3.5 )
ax1.set_ylabel( 'Temperature change' , fontsize = 35)
ax1.tick_params( width = 3 , size = 8 )

ax1.text( 1.4, -0.002 , r'$ E_\mathrm{init} = 0.2$' , fontsize = 30 )
ax1.text( 1.4, -0.006 , r'$ E_\mathrm{end} = 0.0$', fontsize = 30  )
ax1.text( 1.4, -0.01 , '$ E_\mathrm {d} = -0.3$', fontsize = 30  )

#ss = '$ E_\mathrm {d} = $' + str(Pd*3.0*2.0)
#ss1 = Pd*3.0*2.0
#ax1.text(  1.4, -0.008 , ss , fontsize = 30  )

#, label = '$E=0.32 > E_\mathrm{d}=0.3$'
#ax1.legend( frameon=False , loc='lower center', bbox_to_anchor=(0.3, 0.78),  fancybox=True, shadow=True, ncol=4 , prop={'size':30})  
## Define shared x ticks

ax1.tick_params(axis='both', which='major', labelsize = 35)
plt.setp( (ax1) , yticks = [-0.01,0.0,0.01,0.02,0.03,0.04] , ylim = [-0.014 , 0.025] )

# field on case
ax1.text(1.6, 0.018 , 'Field off',bbox={'alpha':0.5, 'pad':10} , fontsize = 30)

Tprime = T0 - (gamma * Pd)**(2.0/3.0) * (27.0*b/4.0)**(1.0/3.0)
print(Tprime)
Tdprime = T0 - (-gamma * Pd+E[0])**(2.0/3.0) * (27.0*b/4.0)**(1.0/3.0)

ax1.plot( [Tprime, Tprime] , [0.02 , -0.03] , color='gray' , linestyle = '--' , linewidth = 2.0 )

ax1.plot( [Tdprime, Tdprime] , [0.02 , -0.03] , color='gray' , linestyle = '--' , linewidth = 2.0 )

ax1.plot( [1.0, 1.0 ] , [0.02 , -0.03] , color='gray' , linestyle = '--' , linewidth = 2.0  )

ax1.text( Tprime-0.03, 0.0215 , r"$ T' $", fontsize = 30  )

ax1.text( Tdprime-0.03, 0.0215 , r"$ T'' $", fontsize = 30  )
ax1.text( 0.98 , 0.0215 , r'$ T_0 $', fontsize = 30  )





















line1 = ax2.plot( Ti[0:(jend-2)] , Wloss[0:(jend-2)], color='green' , marker='<' , label = '$W_\mathrm{loss}$' , markersize = 16 , markevery = 4 , linewidth = 3.5 )
line1 = ax2.plot( Ti[2:(jend-2)] , Stotal[2:(jend-2)], color='blue' , marker='s' , label = '$\Delta S_\mathrm{total}$' , markersize = 16 , markevery = 4 , linewidth = 3.5  )
### only lines
line1 = ax2.plot( Ti[0:(jend-2)] , Stotal[0:(jend-2)], color='blue' , linewidth = 3.5  )


ylim = [-0.41,0.42]

line1 = ax2.plot( Ti[0:(jend-2)] , Sconf[0:(jend-2)], color='black' , marker='o' , label = '$\Delta S_\mathrm{dip}$' , markersize = 16 , linewidth = 3.5  )
line2 = ax2.plot( Ti[0:(jend-2)] , Svib[0:(jend-2)], color='red' , marker='d' , label = '$\Delta S_\mathrm{vib}$' , markersize = 16 , linewidth = 3.5 )

plt.setp( (line1,line2) , markevery = 2 ) 

plt.setp( ax2 , ylim = ylim )


ax2.plot( [Tprime, Tprime] , ylim , color='gray' , linestyle = '--' , linewidth = 2.0 )

ax2.plot( [Tdprime, Tdprime] , ylim , color='gray' , linestyle = '--' , linewidth = 2.0 )

ax2.plot( [1.0, 1.0 ] , ylim , color='gray' , linestyle = '--' , linewidth = 2.0  )




ax2.set_xlabel( 'Temperature' , fontsize = 35)
ax2.set_ylabel( 'Work and entropy' , fontsize = 35)
ax2.tick_params( width = 3 , size = 8 )
ax2.legend(frameon=False ,loc='lower center', bbox_to_anchor=(0.75, 0.63),  fancybox=True, shadow=True, ncol=2 ,prop={'size':30}, columnspacing=0.5 )  


## Define shared x ticks
ax2.tick_params(axis='both', which='major', labelsize = 35 )
suba = plt.setp( (ax2) , xticks = [ 0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0 ]  )


# space between the figures
plt.subplots_adjust(hspace=0.03)
plt.subplots_adjust(wspace=0.02)
## save the figure
plt.savefig('T-S.pdf',bbox_inches='tight', transparent=True)  

plt.close()


