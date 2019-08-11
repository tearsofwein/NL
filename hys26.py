import scipy.optimize as optimize
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
from pylab import rcParams

from scipy.optimize import fsolve

# include the elastic energy implictly
def f(P):
    global E1 , E2 , T1 , a1 , b , Cph , T0 , a , gamma
    return ( 0.5 * a * P**2.0 + 0.25 * b * P**4.0 + gamma * P * Pd - P * E1 )

def equation1(p):
    global E1 , E2 , T1 , a1 , b , Cph , T0 , a , gamma , EE
    PE = (p)
    return( a * PE + b *  PE**3.0 + gamma * Pd - EE  )

#%  %  %  initial temperature


# Define the figure size
rcParams['figure.figsize'] = 6, 6
rcParams['axes.linewidth'] = 1.8

#TT1 = np.linspace(0.1,2.0,20)
#TT1 = np.array([0.5 ])
#TT1 = np.array([0.5 , 0.8 , 1.0 , 1.8])

TT1 = np.array([0.7 , 0.8 , 0.9 , 1.0])
#print TT1
#%  %  %  strength of the coupling between defect and normal polarization
gamma = 3.0
#%  %  %  parameters of P^2
a1 = 1.0 
#%  %  %  curie temperature
T0 = 1.0
#%  %  %  initial temperature
#T1 = 0.8
#%  %  %  parameters of P^4
b = 1.0 / 3.0
#%  %  %  parameters of P^6
c = 1.0 / 3.0
#%  %  %  heat capacity
Cph = 15.0
#######defect dipole
#Pd = 0.01
#Pd = 1.0 / 3.0 * 0.1
Pd = 0.3 / 6.0
#%  %  %  critical field
Ecp = 6.0 * b**2.0 *  np.sqrt( 3.0 * abs(b) / 10.0 / c ) / 25.0 / c


## electric field
E = np.linspace( -0.05 , 0.6 , 500)
## applied electric field 
E_a = np.linspace( 0.2 , 0.0  , 500)
#E = np.array( [0.0,0.4] )
# initialize the polarization
P1 = 0.0 * E
P2 = 0.0 * E
# initialize the polarization (applied)
P1_a = 0.0 * E_a
P2_a = 0.0 * E_a

j = 0








for T1 in TT1 :
  fig = plt.figure()
  ax1 = fig.add_subplot(111)
  a = a1 * (T1 - T0) 
  PF = - np.sqrt(-a/(3*b))          # negative polarization 
  PD = np.sqrt(-a/(3*b))            # positve polarization 
  EF = a * PF + b * PF**3.0 + gamma * Pd   # positive field 
  EB = EF
  ED = a * PD + b * PD**3.0 + gamma * Pd   # negative field
  EE = ED
  PE = fsolve ( equation1, -1 ) 
  PB = -PE   
  # find the minimum value
  i = 0
  for E1 in E:
    P1[i] = optimize.fmin( f , [  1.0 ] )
    P2[i] = optimize.fmin( f , [  -1.0 ] )
    i = i + 1 
  # find the minimum value
  i = 0  
  for E1 in E_a:
    P1_a[i] = optimize.fmin( f , [  1.0 ] )
    P2_a[i] = optimize.fmin( f , [  -1.0 ] )
    i = i + 1
  delta_E =  E_a[1:]-E_a[:-1]
  delta_P1 = P1_a[1:]-P1_a[:-1]
  delta_P2 = P2_a[1:]-P2_a[:-1]    
  delta_P1_n = delta_P1 / np.sqrt(delta_E**2.0 + delta_P1**2.0)
  delta_P2_n = delta_P2 / np.sqrt(delta_E**2.0 + delta_P2**2.0)
  print(delta_E.shape, delta_P1.shape, delta_P2.shape)
  # +++++++++++++++++++++++++++++++++++++++
  # the turning point
  TA1 = T0 - ( gamma * Pd)**(2.0/3.0) * (27.0 * b / 4.0 )**(1.0/3.0)
    ## plot the hysteresis
  line1 = plt.plot( E , P1 , color = 'blue' ,linewidth=3.0 , zorder = 1)
    ## plot the hysteresis
  line2 = plt.plot( E , P2 , color = 'blue' ,linewidth=3.0 , zorder = 1)
  
  # +++++++++++++++++++++++++++++++++++++++ 
  # temperature below TA1
  if T1 <= TA1 :
    delta_E_n = delta_E / np.sqrt(delta_E**2.0 + delta_P1**2.0)
    ax1.quiver(E_a[::250], P1_a[::250], delta_E_n[::250], color='red', angles='xy', edgecolors=('k'), headaxislength=5 , zorder = 2 , scale = 10.0 , headwidth = 3, headlength = 5 , linewidths=(2.0,) , width = 0.015 )
    ax1.plot( min(E_a) , min(P1_a) , marker = 'o' , markersize = 18 , color = 'lime' , zorder = 2 )
    ax1.plot( max(E_a) , max(P1_a) , marker = 'o' , markersize = 18 , color = 'black' , zorder = 2 )    
    ax1.fill_between( E_a , P1_a , max(P1_a) , facecolor='gray' , edgecolor = 'none' )
    
    # +++++++++++++++++++++++++++++++++++++++
    # line of reversible polarization
    P_ir = np.linspace( PD , PB, 300 )
    Pr1 = ( a * P_ir + b* P_ir**3.0 + 2.0 * b * PB**3.0) / ( a + 3.0 * b * PB**2.0 )
    #print a + 3.0 * b * PB**2.0, PB
    # reversible electric field
    Er1 = a * Pr1 + 3.0 * b * PB**2.0 * Pr1 - 2.0 * b * PB**3.0 + gamma * Pd
    # plot line of reversible polarization
    liner = plt.plot( Er1 , Pr1 ,  color = 'gray' ,linewidth=2.0 , zorder = 1)
    #yy2 = [0.0 , 0.0 , min(Pr1) , min(Pr1) ]  
  
    # +++++++++++++++++++++++++++++++++++++++
    # area reversible polarization
    Pr = ( a * P1_a + b* P1_a**3.0 + 2.0 * b * PB**3.0) / ( a + 3.0 * b * PB**2.0 )
    #print a + 3.0 * b * PB**2.0, PB
    # reversible electric field
    Er = a * Pr + 3.0 * b * PB**2.0 * Pr - 2.0 * b * PB**3.0 + gamma * Pd
    ##print 'PE=',PE,'max(P1_a',max(P1_a)
    
    # plot the reversible area
    ax1.fill_between( Er , Pr, max(Pr),  facecolor='green' , edgecolor = 'none' )
    
  # +++++++++++++++++++++++++++++++++++++++
  # temperature above T0
  if T1 >= T0 :
    delta_E_n = delta_E / np.sqrt(delta_E**2.0 + delta_P1**2.0)
    ax1.quiver(E_a[::250], P1_a[::250], delta_E_n[::250], delta_P1_n[::250], color='red', angles='xy', edgecolors=('k'), headaxislength=5 , zorder = 2 , scale = 10.0 , headwidth = 3, headlength = 5 , linewidths=(2.0,) , width = 0.015 )
    ax1.plot( min(E_a) , min(P1_a) , marker = 'o' , markersize = 18 , color = 'lime' , zorder = 2 )
    ax1.plot( max(E_a) , max(P1_a) , marker = 'o' , markersize = 18 , color = 'black' , zorder = 2 )        
    ax1.fill_between( E_a , P1_a , max(P1_a) , facecolor='green' , edgecolor = 'none')
    
  ## +++++++++++++++++++++++++++++++++++++++     
  # both the reversible and irresible process exist, along the bottom line
  if (max(E_a)<=EF) & (max(E_a)>=EE) & (T1>TA1) & (T1<T0)  :  
     delta_E_n = delta_E / np.sqrt(delta_E**2.0 + delta_P2**2.0)     
     ax1.quiver(E_a[::250], P2_a[::250], delta_E_n[::250], delta_P2_n[::250], color='red', angles='xy', edgecolors=('k'), headaxislength=5 , zorder = 2 , scale = 10.0 , headwidth = 3, headlength = 5 , linewidths=(2.0,) , width = 0.015)
     ax1.plot( min(E_a) , min(P2_a) , marker = 'o' , markersize = 18 , color = 'lime' , zorder = 2 )
     ax1.plot( max(E_a) , max(P2_a) , marker = 'o' , markersize = 18 , color = 'black' , zorder = 2 )         
     #print 'TA1=',TA1
     # +++++++++++++++++++++++++++++++++++++++
     # line of reversible polarization
     P_ir = np.linspace( PE , PF, 300 )
     Pr1 = ( a * P_ir + b* P_ir**3.0 + 2.0 * b * PE**3.0) / ( a + 3.0 * b * PE**2.0 )
     #print a + 3.0 * b * PE**2.0, PE
     # reversible electric field
     Er1 = a * Pr1 + 3.0 * b * PE**2.0 * Pr1 - 2.0 * b * PE**3.0 + gamma * Pd
     # plot line of reversible polarization
     liner = plt.plot( Er1 , Pr1 ,  color = 'gray' ,linewidth=2.0 , zorder = 1)
  
     # +++++++++++++++++++++++++++++++++++++++
     # area reversible polarization
     Pr = ( a * P2_a + b* P2_a**3.0 + 2.0 * b * PE**3.0) / ( a + 3.0 * b * PE**2.0 )
     #print a + 3.0 * b * PE**2.0, PE
     # reversible electric field
     Er = a * Pr + 3.0 * b * PE**2.0 * Pr - 2.0 * b * PE**3.0 + gamma * Pd
     
     ## the whole area
     ax1.fill_between( E_a , P2_a , max(P2_a) , facecolor='gray', interpolate=True  , edgecolor = 'none')
     #the reversible area out of the hysteresis
     ax1.fill_between( E_a , P2_a , max(PE) , where = ( (P2_a<=PE)  ) , facecolor='green', interpolate=True , edgecolor = 'none')
     # # the reversible area within the hysteresis
     ax1.fill_between( Er1 , Pr1 , max(Pr) , where = ( (Pr1>=PE) & (Er1<=max(E_a)) ) ,facecolor='green', interpolate=True , edgecolor = 'none')
     
     # rectangle due to the tangent line of hysteresis
     xx = [0.0 , EE  , EE , 0.0 ]
     # area 
     Prr = ( a * PF + b* PF**3.0 + 2.0 * b * PE**3.0) / ( a + 3.0 * b * PE**2.0 )     
     yy1 = [ PE, PE , max(Pr) , max(Pr) ]
     ax1.fill( xx , yy1 , facecolor='green' , edgecolor = 'none' )
     


  # +++++++++++++++++++++++++++++++++++++++
  # only the reversible process existm along the botom line
  if (max(E_a)<=EE) & (T1>TA1) & (T1<T0)  :  
    ax1.fill_between( E_a , P2_a , max(P2_a) , facecolor='green' , edgecolor = 'none' )
    delta_E_n = delta_E / np.sqrt(delta_E**2.0 + delta_P2**2.0)     
    ax1.quiver(E_a[::250], P2_a[::250], delta_E_n[::250], delta_P2_n[::250] , color='red', angles='xy', edgecolors=('k'), headaxislength=5 , zorder = 2 , scale = 10.0 , headwidth = 3, headlength = 5 , linewidths=(2.0,) , width = 0.015 )
    ax1.plot( min(E_a) , min(P2_a) , marker = 'o' , markersize = 18 , color = 'lime' , zorder = 2 )
    ax1.plot( max(E_a) , max(P2_a) , marker = 'o' , markersize = 18 , color = 'black' , zorder = 2 )
  # +++++++++++++++++++++++++++++++++++++++
  # both the reversible and irresible process exist, along the top line
  if (max(E_a)>EF) & (T1>TA1) & (T1<T0)  :  
     print(E_a)
     # +++++++++++++++++++++++++++++++++++++++
     # line of reversible polarization
     P_ir = np.linspace( PB , PD, 300 )
     Pr1 = ( a * P_ir + b* P_ir**3.0 + 2.0 * b * PB**3.0) / ( a + 3.0 * b * PB**2.0 )
     print(a + 3.0 * b * PB**2.0, PB)
     # reversible electric field
     Er1 = a * Pr1 + 3.0 * b * PB**2.0 * Pr1 - 2.0 * b * PB**3.0 + gamma * Pd
     # plot line of reversible polarization
     liner = plt.plot( Er1 , Pr1 ,  color = 'gray' ,linewidth=2.0 , zorder = 1)
     #yy2 = [0.0 , 0.0 , min(Pr1) , min(Pr1) ]  
  
     # +++++++++++++++++++++++++++++++++++++++
     # area reversible polarization
     Pr = ( a * P1_a + b* P1_a**3.0 + 2.0 * b * PB**3.0) / ( a + 3.0 * b * PB**2.0 )
     print(a + 3.0 * b * PB**2.0, PB)
     # reversible electric field
     Er = a * Pr + 3.0 * b * PB**2.0 * Pr - 2.0 * b * PB**3.0 + gamma * Pd
     print('PB=',PB,'max(P1_a',max(P1_a))
    
    
     ax1.fill_between( E_a , P1_a , max(P1_a) , facecolor='gray', interpolate=True  , edgecolor = 'none')
     ax1.fill_between( E_a , P1_a , max(P1_a) , where = ( P1_a>=PB ) , facecolor='green', interpolate=True , edgecolor = 'none' )
     # # from point B to : 
     ax1.fill_between( Er1 , Pr1 , max(P1_a) ,  facecolor='green', interpolate=True , edgecolor = 'none' )
     # # rectangle
     xx = [0.0 , ED  , ED , 0.0 ]
     # area 
     yy1 = [ min(Pr1) , min(Pr1) , max(P1_a) , max(P1_a) ]
     ax1.fill( xx , yy1 , facecolor='green' , edgecolor = 'none' )
     # from point E to final point N
     ax1.fill_between( E_a , P1_a , max(PE), where = ( ( E_a<= EE ) )  , facecolor='green', interpolate=True , edgecolor = 'none' )
     
     delta_E_n = delta_E / np.sqrt(delta_E**2.0 + delta_P1**2.0)     
     ax1.quiver(E_a[::250], P1_a[::250], delta_E_n[::250], delta_P1_n[::250], color='red', angles='xy', edgecolors=('k'), headaxislength=5 , zorder = 2 , scale = 10.0 , headwidth = 3, headlength = 5 , linewidths=(2.0,) , width = 0.015)
     ax1.plot( min(E_a) , min(P1_a) , marker = 'o' , markersize = 18 , color = 'lime' , zorder = 2 )
     ax1.plot( max(E_a) , max(P1_a) , marker = 'o' , markersize = 18 , color = 'black' , zorder = 2 )

  #print 'P1=',P1 
  ## line of reversible
  ax1.plot( 0 , 0 , marker = 'x' , markersize = 10 , color = 'red' , zorder = 2 , markeredgewidth = 2 )
  ax1.plot( Pd*gamma , 0.0  , marker = 'x' , markersize = 10 , color = 'black' , zorder = 2 , markeredgewidth = 2 )
  #ax1.quiver(0 , 0 , Pd*gamma, 0.0 , color='gray', angles='xy', linewidths=(0.05,), edgecolors=('gray'), headaxislength=5 , zorder = 2 , scale = 1.5 , linestyle='dashed')
 

  ###plot figure
  ax1.text( 0.17 , -1.2 , 'T='+str(T1) , fontsize = 30)
  ax1.set_xlabel('Electric field',fontsize=30)
  ax1.set_ylabel('Polarization',fontsize=30)
  # size of the ticklabels
  plt.setp( ax1.get_xticklabels(), fontsize=30)
  plt.setp( ax1.get_yticklabels(), fontsize=30)
  # the ticklabels show
  xlim = [-0.03,0.27] 
  ylim = [-1.3 , 1.3 ]
  plt.setp( ax1 , xticks = [0.0,0.1,0.2] , yticks = [-1.2,-0.6,0.0,0.6,1.2] , xlim= xlim , ylim = ylim )
  # tick thickness
  ax1.tick_params(width=3,size=10)
  # save the figure
  #plot the P0, and E0
  ax1.plot( [0.0,0.0] , ylim , color="gray" , linestyle = "--" )
  ax1.plot( xlim , [0.0,0.0] , color="gray" , linestyle = "--" )      
  j = j + 1
  plt.savefig('hys26-T'+str(j)+'.pdf',bbox_inches='tight', format='pdf')  




