#packages
import numpy as np
'''      discretization   '''
#The CVs of the reciever start from zero at the mid top that face the sun and goes right down
AP_seg_no = 100
GC_seg_no = AP_seg_no*2
AP_seg_ang = 360.0/AP_seg_no 
GC_seg_ang = 360.0/GC_seg_no          #segment size in angle from the center of the pipe
AP_seg_rad = AP_seg_ang/180*np.pi  
GC_seg_rad = GC_seg_ang/180*np.pi     #segement size in radian
dl = 1                                # unit length of the AP and GC

'''           Heat Transfer Fluid   '''
T_fin = 273                           # HTF temp.
T_amb = 300                           #Ambient temp.
T_sky = T_amb - 8                     #Sky temp.
m = 0.01                              # mass flow rate

''' For molten salt '''
density = 882 
k_f = 0.124
Cp = 1711
h_f = 500
visco = 0.00386
segma = 5.6697e-8
###############################################################################
'''      Absorber pipe properties '''
#size
r_ab_i = 0.059/2
r_ab_o = 0.06/2 
dt_ab = r_ab_o - r_ab_i
r_ab_m=(r_ab_o + r_ab_i)/2.0
dA_ab_o=r_ab_o*AP_seg_rad*dl
dA_ab_i=r_ab_i*AP_seg_rad*dl
##############
#Physical properties
k_ab = 16.2
rho_abir = 0.14                      #For bare or HM absorber pipe
epislon_ab = 1 - rho_abir
alph_abvis = 0.86
###############################################################################
'''                 Glass Cover   '''
#Physical properties 
r_g_i= 0.078/2
r_g_o= 0.08/2
dt_g = r_g_o - r_g_i
r_g_m= (r_g_o + r_g_i)/2.0
dA_g_i= r_g_i*GC_seg_rad*dl
dA_g_o= r_g_o*GC_seg_rad*dl
####################
#### Upper half ####
''' IR radiaiton     '''
#glass transmissivity in IR is zero
k_g_up = [220]            #For Alminum outercover
rho_gir_up = [0.90]       #Cavity mirror inside
epislon_gi_up = [0.1]    #Cavity mirror inside
epislon_go_up = [0.1]    #Cavity mirror outside (IR)

# k_g_up = [1.04]          #For HM and glass test
# rho_gir_up = [0.85]      # Hotmirror test
# epislon_gi_up = [0.15]   # Hotmirror test
# epislon_go_up = [0.15]   # Hotmirror test

''' visible radiation '''
tau_g_up = 0               #For outer cavity cover
rho_gvis_up = 0.9         #Alumium outer surface reflectivity

# tau_g_up = 0.875          # HM test
# rho_gvis_up = 0.1         # HM test (It was not activated)

alph_gvis_up = [round(1 - rho_gvis_up - tau_g_up, 4)] 
#alph_gvis_up = [0.5]      # Cavity 
######################################################
'''  It is important part for the opening at which the solar flux inter the system  '''
###lower half ####
''' IR radiaiton '''
k_g_dow = [1.04]           #Cavity opening or HM & glass casees (All have almost the same value)

# rho_gir_dow = [0.85]      # Cavity opening with with HM or HM
# epislon_gi_dow = [0.15]   # Cavity opening with with HM Or HM
# epislon_go_dow = [0.15]   # Cavity opening with with HM or HM

rho_gir_dow = [0.14]       # Cavity opening with with bare glass
epislon_gi_dow = [0.86]    # Cavity opening with with bare glass
epislon_go_dow = [0.86]    # Cavity opening with with bare glass
############
''' Visible radiaiton '''
# tau_g_dow = 0.875         #HM
# rho_gvis_dow = 0.1        #HM

tau_g_dow = 0.935          #bare glass
rho_gvis_dow = 0.045       #bare glass

alph_gvis_dow = [round(1 - rho_gvis_dow - tau_g_dow, 4)]
######################################################
######################################################
 # Outer cover and AP are divided into 200 and 100 cvs, respectively incase of activating Reflection mechansim  

''' Original Case with rim angle 68 degrees, half opening size = 100 CVs means 50 up, 100 down, 50 up '''  
k_g = k_g_up*84 + k_g_dow*33 + k_g_up*83     # cavity opening with 100 CVs (outer cover). AP CVs have the same properties.
rho_G = rho_gir_up*84 + rho_gir_dow*33 + rho_gir_up*83  
tau_g = [tau_g_up]*84 + [tau_g_dow]*33 + [tau_g_up]*83
alph_vis_g = alph_gvis_up*84 + alph_gvis_dow*33 + alph_gvis_up*83 
epislon_Go = epislon_go_up*84 + epislon_go_dow*33 + epislon_go_up*83
epislon_Gi = epislon_gi_up*84 + epislon_gi_dow*33 + epislon_gi_up*83
###############
''' Air outside the glass cover '''
v = 2.6 #m/s the speed of wind
h_g = h_w = (4*v**0.58)*(2*r_g_o)**(-0.42)  #h_w = 4*v**0.58 * d**-0.42 , wher d is the diameter of the glass cover
#######Convection parameters ##################################
#Re = m/(np.pi*visco*r_ab_i**2)
#Pr =  visco*Cp/k_f
#Nu = 0.3+((0.62*Re**0.5*Pr**0.3333)/(1+(0.4/Pr)**0.6667)**0.25)*(1+(Re/282000)**0.625)**0.8
#Rc = 375.3567 
#Beta = 7
#gamma = 0.40278
