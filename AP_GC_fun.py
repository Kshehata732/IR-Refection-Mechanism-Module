''' These calculations include the reflection to the same AP and GC'''
#import fpcetl
import numpy as np
import input as inp
import fun2
#############################################
#############################################
################################################
VF_AG = fun2.VF_AG()
VF_GA = ((inp.r_ab_o**2)/(inp.r_g_i**2))*np.array(fun2.VF_AG())
rho_G = inp.rho_G
#######################################################
###############################################################################################################################################################################################
''' AP functions '''
###################################################### 
#The radiation due to the view factor from AP to GC and the diffused reflection from AP to GC
#Calculating the first reflection first
def Ref1(i, k, Alist):
    Ref1, Ref2, Ref3, Ref4 = 0, 0 ,0 ,0
    for j in range(1, len(VF_AG)):
        if (2*i + j > (inp.GC_seg_no - 1)) and (i + j > (inp.AP_seg_no - 1)):
            Ref1 += (Alist[k-1][(i+j)- inp.AP_seg_no][0]**4)*VF_AG[j]*rho_G[(2*i+j)- inp.GC_seg_no] + (Alist[k-1][i - j][0]**4)*VF_AG[j]*rho_G[2*i-j]
        elif i + j > (inp.AP_seg_no - 1):
            Ref2 += (Alist[k-1][(i+j)- inp.AP_seg_no][0]**4)*VF_AG[j]*rho_G[2*i+j] + (Alist[k-1][i - j][0]**4)*VF_AG[j]*rho_G[2*i-j]
        elif 2*i + j > (inp.GC_seg_no - 1):   
            Ref3 += (Alist[k-1][i+j][0]**4)*VF_AG[j]*rho_G[(2*i+j)- inp.GC_seg_no] + (Alist[k-1][i - j][0]**4)*VF_AG[j]*rho_G[2*i-j]
        else:
            Ref4 += (Alist[k-1][i+j][0]**4)*VF_AG[j]*rho_G[2*i+j] + (Alist[k-1][i - j][0]**4)*VF_AG[j]*rho_G[2*i-j]         
        #xa.append((x + Ref1+ Ref2+ Ref3 + Ref4))
    xa = (Ref1+ Ref2+ Ref3 + Ref4)
    return xa
###############################################
##############################################
def Ref2(i, k, Alist):
    Ref_1R, Ref_2R, Ref_3R, Ref_4R = 0, 0, 0, 0
    Ref_1L, Ref_2L, Ref_3L, Ref_4L = 0, 0, 0, 0
    Ref11, Ref2, Ref3, Ref4 = 0, 0, 0, 0
        #x = AP[i]*(rho_G[2*i]**2)*inp.rho_abir*(F_AG[0]**2)
    for j in range(1, len(VF_AG)):
        if (i + j > (inp.AP_seg_no-1)) and ((2*i + j)>(inp.GC_seg_no-1)):
            Ref_1R += (Alist[k-1][(i+j)- inp.AP_seg_no][0]**4)*VF_AG[j]*rho_G[(2*i+j)-inp.GC_seg_no]*inp.rho_abir*rho_G[2*((i+j) - inp.GC_seg_no)]*VF_AG[0] + VF_AG[j]*rho_G[(2*i+j)-inp.GC_seg_no]*inp.rho_abir*Ref1(((i+j) - inp.AP_seg_no), k, Alist)
            Ref_1L += (Alist[k-1][i-j][0]**4)*VF_AG[j]*rho_G[2*i-j]*inp.rho_abir*rho_G[2*(i-j)]*VF_AG[0] + VF_AG[j]*rho_G[2*i-j]*inp.rho_abir*Ref1((i - j), k, Alist) 
            Ref11 = Ref_1R + Ref_1L
        elif i + j > inp.AP_seg_no-1:
            Ref_2R += (Alist[k-1][(i+j)- inp.AP_seg_no][0]**4)*VF_AG[j]*rho_G[2*i+j]*inp.rho_abir*rho_G[2*((i+j) - inp.GC_seg_no)]*VF_AG[0] + VF_AG[j]*rho_G[2*i+j]*inp.rho_abir*Ref1(((i+j) - inp.AP_seg_no), k, Alist)
            Ref_2L += (Alist[k-1][i-j][0]**4)*VF_AG[j]*rho_G[2*i-j]*inp.rho_abir*rho_G[2*(i-j)]*VF_AG[0] + VF_AG[j]*rho_G[2*i-j]*inp.rho_abir*Ref1((i - j), k, Alist)  
            Ref2 = Ref_2R + Ref_2L
        elif (2*i + j)>inp.GC_seg_no-1:
            Ref_3R += (Alist[k-1][i+j][0]**4)*VF_AG[j]*rho_G[(2*i+j)-inp.GC_seg_no]*inp.rho_abir*rho_G[2*(i+j)]*VF_AG[0] + VF_AG[j]*rho_G[(2*i+j)-inp.GC_seg_no]*inp.rho_abir*Ref1((i + j), k, Alist)
            Ref_3L += (Alist[k-1][i-j][0]**4)*VF_AG[j]*rho_G[2*i-j]*inp.rho_abir*rho_G[2*(i-j)]*VF_AG[0] + VF_AG[j]*rho_G[2*i-j]*inp.rho_abir*Ref1((i - j), k, Alist)  
            Ref3 = Ref_3R + Ref_3L         

        else:
            Ref_4R += (Alist[k-1][i+j][0]**4)*VF_AG[j]*rho_G[2*i+j]*inp.rho_abir*rho_G[2*(i+j)]*VF_AG[0] + VF_AG[j]*rho_G[2*i+j]*inp.rho_abir*Ref1((i + j), k, Alist)
            Ref_4L += (Alist[k-1][i-j][0]**4)*VF_AG[j]*rho_G[2*i-j]*inp.rho_abir*rho_G[2*(i-j)]*VF_AG[0] + VF_AG[j]*rho_G[2*i-j]*inp.rho_abir*Ref1((i - j), k, Alist)      
            Ref4 = Ref_4R + Ref_4L  
        #y.append((x +Ref11 + Ref2 + Ref3 + Ref4))
    y = Ref11 + Ref2 + Ref3 + Ref4
    return y
####################################################################
epislon_Gi = inp.epislon_Gi    
###############################################################################################
def GA_emis(i, k, Glist): #i is the index number for AP. So that I don't care about offset or normal GC when it comes to this calculations.
    GA1, GA2 = 0, 0       
    #if i%2 == 0:    #normal (without off set)
    GA0 =  epislon_Gi[2*i]*VF_GA[0]*(Glist[k-1][2*i][0]**4) 
    for j in range(1,len(VF_AG)): 
        if 2*i+j > inp.GC_seg_no - 1:
            GA1 += epislon_Gi[(2*i+j)-inp.GC_seg_no]*VF_GA[j]*(Glist[k-1][(2*i+j)-inp.GC_seg_no][0]**4) + epislon_Gi[2*i-j]*VF_GA[j]*(Glist[k-1][2*i-j][0]**4)     #GA_n for normal (without offset position)
        else:
            GA2 += epislon_Gi[2*i+j]*VF_GA[j]*(Glist[k-1][2*i+j][0]**4) + epislon_Gi[2*i-j]*VF_GA[j]*(Glist[k-1][2*i-j][0]**4)
    x = GA0 + GA1 + GA2
    return x
########################################################################################################################################################################################################
                      ################################################################################################################################33
''' GC functions '''              
##############################################
################################  
VF_GG = fun2.VF_GG()
###############
def GG_emis(i, k, Glist):
    GG1, GG2 = 0, 0       #
    for j in range(len(VF_GG)): #Emissison from GC to GC CV's. According to VF_GG cal. we don't need to start for loop with 1
        if i+j > inp.GC_seg_no - 1:
            GG1 += epislon_Gi[(i+j) - inp.GC_seg_no]*VF_GG[j]*(Glist[k-1][(i+j) - inp.GC_seg_no][0]**4) + epislon_Gi[i-j]*VF_GG[j]*(Glist[k-1][i - j][0]**4)
        else:
            GG2 += epislon_Gi[i+j]*VF_GG[j]*(Glist[k-1][i+j][0]**4) + epislon_Gi[i-j]*VF_GG[j]*(Glist[k-1][i - j][0]**4)
    x = GG1 + GG2        
    return x        
#######################################################################     
def GC_1ref(i, k, Alist):
    GC_1ref = []
    AG_1ref_1n, AG_1ref_2n, AG_1ref_1off, AG_1ref_2off  = 0, 0, 0, 0
    if i%2 == 0:
        AG_1ref0 = (1-rho_G[i])*(Alist[k-1][int(i/2)][0]**4)*VF_AG[0]
        for j in range(1, int((len(VF_AG)-1)/2)):
            if (i/2 + j) > inp.AP_seg_no - 1:
                AG_1ref_1n += (1-rho_G[i])*((Alist[k-1][int(i/2 + j) - inp.AP_seg_no][0]**4)*VF_AG[2*j] + (Alist[k-1][int(i/2 - j)][0]**4)*VF_AG[2*j])                    
            else:
                AG_1ref_2n += (1-rho_G[i])*((Alist[k-1][int(i/2 + j)][0]**4)*VF_AG[2*j] + (Alist[k-1][int(i/2 - j)][0]**4)*VF_AG[2*j])
        GC_1ref.append(AG_1ref_1n + AG_1ref_2n + AG_1ref0)
            #print(AG_1ref_n, AG_1ref0)
    else:
        for j in range(1, int(len(VF_AG)/2)):
            if ((i+1)/2 + j) > inp.AP_seg_no - 1:
                AG_1ref_1off += (1-rho_G[i])*((Alist[k-1][int((i+1)/2 + j) - inp.AP_seg_no][0]**4)*VF_AG[int(2*j -1)] + (Alist[k-1][int((i+1)/2 - j)][0]**4)*VF_AG[int(2*j -1)])                    
            else:
                AG_1ref_2off += (1-rho_G[i])*((Alist[k-1][int((i+1)/2 + j)][0]**4)*VF_AG[int(2*j -1)] + (Alist[k-1][int((i+1)/2 - j)][0]**4)*VF_AG[int(2*j -1)])
        GC_1ref.append(AG_1ref_1off + AG_1ref_2off)
    GC_ref = sum(GC_1ref)    
    return GC_ref            
###############################################################################
###############################################################################
def GC_2ref(i, k, Alist):
    GC_2ref = []
    AG_2ref_n_1R, AG_2ref_n_1L, AG_2ref_off_1R, AG_2ref_off_1L = 0, 0, 0, 0
    AG_2ref_n_2R, AG_2ref_n_2L, AG_2ref_off_2R, AG_2ref_off_2L = 0, 0, 0, 0
    x, y = 0, 0
    if i%2 == 0:
        AG_2ref0 = (1-rho_G[i])*(Alist[k-1][int(i/2)][0]**4)*(VF_AG[0]**2)*inp.rho_abir*rho_G[i]
        for j in range(1, int((len(VF_AG)-1)/2)):
            if (i/2 + j) > inp.AP_seg_no - 1:
                AG_2ref_n_1R += (1-rho_G[i])*((Alist[k-1][int(i/2 + j) - inp.AP_seg_no][0]**4)*VF_AG[2*j]*inp.rho_abir*rho_G[2*(int(i/2 + j) - inp.AP_seg_no)]*VF_AG[0] + VF_AG[2*j]*inp.rho_abir*Ref1((int(i/2 + j) - inp.AP_seg_no), k, Alist))
                AG_2ref_n_1L += (1-rho_G[i])*((Alist[k-1][int(i/2 - j)][0]**4)*VF_AG[2*j]*inp.rho_abir*rho_G[2*int(i/2 - j)]*VF_AG[0] + VF_AG[2*j]*inp.rho_abir*Ref1((int(i/2 - j)), k, Alist))                   
                    #AG_2ref_n = AG_2ref_n_R + AG_2ref_n_L
                x  = AG_2ref_n_1R + AG_2ref_n_1L
                    #print(x)
            else:
                AG_2ref_n_2R += (1-rho_G[i])*((Alist[k-1][int(i/2 + j)][0]**4)*VF_AG[2*j]*inp.rho_abir*rho_G[2*int(i/2 + j)]*VF_AG[0] + VF_AG[2*j]*inp.rho_abir*Ref1((int(i/2 + j)), k, Alist)) 
                AG_2ref_n_2L += (1-rho_G[i])*((Alist[k-1][int(i/2 - j)][0]**4)*VF_AG[2*j]*inp.rho_abir*rho_G[2*int(i/2 - j)]*VF_AG[0] + VF_AG[2*j]*inp.rho_abir*Ref1((int(i/2 - j)), k, Alist))
                    #AG_2ref_n = AG_2ref_n_R + AG_2ref_n_L
                y =  AG_2ref_n_2R + AG_2ref_n_2L
           # print(x+y)
                 #   break
        GC_2ref.append((x+y+AG_2ref0))            
            #print(AG_1ref_n, AG_1ref0)
    else:
        for j in range(1, int(len(VF_AG)/2)):
            if ((i+1)/2 + j) > inp.AP_seg_no - 1:
                AG_2ref_off_1R += (1-rho_G[i])*((Alist[k-1][int((i+1)/2 + j) - inp.AP_seg_no][0]**4)*VF_AG[int(2*j -1)]*inp.rho_abir*rho_G[2*(int((i+1)/2 + j) - inp.AP_seg_no)]*VF_AG[0] + VF_AG[int(2*j -1)]*inp.rho_abir*Ref1((int((i+1)/2 + j) - inp.AP_seg_no), k, Alist)) 
                AG_2ref_off_1L += (1-rho_G[i])*((Alist[k-1][int((i+1)/2 - j)][0]**4)*VF_AG[int(2*j -1)]*inp.rho_abir*rho_G[2*int((i+1)/2 - j)]*VF_AG[0] + VF_AG[int(2*j -1)]*inp.rho_abir*Ref1((int((i+1)/2 - j)), k, Alist))                    
                x = AG_2ref_off_1R + AG_2ref_off_1L
            else:
                AG_2ref_off_2R += (1-rho_G[i])*((Alist[k-1][int((i+1)/2 + j)][0]**4)*VF_AG[int(2*j -1)]*inp.rho_abir*rho_G[2*int((i+1)/2 + j)]*VF_AG[0] + VF_AG[int(2*j -1)]*inp.rho_abir*Ref1((int((i+1)/2 + j)), k, Alist)) 
                AG_2ref_off_2L += (1-rho_G[i])*((Alist[k-1][int((i+1)/2 - j)][0]**4)*VF_AG[int(2*j -1)]*inp.rho_abir*rho_G[2*int((i+1)/2 - j)]*VF_AG[0] + VF_AG[int(2*j -1)]*inp.rho_abir*Ref1((int((i+1)/2 - j)), k, Alist))                  
                y = AG_2ref_off_2R + AG_2ref_off_2L
        GC_2ref.append(x+y)
    GC_ref = sum(GC_2ref)    
    return GC_ref    
###############################################################################################            
##############################################################################################
    #################################################################################
def AG_emis(i, k, Alist): 
    x = []
    AG_1n, AG_2n, AG_1off, AG_2off  = 0, 0, 0, 0
    if i%2 == 0:
        AG0 = (Alist[k-1][int(i/2)][0]**4)*VF_AG[0]
        for j in range(1, int((len(VF_AG)-1)/2)):
            if (i/2 + j) > inp.AP_seg_no - 1:
                AG_1n += (Alist[k-1][int(i/2 + j) - inp.AP_seg_no][0]**4)*VF_AG[2*j] + (Alist[k-1][int(i/2 - j)][0]**4)*VF_AG[2*j]                   
            else:
                AG_2n += (Alist[k-1][int(i/2 + j)][0]**4)*VF_AG[2*j] + (Alist[k-1][int(i/2 - j)][0]**4)*VF_AG[2*j]
        x.append(AG_1n + AG_2n + AG0)
            #print(AG_1ref_n, AG_1ref0)
    else:
        for j in range(1, int(len(VF_AG)/2)):
            if ((i+1)/2 + j) > inp.AP_seg_no - 1:
                AG_1off += (Alist[k-1][int((i+1)/2 + j) - inp.AP_seg_no][0]**4)*VF_AG[int(2*j -1)] + (Alist[k-1][int((i+1)/2 - j)][0]**4)*VF_AG[int(2*j -1)]                   
            else:
                AG_2off += (Alist[k-1][int((i+1)/2 + j)][0]**4)*VF_AG[int(2*j -1)] + (Alist[k-1][int((i+1)/2 - j)][0]**4)*VF_AG[int(2*j -1)]
        x.append(AG_1off + AG_2off)
    y = sum(x)    
    return y     
        #####################################################################