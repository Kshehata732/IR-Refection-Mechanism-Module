'''Starting the programe''' 
import numpy as np
import time
import input as inp
import pandas as pd
import AP_GC_fun as fun1
import fun2 
#import xlsxwriter
######################################################################################################################
##################################################################################################################
num_unit_length= 5                            # Total length of the Reciever unit 
Iteration = 300                                # Iteration number for accuracy
Excel_file_name = 'test.xlsx'
#Note: Maximum error is 5E-6 for this current version
#e = 0.0000001                                 # Convergence criteria
##########################################################
''' define the output lists for final results '''
T_f = [inp.T_fin]
Q_conv_HTF = []     #HTF heat gain
#Ts = [[] for _ in range(num_unit_length)]     #inner surface temperature of the absorber pipe
''' Intialization '''
#intialize the first row of Ta and Tg 
Ta = [[] for _ in range(num_unit_length+1)]   # The temperature of the outer surface of the absorber pipe
Tg = [[] for _ in range(num_unit_length+1)]   # The temperature of the outer surface of the glass cover
Ta[0].append(np.full(inp.AP_seg_no, inp.T_fin, dtype=np.float64))
Tg[0].append(np.full(inp.GC_seg_no, inp.T_amb, dtype=np.float64)) 
#####################################################
'''The head function for the simulation'''
def starter():     
    time_start = time.time()                          # Time start 
    print ("The calculations have been started...Proceeding...")
    print()
    for It in range(num_unit_length):                 #Iteration along the RU lenght 
        Alist, Glist = Iter_ab(Iteration, It)  # calling the function that iterates the solution of Ta and Tg
        np.seterr(over='raise')
        j = Iteration - 1
        #j = count + 1                                 #Parameter that define the last element of the iteration
## Way to insert the Ta and Tg profile temperatures obtained from functions Glist and Alist at each length of RU in         
        Tgs = []
        Tas = []
        #Ts1 = []
        for i in range(inp.AP_seg_no):                         
            Tas.append(Alist[j][i][0])
        Ta[It+1].append(Tas)
        for k in range(inp.GC_seg_no):
            Tgs.append(Glist[j][k][0])
        Tg[It+1].append(Tgs)               
#Core calculations to find the final values that can be obtained from Ta and Tg        
        Ts1 = (inp.k_ab*(np.array(Ta[It], dtype=np.float64)) + inp.h_f*inp.dt_ab*T_f[It])/(inp.k_ab + inp.h_f*inp.dt_ab)       
        Qconv_CV = inp.h_f*inp.dA_ab_i*(Ts1 - T_f[It])    #It calculates Q HTF gain for each element inside the list [Seg. number]
        Qconv = np.sum(Qconv_CV)                            #Sum the all Q HTF for each segment arounf the cicumference 
        Q_conv_HTF.append(Qconv)                            #Append the Q HTF for the whole dL lenght in one list 
        T_fout = (Qconv/(inp.m*inp.Cp)) + T_f[It]         # Calculating the outlet temperature at each dL
        T_f.append(T_fout)
        #Ts.append(Ts1.tolist())
        print ("Unit length number : "+str(It)+ " Outlet Temperature : "+str(T_f[It])+" HTF Heat Gain : "+str(Q_conv_HTF[It]))
        print()
#Writing down the results in excel file      
    #del T_f[0]                                 # delete the first elemnt T_f list without      
    writer = pd.ExcelWriter(Excel_file_name)
    df1 = pd.DataFrame({'HTF Temp.':T_f})
    df2 = pd.DataFrame({'Q_htf_cv':Q_conv_HTF})

    ####Preparing the AP and GC temperature profiles
    Ta_st = ['Ta'+str(x) for x in range(num_unit_length)]
    Ta_val = [Ta[i][0] for i in range(num_unit_length)]
    Ta_dic = {key:value for key, value in zip(Ta_st, Ta_val)}
    df3 = pd.DataFrame(Ta_dic)
    
    Tg_st = ['Tg'+str(x) for x in range(num_unit_length)]
    Tg_val = [Tg[i][0] for i in range(num_unit_length)]
    Tg_dic = {key:value for key, value in zip(Tg_st, Tg_val)}
    df4 = pd.DataFrame(Tg_dic) 

    df1.to_excel(writer, 'HTF', index = True) 
    df2.to_excel(writer, 'QHTF', index = True)
    df3.to_excel(writer, 'Ta', index = True)
    df4.to_excel(writer, 'Tg', index = True)
    
    writer.save()    
    time_end = time.time()
    duration = int(time_end-time_start)/60 # in minutes
    print("Fluid temperatures "+str(T_f))
    print()
    print("HTF heat gain "+str(Q_conv_HTF))
    print()    
    print("Your calcultion took "+str(duration)+" minutes")
##############################################################  
    ''' Cavity opening sizes are controlled by the properties of the CVs associated with'''
#####################################################  
# Outer cover and AP are divided into 200 and 100 cvs, respectively incase of activating Reflection mechansim     
#############################################################   
'''The main calculation to fill the Iteration matrix  '''
def Iter_ab(Iteration, It):  
    VF_AG = fun2.VF_AG()
    coeff_ab = (inp.k_ab*inp.dt_ab*inp.dl)/(inp.r_ab_m*inp.AP_seg_rad)
    AP_sol = fun2.AP_sol()
    GC_sol = fun2.GC_sol()

#Building Iteration matrix    
    A0 = [[[] for _ in range(inp.AP_seg_no)] for _ in range(1)]             # Each list represents one iteration and include AP segments
    for i in range(inp.AP_seg_no):
        A0[0][i].append(Ta[It][0][i])
    A = [[[] for _ in range(inp.AP_seg_no)] for _ in range(Iteration - 1)]  # Increase the size of the matrix in order to fit the specified number of Iteration. 
    Alist = A0 + A
    
    G0 = [[[] for _ in range(inp.GC_seg_no)] for _ in range(1)]             # Each list represents one iteration and include AP segments
    for i in range(inp.GC_seg_no):
        G0[0][i].append(Tg[It][0][i])
    G = [[[] for _ in range(inp.GC_seg_no)] for _ in range(Iteration - 1)]  # Increase the size of the matrix in order to fit the specified number of Iteration. 
    Glist = G0 + G
    for j in range(1, Iteration):
        for i in range (inp.AP_seg_no):
            Ref1 = fun1.Ref1(i, j, Alist)
            Ref2 = fun1.Ref2(i, j, Alist)
            GA_emis = fun1.GA_emis(i, j, Glist)
           # AP_emis_all = fun1.AP_emis_all(i, j, Alist)
            a_ab = (2*inp.k_ab*inp.dt_ab*inp.dl)/(inp.r_ab_m*inp.AP_seg_rad) + inp.h_f*inp.dA_ab_i \
            + (4*inp.segma*inp.epislon_ab*inp.dA_ab_o*(1 + inp.rho_abir*inp.rho_G[2*i]*VF_AG[0]*(sum(VF_AG) - VF_AG[0])))*(Alist[j-1][i][0]**3) - (1-inp.rho_abir) \
            *(4*inp.segma*inp.epislon_ab*inp.dA_ab_o*(inp.rho_G[2*i]*VF_AG[0] + (inp.rho_G[2*i]**2)*inp.rho_abir*(VF_AG[0]**2))*(Alist[j-1][i][0]**3)) #may be we need to include xa[i]                         

            b_ab = inp.alph_abvis*inp.tau_g[2*i] *AP_sol[i]*inp.dA_ab_o + inp.h_f*inp.dA_ab_i*T_f[It] \
            + (3*inp.epislon_ab*inp.segma*inp.dA_ab_o*(1 + inp.rho_abir*inp.rho_G[2*i]*VF_AG[0]*(sum(VF_AG) - VF_AG[0])))*(Alist[j-1][i][0]**4) - (1-inp.rho_abir)*(3*inp.segma*inp.epislon_ab*inp.dA_ab_o*(inp.rho_G[2*i]*VF_AG[0] + (inp.rho_G[2*i]**2)*inp.rho_abir*(VF_AG[0]**2)))*(Alist[j-1][i][0]**4) \
            + (1-inp.rho_abir)*(inp.epislon_ab*inp.segma*inp.dA_ab_o)*(Ref1 + Ref2) \
            + (1-inp.rho_abir)*inp.segma*inp.dA_g_i*(GA_emis) 
                
            if i == inp.AP_seg_no - 1:
                T_a = (coeff_ab*Alist[j-1][0][0]+coeff_ab*Alist[j-1][i-1][0] + b_ab)/ a_ab
            elif i == 0:
                T_a = (coeff_ab*Alist[j-1][i+1][0]+coeff_ab*Alist[j-1][inp.AP_seg_no - 1][0] + b_ab)/ a_ab
            else:
                T_a = (coeff_ab*Alist[j-1][i+1][0]+coeff_ab*Alist[j-1][i-1][0] + b_ab)/ a_ab

            Alist[j][i].append(T_a)
        
        for i in range (inp.GC_seg_no):
            GC_1Ref = fun1.GC_1ref(i, j, Alist)
            GC_2Ref = fun1.GC_2ref(i, j, Alist)
            GG_emis = fun1.GG_emis(i, j, Glist)
            coeff_GC = (inp.k_g[i]*inp.dt_g*inp.dl)/(inp.r_g_m*inp.GC_seg_rad)
            
            a_g = (2*inp.k_g[i]*inp.dt_g*inp.dl)/(inp.r_g_m*inp.GC_seg_rad) + inp.h_g*inp.dA_g_o + 4*inp.segma*(inp.epislon_Gi[i]*inp.dA_g_i + inp.epislon_Go[i]*inp.dA_g_o)*(Glist[j-1][i][0]**3) 
            # Outer cover and AP are divided into 200 and 100 cvs, respectively incase of activating Reflection mechansim  
            
# Original case  with rim angle 68 degrees, half opening size = 100 CVs. For upper side of circumference i < 50 and i > 149 '''              
            if i < 50 or i > 149:
                b_g = inp.alph_vis_g[i]*GC_sol[i]*inp.dA_g_o + 3*inp.segma*(inp.epislon_Gi[i]*inp.dA_g_i + inp.epislon_Go[i]*inp.dA_g_o)*(Glist[j-1][i][0]**4) \
                + inp.h_g*inp.dA_g_o*inp.T_sky + inp.epislon_Go[i]*inp.segma*inp.dA_g_o*(inp.T_sky**4) + (inp.segma*inp.epislon_ab*inp.dA_ab_o)*(GC_1Ref + GC_2Ref) + inp.epislon_Gi[i]*inp.segma*inp.dA_g_i*GG_emis
            else:
                b_g = inp.alph_vis_g[i]*GC_sol[i]*inp.dA_g_o + 3*inp.segma*(inp.epislon_Gi[i]*inp.dA_g_i + inp.epislon_Go[i]*inp.dA_g_o)*(Glist[j-1][i][0]**4) \
                + inp.h_g*inp.dA_g_o*inp.T_amb + inp.epislon_Go[i]*inp.segma*inp.dA_g_o*(inp.T_amb**4) + (inp.segma*inp.epislon_ab*inp.dA_ab_o)*(GC_1Ref + GC_2Ref) + inp.epislon_Gi[i]*inp.segma*inp.dA_g_i*GG_emis
            
            
            if i == inp.GC_seg_no - 1:
                T_g = (coeff_GC*Glist[j-1][0][0] + coeff_GC*Glist[j-1][i-1][0] + b_g)/ a_g
            elif i == 0:
                T_g = (coeff_GC*Glist[j-1][i+1][0] + coeff_GC*Glist[j-1][(inp.GC_seg_no - 1)][0] + b_g)/ a_g
            else:
                T_g = (coeff_GC*Glist[j-1][i+1][0] + coeff_GC*Glist[j-1][i-1][0] + b_g)/ a_g        

            Glist[j][i].append(T_g)


    return Alist, Glist     
##################################################################################################################################
######################################################################################################################################       