import numpy as np
import input as inp
##################################
'''VF from AP to GC'''
def VF_AG():
    i = 0
    alpha = inp.AP_seg_rad/4
    VF_AG = []
    while (np.pi/2 +alpha + i*(inp.AP_seg_rad/2)) < np.pi:       
        y1 = np.sqrt((inp.r_g_i*np.cos(np.pi/2 + alpha + i*(inp.AP_seg_rad/2)) - inp.r_ab_o*np.cos(np.pi/2 + (inp.AP_seg_rad/2)))**2 + (inp.r_g_i*np.sin(np.pi/2 +alpha + i*(inp.AP_seg_rad/2)) - inp.r_ab_o*np.sin(np.pi/2 + (inp.AP_seg_rad/2)))**2)
        y2 = np.sqrt((inp.r_g_i*np.cos(np.pi/2 - alpha + i*(inp.AP_seg_rad/2)) - inp.r_ab_o*np.cos(np.pi/2 - (inp.AP_seg_rad/2)))**2 + (inp.r_g_i*np.sin(np.pi/2 - alpha + i*(inp.AP_seg_rad/2)) - inp.r_ab_o*np.sin(np.pi/2 - (inp.AP_seg_rad/2)))**2)
        x1 = np.sqrt((inp.r_g_i*np.cos(np.pi/2 + alpha + i*(inp.AP_seg_rad/2)) - inp.r_ab_o*np.cos(np.pi/2 - (inp.AP_seg_rad/2)))**2 + (inp.r_g_i*np.sin(np.pi/2 + alpha + i*(inp.AP_seg_rad/2)) - inp.r_ab_o*np.sin(np.pi/2 - (inp.AP_seg_rad/2)))**2)
        x2 = np.sqrt((inp.r_g_i*np.cos(np.pi/2 - alpha + i*(inp.AP_seg_rad/2)) - inp.r_ab_o*np.cos(np.pi/2 + (inp.AP_seg_rad/2)))**2 + (inp.r_g_i*np.sin(np.pi/2 - alpha + i*(inp.AP_seg_rad/2)) - inp.r_ab_o*np.sin(np.pi/2 + (inp.AP_seg_rad/2)))**2)
        F = ((x1+x2)-(y1+y2))/(2*inp.r_ab_o*inp.AP_seg_rad)
        if F > 0:
            VF_AG.append(F)
            i+=1
        else:
            break
    return VF_AG
#############################################3        
'''VF from GC to GC'''
def VF_GG():
    i = 0
    alpha = inp.AP_seg_rad/4
    VF_GG = []
    while (- alpha + i*(inp.AP_seg_rad/2)) < np.pi/2:
        y1 = np.sqrt((inp.r_g_i*np.cos(-3*alpha + i*(inp.AP_seg_rad/2)) - inp.r_g_i*np.cos(np.pi/2 + alpha))**2 + (inp.r_g_i*np.sin( - 3*alpha + i*(inp.AP_seg_rad/2))- inp.r_g_i*np.sin(np.pi/2 + alpha))**2)
        y2 = np.sqrt((inp.r_g_i*np.cos(- alpha + i*(inp.AP_seg_rad/2)) - inp.r_g_i*np.cos(np.pi/2 - alpha))**2 + (inp.r_g_i*np.sin( - alpha + i*(inp.AP_seg_rad/2))- inp.r_g_i*np.sin(np.pi/2 - alpha))**2)
        x1 = np.sqrt((inp.r_g_i*np.cos(-3*alpha + i*(inp.AP_seg_rad/2)) - inp.r_g_i*np.cos(np.pi/2 - alpha))**2 + (inp.r_g_i*np.sin( - 3*alpha + i*(inp.AP_seg_rad/2))- inp.r_g_i*np.sin(np.pi/2 - alpha))**2)
        x2 = np.sqrt((inp.r_g_i*np.cos(- alpha + i*(inp.AP_seg_rad/2)) - inp.r_g_i*np.cos(np.pi/2 + alpha))**2 + (inp.r_g_i*np.sin( - alpha + i*(inp.AP_seg_rad/2))- inp.r_g_i*np.sin(np.pi/2 + alpha))**2)
        F = ((x1+x2)-(y1+y2))/(2*inp.r_g_i*(inp.GC_seg_rad/2))
        VF_GG.append(F)
#        print (- alpha + i*inp.seg_rad)
        i += 1
    x = VF_GG[::-1]     
    return x
        ###################################################################################################
def VF_sky_GC():
   '''Calculating the view factor from sky to GC'''
   #y = inp.Rc*np.sin(inp.gamma)
   #x = inp.Rc*np.sin(inp.gamma) + (inp.r_g_o*inp.seg_rad)
   #F = (2*x - 2*y)/(inp.Rc*inp.Beta)
   #F = 0.000473
   F = 1/inp.GC_seg_no
   return F
#####################################################################
 ########################################################################33
def AP_sol():  
    
    
    S =[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,203.3501206,271.1334942,271.1334942,413.1500887,555.1666832,555.1666832,1583.288937,2017.724179,2109.452003,2109.452003,2109.452003,4831.891842,7876.982724,8844.935856,8844.935856,12403.59543,15962.255,19577.92306,23193.59111,23922.62039,24651.64966,20074.69234,15497.73502,10872.64923,6247.563434,4471.038459,2694.513483,2694.513483,2260.648749,959.0545482,1125.989561,1626.794599,1349.971821,519.5034875,519.5034875,519.5034875,129.8758719,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    return S 
    
#####################################################################################
def GC_sol():

    S = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,72.90648015,145.8129603,72.90648015,145.8129603,145.8129603,145.8129603,466.1113539,786.4097475,466.1113539,786.4097475,786.4097475,786.4097475,1647.998369,2509.586991,1647.998369,2509.586991,2509.586991,2509.586991,3864.835927,5220.084862,3864.835927,5220.084862,5220.084862,5220.084862,5220.084862,5220.084862,5220.084862,5220.084862,5070.4779,4920.870938,9314.977121,13409.86938,9165.370159,13409.86938,13409.86938,13409.86938,15336.18859,17262.5078,15336.18859,17262.5078,17262.5078,17262.5078,20338.31771,23414.12763,23037.6237,28812.7396,26113.43362,28812.7396,28812.7396,28812.7396,28812.7396,28812.7396,28812.7396,28812.7396,27635.28995,26457.84029,27635.28995,26457.84029,26457.84029,26457.84029,26457.84029,26457.84029,24970.20574,23482.57118,24970.20574,23482.57118,22254.20913,21025.84707,22254.20913,21025.84707,17260.29458,13494.74208,17260.29458,13494.74208,13494.74208,13494.74208,13494.74208,13494.74208,9149.218917,4803.695754,9149.218917,4803.695754,4803.695754,4803.695754,4803.695754,4803.695754,4803.695754,4803.695754,4803.695754,4803.695754,3410.812738,2017.929722,3410.812738,2017.929722,2017.929722,2017.929722,1156.400548,294.8713741,1156.400548,294.8713741,294.8713741,294.8713741,147.4356871,0,147.4356871,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    return S
        