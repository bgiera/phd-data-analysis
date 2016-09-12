# -*- coding: utf-8 -*-

###I want to put all the important libraries in their respective functions
##This function is started from data.py
#	The end product will perform everything indicated in Keck_GM_121511.pptx
#	It would be best to have inputs from the outside world, too...
#import copy
#import gzip
#import sys
import numpy as np #This is only to get the np.array() class
import EDLCalcs
from matplotlib.ticker import MaxNLocator
import colorsys
IfOnCluster = 'NotByDefault'


def ROYGBIV_map(value,max_value,sat=1):
    import colorsys
    return colorsys.hsv_to_rgb(value / max_value , sat, sat)
#    return colorsys.hsv_to_rgb(value / max_value / (1.1), sat, sat)


def get_n_HexCol(n):
    HSV_tuples = [(x*1.0/n, 0.75, 0.35) for x in range(n)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    return RGB_tuples

def floatRgb(mag, cmin, cmax):
       """
       Return a tuple of floats between 0 and 1 for the red, green and
       blue amplitudes.
       """
       import math
       try:
              # normalize to [0,1]
              x = float(mag-cmin)/float(cmax-cmin)
       except:
              # cmax = cmin
              x = 0.5
       blue = min((max((4*(0.75-x), 0.)), 1.))
       red  = min((max((4*(x-0.25), 0.)), 1.))
       green= min((max((4*math.fabs(x-0.5)-1., 0.)), 1.))
       return (red, green, blue)

def Bikerman_Threshold(Phi_bulk_max):
    """This function numerically solves the zeta-potential and surface charge 
    density demanded by the Bikerman model to obtain excess chemical potentials
    that are greater than kBT
Input: 
      None

Output:
	A vector of Sigma_fail(Phis), where Phis is a list of volume fractions
"""
    import matplotlib.pyplot as plt
    Phi_bulkS = np.linspace(1.1*Phi_bulk_max,0.0001,1)
    Phi_star=(np.e-1)/np.e #This is the volume fraction that makes the CP_ex = 1
    zetaFailS = np.log( 0.5*((np.e*Phi_star/Phi_bulkS) + np.sqrt(((np.e*Phi_star/Phi_bulkS)**2)-4))) #These zetas cause CP_ex = 1, too
#    zetaS = np.log( 0.5*((np.e*Phi_star/Phi_bulkS) - np.sqrt((np.e*Phi_star/Phi_bulkS)**2-4)))
    SigmaFailS=np.sqrt(2*(1/Phi_bulkS)*np.log(1 + 2*Phi_bulkS*np.sinh(zetaFailS/2)**2))
#    rho_p = np.exp(-zetaFailS -1)
#    rho_m = rho_p * np.exp(2*zetaFailS)
#    print -0.5*np.log(rho_m*rho_p)
#    print np.log(1 - (Phi_bulkS)*(rho_p+rho_m))
    print 'Analytical Sig = ',SigmaFailS
#    
     #The following loop structure calculates surface charge from the surface potential
    SigmaFailS=[]
    for (PhiB,zeta_f) in zip(Phi_bulkS,zetaFailS):   		
    	V=np.linspace(0.,zeta_f,500000)
    	dV=V[1]-V[0]
    	integrand=0.
    	rhom0_guess=np.min(np.exp(V[0]),Phi_star)
	for v in V:
		max_iter=-1
		PHI = rhom0_guess*(1.+np.exp(-2*v))*PhiB
		per_diff=np.abs(np.log(rhom0_guess)-v - np.log(1-PHI))
		while per_diff>0.000001:
			max_iter+=1
			rhom0_guess+=0.00001 #Not sure about this step...	
			PHI = rhom0_guess*(1.+np.exp(-2*v))*PhiB
			per_diff=np.abs(np.log(rhom0_guess)-v - np.log(1-PHI))/np.log(1-PHI)
			if max_iter==10E6:
				print 'Failure to converge!',fail
		integrand=integrand+rhom0_guess*(1-np.exp(-2*v))*dV
	SigmaFailS.append(np.sqrt(integrand))
    print integrand

    print 'Solved Sig = ',SigmaFailS
    print len(SigmaFailS),len(Phi_bulkS),len(zetaFailS)
    
#    Phi_max=0.65
#    SigmaFail=[]
#    zetaS=[]
#    zeta_guess = 0.0
#    for Phi_bulk in Phi_bulkS: #Evaluating from large to small Phi_bulk is quicker
#    	max_iter = -1
#    	per_dif = -0.5 *(-2 - np.log(np.exp(2*zeta_guess)/(2*np.exp(zeta_guess) + (Phi_bulk/Phi_max)*(1 - 2*np.exp(zeta_guess)+np.exp(2*zeta_guess)))**2))
#    	while per_dif>0.001:
#    		max_iter+=1
#		per_dif = -0.5 *(-2 - np.log(np.exp(2*zeta_guess)/(2*np.exp(zeta_guess) + (Phi_bulk/Phi_max)*(1 - 2*np.exp(zeta_guess)+np.exp(2*zeta_guess)))**2))
#		zeta_guess+=0.001 
#		if max_iter==1E6:
#			print 'Failure to converge!',fail
#	SigmaFail.append(np.sqrt(2*(Phi_max/Phi_bulk)*np.log(1 + 2*(Phi_bulk/Phi_max)*np.sinh(zeta_guess/2)**2)))
#	zetaS.append(zeta_guess)
#    zetaS=np.array(zetaS)
#    rhoPs_0 = 1/(2*np.exp(zetaS) + (Phi_bulkS/Phi_max)*(1-2*np.exp(zetaS) + np.exp(2*zetaS)))
#    rhoMs_0 = rhoPs_0*np.exp(2*zetaS)
#    CP_bik0 = -np.log(1 - (Phi_bulkS)*(rhoPs_0 + rhoMs_0 ))
#    CP_bik00 = -0.5*np.log(rhoPs_0 *rhoMs_0)    
#    print CP_bik0
#    print CP_bik00

    print '\n'
 ###    This plots SigmaFail vs Phi_bulk for Bikerman -- and eventually for data, etc.
    i=-1
    fig=plt.figure()
    fig.subplots_adjust(right=0.85)
    ax1=fig.add_subplot(111)
    graphname = 'P1_F1_Contour' + '.pdf'
    ax1.errorbar(Phi_bulkS,SigmaFailS,yerr=None,color='k',ls='-')
    ax1.set_xlabel(r'$\Phi^{bulk}$',size='x-large') 
    ax1.set_ylabel(r'$\tilde \Sigma_{fail}=4 \pi \Sigma_{fail} \lambda_{D} \lambda_{B} / q$',size='x-large') 
#    ax1.set_ylim(-0.05,0.05)
    plt.savefig(graphname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
  #    plt.show()
    plt.close()
    return Phi_bulkS,SigmaFailS
#Bikerman_Threshold(0.2)

def S_PM0(coeffs,App):
    """
    This calculated S_PM from a fitted function for rho_f(contact).

Input: 
      XXXXfilenameS - A vector of all the filenames to be included in the plots

Output: Plots as meantioned above.
Notes: Would be nice to extende Bikerman theory, too.
"""
    from scipy.integrate import quad

    if len(coeffs)==4:
        inverse_fitted_contact = lambda SIGMA: 1./(coeffs[0]*SIGMA**3+coeffs[1]*SIGMA**2+coeffs[2]*SIGMA+coeffs[3])
    elif len(coeffs)==3:
        inverse_fitted_contact = lambda SIGMA: 1./(coeffs[0]*SIGMA**2+coeffs[1]*SIGMA**1+coeffs[2])
    else:
    	print 'SPM() must be adapted for this fitting function'
    
    fSig = quad(inverse_fitted_contact, 1, App)[0]

#    inverse_GC_contact = lambda SIGMA: 1./ (SIGMA*np.sqrt((SIGMA/2.)**2+1))
#    fSig = quad(inverse_GC_contact, 1, App)[0]
#    print 'Calculating SGC - NOT SPM!!!'
#    print np.exp(-fSig),np.exp(-quad(inverse_GC_contact, 1, App)[0])
    return np.exp(-fSig)

def S_PM0_DH(coeffs,App):
    """
    This calculated S_PM from a fitted function for rho_f(contact).

Input: 
      XXXXfilenameS - A vector of all the filenames to be included in the plots

Output: Plots as meantioned above.
Notes: Would be nice to extende Bikerman theory, too.
"""
    from scipy.integrate import quad
    dlnh = lambda SIGMA: (1./SIGMA)-1./(coeffs[0]*SIGMA**2+coeffs[1]*SIGMA**1+coeffs[2]) if SIGMA>(1-coeffs[1])/coeffs[0] else 0  
    return App*np.exp(quad(dlnh, 0, App)[0])

def S_PM0_GC(coeffs,App):
    """
    This calculates S_PM(z=0) = SPM0 from a fitted function for rho_f(contact) via S_GC.

Input: 
      The applied surface charge density of the system for which SPM(z=0) must be calculated

Output: SPM(z=0), from which all SPM can be determined via SPM = SPM0*exp(-z).
"""
    from scipy.integrate import quad
    from scipy.optimize import fsolve
    
    contact_GC = lambda SIGMA: SIGMA*np.sqrt((SIGMA/2.)**2+1)
    SGC0 = lambda SIGMA: (1/SIGMA)*(np.sqrt(4+SIGMA**2)-2)
    if len(coeffs)==4:
        fitted_contact = lambda SIGMA: coeffs[0]*SIGMA**3+coeffs[1]*SIGMA**2+coeffs[2]*SIGMA+coeffs[3]
    elif len(coeffs)==3:
        fitted_contact = lambda SIGMA: coeffs[0]*SIGMA**2+coeffs[1]*SIGMA**1+coeffs[2]
    else:
    	print 'SPM() must be adapted for this fitting function'

    Integrand = lambda SIGMA:  1./contact_GC(SIGMA) -1./fitted_contact(SIGMA)
    
    First_root = fsolve(Integrand, 0.001)
    dlnh = lambda SIGMA: Integrand(SIGMA) if SIGMA>First_root else 0

    return SGC0(App)*np.exp(quad(dlnh, 0, App)[0])

def S_PM0_or_CS(coeffs,App,SigRhoS,pmorcs):
    """
    This calculates S_PM(z=0) = SPM0 from a fitted function for rho_f(contact) via S_CarnStarling.

Input: 
      The applied surface charge density of the system for which SPM(z=0) must be calculated

Output: SPM(z=0), from which all SPM can be determined via SPM = SPM0*exp(-z).
"""
    from scipy.integrate import quad
    from scipy.optimize import fsolve
    from scipy import interpolate

    if len(coeffs)==4:
        fitted_contact = lambda SIGMA: coeffs[0]*SIGMA**3+coeffs[1]*SIGMA**2+coeffs[2]*SIGMA+coeffs[3]
    elif len(coeffs)==3:
        fitted_contact = lambda SIGMA: coeffs[0]*SIGMA**2+coeffs[1]*SIGMA**1+coeffs[2]
    else:
    	print 'SPM() must be adapted for this fitting function'


    if pmorcs=='CS':
    	interp_SCSofSigma = interpolate.UnivariateSpline(SigRhoS[0],SigRhoS[2],k=5)
    	SCS_MMA = lambda SIGMA: interp_SCSofSigma(SIGMA)[0]
    	return SCS_MMA(App)
    	
    elif pmorcs=='GC':
    	contact_GC = lambda SIGMA: SIGMA*np.sqrt((SIGMA/2.)**2+1)
    	SGC0 = lambda SIGMA: (1/SIGMA)*(np.sqrt(4+SIGMA**2)-2)
    
        Integrand = lambda SIGMA:  1./contact_GC(SIGMA) -1./fitted_contact(SIGMA)

#        print '\n\n\nRoot GCPMCS solve trouble\n'
#    	for x in np.linspace(0.001,4):
#    		print x,'\t',Integrand(x)

        First_root = fsolve(Integrand, 0.001)
        
#        print 'First root,SIG = ',First_root,App
        dlnh = lambda SIGMA: Integrand(SIGMA) if SIGMA>First_root else 0
    
        return SGC0(App)*np.exp(quad(dlnh, 0, App)[0])    	
    else:
    	interp_FreeCSofSigma = interpolate.UnivariateSpline(SigRhoS[0],SigRhoS[1],k=5)
    	contact_CS = lambda SIGMA: interp_FreeCSofSigma(SIGMA)[0]

    	interp_SCSofSigma = interpolate.UnivariateSpline(SigRhoS[0],SigRhoS[2],k=5)
    	SCS_MMA = lambda SIGMA: interp_SCSofSigma(SIGMA)[0]
    
        Integrand = lambda SIGMA:  1./contact_CS(SIGMA) -1./fitted_contact(SIGMA)

        First_root = fsolve(Integrand, 0.001)
#        print 'First root,SIG = ',First_root,App
#        print 'FIX THIS CODE, maybr'

        dlnh = lambda SIGMA: Integrand(SIGMA) if SIGMA>First_root else 0
    
        return SCS_MMA(App)*np.exp(quad(dlnh, 0, App)[0])

def  P1_Plot():
    """
	Write a better description below!
	
    This generates superimposed plots of voltage(z), Electric field(z), n_+(Z), n_-(z), Cap(zeta_meas) all with Guoy-Chapman Theory.

Input: 
      XXXXfilenameS - A vector of all the filenames to be included in the plots

Output: Plots as meantioned above.
Notes: Would be nice to extende Bikerman theory, too.
"""
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from mpl_toolkits.mplot3d import Axes3D
    import sys

    from scipy import optimize
    
    
    #Special groups of filenames
    GC_sig1 = ['Analyzed_GC_800_0.0_1.0_1.0_500.0_44.75.txt','Analyzed_GC_802_0.0125663706144_1.0_1.0_500.0_44.75.txt','Analyzed_GC_802_0.0252584049349_1.0_1.0_500.0_44.75.txt','Analyzed_GC_806_0.0383902622269_1.0_1.0_500.0_44.75.txt','Analyzed_GC_808_0.0520876061965_1.0_1.0_500.0_44.75.txt','Analyzed_GC_816_0.0666645961092_1.0_1.0_500.0_44.75.txt','Analyzed_GC_820_0.082246895671_1.0_1.0_500.0_44.75.txt','Analyzed_GC_830_0.0990858322942_1.0_1.0_500.0_44.75.txt','Analyzed_GC_840_0.117495565244_1.0_1.0_500.0_44.75.txt','Analyzed_GC_864_0.160221225333_1.0_1.0_500.0_44.75.txt','Analyzed_GC_882_0.18510263915_1.0_1.0_500.0_44.75.txt','Analyzed_GC_908_0.21293715006_1.0_1.0_500.0_44.75.txt']
    GC_sig5 = ['Analyzed_GC_800_0.0_1.0_5.0_500.0_44.75.txt','Analyzed_GC_802_0.0125663706144_1.0_5.0_500.0_44.75.txt','Analyzed_GC_802_0.0253212367879_1.0_5.0_500.0_44.75.txt','Analyzed_GC_806_0.038515925933_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.0523389336088_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.0669787553745_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.0826238867894_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.0995884871188_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.118123883775_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.16097520757_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.186045116946_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.214005291563_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.253463695292_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.280418560259_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.364487579669_1.0_5.0_500.0_44.75.txt']#,'Analyzed_GC_808_0.608023842176_1.0_5.0_500.0_44.75.txt'] 
    GC_sig6 = ['Analyzed_GC_800_0.0_1.0_6.0_500.0_44.75.txt','Analyzed_GC_802_0.0125663706144_1.0_6.0_500.0_44.75.txt','Analyzed_GC_802_0.0253212367879_1.0_6.0_500.0_44.75.txt','Analyzed_GC_806_0.038515925933_1.0_6.0_500.0_44.75.txt','Analyzed_GC_808_0.0523389336088_1.0_6.0_500.0_44.75.txt','Analyzed_GC_808_0.0669787553745_1.0_6.0_500.0_44.75.txt','Analyzed_GC_808_0.0826238867894_1.0_6.0_500.0_44.75.txt','Analyzed_GC_808_0.0995884871188_1.0_6.0_500.0_44.75.txt','Analyzed_GC_808_0.118123883775_1.0_6.0_500.0_44.75.txt','Analyzed_GC_808_0.16097520757_1.0_6.0_500.0_44.75.txt','Analyzed_GC_808_0.17775131234_1.0_6.0_500.0_44.75.txt','Analyzed_GC_808_0.186045116946_1.0_6.0_500.0_44.75.txt','Analyzed_GC_808_0.214005291563_1.0_6.0_500.0_44.75.txt','Analyzed_GC_848_0.137790253786_1.0_6.0_500.0_44.75.txt','Analyzed_GC_848_0.14646104951_1.0_6.0_500.0_44.75.txt','Analyzed_GC_862_0.169834498853_1.0_6.0_500.0_44.75.txt','Analyzed_GC_872_0.190443346661_1.0_6.0_500.0_44.75.txt','Analyzed_GC_924_0.280418560259_1.0_6.0_500.0_44.75.txt']
    GC_sig7 = ['Analyzed_GC_800_0.0_1.0_7.0_500.0_44.75.txt','Analyzed_GC_802_0.0125663706144_1.0_7.0_500.0_44.75.txt','Analyzed_GC_802_0.0253212367879_1.0_7.0_500.0_44.75.txt','Analyzed_GC_806_0.038515925933_1.0_7.0_500.0_44.75.txt','Analyzed_GC_808_0.0523389336088_1.0_7.0_500.0_44.75.txt','Analyzed_GC_808_0.0669787553745_1.0_7.0_500.0_44.75.txt','Analyzed_GC_808_0.0826238867894_1.0_7.0_500.0_44.75.txt','Analyzed_GC_808_0.0995884871188_1.0_7.0_500.0_44.75.txt','Analyzed_GC_808_0.118123883775_1.0_7.0_500.0_44.75.txt','Analyzed_GC_808_0.134208838161_1.0_7.0_500.0_44.75.txt','Analyzed_GC_808_0.16097520757_1.0_7.0_500.0_44.75.txt','Analyzed_GC_808_0.186045116946_1.0_7.0_500.0_44.75.txt','Analyzed_GC_808_0.214005291563_1.0_7.0_500.0_44.75.txt','Analyzed_GC_840_0.125412378731_1.0_7.0_500.0_44.75.txt','Analyzed_GC_848_0.142062819795_1.0_7.0_500.0_44.75.txt','Analyzed_GC_862_0.169834498853_1.0_7.0_500.0_44.75.txt','Analyzed_GC_878_0.201438920948_1.0_7.0_500.0_44.75.txt']
    GC_sig8 = ['Analyzed_GC_800_0.0_1.0_8.0_500.0_44.75.txt','Analyzed_GC_802_0.0125663706144_1.0_8.0_500.0_44.75.txt','Analyzed_GC_802_0.0252584049349_1.0_8.0_500.0_44.75.txt','Analyzed_GC_806_0.0383902622269_1.0_8.0_500.0_44.75.txt','Analyzed_GC_808_0.0520876061965_1.0_8.0_500.0_44.75.txt','Analyzed_GC_816_0.0666645961092_1.0_8.0_500.0_44.75.txt','Analyzed_GC_818_0.082246895671_1.0_8.0_500.0_44.75.txt','Analyzed_GC_828_0.107756628018_1.0_8.0_500.0_44.75.txt','Analyzed_GC_838_0.117495565244_1.0_8.0_500.0_44.75.txt','Analyzed_GC_848_0.137790253786_1.0_8.0_500.0_44.75.txt','Analyzed_GC_858_0.160221225333_1.0_8.0_500.0_44.75.txt']
    GC_LBhalf = ['Analyzed_GC_800_0.0_0.5_1.0_500.0_31.64.txt','Analyzed_GC_802_0.00626747734391_0.5_1.0_500.0_31.64.txt','Analyzed_GC_802_0.0126292024674_0.5_1.0_500.0_31.64.txt','Analyzed_GC_806_0.0191951311134_0.5_1.0_500.0_31.64.txt','Analyzed_GC_808_0.0260595110615_0.5_1.0_500.0_31.64.txt','Analyzed_GC_816_0.0333322980546_0.5_1.0_500.0_31.64.txt','Analyzed_GC_818_0.0411077398722_0.5_1.0_500.0_31.64.txt','Analyzed_GC_824_0.0495429161471_0.5_1.0_500.0_31.64.txt','Analyzed_GC_838_0.0587634905854_0.5_1.0_500.0_31.64.txt','Analyzed_GC_858_0.0800949047033_0.5_1.0_500.0_31.64.txt','Analyzed_GC_870_0.0925513195748_0.5_1.0_500.0_31.64.txt','Analyzed_GC_884_0.10646857503_0.5_1.0_500.0_31.64.txt','Analyzed_GC_924_0.139518129746_0.5_1.0_500.0_31.64.txt','Analyzed_GC_980_0.181348435928_0.5_1.0_500.0_31.64.txt']
    GC_LB3 = ['Analyzed_GC_1128_1.81521223524_3.0_1.0_500.0_77.51.txt','Analyzed_GC_800_0.0_3.0_1.0_500.0_77.51.txt','Analyzed_GC_802_0.0378876074023_3.0_1.0_500.0_77.51.txt','Analyzed_GC_802_0.0757752148046_3.0_1.0_500.0_77.51.txt','Analyzed_GC_806_0.11535928224_3.0_1.0_500.0_77.51.txt','Analyzed_GC_808_0.15607432303_3.0_1.0_500.0_77.51.txt','Analyzed_GC_816_0.200182283887_3.0_1.0_500.0_77.51.txt','Analyzed_GC_818_0.246552191454_3.0_1.0_500.0_77.51.txt','Analyzed_GC_824_0.297445992442_3.0_1.0_500.0_77.51.txt','Analyzed_GC_838_0.352298200174_3.0_1.0_500.0_77.51.txt','Analyzed_GC_848_0.413370761359_3.0_1.0_500.0_77.51.txt','Analyzed_GC_858_0.480663675999_3.0_1.0_500.0_77.51.txt','Analyzed_GC_870_0.555307917449_3.0_1.0_500.0_77.51.txt','Analyzed_GC_884_0.63899994574_3.0_1.0_500.0_77.51.txt','Analyzed_GC_902_0.732305247552_3.0_1.0_500.0_77.51.txt','Analyzed_GC_924_0.836920282916_3.0_1.0_500.0_77.51.txt','Analyzed_GC_980_1.08799636779_3.0_1.0_500.0_77.51.txt']
    GC_LB5 = ['Analyzed_GC_800_0.0_5.0_1.0_500.0_103.14.txt','Analyzed_GC_802_0.061261056745_5.0_1.0_500.0_103.14.txt','Analyzed_GC_802_0.12252211349_5.0_1.0_500.0_103.14.txt','Analyzed_GC_806_0.186924762889_5.0_1.0_500.0_103.14.txt','Analyzed_GC_810_0.254469004941_5.0_1.0_500.0_103.14.txt','Analyzed_GC_812_0.325154839647_5.0_1.0_500.0_103.14.txt','Analyzed_GC_812_0.402123859659_5.0_1.0_500.0_103.14.txt','Analyzed_GC_812_0.483805268653_5.0_1.0_500.0_103.14.txt','Analyzed_GC_812_0.574911455607_5.0_1.0_500.0_103.14.txt','Analyzed_GC_812_0.782256570744_5.0_1.0_500.0_103.14.txt','Analyzed_GC_812_0.904778684234_5.0_1.0_500.0_103.14.txt','Analyzed_GC_812_1.03986716834_5.0_1.0_500.0_103.14.txt','Analyzed_GC_812_1.36345121166_5.0_1.0_500.0_103.14.txt','Analyzed_GC_812_1.77185825662_5.0_1.0_500.0_103.14.txt']#,'Analyzed_GC_812_2.95623868703_5.0_1.0_500.0_103.14.txt']
    GC_LB7 = ['Analyzed_GC_1128_4.23637486151_7.0_1.0_500.0_118.4.txt','Analyzed_GC_800_0.0_7.0_1.0_500.0_118.4.txt','Analyzed_GC_802_0.089284063215_7.0_1.0_500.0_118.4.txt','Analyzed_GC_802_0.17548936563_7.0_1.0_500.0_118.4.txt','Analyzed_GC_806_0.267852189645_7.0_1.0_500.0_118.4.txt','Analyzed_GC_808_0.363293774461_7.0_1.0_500.0_118.4.txt','Analyzed_GC_816_0.467971641679_7.0_1.0_500.0_118.4.txt','Analyzed_GC_818_0.575728269697_7.0_1.0_500.0_118.4.txt','Analyzed_GC_824_0.692721180117_7.0_1.0_500.0_118.4.txt','Analyzed_GC_838_0.822029133738_7.0_1.0_500.0_118.4.txt','Analyzed_GC_848_0.963652130562_7.0_1.0_500.0_118.4.txt','Analyzed_GC_858_1.12066893139_7.0_1.0_500.0_118.4.txt','Analyzed_GC_870_1.29615829702_7.0_1.0_500.0_118.4.txt','Analyzed_GC_884_1.49012022745_7.0_1.0_500.0_118.4.txt','Analyzed_GC_902_1.70871224429_7.0_1.0_500.0_118.4.txt','Analyzed_GC_924_1.95193434753_7.0_1.0_500.0_118.4.txt','Analyzed_GC_980_2.53997766043_7.0_1.0_500.0_118.4.txt']
    GC_LB10 = ['Analyzed_GC_800_0.0_10.0_1.0_500.0_141.51.txt','Analyzed_GC_802_0.125663706144_10.0_1.0_500.0_141.51.txt','Analyzed_GC_802_0.251327412287_10.0_1.0_500.0_141.51.txt','Analyzed_GC_806_0.383274303738_10.0_1.0_500.0_141.51.txt','Analyzed_GC_808_0.521504380496_10.0_1.0_500.0_141.51.txt','Analyzed_GC_816_0.666017642561_10.0_1.0_500.0_141.51.txt','Analyzed_GC_818_0.823097275241_10.0_1.0_500.0_141.51.txt','Analyzed_GC_824_0.992743278534_10.0_1.0_500.0_141.51.txt','Analyzed_GC_838_1.17495565244_10.0_1.0_500.0_141.51.txt','Analyzed_GC_848_1.37601758227_10.0_1.0_500.0_141.51.txt','Analyzed_GC_858_1.60221225333_10.0_1.0_500.0_141.51.txt','Analyzed_GC_870_1.85353966562_10.0_1.0_500.0_141.51.txt','Analyzed_GC_884_2.12999981913_10.0_1.0_500.0_141.51.txt','Analyzed_GC_902_2.43787589919_10.0_1.0_500.0_141.51.txt','Analyzed_GC_924_2.78973427639_10.0_1.0_500.0_141.51.txt','Analyzed_GC_980_3.62539792224_10.0_1.0_500.0_141.51.txt','Analyzed_GC_1128_6.05070745081_10.0_1.0_500.0_141.51.txt']
    IDEAL = ['Analyzed_GC_1009_0.431026512073_7.0_7.0_500.0_355.2.txt','Analyzed_GC_1051_0.495680488883_7.0_7.0_500.0_355.2.txt','Analyzed_GC_1178_0.649618528909_7.0_7.0_500.0_355.2.txt','Analyzed_GC_1316_0.846659220142_7.0_7.0_500.0_355.2.txt','Analyzed_GC_1768_1.41315120744_7.0_7.0_500.0_355.2.txt','Analyzed_GC_796_0.0_7.0_7.0_500.0_355.2.txt','Analyzed_GC_797_0.0307876080052_7.0_7.0_500.0_355.2.txt','Analyzed_GC_804_0.0584964552098_7.0_7.0_500.0_355.2.txt','Analyzed_GC_805_0.089284063215_7.0_7.0_500.0_355.2.txt','Analyzed_GC_825_0.12007167122_7.0_7.0_500.0_355.2.txt','Analyzed_GC_836_0.157016800826_7.0_7.0_500.0_355.2.txt','Analyzed_GC_851_0.190883169632_7.0_7.0_500.0_355.2.txt','Analyzed_GC_874_0.230907060039_7.0_7.0_500.0_355.2.txt','Analyzed_GC_898_0.274009711246_7.0_7.0_500.0_355.2.txt','Analyzed_GC_935_0.320191123254_7.0_7.0_500.0_355.2.txt','Analyzed_GC_972_0.372530056863_7.0_7.0_500.0_355.2.txt']
    GC_LB_LD_5_20=['Analyzed_GC_800_0.0_5.0_2.0_500.0_200.13.txt','Analyzed_GC_806_0.0628318530718_5.0_2.0_500.0_200.13.txt','Analyzed_GC_816_0.130376095124_5.0_2.0_500.0_200.13.txt','Analyzed_GC_838_0.20577431881_5.0_2.0_500.0_200.13.txt','Analyzed_GC_870_0.293738913111_5.0_2.0_500.0_200.13.txt','Analyzed_GC_914_0.400553063333_5.0_2.0_500.0_200.13.txt','Analyzed_GC_972_0.532499954783_5.0_2.0_500.0_200.13.txt','Analyzed_GC_1050_0.697433569097_5.0_2.0_500.0_200.13.txt','Analyzed_GC_1454_1.5126768627_5.0_2.0_500.0_200.13.txt']
    LB_LD_Ratio_2 = ['Analyzed_GC_1000_0.0_20.0_1.0_500.0_223.75.txt','Analyzed_GC_1000_0.251327412287_20.0_1.0_500.0_223.75.txt','Analyzed_GC_1000_0.0502654824574_20.0_1.0_500.0_223.75.txt','Analyzed_GC_1000_0.502654824574_20.0_1.0_500.0_223.75.txt','Analyzed_GC_1002_0.125663706144_20.0_1.0_500.0_223.75.txt','Analyzed_GC_1012_1.03044239038_20.0_1.0_500.0_223.75.txt','Analyzed_GC_1024_1.63362817987_20.0_1.0_500.0_223.75.txt','Analyzed_GC_1044_2.3624776755_20.0_1.0_500.0_223.75.txt','Analyzed_GC_1070_3.19185813605_20.0_1.0_500.0_223.75.txt']
	##This is not comprehensive with all of the latest results!
    sig_LD_3_over_15 = ['Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1312_0.0241789536991_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt']
    F3_Data=['Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1312_0.0241789536991_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt']
  
    filenameS=[]

    for f in [[x for x in line.split()] for line in file(sys.argv[1],"r").readlines()]:
    	filenameS.append(f[0])
#    	print 'submit "P1.py %s"\nsleep 2' % str(f[0])
#    filenameS.reverse()
    import random	
    random.shuffle(filenameS)


#    filenameS = GC_sig1

#    print GC_sig1[4:5],GC_LB7[3:5]
#    filenameS = GC_LB_LD_5_20+GC_sig1[4:5]+GC_LB7
#    filenameS = GC_sig1+GC_LB7
#    filenameS = GC_sig1[4:5]

#    filenameS = GC_LB3#  ['Analyzed_GC_1128_1.81521223524_3.0_1.0_500.0_77.51.txt','Analyzed_GC_800_0.0_3.0_1.0_500.0_77.51.txt','Analyzed_GC_802_0.0378876074023_3.0_1.0_500.0_77.51.txt','Analyzed_GC_802_0.0757752148046_3.0_1.0_500.0_77.51.txt','Analyzed_GC_806_0.11535928224_3.0_1.0_500.0_77.51.txt','Analyzed_GC_808_0.15607432303_3.0_1.0_500.0_77.51.txt','Analyzed_GC_816_0.200182283887_3.0_1.0_500.0_77.51.txt','Analyzed_GC_818_0.246552191454_3.0_1.0_500.0_77.51.txt','Analyzed_GC_824_0.297445992442_3.0_1.0_500.0_77.51.txt','Analyzed_GC_838_0.352298200174_3.0_1.0_500.0_77.51.txt','Analyzed_GC_848_0.413370761359_3.0_1.0_500.0_77.51.txt','Analyzed_GC_858_0.480663675999_3.0_1.0_500.0_77.51.txt','Analyzed_GC_870_0.555307917449_3.0_1.0_500.0_77.51.txt','Analyzed_GC_884_0.63899994574_3.0_1.0_500.0_77.51.txt','Analyzed_GC_902_0.732305247552_3.0_1.0_500.0_77.51.txt','Analyzed_GC_924_0.836920282916_3.0_1.0_500.0_77.51.txt','Analyzed_GC_980_1.08799636779_3.0_1.0_500.0_77.51.txt']

#    filenameS = GC_LB10

#    import random
#    test = filenameS[7:-1] 
#    
#    filenameS = filenameS[:7]+test+[filenameS[-1]]
#    for x in filenameS:
#    	print '\'',x,'\''

#    filenameS = GC_sig5[:-1]# + LB_LD_Ratio_2
#    filenameS =['Analyzed_GC_800_0.0_1.0_5.0_500.0_44.75.txt','Analyzed_GC_802_0.0125663706144_1.0_5.0_500.0_44.75.txt','Analyzed_GC_802_0.0253212367879_1.0_5.0_500.0_44.75.txt','Analyzed_GC_806_0.038515925933_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.0523389336088_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.0669787553745_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.0826238867894_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.0995884871188_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.118123883775_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.16097520757_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.186045116946_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.214005291563_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.280418560259_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.364487579669_1.0_5.0_500.0_44.75.txt'] 
    #This is GC_sig5 - select
    
#    filenameS+=LB_LD_Ratio_2
#    filenameS = LB_LD_Ratio_2[:-5]
#    filenameS = GC_sig1
#    filenameS = GC_sig5[:-1] + GC_sig1 + GC_LB5 + LB_LD_Ratio_2
#    filenameS= GC_sig1 + GC_LB5 + LB_LD_Ratio_2  
#    filenameS=GC_LB10
#    filenameS = GC_sig5

#    print '\n\n\n'
#    for f in filenameS:
#    	print f
 
    Volts=[]
    zminS=[]
    zmaxS=[]
    Fields=[]
    N_plusS=[]
    N_minusS=[]
    LDs=[]
    n0s=[]
    SIGMAS=[]
    BjerrumS=[]
    sigWCA_S=[]
    N_totS=[]
    XiS=[]
    L_zS=[]
    areaS=[]
    muexEV_S=[]
    muexEVtot_S=[]
    z_S=[]	#These are the z_pos for plotting is for voltage, field, and muex
    L_binS=[]
    PhiBulkS=[]
    sigHS_S=[]
    
    Bikerman='n0'

   #These variables below are (reasonably) assumed to always equal these values
    z_wall=1.	
    temperature = 1.0
    valency = 1.0
    eps_wall=1.0

    #This loop build all the necessary data to superimpose plots. 
    print 'Loading Data...'
    i=0
    for filename in filenameS:  	
    	print filename
    	i+=1

        #Extract relevant parameters from filename
        j=0
        value=''
        for letter in filename[12:]:
    		if letter!='_':
    			value+=letter
		else:
			j+=1
			if j==1:
				N_tot=int(value)
				N_totS.append(N_tot)
			if j==2:
				Xi=float(value)
				XiS.append(Xi)
			if j==3:
				Bjerrum=float(value)
				BjerrumS.append(Bjerrum)
			if j==4:
				sigWCA=float(value)
				if sigWCA==1.0:
#					print filename
					a=1 #?
				sigWCA_S.append(sigWCA)

				if sigWCA==0.5:
					sig_HS = 0.452400783508				
				elif sigWCA==1.0:
					sig_HS = 0.904801567015
				elif sigWCA==2.0:
					sig_HS = 1.80960313403
				elif sigWCA==3.0:
					sig_HS = 2.714440470105
				elif sigWCA==5.0:
					sig_HS = 4.52400783508
				elif sigWCA==5.4:
					sig_HS = 4.88592846188
				elif sigWCA==6.0:
					sig_HS = 5.42880940209
				elif sigWCA==7.0:
					sig_HS = 6.33361096911
				elif sigWCA==8.0:
					sig_HS = 7.23841253612
				elif sigWCA==9.6:
					sig_HS = 8.68609504334
				else:
					sig_HS = Noro_WCA(1,sigWCA, 10**-20, sigWCA, 10**5) #This should be used when needed!	
					if i==1:
						print '\t**Hard sphere diameters are being calculated. This takes time.'	
					print '(WCA,HS) = (%1.3f,%1.3f)' % (sigWCA,sig_HS)
					print filename,sig_HS
				sigHS_S.append(sig_HS)	
			if j==5:#This length will actually be determined from LAMMPS files
				test2=0 #The L_z that includes wall thickness is more accurately obtained from z_positions below
				
				##Keep these uncommented and eventually delete - 03/07/12 14:34:06 
#				L_z=float(value) #This length DOES NOT include the walls... - 03/07/12 13:02:14 
#				L_zS.append(L_z) 
			value=''
	L_xy=float(value[:-4])
	areaS.append(L_xy*L_xy)

    	z_50=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()])
    
	#Trim incoming data
	z_positions=z_50[:,0]
        SIGMAS.append(z_positions[-1])  
        L_bin=z_positions[1]-z_positions[0]
        L_binS.append(L_bin)
        zpos = z_positions[:-1].tolist()
        z_S.append(zpos)
        L_zS.append(zpos[-1])
        
	Psi_tot=z_50[:,1]
	zminS.append(Psi_tot[-1])
	Volts.append(Psi_tot[:-1].tolist())
	
	E_tot=z_50[:,2]
	zmaxS.append(E_tot[-1])
	Fields.append(E_tot[:-1].tolist())

	N_plus=z_50[:,3]	
	N_plusS.append(N_plus[:-2])

	N_minus=z_50[:,4]
	lam_D=N_minus[-2]
#	print filename,lam_D
	LDs.append(N_minus[-2])
	n0 = (8.*np.pi*Bjerrum*lam_D**2)**-1
	n0s.append(n0)
	N_minusS.append(N_minus[:-2])

	PhiBulkS.append(2.*n0*(np.pi/6)*sig_HS**3)

	muex_EV=z_50[:,5]
	muexEVtot_S.append(muex_EV[-1])
#	print 'muexEVtot_S has been re-defined as average z_inside_wall position\n\t\tNeed to make changes!'
	muexEV_S.append(muex_EV[:-1].tolist())


    ##Turns out Np and Nm notation is somehow switched... This is fixed below:
    temp = N_plusS
    N_plusS = N_minusS
    N_minusS = temp

#    print r"\begin{tabular}{|c|c|c|c|}"
#    print '\hline'
#    print r"$\Sigma_{app}$ & $\lambda_{D}$ & $N_{+}+N{-}$ \zeta_{+}^{measured}\times k_{B}T/q_{\pm} \\ \hline"
##    old_ratio = sigWCA_S[0]/BjerrumS[0]
#    for (Sig_s,lam_D,sig,lam_B,N_p,N_m,V) in zip(SIGMAS,LDs,sigWCA_S,BjerrumS,N_plusS,N_minusS,Volts):
##      if old_ratio !=sig/lam_B:
##      		print '\n'
##      		old_ratio=old_ratio
#      print lam_B,'\t',r"%1.4f & %1.3f & %i & %1.2f \\ \hline" % (Sig_s,lam_D,int(round(np.sum(N_p)+np.sum(N_m))),V[0])


    #for (Sig_s,zL,zR,dp,dp_sq,N_p,N_m,n0,lam_D) in zip(SIGMAS,zeta_dataL,zeta_dataR,dipole_avgs,dipole_avg_sqs,N_plusS,N_minusS,n0s,LDs):
      #print '%1.5f\t\t%1.3f\t\t%1.3f\t\t%1.4f\t\t%1.3f\t\t%i\t\t%1.3f\t\t%1.3f\t\t%1.2f\t\t%1.2f\t\t%i\t\t%1.3f\t\t%1.2f\t\t%i' % (Sig_s,zL,zR,dp,dp_sq,int(round(np.sum(N_p)+np.sum(N_m))),n0,lam_D,Bjerrum,sig,L_z,area,Bjerrum**-1,1)

#    print '*******************\n**********************\nPhi RANGE is: ',np.min(PhiBulkS)*(sigWCA/sig_HS)**3,np.max(PhiBulkS)*(sigWCA/sig_HS)**3

    print '***********\n***********\tPhiBulks(sigma*0.954)\t***********\n',np.array(PhiBulkS)*(sigWCA*0.954/sig_HS)**3,'\n Mean:',np.mean(np.array(PhiBulkS)*(sigWCA*0.954/sig_HS)**3)
    quik= np.array(PhiBulkS)*(sigWCA*0.954/sig_HS)**3
    print np.mean(quik)-np.min(quik),np.mean(quik)-np.max(quik)



##    print PhiBulkS[2:]
#    print len(PhiBulkS)
	
        ##For all plots: 
#    characteristic_lengthS,xlabel=sigWCA_S,r'$z/ \sigma_{WCA}$'
#    characteristic_lengthS,xlabel=L_zS,r'$z/L_{z}$'
    characteristic_lengthS,xlabel=LDs,r'$z/\lambda_{D}$'

#    print '\n\n'
#    for (f,lb,ld,s) in zip(filenameS,BjerrumS,LDs,sigWCA_S):
#    	print f,lb,ld,s

#    print '\n\n'
#    for (f,Sig) in zip(filenameS,SIGMAS):
#    	print Sig,'\t',f
    markers=['v','o','>','s','^','d','<','*','3','D','4','p','h','1','H','2']*15
    colors=[]
    for step in range(0,len(filenameS)):
#      colors.append(ROYGBIV_map(float(step)/len(filenameS),1))
      colors.append(ROYGBIV_map(float(step),len(filenameS)))
#    colors.reverse()

#	##Turn this off for Figure 2 and 4, I think
    colors=[]
    test=[]
    for (Xi,lam_D,Bjerrum,filename) in zip(XiS,LDs,BjerrumS,filenameS):
    		dielectric=(4*np.pi*Bjerrum)**-1
    		print filename,round(Xi/(2.*np.pi*Bjerrum**2)/(dielectric*1.0/(1.0*lam_D)),1),lam_D
#    		if Bjerrum==20.0:
#    	    		test.append([round(Xi/(2.*np.pi*Bjerrum**2)/(dielectric*1.0/(1.0*lam_D)),1),4,1])
#    		else:
    	    	test.append([round(Xi/(2.*np.pi*Bjerrum**2)/(dielectric*1.0/(1.0*lam_D)),1),17,1])
    for t in test:
#    	  print t
          colors.append(ROYGBIV_map(t[0],t[1],t[2]))
#    colors.reverse()


#    print 'F1 Inset...'
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    fig.subplots_adjust(right=0.98) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.11)##Larger adds whitespace
#    fig.subplots_adjust(left=0.12) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.95) ##Smaller adds whitespace to top
#    
#    F1_LD = 1.5
#    F1_Si = 10.
#    ax1.plot(np.linspace(0,10*F1_LD,100),(F1_Si/F1_LD)*np.exp(-np.linspace(0,10*F1_LD,100)/F1_LD),color='k',ls='-',ms=0)  

#    F1_LD = 0.25
#    F1_Si = 10.
#    ax1.plot(np.linspace(0,3*10*F1_LD,100),(0.25/1.5)*(F1_Si/F1_LD)*np.exp(-np.linspace(0,10*F1_LD,100)/F1_LD),color='blue',ls='-',ms=0)  

#    i=-1
#    i_max=10000
##    plt.text(0.2, 1.0, r'$Decreasing$' +' '+ r'$|\Sigma|\longrightarrow$', fontsize='small')
##    for z in np.linspace(0,10*F1_LD,i_max):
##    for colorp in colors:
##    	i+=1
###	ax1.plot(z,(F1_Si/F1_LD)*np.exp(-z/F1_LD),color=ROYGBIV_map((F1_Si/F1_LD**2)*np.exp(-z/F1_LD),(F1_Si/F1_LD**2)*1.5,1),ms=7,marker='.')  
###	colorp = ROYGBIV_map((F1_Si/F1_LD**2)*np.exp(-z/F1_LD),(F1_Si/F1_LD**2)*1.5,1)
###	ax1.plot(z+0.3,1,color=colorp,ms=10,marker='s',mec=colorp)  
###	ax1.plot(z+0.3,1.1,color=colorp,ms=10,marker='s',mec=colorp)    	
###	ax1.plot(z+0.3,1.2,color=colorp,ms=10,marker='s',mec=colorp)    	
###	ax1.plot(z+0.3,1.3,color=colorp,ms=10,marker='s',mec=colorp)    	
###	ax1.plot(z+0.3,1.4,color=colorp,ms=10,marker='s',mec=colorp)    	
##	ax1.plot(z+0.3,1,color=colorp,ms=10,marker='s',mec=colorp)  
##	ax1.plot(z+0.3,1.1,color=colorp,ms=10,marker='s',mec=colorp)    	
##	ax1.plot(z+0.3,1.2,color=colorp,ms=10,marker='s',mec=colorp)    	
##	ax1.plot(z+0.3,1.3,color=colorp,ms=10,marker='s',mec=colorp)    	
##	ax1.plot(z+0.3,1.4,color=colorp,ms=10,marker='s',mec=colorp)    	

##    plt.text(0.2, 0.96, r'$Decreasing$' +' '+ r'$|\Sigma|\longrightarrow0$', fontsize='small',weight='bold')

##    plt.text(0.2, 0.96, r'$Decreasing$' +' '+ r'$\Sigma_{app}\longrightarrow$', fontsize='small',weight='bold')


##    plt.text(0.2, 0.96, r'$Mean-field$' +' '+ r'$predictions$'+ ' '+ r'$for$'+ ' '+ r'$\rho_{f}$'+ ' '+ r'$are$'+ ' ' + r'$similar.$' , fontsize='small',weight='bold')
#    
##    xlabel= r'$\mathrm{Towards}$' +' '+ r'$the$' + ' ' + r'$bulk$' + ' '  r'$z\longrightarrow$'
##    ylabel= r'$Increasing$' +' '+ r'$|\rho_{f}|\longrightarrow$'

#    xlabel= r'$\mathrm{Distance}$' +' '+ r'${\rm from}$' + ' ' + r'${\rm plate,}$' + ' '+ r'$z$'
#    ylabel= r'${\rm Free}$' +' '+r'${\rm charge}$' +' '+ r'${\rm density,}$'+ ' '+ r'$|\rho_{f}|$'
#    
#    ax1.set_xlabel(xlabel,size=14)
#    ax1.set_ylabel(ylabel,size=14)
##    ax1.set_ylim(-0.1*F1_Si/F1_LD,F1_Si/F1_LD)
##    ax1.set_xlim(-0.25*F1_LD,10*F1_LD) 
#    ax1.set_ylim(-0.1*(0.25/1.5)*F1_Si/F1_LD,1.1*(0.25/1.5)*F1_Si/F1_LD)
#    ax1.set_xlim(-0.25*F1_LD,2.5*10*F1_LD) 

#    plt.grid(False)
#    plt.setp(plt.gca(), 'yticklabels', [])
#    plt.setp(plt.gca(), 'xticklabels', [])
#    frame1 = plt.gca()
#    frame1.axes.get_xaxis().set_ticks([])
#    frame1.axes.get_yaxis().set_ticks([])

#    fig.set_size_inches(1.35*2,1.25*2)
#    plt.savefig('F1_Inset_bar.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight') 

##    ax1.legend(loc=0)
##    plt.savefig('F1_Inset_legend.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 
##    plt.show()
#    plt.close()




#    print 'Free shifted...'
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    
#    LD = 1.5
#    blue_zeta = 3.
#    distance = np.linspace(0,10*LD,100)
#    ax1.plot(distance,blue_zeta*np.exp(-distance/LD),color='blue',ls='-',ms=0,lw = 2.)  
#    
#    red_zeta = 1.
#    ax1.plot(distance,red_zeta*np.exp(-distance/LD),color='red',ls='-',ms=0,lw = 2.)  
#    ax1.plot(distance+1.648,red_zeta*np.exp(-distance/LD),color='red',ls='--',ms=0,lw = 2.1)

#    ax1.set_xlabel(r'$z$',size=18)
#    ax1.set_ylabel(r'$|\rho_{f}|$',size=18)

##    ax1.set_ylim(-0.1*(0.25/1.5)*F1_Si/F1_LD,1.1*(0.25/1.5)*F1_Si/F1_LD)
##    ax1.set_xlim(-0.25*F1_LD,2.5*10*F1_LD) 

#    plt.grid(False)
#    plt.setp(plt.gca(), 'yticklabels', [])
#    plt.setp(plt.gca(), 'xticklabels', [])
#    frame1 = plt.gca()
#    frame1.axes.get_xaxis().set_ticks([])
#    frame1.axes.get_yaxis().set_ticks([])
#    ax1.set_xlim(xmax=5*LD)
#    ax1.set_ylim(ymin=-0.05)
#    ax1.set_ylim(ymax=3.1)
#    
#    fig.set_size_inches(3.37,2.45)
#    plt.savefig('F1_shift.png', dpi=1200, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight') 
#    plt.close()

#    print fail




    print 'Plotting Voltages & Assembling Theory...'
    nonDim='yes'  #This is actually important for the RMS deviation calculation below
#    Bikerman='yes'
    Bikerman='no'
#    CarnahanStarling='yes'
    CarnahanStarling='no'    
#    MPBoltz='yes'
    MPBoltz='no'        
    if Bikerman=='yes':
    	print '\t...with Bikerman theory...'
    if CarnahanStarling=='yes':
    	print '\t...with Carnahan-Starling theory...'
    if MPBoltz=='yes':
    	print '\t...with Modified Poisson w/Boltzmann distribution theory...'
    #This plots voltages, measures zeta, and prepares GC Field and GC density theorys
    zeta_dataL,zeta_dataR=[],[]
    E_theorys=[]
    volt_theoryS=[]
    P_Theory=[]
    M_Theory=[]
    rhof_Theory=[]   
    i=-1
    fig=plt.figure()
    fig.subplots_adjust(right=0.85)
    ax1=fig.add_subplot(111)
    for (V,Sigma_s,lam_D,L_z,n0,z_positions,Bjerrum,characteristic_length,filename,Xi) in zip(Volts,SIGMAS,LDs,L_zS,n0s,z_S,BjerrumS,characteristic_lengthS,filenameS,XiS):
    	i+=1
	dielectric=(4*np.pi*Bjerrum)**-1
    	if i==len(colors): #This resets the colors
    		i=0
	if i==0:
		graphname = 'Voltages'
    		
	z_positions=[z/L_z for z in z_positions] #NonDim distances
	
	v_T_cat,v_T_ano=[],[]
	zeta_measL=V[0]
	zeta_measR=V[-1]		

#	print zeta_measL,zeta_measR
	##This code is used to catch erroneously analyzed files
	if zeta_measL>10.:
		print filename
#		print 'submit "P1.py %s.gz"\nsleep 2' % filename[12:-4]

	zeta_dataL.append(zeta_measL)
	zeta_dataR.append(zeta_measR)

	sign=1
	E_theoryL,E_theoryR=[],[]
	
#	###########################
#	## It is unclear if E_field is working. So:
#	##   either the P1__Analysis() code is wrong 
#	##   the theory
#	##   or nothing at all, I haven't needed to look at E_field's yet
#        ##			but they should be good so I'm not worried. - 02/28/12 18:54:06 
#	###########################

	ztheory=np.linspace(0,1,3000)
	for z in ztheory[:len(ztheory)/2]:
    		v_T_cat.append(-np.log(GC_density(z,lam_D/L_z,zeta_measL,n0,sign)/n0))
    		v_T_ano.append(-np.log(GC_density(z,lam_D/L_z,zeta_measR,n0,sign)/n0))
    		E_theoryL.append((-4*np.exp(-z/(lam_D/L_z))*np.tanh(zeta_measL/4))/((lam_D/L_z)*(np.exp(2*z/(lam_D/L_z))-np.tanh(zeta_measL)**2)))
    		E_theoryR.append((-4*np.exp(-z/(lam_D/L_z))*np.tanh(zeta_measR/4))/((lam_D/L_z)*(np.exp(2*z/(lam_D/L_z))-np.tanh(zeta_measR)**2)))
	if len(ztheory) % 2: #If there are an even number of z_positions
		v_T_cat.append(-np.log(GC_density(z_positions[len(z_positions)/2],lam_D/L_z,zeta_measL,n0,sign)/n0))
		z=ztheory[len(ztheory)/2]
		E_theoryL.append((-4*np.exp(-z/(lam_D/L_z))*np.tanh(zeta_measL/4))/((lam_D/L_z)*(np.exp(2*z/(lam_D/L_z))-np.tanh(zeta_measL)**2)))
	v_T_ano.reverse()
	volt_theory=v_T_cat+v_T_ano
	volt_theory=-np.array(volt_theory)
	E_theoryR.reverse()
	E_theorys.append(E_theoryL+E_theoryR)

#	##Original, working theory - 07/19/12 15:21:46 
#	for z in z_positions[:len(z_positions)/2]:
#    		v_T_cat.append(-np.log(GC_density(z,lam_D/L_z,zeta_measL,n0,sign)/n0))
#    		v_T_ano.append(-np.log(GC_density(z,lam_D/L_z,zeta_measR,n0,sign)/n0))
#    		E_theoryL.append((-4*np.exp(-z/(lam_D/L_z))*np.tanh(zeta_measL/4))/((lam_D/L_z)*(np.exp(2*z/(lam_D/L_z))-np.tanh(zeta_measL)**2)))
#    		E_theoryR.append((-4*np.exp(-z/(lam_D/L_z))*np.tanh(zeta_measR/4))/((lam_D/L_z)*(np.exp(2*z/(lam_D/L_z))-np.tanh(zeta_measR)**2)))
#	if len(z_positions) % 2: #If there are an even number of z_positions
#		v_T_cat.append(-np.log(GC_density(z_positions[len(z_positions)/2],lam_D/L_z,zeta_measL,n0,sign)/n0))
#		z=z_positions[len(z_positions)/2]
#		E_theoryL.append((-4*np.exp(-z/(lam_D/L_z))*np.tanh(zeta_measL/4))/((lam_D/L_z)*(np.exp(2*z/(lam_D/L_z))-np.tanh(zeta_measL)**2)))
#	v_T_ano.reverse()
#	volt_theory=v_T_cat+v_T_ano
#	volt_theory=-np.array(volt_theory)
#	E_theoryR.reverse()
#	E_theorys.append(E_theoryL+E_theoryR)

	rhoP=np.exp(-volt_theory)
	rhoM=np.exp(volt_theory)
    	P_Theory.append(rhoP)
    	M_Theory.append(rhoM)
	rhof_Theory.append([(-x+y)/2. for (x,y) in zip(rhoM,rhoP)]) #This 1/2. is required for the updated nondimensionalization scheme of rhoF - 05/30/12 16:40:42 

	#Finally plot the Voltage
#	print filename,np.array(z_positions[:3])*(L_z/characteristic_length),np.array(z_positions[-3:])*L_z
	ax1.errorbar(np.array(z_positions)*(L_z/characteristic_length),V,yerr=None,marker='+',ms=7.0,color=colors[i],ls='None')#,label=r'$\~ \Sigma$'+' = ' + str(round(Xi/(2.*np.pi*Bjerrum**2)/(dielectric*temperature/(valency*lam_D)),3)))

####		##All this code is legit, but the MATLAB files need a new naming convention... easy enough... - 02/21/12 11:55:56 
####	if Bikerman=='yes':
####	   Bik_file='_Bik_zeta_' + ze + '.txt'
####	   x_Bik=[[float(x) for x in line.split()] for line in file('x' + Bik_file,"r").readlines()]
####	   volt_Bik=[[float(x) for x in line.split()] for line in file('volt' + Bik_file,"r").readlines()]
####	   volt_Bik=np.array(volt_Bik[0])
####	   volt_Bik_bulk=np.mean(volt_Bik[len(volt_Bik)/2-2:len(volt_Bik)/2+2])
####	   ax1.plot(np.array(x_Bik[0])*(lam_D/characteristic_length),volt_Bik-volt_Bik_bulk,color=colors[i],lw=1.5,ls='--')
####	if CarnahanStarling=='yes':
####	   CS_file='_CS_zeta_' + ze + '.txt'
####	   x_CS=[[float(x) for x in line.split()] for line in file('x' + CS_file,"r").readlines()]
####	   volt_CS=[[float(x) for x in line.split()] for line in file('volt' + CS_file,"r").readlines()]
####	   volt_CS=np.array(volt_CS[0])
####	   volt_CS_bulk=np.mean(volt_CS[len(volt_CS)/2-2:len(volt_CS)/2+2])
####	   ax1.plot(np.array(x_CS[0])*(lam_D/characteristic_length),volt_CS-volt_CS_bulk,color=colors[i],lw=1.5,ls='-.')
####	if MPBoltz=='yes':
####	   MPBoltz_file='_MPBoltz_zeta_' + ze + '.txt'
####           x_MPBoltz=[[float(x) for x in line.split()] for line in file('x' + MPBoltz_file,"r").readlines()]      		
####	   volt_MPBoltz=[[float(x) for x in line.split()] for line in file('volt' + MPBoltz_file,"r").readlines()]
####	   volt_MPBoltz=np.array(volt_MPBoltz[0])
####	   volt_MPBoltz_bulk=np.mean(volt_MPBoltz[len(volt_MPBoltz)/2-2:len(volt_MPBoltz)/2+2])
####	   ax1.plot(np.array(x_MPBoltz[0])*(lam_D/characteristic_length),volt_MPBoltz-volt_MPBoltz_bulk,color=colors[i],lw=2.0,ls=':')
	volt_theoryS.append(volt_theory)
	ax1.errorbar(np.array(ztheory)*(L_z/characteristic_length),volt_theory,yerr=None,ls='-',color=colors[i],label=r'$\~ \Sigma$'+' = ' + str(round(Xi/(2.*np.pi*Bjerrum**2)/(dielectric*temperature/(valency*lam_D)),3)))

#    ax1.set_xlim(0,10)       
#    ax1.set_ylim(0,6)      
    ax1.set_xlabel(xlabel,size='x-large')
#    ax1.set_ylabel(r'$\~\psi$',size='x-large')
#    ax1.legend(loc=0) 

    print 'Voltage plot was not saved to a file!' 
#    plt.savefig(graphname + '.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 
#    plt.show()
    plt.close()
    zeta_data=[(abs(zL)+abs(zR))/2. for (zL,zR) in zip(zeta_dataL,zeta_dataR)]  #Do NOT use np.mean for these...
	#Done plotting voltages    
    print '\n\n\n********\n',np.min(zeta_data),np.max(zeta_data)
	
#  ###    This plots N+(z) for simulation and GC theory
#    nonDim='no'
##    nonDim='yes'
#    #phi='yes'
#    phi='no'
##    Bikerman='yes'
#    Bikerman='no'
##    CarnahanStarling='yes'
#    CarnahanStarling='no'
##    MPBoltz='yes'
#    MPBoltz='no' 
#    if nonDim=='yes' and phi=='yes':
#      print 'Plotting NonDim N_+(z) with volume fractions...'
#    elif nonDim=='yes' and phi!='yes':
#      print 'Plotting NonDim N_+(z)...'
#    elif nonDim!='yes' and phi!='yes':
#      print 'Plotting Dim N_+(z)...'	
#    else:
#      print 'Plotting Dim N_+(z) with volume fractions...'	
#    if Bikerman=='yes':
#      print '\t\t...with Bikerman theory...'
#    if CarnahanStarling=='yes':
#	    print '\t...with Carnahan-Starling theory...'
#    if MPBoltz=='yes':
#    	print '\t...with Modified Poisson w/Boltzmann distribution theory...'
#    i=-1
#    fig=plt.figure()
#    fig.subplots_adjust(right=0.85)
#    ax1=fig.add_subplot(111)
#    for (Np,P_T,n0,Sigma_s,area,lam_D,L_bin,z_positions,L_z,characteristic_length,Xi) in zip(N_plusS,P_Theory,n0s,SIGMAS,areaS,LDs,L_binS,z_S,L_zS,characteristic_lengthS,XiS):
#	      i+=1
#	      dielectric=(4*np.pi*Bjerrum)**-1
#	      if i==len(markers): #This resets the markers
#		      i=0	      		
#	      if nonDim=='yes':
#		      Np=[npl/(area*L_bin*n0) for npl in Np]
#		      graphname = 'N_plus_NDim' + '.pdf'
#		      ylabel=r'$\~N_{+}$'
#	      else: #This is dimensional
##		      Np=[npl/(area*L_bin) for npl in Np]
#		      P_T=[n0*t for t in P_T]
#		      graphname = 'N_plus_Dim' + '.pdf'
#		      ylabel=r'$N_{+}$'	

#		#This could surely be consolidated within the lines of code below...
#	      z_density = [(x+y)/(L_z*2.) for (x,y) in zip(z_positions[0:len(z_positions)-1],z_positions[1:len(z_positions)])]
#	      z_den_theory=[z/L_z for z in z_positions]  

#	      ax1.errorbar(np.array(z_density)*(L_z/characteristic_length),Np,yerr=None,marker=markers[i],ms=5.0,color=colors[i],ls='None',label=r'$\~ \Sigma$'+' = ' + str(round(Xi/(2.*np.pi*Bjerrum**2)/(dielectric*temperature/(valency*lam_D)),3)))
##	      ax1.errorbar(np.array(z_den_theory)*(L_z/characteristic_length),P_T,yerr=None,ls='-',color=colors[i])#,label=r'$\Sigma$'+' = ' + str(round(Sigma_s,3)))

#		##All this code is legit, but the MATLAB files need a new naming convention... easy enough... - 02/21/12 11:55:56 	      
#	      if Bikerman=='yes':
#		Bik_file='_Bik_zeta_' + ze + '.txt'
#		x_Bik=[[float(x) for x in line.split()] for line in file('x' + Bik_file,"r").readlines()]
#		co_Bik=[[float(x) for x in line.split()] for line in file('co' + Bik_file,"r").readlines()]
#		ax1.plot(np.array(x_Bik[0])*(lam_D/characteristic_length),co_Bik[0],color=colors[i],lw=1.5,ls='--')
#	      if CarnahanStarling=='yes':
#		CS_file='_CS_zeta_' + ze + '.txt'
#	        x_CS=[[float(x) for x in line.split()] for line in file('x' + CS_file,"r").readlines()]
#		co_CS=[[float(x) for x in line.split()] for line in file('co' + CS_file,"r").readlines()]
#		ax1.plot(np.array(x_CS[0])*(lam_D/characteristic_length),co_CS[0],color=colors[i],lw=1.5,ls='-.')	
#	      if MPBoltz=='yes':
#		MPBoltz_file='_MPBoltz_zeta_' + ze + '.txt'
#                x_MPBoltz=[[float(x) for x in line.split()] for line in file('x' + MPBoltz_file,"r").readlines()]      		
#		volt_MPBoltz=[[float(x) for x in line.split()] for line in file('volt' + MPBoltz_file,"r").readlines()]
#		volt_MPBoltz=np.array(volt_MPBoltz[0])
#		volt_MPBoltz_bulk=np.mean(volt_MPBoltz[len(volt_MPBoltz)/2-2:len(volt_MPBoltz)/2+2])
#		if nonDim=='yes':
#		   	ax1.plot(np.array(x_MPBoltz[0])*(lam_D/characteristic_length),np.exp(-volt_MPBoltz+volt_MPBoltz_bulk),color=colors[i],lw=2.0,ls=':')
#		elif nonDim=='no':
#		   	ax1.plot(np.array(x_MPBoltz[0])*(lam_D/characteristic_length),n0*np.exp(-volt_MPBoltz+volt_MPBoltz_bulk),color=colors[i],lw=2.0,ls=':')		
#    locs3=[] # A stupid name, but actually the words that print out on the right y-axis
#    if phi=='yes':	#This most likely will NOT work... - 02/21/12 13:56:36 
#	      locs1,labels=plt.yticks()
#	      for y in [x*(4.*np.pi/3.)*(sig*0.5)**3 for x in locs1]:
#		      s='%1.1E' % y
#		      locs3.append(s)  
#	      locs3.reverse()
#	      ax2=ax1.twinx()  	
#	      locs,labels=plt.yticks()  #These are the location of the new axis
#	      ax2.grid(True)
#	      plt.ylim(ymin=0)
#	      locs = [(x-min(locs1))/(max(locs1)-min(locs1)) for x in locs1] #What a pain in my ass to figure this shit out...
#	      locs.reverse()
#	      plt.yticks(locs, locs3,rotation=15)
#	      ax2.set_ylabel(r'$\Phi_{\sigma}$',size='x-large') 
#    if xlabel==r'$z/L_{z}$':
#      #These must be set manually! 
#      test=0
##      ax1.set_xlim(0,0.15)       
#      ax1.set_ylim(0,2) 
#    elif xlabel==r'$z/ \sigma_{WCA}$':
#      #These must be set manually! 
#      test=0
##      ax1.set_xlim(0,10)       
##      ax1.set_ylim(0,10)      
#    ax1.set_xlabel(xlabel,size='x-large')
#    ax1.set_ylabel(ylabel,size='x-large')
##    ax1.legend(loc=0) 
#    plt.savefig(graphname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    plt.close()
##    plt.show()
##	###Finished plotting N_plus

#  ###    This plots N-(z) for simulation and GC theory
#    nonDim='no'
##    nonDim='yes'
#    phi='no'
##    phi='yes'
##    Bikerman='yes'
#    Bikerman='no'
##    CarnahanStarling='yes'
#    CarnahanStarling='no'
##    MPBoltz='yes'
#    MPBoltz='no' 
#    if nonDim=='yes' and phi=='yes':
#      print 'Plotting NonDim N_-(z) with volume fractions...'
#    elif nonDim=='yes' and phi!='yes':
#      print 'Plotting NonDim N_-(z)...'
#    elif nonDim!='yes' and phi!='yes':
#      print 'Plotting Dim N_-(z)...'	
#    else:
#      print 'Plotting Dim N_-(z) with volume fractions...'	
#    if Bikerman=='yes':
#		print '\t\t...with Bikerman theory...'
#		#Note to user: Bikerman theory must already be evaluated using MATLAB code, which is in ~/sims/MATLAB_Solvers/PBik_Solver_zLzR.m
#    if CarnahanStarling=='yes':
#	      print '\t...with Carnahan-Starling theory...'
#    if MPBoltz=='yes':
#    	print '\t...with Modified Poisson w/Boltzmann distribution theory...'
#    i=-1
#    fig=plt.figure()
#    fig.subplots_adjust(right=0.85)
#    ax1=fig.add_subplot(111)
#    for (Np,P_T,n0,Sigma_s,area,lam_D,L_bin,z_positions,L_z,characteristic_length,filename,Xi) in zip(N_minusS,M_Theory,n0s,SIGMAS,areaS,LDs,L_binS,z_S,L_zS,characteristic_lengthS,filenameS,XiS):
#    		##This is exactly as the N+(z) code above, only with the zip(...,...) changed. It was easy and successful!
#	      i+=1
#	      dielectric=(4*np.pi*Bjerrum)**-1
#	      if i==len(markers): #This resets the markers
#		      i=0	      		

#	      if nonDim=='yes':
#		      Np=[npl/(area*L_bin*n0) for npl in Np]
#		      graphname = 'N_minus_NDim' + '.pdf'
#		      ylabel=r'$\~N_{-}$'
#	      else: #This is dimensional
#		      Np=[npl/(area*L_bin) for npl in Np]
#		      P_T=[n0*t for t in P_T]
#		      graphname = 'N_minus_Dim' + '.pdf'
#		      ylabel=r'$N_{-}$'	

#		#This could surely be consolidated within the lines of code below...
#	      z_density = [(x+y)/(L_z*2.) for (x,y) in zip(z_positions[0:len(z_positions)-1],z_positions[1:len(z_positions)])]
#	      z_den_theory=[z/L_z for z in z_positions]  

#	      ax1.errorbar(np.array(z_density)*(L_z/characteristic_length),Np,yerr=None,marker=markers[i],ms=5.0,color=colors[i],ls='None',label=r'$\~ \Sigma$'+' = ' + str(round(Xi/(2.*np.pi*Bjerrum**2)/(dielectric*temperature/(valency*lam_D)),3)))
##	      ax1.errorbar(np.array(z_den_theory)*(L_z/characteristic_length),P_T,yerr=None,ls='-',color=colors[i])#,label=r'$\Sigma$'+' = ' + str(round(Sigma_s,3)))

#		##All this code is legit, but the MATLAB files need a new naming convention... easy enough... - 02/21/12 11:55:56 	      
#	      if Bikerman=='yes':
#		   Bik_file='_Bik_zeta_' + ze + '.txt'
#		   x_Bik=[[float(x) for x in line.split()] for line in file('x' + Bik_file,"r").readlines()]
#		   c_Bik=[[float(x) for x in line.split()] for line in file('counter' + Bik_file,"r").readlines()]
#		   ax1.plot(np.array(x_Bik[0])*(lam_D/characteristic_length),c_Bik[0],color=colors[i],lw=1.5,ls='--') 
#	      if CarnahanStarling=='yes':
#		   CS_file='_CS_zeta_' + ze + '.txt'
#		   x_CS=[[float(x) for x in line.split()] for line in file('x' + CS_file,"r").readlines()]
#		   c_CS=[[float(x) for x in line.split()] for line in file('counter' + CS_file,"r").readlines()]
#		   ax1.plot(np.array(x_CS[0])*(lam_D/characteristic_length),c_CS[0],color=colors[i],lw=1.5,ls='-.') 
#	      if MPBoltz=='yes':
#		   MPBoltz_file='_MPBoltz_zeta_' + ze + '.txt'
#                   x_MPBoltz=[[float(x) for x in line.split()] for line in file('x' + MPBoltz_file,"r").readlines()]      		
#		   volt_MPBoltz=[[float(x) for x in line.split()] for line in file('volt' + MPBoltz_file,"r").readlines()]
#		   volt_MPBoltz=np.array(volt_MPBoltz[0])
#		   volt_MPBoltz_bulk=np.mean(volt_MPBoltz[len(volt_MPBoltz)/2-2:len(volt_MPBoltz)/2+2])
#		   if nonDim=='yes':
#		   	ax1.plot(np.array(x_MPBoltz[0])*(lam_D/characteristic_length),np.exp(volt_MPBoltz-volt_MPBoltz_bulk),color=colors[i],lw=2.0,ls=':')
#		   elif nonDim=='no':	   	
#		   	ax1.plot(np.array(x_MPBoltz[0])*(lam_D/characteristic_length),n0*np.exp(volt_MPBoltz-volt_MPBoltz_bulk),color=colors[i],lw=2.0,ls=':')
#    locs3=[] # A stupid name, but actually the words that print out on the right y-axis
#    if phi=='yes':	#This most likely will NOT work... - 02/21/12 13:56:36 
#	      locs1,labels=plt.yticks()
#	      for y in [x*(4.*np.pi/3.)*(sig*0.5)**3 for x in locs1]:
#		      s='%1.1E' % y
#		      locs3.append(s)  
#	      locs3.reverse()
#	      ax2=ax1.twinx()  	
#	      locs,labels=plt.yticks()  #These are the location of the new axis
#	      ax2.grid(True)
#	      plt.ylim(ymin=0)
#	      locs = [(x-min(locs1))/(max(locs1)-min(locs1)) for x in locs1] #What a pain in my ass to figure this shit out...
#	      locs.reverse()
#	      plt.yticks(locs, locs3,rotation=15)
#	      ax2.set_ylabel(r'$\Phi_{\sigma}$',size='x-large') 
#    if xlabel==r'$z/L_{z}$':
#      #These must be set manually! 
#      test=0
##      ax1.set_xlim(0,10)       
#      ax1.set_ylim(0,5) 
#    elif xlabel==r'$z/ \sigma_{WCA}$':
#      #These must be set manually! 
#      test=0

##    ax1.set_xlim(0.3,0.5)       
##    ax1.set_ylim(-0.05,0.05)      
#    ax1.set_xlabel(xlabel,size='x-large')
#    ax1.set_ylabel(ylabel,size='x-large')
#    ax1.legend(loc=0) 
#    plt.savefig(graphname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
##    plt.show()
#    plt.close()
#	##Finished plotting N_minus

    ##This plots N_free(z) for simulation and GC theory
    nonDim='no'
    nonDim='yes'
    CarnahanStarling='yes'
    CarnahanStarling='no'
    if nonDim=='yes':
      print 'Plotting NonDim N_free(z)...'
    else:
      print 'Plotting Dim N_free(z)...'	
    if Bikerman=='yes':
		print '\t\t...with Bikerman theory...'
		#Note to user: Bikerman theory must already be evaluated using MATLAB code, which is in ~/sims/PBik_Solver.m
    if CarnahanStarling=='yes':
	      print '\t...with Carnahan-Starling theory...'
    if MPBoltz=='yes':
    	print '\t...with Modified Poisson w/Boltzmann distribution theory...'
    NfS=[]
    GC_NfS=[]
    GC_Nf2S=[]
    SigmaCorrectedS=[]
    zhalfS=[]
    SigEffS=[]
    V_corS=[]
#    Max_effS=[]
    GCs=[]
    BikS=[]
    S_CS=[]
    i=-1
    f3=-1
    f4=-1
    fig=plt.figure()
#    fig.subplots_adjust(right=0.85)
#    fig.subplots_adjust(bottom=0.19)
#    fig.subplots_adjust(left=0.175)

	#F2 Settings
#    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.18)##Larger adds whitespace
#    fig.subplots_adjust(left=0.11) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.96) ##Smaller adds whitespace to top

	#F3 Settings
    fig.subplots_adjust(right=0.96) ##Lower puts more whitespace on the right
    fig.subplots_adjust(bottom=0.20)##Larger adds whitespace
    fig.subplots_adjust(left=0.14) #Larger puts more whitespace on the left
    fig.subplots_adjust(top=0.94) ##Smaller adds whitespace to top

    labelS = [r'$\~ \Sigma_{app}=0$',r'$0.2$',r'$0.5$',r'$0.7$',r'$1.0$',r'$1.3$',r'$1.5$',r'$1.8$',r'$2.2$',r'$2.5$',r'$2.9$',r'$3.3$',r'$3.8$',r'$4.3$',r'$4.8$',r'$5.9$',r'$9.0$']


#    labelS = [r'$12.4$',r'$11.9$',r'$6.5$',r'$5.0$',r'$3.9$',r'$3.4$',r'$3.0$',r'$2.2$',r'$1.9$',r'$1.6$',r'$1.0$',r'$0.5$',r'$0.0$']
    F3_markers = ['^','^','o','o','v','v','s','s','>','>','*','*']*3

    labelS = 10*labelS

    ax1=fig.add_subplot(111)
    SigmaMeasureS=[]
    Nf_for_EffS=[]
    ContactS,VdropS=[],[]
    VEffS=[]
    NtS=[]
    SCS0s=[]
    Keck_Mu_PhiS=[]
#    ax1.errorbar([100,101],[100,101],yerr=None,ls='-',color='k',label=r'$-\~ \rho_{f}^{GC}(\~\Sigma_{app},\lambda_{D}\approx10,\lambda_{B}=3,\sigma=1)$')
#    ax1.errorbar([100,101],[100,101],yerr=None,ls='-',color='k',label=r'$-\~ \rho_{f}^{GC} \left ( \tilde \Sigma_{app},\frac{\sigma}{\lambda_{D}},\frac{\lambda_{B}}{\lambda_{D}} \right )$')
#    ax1.errorbar([100,101],[100,101],yerr=None,ls='-',color='k',label=r'$-\~ \rho_{f}^{GC} \left ( \tilde \Sigma_{app},\sigma/\lambda_{D}\approx 0.1,\lambda_{B}/\lambda_{D}\approx 0.7 \right )$')
    ax1.errorbar([100,101],[100,101],yerr=None,ls='-',color='k',label=r'$-\~ \rho^{GC} $') #\left (\tilde z ,  \tilde \Sigma_{app} \right ) 
#    ax1.errorbar([100,101],[100,101],yerr=None,ls='None',color='white',label=r'$\~ \Sigma_{app}$') #\left (\tilde z ,  \tilde \Sigma_{app} \right ) 
    one_time_iterator=0
    z_forS=[]
    for (Nm,Np,rf_T,n0,Sigma_s,area,lam_D,L_bin,z_positions,L_z,characteristic_length,Bjerrum,filename,sigWCA,Xi,V,muex) in zip(N_plusS,N_minusS,rhof_Theory,n0s,SIGMAS,areaS,LDs,L_binS,z_S,L_zS,characteristic_lengthS,BjerrumS,filenameS,sigWCA_S,XiS,Volts,muexEV_S):
		i+=1
		dielectric=(4*np.pi*Bjerrum)**-1
		if i==len(markers): #This resets the markers
			i=0
			
		if nonDim=='yes':
#			Nf=[(npl-nm)/(area*L_bin*n0) for (npl,nm) in zip(Np,Nm)]
			Nf=[(npl-nm)/(area*L_bin*dielectric*temperature/(valency*lam_D**2)) for (npl,nm) in zip(Np,Nm)]  #This is consistent with other derivatoins
			graphname = 'N_free_NDim'
			ylabel=r'$-\rho \times 4 \pi \lambda_{B} \lambda_{D}^{2}/q_{\pm}$'
#			ylabel=r'$-\tilde \rho_{f}$'
		else: #This is dimensional
			Nf=[(npl-nm)/(area*L_bin) for (npl,nm) in zip(Np,Nm)]
			rf_T=[n0*t for t in rf_T]
			graphname = 'N_free_Dim'
			ylabel=r'$N_{free}$'

		
		##P1_PRL_F2 coding - 08/07/12 21:57:49 
#		one_time_iterator+=1
#		if one_time_iterator==1:
#			special_label = r'$\left (8.1,\frac{5}{11.08},\frac{1}{11.08} \right )$'
##			print special_label,'sig = ',sigWCA,'Bjer=',Bjerrum,'LD= ',lam_D
#		elif one_time_iterator==2:
#			special_label = r'$\left (6.0,\frac{5}{10.73},\frac{1}{10.73} \right )$'
##			print special_label,'sig = ',sigWCA,'Bjer=',Bjerrum,'LD= ',lam_D
#		elif one_time_iterator==3:
##			special_label = r'$(2.3,10,1,1)$'
#			special_label = r'$\left (2.3,\frac{1}{9.92},\frac{1}{9.92} \right )$'
##			print special_label,'sig = ',sigWCA,'Bjer=',Bjerrum,'LD= ',lam_D
#		elif one_time_iterator==4:
##			special_label = r'$(2.0,10,1,1)$'
#			special_label = r'$\left (2.0,\frac{1}{9.97},\frac{1}{9.97} \right )$'
##			print special_label,'sig = ',sigWCA,'Bjer=',Bjerrum,'LD= ',lam_D		
#		elif one_time_iterator==5:
##			special_label = r'$(0.12,10,20,1)$'
#			special_label = r'$\left (0.12,\frac{1}{9.91},\frac{20}{9.91} \right )$'
##			print special_label,'sig = ',sigWCA,'Bjer=',Bjerrum,'LD= ',lam_D
#		elif one_time_iterator==6:
##			special_label = r'$(0.05,10,20,1)$'
#			special_label = r'$\left (0.05,\frac{1}{10.01},\frac{20}{10.01} \right )$'
##			print special_label,'sig = ',sigWCA,'Bjer=',Bjerrum,'LD= ',lam_D
#			colors[i]='purple'

#			Nf=[(npl-nm)/(area*L_bin*dielectric*temperature/(valency*lam_D**2)) for (npl,nm) in zip(Np,Nm)] 


     		markers[i]='s'
	    	if filename in GC_sig5:
			label_str = r'$(\sigma,\lambda_{B},\lambda_{D}) = (5,1,\sim 11)$'
			i=11
	#		name='P1_PRL_F4a.png'
			if filename in ['Analyzed_GC_808_0.280418560259_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.364487579669_1.0_5.0_500.0_44.75.txt']:
				special_line='--'
	    	if filename in GC_sig6:
	  		label_str = r'$\sigma/\lambda_{D}\approx 0.6$'
	  		i=3
	    	if filename in GC_sig7:
	    		label_str = r'$\sigma/\lambda_{D}\approx 0.7$'
	    		i=6
	    	if filename in GC_sig8:
	    		label_str = r'$\sigma/\lambda_{D}\approx 0.8$'
	    		i=9
	    	if filename in GC_sig1:
	    		label_str = r'$\sigma/\lambda_{D}\approx 0.8$'
	    		special_marker = 's'
	    		col = 'orange'
	    	if filename in GC_LBhalf:
	    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.5/10$'
	    		label_str = r'$0.5/10$'
	    		i=0
	    		col = 'orange'
	    		special_marker = 'v'
	    	if filename in GC_LB3:
	    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 3/10$'
	    		label_str = r'$3/10$'
	    		i=3
	    		markers[i]='o'
	    		special_marker = 'o'
	    		col = 'limegreen'
	    	if filename in GC_LB5:
	    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 5/10$'
	    		label_str = r'$5/10$'
	    		i=5
	    		special_marker = 'd'
	    		col = 'cyan'
	    	if filename in GC_LB7:
	    	   	label_str = r'$\lambda_{B}/\lambda_{D}\approx 7/10$'
	    	   	label_str = r'$7/10$'
	    	   	i=7
	    	   	markers[i]='^'
	    	   	special_marker = '^'
	    	   	col = 'blue'
	    	if filename in GC_LB10:
	    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 10/10$'
	    		label_str = r'$10/10$'
	    		i=9
	    		special_marker ='D'
	    		col = 'red'
	    	if filename in IDEAL:
	    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 7/30$'
	    		label_str = r'$7/30$'
	    		i=11
	    	if filename in GC_LB_LD_5_20:
	    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 5/20$'
	    		label_str = r'$5/20$'
	    		i=2
	    		special_marker = '>'
	    		col = 'limegreen'
		if filename in LB_LD_Ratio_2:
			label_str = r'$\lambda_{B}/\lambda_{D}\approx 20/10$'	
			i=17
			markers[i]='*'
			special_ms=5


		#This could surely be consolidated within the lines of code below...
	        z_density = [(x+y)/(L_z*2.) for (x,y) in zip(z_positions[0:len(z_positions)-1],z_positions[1:len(z_positions)])]
	        z_den_theory=ztheory#[z/L_z for z in z_positions]  
		z_density = np.array(z_density)*L_z-0.5*sigWCA #This accounts for all of the wall
		z_density=z_density/L_z

#		special_marker='s'
		special_ls='None'
    		if filename in GC_sig5:
    			special_marker='p'
			if filename in ['Analyzed_GC_808_0.280418560259_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.364487579669_1.0_5.0_500.0_44.75.txt']:
				special_ls='--'
#				print "SIG5 Color = ",colors[i]
    			
    		if filename in LB_LD_Ratio_2:
    			special_marker='*'
			##Restore code with commented out lines below - 08/31/11 15:39:38 
####		col = ROYGBIV_map(np.min(SIGMAS),np.max(SIGMAS))


		#Functional F2 code - 08/08/12 01:04:58 
#		ax1.errorbar(np.array(z_den_theory)*(L_z/characteristic_length),-np.array(rf_T),yerr=None,ls='-',color=colors[i])#label=r'$\~ \Sigma$'+' = ' + str(round(Sigma_s/(dielectric*temperature/(valency*lam_D)),3)))
#		ax1.errorbar(np.array(z_density[1:])*(L_z/characteristic_length),Nf[1:],yerr=None,marker=special_marker,ms=4.5,color=colors[i],ls=special_ls)#,label=special_label)
		ax1.errorbar(np.array(z_density[1:])*(L_z/characteristic_length),np.array(Nf[1:])*2,yerr=None,marker=special_marker,ms=4.5,color=col,ls=special_ls)#,label=special_label)

####		Code dedicated to F3 - 08/08/12 01:25:43 
		f3+=1
		marker_special = '^'
		color_scale = colors[i]
		characteristic_length = lam_D
#		marker_special='h'
#		marker_special = F3_markers[i]
#		color_scale= colors[f3]
#		characteristic_length = sigWCA
#		if filename in ['Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt']:
#			ax1.errorbar(np.array(z_density[1:])*(L_z/characteristic_length),Nf[1:],yerr=None,marker=marker_special,mew=0.1,ms=7.0,color=color_scale,ls=':',label=labelS[i])
#		else:
#			ax1.errorbar(np.array(z_density[1:])*(L_z/characteristic_length),Nf[1:],yerr=None,marker=marker_special,mew=0.1,ms=4.5,color=color_scale,ls='None',label=labelS[i])




		###HERE  HERE  HERE  HERE  HERE  HERE  HERE  HERE  
		##This is the functioning code -- 09/11/13 16:51:59 
#		ax1.errorbar(np.array(z_den_theory)*(L_z/characteristic_length),-np.array(rf_T),yerr=None,ls='-',color=color_scale)#,label=r'$\~ \Sigma$'+' = ' + str(round(Sigma_s/(dielectric*temperature/(valency*lam_D)),3)))
#		ax1.errorbar(np.array(z_density[1:])*(L_z/characteristic_length),Nf[1:],yerr=None,marker=marker_special,mew=0.4,ms=4.5,color=color_scale,ls='None',label=labelS[i])


		##This might be unneeded
#		##Dedicated to code real F3 legend
#		f4+=1
#		marker_special = '^'
##		marker_special='h'
#		color_scale= colors[f3]
##		print filename,labelS[i]
#		ax1.errorbar(np.array(z_den_theory)*(L_z/characteristic_length),-np.array(rf_T),yerr=None,ls='-',color=color_scale)#,label=r'$\~ \Sigma$'+' = ' + str(round(Sigma_s/(dielectric*temperature/(valency*lam_D)),3)))
##		if filename!= 'Analyzed_GC_1128_4.23637486151_7.0_1.0_500.0_118.4.txt':
#		ax1.errorbar(np.array(z_density[1:])*(L_z/characteristic_length),Nf[1:],yerr=None,marker=marker_special,mew=0.1,ms=7.0,color=color_scale,ls='None',label=labelS[i])
##		ax1.errorbar(np.array(z_density[1:])*(L_z/characteristic_length),Nf[1:],yerr=None,marker=F3_markers[i],mew=0.1,ms=4.0,color='k',ls='None',label=labelS[i])

	        Integrand = 0.
	        correction=0.
		Nf=np.array([(npl-nm)/(area*L_bin) for (npl,nm) in zip(Np,Nm)])
		Nf=Nf[:len(Nf)/2]
	        z_restore = np.array([(x+y)/(L_z*2.) for (x,y) in zip(z_positions[0:len(z_positions)-1],z_positions[1:len(z_positions)])]) #SO, this was the original z_density required as input to the code below
		z_pos = z_restore[:len(z_restore)/2]
	        for (y1,y2,z) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos):
	        	if z*L_z<z_wall:
	        	       	Integrand=Integrand + 0.5*(y1+y2)*L_bin
	        	       	correction=correction+Integrand
			else:
			    	Integrand=Integrand + 0.5*(y1+y2)*L_bin	  
        	Sigma_meas = Integrand
		eff=(-Sigma_meas+correction)/(dielectric*temperature/(valency*lam_D))
		SigmaCorrectedS.append(abs(eff))


	        Integrand = eff*(dielectric*temperature/(valency*lam_D))
	        SigEff=[]
	        VEff=[]
	    	volt_correction = Sigma_s/(dielectric*temperature/(valency*lam_D))  - abs(eff)
	    	print 'volt correction equals ',volt_correction
	        Nf_for_Eff=[]
#	        for (y1,y2,z) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos):  #This is the original line of code -- 07/30/13 16:51:36 
	        for (y1,y2,z,volt) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos,V[:len(Nf)/2]):  #This is the original line of code -- 07/30/13 16:51:36 
	        	if z*L_z>z_wall:  #ONe more thing, try this without the equalds
	        		SigEff.append(Integrand) #There's something wrong with the first value assigned to this function
	        		test = volt-volt_correction*z
	        		VEff.append(test)
	        		Integrand=Integrand + 0.5*(y1+y2)*L_bin	
	        		if filename in GC_LB_LD_5_20:
	       				Nf_for_Eff.append(0.5*(y1+y2))	        		    		
	        		else:
	       				Nf_for_Eff.append(y1)
        	SigEffS.append(np.array(SigEff)/(dielectric*temperature/(valency*lam_D)))
        	VEffS.append(VEff)
        	Nf_for_EffS.append(np.array(Nf_for_Eff)/(dielectric*temperature/(valency*lam_D**2)))

        	NfS.append(Nf/(dielectric*temperature/(valency*lam_D**2)))
        	
		GC=[]
		CS=[]
		nf=[]
		nf2=[]
		Bik=[]
		zs=[]
		arr = (sigWCA**3/(24*Bjerrum*lam_D**2)) / 0.14457
		z_shift = (np.array(z_pos)*L_z-0.5*sigWCA-z_wall)/lam_D
		concon=0
		CorSig = abs(eff)
        	for (z,y1,y2,voltage) in zip(z_shift,Nf[0:len(Nf)-1]/(dielectric*temperature/(valency*lam_D**2)),Nf[1:len(Nf)]/(dielectric*temperature/(valency*lam_D**2)),V):		

#        		if filename in GC_LB_LD_5_20:  #Original peice of code
        		if filename in GC_LB_LD_5_20 or filename in F3_Data:  #Last minute  change made on 12/22/12 09:57:45 
				if z>0.: #Next try this with greater than, depending
					if concon==0:
						ContactS.append(y1)
						VdropS.append(voltage)
						concon+=1
					GC.append(np.exp(-z)*(np.sqrt(1+(CorSig/2.)**-2)-(CorSig/2.)**-1))
					zs.append(z)
					es = np.sqrt(4.+CorSig**2)
					Bik.append(np.exp(-z)*np.exp(((es**2-es-2)/(2*es))*arr + ((es**6-14*es**3-6*es**2+72)/(72*es**3))*arr**2)*(np.sqrt(1+(eff/2.)**-2)-(eff/2.)**-1)) #THe original, which is correct...        			
					nf.append(y1) 
					nf2.append(y2)
					
					if CarnahanStarling=='yes':    	 ##NEWW CODE
						Fitted = 'yes'
						#@Fitted = 'no'
						if Fitted == 'yes':
							CS_theory=np.array([[float(x) for x in line.split()] for line in file('MMA_fitWCA_'+filename[11:],"r").readlines()])    				
						else:
							CS_theory=np.array([[float(x) for x in line.split()] for line in file('MMA_WCA_'+filename[11:],"r").readlines()])

						SCS0 = CS_theory[1,4]
						if filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt':
#							print "NOT FORCING S=0"
							SCS0=0
#							print 'Forcibly setting SCS0 for SigApp ~= 0'
						CS.append(np.exp(-z)*SCS0) #THe original, which is correct...					
			else:
				if z>=0.:
					if concon==0:
						ContactS.append(y1)
						VdropS.append(voltage)
						concon+=1
					GC.append(np.exp(-z)*(np.sqrt(1+(CorSig/2.)**-2)-(CorSig/2.)**-1)) #THe original, which is correct...
					zs.append(z)
					es = np.sqrt(4.+CorSig**2)
					Bik.append(np.exp(-z)*np.exp(((es**2-es-2)/(2*es))*arr + ((es**6-14*es**3-6*es**2+72)/(72*es**3))*arr**2)*(np.sqrt(1+(eff/2.)**-2)-(eff/2.)**-1)) #THe original, which is correct...        			
					nf.append(y1)
					nf2.append(y2)
					if CarnahanStarling=='yes':    	
						Fitted = 'yes'
						#@Fitted = 'no'
						if Fitted == 'yes':
							CS_theory=np.array([[float(x) for x in line.split()] for line in file('MMA_fitWCA_'+filename[11:],"r").readlines()])    				
						else:
							CS_theory=np.array([[float(x) for x in line.split()] for line in file('MMA_WCA_'+filename[11:],"r").readlines()])

						SCS0 = CS_theory[1,4]
						if filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt':
#							print "NOT FORCING S=0"
							SCS0=0
#							print 'Forcibly setting SCS0 for SigApp ~= 0'
						CS.append(np.exp(-z)*SCS0) #THe original, which is correct...
		GCs.append(GC)
		z_forS.append(zs)
		BikS.append(Bik)
		GC_NfS.append(nf)
		GC_Nf2S.append(nf2)
		S_CS.append(CS)
    		if CarnahanStarling=='yes':
    			SCS0s.append(SCS0)

	    	volt_correction = Sigma_s/(dielectric*temperature/(valency*lam_D))  - CorSig
	    	keck_i = -1
		for (z,rp,rm,mu,v) in zip(z_positions,[npl/(area*L_bin*dielectric*temperature/(valency*lam_D**2)) for npl in Np] ,[nmi/(area*L_bin*dielectric*temperature/(valency*lam_D**2)) for nmi in Nm],muex,V):
			if z>z_wall and z<(L_z-z_wall):
				keck_i+=1
				v_cor = v - volt_correction*z
				if keck_i==0:
					zetaCor = v_cor
					V_corS.append(zetaCor)


#		if i==0:
#			outname='PMdata_GClike_a.txt'
#		elif i==1:
#			outname='PMdata_GClike_b.txt'
#		elif i==2:
#			outname='PMdata_GClike_c.txt'
#		elif i==3:
#			outname='PMdata_CSlike.txt'
#		elif i==4:
#			outname='PMdata_nonMF_a.txt'
#		elif i==5:
#			outname='PMdata_nonMF_b.txt'

#		if filenameS==GC_sig1:
#			outname = 'GCsig1_Keck_'+str(i)+'.txt'
#		else:
#			outname = 'CSlike_Keck_'+str(i)+'.txt'
##			outname = 'MasterTable_Keck_1_'+str(i)+'.txt'
#				
#	    	Keck_output=file(outname,"w")
#	    	volt_correction = Sigma_s/(dielectric*temperature/(valency*lam_D))  - CorSig
#	    	keck_i = -1
#		for (z,rp,rm,mu,v) in zip(z_positions,[npl/(area*L_bin*dielectric*temperature/(valency*lam_D**2)) for npl in Np] ,[nmi/(area*L_bin*dielectric*temperature/(valency*lam_D**2)) for nmi in Nm],muex,V):
#			if z>z_wall and z<(L_z-z_wall):
#				keck_i+=1
#				v_cor = v - volt_correction*z
#				if keck_i==0:
#					zetaCor = v_cor
#				Keck_output.write("%-1.5f\t%-1.5f\t%-1.5f\t%-1.5f\t%-1.5f\t%-1.5f\n" % (z-z_wall,mu,2.*n0*(rp+rm)*(np.pi/6.)*(0.954*sigWCA)**3,rp,rm,v))
#		Keck_output.write("%-1.5f\t%-1.5f\t%-1.5f\t%-1.5f\t%-1.5f\t%-1.5f" % (lam_D,CorSig,0.954*sigWCA,Bjerrum,n0,zetaCor))
#		Keck_output.close()

	
#		if i==0:
#			print '****************************\n****************************\n****************************\n\n\n'
##		print 'exportname = \"~/Desktop/MMA_WCA_'+filename[11:]+'\"'
##		print 'bulk = ',sigWCA**3/(24*Bjerrum*lam_D**2)
##		print 'effectivefield = ',-eff
##		print '\n\n'
#		print 'exportname = \"~/Desktop/MMA_fitWCA_'+filename[11:]+'\"'
#		print 'bulk = ',(sigWCA*0.954)**3/(24*Bjerrum*lam_D**2)
#		print 'effectivefield = ',-eff
#		print '\n\n'
#		print 'exportname = \"~/Desktop/MMA_kTWCA_'+filename[11:]+'\"'
#		print 'bulk = ',(sigWCA*0.969)**3/(24*Bjerrum*lam_D**2)
#		print 'effectivefield = ',-eff
#		print '\n\n'
		
    		if CarnahanStarling=='yes':    	
			Fitted = 'yes'
#			Fitted = 'no'
			if Fitted == 'yes':
				CS_theory=np.array([[float(x) for x in line.split()] for line in file('MMA_fitWCA_'+filename[11:],"r").readlines()])    				
			else:
#				CS_theory=np.array([[float(x) for x in line.split()] for line in file('MMA_WCA_'+filename[11:],"r").readlines()])
				CS_theory=np.array([[float(x) for x in line.split()] for line in file('MMA_kTWCA_'+filename[11:],"r").readlines()])
#			SCS0,CS = CS_theory[1,4],[]
#			if filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt':
#				SCS0=0
#				print 'Forcibly setting SCS0 for SigApp ~= 0'
#	        	for (z,n) in zip(z_shift,Nf/(dielectric*temperature/(valency*lam_D**2))):
#        			if z>=0.:
#					CS.append(np.exp(-z)*SCS0) #THe original, which is correct...
#			S_CS.append(CS)

			xax,yax = [],[]	
			for (z,rhofCS,effectiveCS) in zip(CS_theory[:,0],-np.array(CS_theory[:,2]),CS_theory[:,3]):
				if effectiveCS<=abs(eff):
					if len(xax)==0:
						zstart=z
					xax.append(z)
					yax.append(rhofCS)
			ax1.errorbar(np.array(xax)-zstart,yax,yerr=None,ls='-.',lw=1.5,color=colors[i],marker='None')
    ax1.set_xlim(0,2)
#    ax1.set_ylim(-0.5,25)
#    ax1.set_ylim(-0.5,15)
    ax1.set_ylim(-0.5,50)
#    xlabel=r'$\~ z = (z-\sigma/2)/\lambda_{D}$'
#    ylabel=r'$-\rho/2 q n^{\infty]$'
#    ax1.set_xlabel(xlabel,size='small')
#    ax1.set_ylabel(ylabel,size='small')
    ax1.xaxis.set_major_locator(MaxNLocator(4))
    plt.setp(ax1.get_xticklabels(), fontsize='small')
    plt.setp(ax1.get_yticklabels(), fontsize='small')
    fig.set_size_inches(3.37,2)
#    plt.show()
    plt.savefig('P1_PRL_F2_raw.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 
    print 'Free charge plot not save to a file!\n'
    fig.set_size_inches(10*.5,7)
##    ax1.legend(loc="best", ncol=2, shadow=False, title=r'$-\~ \rho_{f}^{GC} \left ( \tilde \Sigma_{app},\sigma/\lambda_{D}\approx 0.1,\lambda_{B}/\lambda_{D}\approx 0.7 \right )$')
#    ax1.legend(loc="best", ncol=2, shadow=False, title=r'$-\~ \rho_{f}^{GC} \left ( \tilde \Sigma_{app},\tilde z \right )$')
#           
#    ax1.set_xlim(-2.50,0)   
#    ax1.legend(loc="best", ncol=2, shadow=False)  
#    plt.savefig('P1_PRL_F2_legend2.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 

#    ax1.legend(loc="best", ncol=1, shadow=False, title=r'$-\~ \rho_{f}^{GC} \left ( \tilde \Sigma_{app},\tilde z \right )$')
#    ax1.legend(loc="best", ncol=1, shadow=False)
#    ax1.set_xlim(-10.,-2.)
#    plt.savefig('P1_PRL_F2_legend1.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 

#    ax1.set_xlim(-10.,-2.)
#    ax1.legend(loc="best", ncol=3, shadow=False)
#    plt.savefig('P1_PRL_F2_legend3.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 

#    fig.set_size_inches(3.37,3)
#    plt.savefig(graphname + '.png', dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 
#    plt.show()
    plt.close()
      ##Finished plottting N_free

#    print '***************\n***************\n***************\n***************\n***************\n***************\n'
#    for (filename,Sig) in zip(filenameS,SigmaCorrectedS): 
#    	print Sig
#    print '***************\n***************\n***************\n***************\n***************\n***************\n'

	###LaTex output
#    print r"\begin{tabular}{|c|c|c|c|}"
#    print '\hline'
#    print r"$\Sigma_{app}$ & $\lambda_{D}$ & $N_{+}+N{-}$ \zeta_{+}^{measured}\times k_{B}T/q_{\pm} \\ \hline"
##    old_ratio = sigWCA_S[0]/BjerrumS[0]
#    for (Sig_s,lam_D,sig,lam_B,N_p,N_m,V) in zip(SIGMAS,LDs,sigWCA_S,BjerrumS,N_plusS,N_minusS,V_corS):
##      if old_ratio !=sig/lam_B:
##      		print '\n'
##      		old_ratio=old_ratio
#      print r"%1.4f & %1.3f & %i & %1.2f \\ \hline" % (Sig_s,lam_D,int(round(np.sum(N_p)+np.sum(N_m))),V)



#####    This plots Keck(z) for simulation and GC theory
#    print 'Plotting Keck plot...'
##    VdropS=[]
##    ContactS=[]
#    ND_SigAppS=[]
#    EDL_CountS=[]
#    i=-1
#    f3=-1
#    NtS=[]
#    one_time_iterator=0
#    for (Nm,Np,rf_T,n0,Sigma_s,area,lam_D,L_bin,z_positions,L_z,characteristic_length,Bjerrum,filename,sigWCA,Xi,V) in zip(N_plusS,N_minusS,rhof_Theory,n0s,SIGMAS,areaS,LDs,L_binS,z_S,L_zS,characteristic_lengthS,BjerrumS,filenameS,sigWCA_S,XiS,Volts):
#		i+=1
#		dielectric=(4*np.pi*Bjerrum)**-1
#		if i==len(markers): #This resets the markers
#			i=0

#		special_marker='s'
#		special_ls='None'

#	        Integrand = 0.
#	        correction=0.
#		Nf=np.array([(npl-nm)/(area*L_bin) for (npl,nm) in zip(Np,Nm)])
#		Nf=Nf[:len(Nf)/2]
#		z_restore = np.array([(x+y)/(L_z*2.) for (x,y) in zip(z_positions[0:len(z_positions)-1],z_positions[1:len(z_positions)])]) #SO, this was the original z_density required as input to the code below
#		z_pos = z_restore[:len(z_restore)/2]
#		kk=0
#	        for (y1,y2,z,zeta) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos,V[0:len(V)-1]):
#	        	if z*L_z<=z_wall:
#	        	       	Integrand=Integrand + 0.5*(y1+y2)*L_bin
#	        	       	correction=correction+Integrand
#			else:
#				kk+=1
#				if kk==1:
##					VdropS.append(zeta)
#					print zeta
##					ContactS.append(0.5*(y1+y2)/(dielectric*temperature/(valency*lam_D**2)))
#			    	Integrand=Integrand + 0.5*(y1+y2)*L_bin	    
#        	Sigma_meas = Integrand	       	
#		eff=(-Sigma_meas+correction)/(dielectric*temperature/(valency*lam_D))
#		ND_SigAppS.append(-eff)
#		
#		Nt=np.array([(npl+nm)/(area*L_bin) for (npl,nm) in zip(Np,Nm)])
#		Nt = Nt[:len(Nt)/2]
#		Integrand=0.
#	        for (y1,y2,z) in zip(Nt[0:len(Nt)-1],Nt[1:len(Nt)],z_pos):
#			Integrand=Integrand + 0.5*(y1+y2)*L_bin	 
#		Integrand = Integrand - 2*n0
#		EDL_CountS.append(Integrand/(dielectric*temperature/(valency*lam_D)))
#		
#        	NtS.append(Nt/(0.5*dielectric*temperature/(valency*lam_D**2)))
#    fig=plt.figure()
#    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.19)##Larger adds whitespace
#    fig.subplots_adjust(left=0.14) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.96) ##Smaller adds whitespace to top
#    ax1=fig.add_subplot(111)

#    ax1.errorbar(VdropS,ND_SigAppS,yerr=None,ls='None',lw=1.5,color=colors[1],marker='s')
#    Vtheory = np.linspace(0,3.5,100)
#    ax1.errorbar(Vtheory,2*np.sinh(Vtheory/2.),yerr=None,ls='-',lw=1.5,color='k')

#    ax1.set_ylim(0,4)
#    ax1.set_xlim(0,3.5) 

##    for (x,y,z) in zip(VdropS,ND_SigAppS,EDL_CountS):
##	print x,'\t\t',y,'\t\t',z

#    xlabel=r'$\~ \zeta$'
#    ylabel=r'$-\Sigma_{app} \times 4 \pi \lambda_{B} \lambda_{D}^{2}/q$' 
#    ax1.set_xlabel(xlabel,size='x-small')
#    ax1.set_ylabel(ylabel,size='x-small')
#    plt.setp(ax1.get_xticklabels(), fontsize='x-small')
#    plt.setp(ax1.get_yticklabels(), fontsize='x-small')
#    fig.set_size_inches(3.37,2)
##    plt.show()
##    plt.savefig('Sig_vs_zeta.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 
#    plt.close()

#    quik= np.array(PhiBulkS)*(sigWCA*0.954/sig_HS)**3
#    print 'attempting volume frctions'

#    print L_bin
#    for (t,phib,filename) in zip(NtS,quik,filenameS):
#    		
#        length=0
#	print filename,'\n'
#	print np.array(t)*phib/2.
#	n=-0
#	print 'Phi(z),z'
#	while length<10:
#		n+=1
#		length+=L_bin
#		print t[n]*phib/2.,length
    
#    print 'Uncomment code above to get volume fraction informations'
#    fig=plt.figure()
#    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.19)##Larger adds whitespace
#    fig.subplots_adjust(left=0.20) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.96) ##Smaller adds whitespace to top
#    ax1=fig.add_subplot(111)

#    ax1.errorbar(VdropS,EDL_CountS,yerr=None,ls='None',lw=1.5,color=colors[1],marker='s')

##    ax1.set_ylim(-0.5,15)
#    ax1.set_xlim(0,3.5) 

#    xlabel=r'$\~ \zeta$'
#    ylabel=r'$  4 \pi \lambda_{B} \lambda_{D}^{2} \times \int \rho_{tot}\mathrm{d}z-2 \rho^{bulk}$' 
#    ax1.set_xlabel(xlabel,size='x-small')
#    ax1.set_ylabel(ylabel,size='x-small')
#    plt.setp(ax1.get_xticklabels(), fontsize='x-small')
#    plt.setp(ax1.get_yticklabels(), fontsize='x-small')
#    fig.set_size_inches(3.37,2)
##    plt.show()
#    plt.savefig('Count_vs_zeta.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 
#    plt.close()

#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    fig.subplots_adjust(right=0.98) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.18)##Larger adds whitespace
#    fig.subplots_adjust(left=0.28) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top

#	##Horizontal cconfiguration
##    fig.subplots_adjust(right=0.98) ##Lower puts more whitespace on the right
##    fig.subplots_adjust(bottom=0.27)##Larger adds whitespace
##    fig.subplots_adjust(left=0.14) #Larger puts more whitespace on the left
##    fig.subplots_adjust(top=0.95) ##Smaller adds whitespace to top

#    x_fit = np.linspace(0,13,100)
#    ax1.plot(x_fit,0.5*x_fit*np.sqrt(x_fit**2+4),color='k',lw=1.5,ls='-')  
#    for (x,y,c) in zip(SigmaCorrectedS,ContactS,colors):
#    	ax1.errorbar(x,y,yerr=None,marker='^',mew=0.4,ms=6.5,color=c)
#    ax1.set_ylim(0,50)
#    ax1.set_xlim(0,9.5) 
#    xlabel=r'$\Sigma_{app} \times 4 \pi \lambda_{B} \lambda_{D}/q_{\pm} $'
#    ylabel=r'$-\rho_{f}(0,\Sigma_{app}) \times 8 \pi \lambda_{B} \lambda_{D}^{2}/q_{\pm}$'
#    ylabel=r'$-\~\rho_{f}(0,\~ \Sigma_{app})$'
#    ax1.set_xlabel(xlabel,size='medium')
#    ax1.set_ylabel(ylabel,size='medium')
#    
#    plt.setp(ax1.get_xticklabels(), fontsize='small')
#    plt.setp(ax1.get_yticklabels(), fontsize='small')
#	##Horizontal
##    fig.set_size_inches(1.55*2,0.65*2)
#	##Vertical
#    fig.set_size_inches(0.92*2,1.25*2)

##    plt.show()
#    plt.savefig('Contact_vs_App.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 
#    plt.close()

#    ratioS,legend_key=np.array(sigWCA_S)/np.array(LDs),'sigs'
    ratioS,legend_key=np.array(BjerrumS)/np.array(LDs),'ele'

    print 'F3 Plotting rhofGC,rhofCS,rhofPM vs. S...'
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    i,cc=-1,-1  #Rando iteraters
    labels=[]
    F3_colorS=[]
    F3_markerS=[]    
    maxX = -1E9
    maxY = maxX
    print 'Investigate adding more files into this figure, especially for the GC figure!'

#    CS_PhiS = ['0.01','0.02','0.03','0.10','0.20','0.50']  
    CS_PhiS = ['0.01','0.02','0.07','0.10','0.50']  
    for P in CS_PhiS:
	Tname = 'MMA_CS_PhiB_'+P+'.txt'
	P=float(P)
    	CS_theory=np.array([[float(x) for x in line.split()] for line in file(Tname,"r").readlines()])
#    	print P,cols[dd]
	ax1.errorbar(4*np.array(CS_theory[:,1]),-np.array(CS_theory[:,2]),yerr=None,ls='-',lw=0.75/(1.-P),color='grey',marker='None')     	
#    ax1.errorbar(4*np.array(CS_theory[:,1]),-np.array(CS_theory[:,2]),yerr=None,ls='-',lw=0.75/(1.-P),color='blue',marker='None')     		
#	ax1.errorbar(CS_theory[:,1],-np.array(CS_theory[:,2]),yerr=None,ls='-',lw=0.75/(1.-P),color=cols[dd],marker='None')     	
    y=np.linspace(0,0.97)
    ax1.errorbar(4*y,4*y*(y**2+1)/(1-y**2)**2,yerr=None,ls='-',lw=0.55,color='k',marker='None')

    ax1.errorbar([100,101],[100,101],yerr=None,ls='None',color='white',label=r'$  $') #\left (\tilde z ,  \tilde \Sigma_{app} \right ) 
    xy_data=[[],[]]
    
    theory_col='blue'
    for (Nf,cs,gc,filename,ratio,bulkvc,hs,wca,zpos,Sig) in zip(GC_NfS,S_CS,GCs,filenameS,ratioS,PhiBulkS,sigHS_S,sigWCA_S,z_forS,SigEffS):  #GC_NfS is newest change
        cc+=1
        i+=1
    	special_line='None'
    	special_ms=3.5
    	
	lw_special = 0.75/(1.-bulkvc*(wca*0.95/hs)**3)
	edge_width=0.1
	if i>=len(markers):
		i=0
	special_ms=markers[i] 
	ls_x,lw_x='',0.
	ms_x = 3.0
	label_str=''
	legS=[]
	alabel=''
	if filename == 'Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		ls_x,lw_x=':',0.5
		special_ms='*'
		ms_x = 3.5
		label_str=r'$11.6$'		
	elif filename == 'Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		ls_x,lw_x=':',0.5
		special_ms='h'
		label_str=r'$11.1$'
	elif filename == 'Analyzed_GC_1312_0.0241789536991_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='h'
		label_str=r'$6.0$'
		alabel=r'$0.92$'
	elif filename == 'Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='o'
		label_str=r'$5.0$'
		alabel=r'$0.75$'
	elif filename == 'Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='o'	
		label_str=r'$3.9$'
		alabel=r'$0.59$'
	elif filename == 'Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='s'	
		label_str=r'$3.4$'
		alabel=r'$0.51$'
	elif filename == 'Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='s'
		ms_x=3.5
		label_str=r'$3.0$'
		alabel=r'$0.44$'
	elif filename == 'Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='p'	
		ms_x=3.5	
		label_str=r'$2.2$'
		alabel=r'$0.32$'
	elif filename == 'Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='p'	
		ms_x=3.5
		label_str=r'$1.9$'
		alabel=r'$0.26$'
	elif filename == 'Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='H'
		ms_x=3.5
		label_str=r'$1.6$'
		alabel=r'$0.22$'
	elif filename == 'Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='H'	
		label_str=r'$1.0$'
		alabel=r'$0.13$'
	elif filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='D'	
		label_str=r'$0.5$'
		alabel=r'$0.06$'
	elif filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='D'
		label_str=r'$0.0$'
		alabel=r'$0.00$'
#	label_str=alabel

#	print 'This must agree with the legend - CHECK THIS EVENTUALLY'
#	print filename,'\n',Sig[0],Sig[1]


    	F3_colorS.append(colors[cc])
    	F3_markerS.append([special_ms,edge_width,ms_x])
    		#This is the way the code way

    	label_str=label_str+r'$,$'+' '+alabel
	ax1.errorbar(4*np.array(cs),Nf,yerr=None,color=colors[cc],marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)	
	legS.append(label_str)


###	print zpos
###	cs_2 = [SCS_MMA(effective) for effective in Sig[1:]]
###	print len(cs_2),len(Sig[1:]),SCS_MMA(effective),len(Nf[1:])
###	ax1.errorbar(cs_2,Nf,yerr=None,color=colors[cc],marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)	
#####	
	if filename in ['Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt']:#,'Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt']:#or True:
#########		ax1.errorbar(cs,Nf,yerr=None,color=colors[cc],marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x)#,label=labelS[cc])
   	    	CS_theory=np.array([[float(x) for x in line.split()] for line in file('MMA_fitWCA_'+filename[11:],"r").readlines()])
   	 	ax1.errorbar(4*np.array(CS_theory[:,1]),-np.array(CS_theory[:,2]),yerr=None,ls='-',lw=lw_special,color=colors[cc],marker='None')	

#    	maxY = np.max([maxY,np.max(np.array(Nf))])
#    	maxX = np.max([maxX,np.max(4*np.array(cs))])

#    P = '0.0436891008773'  #AvERAGE
#    Tname = 'MMA_CS_PhiB_'+P+'.txt'
#    P=float(P)
#    CS_theory=np.array([[float(x) for x in line.split()] for line in file(Tname,"r").readlines()])
#    ax1.errorbar(CS_theory[:,1],-np.array(CS_theory[:,2]),yerr=None,ls='-',lw=0.75/(1.-P),color=theory_col,marker='None')  

    P = '0.0438'  #AvERAGE
    Tname = 'MMA_CS_PhiB_'+P+'.txt'
    P=float(P)
    CS_theory=np.array([[float(x) for x in line.split()] for line in file(Tname,"r").readlines()])
    ax1.errorbar(4*np.array(CS_theory[:,1]),-np.array(CS_theory[:,2]),yerr=None,ls='-',lw=0.75/(1.-P),color=theory_col,marker='None')  
#    ax1.errorbar(4*np.array(CS_theory[:,1]),-np.array(CS_theory[:,2]),yerr=None,ls='-',lw=0.75/(1.-P),color='grey',marker='None') 

    ax1.errorbar(4*y,4*y*(y**2+1)/(1-y**2)**2,yerr=None,ls='-',lw=0.55,color='k',marker='None')

    
    ax1.set_xlabel(r'$\~ S_{CS}$',fontsize=8.)
    ax1.set_ylabel(r'$-\rho /2 qe n^{\mathrm{B}}$',fontsize=8.)
    plt.setp(ax1.get_xticklabels(), fontsize='small')
    plt.setp(ax1.get_yticklabels(), fontsize='small')
    ax1.set_xlim(0,6) 
    ax1.xaxis.set_major_locator(MaxNLocator(6))
#    ax1.set_xlim(0,3.0) 
    ax1.set_ylim(0,6.5)
    plt.setp(ax1.get_xticklabels(), fontsize=8.)
    plt.setp(ax1.get_yticklabels(), fontsize=8.)
    fig.subplots_adjust(right=0.99) ##Lower puts more whitespace on the right
    fig.subplots_adjust(bottom=0.10)##Larger adds whitespace
    fig.subplots_adjust(left=0.10) #Larger puts more whitespace on the left
    fig.subplots_adjust(top=0.98) ##Smaller adds whitespace to top
    fig.set_size_inches(3.37,1.80)
#    print 'Not saveing rho vs CS'
    plt.savefig('F3_B_Free_vs_CS.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)  

#    ax1.set_xlim(-100,-50) 
##    ax1.set_xlim(0,3.0) 
#    ax1.set_ylim(-100,-50)     
#    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.05,borderpad=0.25,labelspacing=0.05,ncol=1,handletextpad=0.05)
#    plt.savefig('F3a_Legend_1col_new.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
###    plt.show()
##    fig.set_size_inches(6.,3.5)
##    plt.savefig('Pres_rhof_vs_meanfieldS_CS.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    plt.close()

#    print 'Plotting rhof(z) vs. SigmaEff(z) for F3...'
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    #For square inset
#    i,cc=-1,-1  #Rando iteraters
#    labels=[]
#    maxX = -1E9
#    maxY = maxX
#    xy_data=[[],[]]   
#    F3_SigmaCorrectedS=[] 

#    x_fit = np.linspace(0,17,100)		
#    ax1.plot(x_fit,x_fit*np.sqrt((x_fit/2.)**2+1),color='k',lw=0.30,ls='-',label=r'$ \rho_{\rm GC}(0,\Sigma)$') 
#    for (Nf,Sig,filename,ratio,bulkvc,hs,wca,lbin,lam_D) in zip(Nf_for_EffS,SigEffS,filenameS,ratioS,PhiBulkS,sigHS_S,sigWCA_S,L_binS,LDs):  #GC_NfS is newest change
#        cc+=1
#        i+=1
#    	special_line='None'
#    	special_ms=3.5

#	lw_special = 0.75/(1.-bulkvc*(wca*0.95/hs)**3)
#	edge_width=0.1
#	if i>=len(markers):
#		i=0
#	special_ms=markers[i]
#	ls_x,lw_x='',0.
#	ms_x = 3.0
#	label_str=''
#	alabel=''
#	if filename == 'Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		ls_x,lw_x=':',0.5
#		special_ms='*'
#		ms_x = 3.5
#		label_str=r'$11.6$'		
#	elif filename == 'Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		ls_x,lw_x=':',0.5
#		special_ms='h'
#		label_str=r'$11.1$'
#	elif filename == 'Analyzed_GC_1312_0.0241789536991_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='h'
#		label_str=r'$6.0$'
#		alabel=r'$0.92$'
#	elif filename == 'Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='o'
#		label_str=r'$5.0$'
#		alabel=r'$0.8$'
#	elif filename == 'Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='o'	
#		label_str=r'$3.9$'
#		alabel=r'$0.6$'
#	elif filename == 'Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='s'	
#		label_str=r'$3.4$'
#		alabel=r'$0.5$'
#	elif filename == 'Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='s'
#		ms_x=3.5
#		label_str=r'$3.0$'
#		alabel=r'$0.4$'
#	elif filename == 'Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='p'	
#		ms_x=3.5	
#		label_str=r'$2.2$'
#		alabel=r'$0.3$'
#	elif filename == 'Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='p'	
#		ms_x=3.5
#		label_str=r'$1.9$'
#		alabel=r'$0.5$'
#	elif filename == 'Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='H'
#		ms_x=3.5
#		label_str=r'$1.6$'
#		alabel=r'$0.22$'
#	elif filename == 'Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='H'	
#		label_str=r'$1.0$'
#		alabel=r'$0.22$'
#	elif filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='D'	
#		label_str=r'$0.5$'
#		alabel=r'$0.13$'
#	elif filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='D'
#		label_str=r'$0.0$'
#		alabel=r'$0.06$'
#    	maxX = np.max([maxX,np.max(-np.array(Sig))])
#    	maxY = np.max([maxY,np.max(np.array(Nf))])

#	FreeTrim = Nf[2:]  ## 2 ==int(round(wca/lbin))-1
#	SigEffTrim = Sig[2:]
#	for (x,y) in zip(SigEffTrim,FreeTrim): #This trims off the first layer of ions, to be explained in computer
#		if filename not in ['Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt']:
#			xy_data[0].append(abs(x))
#			xy_data[1].append(y)
#	F3_SigmaCorrectedS.append(-SigEffTrim[0])

#	print filename
#	print label_str,Sig[0]-SigEffTrim[0]


#	if filename not in ['Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt']:
#		ax1.errorbar(-np.array(SigEffTrim),FreeTrim,yerr=None,color=colors[cc],marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x)#,label=labelS[cc])

##	monolayer_sig = -np.array(SigEffTrim[-int(round(wca/lbin))-1:])
##	monolayer_rho = FreeTrim[-int(round(wca/lbin))-1:]
##	print 'Monolayer values?'
##	print monolayer_sig
##	print monolayer_rho
#	if filename in ['Analyzed_GC_1312_0.0241789536991_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt']:
#		ax1.errorbar(-Sig[2:4],Nf[2:4],yerr=None,color='white',marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x)#,label=labelS[cc])


#    P = '0.0438'  #AvERAGE
##    P = '0.0436891008773'
#    Tname = 'MMA_CS_PhiB_'+P+'.txt'
#    P=float(P)
#    CS_theory=np.array([[float(x) for x in line.split()] for line in file(Tname,"r").readlines()])
#    ax1.errorbar(CS_theory[:,3],-np.array(CS_theory[:,2]),yerr=None,ls='-',lw=0.75/(1.-P),color=theory_col,marker='None',label=r'$\rho_{\mathrm{CS}}(0, \Sigma)$')#

#    coeffs = np.polyfit(xy_data[0],xy_data[1],deg=2)
#    print coeffs
#    coeffs[2]=0. 
#    ax1.plot(x_fit,np.polyval(np.poly1d(coeffs),x_fit),color=theory_col,lw=1.0,ls=':',label=r'$\rho_{\mathrm{PM}}^{\prime} (\sigma, \Sigma^{\prime})$')
#    ymean = np.mean(xy_data[1])
#    R_sq = 1. - (sum([(yi-fi)**2 for (yi,fi) in zip(xy_data[1],np.polyval(np.poly1d(coeffs),xy_data[0]))]) / sum([(yi-ymean)**2 for yi in xy_data[1]]))
#    print '2nd order forced intercept R_sq = %1.5f' % R_sq

#    ax1.plot(x_fit,x_fit*np.sqrt((x_fit/2.)**2+1),color='k',lw=0.30,ls='-') 

##    ax1.set_xlabel(r'$ \int_{\infty}^{z} \rho \mathrm{d}\hat{z} / 2 \lambda_D q n^\infty $',fontsize=8.)
##    ax1.set_ylabel(r'$-\rho(z) /2 q n_{\pm}^{\mathrm{B}}$',fontsize=6.)
#    plt.setp(ax1.get_xticklabels(), fontsize=6.)
#    plt.setp(ax1.get_yticklabels(), fontsize=6.)
##    ax1.set_xlim(0,maxX*1.1) 
##    ax1.set_ylim(0,maxY*1.1)
#    ax1.set_xlim(0,7) 
#    ax1.set_ylim(0,6.5)
#    fig.set_size_inches(1.27,0.80) 
#    fig.subplots_adjust(right=0.98) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.15)##Larger adds whitespace
#    fig.subplots_adjust(left=0.08) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.99) ##Smaller adds whitespace to tops    
#    plt.savefig('F3_Ainset_rhof_vs_SigEff.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
##    plt.show()
#    plt.close()



###### rhof(z) vs. SigmaEff(z) for BikSucks...'
    print 'Plotting contact expressions'
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    #For square inset
    i,cc=-1,-1  #Rando iteraters
    labels=[]
    maxX = -1E9
    maxY = maxX
    xy_data=[[],[]]   
    F3_SigmaCorrectedS=[] 
    
    	
    P = '0.0438'  #AvERAGE
#    P = '0.0436891008773'
    Tname = 'MMA_CS_PhiB_'+P+'.txt'
    P=float(P)
    CS_theory=np.array([[float(x) for x in line.split()] for line in file(Tname,"r").readlines()])
    ax1.errorbar(CS_theory[:,3],-np.array(CS_theory[:,2]),yerr=None,ls='-',lw=0.75/(1.-P),color=theory_col,marker='None',label=r'$\mathrm{CS}$')#r'$\rho_{\mathrm{CS}}(0, \tilde \Sigma)$'

#    P = '0.0438'  #AvERAGE
##    P = '0.0436891008773'
    Tname1 = 'Bik_zeta_phi_3.9_0.438_rho.txt'
    Tname2 = 'Bik_zeta_phi_3.9_0.438_Sig.txt'
    P=float(P)
    Bik_theory_rho=np.array([[float(x) for x in line.split()] for line in file(Tname1,"r").readlines()])
    Bik_theory_Sig=np.array([[float(x) for x in line.split()] for line in file(Tname2,"r").readlines()])
#    Bik_theory = Bik_theory
#    print Bik_theory_Sig

    ax1.errorbar(Bik_theory_Sig[0],np.array(Bik_theory_rho[0])*0.5,yerr=None,ls='-.',lw=0.75/(1.-P),color='k',marker='None',label=r'$\mathrm{Bik}$')#
    x_fit = np.linspace(0,17,100)	
	#Gouy-Chapaman
    ax1.plot(x_fit,x_fit*np.sqrt((x_fit/2.)**2+1),color='k',lw=0.30,ls='-',label=r'$\mathrm{GC}$')
    col_etc=[]
#    colors = get_n_HexCol(len(filenameS))
    
    for (Nf,Sig,filename,ratio,bulkvc,hs,wca,lbin,lam_D) in zip(Nf_for_EffS,SigEffS,filenameS,ratioS,PhiBulkS,sigHS_S,sigWCA_S,L_binS,LDs):  #GC_NfS is newest change
        cc+=1
        i+=1
    	special_line='None'
    	special_ms=3.5

	lw_special = 0.75/(1.-bulkvc*(wca*0.95/hs)**3)
	edge_width=0.1
	if i>=len(markers):
		i=0
	special_ms=markers[i]
	ls_x,lw_x='',0.
	ms_x = 4.0
	label_str=''
	alabel=''
	if filename == 'Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		ls_x,lw_x=':',0.5
		special_ms='*'
		ms_x = 3.5
		label_str=r'$11.6$'		
	elif filename == 'Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		ls_x,lw_x=':',0.5
		special_ms='h'
		label_str=r'$11.1$'
	elif filename == 'Analyzed_GC_1312_0.0241789536991_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='h'
		label_str=r'$\tilde \Sigma = 6.0$'
		alabel=r'$0.92$'	
		## NEW
		colors[cc] = ROYGBIV_map(0.,10)	
	elif filename == 'Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='o'
		label_str=r'$5.0$'
		alabel=r'$0.8$'
		## NEW P2 Markers
		special_ms='*'	
		colors[cc] = ROYGBIV_map(1.,10)		
	elif filename == 'Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='o'	
		label_str=r'$3.9$'
		alabel=r'$0.6$'
		colors[cc] = ROYGBIV_map(2.,10)	
	elif filename == 'Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='s'	
		label_str=r'$3.4$'
		alabel=r'$0.5$'
		## NEW P2 Markers
		special_ms='d'	
		colors[cc] = ROYGBIV_map(3.,10)	
	elif filename == 'Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='s'
		ms_x=3.5
		label_str=r'$3.0$'
		alabel=r'$0.4$'
		colors[cc] = ROYGBIV_map(4.,10)	
	elif filename == 'Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='p'	
		ms_x=3.5	
		label_str=r'$2.2$'
		alabel=r'$0.3$'
		colors[cc] = ROYGBIV_map(5.,10)	
	elif filename == 'Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='p'	
		ms_x=3.5
		label_str=r'$1.9$'
		alabel=r'$0.5$'
		## NEW P2 Markers
		special_ms='v'	
		colors[cc] = ROYGBIV_map(6.,10)	
	elif filename == 'Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='H'
		ms_x=3.5
		label_str=r'$1.6$'
		alabel=r'$0.22$'
		colors[cc] = ROYGBIV_map(7.,10)	
	elif filename == 'Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='H'	
		label_str=r'$1.0$'
		alabel=r'$0.22$'
		## NEW P2 Markers
		special_ms='>'	
		colors[cc] = ROYGBIV_map(8.,10)	
	elif filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='D'	
		label_str=r'$0.5$'
		alabel=r'$0.13$'
		colors[cc] = ROYGBIV_map(9.,10)	
	elif filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='D'
		label_str=r'$0.0$'
		alabel=r'$0.06$'
    	maxX = np.max([maxX,np.max(-np.array(Sig))])
    	maxY = np.max([maxY,np.max(np.array(Nf))])
    	ms_x = 4.

	FreeTrim = Nf[2:]  ## 2 ==int(round(wca/lbin))-1
	SigEffTrim = Sig[2:]
	for (x,y) in zip(SigEffTrim,FreeTrim): #This trims off the first layer of ions, to be explained in computer
		if filename not in ['Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt']:
			xy_data[0].append(abs(x))
			xy_data[1].append(y)
	F3_SigmaCorrectedS.append(-SigEffTrim[0])

#	print filename
#	print label_str,Sig[0]-SigEffTrim[0]


#	ax1.errorbar(-np.array(Sig),Nf,yerr=None,color=colors[cc],marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x)#,label=labelS[cc])

	if filename not in ['Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt']:
		col_etc.append([colors[cc],special_ms,edge_width,ms_x,ls_x,lw_x])
		ax1.errorbar(-np.array(SigEffTrim),FreeTrim,yerr=None,color=colors[cc],marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)

#	if filename in ['Analyzed_GC_1312_0.0241789536991_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt']:
#		ax1.errorbar(-Sig[2:4],Nf[2:4],yerr=None,color='white',marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x)#,label=labelS[cc])


    P = '0.0438'  #AvERAGE
#    P = '0.0436891008773'
    Tname = 'MMA_CS_PhiB_'+P+'.txt'
    P=float(P)
    CS_theory=np.array([[float(x) for x in line.split()] for line in file(Tname,"r").readlines()])
    ax1.errorbar(CS_theory[:,3],-np.array(CS_theory[:,2]),yerr=None,ls='-',lw=0.75/(1.-P),color=theory_col,marker='None')#,label=r'$\rho_{\mathrm{CS}}(0, \tilde \Sigma)$')#

#    P = '0.0438'  #AvERAGE
##    P = '0.0436891008773'
    Tname1 = 'Bik_zeta_phi_3.9_0.438_rho.txt'
    Tname2 = 'Bik_zeta_phi_3.9_0.438_Sig.txt'
    P=float(P)
    Bik_theory_rho=np.array([[float(x) for x in line.split()] for line in file(Tname1,"r").readlines()])
    Bik_theory_Sig=np.array([[float(x) for x in line.split()] for line in file(Tname2,"r").readlines()])
#    Bik_theory = Bik_theory
#    print Bik_theory_Sig

    ax1.errorbar(Bik_theory_Sig[0],np.array(Bik_theory_rho[0])*0.5,yerr=None,ls='-.',lw=0.75/(1.-P),color='k',marker='None')#,label=r'$\rho_{\mathrm{Bik}}(0, \tilde \Sigma)$')#
	#Gouy-Chapaman
    ax1.plot(x_fit,x_fit*np.sqrt((x_fit/2.)**2+1),color='k',lw=0.30,ls='-')#,label=r'$\rho_{\mathrm{GC}}(0, \tilde \Sigma)$')#

#    coeffs = np.polyfit(xy_data[0],xy_data[1],deg=2)
#    print coeffs
#    coeffs[2]=0. 
#    ax1.plot(x_fit,np.polyval(np.poly1d(coeffs),x_fit),color=theory_col,lw=1.0,ls=':',label=r'$\rho_{\mathrm{PM}}^{\prime} (\sigma, \Sigma^{\prime})$')
#    ymean = np.mean(xy_data[1])
#    R_sq = 1. - (sum([(yi-fi)**2 for (yi,fi) in zip(xy_data[1],np.polyval(np.poly1d(coeffs),xy_data[0]))]) / sum([(yi-ymean)**2 for yi in xy_data[1]]))
#    print '2nd order forced intercept R_sq = %1.5f' % R_sq

    ax1.set_xlabel(r'$ \Sigma_{\mathrm{eff}}[\tilde z] / \Sigma_{\mathrm{ref}} $',fontsize=10.)
    ax1.set_ylabel(r'$-\rho[z/ \lambda_{\mathrm{D}} ] /2 qe n^{\mathrm{B}} $',fontsize=10.)

    plt.setp(ax1.get_xticklabels(), fontsize=8.)
    plt.setp(ax1.get_yticklabels(), fontsize=8.)
#    ax1.set_xlim(0,maxX*1.1) 
#    ax1.set_ylim(0,maxY*1.1)
    ax1.set_xlim(0,7) 
    ax1.set_ylim(0,6.5)

#Biksucks
    fig.set_size_inches(3.37,3.5)
    ax1.legend(loc='best',numpoints=1,prop={"size":10},columnspacing=0.07,borderpad=0.18,labelspacing=0.10,handletextpad=0.07)
#    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.05,borderpad=0.25,labelspacing=0.05,ncol=1,handletextpad=0.05)    
    plt.savefig('F2_Contact_BikSucks.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight') 
#    plt.show()
    plt.close()


###### Capacitance for BikSucks...'

#    print 'Plotting Charge-voltage expressions' 
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    #For square inset
#    i,cc=-1,-1  #Rando iteraters
#    labels=[]
#    maxX = -1E9
#    maxY = maxX
#    xy_data=[[],[]]   
#    F3_SigmaCorrectedS=[] 

#    P = '0.0438'  #AvERAGE
##    P = '0.0436891008773'
#    Tname = 'CS_0.0438_5.txt'
#    P=float(P)
#    QV_CS=np.array([[float(x) for x in line.split()] for line in file(Tname,"r").readlines()])
#    ax1.errorbar(QV_CS[:,4],QV_CS[:,5],yerr=None,ls='-',lw=0.75/(1.-P),color=theory_col,marker='None',label=r'${\rm CS}$')#

#    Tname = 'Bik_0.0438_5.txt'
#    P=float(P)
#    QV_CS=np.array([[float(x) for x in line.split()] for line in file(Tname,"r").readlines()])
#    ax1.errorbar(QV_CS[:,4],QV_CS[:,5],yerr=None,ls='--',lw=0.75/(1.-P),color=theory_col,marker='None',label=r'${\rm Bik}$')#
#    
#    theory = np.linspace(0,5)
#    ax1.errorbar(theory,2*np.sinh(theory/2),yerr=None,ls='-',lw=0.3,color='k',marker='None',label=r'${\rm GC}$')#
#    
##dimensional_Sigma = np.array([0.3848,0.2961,0.2259,0.1964,0.1700,0.1247,0.1051,0.0872,0.0553,0.0268 ,0])
#    for (Q,V,corr,plotables) in zip([6,5,3.9,3.4,3,2.2,1.9,1.6,1,0.5,0],[4.55,3.57,2.9,2.56,2.25,1.73,1.53,1.27,0.83,0.39,0],[0.92,0.75,0.59,0.51,0.44,0.32,0.26,0.22,0.13,0.06,0.0],col_etc):
#	ax1.errorbar(V,Q,yerr=None,color=plotables[0],marker=plotables[1],mew=plotables[2],ms=plotables[3],ls=plotables[4],lw=plotables[5])
#	ax1.errorbar(V,Q-corr,yerr=None,color='k',marker='*',mew=plotables[2],ms=plotables[3],ls=plotables[4],lw=plotables[5])

#    ax1.set_ylabel(r'$ -\Sigma / \Sigma_{\rm ref}$',fontsize=8.)
#    ax1.set_xlabel(r'$ \zeta qe / k_{\rm B}T $',fontsize=8.)

#    plt.setp(ax1.get_xticklabels(), fontsize=8.)
#    plt.setp(ax1.get_yticklabels(), fontsize=8.)
##    ax1.set_xlim(0,maxX*1.1) 
##    ax1.set_ylim(0,maxY*1.1)
#    ax1.set_xlim(0,5) 
#    ax1.set_ylim(0,7)
##Biksucks
#    fig.set_size_inches(3.37,3.5)
#    fig.subplots_adjust(right=0.98) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.10)##Larger adds whitespace
#    fig.subplots_adjust(left=0.13) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.98) ##Smaller adds whitespace to tops
#    ax1.legend(loc='best',numpoints=1,prop={"size":8},columnspacing=0.07,borderpad=0.18,labelspacing=0.10,handletextpad=0.07)

##    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.05,borderpad=0.25,labelspacing=0.05,ncol=1,handletextpad=0.05)    
#    plt.savefig('F5_ChargeVoltage_BikSucks_corr.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
##    plt.show()
#    plt.close()



    print 'Plotting Charge-voltage expressions'
    
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    #For square inset
    i,cc=-1,-1  #Rando iteraters
    labels=[]
    maxX = -1E9
    maxY = maxX
    xy_data=[[],[]]   
    F3_SigmaCorrectedS=[] 

    for (Sigma,V,plotables,filename) in zip(SigEffS,VEffS,col_etc,filenameS):
    	print filename,np.max(V)
	ax1.errorbar(V,-np.array(Sigma),yerr=None,color=plotables[0],marker=plotables[1],mew=plotables[2],ms=plotables[3],ls=plotables[4],lw=plotables[5])


    P = '0.0438'  #AvERAGE
#    P = '0.0436891008773'
    Tname = 'CS_0.0438_7.txt'
    P=float(P)
    QV_CS=np.array([[float(x) for x in line.split()] for line in file(Tname,"r").readlines()])
    ax1.errorbar(QV_CS[:,4],QV_CS[:,5],yerr=None,ls='-',lw=0.75/(1.-P),color=theory_col,marker='None',label=r'${\rm CS}$')#

    Tname = 'Bik_0.0438_5.txt'
    P=float(P)
    QV_CS=np.array([[float(x) for x in line.split()] for line in file(Tname,"r").readlines()])
    ax1.errorbar(QV_CS[:,4],QV_CS[:,5],yerr=None,ls='--',lw=0.75/(1.-P),color='k',marker='None',label=r'${\rm Bik}$')#
    
    theory = np.linspace(0,5)
    ax1.errorbar(theory,2*np.sinh(theory/2),yerr=None,ls='-',lw=0.3,color='k',marker='None',label=r'${\rm GC}$')#

    ax1.set_ylabel(r'$ -\Sigma[z / \lambda_{\rm D}] / \Sigma_{\rm ref}$',fontsize=10.)
    ax1.set_xlabel(r'$ (\phi[\tilde z] - \phi^{\rm B})  qe / k_{\rm B}T $',fontsize=10.)

    plt.setp(ax1.get_xticklabels(), fontsize=8.)
    plt.setp(ax1.get_yticklabels(), fontsize=8.)
#    ax1.set_xlim(0,maxX*1.1) 
#    ax1.set_ylim(0,maxY*1.1)
    ax1.set_xlim(0,7) 
    ax1.set_ylim(0,7)
#Biksucks
    fig.set_size_inches(3.37,3.5)
    fig.subplots_adjust(right=0.98) ##Lower puts more whitespace on the right
    fig.subplots_adjust(bottom=0.10)##Larger adds whitespace
    fig.subplots_adjust(left=0.13) #Larger puts more whitespace on the left
    fig.subplots_adjust(top=0.98) ##Smaller adds whitespace to tops
#    ax1.legend(loc='best',numpoints=1,prop={"size":10},columnspacing=0.07,borderpad=0.18,labelspacing=0.10,handletextpad=0.07)

#    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.05,borderpad=0.25,labelspacing=0.05,ncol=1,handletextpad=0.05)    
    plt.savefig('F5_ChargeVoltage_BikSucks.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight') 
#    plt.show()
    plt.close()    






#    print 'Plotting rhof vs. SPM for F3...'
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    SPM0S=[]
#    SCS0S=[] 
##    P = '0.0436891008773'  #AvERAGE
#    P = '0.0438'  #AvERAGE
#    Tname = 'MMA_CS_PhiB_'+P+'.txt'
#    SigCS,rhoCS,SCS_MMA=[],[],[]
#    CS_theory=np.array([[float(x) for x in line.split()] for line in file(Tname,"r").readlines()])
#    for (x,y,z) in zip(CS_theory[:,3],-np.array(CS_theory[:,2]),CS_theory[:,1]):
#	if x>=0 and y>=0:
#		SigCS.append(x)
#		rhoCS.append(y)
#		SCS_MMA.append(z)

#    print len(F3_SigmaCorrectedS),len(GC_NfS),len(z_forS),len(F3_colorS),len(F3_markerS),len(filenameS)
#    for (SigApp,Nf,zs,col,mark,filename) in zip(F3_SigmaCorrectedS,GC_NfS,z_forS,F3_colorS,F3_markerS,filenameS):
#	SPM0 = 4*S_PM0_or_CS(coeffs,SigApp,[SigCS,rhoCS,SCS_MMA],'PM') 
#	SPM0S.append(SPM0)
#	scs0 = 4*S_PM0_or_CS(coeffs,SigApp,[SigCS,rhoCS,SCS_MMA],'CS')
#	SCS0S.append(scs0)
#	leg=''
#	if filename  not in ['Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt']:    		
#	    	ax1.errorbar(np.exp(-np.array(zs))*SPM0,Nf,yerr=None,color=col,marker=mark[0],ms=mark[2],mew=mark[1],ls='None',label=leg)#marker=markers[i],color=colors[i]      	
#	else:
#		ax1.errorbar(-100*np.ones(len(Nf)),Nf,yerr=None,color=col,marker=mark[0],ms=mark[2],mew=mark[1],ls='None',label=leg)#marker=markers[i],color=colors[i]      			

#	if filename in ['Analyzed_GC_1312_0.0241789536991_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt']:
##		ax1.errorbar(-Sig[1:3],Nf[1:3],yerr=None,color='white',marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x)#,label=labelS[cc])

#	    	ax1.errorbar(np.exp(-np.array(zs[:2]))*SPM0,Nf[:2],yerr=None,color='white',marker=mark[0],ms=mark[2],mew=mark[1],ls='None',label=leg)#marker=markers[i],color=colors[i]      	


#    SigApps=np.linspace(0,np.max(SigmaCorrectedS)*1.2,1000)
#    Free_from_App = np.polyval(np.poly1d(coeffs),SigApps)
#    S_from_App = [S_PM0_or_CS(coeffs,SigApp,[SigCS,rhoCS,SCS_MMA],'PM') for SigApp in SigApps] 
#    ax1.errorbar(4*np.array(S_from_App),Free_from_App,yerr=None,ls=':',lw=1.1,color=theory_col,marker='None')#,label=r'$\rho_{\mathrm{PM}}(\bar{\Phi}\approx0.044)$')	##Old
#    ax1.set_ylabel(r'$-\rho /2 qe n^{\mathrm{B}}$',fontsize=8.)
#    plt.setp(ax1.get_xticklabels(), fontsize=8.)
#    plt.setp(ax1.get_yticklabels(), fontsize=8.)
#    fig.set_size_inches(3.37,1.75)
#    fig.subplots_adjust(right=0.99) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.11)##Larger adds whitespace
#    fig.subplots_adjust(left=0.10) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.99) ##Smaller adds whitespace to top
#    ax1.set_xlim(0,4*1.5) 
#    ax1.xaxis.set_major_locator(MaxNLocator(7))
##    ax1.set_ylim(-5,25)  
##    ax1.set_ylim(0,maxY*1.1)
#    ax1.set_ylim(0,6.5)
# #    ax1.legend(loc='upper left',numpoints=1,prop=dict(size=6.),columnspacing=0.05,borderpad=0.25,labelspacing=0.05,ncol=3,handletextpad=0.05)	    	
#    plt.savefig('F3_A_rhof_vs_SPM.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    plt.close()

##    fig.subplots_adjust(right=0.99) ##Lower puts more whitespace on the right
##    fig.subplots_adjust(bottom=0.11)##Larger adds whitespace
##    fig.subplots_adjust(left=0.12) #Larger puts more whitespace on the left
##    fig.subplots_adjust(top=0.965) ##Smaller adds whitespace to top

#    print 'Plotting S vs. S for F3...'
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    for (cs0,pm,col,mark,filename,test) in zip(SCS0S,SPM0S,F3_colorS,F3_markerS,filenameS,S_CS):
#	if filename  not in ['Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt']:
##	    	ax1.errorbar(cs0,pm/cs0,yerr=None,color=col,marker=mark[0],ms=mark[2],mew=mark[1],ls='None')
#		rat=pm/test[0]
#		if test[0]==0.0:
#			rat=1.0

#		if filename  not in ['Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt']:
#		    	ax1.errorbar(4*test[0],rat/4,yerr=None,color=col,marker=mark[0],ms=mark[2],mew=mark[1],ls='None')
##    ax1.set_xlabel(r'$S_{\mathrm{CS}}$',fontsize=8.)
##    ax1.set_ylabel(r'$S_{\mathrm{PM}}/S_{\mathrm{CS}}$',fontsize=6.)
#    plt.setp(ax1.get_xticklabels(), fontsize=6.)
#    plt.setp(ax1.get_yticklabels(), fontsize=6.)
#    ax1.xaxis.set_major_locator(MaxNLocator(4))
#    ax1.yaxis.set_major_locator(MaxNLocator(5))
#    ax1.set_xlim(0,5.95)#np.max(SGC0S)*1.05) 
##    ax1.set_ylim(0.99,1.05)  
#    fig.set_size_inches(0.90,0.75)
#    fig.subplots_adjust(right=0.99) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.168)##Larger adds whitespace
#    fig.subplots_adjust(left=0.255) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.955) ##Smaller adds whitespace to top
#    plt.savefig('F3_Binset_SPM0_vs_SCS0.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    plt.close()



#	##Figure 2 - homestretch on 12/19/12 00:55:27 
#    print '\nPlotting rhof(GC) vs. GC(z,Sigma)...'
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    i=-1
#    cc=-1
#    labels=[]
#    maxX = -1E9
#    xy_data=[[],[]]
#    y=np.linspace(0,0.95)
##    ax1.errorbar(y,4*y*(y**2+1)/(1-y**2)**2,yerr=None,ls='-',lw=1.0,color='k',marker='None',label=r'$- \rho_{\mathrm{GC}}$')	##Old
#    ax1.errorbar(4*y,4*y*(y**2+1)/(1-y**2)**2,yerr=None,ls='-',color='white',marker='None',label=r'$\lambda_{B}/\lambda_{D}$')	##Old
#    name=''
#    for (Nf,S,filename,ratio,Sig) in zip(GC_NfS,GCs,filenameS,ratioS,SigEffS):
##     	cc+=1
#    	label_str=''
#    	special_line='None'
#    	special_ms=3.5
#    	if filename in GC_sig1:
#    		if legend_key=='sigs':
#    			  label_str = r'$\sigma/\lambda_{D}\approx 0.1$'
# 		elif legend_key=='ele':
# 		   	label_str = r'$\lambda_{B}/\lambda_{D}\approx 1/10$'
# 		   	label_str = r'$1/10$'
## 		   	label_str = r'$\lambda_{B}=1$'
#     		i=1
#     		
#     		markers[i]='s'
#    	if filename in GC_sig5:
#		label_str = r'$(\sigma,\lambda_{B},\lambda_{D}) = (5,1,\sim 11)$'
#		i=11
##		name='P1_PRL_F4a.png'
#		if filename in ['Analyzed_GC_808_0.280418560259_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.364487579669_1.0_5.0_500.0_44.75.txt']:
#			special_line='--'
#    	if filename in GC_sig6:
#  		label_str = r'$\sigma/\lambda_{D}\approx 0.6$'
#  		i=3
#    	if filename in GC_sig7:
#    		label_str = r'$\sigma/\lambda_{D}\approx 0.7$'
#    		i=6
#    	if filename in GC_sig8:
#    		label_str = r'$\sigma/\lambda_{D}\approx 0.8$'
#    		i=9
#    	if filename in GC_LBhalf:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.5/10$'
#    		label_str = r'$0.5/10$'
#    		i=0
#    	if filename in GC_LB3:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 3/10$'
#    		label_str = r'$3/10$'
#    		i=3
#    		markers[i]='o'
#    	if filename in GC_LB5:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 5/10$'
#    		label_str = r'$5/10$'
#    		i=5
#    	if filename in GC_LB7:
#    	   	label_str = r'$\lambda_{B}/\lambda_{D}\approx 7/10$'
#    	   	label_str = r'$7/10$'
#    	   	i=7
#    	   	markers[i]='^'
#    	if filename in GC_LB10:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 10/10$'
#    		label_str = r'$10/10$'
#    		i=9
#    	if filename in IDEAL:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 7/30$'
#    		label_str = r'$7/30$'
#    		i=11
#    	if filename in GC_LB_LD_5_20:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 5/20$'
#    		label_str = r'$5/20$'
#    		i=2
#	if filename in LB_LD_Ratio_2:
#		label_str = r'$\lambda_{B}/\lambda_{D}\approx 20/10$'	
#		i=17
#		markers[i]='*'
#		special_ms=5

#	special_ms = 3.0
##	special_ms = 7.
#	if label_str in labels:
#		ax1.errorbar(4*np.array(S),Nf,yerr=None,color=ROYGBIV_map(ratio,np.max(ratioS)),marker=markers[i],ms=special_ms,mew=0.1,ls='None')#,label=label_str)#,label=r'$CS$')#marker=markers[i],color=colors[i]    
##		for (x,y,col) in zip(S,Nf,Sig):
##			special_color = ROYGBIV_map(abs(col),17,1)
##			ax1.errorbar(x,y,yerr=None,color=special_color,marker=markers[i],ms=special_ms,mew=0.1,ls='None')#,label=label_str)#,label=r'$CS$')#marker=markers[i],color=colors[i]    			
#		
#    	else:
#		ax1.errorbar(4*np.array(S),Nf,yerr=None,color=ROYGBIV_map(ratio,np.max(ratioS)),marker=markers[i],ms=special_ms,mew=0.1,ls='None',label=label_str)#,label=r'$CS$')#marker=markers[i],color=colors[i]    
##		for (x,y,col) in zip(S,Nf,Sig):
##			special_color = ROYGBIV_map(abs(col),17,1)
##			ax1.errorbar(x,y,yerr=None,color=special_color,marker=markers[i],ms=special_ms,mew=0.1,ls='None')#,label=label_str)#,label=r'$CS$')#marker=markers[i],color=colors[i]    			
##		ax1.errorbar(x,y,yerr=None,color=special_color,marker=markers[i],ms=special_ms,mew=0.1,ls='None',label=label_str)#,label=r'$CS$')#marker=markers[i],color=colors[i]    

#		labels.append(label_str)
##		cc=temp
#	for (x,y) in zip(S,Nf):
#			xy_data[0].append(x)
#			xy_data[1].append(y)
#    	maxX = np.max([maxX,np.max(S)])
#    y=np.linspace(0,0.95)
#    ax1.errorbar(4*y,4*y*(y**2+1)/(1-y**2)**2,yerr=None,ls='-',lw=1.0,color='k',marker='None')
#	#This ONLY works for F3 
##    ax1.errorbar(S,Nf,yerr=None,color='blue',marker=markers[i],ms=special_ms,ls=':')#,label=label_str)#,label=r'$CS$')#marker=markers[i],color=colors[i]    

##    ax1.set_xlabel(r'$S_{\mathrm{GC}}=\exp(-z) (\sqrt{4 + \Sigma} - 2 ) / \Sigma$',fontsize=8.)
#    ax1.set_ylabel(r'$-\rho / 2 q e n^{\mathrm{B}}$',fontsize=8.)
#    ax1.set_ylim(-5,50)  
#    ax1.xaxis.set_major_locator(MaxNLocator(4))
#    plt.setp(ax1.get_xticklabels(), fontsize=8.)
#    plt.setp(ax1.get_yticklabels(), fontsize=8.)
#    fig.subplots_adjust(right=0.97) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.12)##Larger adds whitespace
#    fig.subplots_adjust(left=0.12) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.965) ##Smaller adds whitespace to top
#    fig.set_size_inches(3.37,1.70)
#    ax1.legend(loc='lower right',numpoints=1,prop=dict(size=6.),columnspacing=0.05,borderpad=0.25,labelspacing=0.10,handletextpad=0.1)
#    plt.savefig('F2_B_rhof_vs_SGC.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 
##    plt.show()
#    plt.close()





	##Figure 2 - homestretch on 12/19/12 00:55:27 
#    print '\nPlotting **** BONUS ROUND **** rhof(GC) vs. GC(z,Sigma)...'
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    i=-1
#    cc=-1
#    labels=[]
#    maxX = -1E9
#    xy_data=[[],[]]
#    y=np.linspace(0,0.95)
##    ax1.errorbar(y,4*y*(y**2+1)/(1-y**2)**2,yerr=None,ls='-',lw=1.0,color='k',marker='None',label=r'$- \rho_{\mathrm{GC}}$')	##Old
##    ax1.errorbar(y,4*y*(y**2+1)/(1-y**2)**2,yerr=None,ls='-',color='white',marker='None',label=r'$\lambda_{B}/\lambda_{D}$')	##Old
#    name=''
#    bonus_colors=colors
#    xxx=-1
##    bonus_colors=[]
##    test=[]

#    for (Sig,lam_D,Bjerrum,filename) in zip(SigEffS,LDs,BjerrumS,filenameS):
##    		dielectric=(4*np.pi*Bjerrum)**-1
##    		print filename,round(Xi/(2.*np.pi*Bjerrum**2)/(dielectric*1.0/(1.0*lam_D)),1),lam_D
##    		if Bjerrum==20.0:
##    	    		test.append([round(Xi/(2.*np.pi*Bjerrum**2)/(dielectric*1.0/(1.0*lam_D)),1),4,1])
##    		else:
#    	    	test.append([Sig[0],17,1])
#    for t in test:
##    	  print t
#          bonus_colors.append(ROYGBIV_map(t[0],t[1],t[2]))

#    for (Nf,S,filename,ratio,Sig) in zip(GC_NfS,GCs,filenameS,ratioS,SigEffS):
#    	xxx+=1
##     	cc+=1
#    	label_str=''
#    	special_line='None'
#    	special_ms=3.5
#    	if filename in GC_sig1:
#    		if legend_key=='sigs':
#    			  label_str = r'$\sigma/\lambda_{D}\approx 0.1$'
# 		elif legend_key=='ele':
# 		   	label_str = r'$\lambda_{B}/\lambda_{D}\approx 1/10$'
# 		   	label_str = r'$1/10$'
## 		   	label_str = r'$\lambda_{B}=1$'
#     		i=1
#     		
#     		markers[i]='s'
#    	if filename in GC_sig5:
#		label_str = r'$(\sigma,\lambda_{B},\lambda_{D}) = (5,1,\sim 11)$'
#		i=11
##		name='P1_PRL_F4a.png'
#		if filename in ['Analyzed_GC_808_0.280418560259_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.364487579669_1.0_5.0_500.0_44.75.txt']:
#			special_line='--'
#    	if filename in GC_sig6:
#  		label_str = r'$\sigma/\lambda_{D}\approx 0.6$'
#  		i=3
#    	if filename in GC_sig7:
#    		label_str = r'$\sigma/\lambda_{D}\approx 0.7$'
#    		i=6
#    	if filename in GC_sig8:
#    		label_str = r'$\sigma/\lambda_{D}\approx 0.8$'
#    		i=9
#    	if filename in GC_LBhalf:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.5/10$'
#    		label_str = r'$0.5/10$'
#    		i=0
#    	if filename in GC_LB3:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 3/10$'
#    		label_str = r'$3/10$'
#    		i=3
#    		markers[i]='o'
#    	if filename in GC_LB5:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 5/10$'
#    		label_str = r'$5/10$'
#    		i=5
#    	if filename in GC_LB7:
#    	   	label_str = r'$\lambda_{B}/\lambda_{D}\approx 7/10$'
#    	   	label_str = r'$7/10$'
#    	   	i=7
#    	   	markers[i]='^'
#    	if filename in GC_LB10:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 10/10$'
#    		label_str = r'$10/10$'
#    		i=9
#    	if filename in IDEAL:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 7/30$'
#    		label_str = r'$7/30$'
#    		i=11
#    	if filename in GC_LB_LD_5_20:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 5/20$'
#    		label_str = r'$5/20$'
#    		i=2
#	if filename in LB_LD_Ratio_2:
#		label_str = r'$\lambda_{B}/\lambda_{D}\approx 20/10$'	
#		i=17
#		markers[i]='*'
#		special_ms=5
#	special_ms = 3.0
##	special_ms = 7.

#	print 'here in progress'

#	ax1.errorbar(4*np.array(S),Nf,yerr=None,color=bonus_colors[xxx],marker=markers[i],ms=special_ms,mew=0.1,ls='None')#,label=label_str)#,label=r'$CS$')#marker=markers[i],color=colors[i]    

##	s_col = np.linspace(0,np.max(S),300)
##	ax1.errorbar(s_col,np.ones(len(s_col)),yerr=None,color=bonus_colors[xxx],marker='s',ms=10.,mew=0.0,ls='None')
##	ax1.errorbar(s_col,np.ones(len(s_col))*1.1,yerr=None,color=bonus_colors[xxx],marker='s',ms=10.,mew=0.0,ls='None')
##	ax1.errorbar(s_col,np.ones(len(s_col))*1.2,yerr=None,color=bonus_colors[xxx],marker='s',ms=10.,mew=0.0,ls='None')
##	ax1.errorbar(s_col,np.ones(len(s_col))*0.9,yerr=None,color=bonus_colors[xxx],marker='s',ms=10.,mew=0.0,ls='None')



##	if label_str in labels:
##		ax1.errorbar(S,Nf,yerr=None,color=ROYGBIV_map(ratio,np.max(ratioS)),marker=markers[i],ms=special_ms,mew=0.1,ls='None')#,label=label_str)#,label=r'$CS$')#marker=markers[i],color=colors[i]    
###		for (x,y,col) in zip(S,Nf,Sig):
###			special_color = ROYGBIV_map(abs(col),17,1)
###			ax1.errorbar(x,y,yerr=None,color=special_color,marker=markers[i],ms=special_ms,mew=0.1,ls='None')#,label=label_str)#,label=r'$CS$')#marker=markers[i],color=colors[i]    			
##		
##    	else:
##		ax1.errorbar(S,Nf,yerr=None,color=ROYGBIV_map(ratio,np.max(ratioS)),marker=markers[i],ms=special_ms,mew=0.1,ls='None',label=label_str)#,label=r'$CS$')#marker=markers[i],color=colors[i]    
###		for (x,y,col) in zip(S,Nf,Sig):
###			special_color = ROYGBIV_map(abs(col),17,1)
###			ax1.errorbar(x,y,yerr=None,color=special_color,marker=markers[i],ms=special_ms,mew=0.1,ls='None')#,label=label_str)#,label=r'$CS$')#marker=markers[i],color=colors[i]    			
###		ax1.errorbar(x,y,yerr=None,color=special_color,marker=markers[i],ms=special_ms,mew=0.1,ls='None',label=label_str)#,label=r'$CS$')#marker=markers[i],color=colors[i]    

##		labels.append(label_str)
###		cc=temp
#	for (x,y) in zip(S,Nf):
#			xy_data[0].append(x)
#			xy_data[1].append(y)
#    	maxX = np.max([maxX,np.max(S)])
#    y=np.linspace(0,0.95)
#    ax1.errorbar(4*y,4*y*(y**2+1)/(1-y**2)**2,yerr=None,ls='-',lw=1.0,color='k',marker='None')
#    
#	#This ONLY works for F3 
##    ax1.errorbar(S,Nf,yerr=None,color='blue',marker=markers[i],ms=special_ms,ls=':')#,label=label_str)#,label=r'$CS$')#marker=markers[i],color=colors[i]    

##    ax1.set_xlabel(r'$S_{\mathrm{GC}}=\exp(-z) (\sqrt{4 + \Sigma} - 2 ) / \Sigma$',fontsize=8.)
#    ax1.set_ylabel(r'$-\rho / 2 qe  n^{\mathrm{B}}$',fontsize=8.)
#    ax1.set_ylim(-5,50)  
#    ax1.set_xlim(0,4)  
#    ax1.xaxis.set_major_locator(MaxNLocator(4))
#    plt.setp(ax1.get_xticklabels(), fontsize=8.)
#    plt.setp(ax1.get_yticklabels(), fontsize=8.)
#    fig.subplots_adjust(right=0.97) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.12)##Larger adds whitespace
#    fig.subplots_adjust(left=0.12) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.965) ##Smaller adds whitespace to top
#    fig.set_size_inches(3.37,1.70)
##    ax1.legend(loc='lower right',numpoints=1,prop=dict(size=6.),columnspacing=0.05,borderpad=0.25,labelspacing=0.10,handletextpad=0.1)
#    plt.savefig('F2_B_rhof_vs_SGC_bonus.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 
##    plt.savefig('F2_B_rhof_vs_SGC_bonus_scale.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 
##    plt.show()
#    plt.close()

#    print 'Plotting rhof(z) vs. SigmaEff(z) for F2...'
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    i,cc=-1,-1  #Rando iteraters
#    labels=[]
#    maxX = -1E9
#    xy_data=[[],[]]
#    xy_GC_LBhalf=[[],[]]
#    xy_GC_sig1=[[],[]]
#    xy_GC_LB_LD_5_20=[[],[]]
#    xy_GC_LB3=[[],[]]
#    xy_GC_LB5=[[],[]]
#    xy_GC_LB7=[[],[]]
#    xy_GC_LB10=[[],[]]
#    F2_colorS=[]
#    F2_markerS=[]
##    SigEff_output=file("SigEff_Output_LB7.txt","w")
#    for (Nf,Sig,filename,ratio,bulkvc,App) in zip(Nf_for_EffS,SigEffS,filenameS,ratioS,PhiBulkS,SigmaCorrectedS):  #GC_NfS is newest change
#        cc+=1
#    	label_str=''
#    	special_line='None'
#    	special_ms=3.5
#    	if filename in GC_sig1:
#    		if legend_key=='sigs':
#    			  label_str = r'$\sigma/\lambda_{D}\approx 0.1$'
# 		elif legend_key=='ele':
# 		   	label_str = r'$\lambda_{B}/\lambda_{D}\approx 1/10$'
## 		   	label_str = r'$\lambda_{B}=1$'
#     		i=1
#		for (x,y) in zip(Sig,Nf):
#				xy_GC_sig1[0].append(-x)
#				xy_GC_sig1[1].append(y)     		
#     		markers[i]='s'
#    	if filename in GC_sig5:
#		label_str = r'$(\sigma,\lambda_{B},\lambda_{D}) = (5,1,\sim 11)$'
#		i=11
#		if filename in ['Analyzed_GC_808_0.280418560259_1.0_5.0_500.0_44.75.txt','Analyzed_GC_808_0.364487579669_1.0_5.0_500.0_44.75.txt']:
#			special_line='--'
#    	if filename in GC_sig6:
#  		label_str = r'$\sigma/\lambda_{D}\approx 0.6$'
#  		i=3
#    	if filename in GC_sig7:
#    		label_str = r'$\sigma/\lambda_{D}\approx 0.7$'
#    		i=6
#    	if filename in GC_sig8:
#    		label_str = r'$\sigma/\lambda_{D}\approx 0.8$'
#    		i=9
#    	if filename in GC_LBhalf:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.5/10$'
#		for (x,y) in zip(Sig,Nf):
#				xy_GC_LBhalf[0].append(-x)
#				xy_GC_LBhalf[1].append(y)     
#    		i=0
#    	if filename in GC_LB3:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 3/10$'
#		for (x,y) in zip(Sig,Nf):
#				xy_GC_LB3[0].append(-x)
#				xy_GC_LB3[1].append(y)   
#    		i=3
#    		markers[i]='o'
#    	if filename in GC_LB5:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 5/10$'
#		for (x,y) in zip(Sig,Nf):
#				xy_GC_LB5[0].append(-x)
#				xy_GC_LB5[1].append(y)   
#    		i=5
#    	if filename in GC_LB7:
#    	   	label_str = r'$\lambda_{B}/\lambda_{D}\approx 7/10$'
#		for (x,y) in zip(Sig,Nf):
#				xy_GC_LB7[0].append(-x)
#				xy_GC_LB7[1].append(y)   
#    	   	i=7
#    	   	markers[i]='^'
#    	if filename in GC_LB10:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 10/10$'
#		for (x,y) in zip(Sig,Nf):
#				xy_GC_LB10[0].append(-x)
#				xy_GC_LB10[1].append(y)   
#    		i=9
#    		special_line=':'
#    	if filename in IDEAL:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 7/30$'
#    		i=11
#    	if filename in GC_LB_LD_5_20:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 5/20$'
#		for (x,y) in zip(Sig,Nf):
#				xy_GC_LB_LD_5_20[0].append(-x)
#				xy_GC_LB_LD_5_20[1].append(y)   
#    		i=2
#	if filename in LB_LD_Ratio_2:
#		label_str = r'$\lambda_{B}/\lambda_{D}\approx 20/10$'
##	   	label_str = r'$\lambda_{B}=20$'		
#		i=17
#		markers[i]='*'
#		special_ms=5.
#	if filename in sig_LD_3_over_15:
#		label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.1/15$'
##	   	label_str = r'$\lambda_{B}=20$'		
#		i=12
##		markers[i]='*'

#	special_ms = 3.0
#	F2_colorS.append(ROYGBIV_map(ratio,np.max(ratioS)))
#	F2_markerS.append(markers[i])

#	if label_str in labels:
#	    	ax1.errorbar(-np.array(Sig),np.array(Nf),yerr=None,color=ROYGBIV_map(ratio,np.max(ratioS)),marker=markers[i],ms=special_ms,mew=0.1,ls='None')#,label=r'$CS$')#marker=markers[i],color=colors[i]    
#    	else:
#		ax1.errorbar(-np.array(Sig),np.array(Nf),yerr=None,color=ROYGBIV_map(ratio,np.max(ratioS)),marker=markers[i],ms=special_ms,mew=0.1,ls='None')#,label=label_str)
#		labels.append(label_str)
##	if App>1:
#	for (x,y) in zip(Sig,Nf):
#			xy_data[0].append(-x)
#			xy_data[1].append(y)
##			SigEff_output.write("%1.5f\t\t%1.5f\n" % (-x,y))			
#    	maxX = np.max([maxX,np.max(-np.array(Sig))])
##    SigEff_output.close()  
#    x_fit = np.linspace(0,maxX*1.3,100)		
#    ax1.plot(x_fit,x_fit*np.sqrt((x_fit/2.)**2+1),color='k',lw=1.0,ls='-',label=r'$\rho_{\mathrm{GC}}(0,\Sigma)$')#,label=r'$\tilde \rho_{f}^{GC} = -\~\Sigma \sqrt{(\~\Sigma/2)^{2}+1}$')  

#    coeffs = np.polyfit(xy_data[0],xy_data[1],deg=2)
##    print coeffs
#    coeffs[2]=0. 
#    ax1.plot(x_fit,np.polyval(np.poly1d(coeffs),x_fit),color='grey',lw=1.1,ls=':',label=r'$\rho_{\mathrm{PM}}(0,\Sigma)$')#,label=leg_info)
#    ymean = np.mean(xy_data[1])
#    R_sq = 1. - (sum([(yi-fi)**2 for (yi,fi) in zip(xy_data[1],np.polyval(np.poly1d(coeffs),xy_data[0]))]) / sum([(yi-ymean)**2 for yi in xy_data[1]]))
#    print '2nd order forced intercept R_sq = %1.5f' % R_sq

##    GC_contact = lambda SIGMA: (SIGMA*np.sqrt((SIGMA/2.)**2+1))
##    R_sq = 1. - (sum([(yi-fi)**2 for (yi,fi) in zip(xy_data[1],[GC_contact(x) for x in xy_data[0]])]) / sum([(yi-ymean)**2 for yi in xy_data[1]]))
##    print 'Using GC contact R_sq = %1.5f' % R_sq
##    ax1.set_xlabel(r'$ \int_{\infty}^{z} \rho \mathrm{d}\hat{z} / 2 \lambda_D q n^\infty $',fontsize=8.)
#    ax1.set_ylabel(r'$-\tilde \rho( z / \lambda_{\rm D})$',fontsize=6.)
#    plt.setp(ax1.get_xticklabels(), fontsize=6.)
#    plt.setp(ax1.get_yticklabels(), fontsize=6.)
#    ax1.set_xlim(0,10.6) 
#    ax1.set_ylim(0,50)  
#    #This is for an inset
#    fig.set_size_inches(1.7,1.0) 
#    fig.subplots_adjust(right=0.99) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.14)##Larger adds whitespace
#    fig.subplots_adjust(left=0.19) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.95) ##Smaller adds whitespace to tops
#    ax1.legend(loc='best',numpoints=2,prop=dict(size=6.),columnspacing=0.07,borderpad=0.20,labelspacing=0.10,handletextpad=0.1)
#    plt.savefig('F2_Ainset_rhof_vs_SigEff.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
##    plt.show()
##    ax1.legend(loc=0)
##    fig.set_size_inches(11,7)
##    plt.savefig('rhof_vs_SigEff_legend.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 
#    plt.close()

#    print 'Plotting rhof vs. SPM for F2...'
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    SGC0S=[]
#    SPM0S=[]
#    for (SigApp,Nf,zs,col,mark,filename) in zip(SigmaCorrectedS,GC_NfS,z_forS,F2_colorS,F2_markerS,filenameS):
#    		SPM0 = 4*S_PM0_GC(coeffs,SigApp)
#    		SGC0S.append((np.sqrt(1+(SigApp/2.)**-2)-(SigApp/2.)**-1))
#		SPM0S.append(SPM0)
#	    	ax1.errorbar(np.exp(-np.array(zs))*SPM0,Nf,yerr=None,color=col,marker=mark,ms=3.0,mew=0.1,ls='None')#,label=r'$CS$')#marker=markers[i],color=colors[i]    

#    SigApps=np.linspace(0,np.max(SigmaCorrectedS)*3,1000)
#    Free_from_App = np.polyval(np.poly1d(coeffs),SigApps)
#    S_from_App = [S_PM0_GC(coeffs,SigApp) for SigApp in SigApps]
#    ax1.errorbar(4*np.array(S_from_App),Free_from_App,yerr=None,ls=':',lw=1.1,color='grey',marker='None',label=r'$\rho_{\mathrm{PM}}$')	##Old
##    ax1.errorbar(y,100*y,yerr=None,ls='-',lw=1.0,color='k',marker='None',label=r'$\rho_{\mathrm{GC}}$')
##    ax1.set_xlabel(r'$S_{\mathrm{PM}}$',fontsize=8.)
#    ax1.set_ylabel(r'$-\rho / 2 qe n^{\mathrm{B}}$',fontsize=8.)
#    ax1.xaxis.set_major_locator(MaxNLocator(4))
#    plt.setp(ax1.get_xticklabels(), fontsize=8.)
#    plt.setp(ax1.get_yticklabels(), fontsize=8.)

#    fig.set_size_inches(3.37,1.70)
#    fig.subplots_adjust(right=0.97) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.12)##Larger adds whitespace
#    fig.subplots_adjust(left=0.12) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.965) ##Smaller adds whitespace to top
#    print 'Size xlim based upon this!'
#    print np.max(SPM0S)
##    ax1.set_xlim(0,np.max(SPM0S)*1.01) 
#    ax1.set_ylim(-5,50)  
##    ax1.xaxis.set_major_locator(MaxNLocator(6))
##    ax1.legend(loc='lower right',numpoints=2,prop=dict(size=6.),columnspacing=0.10,borderpad=0.25,labelspacing=0.10,handletextpad=0.1)
#    plt.savefig('F2_A_rhof_vs_SPM.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    plt.close()

#    print 'Plotting S vs. S for F2...'
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    for (gc,pm,col,mark) in zip(SGC0S,SPM0S,F2_colorS,F2_markerS):
#	    	ax1.errorbar(4*gc,pm/(4*gc),yerr=None,color=col,marker=mark,ms=3.0,mew=0.1,ls='None')#,label=r'$CS$')#marker=markers[i],color=colors[i]    
##    ax1.errorbar(np.linspace(0,5),np.linspace(0,5),yerr=None,color='k',ls='-')#,label=r'$CS$')#marker=markers[i],color=colors[i]   
##    ax1.set_xlabel(r'$S_{\mathrm{GC}}$',fontsize=8.)
#    ax1.set_ylabel(r'$S_{\mathrm{PM}}/S_{\mathrm{GC}}$',fontsize=6.)
#    plt.setp(ax1.get_xticklabels(), fontsize=6.)
#    plt.setp(ax1.get_yticklabels(), fontsize=6.)
#    ax1.xaxis.set_major_locator(MaxNLocator(5))
##    ax1.yaxis.set_major_locator(MaxNLocator(5))
##    ax1.set_xlim(0,1)#np.max(SGC0S)*1.05) 
##    ax1.set_ylim(0.99,1.05)  
#    fig.set_size_inches(1.58,1)
#    fig.subplots_adjust(right=0.99) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.15)##Larger adds whitespace
#    fig.subplots_adjust(left=0.26) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.95) ##Smaller adds whitespace to top
#    plt.savefig('F2_Binset_SPM0_vs_SGC0_b.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    plt.close()

#    print 'Plotting S vs. S for F2...'
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    for (gc,pm,col,mark) in zip(SGC0S,SPM0S,F2_colorS,F2_markerS):
#	    	ax1.errorbar(pm,gc/pm,yerr=None,color=col,marker=mark,ms=3.0,mew=0.1,ls='None')#,label=r'$CS$')#marker=markers[i],color=colors[i]    
##    ax1.errorbar(np.linspace(0,5),np.linspace(0,5),yerr=None,color='k',ls='-')#,label=r'$CS$')#marker=markers[i],color=colors[i]   
##    ax1.set_xlabel(r'$S_{\mathrm{GC}}$',fontsize=8.)
#    ax1.set_ylabel(r'$S_{\mathrm{GC}}/S_{\mathrm{PM}}$',fontsize=6.)
#    plt.setp(ax1.get_xticklabels(), fontsize=6.)
#    plt.setp(ax1.get_yticklabels(), fontsize=6.)
#    ax1.xaxis.set_major_locator(MaxNLocator(5))
##    ax1.yaxis.set_major_locator(MaxNLocator(5))
##    ax1.set_xlim(0,1)#np.max(SGC0S)*1.05) 
##    ax1.set_ylim(0.99,1.05)  
#    fig.set_size_inches(1.5,1)
#    fig.subplots_adjust(right=0.96) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.15)##Larger adds whitespace
#    fig.subplots_adjust(left=0.26) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.95) ##Smaller adds whitespace to top
#    plt.savefig('F2_Binset_SPM0_vs_SGC0_SWITCHAXIS.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    plt.close()


#    print 'Plotting rhof_CONTACT vs. GC (SigmaEff(z))...'
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    fig.subplots_adjust(right=0.90)
#    i=-1
#    labels=[]
#    maxX = -1E9
#    for (Nf,Sig,filename,ratio) in zip(NfS,SigEffS,filenameS,ratioS):
#    	Sig = np.array(Sig)
#    	label_str=''
#    	if filename in GC_sig1:
#    		if legend_key=='sigs':
#    			  label_str = r'$\sigma/\lambda_{D}\approx 0.1$'
# 		elif legend_key=='ele':
# 		   	label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.1$'
#     		i=1
#    	if filename in GC_sig5:
#		label_str = r'$\sigma/\lambda_{D}\approx 0.5$'
#		i=11
#    	if filename in GC_sig6:
#  		label_str = r'$\sigma/\lambda_{D}\approx 0.6$'
#  		i=3
#    	if filename in GC_sig7:
#    		label_str = r'$\sigma/\lambda_{D}\approx 0.7$'
#    		i=6
#    	if filename in GC_sig8:
#    		label_str = r'$\sigma/\lambda_{D}\approx 0.8$'
#    		i=9
#    	if filename in GC_LBhalf:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.05$'
#    		i=0
#    	if filename in GC_LB3:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.3$'
#    		i=3
#    	if filename in GC_LB5:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.5$'
#    		i=5
#    	if filename in GC_LB7:
#    	   	label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.7$'
#    	   	i=7
#    	if filename in GC_LB10:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 1.0$'
#    		i=9
#	if label_str in labels:
#	    	ax1.errorbar(Sig[1:]*np.sqrt((Sig[1:]*0.5)**2+1),Nf[1:-1],yerr=None,color=ROYGBIV_map(ratio,np.max(ratioS)),marker=markers[i],ms=5.0,ls='None')#,label=r'$CS$')#marker=markers[i],color=colors[i]    
#    	else:
#		ax1.errorbar(Sig[1:]*np.sqrt((Sig[1:]*0.5)**2+1),Nf[1:-1],yerr=None,color=ROYGBIV_map(ratio,np.max(ratioS)),marker=markers[i],ms=5.0,ls='None',label=label_str)
#		labels.append(label_str)
#    	maxX = np.max([maxX,np.max(Sig)])    
#    ax1.errorbar(np.linspace(0,50),np.linspace(0,50),yerr=None,color='k',ls='-',label=r'$GC$')
#    ax1.set_xlabel(r'$ \~ \Sigma_{eff} \sqrt{ (\~ \Sigma_{eff}/2)^{2}+1}$',size='x-large')
#    ax1.set_ylabel(r'$\~ \rho_{f}(z \rightarrow 0)$',size='x-large')
#    ax1.legend(loc=0)
#    ax1.set_xlim(0,10) 
#    ax1.set_ylim(0,10) 
##    plt.xlim(xmin=0)
##    plt.ylim(ymin=-0.5)
#    plt.savefig('GC_Contact.pdf', dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
##    plt.show()
#    plt.close()

#    print 'Plotting rhof_CONTACT vs. Bik (SigmaEff(z))...'
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    fig.subplots_adjust(right=0.90)
#    i=-1
#    labels=[]
#    maxX = -1E9
#    for (Nf,Sig,filename,PhiB,ratio) in zip(NfS,SigEffS,filenameS,PhiBulkS,ratioS):
#    	Sig = np.array(Sig)
#    	r = PhiB/0.65
#    	label_str=''
#    	if filename in GC_sig1:
#    		if legend_key=='sigs':
#    			  label_str = r'$\sigma/\lambda_{D}\approx 0.1$'
# 		elif legend_key=='ele':
# 		   	label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.1$'
#     		i=1
#    	if filename in GC_sig5:
#		label_str = r'$\sigma/\lambda_{D}\approx 0.5$'
#		i=11
#    	if filename in GC_sig6:
#  		label_str = r'$\sigma/\lambda_{D}\approx 0.6$'
#  		i=3
#    	if filename in GC_sig7:
#    		label_str = r'$\sigma/\lambda_{D}\approx 0.7$'
#    		i=6
#    	if filename in GC_sig8:
#    		label_str = r'$\sigma/\lambda_{D}\approx 0.8$'
#    		i=9
#    	if filename in GC_LBhalf:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.05$'
#    		i=0
#    	if filename in GC_LB3:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.3$'
#    		i=3
#    	if filename in GC_LB5:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.5$'
#    		i=5
#    	if filename in GC_LB7:
#    	   	label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.7$'
#    	   	i=7
#    	if filename in GC_LB10:
#    		label_str = r'$\lambda_{B}/\lambda_{D}\approx 1.0$'
#    		i=9
#	if label_str in labels:
##		(1./r)*np.sqrt((1.-np.exp(-r*0.5*Sig[1:]**2))*(1-np.exp(-r*0.5*Sig[1:]**2)*(1.-2*r)))
#	    	ax1.errorbar((1./r)*np.sqrt((1.-np.exp(-r*0.5*Sig[1:]**2))*(1-np.exp(-r*0.5*Sig[1:]**2)*(1.-2*r))),Nf[1:-1],yerr=None,color=ROYGBIV_map(ratio,np.max(ratioS)),marker=markers[i],ms=5.0,ls='None')#,label=r'$CS$')#marker=markers[i],color=colors[i]    
#    	else:
#		ax1.errorbar((1./r)*np.sqrt((1.-np.exp(-r*0.5*Sig[1:]**2))*(1-np.exp(-r*0.5*Sig[1:]**2)*(1.-2*r))),Nf[1:-1],yerr=None,color=ROYGBIV_map(ratio,np.max(ratioS)),marker=markers[i],ms=5.0,ls='None',label=label_str)
#		labels.append(label_str)
#    	maxX = np.max([maxX,np.max(Nf)])    
#    ax1.errorbar(np.linspace(0,50),np.linspace(0,50),yerr=None,color='k',ls='--',label=r'$Bikerman$')
#    ax1.set_xlabel(r'$ (1/r) \sqrt{(1-\exp(-r \~ \Sigma_{eff}^{2}/2)[1-\exp(-r \~ \Sigma_{eff}^{2}/2)(1-2r)]}$',size='x-large')
#    ax1.set_ylabel(r'$\~ \rho_{f}(z \rightarrow 0)$',size='x-large')
#    ax1.legend(loc=0)
#    ax1.set_xlim(0,10) 
#    ax1.set_ylim(0,10) 
#    plt.savefig('Bik_Contact.pdf', dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
##    plt.show()
#    plt.close()

#    print 'Plotting similar rhof(z) curves... for GC data only, primarily'
#    #When plotting only GC data, always operating in ND distances is fine...
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    fig.subplots_adjust(right=0.90)
#    i=-1
#    for (rhofi,Sigma_s,z_pos,filename,ratio) in zip(NfS,SIGMAS,zhalfS,filenameS,ratioS):
#    	label_str=''
#    	max_info=[]
#    	if filename in GC_sig1:
#    		for x in Max_effS:
#    			if x[0] in GC_sig1:
#    				max_info=x
#    				i=1
#    				if max_info[3]==Sigma_s:
#			    		if legend_key=='sigs':
#			    			  label_str = r'$\sigma/\lambda_{D}\approx 0.1$'
#			 		elif legend_key=='ele':
#			 		   	label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.1$'
#    				break
#    	if filename in GC_sig5:
#    		for x in Max_effS:
#    			if x[0] in GC_sig5:
#    				max_info=x
#    				i=11
#    				if max_info[3]==Sigma_s:
#    					label_str = r'$\sigma/\lambda_{D}\approx 0.5$'
#    				break
#    	if filename in GC_sig6:
#    		for x in Max_effS:
#    			if x[0] in GC_sig6:
#    				max_info=x
#    				i=3
#    				if max_info[3]==Sigma_s:
#    					label_str = r'$\sigma/\lambda_{D}\approx 0.6$'
#    				break
#    	if filename in GC_sig7:
#    		for x in Max_effS:
#    			if x[0] in GC_sig7:
#    				max_info=x
#    				i=6
#    				if max_info[3]==Sigma_s:
#    					label_str = r'$\sigma/\lambda_{D}\approx 0.7$'
#    				break
#    	if filename in GC_sig8:
#    		for x in Max_effS:
#    			if x[0] in GC_sig8:
#    				max_info=x
#    				i=9
#    				if max_info[3]==Sigma_s:
#    					label_str = r'$\sigma/\lambda_{D}\approx 0.8$'
#    				break
#    	if filename in GC_LBhalf:
#    		for x in Max_effS:
#    			if x[0] in GC_LBhalf:
#    				max_info=x
#    				i=0
#    				if max_info[3]==Sigma_s:
#    					label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.05$'
#    				break
#    	if filename in GC_LB3:
#    		for x in Max_effS:
#    			if x[0] in GC_LB3:
#    				max_info=x
#    				i=3
#    				if max_info[3]==Sigma_s:
#    					label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.3$'
#    				break
#    	if filename in GC_LB5:
#    		for x in Max_effS:
#    			if x[0] in GC_LB5:
#    				max_info=x
#    				i=5
#    				if max_info[3]==Sigma_s:
#    					label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.5$'
#    				break	
#    	if filename in GC_LB7:
#    		for x in Max_effS:
#    			if x[0] in GC_LB7:
#    				max_info=x
#    				i=7
#    				if max_info[3]==Sigma_s:
#    					label_str = r'$\lambda_{B}/\lambda_{D}\approx 0.7$'
#    				break	
#    	if filename in GC_LB10:
#    		for x in Max_effS:
#    			if x[0] in GC_LB10:
#    				max_info=x
#    				i=9
#    				if max_info[3]==Sigma_s:
#    					label_str = r'$\lambda_{B}/\lambda_{D}\approx 1.0$'
#    				break	
#	if len(max_info)!=0:
#		for (z_shift,SigEff) in zip(max_info[1],max_info[2]):
#			if SigEff<=Sigma_s:
#				if label_str=='':
#					ax1.errorbar(np.array(z_pos[1:])+z_shift,np.array(rhofi[1:])+i,yerr=None,color=ROYGBIV_map(ratio,np.max(ratioS)),marker=markers[i],ms=5.0,ls='None')#,label=r'$\Sigma$'+' = ' + str(round(Sigma_s,3)))
#				else:
#					ax1.errorbar(np.array(z_pos[1:])+z_shift,np.array(rhofi[1:])+i,yerr=None,color=ROYGBIV_map(ratio,np.max(ratioS)),marker=markers[i],ms=5.0,ls='--',label=label_str)
#				break	
#	
#    ax1.set_ylabel('Shifted ' + r'$\~ \rho_{f}(z)$',size='x-large')
#    ax1.set_xlabel(r'$z/L_{z}$',size='x-large')
#    ax1.legend(loc=0)
#    ax1.set_xlim(0,0.15)       
#    ax1.set_ylim(-0.1,25)
#    plt.savefig('rhof_z_superimpose.pdf', dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
##    plt.show()
#    plt.close()    		
    
##    ##This plots mu_ex^EV(z) for simulation and modified theories
#    muexEV_bulkS=[]
#    nonDim='yes'
##    Bikerman='yes'
#    Bikerman='no'
##    CarnahanStarling='yes'
#    CarnahanStarling='no'
#    if nonDim=='yes':
#      print 'Plotting NonDim mu_ex^EV(z)...'
#    else:
#      print 'Plotting Dim mu_ex^EV(z)...'	
#    if Bikerman=='yes':
#		print '\t\t...with Bikerman theory...'
#		#Note to user: Bikerman theory must already be evaluated using MATLAB code, which is in ~/sims/PBik_Solver.m
#    if CarnahanStarling=='yes':
#	      print '\t...with Carnahan-Starling theory...'
#    i=-1
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    fig.subplots_adjust(right=0.90)
#    for (z_positions,characteristic_length,muex_EV,L_z) in zip(z_S,characteristic_lengthS,muexEV_S,L_zS):
#    	i+=1

#	if i==len(markers): #This resets the markers
#		i=0

#	if i==0:
#		graphname = 'muex_EV_z' + '.pdf'
#		ylabel=r'$\~\mu_{ex}^{EV}$'
#		
#	z_positions=[z/L_z for z in z_positions] #NonDim distances

#        muexEV_bulk=np.mean(muex_EV[len(muex_EV)/2-2:len(muex_EV)/2+2])
#        muexEV_bulkS.append(muexEV_bulk)

#	##The code below includes bulk values:
#	ax1.errorbar(np.array(z_positions)*(L_z/characteristic_length),np.array(muex_EV),yerr=None,marker=markers[i],ms=3.0,color=colors[i],ls='None')
#	
##	#The code below subtracts off bulk values:
##	ax1.errorbar(np.array(z_positions)*(L_z/characteristic_length),np.array(muex_EV)-muexEV_bulk,yerr=None,marker=markers[i],ms=3.0,color=colors[i],ls='None')
#	
#	if Bikerman=='yes':
#	  Bik_file='_Bik_zeta_' + ze + '.txt'
#	  x_Bik=[[float(x) for x in line.split()] for line in file('x' + Bik_file,"r").readlines()]
#	  x_Bik=np.array(x_Bik)*lam_D/characteristic_length
#	  c_Bik=[[float(x) for x in line.split()] for line in file('counter' + Bik_file,"r").readlines()]
#	  co_Bik=[[float(x) for x in line.split()] for line in file('co' + Bik_file,"r").readlines()]
#	  ax1.plot(x_Bik[0],-0.5*np.log(np.array(co_Bik[0])*np.array(c_Bik[0])),color=colors[i],lw=1.5,ls='--') 
#	if CarnahanStarling=='yes':
#	  CS_file='_CS_zeta_' + ze + '.txt'
#          x_CS=[[float(x) for x in line.split()] for line in file('x'+CS_file,"r").readlines()]
#          x_CS=np.array(x_CS)*lam_D/characteristic_length
#	  c_CS=[[float(x) for x in line.split()] for line in file('counter' + CS_file,"r").readlines()]
#	  co_CS=[[float(x) for x in line.split()] for line in file('co' + CS_file,"r").readlines()]	
#	  ax1.plot(x_CS[0],-0.5*np.log(np.array(co_CS[0])*np.array(c_CS[0])),color=colors[i],lw=1.5,ls='-.')      
#  		##Keep this code for sometime - 02/22/12 14:43:58 
##	  if sig==1:
##	  	ax1.plot(x_CS[0],0.5*np.log(np.array(co_CS[0])*np.array(c_CS[0])),color=colors[i],lw=1.5,ls='-.')
##	  else:
##	  	ax1.plot(x_CS[0],-0.5*np.log(np.array(co_CS[0])*np.array(c_CS[0])),color=colors[i],lw=1.5,ls='-.')
#    if xlabel==r'$z/L_{z}$':
#      #These must be set manually! 
#      test=0
##      ax1.set_xlim(0,0.2)       
##      ax1.set_ylim(0,10) 
#    elif xlabel==r'$z/ \sigma_{WCA}$':
#      #These must be set manually! 
#      test=0
##      ax1.set_xlim(0,0.2)       
##      ax1.set_ylim(0,10)      
#    ax1.set_xlabel(xlabel,size='x-large')
#    ax1.set_ylabel(ylabel,size='x-large')
##    ax1.legend(loc=0) 
#    plt.savefig(graphname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
##    plt.show()
#    plt.close()
##    #Done plotting mu_excess^EV(z)

#    print 'Plotting muexEV = muexEV(Phi_bulk)...'
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
#    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
#    i=-1
#    j=-1
#    PhiWCAs=[]
#    Phi_shifted=[]
#    Phi_HS=[]
#    MAX = -1e9
#    for (mexb,Xi,n0,sig_WCA,sigHS) in zip(muexEV_bulkS,XiS,n0s,sigWCA_S,sigHS_S):
#    	i+=1
#    	j+=1
#	if i==len(markers): #This resets the markers
#			i=0
#    	Phi_WCA = 2.*n0*(np.pi/6)*sig_WCA**3
##    	print '%1.5f\t\t%1.5f' % (Phi_WCA,mexb)
#    	MAX = np.max([MAX,Phi_WCA])
#    	PhiWCAs.append(Phi_WCA)
#    	Phi_shifted.append(Phi_WCA*0.953259114**3)  #Unweighted = 0.981989165  0.954028946 
#    	Phi_HS.append(2.*n0*(np.pi/6)*sigHS**3)
#    	if j==0:
#    		qqq=1
#    		ax1.errorbar(Phi_WCA,mexb,yerr=None,mfc='white',mec='blue',marker='o',ms=3.0,color='blue',ls='None',label=r'$\~ \mu_{ex}^{EV}(\sigma_{WCA})$')
######    		ax1.errorbar(Phi_WCA,mexb,yerr=None,color='blue',ls='-',label=r'$\~ \mu_{ex}^{EV}(\sigma_{WCA})$')
#	else:
#		qqq=1
#    		ax1.errorbar(Phi_WCA,mexb,yerr=None,mfc='white',mec='blue',marker='o',ms=3.0,color='blue',ls='None')
######    		ax1.errorbar(Phi_WCA,mexb,yerr=None,color='blue',ls='-')
#    Phi_theory = np.linspace(0,MAX,1000)
#    ax1.errorbar(Phi_shifted,muexEV_bulkS,yerr=None,color='blue',ls='None',marker='*',ms=6.0,label=r'$\~ \mu_{ex}^{EV}(0.953*\sigma_{WCA})$') 

#    ax1.errorbar(np.array(Phi_shifted)*(0.984796597/0.953259114)**3,muexEV_bulkS,yerr=None,color='orange',ls='None',marker='^',ms=3.0,label=r'$\~ \mu_{ex}^{EV}(0.955*\sigma_{WCA})$') 
#    
#    ax1.errorbar(Phi_HS,muexEV_bulkS,yerr=None,color='blue',marker='o',ms=3.0,ls='None',label=r'$\~ \mu_{ex}^{EV}(\sigma_{HS})$') 
#    ax1.errorbar(Phi_theory,Phi_theory*(8-9*Phi_theory+3*Phi_theory**2)/(1-Phi_theory)**3,yerr=None,color='k',ls='-.',lw=1.5,label=r'$\mu_{ex}^{CS}(\Phi)$')#marker=markers[i],color=colors[i]    
#    ax1.errorbar(Phi_theory,-np.log(1-Phi_theory/0.65),yerr=None,color='r',ls='--',label=r'$-\ln(1-\Phi/0.65)$')#marker=markers[i],color=colors[i]    # /0.65)
#    ax1.errorbar(Phi_theory,-np.log(1-Phi_theory),yerr=None,color='k',ls='--',label=r'$-\ln(1-\Phi)$')#marker=markers[i],color=colors[i]    # /0.65)

#    ax1.set_xlabel(r'$\Phi^{bulk}$',size='x-small')
#    ax1.set_ylabel(r'$\~\mu_{ex}^{EV}(\Phi^{bulk})$',size='x-small')
#	##Important fitting text below
##    coeffs = np.polyfit(PhiWCAs,np.exp(-np.array(muexEV_bulkS)),deg=3)
##    x_fit = np.linspace(0,MAX*1.1,100)
##    print coeffs
##    leg_info = r'$-\ln(%s \Phi_{WCA}^{3} + %s \Phi_{WCA}^{2}  %s\Phi_{WCA}+%s)$' % (str(round(coeffs[0],1)),str(round(coeffs[1],1)),str(round(coeffs[2],2)),str(round(coeffs[3],2)))
##    ax1.plot(x_fit,-np.log(np.polyval(np.poly1d(coeffs),x_fit)),color='b',lw=2.0,ls='-',label=r'$-\ln(1-f[\Phi])$')
##    ymean = np.mean(np.exp(-np.array(muexEV_bulkS)))
##    R_sq = 1. - (sum([(yi-fi)**2 for (yi,fi) in zip(np.exp(-np.array(muexEV_bulkS)),np.polyval(np.poly1d(coeffs),PhiWCAs))]) / sum([(yi-ymean)**2 for yi in np.exp(-np.array(muexEV_bulkS))]))
##    print 'R_sq = ',R_sq
##    ax1.text(1,12,'Fitting \~\mu_{ex}^{EV}(\sigma_{WCA}):\nR^2 = %s' % str(R_sq))

#    plt.setp(ax1.get_xticklabels(), fontsize='x-small')
#    plt.setp(ax1.get_yticklabels(), fontsize='x-small')

#    ax1.set_xlim(0,MAX)       
#    ax1.set_ylim(0,np.max(muexEV_bulkS))
#    
#    fig.set_size_inches(3.37,3.5)
#    plt.savefig('muexEV_vs_Phi_B.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
###    plt.show()
##    ax1.legend(loc=0)#    fig.set_size_inches(11,7)

#    ax1.set_xlim(0,0.05805)
#    ax1.set_ylim(0,0.5) 
#    plt.savefig('muexEV_vs_Phi_B_LT1.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 
#    plt.close()

#    print 'Ratio of sigma/LD'
#    for (x,y) in zip(sigHS_S,LDs):
#    	print 'sig,LD, sig/LD = ',x,y,x/y

###    This plots mu_excess^EC(z) for simulation and GC theory
#    print 'Plotting muex^EC_pm(z)...\n\t***This code must be executed in order to make countour plots!'
#    i=-1
#    fig=plt.figure()
##    ax1=fig.add_subplot(111)
##    fig.subplots_adjust(right=0.90)
#    ax1=fig.add_subplot(211)
#    ax2=fig.add_subplot(212)
#    muex_EC_pS,muex_EC_mS=[],[]
#    Con_ECp,Con_ECm,Con_EV = [[],[]],[[],[]],[[],[]] #Contour plot information, to be determined below
#    ThreeD_Plot=[]
#    muexEC_test=[]
#    for (Nm,Np,characteristic_length,z_positions,Psi,EV,L_z,n0,area,PB,muexEV_bulk,Nf_half,Bjerrum,sigWCA,filename) in zip(N_minusS,N_plusS,characteristic_lengthS,z_S,Volts,muexEV_S,L_zS,n0s,areaS,PhiBulkS,muexEV_bulkS,NfS,BjerrumS,sigWCA_S,filenameS):
#	i+=1
#	
#	if i==0:
#		graphname = 'muex_EC_z' + '.pdf'
#	if i==1:
#		ax1.errorbar(np.array(z_density)*(L_z/characteristic_length),np.zeros(len(z_density)),yerr=None,color='k',ls=':')
#		ax2.errorbar(np.array(z_density)*(L_z/characteristic_length),np.zeros(len(z_density)),yerr=None,color='k',ls=':')
#		ax1.errorbar(np.array(z_density)*(L_z/characteristic_length),np.ones(len(z_density)),yerr=None,color='k',ls='--')
#		ax2.errorbar(np.array(z_density)*(L_z/characteristic_length),-np.ones(len(z_density)),yerr=None,color='k',ls='--')
#		
#	if i==len(markers): #This resets the markers
#		i=0

#	##ND concentrations are required
#	Nm=np.array([nm/(area*L_bin*n0) for nm in Nm])
#	Np=np.array([np/(area*L_bin*n0) for np in Np])
#	Nt=Nm+Np

#	##Need bin centers because that's how Nm and Np are calculated
#	z_density = [(x+y)/(L_z*2.) for (x,y) in zip(z_positions[0:len(z_positions)-1],z_positions[1:len(z_positions)])]

#	import numpy as np #Why the fuck do I need to do this again? - 02/29/12 13:12:37 
#	##Voltages have always been calculated as (Volt - Volt_bulk), but Volt_bulk is strictly zero for acceptable simulation data
#	Psi = np.array([(x+y)/2. for (x,y) in zip(Psi[0:len(Psi)-1],Psi[1:len(Psi)])])  ##Determine voltage in the middle of the bin... since density and voltage have different bins

#	##muexEV - bulk values:
#	EV = np.array([(exV-muexEV_bulk) for exV in EV[:-1]])

#	mu_ex_wall=np.zeros(len(z_density),dtype=float)
#	k=-1
#	for z in [value*L_z for value in z_density]: #z is a dimensional variable of bin centers
#		k+=1
#		if z<=z_wall:
#			zt=z
#			mu_ex_wall[k]=eps_wall*((2./15.)*(1.165/zt)**9.-(1.165/zt)**3.-(2./15.)*((2./5.)**(-1./6.)*1.)**9.+((2./5.)**(-1./6.)*1.)**3.)
#		if z>=(L_z-z_wall):
#		   	zt=z-(L_z-z_wall)
#			mu_ex_wall[k]=eps_wall*((2./15.)*(1.165/zt)**9.-(1.165/zt)**3.-(2./15.)*((2./5.)**(-1./6.)*1.)**9.+((2./5.)**(-1./6.)*1.)**3.)

#	##This may perhaps be made permanent...
##	mu_ex_wall=np.zeros(len(z_density),dtype=float)
#	muex_EC_p=-(np.log(Np)+np.array(Psi)+EV+mu_ex_wall)
#	muex_EC_m=-(np.log(Nm)-np.array(Psi)+EV+mu_ex_wall)

#	muex_EC_p_bulk=np.mean(muex_EC_p[len(muex_EC_p)/2-2:len(muex_EC_p)/2+2])
#	muex_EC_m_bulk=np.mean(muex_EC_m[len(muex_EC_m)/2-2:len(muex_EC_m)/2+2])
#	
#	ax1.errorbar(np.array(z_density)*(L_z/characteristic_length),muex_EC_p-muex_EC_p_bulk,yerr=None,marker='+',ms=7.0,color=colors[i])#,ls='-')#,label=r'$\Sigma$'+' = ' + str(round(Sigma_s,4)))
#	ax2.errorbar(np.array(z_density)*(L_z/characteristic_length),muex_EC_m-muex_EC_m_bulk,yerr=None,marker='_',ms=7.0,color=colors[i])#,ls='--')#,label=r'$\~ \zeta$'+' = ' + str(round(zd,3)))
#	muex_EC_pS.append(muex_EC_p)
#	muex_EC_mS.append(muex_EC_m)
#	
#	#This indicated a wall marker
#	ax1.errorbar(np.ones(100)*(z_wall/characteristic_length),np.linspace(-2,0,100),yerr=None,marker='|',ms=7.0,color='k')#,ls='-')#,label=r'$\Sigma$'+' = ' + str(round(Sigma_s,4)))
#	ax2.errorbar(np.ones(100)*(z_wall/characteristic_length),np.linspace(0,2,100),yerr=None,marker='|',ms=7.0,color='k')#,ls='-')#,label=r'$\Sigma$'+' = ' + str(round(Sigma_s,4)))	

#	##Take average of all profiles, keeping only from electrode to bulk
#	EV = EV.tolist()
#	ECp = muex_EC_p.tolist()
#	ECm = muex_EC_m.tolist()
#	Volt = Psi.tolist()
#	Nt = Nt.tolist()
#	
#	z_pos = [x*L_z for x in z_density[:len(z_density)/2]]
#	z_pos.reverse()

#	EV_right = EV[len(EV)/2:]
#	EV_right.reverse()
#	EV_average = [((x+y)*0.5)+muexEV_bulk for (x,y) in zip(EV[:len(EV)/2],EV_right)]
##	EV_average = [((x+y)*0.5) for (x,y) in zip(EV[:len(EV)/2],EV_right)]
#	EV_average.reverse()

#	Volt_right = Volt[len(Volt)/2:]
#	Volt_right.reverse()
#	Volt_average = [(x-y)*0.5 for (x,y) in zip(Volt[:len(Volt)/2],Volt_right)]	
#	Volt_average.reverse()

#	ECp_right = ECp[len(ECp)/2:]
#	ECp_right.reverse() #This is needed for ECm average
#	ECm_average = [(x+y)*0.5 for (x,y) in zip(ECm[:len(ECm)/2],ECp_right)]
#	ECm_average.reverse()

#	ECm_right = ECm[len(ECm)/2:]
#	ECm_right.reverse() #This is needed for ECp average
#	ECp_average = [(x+y)*0.5 for (x,y) in zip(ECp[:len(ECp)/2],ECm_right)]
#	ECp_average.reverse()

#	Nt_right = Nt[len(Nt)/2:]
#	Nt_right.reverse() #This is needed for ECp average
#	Nt_average = [(x+y)*0.5 for (x,y) in zip(Nt[:len(Nt)/2],Nt_right)]
#	Nt_average.reverse()

#	xaxis,yaxis,zaxis=[],[],[]
#	for (z,mex,phi) in zip(z_pos,[(qq+muexEV_bulk) for qq in EV_average],[ww*n0*(np.pi/6)*(sigWCA)**3 for ww in Nt_average]):
#		z=2*(z-z_wall)/(L_z-2*z_wall)
#		if z>=0:
#			xaxis.append(z)
#			yaxis.append(phi)
#			zaxis.append(mex)
#	ThreeD_Plot.append([xaxis,yaxis,zaxis])

#	if muexEV_bulk>=2:
#		print filename

##	Nf_half = Nf_half.tolist()
##	Nf_half.reverse()
##	for (z,ecP,ecM,nfree) in zip(z_pos,ECp_average,ECm_average,Nf_half):
##		if z>z_wall:
##			muexEC_test.append([nfree,ecP,ecM,Bjerrum])

#	a,b,c=0,0,0
#	setpoint = 0.05
##	print '\n'
#	for (z,ev,ecP,ecM,volt,nt) in zip(z_pos,EV_average,ECp_average,ECm_average,Volt_average,Nt_average):
#		if z>z_wall:

#			if a==0 and ev>=setpoint:
#				print 'z,ev,volt,phi = ',z,ev,volt,nt*n0*(np.pi/6)*(sigWCA*0.953)**3,filename
#				Con_EV[0].append((n0*nt*(np.pi/6)*(sigWCA*0.923)**3)) #originally (n0*nt*(np.pi/6)*(sigWCA*0.953)**3)
#				Con_EV[1].append(volt)
#				a+=1
#			if b==0 and ecP<=-setpoint:
#				Con_ECp[0].append(PB)
#				Con_ECp[1].append(volt)
#				b+=1
#			if c==0 and ecM>=setpoint:
#				Con_ECm[0].append(PB)
#				Con_ECm[1].append(volt)
#				c+=1
#    ax1.set_ylabel(r'$\~\mu_{ex,+}^{nMF}$',size='x-large') 
#    ax2.set_ylabel(r'$\~\mu_{ex,-}^{nMF}$',size='x-large') 
#    if xlabel==r'$z/L_{z}$':
#      #These must be set manually! 
#      test=0
##      ax1.set_xlim(0,0.15)       
#      ax1.set_ylim(-10,10) 
#      ax2.set_xlim(0,0.15)       
#      ax2.set_ylim(-10,10) 
#    elif xlabel==r'$z/ \sigma_{WCA}$':
#      #These must be set manually! 
#      test=0
##      ax1.set_xlim(0,10)       
##      ax1.set_ylim(0,10)      
#    ax2.set_xlabel(xlabel,size='x-large')
##    ax1.legend(loc=0) 
#    plt.savefig(graphname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
##    plt.show()
#    plt.close()
#      #Finished plottting excess chemical potential DUE TO ELECTROSTATICS = f(z)


#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    fig.subplots_adjust(right=0.90)
#    graphname = 'EV_vs_Phi_ofZ.pdf'
#    for (x,col) in zip(ThreeD_Plot,F2_colorS):
#    	ax1.errorbar(x[1],x[2],yerr=None,color=col,ls='-',marker='x')#,label=r'$\tilde \mu_{ex,+}^{nMF}$')
#    ax1.set_ylabel(r'$\mu_{ex}^{EV}(z)$',size='x-large') 
#    ax1.set_xlabel(r'$\Phi (z)$',size='x-large') 
#    plt.savefig(graphname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
##    plt.show()
#    plt.close()


##    ##This generates a 3D plot of (x,y,z) = (distance from electrode,phi,muexEV)
#    print 'Plotting 3D graph of 3D plot of (x,y,z) = (distance from electrode,phi,muexEV)...'
#    fig=plt.figure()
#    ax = Axes3D(fig)
##    fig.subplots_adjust(bottom=0.17)
#    graphname = '3D_muexEV.pdf'	
#    RatioS = [(sig/LB) for (sig,LB) in zip(sigWCA_S,BjerrumS)]
#    i=-1
#    aa,bb=0,0
#    for (sheet,SM,filename) in zip(ThreeD_Plot,SigmaMeasureS,filenameS):
#    	i+=1
##	ax.scatter(sheet[0],sheet[1],sheet[2],color=ROYGBIV_map(np.abs(SM),np.max(SigmaMeasureS)),marker='o')   #  
#	ax.scatter(sheet[0],sheet[1],sheet[2],color='r',marker='o')   #  
#	if 'P1_1000_7' in filename:
#		if aa==0:
#			ax.scatter(sheet[0],sheet[1],sheet[2],color='r',marker='o',label=r'$\~ \Sigma = 0.0888$')   #  
#			aa+=1
#		else:
#			ax.scatter(sheet[0],sheet[1],sheet[2],color='r',marker='o')
#	elif 'P1_1000_5' in filename:
#		if bb==0:
#			ax.scatter(sheet[0],sheet[1],sheet[2],color='r',marker='o',label=r'$\~ \Sigma = 0.0633$')   #  
#			bb+=1
#		else:
#			ax.scatter(sheet[0],sheet[1],sheet[2],color='r',marker='o')
##    for z in np.linspace(0,1,50):
#    ax.scatter(0.954**3*np.ones(len(Phi_theory)),Phi_theory,Phi_theory*(8-9*Phi_theory+3*Phi_theory**2)/(1-Phi_theory)**3,color='k',marker='x')#,label=r'$\mu_{ex}^{CS}(\Phi^{\infty})$')    	
#    ax.set_xlabel(r'$(z-z_{wall})/z_{\infty}$',size='x-large') 
##    ax.set_xlim(0,1.)
#    ax.set_ylabel(r'$\Phi$',size='x-large') 
##    ax.set_ylim(0,0.5)
#    ax.set_zlabel(r'$\mu_{ex}^{EV}/k_{B}T$',size='x-large') 
##    ax.set_zlim(0,10.)
##    ax.legend(loc=0) 
##    plt.show()
#    plt.savefig(graphname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    plt.close()	

####    ##This plots muexEC vs. rho_free... why not? The code works, but the plot does NOT
###    print 'Plotting muexEC vs. rho_free...'
###    print '\tThis ROUGH plot needs refining in terms of legend information etc.'
###    i=0
###    fig=plt.figure()
###    ax1=fig.add_subplot(111)
###    fig.subplots_adjust(right=0.90)
###    graphname = 'muexEC_vs_rhoF.pdf'
###    for x in muexEC_test:
###    	ax1.errorbar(x[0],x[1]/x[3],yerr=None,color=ROYGBIV_map(x[3],np.max(BjerrumS)),ls='None',marker='+',label=r'$\tilde \mu_{ex,+}^{nMF}$')
####	ax1.errorbar(muexEC_test[0],muexEC_test[2],yerr=None,color='r',ls='None',marker='+',label=r'$\tilde \mu_{ex,-}^{nMF}$')
###    ax1.set_ylabel(r'$\~\mu_{ex,\pm}^{nMF}(z)$',size='x-large') 
###    ax1.set_xlabel(r'$\~\rho_{f}(z)$',size='x-large') 
###    plt.savefig(graphname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
####    plt.show()
###    plt.close()

##    ##This plots contours of zeta_threshold vs PhiBulk
#    print 'Plotting countor plot: zeta_fail vs PhiBulk'
#    print '\tThis plot needs refining in terms of legend information etc.\n\tAlso, supporting theory needs to be developed that correctly accounts for bulk values!'
#    i=0
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    fig.subplots_adjust(right=0.90)
#    graphname = 'Contour_zetaF_vs_PhiB.pdf'
#    PhiB_theory = np.linspace(np.min(PhiBulkS)*0.001,np.max(PhiBulkS)*1.1,100)
##    PhiB_theory = np.linspace(0.001,0.65,100)
#    ax1.errorbar(PhiB_theory,np.zeros(len(PhiB_theory)),yerr=None,color='k',ls=':')

#    if len(Con_ECp[0])!=0:
#        ax1.errorbar(Con_ECp[0],Con_ECp[1],yerr=None,marker='+',ms=7.0,color='orange',ls='None')#marker=markers[i],color=colors[i]
#        i+=3
#        print 'Change this sometime'
#    if len(Con_ECm[0])!=0:
#    	ax1.errorbar(Con_ECm[0],Con_ECm[1],yerr=None,marker='_',ms=7.0,color='blue',ls='None')#marker=markers[i],color=colors[i]
#    	i+=1
#    if len(Con_EV[0])!=0:
#    	ax1.errorbar(Con_EV[0],Con_EV[1],yerr=None,marker='s',ms=3.0,color=colors[i],ls='None')
##    	print markers[i]
##    	Phi_star = 0.101839 #CS Fail
##    	ax1.plot(PhiB_theory,np.log((np.e*Phi_star/PhiB_theory)*0.5 + np.sqrt((np.e*Phi_star/PhiB_theory)**2 - 4)*0.5),ls='-.',marker='None',color=colors[i])
##    	Phi_star = (np.e-1)/np.e  #Bik fail
##    	ax1.plot(PhiB_theory,np.log((np.e*Phi_star/PhiB_theory)*0.5 + np.sqrt((np.e*Phi_star/PhiB_theory)**2 - 4)*0.5),ls='--',marker='None',color=colors[i])

#	if setpoint ==1.0:
#    		Phi_star = 0.101839 #CS Fail
#	elif setpoint == 0.5:
#		Phi_star = 0.0560375
#	elif setpoint == 0.25:
#		Phi_star = 0.0295338
#	elif setpoint == 0.10:
#		Phi_star = 0.0122147
#	elif setpoint == 0.05:
#		Phi_star = 0.00617773
#	else:
#		print 'Need to calculated Phi_star for CS! Default value used.'
#		Phi_star = 0.0295338

#	print 'huh',setpoint,Phi_star
#    	ax1.plot(PhiB_theory,np.log((np.exp(setpoint)*Phi_star/PhiB_theory)*0.5 + np.sqrt((np.exp(setpoint)*Phi_star/PhiB_theory)**2 - 4)*0.5),ls='-.',marker='None',color=colors[i])
#    	Phi_star = 1 - np.exp(-setpoint)
#    	ax1.plot(PhiB_theory,np.log((np.exp(setpoint)*Phi_star/PhiB_theory)*0.5 + np.sqrt((np.exp(setpoint)*Phi_star/PhiB_theory)**2 - 4)*0.5),ls='--',marker='None',color=colors[i])

#    	i+=1   
#    ax1.set_xlabel(r'$\Phi$',size='x-large')
#    ax1.set_ylabel(r'$\~\zeta^{*}(\tilde \mu_{ex}^{i}=$'+str(setpoint)+')',size='x-large')
##    ax1.legend(loc=0) 
##    ax1.set_xlim(0,np.max(PhiBulkS)*1.1)       
#    ax1.set_ylim(-0.1,2)  
#    plt.savefig(graphname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
##    plt.show()
#    plt.close()
##    #Done plotting contours


##	##This code creates a files with all necessary MATLAB inputs
#    print 'MATLAB & MMA preparation code commence'
#    print 'x_BIk should be plotted as xBik*(L_z-2*z_wall+z_wall)\n\twith same y_vals'
#    i=-1
#    MATLAB_TEXT = file("MATLAB_COMMANDS.txt","w")
#    MMA_TEXT = file("MMA_COMMANDS.txt","w")
#    for (Volt,sig_HS,n0,ND_Sig,lam_D,L_bin,z_positions,filename,L_z) in zip(Volts,sigHS_S,n0s,ND_SigS,LDs,L_binS,z_S,filenameS,L_zS):
#    	i+=1

#    	Phi_Bik = 2.*n0*sig_HS**3
#    	zL=[]
#    	MMA_string=''
#    	ll=-1
#     	for (V,z) in zip(Volt,z_positions):
#     		ll+=1
#     		if ll==0:
#			MMA_string+='\nSimData={'

#		if z<=L_z*0.5:
#			MMA_string+='{'+str(Phi_Bik/0.65)+','+str(ND_Sig)+','+str(z/lam_D)+'},'
#	
#     		if z>=z_wall and type(zL)==type([]): #The first time z>z_wall, grab voltage
#     			zL=V
#        MMA_string=MMA_string[:-1]+'}'
#	Volt.reverse()
#	z_positions.reverse()
#    	zR=[]
#     	for (V,z) in zip(Volt,z_positions):
#     		if z<=L_z-z_wall and type(zR)==type([]): #The first time z>z_wall, grab voltage
#     			zR=V
#     			break

#	MMA_TEXT.write('%s' % MMA_string)
#	MMA_TEXT.write('\nresult = Map[ynew, SimData]\nExport["~/Desktop/y_Bik/MMA_%s.txt", result, "Table"]' % filename[12:-4])
#		
#	MATLAB_TEXT.write('P1_PBik(%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,\'MATLAB_%s.txt\');\n' % (zL,zR,lam_D/(L_z-2*z_wall),Phi_Bik/0.65,L_bin/(L_z-2*z_wall),filename[12:-4]))
#    MMA_TEXT.close()
#    MATLAB_TEXT.close()		

#    print '\t\tDebugging rho_f(y_Bik) code!'
#    #This plots N_free(z) vs y_Shell_Bik for simulation and Bikerman theory
#    print 'Plotting Nf(ShellY_Bik)...'
#    i=-1
#    fig=plt.figure()
##    fig.subplots_adjust(left=0.16)
#    fig.subplots_adjust(bottom=0.14)
#    ax1=fig.add_subplot(111)
#    for (Volt,Nm,Np,n0,Sigma_s,area,lam_D,rho_fT,L_bin,Bjerrum,z_positions,filename,L_z) in zip(Volts,N_minusS,N_plusS,n0s,SIGMAS,areaS,LDs,rhof_Theory,L_binS,BjerrumS,z_S,filenameS,L_zS):
#		i+=1
#		if i==len(markers): #This resets the markers
#			i=0

#		if i==0:
#			graphname = 'Nf_v_ShellY_Bik'
#			ylabel=r'$\~ \rho_{f}(\~ y^{Bik})$'
#			z_den_T = np.array(np.linspace(0.+z_wall/L_z,1.-z_wall/L_z,2*5000))*L_z

#                dielectric=(4*np.pi*Bjerrum)**-1
#	        z_density = [(x+y)/(L_z*2.) for (x,y) in zip(z_positions[0:len(z_positions)-1],z_positions[1:len(z_positions)])]
#                
##		##This plots my data
#		Nf=np.array([(-npl+nm)/(area*L_bin) for (npl,nm) in zip(Np,Nm)])
#		Nf=0.5*Nf/(dielectric*temperature/(valency*lam_D**2)) ##Assumes temp and valency are 1		

#		
#		ND_Sig = Sigma_s / (dielectric*temperature / (valency*lam_D))			
#		rho_f,y_Shell=[],[]
#		rho_fi=[]
#		f=0
#		index=-1
#		for (nfi,nfree,z) in zip(Nf_idealized,Nf,np.array(z_density)*L_z):
#			index+=1
##			if z<L_z*0.5:
#			if index<len(Nf)*0.5: #Nf better have an even number of elements...
#				rho_f.append(nfree)
#				rho_fi.append(nfi)
#				y_Shell.append(np.exp(-z/lam_D)*(np.sqrt(1+(ND_Sig/2.)**-2)-(ND_Sig/2.)**-1))
#			else:
#				f+=1
#				rho_f[-f]=0.5*(rho_f[-f]-nfree) ##Take average value of these
#				rho_fi[-f]=0.5*(rho_fi[-f]-nfi)				
#		datamin=np.min([datamin,np.min(rho_f)])
#		Sigma_s=ND_Sig 
#		if '_GC_' in filename:
#			ax1.errorbar(y_Shell[1:],rho_f[1:],ls='None',yerr=None,color=colors[i],marker=markers[i],label=r'$\tilde \Sigma$'+' = ' + str(round(Sigma_s,3)))
#		else:
#			ax1.errorbar(y_Shell,rho_f,ls='None',yerr=None,color=colors[i],marker=markers[i],label=r'$\tilde \Sigma$'+' = ' + str(round(Sigma_s,3)))


#		#The code below generates theory lines
#		z_ND = z_den_T[z_den_T<L_z*0.5]/lam_D
#		y = np.exp(-z_ND)*(np.sqrt(1+(ND_Sig/2.)**-2)-(ND_Sig/2.)**-1) #This is GC
#		ax1.errorbar(y,-2*y*(y**2+1)/(1-y**2)**2,yerr=None,ls='-',lw=2.5,color='k',marker='None')
##			print 'Testing to see if there is a missing factor'
##			ax1.errorbar(y,-(np.pi/2)*2*y*(y**2+1)/(1-y**2)**2,yerr=None,ls='-',lw=2.5,color='red',marker='None')
#		xmax=np.max(y)

##		if i2==len(filenameS)-1:
#			#The code below gives Bikerman predictions, which are about the same as GC
##			r = 0.5*(dielectric*temperature/(valency*lam_D**2))*((np.pi/6)*sig_HS**3)/0.74048
##			y2 = (np.sqrt(1+(ND_Sig/2.)**-2)-(ND_Sig/2.)**-1)*np.exp(r*(2.+ND_Sig**2)/(2.*np.sqrt(4.+ND_Sig**2)) - z_ND)
##			Psi_Bik = 2.*np.arcsinh(np.sqrt((np.exp(0.5*r*ND_Sig**2)-1.)/(2.*r)))- ND_Sig*z_ND + (1./(2.*r)) * np.sqrt((1.-np.exp(-0.5*r*ND_Sig**2))*(1.-np.exp(-0.5*r*ND_Sig**2)*(1.-2.*r)))*z_ND**2
##			rf_Bik=-np.sinh(Psi_Bik) / ( 1. + 2.*r*np.sinh(Psi_Bik/2.)**2)
##			y = y2
##			ax1.errorbar(y,-2*y*(y**2+1)/(1-y**2)**2,yerr=None,ls='-',lw=2.5,color='red',marker='None')
##			xmax=np.max([np.max(y2),xmax])
#    ax1.set_xlabel(r'$\~y^{GC}=\exp(-z/\lambda_{D}) \left [ \sqrt{2\tilde \Sigma_{s}^{-2} + 1} - 2\tilde \Sigma_{s}^{-1} \right ]$',size='x-large')
#    ax1.set_xlim(0,xmax*1.5)
#    ax1.set_ylim(ymax=0)
#    ax1.set_ylabel(ylabel,size='x-large')
#    plt.savefig(graphname+'_nl.pdf', dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
##    ax1.legend(loc=0) 
##    plt.savefig(graphname+'_l.pdf', dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)          
##    plt.show()
#    plt.close()    
#	###Done with rho_f(Shell_yBik)
#    print 'Done debugging rhof(Shell_Y)'

#    print 'Plotting lateral distribution functions...'
#    graphname = 'LDFs' + '.pdf'
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    i=-1
#    for (filename,HS) in zip(filenameS,sigHS_S):  
#    	CDF_file = 'LDF_50_' + filename[12:] 	
#    	i+=1
#    	z_50=np.array([[x for x in line.split()] for line in file(CDF_file,"r").readlines()])

#    	if i==len(colors): #This resets the colors
#    		i=0

#	r_bins=z_50[:,0]
#	r_bins=[float(x)/HS for x in r_bins[1:]]
#        
#	g_AA=z_50[:,1]
#	g_AA=[float(x) for x in g_AA[1:]]

#	g_AC=z_50[:,2]
#	g_AC=[float(x) for x in g_AC[1:]]

#	g_CC=z_50[:,3]
#	g_CC=[float(x) for x in g_CC[1:]]
#	ax1.errorbar(r_bins,g_CC,yerr=None,ls='-',lw=2.0,color=colors[i])#,label=r'$\Sigma$'+' = ' + str(round(Sigma_s,4)))
##	ax1.errorbar(r_bins,g_AA,yerr=None,ls='-',lw=2.0,color='k')#,label=r'$\Sigma$'+' = ' + str(round(Sigma_s,4)))
##	ax1.errorbar(r_bins,g_AC,yerr=None,ls='--',lw=2.0,color='k')#,label=r'$\Sigma$'+' = ' + str(round(Sigma_s,4)))	
#    ax1.plot(r_bins,np.ones(len(r_bins)),ls=':',color='k')
#    ax1.set_xlabel(r'$r/\sigma_{HS}$',size='x-large') 
#    ax1.set_ylabel(r'$g_{++}$',size='x-large')
##    plt.title(r'$\lambda_B = $'+str(Bjerrum))
#    plt.xlim(xmin=0.0)
##    plt.xlim(xmax=np.max(r_bins))
#    #plt.ylim(ymin=-0.05)
##    plt.ylim(ymax=2.)  
##    ax1.legend(loc=0) 
#    plt.savefig(graphname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
##    plt.show()
#    plt.close()


#    print 'Determining k-analysis on muex^EC_pm(z)...'
#    i=-1
#    fig=plt.figure()
##    ax1=fig.add_subplot(111)
##    fig.subplots_adjust(right=0.90)
#    ax1=fig.add_subplot(211)
#    ax2=fig.add_subplot(212)
#    for (ex_p,ex_m,Sigma_s,zd) in zip(muex_EC_pS,muex_EC_mS,SIGMAS,zeta_data):
#	i+=1
#	
#	if i==0:
#		graphname = 'k_of_z' + '.pdf'
#		exEC_max_lp=np.array(muex_EC_pS[-1])
#		exEC_max_lm=np.array(muex_EC_mS[-1])
#		ihalf=len(exEC_max_lp)/2
#		print ihalf ##Need this to only perform analysis for half the box!
#		

#		
##		ax1.errorbar(np.array(z_plot)*(L_z/characteristic_length),np.zeros(len(z_plot)),yerr=None,color='k',ls=':')
##		ax2.errorbar(np.array(z_plot)*(L_z/characteristic_length),np.zeros(len(z_plot)),yerr=None,color='k',ls=':')
#		
#	if i==len(markers): #This resets the markers
#		i=0

#	muex_EC_p_bulk=np.mean(ex_p[len(ex_p)/2-2:len(ex_p)/2+2])
#	muex_EC_m_bulk=np.mean(ex_m[len(ex_m)/2-2:len(ex_m)/2+2])

#	print muex_EC_p_bulk,muex_EC_m_bulk

#	ax1.errorbar(np.array(z_density)*(L_z/characteristic_length),(ex_p-muex_EC_p_bulk)/exEC_max_lp,yerr=None,marker='+',ms=7.0,color=colors[i],ls='-',label=r'$\Sigma$'+' = ' + str(round(Sigma_s,4)))
#	ax2.errorbar(np.array(z_density)*(L_z/characteristic_length),(ex_m-muex_EC_m_bulk)/exEC_max_lm,yerr=None,marker='_',ms=7.0,color=colors[i],ls='--',label=r'$\~ \zeta$'+' = ' + str(round(zd,3)))
#    ax1.set_ylabel(r'$CHANGE THESE \~\mu_{ex,+}^{EC}$',size='x-large') 
#    ax2.set_ylabel(r'$CHANGE THESE \~\mu_{ex,-}^{EC}$',size='x-large') 
#    if characteristic_length==L_z:
#      ax1.set_xlabel(r'$z/L_{box}$',size='x-large') 
#      plt.xlim(xmin=0.0)
#      plt.xlim(xmax=0.25)
##      plt.ylim(ymin=-0.75)
##      plt.ylim(ymax=3.0)
#    elif characteristic_length==sig:
#      ax1.set_xlabel(r'$z/ \sigma_{WCA}$',size='x-large')
#      ax2.set_xlabel(r'$z/ \sigma_{WCA}$',size='x-large')
##      ax1.set_xlim(0,40)
##      ax1.set_ylim(-4,1)
##      ax2.set_xlim(0,40)
##      ax2.set_ylim(-1,2)
#      if sig==8:
#      	ax1.set_xlim(0,10)
##      	ax2.set_xlim(0,10)
##    ax1.legend(loc=0) 
##    ax2.legend(loc=0) 
#    plt.savefig(graphname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    #plt.show()
#    plt.close()
#      #Finished plottting excess chemical potential DUE TO ELECTROSTATICS = f(z)

###   ##Bikerman MATLAB commands.
##    for (ze,zL,zR,lam_D,n0) in zip(zetas,zeta_dataL,zeta_dataR,LDs,n0s):
##	  co_text='\'co_Bik_at_simZ_zeta_'+ze+'.txt\''
##	  c_text='\'counter_Bik_at_simZ_zeta_'+ze+'.txt\''
##	  volt_text='\'volt_Bik_at_simZ_zeta_'+ze+'.txt\''
##	  print 'PBik_Solver_zLzR_RMS(%1.5f,%1.5f,%1.5f,%1.5f,%s,%s,%s);\n' % (zL,zR,lam_D,2*n0*sig_HS**3,co_text,c_text,volt_text)

###    #Carnahan-Starling commands
###    i=-1
##    for (ze,zL,zR,lam_D,n0) in zip(zetas,zeta_dataL,zeta_dataR,LDs,n0s):
##          i+=1
##	  co_text='\'co_CS_at_simZ_zeta_'+ze+'.txt\''
##	  c_text='\'counter_CS_at_simZ_zeta_'+ze+'.txt\''
##	  volt_text='\'volt_CS_at_simZ_zeta_'+ze+'.txt\''
##	  if i==0:
##	  	  print 'MPB_CS_RMS(%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%s,%s,%s);\n' % (zL,zR,500.,lam_D, 2*n0*(np.pi/6.)*sig_HS**3,co_text,c_text,volt_text)
##	  else:
##	  	  print 'MPB_CS_RMS(%1.5f,%1.5f,%1.5f,%1.5f,%1.5f,%s,%s,%s);\n' % (zL,zR,500.,lam_D, 2*n0*(np.pi/6.)*sig_HS**3,co_text,c_text,volt_text)

##    ##MF Carnahan-Starling commands
##    phi_bulkS=[]
##    for (ze,n0,LD) in zip(zetas,n0s,LDs):
##        phi_bulkS.append(2*n0*(np.pi/6)*sig_HS**3)
##    print 'MF_CS_Compare(%1.5f,%1.5f,%1.1f,%1.1f);\n' % (np.mean(phi_bulkS),np.mean(LDs),Bjerrum,sig)

    return
#P1_Plot()


def P1_CDF():
    """Calculates lateral distribution functions from LAMMPS dump files.
    This demands that L_x = L_y
    It also is ONLY for CDF of a sheet that is up against the wall.
Input:
    None:
    	A filename argument is provided from the command line
Output:
    None, a file of the results is created the directory the code was ran.
""" 
#    Hard-coded values that cannot be determined from file
    r_wall=1.0 #This actually can be determined from the filename r_wall = (L_z - LZ from LAMMPS file)/2
#    z_lo=-251.
#    z_hi=-250+0.55*2*r_ion
    r_AA,r_AC,r_CC=[],[],[]
    NA,NC=0,0
    A,C=[],[]


   #These are the necessary libraries for this function.
    	#Check to see if these are on the machines this file is being run on...
    import gzip
    import sys
    import math
    filename = sys.argv[1]
    filename_orignal=filename

    #Extract relevant parameters from filename
    j=0
    value=''
    for letter in filename:
    	if letter!='_':
    		value+=letter
	else:
		j+=1
		if j==1:
			N_tot=int(value)
		if j==2:
			Xi=float(value)
		if j==3:
			Bjerrum=float(value)
		if j==4:
			sig_WCA=float(value)
			r_ion=0.5*sig_WCA
		if j==5:#This length will actually be determined from LAMMPS files
			L_z=float(value) 
			z_lo=-0.5*L_z-r_wall
			z_hi=-0.5*L_z+0.55*2*r_ion
			print z_lo,z_hi
		value=''

    if filename[-3:]=='.gz': #This length will actually be determined from LAMMPS files
        L_xy=float(value[:-3])
    else:
        L_xy=float(value[:-4])

##	Except for r_ion, z_lo, and z_hi none of the above variables are kept.

    Numbins=50
    print "Number of bins = %i\n\t**This may need to be manually changed!" % Numbins

    if filename[-3:]=='.gz':
        Data=gzip.GzipFile(fileobj=open(filename,'rb'))
	filename=filename[0:len(filename)-2]+'txt'
    else:
        Data=file(filename,"r")

    print '\n'
    print filename
    graph_prefix=filename[0:len(filename)-4]
    print graph_prefix,' This needs to be changed?'

    nth=1     # new addition, untested - 02/10/10 13:17:25 , also unused - 11/18/10 11:47:02    
    n_obs=0
    i,k=0,0
    get_box_dims=0
    for y in Data:
      x=y.split()
           
      #Acquire box dimensions
      get_box_dims+=1
      if get_box_dims<=8:
	if get_box_dims==6:
	  L_x=float(x[1])-float(x[0])
	elif get_box_dims==7:
	  L_y=float(x[1])-float(x[0])
	  area=L_x*L_y
	elif get_box_dims==8:
	  L_box=float(x[1])-float(x[0])
	  z_wall=-0.5*L_box+r_wall
	  L_xy=L_x
	  r_max=0.5*L_xy ##99% sure this is the old way <- Huh? - 02/28/12 17:11:25 
	  V=(z_hi-z_lo-r_wall)*L_xy**2 #The volume of the "control volume" minus the wall volume
	  
      if not x[0]=="ITEM:":
	x=[float(xx) for xx in x]

	if len(x)<=2:
	  i+=1
	  if i==6: #This means an entire timestep worth of data is had
	    n_obs+=1
	    if (n_obs % 1000) == 0:
	    	print 'Calculated %i time steps, percent to completion = %1.2f' % (n_obs,1-float(n_obs)/50000)
	    
	    for (i,ri) in enumerate(A):
    		for (j,rj) in enumerate(A): 
    			if i!=j:
   				r=ri[0:2]-rj[0:2]
   				r=r-L_xy*np.round_(r/L_xy)
   				r=np.sqrt(r[0]**2+r[1]**2)
   				if r<=r_max:
    					r_AA.append(r) #This gives AA and must be changed				
    		for (j,rj) in enumerate(C):
			r=ri[0:2]-rj[0:2]
			r=r-L_xy*np.round_(r/L_xy)
			r=np.sqrt(r[0]**2+r[1]**2)
   			if r<=r_max:
				r_AC.append(r) 	 #This gives AC and must be changed			
	    for (i,ri) in enumerate(C):
    		for (j,rj) in enumerate(C):
    			if i!=j:
   				r=ri[0:2]-rj[0:2]
   				r=r-L_xy*np.round_(r/L_xy)
   				r=np.sqrt(r[0]**2+r[1]**2)
   				if r<=r_max:
					r_CC.append(r) #This gives CC and must be changed	
	    i=1
	    A,C=[],[]
	elif x[1]==1.0:
	      if x[4]<=z_hi and x[4]>z_lo:
		xyz=[x[2],x[3]]
		A.append(np.array(xyz))
		NA+=len(A)
	elif x[1]==2.0:
	      if x[4]<=z_hi and x[4]>z_lo:
		xyz=[x[2],x[3]]
		C.append(np.array(xyz))	  
		NC+=len(C)
			    				
    CDF_output=file('LDF'+'_'+ str(Numbins) + '_' + graph_prefix + ".txt","w")

    NA=NA/float(n_obs)
    NC=NC/float(n_obs)

    g_AA,g_AC,g_CC=[],[],[]
    hi,lo=z_hi,z_lo
    hist,bin=np.histogram(np.array(r_AA), Numbins, (0,r_max))
    bin_width=(bin[1]-bin[0])
    for (i,ci) in zip(np.arange(len(hist)),hist):   	
	g_AA.append((ci/(n_obs*NA*(NA-1.)))*(V/(np.pi*(hi-lo)*(2.*i+1.)*bin_width**2.)))    	
    hist,bin=np.histogram(np.array(r_CC), Numbins, (0,r_max))
    for (i,ci) in zip(np.arange(len(hist)),hist):   	
    	g_CC.append((ci/(n_obs*NC*(NC-1.)))*(V/(np.pi*(hi-lo)*(2.*i+1.)*bin_width**2.))) 
    hist,bin=np.histogram(np.array(r_AC), Numbins, (0,r_max))
    for (i,ci) in zip(np.arange(len(hist)),hist):
    	g_AC.append((ci/(n_obs*NA*NC))*(V/(np.pi*(hi-lo)*(2.*i+1.)*bin_width**2.)))

    g_AA_bulk=np.mean(g_AA[-5:len(g_AA)-1])
    g_CC_bulk=np.mean(g_CC[-5:len(g_CC)-1])
    g_AC_bulk=np.mean(g_AC[-5:len(g_AC)-1])

    g_AA=[x/g_AA_bulk for x in g_AA]
    g_CC=[x/g_CC_bulk for x in g_CC]
    g_AC=[x/g_AC_bulk for x in g_AC]
    bin=[(x+y)/2. for (x,y) in zip(bin[0:len(bin)-1],bin[1:len(bin)])] #Not sure what I feel about this...

    CDF_output.write("Bin		g_AA		g_AC		g_CC\n") 
    for (r,AA,AC,CC) in zip(bin,g_AA,g_AC,g_CC):
      CDF_output.write("%1.3f		%1.3f		%1.3f		%1.3f\n" % (r,AA,AC,CC))
    CDF_output.close()
             
    return

    
def MyHist(Pos,bin):
    """Returns a histogram of the z-positions for ions.
    This does the exact thing as np.histogram() before numpy was updated!
Input:
    Pos:  (N,many) sorted array of atomic positions
    bin:  vector of bin dimensions
Output:
	counts:	a histogram of counts for z-values within their respective bins
"""
    #Start with zeroed bins
    counts= np.zeros(len(bin)-1, int)
    k_store=0
    for x in sorted(Pos): #Scroll through all the ion positions of Pos
        k,get_next_x=k_store,False
        while (k<(len(bin)-1) and not(get_next_x)): #iterate through all the bin sets
	  if x>=bin[k] and x<bin[k+1]: #check to see if the z-position is within a particular bin
	    counts[k]+=1
	    k=k_store
	    get_next_x=True
	  k+=1 #iterare the bin index
    return counts

def Calc_Volt(z_positions,z_plus,z_minus,Sigma_s,epsilon,A_xy):
    """
	Unclear if I need this. Keep it, but comment out the function in the Interogation function
    
Input:
    z_positions:   Bin discretation
    z_plus:	   array of z positions for cations (attracted to right plate)
    z_minus:	   array of z positions for anions (attracted to left plate)
    Sigma_s:	   the applied surface charge density
    epsilon:	   dialectric, or inverse bjerrum
    A_xy:	   Cross-sectional area of simulation box
Output:
    Psi_tot:	   Final calculation of the voltage
    E_tot:	   Final calculation of Efield
"""
  #Initialize array of potential and field   
    E_tot=np.ones(len(z_positions), dtype=float)*(Sigma_s/epsilon)
    Psi_tot=[-Sigma_s/epsilon*z for z in z_positions]
    charge=-1.	##This is a trial switch! - 02/21/12 14:43:20 
    for (i,z) in enumerate(z_positions):
	for z_charge in z_plus[z_plus<=z]:
	    Psi_tot[i]=Psi_tot[i]-charge*(z-z_charge)/(epsilon*A_xy)
	    E_tot[i]=E_tot[i]-charge/(epsilon*A_xy)
    charge=1.
    for (i,z) in enumerate(z_positions):
	for z_charge in z_minus[z_minus<=z]:
	    Psi_tot[i]=Psi_tot[i]-charge*(z-z_charge)/(epsilon*A_xy)
	    E_tot[i]=E_tot[i]-charge/(epsilon*A_xy)
	    
    return np.array(Psi_tot),np.array(E_tot)

def quickCalc_VoltandField(z_positions,z_plus,z_minus,Sigma_s,epsilon,A_xy,valence):
    """
    This calculates the following this using Fortran (f2py)
Input:
    z_positions:   Bin discretation
    z_plus:	   array of z positions for cations (attracted to right plate)
    z_minus:	   array of z positions for anions (attracted to left plate)
    Sigma_s:	   the applied surface charge density
    epsilon:	   dialectric, or inverse bjerrum
    A_xy:	   Cross-sectional area of simulation box
Output:
    Psi_tot:	   Final calculation of the voltage
    E_tot:	   Final calculation of Efield
"""

  #Initialize array of potential and field   
    E_tot=np.ones(len(z_positions), dtype=float)*(Sigma_s/epsilon)
    Psi_tot=np.array([-Sigma_s/epsilon*z for z in z_positions])

    Psi_tot,E_tot = EDLCalcs.calcvoltandfield(np.array(z_positions), np.array(z_plus), np.array(z_minus), epsilon, A_xy, valence, Psi_tot, E_tot, len(z_positions), len(z_plus))	    
    return np.array(Psi_tot),np.array(E_tot)


def Calc_muex_EV(z_positions,Pos_A,Pos_C,M,r_ion,eps_WCA,A_xy):
	    """

	Unclear if I need this. Keep it, but comment out the function in the Interogation function
	    
	Input:
	    z_positions:   Bin discretation
	    Pos_A:	   array of (x,y,z) positions for ions (attracted to right plate)
	    Pos_C:	   array of (x,y,z) positions for ions (attracted to left plate)
	    M:		   the number of random insertions per bin
	    r_ion:	   WCA radius for ions
	    eps_WCA:	   Energy scale for WCA potential
	    A_xy:	   Cross-sectional area of simulation box
	Output:
	    exp_muex_EV:	   Excluded volume contributions to mu_ex 
	    			NOTE: muex_EV[z] = -ln(exp_muex_EV[z]/n_obs) must be taken before printing 
	"""
	    import random
	    L_xy = A_xy**0.5
	    hL_xy=L_xy*0.5
	    iL_xy=1./L_xy
	    sig_WCA=2*r_ion
	    sig_WCA_sq=sig_WCA**2
	    
	  #Initialize array of exp(muex_EV), which we later average 
	  #	and subesequently take the natural log of to get muex_EV(z)
	    exp_muex_EV=np.zeros(len(z_positions), dtype=float)

	    i=-1 #This index is used to fill muex_EV[i]
	    for (zlo,zhi) in zip(z_positions[0:-1],z_positions[1:len(z_positions)]):
		i+=1
	    	m=0 #This index and while loop to ensure that 
	    	
		#Trim away ions whose z-coordinates are +/- sig_WCA from bin edges	
	    	Pos_near_A = Pos_A[np.logical_and(Pos_A[:,2]>(zlo-sig_WCA),Pos_A[:,2]<(zhi+sig_WCA))]
		Pos_near_C = Pos_C[np.logical_and(Pos_C[:,2]>(zlo-sig_WCA),Pos_C[:,2]<(zhi+sig_WCA))]
		#This can't be done apriori for (x,y), since these are periodic dimensions
	
		if (len(Pos_near_A)+len(Pos_near_C))!=0: #Perform M random 
				#insertions if simulation particles have
				#nearby z-coordinated.
			while m < M:
				m+=1
				#Reset energies
				del_Energies_A,del_Energies_C=[],[] 
				
				#Randomly choose (x,y,z) of test ion with a z
				#	position inside the bin of interest
				r0=np.array([random.uniform(-hL_xy,hL_xy),random.uniform(-hL_xy,hL_xy),random.uniform(zlo,zhi)])

				#Calculate the distance between the test ion and Anions
				if len(Pos_near_A)!=0: #provided anions are near
					r_0j=r0-Pos_near_A
					r_0j[:,:2]=r_0j[:,:2]-L_xy*np.round_(r_0j[:,:2]*iL_xy)
					r_0j_sq=np.sum(r_0j*r_0j,axis=1)
					sig_over_r_p6=(sig_WCA_sq/r_0j_sq[r_0j_sq<sig_WCA_sq])**3
					#Calculate energy due to inserting test ion
					del_Energies_A = 4.*eps_WCA*(sig_over_r_p6*(sig_over_r_p6-1)+0.25)
				else:
					del_Energies_A=[]

				if len(Pos_near_C)!=0:
					#Do the same thing for Cations...
					r_0j=r0-Pos_near_C
					r_0j[:,:2]=r_0j[:,:2]-L_xy*np.round_(r_0j[:,:2]*iL_xy)
					r_0j_sq=np.sum(r_0j*r_0j,axis=1)
					sig_over_r_p6=(sig_WCA_sq/r_0j_sq[r_0j_sq<sig_WCA_sq])**3
					del_Energies_C = 4.*eps_WCA*(sig_over_r_p6*(sig_over_r_p6-1)+0.25)
				else:
					del_Energies_C=[]

		    		exp_muex_EV[i]+=np.exp(-(sum(del_Energies_A)+sum(del_Energies_C)))
		else:
			del_Energies_A,del_Energies_C=[],[]
			exp_muex_EV[i]+=M*np.exp(-(sum(del_Energies_A)+sum(del_Energies_C)))
	    return exp_muex_EV/M #Calculate average after M insertions

def quick_Calc_muex_EV(z_positions,Pos_A,Pos_C,M,r_ion,eps_WCA,A_xy):
	    """
	Input:
	    z_positions:   Bin discretation
	    Pos_A:	   array of (x,y,z) positions for ions (attracted to right plate)
	    Pos_C:	   array of (x,y,z) positions for ions (attracted to left plate)
	    M:		   the number of random insertions per bin
	    r_ion:	   WCA radius for ions
	    eps_WCA:	   Energy scale for WCA potential
	    A_xy:	   Cross-sectional area of simulation box
	Output:
	    exp_muex_EV:	   Excluded volume contributions to mu_ex 
	    			NOTE: muex_EV[z] = -ln(exp_muex_EV[z]/n_obs) must be taken before printing 
	"""
	    Pos = Pos_A + Pos_C 
	    exp_EV = np.zeros(len(z_positions), dtype=float)
	    exp_EV_HS = np.zeros(len(z_positions), dtype=float)


		##Transfer this to P1 for P1_data...
	    return EDLCalcs.calc_ev(z_positions, np.array(Pos), exp_EV, exp_EV_HS, M, np.sqrt(A_xy), eps_WCA, 2.*r_ion, len(z_positions), len(Pos))


def Calc_2_muex_EV(z_positions,Pos_A,Pos_C,M,r_ion,eps_WCA,A_xy):
	    """
	    
	Unclear if I need this. Keep it, but comment out the function in the Interogation function
		    
	Input:
	    z_positions:   Bin discretation
	    Pos_A:	   array of (x,y,z) positions for ions (attracted to right plate)
	    Pos_C:	   array of (x,y,z) positions for ions (attracted to left plate)
	    M:		   the number of random insertions per bin
	    r_ion:	   list() of WCA radii for ions r_ion = [smaller (type 1 or A),larger]
	    eps_WCA:	   Energy scale for WCA potential
	    A_xy:	   Cross-sectional area of simulation box
	Output:
	    exp_muex_EV:	   Excluded volume contributions to mu_ex 
	    			NOTE: muex_EV[z] = -ln(exp_muex_EV[z]/n_obs) must be taken before printing 
	"""
	    import random
	    L_xy = A_xy**0.5
	    hL_xy=L_xy*0.5
	    iL_xy=1./L_xy
	    sig_WCA=2*np.array(r_ion)
	    sig_WCA_sq=sig_WCA**2
	    sig_WCA_12 = 0.5*r_ion[0] + 0.5*r_ion[1]
	    sig_WCA_12_sq = sig_WCA_12**2
	   
	  #Initialize array of exp(muex_EV), which we later average 
	  #	and subesequently take the natural log of to get muex_EV(z)
	    exp_1muex_EV=np.zeros(len(z_positions), dtype=float)
	    exp_2muex_EV=np.zeros(len(z_positions), dtype=float)

	    i=-1 #This index is used to fill muex_EV[i]
	    for (zlo,zhi) in zip(z_positions[0:-1],z_positions[1:len(z_positions)]):
		i+=1
	    	m=0 #This index and while loop to ensure that 
	    	
		#Trim away ions whose z-coordinates are +/- larger sig_WCA from bin edges	
	    	Pos_near_A = Pos_A[np.logical_and(Pos_A[:,2]>(zlo-sig_WCA[1]),Pos_A[:,2]<(zhi+sig_WCA[1]))]
		Pos_near_C = Pos_C[np.logical_and(Pos_C[:,2]>(zlo-sig_WCA[1]),Pos_C[:,2]<(zhi+sig_WCA[1]))]
		#This can't be done apriori for (x,y), since these are periodic dimensions
	
		if (len(Pos_near_A)+len(Pos_near_C))!=0: #Perform M random 
				#insertions if simulation particles have
				#nearby z-coordinated.
			while m < M:
				m+=1
				#Reset energies
				del_Energies_A1,del_Energies_C1=[],[] #Small test ions
				del_Energies_A2,del_Energies_C2=[],[] #Larger test ions

				#Randomly choose (x,y,z) of test ion with a z
				#	position inside the bin of interest
				r0=np.array([random.uniform(-hL_xy,hL_xy),random.uniform(-hL_xy,hL_xy),random.uniform(zlo,zhi)])

				#Calculate the distance between the test ions and small Anions
				if len(Pos_near_A)!=0: #provided anions are near
					r_0j=r0-Pos_near_A
					r_0j[:,:2]=r_0j[:,:2]-L_xy*np.round_(r_0j[:,:2]*iL_xy)
					r_0j_sq=np.sum(r_0j*r_0j,axis=1)

					#Calculate energy due to inserting smaller test ion atop small ion					
					sig_over_r_p6=(sig_WCA_sq[0]/r_0j_sq[r_0j_sq<sig_WCA_sq[0]])**3
					del_Energies_A1 = 4.*eps_WCA*(sig_over_r_p6*(sig_over_r_p6-1)+0.25)

					#Calculate energy due to inserting larger test ion atop small ion
					sig_over_r_p6=(sig_WCA_12_sq/r_0j_sq[r_0j_sq<sig_WCA_12_sq])**3
					del_Energies_A2 = 4.*eps_WCA*(sig_over_r_p6*(sig_over_r_p6-1)+0.25)
				else:
					del_Energies_A1,del_Energies_A2=[],[]

				if len(Pos_near_C)!=0:
					#Do the same thing for Cations...
					r_0j=r0-Pos_near_C
					r_0j[:,:2]=r_0j[:,:2]-L_xy*np.round_(r_0j[:,:2]*iL_xy)
					r_0j_sq=np.sum(r_0j*r_0j,axis=1)

					#Calculate energy due to inserting smaller test ion atop larger ion					
					sig_over_r_p6=(sig_WCA_12_sq/r_0j_sq[r_0j_sq<sig_WCA_12_sq])**3
					del_Energies_C1 = 4.*eps_WCA*(sig_over_r_p6*(sig_over_r_p6-1)+0.25)

					#Calculate energy due to inserting larger test ion atop small ion
					sig_over_r_p6=(sig_WCA_sq[1]/r_0j_sq[r_0j_sq<sig_WCA_sq[1]])**3
					del_Energies_C2 = 4.*eps_WCA*(sig_over_r_p6*(sig_over_r_p6-1)+0.25)
				else:
					del_Energies_C1,del_Energies_C2=[],[]

				##MOdify this now~!!!
		    		exp_1muex_EV[i]+=np.exp(-(sum(del_Energies_A1)+sum(del_Energies_C1)))
		    		exp_2muex_EV[i]+=np.exp(-(sum(del_Energies_A2)+sum(del_Energies_C2)))
		else:
			del_Energies_A1,del_Energies_C1=[],[]
			del_Energies_A2,del_Energies_C2=[],[]
			
	    		exp_1muex_EV[i]+=np.exp(-(sum(del_Energies_A1)+sum(del_Energies_C1)))
	    		exp_2muex_EV[i]+=np.exp(-(sum(del_Energies_A2)+sum(del_Energies_C2)))
	    		
	    return exp_1muex_EV/M,exp_2muex_EV/M #Calculate average after M insertions

def GC_density(z,lam_D,zeta,rho_bulk,sign):
    """SAVE - Gives the GCT density as used for InitialPositiions(...) [and possibly GC_System(...) in the future]
Input:
    z:		z from the wall
    lam_D:      screening length
    zeta:       the applied zeta potential of the system
    sign:       is positive for counter-ions.
Output:
    none
"""    
    psi_GCT=2*np.log((1+np.exp(-(z)/lam_D)*np.tanh(zeta/4.))/(1-np.exp(-(z)/lam_D)*np.tanh(zeta/4.)))
    return rho_bulk*np.exp(sign*psi_GCT)

def simpson(lam_D_spec,zeta,rho_bulk,sign, a, b, N=1000): # N must be even on entry
    "SAVE - Approximate the definite integral of f from a to b, using Composite Simpson's rule - thank you Wikipedia!"
    h = float(b-a)/N   # this is factored out of the loop just for efficiency
    return  h/3 * ( GC_density(a,lam_D_spec,zeta,rho_bulk,sign) + 2*sum((k%2+1) * GC_density(a + h*k,lam_D_spec,zeta,rho_bulk,sign) for k in range(1,N)) + GC_density(b,lam_D_spec,zeta,rho_bulk,sign) )

def WCA(eps,r_cut,r):
    """WCA potential
Input:
    eps:	characteristic energy, units of kT
    r_cut:      intermolecular distance beyond which potential = 0
    		I  report this as \sigma in my results
    r:       	intermolecular separation
Output:
    Evaluation of the WCA potential at r
"""    
    if r>r_cut:
    	return 0
    else:
    	return 4*eps*(0.25*(r_cut/r)**12 - 0.5*(r_cut/r)**6+0.25)
#    	return 4*eps*(0.25*(1/r)**12 - 0.5*(1/r)**6+0.25)

def B2_integrand(eps,r_cut,r):
    return 2*np.pi*(1-np.exp(-WCA(eps,r_cut,r)))*r**2

def B2_WCA(eps,r_cut, a, b, N): # N must be even on entry
    """Approximate the definite integral of f from a to b, using Composite Simpson's rule - thank you Wikipedia!
    http://en.wikipedia.org/wiki/Simpson's_rule#Composite_Simpson.27s_rule"""
    
    h = float(b-a)/N   # this is factored out of the loop just for efficiency
    return h/3 * ( B2_integrand(eps,r_cut,a) + 2*sum((k%2+1) * B2_integrand(eps,r_cut,a + h*k) for k in range(1,N)) + B2_integrand(eps,r_cut,b) )

def Noro_integrand(eps,r_cut,r):
    return (1-np.exp(-WCA(eps,r_cut,r)))

def Noro_WCA(eps,r_cut, a, b, N): # N must be even on entry
    """Approximate the definite integral of f from a to b, using Composite Simpson's rule - thank you Wikipedia!
    http://en.wikipedia.org/wiki/Simpson's_rule#Composite_Simpson.27s_rule"""
    
    h = float(b-a)/N   # this is factored out of the loop just for efficiency
    return h/3 * ( Noro_integrand(eps,r_cut,a) + 2*sum((k%2+1) * Noro_integrand(eps,r_cut,a + h*k) for k in range(1,N)) + Noro_integrand(eps,r_cut,b) )

def P1_LAMMPS_GC_Initialization(Xi_target,Bjerrum,sig_WCA,N_tot,input_files,Max_vol_frac):
    """This will create appropriately named and ready-to-go:
-  input files for use in LAMMPS for the implicit solvent model for P1.
-  associated restart files 
-  required qsub files
    Input: (always dimensional!)
    Xi_target = set point for the coupling parameter. This may be lessened to accomodate the number of specified particles.
    Bjerrum = Bjerrum length 
    sig_WCA = Weeks-Chandler-Andersen particle diameter
    N_tot = set point for the total number of particles placed. This may be increased to ensure there is a bulk, specifically minimum(lamD/L_geo)=1/20
    input_files = number of files to create, keep this lower than 4
    Max_vol_frac = Each input file will run a system with increasing volume fraction and the largest volume fraction equals this parameter
Output:
    None
"""
    import numpy as np
    import random

    pi = np.pi  #This is because I didn't import numpy for the first bit of code I wrote...
    fudge = 1.2 	#Who really knows what this should be?

    Initial_Inputs = '(N_tot,Xi_tar,Bjerrum,sig_WCA) = ('+str(N_tot)+','+str(Xi_target)+','+str(Bjerrum)+','+str(sig_WCA)+')'
    
    Sigma_s = Xi_target / (2.*pi*Bjerrum**2)
    l_geo = Bjerrum * (fudge * N_tot * pi / Xi_target)**0.5 #This x-y (& z) length ensures electroneutrality, 2*Sigma_s*Area = fudge*N_tot, Area = l_geo^2
    dielectric=Bjerrum**-1

	##Trying something new... - 03/02/12 13:15:17 
	#This loop maximized Xi for a total volume fraction less than ~pi/6 (or maximium grid packing of ions)
    while (0.95*l_geo) <= sig_WCA*int(N_tot**(1./3) + 1): #	95% because ions are inserted into a slightly smaller grid than the box - to prevent overlaps
    	Xi_target = Xi_target - 1
	if Xi_target<0:
		print 'Need more ions!',fail
	else:
#		print Xi_target
	    	l_geo = Bjerrum * (fudge * N_tot * pi / Xi_target)**0.5
    Sigma_s = Xi_target / (2.*pi*Bjerrum**2)
    Xi = Xi_target

    #Not 100% sure about the following. It's OK to do, but necessary? Is 20 enough?
    j=0
    while (l_geo)**0.5 >= (1./20)*(16*pi*Bjerrum*N_tot)**0.5: 
    	N_tot = N_tot + 2 #It's always OK to add more particles
    	j+=2
    if j!=0:
    	print '\nHad to add %i ions to ensure sufficient screening.' % j
    
    Phi_totS = np.linspace((N_tot*pi*sig_WCA**3)/(6.*l_geo**3),Max_vol_frac,input_files)
    L_xyS = np.sqrt((N_tot*pi*sig_WCA**3) / (6.*Phi_totS*l_geo))

    if Phi_totS[0]>=Max_vol_frac:
    	print 'Warning! This series of input files does NOT increase in total volume fraction???'


    
    wall_thickness = 1.0 #There is no reason to ever change this.
    L_z = l_geo + 2*wall_thickness
	#This includes all information from the initial, largest system

    if Xi==10000:
    	print '***Special modification to N_tot!!!'
    	N_tot=N_tot*100 #This accounts for overcharging

    print '\nAll initial scales are set!\nN_tot = %i\nXi = %1.2f\nPhi_initial = %1.4e\nBjerrum = %1.2f\nsig_WCA = %1.2f\nl_geo = %1.2f' % (N_tot,Xi,(N_tot*pi*sig_WCA**3)/(6.*l_geo**3),Bjerrum,sig_WCA,l_geo)
    	
    root_name = str(N_tot) + '_' + str(Xi)+'_'+str(Bjerrum)+'_'+str(sig_WCA)+'_'+str(round(l_geo,1))

    k=0
    for L_xy in L_xyS:
    	k+=1
	name = root_name +'_'+str(round(L_xy,1))

    	filename = 'in.N_Xi_lB_WCA_lG_'+name

    	total_steps=0
    	input_file=file(filename,"w")
    	OPT_file=file(filename+'_OPTIMIZATION',"w")
    	
    	input_file.write("""
#################################################
##              P1 Simulations                ##""")
    	input_file.write('\n##\tINPUT: '+Initial_Inputs)
    	input_file.write('\n##\t(N,Xi,lB,WCA,Lz,L_xy) = ('+name+'\n')
    	input_file.write("""#################################################
clear
# Initialize simulation box
dimension	3
boundary	p p f
units		lj
atom_style	charge

# Create geometry""")
    	input_file.write('\nregion  \tsimbox block %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f units box\ncreate_box\t2 simbox\n\n#Insert Particles: Insert ions randomly onto a lattice' % (-0.5*L_xy,0.5*L_xy,-0.5*L_xy,0.5*L_xy,-0.5*L_z,0.5*L_z))
    	OPT_file.write("""
#################################################
##              P1 Simulations                ##""")
    	OPT_file.write('\n##\tINPUT: '+Initial_Inputs)
    	OPT_file.write('\n##\t(N,Xi,lB,WCA,Lz,L_xy) = ('+name+'\n')
    	OPT_file.write("""#################################################
clear
# Initialize simulation box
dimension	3
boundary	p p f
units		lj
atom_style	charge

# Create geometry""")
    	OPT_file.write('\nregion  \tsimbox block %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f units box\ncreate_box\t2 simbox\n\n#Insert Particles: Insert ions randomly onto a lattice' % (-0.5*L_xy,0.5*L_xy,-0.5*L_xy,0.5*L_xy,-0.5*L_z,0.5*L_z))

	#This inserts ions onto a 3D grid with spacing greater than the ion size 
	spacer = l_geo  #"spacer" will be reduced as needed to ensure a minimum number of lattice points
    	NLat_z = l_geo / spacer
   	r_z = np.linspace(-0.5*(L_z-2*wall_thickness)+0.5*sig_WCA,0.5*(L_z-2*wall_thickness)-0.5*sig_WCA,NLat_z)
    	NLat_xy = L_xy / spacer
    	r_xy = np.linspace(-0.5*L_xy+0.5*sig_WCA,0.5*L_xy,NLat_xy)
    	positions=len(r_z)*len(r_xy)**2

	while positions < N_tot:
		spacer=spacer - min(1,sig_WCA)
	    	NLat_z = l_geo / spacer
	   	r_z = np.linspace(-0.5*(L_z-2*wall_thickness)+0.5*sig_WCA,0.5*(L_z-2*wall_thickness)-0.5*sig_WCA,NLat_z)
	    	NLat_xy = L_xy / spacer
	    	r_xy = np.linspace(-0.5*L_xy*0.9+0.5*sig_WCA,0.5*L_xy*0.9,NLat_xy) #The extra 0.9 is there to shrink down the size, again
	    	positions=len(r_z)*len(r_xy)**2
	    	if spacer<=0:
	    		print 'This configuration is too closely packed. Need smaller ions? Disproportionate spacing?',fail

    	all_pos=range(positions)    
    	co_insert=random.sample(all_pos,int(0.5*N_tot))
    	c_insert=random.sample([x for x in set(all_pos).difference(set(co_insert))],int(0.5*N_tot))    
    	i=-1
    	for z in r_z:
		for x in r_xy:
	    		for y in r_xy:
	        		i+=1
	        		if i in co_insert: #Type 1 is co to the LEFT wall
	        			input_file.write("\ncreate_atoms	%i single %1.2f %1.2f %1.2f units box" % (1,x,y,z))
	        			OPT_file.write("\ncreate_atoms		%i single %1.2f %1.2f %1.2f units box" % (1,x,y,z))
				if i in c_insert: #Type 1 is counter to the RIGHT wall
	        			input_file.write("\ncreate_atoms	%i single %1.2f %1.2f %1.2f units box" % (2,x,y,z))	
	        			OPT_file.write("\ncreate_atoms		%i single %1.2f %1.2f %1.2f units box" % (2,x,y,z))	
	
        input_file.write('\n\n##System information:\n# N_tot = %i\n# Xi = %1.2f\n# Phi_tot = %1.4e\n# Bjerrum = %1.2f\n# sig_WCA = %1.2f\n# L_z = %1.2f (includes wall)\n# L_x = L_y = %1.2f\n# Initial Closest Distance = %1.2f' % (N_tot,Xi,(N_tot*pi*sig_WCA**3)/(6.*l_geo*L_xy**2),Bjerrum,sig_WCA,L_z,L_xy,spacer))
        OPT_file.write('\n\n##System information:\n# N_tot = %i\n# Xi = %1.2f\n# Phi_tot = %1.4e\n# Bjerrum = %1.2f\n# sig_WCA = %1.2f\n# L_z = %1.2f (includes wall)\n# L_x = L_y = %1.2f\n# Initial Closest Distance = %1.2f' % (N_tot,Xi,(N_tot*pi*sig_WCA**3)/(6.*l_geo*L_xy**2),Bjerrum,sig_WCA,L_z,L_xy,spacer))
    	input_file.write("""

# Create groups
group		A type 1
group		C type 2
group		ions type 1 2

# Set masses   
mass		1 1.0
mass		2 1.0

# Set charges
set		group A charge -1
set		group C charge 1

# Initilize velocities
velocity	all create 1.0 3""")
    	OPT_file.write("""

# Create groups
group		A type 1
group		C type 2
group		ions type 1 2

# Set masses   
mass		1 1.0
mass		2 1.0

# Set charges
set		group A charge -1
set		group C charge 1

# Initilize velocities
velocity	all create 1.0 3""")

	opt_cut = 'OPT_CUT'
	    	
    	input_file.write("""
    	
# Thermostat & time integration
fix		thermostat all langevin 1.0 1.0 25 3
fix		timeintegration all nve""")

    	input_file.write('\n\ndielectric	%1.4f #=1/Bjerrum' % (dielectric))	
    	input_file.write('\n\nfix		anode all wall/lj93 zlo %1.4f 1. 1.165 1.0 units box' % (-0.5*L_z))
    	input_file.write('\nfix		cathode all wall/lj93 zhi %1.4f 1. 1.165 1.0 units box' % (0.5*L_z))	 
    	input_file.write('\n\nfix		A_field A addforce 0.0 0.0 %1.5f # F = q*E = 1*(Sigma_s/dielectric) = %1.5f/dielectric' % (-Sigma_s/dielectric,-Sigma_s))
    	input_file.write('\nfix		C_field C addforce 0.0 0.0 %1.5f\n' % (Sigma_s/dielectric))	
    	input_file.write("""
timestep 	0.001 #Who knows?
thermo_style	custom step temp etotal pe ecoul evdwl press cpu
thermo		1000

neigh_modify	delay 0 every 1 check yes

# Set potentials""")

    	OPT_file.write("""
    	
# Thermostat & time integration
fix		thermostat all langevin 1.0 1.0 25 3
fix		timeintegration all nve""")

    	OPT_file.write('\n\ndielectric	%1.4f #=1/Bjerrum' % (dielectric))	
    	OPT_file.write('\n\nfix		anode all wall/lj93 zlo %1.4f 1. 1.165 1.0 units box' % (-0.5*L_z))
    	OPT_file.write('\nfix		cathode all wall/lj93 zhi %1.4f 1. 1.165 1.0 units box' % (0.5*L_z))	 
    	OPT_file.write('\n\nfix		A_field A addforce 0.0 0.0 %1.5f # F = q*E = 1*(Sigma_s/dielectric) = %1.5f/dielectric' % (-Sigma_s/dielectric,-Sigma_s))
    	OPT_file.write('\nfix		C_field C addforce 0.0 0.0 %1.5f\n' % (Sigma_s/dielectric))	
    	OPT_file.write("""
timestep 	0.001 #Who knows?
thermo_style	custom step temp etotal pe ecoul evdwl press cpu
thermo		1000

neigh_modify	delay 0 every 1 check yes

# Set potentials""")

		###Potentially useful text below... - 01/18/12 13:30:15 
#    	input_file.write('\n#This is only to prevent particles from overlapping each other.\npair_style	soft %1.4f\npair_coeff	* * 60.0 %1.4f\nrun		35000\n' % (sig_WCA*1.0,sig_WCA*1.0))   	    
	
	
    	input_file.write('\npair_style	lj/cut/coul/long %1.4f %s #The last number in this line must be mannually optimized!' % (sig_WCA,opt_cut))
    	input_file.write('\npair_coeff	* * 1. %1.5f' % (sig_WCA/2.**(1./6.)))
    	input_file.write("""
kspace_style	pppm 1E-4
kspace_modify	slab 3.0
pair_modify     shift yes""")

    	input_file.write('\n\n#dump	        DEBUG all custom 1000 '+name+'_debug.txt id type x y z xu yu zu vx vy vz fx fy fz')
    	input_file.write('\n#run	\t10000\n#write_restart	DEBUG crash\n\nprint "_"\nprint "_"\nprint "Equilibration, (N,Xi,lB,WCA,Lgeo,Lxy) = ('+name+')"\nprint "_"\nprint "_"')
    	input_file.write('\n#dump	        dump0 all custom 1000 '+name+'_equil.txt id type x y z xu yu zu vx vy vz fx fy fz')

    	equilibration=5000000
	production=50000000

    	input_file.write('\nrestart  \t1000 P1_'+name+'_RESa P1_'+name+'_RESb\nrun  \t\t%i' % equilibration)
    	input_file.write('\n#undump\t\tdump0')   
    	input_file.write('\n\nprint "_"\nprint "_"\nprint "Production, (N,Xi,lB,WCA,Lgeo,Lxy) = ('+name+')"\nprint "_"\nprint "_"')
    	input_file.write('\ndump	        dump1 all custom 1000 '+name+'.txt id type x y z xu yu zu vx vy vz fx fy fz')
    	
    	input_file.write('\nrun  \t\t%i' % production)
    	input_file.write('\nwrite_restart   P1_'+name+'_COMPLETE_*')
    	
    	total_steps+=(equilibration+production)
    	input_file.write('\n#################################################')
    	input_file.write('\n# The total number of MD steps for this run file is %i' % total_steps)
    	input_file.close()
    	
	cut_range = np.linspace(0.35*L_xy,L_xy*1.5,45)
	l=-1
	for cut in cut_range:
		l+=1
		OPT_file.write('\npair_style	lj/cut/coul/long %1.4f %1.3f #The last number in this line must be mannually optimized!' % (sig_WCA,cut))
    		if l==0:
			OPT_file.write('\nthermo		100')
    			OPT_file.write('\npair_coeff	* * 1. %1.5f' % (sig_WCA/2.**(1./6.)))
	    		OPT_file.write("""
kspace_style	pppm 1E-4
kspace_modify	slab 3.0
pair_modify     shift yes""")

    		OPT_file.write('\nprint "_"\nprint "_"\nprint "Cut off  = %1.2f"\nprint "_"\nprint "_"' % cut)
	    	OPT_file.write('\nrun  \t\t500\n')
  	OPT_file.close()


	restart_file=file(filename+'_RES',"w")
	restart_file.write('clear\n##INPUT: '+Initial_Inputs+'\nread_restart 	P1_'+name+'RES_aORb')
	restart_file.write("""
group		A type 1
group		C type 2
group		ions type 1 2

fix		thermostat all langevin 1.0 1.0 25 3
fix		timeintegration all nve""")
	restart_file.write('\ndielectric	%1.4f #=1/Bjerrum\n' % (dielectric))
	restart_file.write('\nfix		anode all wall/lj93 zlo %1.4f 1. 1.165 1.0 units box' % (-0.5*L_z))
    	restart_file.write('\nfix		cathode all wall/lj93 zhi %1.4f 1. 1.165 1.0 units box' % (0.5*L_z))
	restart_file.write('\n\nfix		A_field A addforce 0.0 0.0 %1.5f # F = q*E = 1*(Sigma_s/dielectric) = %1.5f/dielectric' % (-Sigma_s/dielectric,-Sigma_s))
    	restart_file.write('\nfix		C_field C addforce 0.0 0.0 %1.5f\n' % (Sigma_s/dielectric))	
	restart_file.write("""
timestep 	0.001 ##Who knows?
thermo_style	custom step temp etotal pe ecoul evdwl press cpu
thermo		1000
""")
	restart_file.write('\npair_style	lj/cut/coul/long %1.4f OPT_CUT #The last number in this line must be mannually optimized!' % sig_WCA)
	restart_file.write('\npair_coeff	* * 1. %1.5f' % (sig_WCA/2.**(1./6.)))
	restart_file.write("""
kspace_style	pppm 1E-4
kspace_modify	slab 3.0
pair_modify     shift yes""")
    	restart_file.write('\n\ndump	        dump1 all custom 1000 '+name+'.gz id type x y z xu yu zu vx vy vz fx fy fz\ndump_modify	dump1 append yes')
    	restart_file.write('\n\nprint "_"\nprint "_"\nprint "Simulation Continued, (N,Xi,lB,WCA,Lgeo,Lxy) = ('+name+')"\nprint "_"\nprint "_"') 

    	restart_file.write('\nrestart  \t1000 P1_'+name+'_RESa P1_'+name+'_RESb\nrun  \t\t%i upto' % total_steps)
    	restart_file.write('\nwrite_restart   P1_'+name+'_COMPLETE_*')
    	restart_file.close()

#	qsub_file.write("""

##$ -o $JOB_NAME.o$JOB_ID
##$ -pe 16way 16 # the last number is (16*#ofNodes), or number of cores?
##$ -q long #"long" for 48 hours, "development" for high priority & 2-hour, or "normal" for 24 hour limit
##$ -l h_rt=48:00:00 # Run time (hh:mm:ss) - 10 minutes (For now, while debugging)
##$ -M bgiera@umail.ucsb.edu
##$ -m be

#module swap pgi intel
#module load mvapich
#module load mkl
#module load fftw2
#module load lammps
#""")	

#	qsubname+=1
	qsubname=name
	qsub_file=file(qsubname+'_qsub.sh',"w")

	qsub_file.write("""
#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
""")  
	qsub_file.write('#$ -N LAMMPS_%s\n\ndate' % name)
   	qsub_file.write('\nhostname\npwd\nls -lh *%s*\n' % name)	

	
	qsub_file.write('TMPDIR=/state/partition1/tmp/$JOB_NAME\nmkdir $TMPDIR')
	qsub_file.write('\ncp -rf $SGE_O_WORKDIR/*%s* $TMPDIR\n' % name)
	qsub_file.write('cp -rf $SGE_O_WORKDIR/P1.py $TMPDIR\n')
	   	
   	#1-proc
#	qsub_file.write('\nLAMMPS_HOME=/share/apps/x86_64/lammps\n')
##	#Uncomment below for initial optimization - 1 processor	 
#	qsub_file.write('\n$LAMMPS_HOME/lmp < %s_OPTIMIZATION > THERMO_%s_OPTIMIZATION.txt \n' % (filename,name))
#   	qsub_file.write('\nhostname\npwd\nls -lh\nsleep 5')
#   	qsub_file.write('\ndate')

#	Uncomment below for prime time - 1 proc

#	qsub_file.write('\n$LAMMPS_HOME/lmp < %s > THERMO_%s.txt\n' % (filename,name))
#	qsub_file.write('\npython P1.py %s.gz\nwait\nmv %s.gz COMPLETED/DATA/.\nwait\nmv THERMO*%s* COMPLETED/THERMO/.\nmv in*%s COMPLETED/LAMMPS_Files/.\nwait\nmv *RES* COMPLETED/RESTARTS/.\n' % (name,name,name,name,name,name,name))
#   	qsub_file.write('\nhostname\npwd\nls -lh\nsleep 5')
#   	qsub_file.write('\ndate')

	numProc = 4
#	#multi-proc
#	print 'There may be OPT qsub problems associated with creating temporary directories...'
#	qsub_file.write('\nLAMMPS_HOME=/home/cask0/home/bgiera/LAMMPS\n#$ -pe mpich %i\n' % numProc)
##	#Uncomment below for initial optimization - multiple processor	 5
#	qsub_file.write('\n#$ -l h_rt=02:00:01\nmpiexec -np %s $LAMMPS_HOME/lmp_ZIN < %s_OPTIMIZATION > $SGE_O_WORKDIR/THERMO_%s_OPTIMIZATION_%s.txt \n' % (str(numProc),filename,name,str(numProc)))
#	qsub_file.write('\nmv *%s* $SGE_O_WORKDIR/.\n' % (name))
#   	qsub_file.write('\nhostname\npwd\nls -lh\nrm -rf $TMPDIR\ndate')

#	print 'new code, check here if there are problems.'
#	#Uncomment below for prime time - multiple processor	 
	qsub_file.write('\nLAMMPS_HOME=/home/cask0/home/bgiera/LAMMPS\n#$ -pe mpich %i\ncd $TMPDIR\n' % numProc)
	qsub_file.write('\nmpiexec -np %s $LAMMPS_HOME/lmp_ZIN < %s > $SGE_O_WORKDIR/THERMO_%s_%s.txt \n' % (str(numProc),filename,name,str(numProc)))
	qsub_file.write('\ngzip %s.txt\nwait\nmv %s.txt.gz %s.gz\nwait\npython P1.py %s.gz\nwait\nmv %s.gz /home/cask0/home/bgiera/P1_Data/COMPLETED/DATA/.\nwait\nmv THERMO*%s* /home/cask0/home/bgiera/P1_Data/COMPLETED/THERMO/.\nwait\nmv in*%s /home/cask0/home/bgiera/P1_Data/COMPLETED/LAMMPS_Files/.\nwait\nmv *RES* /home/cask0/home/bgiera/P1_Data/COMPLETED/RESTARTS/.\n' % (name,name,name,name,name,name,name))
   	qsub_file.write('\nhostname\npwd\nls -lh\nrm -rf $TMPDIR\nwait\ndate')

    	qsub_file.close()
    	print 'qsub %s_qsub.sh\nsleep 2' % str(qsubname)
    return
#IfOnCluster = P1_Analysis()
#P1_CDF()
#P1_Plot()


def P1_LAMMPS_WCA_Initialization(Xi_target,Bjerrum,sig_WCA,N_tot,input_files,Max_vol_frac):
    """This will create appropriately named and ready-to-go:
-  input files for use in LAMMPS for UNCHARGED WCA PARTICLES ONLY
-  associated restart files 
-  required qsub files
    Input: (always dimensional!)
    Xi_target = set point for the coupling parameter. This may be lessened to accomodate the number of specified particles.
    Bjerrum = Bjerrum length 
    sig_WCA = Weeks-Chandler-Andersen particle diameter
    N_tot = set point for the total number of particles placed. This may be increased to ensure there is a bulk, specifically minimum(lamD/L_geo)=1/20
    input_files = number of files to create, keep this lower than 4
    Max_vol_frac = Each input file will run a system with increasing volume fraction and the largest volume fraction equals this parameter
Output:
    None
"""
    import numpy as np
    import random

    pi = np.pi  #This is because I didn't import numpy for the first bit of code I wrote...
    fudge = 1.2 	#Who really knows what this should be?

    Initial_Inputs = '(N_tot,Xi_tar,Bjerrum,sig_WCA) = ('+str(N_tot)+','+str(Xi_target)+','+str(Bjerrum)+','+str(sig_WCA)+')'
    
    Sigma_s = Xi_target / (2.*pi*Bjerrum**2)
    l_geo = Bjerrum * (fudge * N_tot * pi / Xi_target)**0.5 #This x-y (& z) length ensures electroneutrality, 2*Sigma_s*Area = fudge*N_tot, Area = l_geo^2
    dielectric=Bjerrum**-1

	##Trying something new... - 03/02/12 13:15:17 
	#This loop maximized Xi for a total volume fraction less than ~pi/6 (or maximium grid packing of ions)
    while (0.95*l_geo) <= sig_WCA*int(N_tot**(1./3) + 1): #	95% because ions are inserted into a slightly smaller grid than the box - to prevent overlaps
    	Xi_target = Xi_target - 1
	if Xi_target<0:
		print 'Need more ions!',fail
	else:
#		print Xi_target
	    	l_geo = Bjerrum * (fudge * N_tot * pi / Xi_target)**0.5
    Sigma_s = Xi_target / (2.*pi*Bjerrum**2)
    Xi = Xi_target

    #Not 100% sure about the following. It's OK to do, but necessary? Is 20 enough?
    j=0
    while (l_geo)**0.5 >= (1./20)*(16*pi*Bjerrum*N_tot)**0.5: 
    	N_tot = N_tot + 2 #It's always OK to add more particles
    	j+=2
    if j!=0:
    	print '\nHad to add %i ions to ensure sufficient screening.' % j
    
    Phi_totS = np.linspace((N_tot*pi*sig_WCA**3)/(6.*l_geo**3),Max_vol_frac,input_files)
    L_xyS = np.sqrt((N_tot*pi*sig_WCA**3) / (6.*Phi_totS*l_geo))

    if Phi_totS[0]>=Max_vol_frac:
    	print 'Warning! This series of input files does NOT increase in total volume fraction???'


    
    wall_thickness = 1.0 #There is no reason to ever change this.
    L_z = l_geo + 2*wall_thickness
	#This includes all information from the initial, largest system

    if Xi==10000:
    	print '***Special modification to N_tot!!!'
    	N_tot=N_tot*100 #This accounts for overcharging

    print '\nAll initial scales are set!\nN_tot = %i\nXi = %1.2f\nPhi_initial = %1.4e\nBjerrum = %1.2f\nsig_WCA = %1.2f\nl_geo = %1.2f' % (N_tot,Xi,(N_tot*pi*sig_WCA**3)/(6.*l_geo**3),Bjerrum,sig_WCA,l_geo)
    	
    root_name = str(N_tot) + '_' + str(Xi)+'_'+str(Bjerrum)+'_'+str(sig_WCA)+'_'+str(round(l_geo,1))

    k=0
    for L_xy in L_xyS:
    	k+=1
	name = root_name +'_'+str(round(L_xy,1))

    	filename = 'in.N_Xi_lB_WCA_lG_'+name

    	total_steps=0
    	input_file=file(filename,"w")
    	
    	input_file.write("""
#################################################
##              WCA Simulations                ##""")
    	input_file.write('\n##\tINPUT: '+Initial_Inputs)
    	input_file.write('\n##\t(N,Xi,lB,WCA,Lz,L_xy) = ('+name+'\n')
    	input_file.write("""#################################################
clear
# Initialize simulation box
dimension	3
boundary	p p f
units		lj
atom_style	atomic

# Create geometry""")
    	input_file.write('\nregion  \tsimbox block %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f units box\ncreate_box\t2 simbox\n\n#Insert Particles: Insert ions randomly onto a lattice' % (-0.5*L_xy,0.5*L_xy,-0.5*L_xy,0.5*L_xy,-0.5*L_z,0.5*L_z))

	#This inserts ions onto a 3D grid with spacing greater than the ion size 
	spacer = l_geo  #"spacer" will be reduced as needed to ensure a minimum number of lattice points
    	NLat_z = l_geo / spacer
   	r_z = np.linspace(-0.5*(L_z-1.5*wall_thickness)+0.5*sig_WCA,0.5*(L_z-1.5*wall_thickness)-0.5*sig_WCA,NLat_z)
    	NLat_xy = L_xy / spacer
    	r_xy = np.linspace(-0.5*L_xy+0.5*sig_WCA,0.5*L_xy,NLat_xy)
    	positions=len(r_z)*len(r_xy)**2

	while positions < N_tot:
		spacer=spacer - min(1,sig_WCA)
	    	NLat_z = l_geo / spacer
	   	r_z = np.linspace(-0.5*(L_z-1.5*wall_thickness)+0.5*sig_WCA,0.5*(L_z-1.5*wall_thickness)-0.5*sig_WCA,NLat_z)
	    	NLat_xy = L_xy / spacer
	    	r_xy = np.linspace(-0.5*L_xy*0.95+0.5*sig_WCA,0.5*L_xy*0.95,NLat_xy) #The extra 0.9 is there to shrink down the size, again
	    	positions=len(r_z)*len(r_xy)**2
	    	if spacer<=0:
	    		print 'This configuration is too closely packed. Need smaller ions? Disproportionate spacing?',fail

    	all_pos=range(positions)    
    	co_insert=random.sample(all_pos,int(0.5*N_tot))
    	c_insert=random.sample([x for x in set(all_pos).difference(set(co_insert))],int(0.5*N_tot))    
    	i=-1
    	for z in r_z:
		for x in r_xy:
	    		for y in r_xy:
	        		i+=1
	        		if i in co_insert: #Type 1 is co to the LEFT wall
	        			input_file.write("\ncreate_atoms	%i single %1.2f %1.2f %1.2f units box" % (1,x,y,z))
				if i in c_insert: #Type 1 is counter to the RIGHT wall
	        			input_file.write("\ncreate_atoms	%i single %1.2f %1.2f %1.2f units box" % (2,x,y,z))	
	
        input_file.write('\n\n##System information:\n# N_tot = %i\n# Xi = %1.2f\n# Phi_tot = %1.4e\n# Bjerrum = %1.2f\n# sig_WCA = %1.2f\n# L_z = %1.2f (includes wall)\n# L_x = L_y = %1.2f\n# Initial Closest Distance = %1.2f' % (N_tot,Xi,(N_tot*pi*sig_WCA**3)/(6.*l_geo*L_xy**2),Bjerrum,sig_WCA,L_z,L_xy,spacer))
    	input_file.write("""

# Create groups
group		A type 1
group		C type 2
group		ions type 1 2

# Set masses   
mass		1 1.0
mass		2 1.0

# Set charges
##ONLY WCA PARTICLES

# Initilize velocities
velocity	all create 1.0 3""")

	    	
    	input_file.write("""
    	
# Thermostat & time integration
fix		thermostat all langevin 1.0 1.0 25 3
fix		timeintegration all nve""")
	
    	input_file.write('\n\nfix		anode all wall/lj93 zlo %1.4f 1. 1.165 1.0 units box' % (-0.5*L_z))
    	input_file.write('\nfix		cathode all wall/lj93 zhi %1.4f 1. 1.165 1.0 units box' % (0.5*L_z))	 
    	input_file.write("""
timestep 	0.0001 #Who knows?
thermo_style	custom step temp etotal pe ecoul evdwl press cpu
thermo		10000

neigh_modify	delay 0 every 1 check yes

# Set potentials""")

	
    	input_file.write('\npair_style	lj/cut %1.4f ' % (sig_WCA))	
    	input_file.write('\npair_coeff	* * 1. %1.5f' % (sig_WCA/2.**(1./6.)))

    	input_file.write('\n\n#dump	        DEBUG all custom 1000 '+name+'_debug.txt id type x y z xu yu zu vx vy vz fx fy fz')
    	input_file.write('\n#run	\t10000\n#write_restart	DEBUG crash\n\nprint "_"\nprint "_"\nprint "Equilibration, (N,Xi,lB,WCA,Lgeo,Lxy) = ('+name+')"\nprint "_"\nprint "_"')
    	input_file.write('\n#dump	        dump0 all custom 1000 '+name+'_equil.txt id type x y z xu yu zu vx vy vz fx fy fz')

    	equilibration=50000000
	production=500000000

    	input_file.write('\nrestart  \t10000 P1_'+name+'_RESa P1_'+name+'_RESb\nrun  \t\t%i' % equilibration)
    	input_file.write('\n#undump\t\tdump0')   
    	input_file.write('\n\nprint "_"\nprint "_"\nprint "Production, (N,Xi,lB,WCA,Lgeo,Lxy) = ('+name+')"\nprint "_"\nprint "_"')
    	input_file.write('\ndump	        dump1 all custom 10000 '+name+'.gz id type x y z xu yu zu vx vy vz fx fy fz')
    	
    	input_file.write('\nrun  \t\t%i' % production)
    	input_file.write('\nwrite_restart   P1_'+name+'_COMPLETE_*')
    	
    	total_steps+=(equilibration+production)
    	input_file.write('\n#################################################')
    	input_file.write('\n# The total number of MD steps for this run file is %i' % total_steps)
    	input_file.close()



#	qsub_file.write("""

##$ -o $JOB_NAME.o$JOB_ID
##$ -pe 16way 16 # the last number is (16*#ofNodes), or number of cores?
##$ -q long #"long" for 48 hours, "development" for high priority & 2-hour, or "normal" for 24 hour limit
##$ -l h_rt=48:00:00 # Run time (hh:mm:ss) - 10 minutes (For now, while debugging)
##$ -M bgiera@umail.ucsb.edu
##$ -m be

#module swap pgi intel
#module load mvapich
#module load mkl
#module load fftw2
#module load lammps
#""")	

#	qsubname+=1
#	qsubname=name
#	qsub_file=file(qsubname+'_qsub.sh',"w")

#	qsub_file.write("""
##!/bin/bash
##
##$ -V
##$ -cwd
##$ -j y
##$ -S /bin/bash
#""")  
#	qsub_file.write('#$ -N LAMMPS_%s\n\ndate' % name)
#   	qsub_file.write('\nhostname\npwd\nls -lh *%s*\n' % name)	

#	
#	qsub_file.write('TMPDIR=/state/partition1/tmp/$JOB_NAME\nmkdir $TMPDIR')
#	qsub_file.write('\ncp -rf $SGE_O_WORKDIR/*%s* $TMPDIR\n' % name)
#	qsub_file.write('cp -rf $SGE_O_WORKDIR/P1.py $TMPDIR\n')
#	   	
#   	#1-proc
##	qsub_file.write('\nLAMMPS_HOME=/share/apps/x86_64/lammps\n')
###	#Uncomment below for initial optimization - 1 processor	 
##	qsub_file.write('\n$LAMMPS_HOME/lmp < %s_OPTIMIZATION > THERMO_%s_OPTIMIZATION.txt \n' % (filename,name))
##   	qsub_file.write('\nhostname\npwd\nls -lh\nsleep 5')
##   	qsub_file.write('\ndate')

##	Uncomment below for prime time - 1 proc

##	qsub_file.write('\n$LAMMPS_HOME/lmp < %s > THERMO_%s.txt\n' % (filename,name))
##	qsub_file.write('\npython P1.py %s.gz\nwait\nmv %s.gz COMPLETED/DATA/.\nwait\nmv THERMO*%s* COMPLETED/THERMO/.\nmv in*%s COMPLETED/LAMMPS_Files/.\nwait\nmv *RES* COMPLETED/RESTARTS/.\n' % (name,name,name,name,name,name,name))
##   	qsub_file.write('\nhostname\npwd\nls -lh\nsleep 5')
##   	qsub_file.write('\ndate')

#	numProc = 2
##	#multi-proc
##	print 'There may be OPT qsub problems associated with creating temporary directories...'
##	qsub_file.write('\nLAMMPS_HOME=/home/cask0/home/bgiera/LAMMPS\n#$ -pe mpich %i\n' % numProc)
###	#Uncomment below for initial optimization - multiple processor	 5
##	qsub_file.write('\n#$ -l h_rt=02:00:01\nmpiexec -np %s $LAMMPS_HOME/lmp_ZIN < %s_OPTIMIZATION > $SGE_O_WORKDIR/THERMO_%s_OPTIMIZATION_%s.txt \n' % (str(numProc),filename,name,str(numProc)))
##	qsub_file.write('\nmv *%s* $SGE_O_WORKDIR/.\n' % (name))
##   	qsub_file.write('\nhostname\npwd\nls -lh\nrm -rf $TMPDIR\ndate')

##	print 'new code, check here if there are problems.'
##	#Uncomment below for prime time - multiple processor	 
#	qsub_file.write('\nLAMMPS_HOME=/home/cask0/home/bgiera/LAMMPS\n#$ -pe mpich %i\ncd $TMPDIR\n' % numProc)
#	qsub_file.write('\nmpiexec -np %s $LAMMPS_HOME/lmp_ZIN < %s > $SGE_O_WORKDIR/THERMO_%s_%s.txt \n' % (str(numProc),filename,name,str(numProc)))
#	qsub_file.write('\ngzip %s.txt\nwait\nmv %s.txt.gz %s.gz\nwait\npython P1.py %s.gz\nwait\nmv %s.gz /home/cask0/home/bgiera/P1_Data/COMPLETED/DATA/.\nwait\nmv THERMO*%s* /home/cask0/home/bgiera/P1_Data/COMPLETED/THERMO/.\nwait\nmv in*%s /home/cask0/home/bgiera/P1_Data/COMPLETED/LAMMPS_Files/.\nwait\nmv *RES* /home/cask0/home/bgiera/P1_Data/COMPLETED/RESTARTS/.\n' % (name,name,name,name,name,name,name))
#   	qsub_file.write('\nhostname\npwd\nls -lh\nrm -rf $TMPDIR\nwait\ndate')

#    	qsub_file.close()
#    	print 'qsub %s_qsub.sh\nsleep 2' % str(qsubname)
    return
#print 'before function'
#### P1_LAMMPS_WCA_Initialization(Xi_target,Bjerrum,sig_WCA,N_tot,input_files,Max_vol_frac):
#P1_LAMMPS_WCA_Initialization(1.,1.,1,1000,16,0.50)
#print 'now done'







def ExpWall_Sexp():
    import matplotlib.pyplot as plt
    import numpy as np
#    from matplotlib.colors import LinearSegmentedColormap
#    from mpl_toolkits.mplot3d import Axes3D
    import sys

    print 'Plotting along Swall for exponential wall'
    filenameS=[]
    for f in [[x for x in line.split()] for line in file(sys.argv[1],"r").readlines()]:
    	filenameS.append(f[0])

#    filenameS=['MFS_0.05_0.2.txt']
    
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    colS=[]
    i=-1
#    sigS = [0.078,0.236,0.394,1.235,1.325,1.511,0.802,1.704]
    sigS = [0.802,1.235,1.325,1.511,1.704]
    labelX = [r'$\tilde \Sigma = 0.80$',r'$1.24$',r'$1.33$',r'$1.51$',r'$1.70$']
    y=np.linspace(0,0.95)
    ax1.errorbar(4*y,4*y*(y**2+1)/(1-y**2)**2,yerr=None,ls='-',lw=1.3,color='k',marker='None',label=r'$\tilde \rho_{\rm GC}$')  

    rhoatwallS=[]
    for (sig,filename,label_x) in zip(sigS,filenameS,labelX):
    	i+=1
    	print filename

       	MMA=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()])
       	rho = -np.array(MMA[:,0])
       	SGC = MMA[:,1]
       	Sexp = MMA[:,2]

       	col = ROYGBIV_map(sig,2.0,1)
       	colS.append(col)

	rhoatwallS.append(rho[99])
       	ax1.errorbar(Sexp[:99],rho[:99],yerr=None,ls='--',lw=1.0,color=col,marker='None')#,label=label_x)
       	ax1.errorbar(Sexp[99],rho[99],yerr=None,marker = 'o',color=col)#,label=label_x)
       	ax1.errorbar(Sexp[99:],rho[99:],yerr=None,ls='-',lw=1.0,color=col,marker='None',label=label_x)


#    ax1.errorbar(Sexp[49:],rho[49:],yerr=None,ls='--',lw=1.2,color='k',marker='None')#,label=label_x)
#    y=np.linspace(0,0.95)
#    ax1.errorbar(4*y,4*y*(y**2+1)/(1-y**2)**2,yerr=None,ls='-',lw=1.0,color='k',marker='None')  
       	
    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
    ax1.set_xlabel(r'$S_{\rm GC}$',size='small')
    ax1.set_ylabel(r'$-\tilde \rho$',size='small')
    plt.setp(ax1.get_xticklabels(), fontsize='small')
    plt.setp(ax1.get_yticklabels(), fontsize='small')
#    ax1.xaxis.set_major_locator(MaxNLocator(6))
#    ax1.yaxis.set_major_locator(MaxNLocator(5))
    ax1.set_xlim(0,2)     
#    plt.xlim(xmin=0)
#    plt.ylim(ymin=0)
    ax1.set_ylim(0,2)
    fig.set_size_inches(3.37,3.5)
    graphname = 'ExpWall_AlongSGC.png'
    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    plt.close() 

    return rhoatwallS
#rhoatwallS = ExpWall_Sexp()

def ExpWall_contact(rhoatwallS):
    import matplotlib.pyplot as plt
    import numpy as np
#    from matplotlib.colors import LinearSegmentedColormap
#    from mpl_toolkits.mplot3d import Axes3D
    import sys

    print 'Plotting contact for exponential wall'
    filenameS=[]
    for f in [[x for x in line.split()] for line in file(sys.argv[2],"r").readlines()]:
    	filenameS.append(f[0])

#    filenameS=['MFS_0.05_0.2.txt']
    
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    colS=[]
    i=-1

    labelX = [r'$\tilde \Sigma = 0.80$',r'$1.24$',r'$1.33$',r'$1.51$',r'$1.70$']
    x_fit = np.linspace(0,2)
    ax1.errorbar(x_fit,0.5*x_fit*np.sqrt(x_fit**2+4),yerr=None,ls='-',lw=1.3,color='k',marker='None',label=r'$\tilde \rho_{\rm GC}$')
    for (filename,setpoint,label_x) in zip(filenameS,rhoatwallS,labelX):
    	i+=1
    	print filename

       	MMA=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()])
       	Sigma = MMA[:,0]
       	rhocontact = -np.array(MMA[:,1])

#	if i==0:
#		SigMax = Sigma[0]

       	col = ROYGBIV_map(Sigma[len(Sigma)-1],2.0,1)
       	print Sigma[len(Sigma)-1]
       	colS.append(col)


	indexofwall=-1
	for rho in rhocontact:
		indexofwall+=1
		if rho>=setpoint:
			print indexofwall,rho,setpoint
			break

       	ax1.errorbar(Sigma[indexofwall:],rhocontact[indexofwall:],yerr=None,ls='--',lw=1.0,color=col,marker='None')
       	ax1.errorbar(Sigma[indexofwall],rhocontact[indexofwall],yerr=None,marker = 'o',color=col)#,label=label_x)
       	ax1.errorbar(Sigma[:indexofwall],rhocontact[:indexofwall],yerr=None,ls='-',lw=1.0,color=col,marker='None',label=label_x)
		
#       	ax1.errorbar(Sigma,rhocontact,yerr=None,ls='-',lw=1.0,color=col,marker='None')#,label=label_x)
	    
    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
    ax1.set_xlabel(r'$\tilde \Sigma$',size='small')
    ax1.set_ylabel(r'$-\tilde \rho$',size='small')
    plt.setp(ax1.get_xticklabels(), fontsize='small')
    plt.setp(ax1.get_yticklabels(), fontsize='small')
#    ax1.xaxis.set_major_locator(MaxNLocator(6))
#    ax1.yaxis.set_major_locator(MaxNLocator(5))
    ax1.set_xlim(0,2)     
    plt.xlim(xmin=0)
#    plt.ylim(ymin=0)
    ax1.set_ylim(0,2)
    fig.set_size_inches(3.37,3.5)
    graphname = 'ExpWall_Contact.png'
    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    plt.close() 

    return
#ExpWall_contact(rhoatwallS)



#def ExpWall_SGC():
#    import matplotlib.pyplot as plt
#    import numpy as np
##    from matplotlib.colors import LinearSegmentedColormap
##    from mpl_toolkits.mplot3d import Axes3D
#    import sys

#    print 'Plotting along SGC for exponential wall -- THIS IS FUCKED...'
#    filenameS=[]
#    for f in [[x for x in line.split()] for line in file(sys.argv[1],"r").readlines()]:
#    	filenameS.append(f[0])

##    filenameS=['MFS_0.05_0.2.txt']
#    
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
#    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
#    colS=[]
#    i=-1
#    sigS = [0.078,0.236,0.394,1.235,1.325,1.511,0.802,1.704]
#    y=np.linspace(0,0.95)
#    ax1.errorbar(4*y,4*y*(y**2+1)/(1-y**2)**2,yerr=None,ls='-',lw=1.0,color='k',marker='None')    
#    for (sig,filename) in zip(sigS,filenameS):
#    	i+=1
#    	print filename

#       	MMA=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()])
#       	rho = -np.array(MMA[:,0])
#       	SGC = MMA[:,1]
#       	Sexp = MMA[:,2]

#       	col = ROYGBIV_map(sig,2.0,1)
#       	colS.append(col)

#       	ax1.errorbar(SGC[:99],rho[:99],yerr=None,ls='--',lw=1.0,color=col,marker='None')#,label=label_x)
#       	ax1.errorbar(SGC[99],rho[99],yerr=None,marker = 's',color=col)#,label=label_x)
#       	ax1.errorbar(SGC[99:],rho[99:],yerr=None,ls='-',lw=1.0,color=col,marker='None')#,label=label_x)
#       	
##    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=2,handletextpad=0.1) 
#    ax1.set_xlabel(r'$S_{\rm GC}$',size='small')
#    ax1.set_ylabel(r'$-\tilde \rho$',size='small')
#    plt.setp(ax1.get_xticklabels(), fontsize='small')
#    plt.setp(ax1.get_yticklabels(), fontsize='small')
##    ax1.xaxis.set_major_locator(MaxNLocator(6))
##    ax1.yaxis.set_major_locator(MaxNLocator(5))
#    ax1.set_xlim(0,2)     
##    plt.xlim(xmin=0)
##    plt.ylim(ymin=0)
#    ax1.set_ylim(0,2)
#    fig.set_size_inches(3.37,3.5)
#    graphname = 'ExpWall_AlongSGC.png'
#    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    plt.close() 

#    return
##ExpWall_SGC()


def P1_Analysis():
    """
    This function analyzes MD trajectory data from simulations designed according
    	to P1. It will calculate Voltage, density, excluded volume contributions
    	to the excess chemical potential - all needed for P1. All relevant information
    	is printed to *_DATA.txt files to be used for plotting purposes.
    Error analysis has been commented out.
Input:
    The MD trajectory file to be analyzed is entered via the command line. These 
    	trajectory filenames have been designed to contain all necessary simulation info.
Output:
    None, only printed statements that are written to a .txt file from the command line
"""
    #These are hard coded values that cannot be determined from Data[]
    valency=1.0
    valence=1
    r_wall=1.0
    eps_wall=1.0
    eps_WCA=1.0
    SymmetricIons = 'yes'

    
    #These are the necessary libraries for this function.
    	#Check to see if these are on the machines this file is being run on...
    import gzip
    import sys
    import math
    filename = sys.argv[1]
    filename_orignal=filename

    #Extract relevant parameters from filename
    j=0
    value=''
    for letter in filename:
    	if letter!='_':
    		value+=letter
	else:
		j+=1
		if j==1:
			N_tot=int(value)
		if j==2:
			Xi=float(value)
		if j==3:
			Bjerrum=float(value)
		if j==4:
			sig_WCA=float(value)
			r_ion=0.5*sig_WCA
		if j==5:#This length will actually be determined from LAMMPS files
			L_z=float(value) 
		value=''

    if filename[-3:]=='.gz': #This length will actually be determined from LAMMPS files
        L_xy=float(value[:-3])
    else:
        L_xy=float(value[:-4])

    print '\nBeginning data analysis for %s, which has the following parameters:' % filename_orignal
    print 'N_tot = %i\nXi = %1.2f\nBjerrum = %1.2f\nsig_WCA = %1.2f\nPhi_tot = %1.3e' % (N_tot,Xi,Bjerrum,sig_WCA,(N_tot*sig_WCA**3*np.pi)/(6*L_z*L_xy**2))

    Numbins = 500
    bin=np.array([((-L_z/2.)+i*(L_z/float(Numbins))) for i in range(Numbins)]+[L_z*0.5])
    while (bin[1]-bin[0])>sig_WCA:	#This loops makes sure the bins are thinner than the particle diameter
    	Numbins+=1
        bin=np.array([((-L_z/2.)+i*(L_z/float(Numbins))) for i in range(Numbins)]+[L_z*0.5])
    print 'Number of bins = ',len(bin)	

    #Numbins*M_per_slice*number_of_snapshots needs to be a big number! Default is 250E6
    M_per_slice = 10 

    M_per_slice=math.ceil(250E6/(50E3*Numbins)) #This (reasonably) assumes 50E3 snapshots will be observed... - 01/31/12 13:24:45 
    print 'Number of Widom insertions = ',Numbins*M_per_slice*50E3

    #z='3.776'
    #filename='1000_ss_zeta_'+z+'.txt'
    if filename[-3:]=='.gz':
        Data=gzip.GzipFile(fileobj=open(filename,'rb'))
	filename=filename[0:len(filename)-2]+'txt'
    else:
        Data=file(filename,"r")
        
    dielectric=(4*np.pi*Bjerrum)**-1 

    ###I have doubts the "dielectric" should be there. Maybe it's just the 4pi that matters... - 04/07/12 11:31:23 
#    Sigma_s = dielectric*Xi / (2.*np.pi*Bjerrum**2)  #This is the only modification made

	##This is a new thing based upon solid investigation. - 04/07/12 19:30:17 
    print 'This is a new thing!'
    Sigma_s = (1/(4*np.pi))*Xi / (2.*np.pi*Bjerrum**2)  #This is the only modification made

    print filename,Sigma_s
    
    #Create empty vectors to fill
    BA_Data,BA_z=[],[]
    BC_Data,BC_z=[],[]
    Pos_A,Pos_C=[],[]
    #Might as well leave these embeded into the code, commented out... - 01/31/12 11:59:45 
    p_squared,p_avg=0,0 #Am I still going to calculate these?
    			
    wallZ_L,wallZ_R=[],[]
    
    nth=1     # new addition, untested (prob need to add a k+=1 near qz_avg) - 02/10/10 13:17:25 , also never actually used... - 11/18/10 11:47:02    
    n_obs=0
    i,k=0,0 	
    get_box_dims=0
    for y in Data:
      x=y.split()
      ##Format for x is:
#     x= [id type x y z xu yu zu vx vy vz fx fy fz]
      
      get_box_dims+=1
      if get_box_dims<=8:
	if get_box_dims==6:
	  L_x=float(x[1])-float(x[0])
	elif get_box_dims==7:
	  L_y=float(x[1])-float(x[0])
	  area=L_x*L_y
	elif get_box_dims==8:
	  L_box=float(x[1])-float(x[0])
	  z_wall=-0.5*L_box+r_wall
	  
      if not x[0]=="ITEM:":
	x=[float(xx) for xx in x]
	if len(x)<=2:
	  i+=1
	  if i==6: #This means an entire timestep worth of data is had
	    if not (k % nth): #This means data is collected every nth timestep.
		    BA_z=np.array(sorted(BA_z))
	            BC_z=np.array(sorted(BC_z))

		    #Calculate the average z-positions of all ions penetrating the wall
		    wallL,wallR=[],[]
		    for z in BA_z:
		    	if z<(-L_box*0.5+r_wall):
		    		wallL.append(z)
	    		if z>(L_box*0.5-r_wall):
	    			wallR.append(z)
		    for z in BC_z:
		    	if z<(-L_box*0.5+r_wall):
		    		wallL.append(z)
	    		if z>(L_box*0.5-r_wall):
	    			wallR.append(z)		
		    if len(wallL)!=0:
		    	wallZ_L.append(sum(wallL)/len(wallL))
	    	    if len(wallR)!=0:
	    	        wallZ_R.append(sum(wallR)/len(wallR))
	            
	            n_obs+=1
		    if k==0: #This indicates that the first timestep worth of data was collected
		      BA_qty=len(BA_z)
		      BC_qty=len(BC_z)
		      N_tot=BA_qty+BC_qty

      		      bin=np.array([((-L_box/2.)+i*(L_box/float(Numbins))) for i in range(Numbins)]+[L_box*0.5])

		      hist_A = MyHist(BA_z,bin)      
		      
		      Sum_A=hist_A
		      SumSq_A=Sum_A**2
		      
		      hist_C = MyHist(BC_z,bin)   
		      Sum_C=hist_C
		      SumSq_C=Sum_C**2

		      bin_plot=[(x+y)/2. for (x,y) in zip(bin[0:len(bin)-1],bin[1:len(bin)])]
		      L_bin=(bin[1]-bin[0])

		      #Voltage initialization
		      Shell_bin=np.array([x+bin[len(bin)-1] for x in bin])
		      #bin[len(bin)-1] is added to BA_z and BC_z since this voltage calculation assumes the walls go from 0 - L (not -0.5L -> 0.5L as all else)
#		      Psi,Efield=Calc_Volt(Shell_bin,BA_z+bin[len(bin)-1],BC_z+bin[len(bin)-1],Sigma_s,dielectric,area)
			##New
      		      Psi,Efield=quickCalc_VoltandField(Shell_bin,BA_z+bin[len(bin)-1],BC_z+bin[len(bin)-1],Sigma_s,dielectric,area,valence)	
		      Psi_tot=Psi
		      E_tot=Efield
		      Psi_sq=Psi**2
		      E_sq=Efield**2

		      ##Use modified Widom to calculate instantaneous exp_muex_EV = exp(-muex_EV) to add to a total for later averaging
		      		#muex_EV = -np.ln(exp_muex_EV/n_obs)
	    	      if SymmetricIons == 'yes':
			      exp_muex_EV,exp_muex_HS = quick_Calc_muex_EV(bin,Pos_A,Pos_C,M_per_slice,r_ion,eps_WCA,area)
			      exp_muex_EV_tot=exp_muex_EV#+exp_muex_EV_tot
			      exp_muex_EV_sq=exp_muex_EV**2#+exp_muex_EV_sq

			      exp_muex_HS_tot=exp_muex_HS#+exp_muex_HS_tot
			      exp_muex_HS_sq=exp_muex_HS**2#+exp_muex_HS_sq
		      elif SymmetricIons == 'no':
			      exp_1muex_EV,exp_2muex_EV = Calc_2_muex_EV(bin,np.array(Pos_A),np.array(Pos_C),M_per_slice,r_ion,eps_WCA,area)
			      exp_1muex_EV_tot,exp_2muex_EV_tot = exp_1muex_EV,exp_2muex_EV
			      exp_1muex_EV_sq,exp_2muex_EV_sq = exp_1muex_EV**2,exp_2muex_EV**2			      

			##Who knows if this should actually be figured out, but here's the code to start implementing it...
#		      #Correlation time and error analysis, for data sets of 50,000 snapshots!
#		      Psi_5000,Psi_sq_5000=np.zeros(len(Psi),dtype=float),np.zeros(len(Psi),dtype=float)
#		      Psi_7500,Psi_sq_7500=np.zeros(len(Psi),dtype=float),np.zeros(len(Psi),dtype=float)
#		      Psi_10000,Psi_sq_10000=np.zeros(len(Psi),dtype=float),np.zeros(len(Psi),dtype=float)
#		      Psi_22500,Psi_sq_22500=np.zeros(len(Psi),dtype=float),np.zeros(len(Psi),dtype=float)

#		      E_5000,E_sq_5000=np.zeros(len(Efield),dtype=float),np.zeros(len(Efield),dtype=float)
#		      E_7500,E_sq_7500=np.zeros(len(Efield),dtype=float),np.zeros(len(Efield),dtype=float)
#		      E_10000,E_sq_10000=np.zeros(len(Efield),dtype=float),np.zeros(len(Efield),dtype=float)
#		      E_22500,E_sq_22500=np.zeros(len(Efield),dtype=float),np.zeros(len(Efield),dtype=float)

#		      A_count_5000,A_count_sq_5000=np.zeros(len(Sum_A),dtype=float),np.zeros(len(Sum_A),dtype=float)
#		      A_count_7500,A_count_sq_7500=np.zeros(len(Sum_A),dtype=float),np.zeros(len(Sum_A),dtype=float)
#		      A_count_10000,A_count_sq_10000=np.zeros(len(Sum_A),dtype=float),np.zeros(len(Sum_A),dtype=float)
#		      A_count_22500,A_count_sq_22500=np.zeros(len(Sum_A),dtype=float),np.zeros(len(Sum_A),dtype=float)

#		      C_count_5000,C_count_sq_5000=np.zeros(len(Sum_C),dtype=float),np.zeros(len(Sum_C),dtype=float)
#		      C_count_7500,C_count_sq_7500=np.zeros(len(Sum_C),dtype=float),np.zeros(len(Sum_C),dtype=float)
#		      C_count_10000,C_count_sq_10000=np.zeros(len(Sum_C),dtype=float),np.zeros(len(Sum_C),dtype=float)
#		      C_count_22500,C_count_sq_22500=np.zeros(len(Sum_C),dtype=float),np.zeros(len(Sum_C),dtype=float)
		      
		      #Using built in std_dev functions for these
#		      LD_5000,zeta_5000=[],[]
#		      LD_7500,zeta_7500=[],[]
#		      LD_10000,zeta_10000=[],[]
#		      LD_22500,zeta_22500=[],[]
#		      LD_tot=0

#		      LD=(8*np.pi*((np.mean(hist_A[Numbins/2-2:Numbins/2+2])+np.mean(hist_C[Numbins/2-2:Numbins/2+2]))/(2.*area*L_bin))*Bjerrum)**-0.5
#		      LD_5000.append(LD)
#		      LD_tot+=LD

		      z_min=np.min([np.min(BA_z),np.min(BC_z)])
		      z_max=np.max([np.max(BA_z),np.max(BC_z)])
		    #At some point error analysis is needed on voltage and density profiles
		      k+=1
		    else:
		      hist_A = MyHist(BA_z,bin)
		      Sum_A=Sum_A+hist_A
		      SumSq_A=SumSq_A+hist_A**2

		      hist_C = MyHist(BC_z,bin)
		      Sum_C=Sum_C+hist_C
		      SumSq_C=SumSq_C+hist_C**2      
		  
#		      Psi,Efield=Calc_Volt(Shell_bin,BA_z+bin[len(bin)-1],BC_z+bin[len(bin)-1],Sigma_s,dielectric,area)
			#New
      		      Psi,Efield=quickCalc_VoltandField(Shell_bin,BA_z+bin[len(bin)-1],BC_z+bin[len(bin)-1],Sigma_s,dielectric,area,valence)	
		      Psi_tot=Psi_tot+Psi
		      E_tot=E_tot+Efield		      
		      Psi_sq=Psi_sq+Psi**2
		      E_sq=E_sq+Efield**2  
 
		      ##Use modified Widom to calculate instantaneous exp_muex_EV = exp(-muex_EV) to add to a total for later averaging
		      		#muex_EV = -np.ln(exp_muex_EV/n_obs)
      	    	      if SymmetricIons == 'yes':
			      exp_muex_EV,exp_muex_HS = quick_Calc_muex_EV(bin,Pos_A,Pos_C,M_per_slice,r_ion,eps_WCA,area)
			      exp_muex_EV_tot=exp_muex_EV+exp_muex_EV_tot
			      exp_muex_EV_sq=exp_muex_EV**2+exp_muex_EV_sq

			      exp_muex_HS_tot=exp_muex_HS+exp_muex_HS_tot
			      exp_muex_HS_sq=exp_muex_HS**2+exp_muex_HS_sq			      
		      elif SymmetricIons == 'no':
			      exp_1muex_EV,exp_2muex_EV = Calc_2_muex_EV(bin,np.array(Pos_A),np.array(Pos_C),M_per_slice,r_ion,eps_WCA,area)
			      exp_1muex_EV_tot+=exp_1muex_EV
			      exp_2muex_EV_tot+=exp_2muex_EV
			      exp_1muex_EV_sq = exp_1muex_EV_sq + exp_1muex_EV**2
			      exp_2muex_EV_sq = exp_2muex_EV_sq + exp_2muex_EV**2	
			      
		      z_min=np.min([np.min(BA_z),np.min(BC_z),z_min])
		      z_max=np.max([np.max(BA_z),np.max(BC_z),z_max])

  ######	Capacitance information to go here!!! - 11/05/10 12:14:29
  			#This is actually worthless, but perhaps a useful future thing to calculate... 
		    z_plus=np.mean(BC_z)
		    z_minus=np.mean(BA_z)
		    qz_avg=(z_plus-z_minus)/2.

		    #These quantities still need to be averaged after the loop is through
		    p_avg+=(N_tot*qz_avg)
		    p_squared+=(N_tot*qz_avg)**2

		    #The following is a progress meter
		    if not n_obs % 1000:
		    	print filename,n_obs,n_obs/float(50000)
		    
#		    #Correlation time and error analysis, these still must be averaged after the loop is through! 
#		    if n_obs<=5000:
#			Psi_5000+=Psi
#			Psi_sq_5000+=Psi**2
#			E_5000+=Efield
#			E_sq_5000+=Efield**2
#			A_count_5000+=hist_A
#			A_count_sq_5000+=hist_A**2
#			C_count_5000+=hist_C
#			C_count_sq_5000+=hist_C**2
#			LD_5000.append((8*np.pi*((np.mean(hist_A[Numbins/2-2:Numbins/2+2])+np.mean(hist_C[Numbins/2-2:Numbins/2+2]))/(2.*area*L_bin))*Bjerrum)**-0.5)
#			zeta_5000.append(np.mean(abs(Psi[1]),abs(Psi[-1])))    
#		    elif n_obs>5000 and n_obs<=12500:
#			Psi_7500+=Psi
#			Psi_sq_7500+=Psi**2
#			E_7500+=Efield
#			E_sq_7500+=Efield**2
#			A_count_7500+=hist_A
#			A_count_sq_7500+=hist_A**2
#			C_count_7500+=hist_C
#			C_count_sq_7500+=hist_C**2
#			LD_7500.append((8*np.pi*((np.mean(hist_A[Numbins/2-2:Numbins/2+2])+np.mean(hist_C[Numbins/2-2:Numbins/2+2]))/(2.*area*L_bin))*Bjerrum)**-0.5)
#			zeta_7500.append(np.mean(abs(Psi[1]),abs(Psi[-1])))    
#		    elif n_obs>12500 and n_obs<=22500:
#			Psi_10000+=Psi
#			Psi_sq_10000+=Psi**2
#			E_10000+=Efield
#			E_sq_10000+=Efield**2
#			A_count_10000+=hist_A
#			A_count_sq_10000+=hist_A**2
#			C_count_10000+=hist_C
#			C_count_sq_10000+=hist_C**2
#			LD_10000.append((8*np.pi*((np.mean(hist_A[Numbins/2-2:Numbins/2+2])+np.mean(hist_C[Numbins/2-2:Numbins/2+2]))/(2.*area*L_bin))*Bjerrum)**-0.5)
#			zeta_10000.append(np.mean(abs(Psi[1]),abs(Psi[-1]))) 
#		    elif n_obs>22500:
#			Psi_22500+=Psi
#			Psi_sq_22500+=Psi**2
#			E_22500+=Efield
#			E_sq_22500+=Efield**2
#			A_count_22500+=hist_A
#			A_count_sq_22500+=hist_A**2
#			C_count_22500+=hist_C
#			C_count_sq_22500+=hist_C**2
#			LD_22500.append((8*np.pi*((np.mean(hist_A[Numbins/2-2:Numbins/2+2])+np.mean(hist_C[Numbins/2-2:Numbins/2+2]))/(2.*area*L_bin))*Bjerrum)**-0.5)
#			zeta_22500.append(np.mean(abs(Psi[1]),abs(Psi[-1])))  
			
		    BA_z=[]
		    BC_z=[]
		    Pos_A=[]
		    Pos_C=[]
		    i=1
	elif x[1]==1:
	  BA_z.append(x[4])
	  Pos_A.append([x[2],x[3],x[4]])
	elif x[1]==2:
	  BC_z.append(x[4])
  	  Pos_C.append([x[2],x[3],x[4]])

    #And to get the last (tricky) snapshot... - 4/27/11 @ 3.30p
    n_obs+=1
    
    hist_A = MyHist(BA_z,bin)
    Sum_A=Sum_A+hist_A
    SumSq_A=SumSq_A+hist_A**2

    hist_C = MyHist(BC_z,bin)
    Sum_C=Sum_C+hist_C
    SumSq_C=SumSq_C+hist_C**2     
    
#    Psi,Efield=Calc_Volt(Shell_bin,BA_z+bin[len(bin)-1],BC_z+bin[len(bin)-1],Sigma_s,dielectric,area)
	##New
    Psi,Efield=quickCalc_VoltandField(Shell_bin,BA_z+bin[len(bin)-1],BC_z+bin[len(bin)-1],Sigma_s,dielectric,area,valence)	 
    Psi_tot=Psi_tot+Psi
    E_tot=E_tot+Efield		      
    Psi_sq=Psi_sq+Psi**2
    E_sq=E_sq+Efield**2    

    z_plus=np.mean(BC_z)
    z_minus=np.mean(BA_z)
    qz_avg=(z_plus-z_minus)/2.

    p_avg+=(N_tot*qz_avg)
    p_squared+=(N_tot*qz_avg)**2

#    #This must be _LARGEST_BLOCK_SIZE
#    Psi_22500+=Psi
#    Psi_sq_22500+=Psi**2
#    E_22500+=Efield
#    E_sq_22500+=Efield**2
#    A_count_22500+=hist_A
#    A_count_sq_22500+=hist_A**2
#    C_count_22500+=hist_C
#    C_count_sq_22500+=hist_C**2
#    LD_22500.append((8*np.pi*((np.mean(hist_A[Numbins/2-2:Numbins/2+2])+np.mean(hist_C[Numbins/2-2:Numbins/2+2]))/(2.*area*L_bin))*Bjerrum)**-0.5)
#    zeta_22500.append(np.mean(abs(Psi[1]),abs(Psi[-1])))  
   
    print "Sorted data..."

    L_z=L_box-2*r_wall
    pot_drop=Sigma_s*L_z/dielectric
    dipole_avg=p_avg/n_obs
    dipole_avg_sq=p_squared/n_obs

    kBT=1.
    beta=1./kBT
    dipole_var=beta*(dipole_avg_sq-dipole_avg**2.)

	##As is, this is fairly wortheless information.
    #Per conversation with Shell, to match up with theory, remove the 1 in (area*dielectric/L_z)*(1-dipole_avg/(dielectric*pot_drop))**-1
    C_tot=(area*dielectric/L_z)*(1-dipole_avg/(dielectric*pot_drop))**-1
    #Per conversation with Shell, to match up with theory, remove the 1 in (area*dielectric/L_z)*(1-dipole_var/(dielectric*area*L_z))**-1
    C_dif=(area*dielectric/L_z)*(1-dipole_var/(dielectric*area*L_z))**-1

    A_count=Sum_A/float(n_obs)
    C_count=Sum_C/float(n_obs)

    BA_den=Sum_A/(float(n_obs)*area*L_bin) #Unnormalized
    BC_den=Sum_C/(float(n_obs)*area*L_bin) #Unnormalized
      
    n0_BA=np.mean(BA_den[Numbins/2-2:Numbins/2+2])
    n0_BC=np.mean(BC_den[Numbins/2-2:Numbins/2+2])
       
    n0=(n0_BA+n0_BC)/2. #I think this is plenty fine as long as cations and anions differ only in +/- charge
#    print "Bulk concentration is %1.8f" % n0

    lambda_D=(8.*np.pi*n0*Bjerrum)**(-0.5) # added valency on 02/18/10 14:20:23 
    print "Lambda_D is %1.5f" % lambda_D
    if (L_z/lambda_D)<20:
    	print '\nWarning: The electrode double layers may be overlapping!!!\n'
       
#    Psi_var=Psi_sq/n_obs - (Psi_tot/n_obs)**2
#    E_var=E_sq/n_obs - (E_tot/n_obs)**2
#    A_var=SumSq_A/n_obs - (Sum_A/n_obs)**2
#    C_var=SumSq_C/n_obs - (Sum_C/n_obs)**2
#    LD_var,zeta_var=[],[]
#    for (L,z) in zip(LD_5000,zeta_5000):
#    	LD_var.append(L)
#    	zeta_var.append(z)  
#    for (L,z) in zip(LD_7500,zeta_7500):
#    	LD_var.append(L)
#    	zeta_var.append(z)  
#    for (L,z) in zip(LD_10000,zeta_10000):
#    	LD_var.append(L)
#    	zeta_var.append(z)  
#    for (L,z) in zip(LD_22500,zeta_22500):
#    	LD_var.append(L)
#    	zeta_var.append(z)  

#    print 'Other lam_D = %1.5f' % np.mean(LD_var)	
#    LD_var=np.var(LD_var)
#    zeta_var=np.var(zeta_var)
    
#    Psi_5000_var=Psi_sq_5000/5000 - (Psi_5000/5000)**2
#    E_5000_var=E_sq_5000/5000 - (E_5000/5000)**2
#    A_5000_var=A_count_sq_5000/5000 - (A_count/5000)**2
#    C_5000_var=C_count_sq_5000/5000 - (C_count/5000)**2
#    LD_5000_var=np.var(LD_5000)
#    zeta_5000_var=np.var(zeta_5000)

#    Psi_7500_var=Psi_sq_7500/7500 - (Psi_7500/7500)**2
#    E_7500_var=E_sq_7500/7500 - (E_7500/7500)**2
#    A_7500_var=A_count_sq_7500/7500 - (A_count/7500)**2
#    C_7500_var=C_count_sq_7500/7500 - (C_count/7500)**2
#    LD_7500_var=np.var(LD_7500)
#    zeta_7500_var=np.var(zeta_7500)
#    
#    Psi_10000_var=Psi_sq_10000/10000 - (Psi_10000/10000)**2
#    E_10000_var=E_sq_10000/10000 - (E_10000/10000)**2
#    A_10000_var=A_count_sq_10000/10000 - (A_count/10000)**2
#    C_10000_var=C_count_sq_10000/10000 - (C_count/10000)**2
#    LD_10000_var=np.var(LD_10000)
#    zeta_10000_var=np.var(zeta_10000)

#    Psi_22500_var=Psi_sq_22500/22500 - (Psi_22500/22500)**2
#    E_22500_var=E_sq_22500/22500 - (E_22500/22500)**2
#    A_22500_var=A_count_sq_22500/22500 - (A_count/22500)**2
#    C_22500_var=C_count_sq_22500/22500 - (C_count/22500)**2
#    LD_22500_var=np.var(LD_22500)
#    zeta_22500_var=np.var(zeta_22500)

#    print '\nPsi error ' + filename
#    Psi_err=[]
#    o=-1
#    for (y_5000,y_7500,y_10000,y_22500,constant) in zip(Psi_5000_var,Psi_7500_var,Psi_10000_var,Psi_22500_var,Psi_var):
#	o+=1
#    	y=[y_5000,y_7500,y_10000,y_22500]
#    	x=[constant/5000,constant/7500,constant/10000,constant/22500]
#	two_tau,intercept=np.polyfit(x,y,1)
#	Psi_err.append(constant*2/(n_obs*two_tau*0.5))
#	print bin[o]
#	for (yp,xp) in zip(y,x):
#	  print xp,yp
#	print '\n'  
#	
#    E_err=[]
#    for (y_5000,y_7500,y_10000,y_22500,constant) in zip(E_5000_var,E_7500_var,E_10000_var,E_22500_var,E_var):
#    	y=[y_5000,y_7500,y_10000,y_22500]
#    	x=[constant/5000,constant/7500,constant/10000,constant/22500]
#	two_tau,intercept=np.polyfit(x,y,1)
#	E_err.append(constant*2/(n_obs*two_tau*0.5))

#    A_err=[]
#    for (y_5000,y_7500,y_10000,y_22500,constant) in zip(A_5000_var,A_7500_var,A_10000_var,A_22500_var,A_var):
#    	y=[y_5000,y_7500,y_10000,y_22500]
#    	x=[constant/5000,constant/7500,constant/10000,constant/22500]
#	two_tau,intercept=np.polyfit(x,y,1)
#	A_err.append(constant*2/(n_obs*two_tau*0.5))

#    C_err=[]
#    for (y_5000,y_7500,y_10000,y_22500,constant) in zip(C_5000_var,C_7500_var,C_10000_var,C_22500_var,C_var):
#    	y=[y_5000,y_7500,y_10000,y_22500]
#    	x=[constant/5000,constant/7500,constant/10000,constant/22500]
#	two_tau,intercept=np.polyfit(x,y,1)
#	C_err.append(constant*2/(n_obs*two_tau*0.5))

#    y=[LD_5000_var,LD_7500_var,LD_10000_var,LD_22500_var]
#    x=[LD_var/5000,LD_var/7500,LD_var/10000,LD_var/22500]
#    two_tau,intercept=np.polyfit(x,y,1)
#    LD_err=LD_var*2/(n_obs*two_tau*0.5)

#    y=[zeta_5000_var,zeta_7500_var,zeta_10000_var,zeta_22500_var]
#    x=[zeta_var/5000,zeta_var/7500,zeta_var/10000,zeta_var/22500]
#    two_tau,intercept=np.polyfit(x,y,1)
#    zeta_err=zeta_var*2/(n_obs*two_tau*0.5)

    Psi_tot=Psi_tot/float(n_obs)  
    E_tot=E_tot/float(n_obs)
    Psi_tot=Psi_tot-np.mean(Psi_tot)
    print 'Voltage is is calculated the traditional way'

    AvInWallZ = r_wall - 0.5*(abs(np.mean(wallZ_L)+L_z*0.5)+(np.mean(wallZ_R)-L_z*0.5))

    if SymmetricIons == 'yes':
	    exp_muex_EV,exp_muex_HS = quick_Calc_muex_EV(bin,Pos_A,Pos_C,M_per_slice,r_ion,eps_WCA,area)
	    exp_muex_EV_tot=exp_muex_EV+exp_muex_EV_tot
	    exp_muex_EV_sq=exp_muex_EV**2+exp_muex_EV_sq
	    exp_muex_HS_tot=exp_muex_HS+exp_muex_HS_tot
	    exp_muex_HS_sq=exp_muex_HS**2+exp_muex_HS_sq
    
	    muex_EV=-np.log(exp_muex_EV_tot/float(n_obs))
	    muex_EV_total=-np.log(sum(exp_muex_EV_tot)/float(n_obs*len(exp_muex_EV_tot)))

	    muex_HS=-np.log(exp_muex_HS_tot/float(n_obs))
	    muex_HS_total=-np.log(sum(exp_muex_HS_tot)/float(n_obs*len(exp_muex_HS_tot)))


#	    total_prefix='Analyzed_error' + filename[4:len(filename)-4]
#	    total_output=file(total_prefix+"_DATA.txt","w")
#	    for (z,V,E,Np,Nm,V_err,Ef_err,A_err,C_err,Ovar_Psi,Ovar_E,Ovar_A,Ovar_C) in zip(Shell_bin.tolist()+[Sigma_s],Psi_tot.tolist()+[dipole_avg],E_tot.tolist()+[dipole_avg_sq],A_count.tolist()+[area,Bjerrum],C_count.tolist()+[lambda_D,r_ion],Psi_err+[LD_err],E_err+[zeta_err],A_err+[LD_var,zeta_var],C_err+[np.pi,np.pi],Psi_var.tolist()+[np.pi],E_var.tolist()+[np.pi],A_var.tolist()+[np.pi,np.pi],C_var.tolist()+[np.pi,np.pi]):
#	    	total_output.write("%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\n" % (z,V,E,Np,Nm,V_err,Ef_err,A_err,C_err,Ovar_Psi,Ovar_E,Ovar_A,Ovar_C))   
#	    total_output.close()

#	    total_prefix='Analyzed_' + filename[4:len(filename)-4]
#	    total_output=file(total_prefix+"_DATA.txt","w")
#	    for (z,V,E,Np,Nm) in zip(Shell_bin.tolist()+[Sigma_s],Psi_tot.tolist()+[dipole_avg],E_tot.tolist()+[dipole_avg_sq],A_count.tolist()+[area,Bjerrum],C_count.tolist()+[lambda_D,r_ion]):
#	    	total_output.write("%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\n" % (z,V,E,Np,Nm))   
#	    if (len(wallZ_L)!=0) and (len(wallZ_R)!=0):
#		total_output.write("%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\n" % (0,0,0,np.mean(wallZ_L),np.mean(wallZ_R))) 
#	    total_output.close()
#2**3/(24*0.96009**2)



	    total_prefix='Analyzed_P1_' + filename[:len(filename)-4]
	    total_output=file(total_prefix+".txt","w")
	    for (z,V,E,Np,Nm,muexEV,muexHS) in zip(Shell_bin.tolist()+[Sigma_s],Psi_tot.tolist()+[z_min],E_tot.tolist()+[z_max],A_count.tolist()+[area,Bjerrum],C_count.tolist()+[lambda_D,r_ion],muex_EV.tolist()+[AvInWallZ],muex_HS.tolist()+[AvInWallZ]):
	    	total_output.write("%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\n" % (z,V,E,Np,Nm,muexEV,muexHS))   
	    total_output.close()    
    elif SymmetricIons == 'no':
	    muex_EV1=-np.log(exp_1muex_EV_tot/float(n_obs))
	    muex_EV2=-np.log(exp_2muex_EV_tot/float(n_obs))

# 	    total_prefix='Analyzed_sigr' + filename[4:len(filename)-4]
#	    total_output=file(total_prefix+"_DATA.txt","w")
#	    for (z,V,E,Np,Nm) in zip(Shell_bin.tolist()+[Sigma_s],Psi_tot.tolist()+[dipole_avg],E_tot.tolist()+[dipole_avg_sq],A_count.tolist()+[area,Bjerrum],C_count.tolist()+[lambda_D,r_ion[0]]):
#	    	total_output.write("%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\n" % (z,V,E,Np,Nm))   
#	    if (len(wallZ_L)!=0) and (len(wallZ_R)!=0):
#		total_output.write("%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\n" % (0,0,0,np.mean(wallZ_L),np.mean(wallZ_R))) 
#	    total_output.close()
		      	
	    total_prefix='Analyzed_P1_sigr_' + str(M_per_slice) + filename[4:len(filename)-4]
	    total_output=file(total_prefix+".txt","w")
	    for (z,V,E,Np,Nm,muexEV1,muexEV2) in zip(Shell_bin.tolist()+[Sigma_s],Psi_tot.tolist()+[dipole_avg],E_tot.tolist()+[dipole_avg_sq],A_count.tolist()+[area,Bjerrum],C_count.tolist()+[lambda_D,r_ion[1]],muex_EV1.tolist()+[np.mean(wallZ_L)],muex_EV2.tolist()+[np.mean(wallZ_R)]):
	    	total_output.write("%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\t\t\t\t%-1.5f\n" % (z,V,E,Np,Nm,muexEV1,muexEV2))   
	    total_output.close()

    print 'Finished ' + filename_orignal
    return 'OnlyAnalyze'
P1_Analysis()
#P1_Plot()






#if IfOnCluster!='OnlyAnalyze':
######	P1_LAMMPS_GC_Initialization(Xi_target,Bjerrum,sig_WCA,N_tot,input_files,Max_vol_frac)
###	#
#	##Successfully run code below. These are in zin:~/P1_Data/COMPLETED
###	P1_LAMMPS_GC_Initialization(5,1,1,1000,32,0.20) #21 successfully running (L_xy<=11.5 did not run) - 02/07/12 10:31:50
###	P1_LAMMPS_GC_Initialization(7,1,1,1000,32,0.20) #All submitted (but didn't immediately start running) underway - 02/09/12 18:34:49 
###	P1_LAMMPS_GC_Initialization(5,1,2,1000,32,0.50) # - 02/07/12 10:31:50
###	P1_LAMMPS_GC_Initialization(7,1,2,1000,32,0.50) #L_xy<21.7 didn't run - 02/07/12 10:31:50

######	Third Fleet, attempted on Mellon
###	P1_LAMMPS_GC_Initialization(10,1,1,1000,4,0.05) #Lowest Phi_bulk Submitted on Mellon - 02/09/12 18:34:49 
###	P1_LAMMPS_GC_Initialization(100,10,1,1000,4,0.1) #Lowest Phi_bulk Submitted on Mellon- 02/09/12 18:34:49 
###	P1_LAMMPS_GC_Initialization(10000,100,1,1000,4,0.1) #Lowest Phi_bulk Submitted on Mellon- 02/13/12 12:22:38 
##		###The following Wigner attempt does N_tot=5*N_tot right before LAMMPS creation!
###	P1_LAMMPS_GC_Initialization(10000,100,0.5,1000/5,4,0.1) #Lowest Phi_bulk Submitted on Mellon -~02/23/12 15:30:46 
##	P1_LAMMPS_GC_Initialization(10000,100,0.5,1000/10,1,0.1) #Lowest Phi_bulk Submitted on Mellon -~03/02/12 11:10:36 

##	P1_LAMMPS_GC_Initialization(10000,100,0.5,1000/100,1,0.1) #Lowest Phi_bulk Submitted on Mellon -~3/21/12
#	P1_LAMMPS_GC_Initialization(10000,100,0.5,5000/100,1,0.1) #Lowest Phi_bulk Submitted on Mellon -~3/21/12


#	P1_LAMMPS_GC_Initialization(0.1,0.1,1.0,1000,32,0.50) ## - 03/06/12 13:35:42 
#	P1_LAMMPS_GC_Initialization(0.1,0.1,0.5,1000,32,0.50) ## - 03/09/12 11:33:28 	
#	P1_LAMMPS_GC_Initialization(0.5,0.1,0.5,1500,32,0.50) ## - 03/06/12 13:35:42

#	a=1
#	P1_LAMMPS_GC_Initialization(0.1,0.1,0.5,1000,10,0.10) ## - ~05/02/12 15:27:49 
#	P1_LAMMPS_GC_Initialization(0.25,0.1,0.5,1000,10,0.10) ## - ~05/02/12 15:27:49 
#	P1_LAMMPS_GC_Initialization(0.4,0.1,0.5,1000,10,0.10) ## - ~05/02/12 15:27:49 
	
##IfOnCluster = P1_Analysis()re()
#    ax1=fig.add_subplot(111)
#    fig.subplots_adjust(right=0.98) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.17)##Larger adds whitespace
#    fig.subplots_adjust(left=0.26) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.95) ##Smaller adds whitespace to top

#	##Horizontal cconfiguration
##    fig.subplots_adjust(right=0.98) ##Lower puts more whitespace on the right
##    fig.subplots_adjust(bottom=0.27)##Larger adds whitespace
##    fig.subplots_adjust(left=0.14) #Larger puts more whitespace on the left
##    fig.subplots_adjust(top=0.95) ##Smaller adds whitespace to top

#    x_fit = np.linspace(0,13,100)
#    ax1.plot(x_fit,0.5*x_fit*np.sqrt(x_fit**2+4),color='k',lw=1.5,ls='-')  
#    for (x,y,c) in zip(SigmaCorrectedS,ContactS,colors):
#    	ax1.errorbar(x,y,yerr=None,marker='^',ms=5.5,color=c)
#    ax1.set_ylim(0,50)
#    ax1.set_xlim(0,9.5) 
#    xlabel=r'$\Sigma_{app} \times 4 \pi \lambda_{B} \lambda_{D}/q_{\pm} $'
#    ylabel=r'$-\rho_{f}(0) \times 8 \pi \lambda_{B} \lambda_{D}^{2}/q_{\pm}$'
#    ax1.set_xlabel(xlabel,size='small')
#    ax1.set_ylabel(ylabel,size='small')
#    
#    plt.setp(ax1.get_xticklabels(), fontsize='small')
#    plt.setp(ax1.get_yticklabels(), fontsize='small')
#	##Horizontal
##    fig.set_size_inches(1.55*2,0.65*2)
#	##Vertical
#    fig.set_size_inches(0.92*2,1.25*2)

##    plt.show()
#    plt.savefig('Contact_vs_App.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 
#    plt.close()
#'
#		print 'exportname = \"~/Desktop/MMA_fitWCA_'+filename[11:]+'\"'
#		print 'bulk = ',(sigWCA*0.954)**3/(24*Bjerrum*lam_D**2)
#		print 'effectivefield = ',-eff
#		print '\n\n'
