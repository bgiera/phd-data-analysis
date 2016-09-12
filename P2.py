# -*- coding: utf-8 -*-
#!/usr/bin/python
import copy
import numpy as np
import gzip
import EDLCalcs


def CSInterpolation_ntot(MMAz,MMAnp,MMAnm,MMAfree,MMAvolt,MMASig,zeta_target):
    """This takes a working solution from a larger CS zeta file to 
    	calculate the CS result for a smaller zeta file..
It's general enough to be modified to extract other things, like interpolated free charge
Input:
    MMAz:	CS distance data
    MMAvolt:    CS voltage data
    MMA_phi:    CS volume fraction data
    zeta_target: target zeta for new theory curve, zeta_target<=MMAvolt[0]
Output:
	this can put out anything, I suppose
    
""" 
    from scipy import interpolate
    print 'WARNING: this does not appear to work for ntot quantities...'
     
    t1 = np.linspace(MMAz[0],MMAz[-1],20)[1:-1]  #Check out "20" if this is giving issues
    CS_volt = interpolate.LSQUnivariateSpline(MMAz,MMAvolt,t1,k=5)
    fitCS_volt = lambda z: CS_volt(z)[0]

    dif = 1
    i=-1
    finer_z = np.linspace(MMAz[0],MMAz[-1],len(MMAz)*20)
    while dif>=0.001:
    	i+=1
#    	print i,finer_z[i],fitCS_volt(finer_z[i]),dif
    	dif = fitCS_volt(finer_z[i])-zeta_target
    	
    zset = finer_z[i]
    
    t1 = np.linspace(MMAz[0],MMAz[-1],20)[1:-1]  #Check out "20" if this is giving issues
    ntot = np.array(MMAnp) + np.array(MMAnm)
    CS_tot = interpolate.LSQUnivariateSpline(MMAz,ntot,t1,k=5)
    fitCS_tot = lambda z: CS_tot(z)[0]
    
    interpCS_z = np.linspace(zset,MMAz[-1],len(MMAz))
    interpCS_tot = [fitCS_tot(z) for z in interpCS_z]

    return interpCS_z-zset,interpCS_tot

def ROYGBIV_map(value,max_value,sat=1):
    import colorsys
    return colorsys.hsv_to_rgb(value / max_value , sat, sat)
#    return colorsys.hsv_to_rgb(value / max_value / (1.1), sat, sat)

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
    charge=1.
    for (i,z) in enumerate(z_positions):
	for z_charge in z_plus[z_plus<=z]:
	    Psi_tot[i]=Psi_tot[i]-charge*(z-z_charge)/(epsilon*A_xy)
	    E_tot[i]=E_tot[i]-charge/(epsilon*A_xy)
    charge=-1.
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

#		###The next two lines of code ensure calculations are made at the
#		#######	same value of z (not random z's). This is for debugging.
#		#08/25/11 16:21:38  
#		zmid = 0.5*(zlo+zhi)
#		zlo,zhi=zmid,zmid
	    	
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

		##Erase the following commented out sections when you have some time... - 03/18/13 15:59:35 
#def Calc_muex_EV(z_positions,Pos_A,Pos_C,M,r_ion,eps_WCA,A_xy):
#	    """
#	Input:
#	    z_positions:   Bin discretation
#	    Pos_A:	   array of (x,y,z) positions for ions (attracted to right plate)
#	    Pos_C:	   array of (x,y,z) positions for ions (attracted to left plate)
#	    M:		   the number of random insertions per bin
#	    r_ion:	   WCA radius for ions
#	    eps_WCA:	   Energy scale for WCA potential
#	    A_xy:	   Cross-sectional area of simulation box
#	Output:
#	    exp_muex_EV:	   Excluded volume contributions to mu_ex 
#	    			NOTE: muex_EV[z] = -ln(exp_muex_EV[z]/n_obs) must be taken before printing 
#	"""
#	    import random
#	    L_xy = A_xy**0.5
#	    hL_xy=L_xy*0.5
#	    iL_xy=1./L_xy
#	    sig_WCA=2*r_ion
#	    sig_WCA_sq=sig_WCA**2


#	    print 'erase this after done testing'
#	   
#	  #Initialize array of exp(muex_EV), which we later average 
#	  #	and subesequently take the natural log of to get muex_EV(z)
#	    exp_muex_EV=np.zeros(len(z_positions), dtype=float)

#	    i=-1 #This index is used to fill muex_EV[i]
#	    for (zlo,zhi) in zip(z_positions[0:-1],z_positions[1:len(z_positions)]):
##	    	print zlo,zhi

#		i+=1
#	    	m=0 #This index and while loop to ensure that 
#	    	
#		#Trim away ions whose z-coordinates are +/- sig_WCA from bin edges	
#	    	Pos_near_A = Pos_A[np.logical_and(Pos_A[:,2]>(zlo-sig_WCA),Pos_A[:,2]<(zhi+sig_WCA))]
#		Pos_near_C = Pos_C[np.logical_and(Pos_C[:,2]>(zlo-sig_WCA),Pos_C[:,2]<(zhi+sig_WCA))]
#		#This can't be done apriori for (x,y), since these are periodic dimensions


#######		print 'erase this after done testing'
#		Pos_near_A = Pos_A
#		Pos_near_C = Pos_C
#		
#	
#		if (len(Pos_near_A)+len(Pos_near_C))!=0: #Perform M random 
#				#insertions if simulation particles have
#				#nearby z-coordinated.
#			while m < M:
#				m+=1
#				#Reset energies
#				del_Energies_A,del_Energies_C=[],[] 
#				
#				#Randomly choose (x,y,z) of test ion with a z
#				#	position inside the bin of interest
#				r0=np.array([random.uniform(-hL_xy,hL_xy),random.uniform(-hL_xy,hL_xy),random.uniform(zlo,zhi)])

#				r0 = np.array([-8.5, -5.9, -248.2])

#				#Calculate the distance between the test ion and Anions
#				if len(Pos_near_A)!=0: #provided anions are near
#					r_0j=r0-Pos_near_A
#					r_0j[:,:2]=r_0j[:,:2]-L_xy*np.round_(r_0j[:,:2]*iL_xy)
#					r_0j_sq=np.sum(r_0j*r_0j,axis=1)

#					sig_over_r_p6=(sig_WCA_sq/r_0j_sq[r_0j_sq<sig_WCA_sq])**3
#					#Calculate energy due to inserting test ion
#					del_Energies_A = 4.*eps_WCA*(sig_over_r_p6*(sig_over_r_p6-1)+0.25)
#				else:
#					del_Energies_A=[]
#				
#				if len(Pos_near_C)!=0:
#					#Do the same thing for Cations...
#					r_0j=r0-Pos_near_C
#					r_0j[:,:2]=r_0j[:,:2]-L_xy*np.round_(r_0j[:,:2]*iL_xy)
#					r_0j_sq=np.sum(r_0j*r_0j,axis=1)
#					sig_over_r_p6=(sig_WCA_sq/r_0j_sq[r_0j_sq<sig_WCA_sq])**3
#					del_Energies_C = 4.*eps_WCA*(sig_over_r_p6*(sig_over_r_p6-1)+0.25)
#				else:
#					del_Energies_C=[]
#					
#		    		exp_muex_EV[i]+=np.exp(-(sum(del_Energies_A)+sum(del_Energies_C)))
#		else:
#			del_Energies_A,del_Energies_C=[],[]
#			exp_muex_EV[i]+=M*np.exp(-(sum(del_Energies_A)+sum(del_Energies_C)))

#	    return exp_muex_EV/M #Calculate average after M insertions

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

def CountsandVelocity(zVz,bin):
    """Returns a histogram of the z-positions for ions
Input:
    zVz:  (N,many) sorted array of atomic positions
    bin:  vector of bin dimensions
Output:
	counts:	a histogram of counts for z-values within their respective bins
"""
    #Start with zeroed bins
    counts= np.zeros(len(bin)-1, int)
    vz= np.zeros(len(bin)-1, int)
    k_store=0
  
    for x in zVz:  
        k,get_next_x=k_store,False
        while (k<(len(bin)-1) and not(get_next_x)): #iterate through all the bin sets
	  if x[0]>=bin[k] and x[0]<bin[k+1]: #check to see if the z-position is within a particular bin
	    vz[k]+=x[1]
	    counts[k]+=1
	    k=k_store
	    get_next_x=True
	  k+=1 #iterare the bin index
	  
    mean_vz=[]
    for (count,z) in zip(counts,vz):
    	if count!=0:
    	     count=float(count)   
    	     mean_vz.append(z/count)
    	elif count==0:
    	     mean_vz.append(0.)
       
    return mean_vz,counts

def TransientAnalysis(filename,Bjerrum,r_ion):#,fnam):
    """As is, this function must be double checked for functionality.
	It should be able to calculate Q_EDL=Q_EDL(time).
	This will serve as the template to calculate other spatially and
		temporally varying quanties: voltage(t,z), concentration(t,z), etc. 

Input
    filename: str() the name of the file with MD trajectory data
    Bjerrum:  float() the Bjerrum length that was set for the system
    r_ion:    float() the WCA radius of the ion OR
    		list() of WCA radii r_ion =[smaller,larger]
    M_per_slice: number of insertions per bin for muex_EV calculations
    ##CONSIDER ADDING VALENCY AS A FUCNTION INPUT!
Output:
    None, only printed statements that are written to a .txt file from the command line
"""
    from operator import itemgetter
    #These are hard coded values that cannot be determined from Data[]
    valence=1
    r_wall=1.0
    eps_wall=1.0
    eps_WCA=1.0
    Numbins = 502
#    print 'Adjusting bin size to ',Numbins

    print 'This code calculates excess chemical potential of hard sphere, but does not yet store running averages nor print this information to analysis output file!'

    M_per_slice = 10

    if type(r_ion) != type(0.5):
    	print 'Output code must be adjusted before analyzing this data!',fail
    	

    filename_orignal=filename

    if filename[-3:]=='.gz':
        Data=gzip.GzipFile(fileobj=open(filename,'rb'))
	filename=filename[0:len(filename)-2]+'txt'
    else:
        Data=file(filename,"r")


    if '1000' in filename:
    	vz_index = 10
    else:
	print 'Assuming efficient data output'
    	vz_index = 5    	

    print '\n'
    print filename,1

    

    dielectric=(4*np.pi*Bjerrum)**-1 #This used to be dielectric=Bjerrum**-1

    ze=str(filename[-9:-4])
#    print 'The goal as of 11/05/10 16:08:46 should be to feed Sigma_s via the LAMMPS filename.'
	
    if	Bjerrum==0.1:
		if ze=='0.000':
			Sigma_s=0.0*dielectric
		elif ze=='0.500':
			Sigma_s=0.03368*dielectric
		elif ze=='1.000':
			Sigma_s=0.06948*dielectric
		elif ze=='1.500':
			Sigma_s=0.10964*dielectric
		elif ze=='1.750':
			Sigma_s=0.13213*dielectric
		elif ze=='2.000':
			Sigma_s=0.15669*dielectric
		elif ze=='2.500':
			Sigma_s=0.21359*dielectric
		elif ze=='2.750':
			Sigma_s=0.24682*dielectric
		elif ze=='3.000':
			Sigma_s=0.28390*dielectric
		elif ze=='3.500':
			Sigma_s=0.37206*dielectric						
		elif ze=='4.000':
			Sigma_s=0.48358*dielectric
		elif ze=='5.000':
			Sigma_s=0.80669*dielectric
		elif ze=='6.000':
			Sigma_s=1.33572*dielectric		
		elif ze=='7.000':
			Sigma_s=2.20568*dielectric
		else:
			Sigma_s=0.0000				
    else:
	Sigma_s=0.0000
	print 'Sigma_s was set to a default value of 0.000'
	#This code is as new as 11/05/10 16:07:32 

    nth=1     # new addition, untested (prob need to add a k+=1 near qz_avg) - 02/10/10 13:17:25 , also never actually used... - 11/18/10 11:47:02    
    n_obs=0
    i=0 
    get_time=0	
    get_box_dims=0
    time_old=-1000
    time_new=time_old
    charges=1
    
    #Create empty vectors to fill
    BA_Data,BA_z=[],[]
    BC_Data,BC_z=[],[]
    Pos_A,Pos_C=[],[]
    Q_DL=[]
    Psi_tr=[]
    Efield_tr=[]
    Sum_A_tr=[]
    Sum_C_tr=[]
    exp_muex_EV_tr=[]
    exp_1muex_EV_tr=[]
    exp_2muex_EV_tr=[]
    A_vz,C_vz=[],[]
    A_vz_tr,C_vz_tr=[],[]
    Q_L=[]
    Q_R=[]
    
        
    t_index=0

    jj=0

    for y in Data:    
      x=y.split()
      
      ##This gets the current time_step
      if get_time==1:
        time_old=time_new
	time_new=int(x[0])
      	get_time=0	
      if len(x)==2:
      	if x[1]=='TIMESTEP':
     		get_time=1

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
    	  bin=np.array([((-L_box/2.)+qpp*(L_box/float(Numbins))) for qpp in range(Numbins)]+[L_box*0.5])
	  Shell_bin=np.array([value+bin[len(bin)-1] for value in bin])
	  L_bin = Shell_bin[1]-Shell_bin[0]
	  z_wall=-0.5*L_box+r_wall

      if not x[0]=="ITEM:":
	x=[float(xx) for xx in x]
	if len(x)<=2:
	  i+=1
	  if i==6: #This means an entire timestep worth of data is had
	    n_obs+=1

	    BA_z=np.array(BA_z)
	    BC_z=np.array(BC_z)

	    A_vz=sorted(np.array(A_vz), key=itemgetter(0))
	    C_vz=sorted(C_vz, key=itemgetter(0))
	    NA_L= len(BA_z[BA_z<0.])
	    NA_R= len(BA_z[BA_z>=0.])
	    NC_L= len(BC_z[BC_z<0.])
	    NC_R= len(BC_z[BC_z>=0.])
	    mean_qdl = ((NA_L-NC_L) + (NC_R-NA_R))*0.5
	    q_r = (NC_R-NA_R)
	    q_l = (NA_L-NC_L)

	    if charges==1:  #This indicates the first charge cycle data has been gathered.
	    	Q_DL.append(mean_qdl)
	    	Q_L.append(q_l)
	    	Q_R.append(q_r)

	    	##Counts (for densities) initilization 
		A_vz,hist_A = CountsandVelocity(A_vz,bin)
#		hist_A = MyHist(BA_z,bin)       
		Sum_A_tr.append(np.array(hist_A))
		A_vz_tr.append(A_vz)

		C_vz,hist_C = CountsandVelocity(C_vz,bin)
#		hist_C = MyHist(BC_z,bin)   
		Sum_C_tr.append(np.array(hist_C))
		C_vz_tr.append(C_vz)

		#Voltage & Efield initialization
		##Old, slow code
#		Psi,Efield=Calc_Volt(Shell_bin,BA_z+bin[len(bin)-1],BC_z+bin[len(bin)-1],Sigma_s,dielectric,area)


		Psi,Efield=quickCalc_VoltandField(Shell_bin,BA_z+bin[len(bin)-1],BC_z+bin[len(bin)-1],Sigma_s,dielectric,area,valence)	
		Psi_tr.append(np.array(Psi))
		Efield_tr.append(np.array(Efield))
	        N_tot = sum(hist_A)+sum(hist_C)
    	        if type(r_ion) == type(0.5):
#			      exp_muex_EV= Calc_muex_EV(bin,np.array(Pos_A),np.array(Pos_C),M_per_slice,r_ion,eps_WCA,area)
			      exp_muex_EV,exp_muex_HS = quick_Calc_muex_EV(bin,Pos_A,Pos_C,M_per_slice,r_ion,eps_WCA,area)
			      exp_muex_EV_tr.append(np.array(exp_muex_EV))

		elif len(r_ion) == 2:
			      print "Warning: These calculations are much slower as they are not coded in Fortran"
			      exp_1muex_EV,exp_2muex_EV = Calc_2_muex_EV(bin,np.array(Pos_A),np.array(Pos_C),M_per_slice,r_ion,eps_WCA,area)
			      exp_1muex_EV_tr.append(np.array(exp_1muex_EV))
			      exp_2muex_EV_tr.append(np.array(exp_2muex_EV))	    	
	    else:

	    	Q_DL[t_index]+=mean_qdl 
	    	Q_L[t_index]+=q_l
	    	Q_R[t_index]+=q_r

		A_vz,hist_A = CountsandVelocity(A_vz,bin)
#		hist_A = MyHist(BA_z,bin)  
		Sum_A_tr[t_index]+=np.array(hist_A)
		A_vz_tr[t_index]+=np.array(A_vz)

		C_vz,hist_C = CountsandVelocity(C_vz,bin)
#		hist_C = MyHist(BC_z,bin)   
		Sum_C_tr[t_index]+=np.array(hist_C)
		C_vz_tr[t_index]+=np.array(C_vz)

		Psi,Efield=quickCalc_VoltandField(Shell_bin,BA_z+bin[len(bin)-1],BC_z+bin[len(bin)-1],Sigma_s,dielectric,area,valence)
		Psi_tr[t_index]+=np.array(Psi)
		Efield_tr[t_index]+=np.array(Efield)

    	        if type(r_ion) == type(0.5):
#			      exp_muex_EV= Calc_muex_EV(bin,np.array(Pos_A),np.array(Pos_C),M_per_slice,r_ion,eps_WCA,area)
			      exp_muex_EV,exp_muex_HS = quick_Calc_muex_EV(bin,Pos_A,Pos_C,M_per_slice,r_ion,eps_WCA,area)
			      exp_muex_EV_tr[t_index]+=np.array(exp_muex_EV)
		elif len(r_ion) == 2: ##This is for when ions have different sizes
			      exp_1muex_EV,exp_2muex_EV = Calc_2_muex_EV(bin,np.array(Pos_A),np.array(Pos_C),M_per_slice,r_ion,eps_WCA,area)
			      exp_1muex_EV_tr[t_index]+=np.array(exp_1muex_EV)
			      exp_2muex_EV_tr[t_index]+=np.array(exp_2muex_EV)	


		##This is old code that works only for non-parallelized Transient codes
#	    if time_old>time_new:
#		charges+=1
#		print filename,charges,t_index
#		t_index=0 	
#            else:
#		t_index+=1	

		##Added 03/14/13 10:24:52 
		##This bit of code is robust and can handle any type of transient data input
	    if charges>1:
		if t_index<charge_time:
			t_index+=1
		else:
			charges+=1
			print filename,charges,t_index
			t_index=0 
			
	    else:	
		    if time_old>time_new:
			charge_time = t_index
			charges+=1
			print filename,charges,t_index
			t_index=0 	
		    else:
			t_index+=1
	    BA_z=[]
	    BC_z=[]
	    Pos_A=[]
	    Pos_C=[]
	    A_vz=[]
	    C_vz=[]
	    i=1
	elif x[1]==1:
	  BA_z.append(x[4])
	  Pos_A.append([x[2],x[3],x[4]])
	  A_vz.append([x[4],x[vz_index]])
	elif x[1]==2:
	  BC_z.append(x[4])
  	  Pos_C.append([x[2],x[3],x[4]])
  	  C_vz.append([x[4],x[vz_index]])

    #And to get the last (tricky) snapshot... - 4/27/11 @ 3.30p
    n_obs+=1
    BA_z=np.array(BA_z)
    BC_z=np.array(BC_z)
    N_tot = len(BA_z)+len(BC_z)
    NA_L= len(BA_z[BA_z<0.])
    NA_R= len(BA_z[BA_z>=0.])
    NC_L= len(BC_z[BC_z<0.])
    NC_R= len(BC_z[BC_z>=0.])
    mean_qdl = ((NA_L-NC_L) + (NC_R-NC_R))*0.5
    q_r = (NC_R-NA_R)
    q_l = (NA_L-NC_L)

		#UNcomment all of this stuff, it's still good
    if charges!=1:
	    Q_DL[t_index]+=mean_qdl
	    Q_L[t_index]+=q_l
	    Q_R[t_index]+=q_r

	    A_vz,hist_A = CountsandVelocity(A_vz,bin)
	#    hist_A = MyHist(BA_z,bin)  
	    Sum_A_tr[t_index]+=np.array(hist_A)
	    A_vz_tr[t_index]+=np.array(A_vz)

	    C_vz,hist_C = CountsandVelocity(C_vz,bin)
	#	hist_C = MyHist(BC_z,bin)   
	    Sum_C_tr[t_index]+=np.array(hist_C)
	    C_vz_tr[t_index]+=np.array(C_vz)

	    Psi,Efield=quickCalc_VoltandField(Shell_bin,BA_z+bin[len(bin)-1],BC_z+bin[len(bin)-1],Sigma_s,dielectric,area,valence)
	    Psi_tr[t_index]+=np.array(Psi)
	    Efield_tr[t_index]+=np.array(Efield)

	    if type(r_ion) == type(0.5):
	#	      exp_muex_EV= Calc_muex_EV(bin,np.array(Pos_A),np.array(Pos_C),M_per_slice,r_ion,eps_WCA,area)
		      exp_muex_EV,exp_muex_HS = quick_Calc_muex_EV(bin,Pos_A,Pos_C,M_per_slice,r_ion,eps_WCA,area)
		      exp_muex_EV_tr[t_index]+=np.array(exp_muex_EV)
	    elif len(r_ion) == 2:
		      exp_1muex_EV,exp_2muex_EV = Calc_2_muex_EV(bin,np.array(Pos_A),np.array(Pos_C),M_per_slice,r_ion,eps_WCA,area)
	 	      exp_1muex_EV_tr[t_index]+=np.array(exp_1muex_EV)
	 	      exp_2muex_EV_tr[t_index]+=np.array(exp_2muex_EV)	
    	
###    #Average the results (negative sign just to flip curves over x-axis):

    Q_DL=-np.array(Q_DL)/float(charges)
    Q_L=-np.array(Q_L)/float(charges) # (negative sign just to flip curves over x-axis)
    Q_R=-np.array(Q_R)/float(charges) # (negative sign just to flip curves over x-axis)	
    Sum_A_tr = np.array(Sum_A_tr) / float(charges)
    Sum_C_tr = np.array(Sum_C_tr) / float(charges)
    A_vz_tr = np.array(A_vz_tr ) / float(charges)
    C_vz_tr = np.array(C_vz_tr ) / float(charges)    
    Psi_tr = np.array(Psi_tr) / float(charges)
    Efield_tr = np.array(Efield_tr) / float(charges)
    
    if type(r_ion) == type(0.5):
  	    muex_EV_tr=-np.log(np.array(exp_muex_EV_tr)/float(charges))
    elif len(r_ion) == 2:
  	    muex_EV1_tr=-np.log(np.array(exp_1muex_EV_tr)/float(charges))
  	    muex_EV2_tr=-np.log(np.array(exp_2muex_EV_tr)/float(charges))
  
    print "Writing analyzed data to file..."
    total_prefix='Transient_' + filename[4:len(filename)-4]
    total_output=file(total_prefix+"_DATA.txt","w")
    for q in Q_DL:
    	total_output.write("%-1.5f\n" % (q))    
    total_output.close()

		##Convert htis to proper outpit
    Xi = Sigma_s*2.*np.pi*Bjerrum**2
#    newname = 'Analyzed_tV_'+str(N_tot)+'_'+str(Xi)+'_'+str(Bjerrum)+'_'+str(2*r_ion)+'_'+str(L_box)+'_'+str(round(np.sqrt(area),2))

    newname = 'Analy'+str(Numbins)+'_Dr_'+str(N_tot)+'_'+str(Xi)+'_'+str(Bjerrum)+'_'+str(2*r_ion)+'_'+str(L_box)+'_'+str(round(np.sqrt(area),2))

#    newname = str(fnam)+'Analy'+str(Numbins)+'_tV_'+str(N_tot)+'_'+str(Xi)+'_'+str(Bjerrum)+'_'+str(2*r_ion)+'_'+str(L_box)+'_'+str(round(np.sqrt(area),2))
    
    total_output= file(newname+".txt","w") 

    for tcount in np.arange(len(Q_DL)):
    	A_count = Sum_A_tr[tcount]
    	C_count = Sum_C_tr[tcount]

    	vzA = A_vz_tr[tcount]
    	vzC = C_vz_tr[tcount]

    	ql = Q_L[tcount]
    	qr = Q_R[tcount]
    	
    	A_den = np.mean(A_count[Numbins/2-2:Numbins/2+2])/(area*L_bin)
    	C_den = np.mean(C_count[Numbins/2-2:Numbins/2+2])/(area*L_bin)
    	n0 = (A_den + C_den)*0.5
    	lambda_D=(8.*np.pi*n0*Bjerrum)**(-0.5)
    	    	    	
    	for (z,V,E,Np,Nm,muexEV,vA,vC) in zip(Shell_bin,Psi_tr[tcount],Efield_tr[tcount],A_count.tolist()+[tcount],C_count.tolist()+[lambda_D],muex_EV_tr[tcount],vzA.tolist()+[ql],vzC.tolist()+[qr]):
		total_output.write("%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\n" % (z,V,E,Nm,Np,muexEV,vA,vC))
    total_output.close()

    print 'Finished ' + filename_orignal
    return
####    def TransientAnalysis(filename,Bjerrum,r_ion):
temp_rion =0.15*10
#temp_rion = 4.8000*0.5
print 'This should be programmed to work on diffusivity ratios'
fnam = 50

#TransientAnalysis('100C_1000_tr_zeta_3.500.txt',0.1,0.15*10)

#TransientAnalysis(str(fnam)+'_1000_tr_zeta_3.500.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_3.500.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_3.500.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_3.500.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_3.500.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_3.500.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_3.500.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_3.500.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_3.500.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_3.500.txt',0.1,temp_rion,fnam)
#fnam+=1


#TransientAnalysis(str(fnam)+'_1000_tr_zeta_2.000.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_2.000.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_2.000.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_2.000.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_2.000.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_2.000.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_2.000.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_2.000.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_2.000.txt',0.1,temp_rion,fnam)
#fnam+=1
#TransientAnalysis(str(fnam)+'_1000_tr_zeta_2.000.txt',0.1,temp_rion,fnam)
#fnam+=1






#TransientAnalysis('1000_tr_zeta_4.000.gz',0.1,temp_rion)
#TransientAnalysis('1000_tr_zeta_3.500.gz',0.1,temp_rion)
#TransientAnalysis('1000_tr_zeta_3.000.gz',0.1,temp_rion)

#TransientAnalysis('1000_tr_zeta_2.750.gz',0.1,temp_rion)
#TransientAnalysis('1000_tr_zeta_2.500.gz',0.1,temp_rion)
#TransientAnalysis('1000_tr_zeta_2.000.gz',0.1,temp_rion)

#TransientAnalysis('1000_tr_zeta_1.750.gz',0.1,temp_rion)
#TransientAnalysis('1000_tr_zeta_1.500.gz',0.1,temp_rion)
#TransientAnalysis('1000_tr_zeta_1.000.gz',0.1,temp_rion)
#TransientAnalysis('1000_tr_zeta_0.500.gz',0.1,temp_rion)


#TransientAnalysis('500_tr_zeta_3.500.gz',0.1,temp_rion)
#TransientAnalysis('500_tr_zeta_1.750.gz',0.1,temp_rion)

#print 'here'
#TransientAnalysis('1000_tr_zeta_3.776.txt',0.1,temp_rion)
#print 'now here'

#TransientAnalysis('1000_tr_zeta_3.500.txt',0.1,temp_rion)


def SSDiffusion(filename,Bjerrum):
    """
    This function finds an the diffusion constant.
    This code could easily be modified to find multiple diffusion constants in order to get an
    error on the value
Input:
    Data:     (N,3) array of atomic positions
    lambda_D: Dimless Debye-Length
Output:
    None, only printed statements that are written to a .txt file from the command line
"""
    #This is a one-time deal function. It wasn't designed to be used for every data set, but could
    # be modfied to allow for this kind of functionality
    from operator import itemgetter

    LAMMPS_dt=0.001 #This must be manually entered!
    r_wall=1.0
    Numbins = 502


    
	##THis is the old code
#    print 'This function is working...'
#    if filename[-3:]=='.gz':
#        Data=[[float(x) for x in line.split()] for line in gzip.GzipFile(fileobj=open(filename,'rb')).readlines() if not line.startswith("ITEM")]
#	filename=filename[0:len(filename)-2]+'txt'
#    else:
#        Data=[[float(x) for x in line.split()] for line in file(filename,"r").readlines() if not line.startswith("ITEM")]
#    
#    print '\n'
#    print filename
#    print "Loaded data..."
    
	##This is the code that must be run::


    #z='3.776'
    #filename='1000_ss_zeta_'+z+'.txt'
    filename_orignal=filename

    if filename[-3:]=='.gz':
        Data=gzip.GzipFile(fileobj=open(filename,'rb'))
	filename=filename[0:len(filename)-2]+'txt'
    else:
        Data=file(filename,"r")
    print '\n'
    print filename
    
#    #"Data" for this function must be a system that is equilibrated 1MM+ runs or more

#    L=np.sum(abs(Data[4][0]-Data[4][1]))-2*r_wall
#                  
#    #Various iterators
#    i,j,k,l=0,0,0,0 #l can go
#    a,b,c,d=0,0,-1,-1 #not sure why these are here

#    #Create empty vectors to fill
#    BA_xyz,BA_0,MSD_BA=[],[],[]
#    BC_xyz,BC_0,MSD_BC=[],[],[]
#    
#    del_time=int(np.sum(Data[int(np.sum(Data[1]))+5])-np.sum(Data[0]))      
#    
#    for x in Data:
#      if len(x)!=14:
#	i+=1
#	#Sorting keeps atoms' XYZ data in order
#	BA_xyz=sorted(BA_xyz)
#	BC_xyz=sorted(BC_xyz)	
#	if i==6: #This means an entire timestep worth of data is had
#	  if k==0: #This indicates that it's the first timestep worth of data having been collected, Pos0's need to be generated
#	    l+=1 #Not sure what this l's doing, 99.9% sure it can be deleted   
#	    BA_xyz0=copy.deepcopy(BA_xyz)
#	    BC_xyz0=copy.deepcopy(BC_xyz)
#	    k+=1 #This stops XX_0 positions to be copied again
#	    
##	  print np.mean([((xyz[1][0]-xyz0[1][0])**2+(xyz[1][1]-xyz0[1][1])**2+(xyz[1][2]-xyz0[1][2])**2) for (xyz,xyz0) in zip(BA_xyz,BA_xyz0)])
#	  MSD_BA.append(np.mean([((xyz[1][0]-xyz0[1][0])**2+(xyz[1][1]-xyz0[1][1])**2+(xyz[1][2]-xyz0[1][2])**2) for (xyz,xyz0) in zip(BA_xyz,BA_xyz0)]))
#	  MSD_BC.append(np.mean([((xyz[1][0]-xyz0[1][0])**2+(xyz[1][1]-xyz0[1][1])**2+(xyz[1][2]-xyz0[1][2])**2) for (xyz,xyz0) in zip(BC_xyz,BC_xyz0)]))

#	  BA_xyz=[]
#	  BC_xyz=[]

#	  i=1
#      #XYZ data is stored for every timestep for each particle type
#      elif x[1]==1:  #this is only for the two type model
#	point=[x[0],x[5:8]]
#	BA_xyz.append(point)
#      elif x[1]==2: #this is only for the two type model
#	point=[x[0],x[5:8]]
#	BC_xyz.append(point)


    BA_xyz,BC_xyz=[],[]
    MSDs,MSDAs,MSDCs=[],[],[]
    BA_z,BC_z=[],[]
#    Sum_A_tr=[]
#    Sum_C_tr=[]
    
    nth=1     # new addition, untested (prob need to add a k+=1 near qz_avg) - 02/10/10 13:17:25 , also never actually used... - 11/18/10 11:47:02    
    n_obs=0
    i,k=0,0 	
    get_time=0	
    time_old=-1000
    time_new=time_old
    get_box_dims=0
    for y in Data:
      x=y.split()
      ##Format for x is:
#     x= [id type x y z xu yu zu vx vy vz fx fy fz]

      ##This gets the current time_step
      if get_time==1:
        time_old=time_new
	time_new=int(x[0])
      	get_time=0	
      if len(x)==2:
      	if x[1]=='TIMESTEP':
     		get_time=1
      
      get_box_dims+=1
      if get_box_dims<=8:
	if get_box_dims==6:
	  L_x=float(x[1])-float(x[0])
	elif get_box_dims==7:
	  L_y=float(x[1])-float(x[0])
	  area=L_x*L_y
	elif get_box_dims==8:
	  L_box=float(x[1])-float(x[0])
      	  bin=np.array([((-L_box/2.)+qpp*(L_box/float(Numbins))) for qpp in range(Numbins)]+[L_box*0.5])
	  L_bin = bin[1]-bin[0]
	  z_wall=-0.5*L_box+r_wall
	  
      if not x[0]=="ITEM:":
	x=[float(xx) for xx in x]
	if len(x)<=2:
	  i+=1
	  if i==6: #This means an entire timestep worth of data is had
	    if not (k % nth): #This means data is collected every nth timestep.
	    
		    BA_xyz=sorted(BA_xyz, key=itemgetter(0))
		    BC_xyz=sorted(BC_xyz, key=itemgetter(0))

		    BA_z=np.array(BA_z)
	            BC_z=np.array(BC_z)

		    PosA,PosC = [],[]
		    for (x,y) in zip(BA_xyz,BC_xyz):
		    	PosA.append(x[1])
		    	PosC.append(y[1])
    		    Pos = np.array(PosA+PosC)
    		    PosA = np.array(PosA)
    		    PosC = np.array(PosC)
	            
	            n_obs+=1
		    if k==0: #This indicates that the first timestep worth of data was collected
			
		    ##Collect Pos0 here
		      Pos0 = Pos.copy()
		      MSD = np.sum((Pos - Pos0)**2,axis = 1).mean()
		      MSDs.append(MSD)   

		      Pos0A = PosA.copy()	
		      MSD = np.sum((PosA - Pos0A)**2,axis = 1).mean()
		      MSDAs.append(MSD)   	    

		      Pos0C = PosC.copy()	
		      MSD = np.sum((PosC - Pos0C)**2,axis = 1).mean()
		      MSDCs.append(MSD)   

      		      hist_A = MyHist(BA_z,bin)
		      Sum_A=hist_A	

      		      hist_C = MyHist(BC_z,bin)
		      Sum_C=hist_C

		      del_time = time_new-time_old
		      
		      k+=1
		    else:
#		      if k==1:

		      
		      MSD = np.sum((Pos - Pos0)**2,axis = 1).mean()
		      MSDs.append(MSD)    

		      MSD = np.sum((PosA - Pos0A)**2,axis = 1).mean()
		      MSDAs.append(MSD)    	

		      MSD = np.sum((PosC - Pos0C)**2,axis = 1).mean()
		      MSDCs.append(MSD)   

      		      hist_A = MyHist(BA_z,bin)
		      Sum_A+=hist_A	

      		      hist_C = MyHist(BC_z,bin)
		      Sum_C+=hist_C

		    #The following is a progress meter
		    if not n_obs % 1000:
		    	print filename,n_obs,n_obs/float(50000)
			
		    BA_xyz=[]
		    BC_xyz=[]
		    BA_z,BC_z=[],[]
		    i=1
        elif x[1]==1:  #this is only for the two type model
        	BA_z.append(x[4])
		point=[x[0],x[5:8]]
		BA_xyz.append(point)
        elif x[1]==2: #this is only for the two type model
        	BC_z.append(x[4])
		point=[x[0],x[5:8]]
		BC_xyz.append(point)


    #And to get the last (tricky) snapshot... - 4/27/11 @ 3.30p
    n_obs+=1
    MSD = np.sum((Pos - Pos0)**2,axis = 1).mean()
    MSDs.append(MSD)    

    MSD = np.sum((PosA - Pos0A)**2,axis = 1).mean()
    MSDAs.append(MSD)    	

    MSD = np.sum((PosC - Pos0C)**2,axis = 1).mean()
    MSDCs.append(MSD)  

    hist_A = MyHist(BA_z,bin)
    Sum_A+=hist_A	

    hist_C = MyHist(BC_z,bin)
    Sum_C+=hist_C


    A_count=Sum_A/float(n_obs)
    C_count=Sum_C/float(n_obs)

    BA_den=Sum_A/(float(n_obs)*area*L_bin) #Unnormalized
    BC_den=Sum_C/(float(n_obs)*area*L_bin) #Unnormalized
    n0_BA=np.mean(BA_den[Numbins/2-2:Numbins/2+2])
    n0_BC=np.mean(BC_den[Numbins/2-2:Numbins/2+2])
    n0=(n0_BA+n0_BC)/2. #I think this is plenty fine as long as cations and anions differ only in +/- charge
    lambda_D=(8.*np.pi*n0*Bjerrum)**(-0.5) # added valency on 02/18/10 14:20:23 

    Diff_output=file("Diffusion_Results.txt","a")   
    Diff_output.write(filename)   

    Diff_output.write("\n(lamD,Lbox,dt)= (%1.5f,%1.5f,%1.5f)" % (lambda_D,L_box,LAMMPS_dt))
    for (t,both,aa,bb) in zip([del_time*x for x in range(len(MSDs))],MSDs,MSDAs,MSDCs):
      Diff_output.write("\n%1.1f\t\t%1.8f\t\t%1.8f\t\t%1.8f" % (t,both,aa,bb))
    Diff_output.close()

    return
#SSDiffusion('1000_tr_DEBU_3.776.txt',15.)
#SSDiffusion('1000_tr_zeta_3.776.txt',15.)
#SSDiffusion('1000_Dab_zeta_0.000.gz',0.1)
#SSDiffusion('1000_SD_zeta_0.500.gz',0.1)


#The following function si redundant and should eventually be erased (Same with Q_EDL)
def TransientFree(filename,Bjerrum,r_ion):
    """As is, this function must be double checked for functionality.
	It should be able to calculate Q_EDL=Q_EDL(time).

	***This just gets free charge desnity

Input
    filename: str() the name of the file with MD trajectory data
    Bjerrum:  float() the Bjerrum length that was set for the system
    r_ion:    float() the WCA radius of the ion OR
    		list() of WCA radii r_ion =[smaller,larger]
    M_per_slice: number of insertions per bin for muex_EV calculations
    ##CONSIDER ADDING VALENCY AS A FUCNTION INPUT!
Output:
    None, only printed statements that are written to a .txt file from the command line
"""
    #These are hard coded values that cannot be determined from Data[]
    valency=1.0
    r_wall=1.0
    eps_wall=1.0
    eps_WCA=1.0
    Numbins = 502
    M_per_slice = 10

    if type(r_ion) != type(0.5):
    	print 'Output code must be adjusted before analyzing this data!',fail
    	

    filename_orignal=filename

    if filename[-3:]=='.gz':
        Data=gzip.GzipFile(fileobj=open(filename,'rb'))
	filename=filename[0:len(filename)-2]+'txt'
    else:
        Data=file(filename,"r")

    print '\n'
    print filename,1

    dielectric=(4*np.pi*Bjerrum)**-1 #This used to be dielectric=Bjerrum**-1

    ze=str(filename[-9:-4])
#    print 'The goal as of 11/05/10 16:08:46 should be to feed Sigma_s via the LAMMPS filename.'
	
    if	Bjerrum==0.1:
		if ze=='0.000':
			Sigma_s=0.0*dielectric
		if ze=='0.500':
			Sigma_s=0.03368*dielectric
		if ze=='1.000':
			Sigma_s=0.06948*dielectric
		if ze=='1.500':
			Sigma_s=0.10964*dielectric
		if ze=='1.750':
			Sigma_s=0.13213*dielectric
		if ze=='2.000':
			Sigma_s=0.15669*dielectric
		if ze=='2.500':
			Sigma_s=0.21359*dielectric
		if ze=='2.750':
			Sigma_s=0.24682*dielectric
		if ze=='3.000':
			Sigma_s=0.28390*dielectric
		if ze=='3.500':
			Sigma_s=0.37206*dielectric						
		if ze=='4.000':
			Sigma_s=0.48358*dielectric
		if ze=='5.000':
			Sigma_s=0.80669*dielectric
		if ze=='6.000':
			Sigma_s=1.33572*dielectric		
		if ze=='7.000':
			Sigma_s=2.20568*dielectric
    else:
	Sigma_s=0.0000
	print 'Sigma_s was set to a default value of 0.000'
	#This code is as new as 11/05/10 16:07:32 
    
    nth=1     # new addition, untested (prob need to add a k+=1 near qz_avg) - 02/10/10 13:17:25 , also never actually used... - 11/18/10 11:47:02    
    n_obs=0
    i=0 
    get_time=0	
    get_box_dims=0
    time_old=-1000
    time_new=time_old
    charges=1
##    cycle_time = -1000
    
    #Create empty vectors to fill
    BA_Data,BA_z=[],[]
    BC_Data,BC_z=[],[]
    Pos_A,Pos_C=[],[]
    Q_DL=[]
    Q_L=[]
    Q_R=[]
    Psi_tr=[]
    Efield_tr=[]
    Sum_A_tr=[]
    Sum_C_tr=[]
    NetFree_tr=[]
    exp_muex_EV_tr=[]
    exp_1muex_EV_tr=[]
    exp_2muex_EV_tr=[]
        
    t_index=0

    jj=0

    for y in Data:    
      x=y.split()
      
      ##This gets the current time_step
      if get_time==1:
        time_old=time_new
	time_new=int(x[0])
      	get_time=0	
      if len(x)==2:
      	if x[1]=='TIMESTEP':
     		get_time=1

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
    	  bin=np.array([((-L_box/2.)+qpp*(L_box/float(Numbins))) for qpp in range(Numbins)]+[L_box*0.5])
	  Shell_bin=np.array([value+bin[len(bin)-1] for value in bin])
	  L_bin = Shell_bin[1]-Shell_bin[0]
	  z_wall=-0.5*L_box+r_wall

      if not x[0]=="ITEM:":
	x=[float(xx) for xx in x]
	if len(x)<=2:
	  i+=1
	  if i==6: #This means an entire timestep worth of data is had
	    n_obs+=1

	    BA_z=np.array(BA_z)
	    BC_z=np.array(BC_z)
	    NA_L= len(BA_z[BA_z<0.])
	    NA_R= len(BA_z[BA_z>=0.])
	    NC_L= len(BC_z[BC_z<0.])
	    NC_R= len(BC_z[BC_z>=0.])
	    mean_qdl = ((NA_L-NC_L) + (NC_R-NA_R))*0.5
	    q_r = (NC_R-NA_R)
	    q_l = (NA_L-NC_L)

	    if charges==1:  #This indicates only 1 charge cycle is completed.
	    	Q_DL.append(mean_qdl)
	    	Q_L.append(q_l)
	    	Q_R.append(q_r)

	    	##Counts (for densities) initilization 
		hist_A = MyHist(BA_z,bin)            
		Sum_A_tr.append(np.array(hist_A))

		hist_C = MyHist(BC_z,bin)   
		Sum_C_tr.append(np.array(hist_C))

		freet=np.array(hist_C)-np.array(hist_A)
		NetFree_tr.append(freet)
	    else:
	    	Q_DL[t_index]+=mean_qdl 
	    	Q_L[t_index]+=q_l
	    	Q_R[t_index]+=q_r

		hist_A = MyHist(BA_z,bin)            
		Sum_A_tr[t_index]+=np.array(hist_A)

		hist_C = MyHist(BC_z,bin)   
		Sum_C_tr[t_index]+=np.array(hist_C)

		freet=np.array(hist_C)-np.array(hist_A)
		NetFree_tr[t_index]+=freet


		##This is old code that works only for non-parallelized Transient codes
#	    if time_old>time_new:
#		charges+=1
#		print filename,charges,t_index
#		t_index=0 	
#            else:
#		t_index+=1	


		##This bit of code is robust and can handle any type of transient data input
	    if charges>1:
		if t_index<charge_time:
			t_index+=1
		else:
			charges+=1
			print filename,charges,t_index
			t_index=0 
			
	    else:	
		    if time_old>time_new:
			charge_time = t_index
			charges+=1
			print filename,charges,t_index
			print 'charge time is ',charge_time
			t_index=0 	
		    else:
			t_index+=1
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
    BA_z=np.array(BA_z)
    BC_z=np.array(BC_z)
    N_tot = len(BA_z)+len(BC_z)
    NA_L= len(BA_z[BA_z<0.])
    NA_R= len(BA_z[BA_z>=0.])
    NC_L= len(BC_z[BC_z<0.])
    NC_R= len(BC_z[BC_z>=0.])
    mean_qdl = ((NA_L-NC_L) + (NC_R-NC_R))*0.5
    q_r = (NC_R-NA_R)
    q_l = (NA_L-NC_L)

    Q_DL[t_index]+=mean_qdl
    Q_L[t_index]+=q_l
    Q_R[t_index]+=q_r

    hist_A = MyHist(BA_z,bin)            
    Sum_A_tr[t_index]+=np.array(hist_A)

    hist_C = MyHist(BC_z,bin)   
    Sum_C_tr[t_index]+=np.array(hist_C)

    freet=np.array(hist_C)-np.array(hist_A)
    NetFree_tr[t_index]+=freet

###    #Average the results:
    Q_DL=-np.array(Q_DL)/float(charges) # (negative sign just to flip curves over x-axis)
    Q_L=-np.array(Q_L)/float(charges) # (negative sign just to flip curves over x-axis)
    Q_R=-np.array(Q_R)/float(charges) # (negative sign just to flip curves over x-axis)	
    Sum_A_tr = np.array(Sum_A_tr) / float(charges)
    Sum_C_tr = np.array(Sum_C_tr) / float(charges)
    NetFree_tr = np.array(NetFree_tr) / float(charges)

  
    print "Writing analyzed data to file..."
    total_prefix='Q_meanLR_' + filename[4:len(filename)-4]
    total_output=file(total_prefix+"_DATA.txt","w")
    for (q,l,r) in zip(Q_DL,Q_L,Q_R):
    	total_output.write("%-1.5f\t\t%-1.5f\t\t%-1.5f\n" % (q,l,r))    
    total_output.close()

		##Convert htis to proper outpit
    Xi = Sigma_s*2.*np.pi*Bjerrum**2
    newname = 'Analyzed_trFREE_'+str(N_tot)+'_'+str(Xi)+'_'+str(Bjerrum)+'_'+str(2*r_ion)+'_'+str(L_box)+'_'+str(round(np.sqrt(area),2))
    total_output= file(newname+".txt","w") 

    for tcount in np.arange(len(Q_DL)):
    	A_count = Sum_A_tr[tcount]
    	C_count = Sum_C_tr[tcount]
    	F_count = NetFree_tr[tcount]

    	ql = Q_L[tcount]
    	qr = Q_R[tcount]

    	A_den = np.mean(A_count[Numbins/2-2:Numbins/2+2])/(area*L_bin)
    	C_den = np.mean(C_count[Numbins/2-2:Numbins/2+2])/(area*L_bin)
    	n0 = (A_den + C_den)*0.5
    	lambda_D=(8.*np.pi*n0*Bjerrum)**(-0.5)
    	    	    	
    	for (z,Np,Nm,f) in zip(Shell_bin,A_count.tolist()+[ql],C_count.tolist()+[lambda_D],F_count.tolist()+[qr]):
		total_output.write("%-1.5f\t\t%-1.5f\t\t%-1.5f\t\t%-1.5f\n" % (z,Nm,Np,f))
    total_output.close()

    print 'Finished ' + filename_orignal
    return
temp_rion =0.15#*10
#TransientFree('1000_tr_zeta_4.000.gz',0.1,temp_rion)
#TransientFree('1000_tr_zeta_3.500.gz',0.1,temp_rion)
#TransientFree('1000_tr_zeta_3.000.gz',0.1,temp_rion)
#TransientFree('1000_tr_zeta_2.750.gz',0.1,temp_rion)
#TransientFree('1000_tr_zeta_2.500.gz',0.1,temp_rion)
#TransientFree('1000_tr_zeta_2.000.gz',0.1,temp_rion)
#TransientFree('1000_tr_zeta_1.750.gz',0.1,temp_rion)
#TransientFree('1000_tr_zeta_1.500.gz',0.1,temp_rion)
#TransientFree('1000_tr_zeta_1.000.gz',0.1,temp_rion)
#TransientFree('1000_tr_zeta_0.500.gz',0.1,temp_rion)
#TransientFree('1000_tr_zeta_0.500.txt',0.1,temp_rion)

def Transient_Q_EDL(filename,Bjerrum,r_ion):
    """As is, this function must be double checked for functionality.
	It should be able to calculate Q_EDL=Q_EDL(time).
	This will serve as the template to calculate other spatially and
		temporally varying quanties: voltage(t,z), concentration(t,z), etc. 

Input
    filename: str() the name of the file with MD trajectory data
    Bjerrum:  float() the Bjerrum length that was set for the system
    r_ion:    float() the WCA radius of the ion OR
    		list() of WCA radii r_ion =[smaller,larger]
    M_per_slice: number of insertions per bin for muex_EV calculations
    ##CONSIDER ADDING VALENCY AS A FUCNTION INPUT!
Output:
    None, only printed statements that are written to a .txt file from the command line
"""
    #These are hard coded values that cannot be determined from Data[]
    valency=1.0
    r_wall=1.0
    eps_wall=1.0
    eps_WCA=1.0

    filename_orignal=filename

    if filename[-3:]=='.gz':
        Data=gzip.GzipFile(fileobj=open(filename,'rb'))
	filename=filename[0:len(filename)-2]+'txt'
    else:
        Data=file(filename,"r")

    print '\n'
    print filename

    dielectric=(4*np.pi*Bjerrum)**-1 #This used to be dielectric=Bjerrum**-1

#    graph_prefix=str(Numbins) + filename[4:len(filename)-4]

    #Create empty vectors to fill
    BA_Data,BA_z=[],[]
    BC_Data,BC_z=[],[]
#    Pos_A,Pos_C=[],[]

    ze=str(filename[-9:-4])
#    print 'The goal as of 11/05/10 16:08:46 should be to feed Sigma_s via the LAMMPS filename.'
	
    if	Bjerrum==0.1:
		if ze=='0.000':
			Sigma_s=0.0*dielectric
		if ze=='0.500':
			Sigma_s=0.03368*dielectric
		if ze=='1.000':
			Sigma_s=0.06948*dielectric
		if ze=='1.500':
			Sigma_s=0.10964*dielectric
		if ze=='1.750':
			Sigma_s=0.13213*dielectric
		if ze=='2.000':
			Sigma_s=0.15669*dielectric
		if ze=='2.500':
			Sigma_s=0.21359*dielectric
		if ze=='2.750':
			Sigma_s=0.24682*dielectric
		if ze=='3.000':
			Sigma_s=0.28390*dielectric
		if ze=='3.500':
			Sigma_s=0.37206*dielectric						
		if ze=='4.000':
			Sigma_s=0.48358*dielectric
		if ze=='5.000':
			Sigma_s=0.80669*dielectric
		if ze=='6.000':
			Sigma_s=1.33572*dielectric		
		if ze=='7.000':
			Sigma_s=2.20568*dielectric
    else:
	Sigma_s=0.0000
	print 'Sigma_s was set to a default value of 0.000'
	#This code is as new as 11/05/10 16:07:32 
    
    nth=1     # new addition, untested (prob need to add a k+=1 near qz_avg) - 02/10/10 13:17:25 , also never actually used... - 11/18/10 11:47:02    
    n_obs=0
    i=0 
    get_time=0	
    get_box_dims=0
    time_old=-1000
    time_new=time_old
    charges=1
    Q_DL=[]
    t_index=0

    jj=0
    for y in Data:    
      x=y.split()
      
      ##This gets the current time_step
      if get_time==1:
        time_old=time_new
	time_new=int(x[0])
      	get_time=0	
      if len(x)==2:
      	if x[1]=='TIMESTEP':
     		get_time=1

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

#      print 'QEDL'
#      print x
      if not x[0]=="ITEM:":
	x=[float(xx) for xx in x]
	if len(x)<=2:
	  i+=1
	  if i==6: #This means an entire timestep worth of data is had
	    n_obs+=1

	    BA_z=np.array(BA_z)
	    BC_z=np.array(BC_z)
	    NA_L= len(BA_z[BA_z<0.])
	    NA_R= len(BA_z[BA_z>=0.])
	    NC_L= len(BC_z[BC_z<0.])
	    NC_R= len(BC_z[BC_z>=0.])

	    mean_qdl = ((NA_L-NC_L) + (NC_R-NA_R))*0.5

	    if charges==1:  #This indicates only 1 charge cycle is completed.
	    	Q_DL.append(mean_qdl)
	    	#Insert other quantities here if they are to be averaged...
	    else:
	    	Q_DL[t_index]+=mean_qdl 
	    	
#	    print Q_DL,t_index,time_old,BA_z


	    
	    if time_old>time_new:
		charges+=1
		t_index=0 	
	        print 'Just Q_EDL',filename,charges
            else:
		t_index+=1

	    BA_z=[]
	    BC_z=[]
	    i=1
	elif x[1]==1:
	  BA_z.append(x[4])
	elif x[1]==2:
	  BC_z.append(x[4])

    #And to get the last (tricky) snapshot... - 4/27/11 @ 3.30p

    print jj
    n_obs+=1
    BA_z=np.array(BA_z)
    BC_z=np.array(BC_z)
    NA_L= len(BA_z[BA_z<0.])
    NA_R= len(BA_z[BA_z>=0.])
    NC_L= len(BC_z[BC_z<0.])
    NC_R= len(BC_z[BC_z>=0.])
    mean_qdl = ((NA_L-NC_L) + (NC_R-NC_R))*0.5
    Q_DL[t_index]+=mean_qdl	
    #Average the results (negative sign just to flip curves over x-axis):
    Q_DL=-np.array(Q_DL)/float(charges)
   
    print "Writing analyzed data to file..."

    	##Do averages here!
    total_prefix='Transient_' + filename[4:len(filename)-4]
    total_output=file(total_prefix+"_DATA.txt","w")
    for q in Q_DL:
    	total_output.write("%-1.5f\n" % (q))    
    total_output.close()

    print 'Finished ' + filename_orignal
    return
#temp_rion = 0.15
#Transient_Q_EDL("1000_tr_zeta_0.500.gz",0.1,temp_rion)
#Transient_Q_EDL("1000_tr_zeta_1.000.gz",0.1,temp_rion)
#Transient_Q_EDL("1000_tr_zeta_1.500.gz",0.1,temp_rion)
#Transient_Q_EDL("1000_tr_zeta_1.750.gz",0.1,temp_rion)
#Transient_Q_EDL("1000_tr_zeta_2.000.gz",0.1,temp_rion)
#Transient_Q_EDL("1000_tr_zeta_2.500.gz",0.1,temp_rion)
#Transient_Q_EDL("1000_tr_zeta_2.750.gz",0.1,temp_rion)
#Transient_Q_EDL("1000_tr_zeta_3.000.gz",0.1,temp_rion)
#Transient_Q_EDL("1000_tr_zeta_3.500.gz",0.1,temp_rion)
#Transient_Q_EDL("1000_tr_zeta_4.000.gz",0.1,temp_rion)

def ChargingCycles(filename):
    """The following commented out text must be printed 100X (or however many...)
    so that it can be appended to the run files.


Input
    filename: str() the name of the file with MD trajectory data
    Bjerrum:  float() the Bjerrum length that was set for the system
    r_ion:    float() the WCA radius of the ion OR
    		list() of WCA radii r_ion =[smaller,larger]
    M_per_slice: number of insertions per bin for muex_EV calculations
    ##CONSIDER ADDING VALENCY AS A FUCNTION INPUT!
Output:
    None, only printed statements that are written to a .txt file from the command line
"""


    zeta = filename[-9:-4]
    print 'This ONLY works for CS like, slight modifications required to make this more general\n'
    ##To increase generality, include the following hard set variables as function inputs
    eps = str(10.0000)+'   #=1/Bjerrum'


    if zeta=='0.000':
	Fz=0.0
    if zeta=='0.500':
	Fz=0.03368
    if zeta=='1.000':
	Fz=0.06948
    if zeta=='1.500':
	Fz=0.10964 
    if zeta=='1.750':
	Fz=0.13213 
    if zeta=='2.000':
	Fz=0.15669 
    if zeta=='2.500':
	Fz=0.21359 
    if zeta=='2.750':
	Fz=0.24682 
    if zeta=='3.000':
	Fz=0.28390 
    if zeta=='3.500':
	Fz=0.37206 					
    if zeta=='4.000':
	Fz=0.48358 
    if zeta=='5.000':
	Fz=0.80669
    if zeta=='6.000':
	Fz=1.33572 
    if zeta=='7.000':
	Fz=2.20568 


#    for i in range(1,101):
    for i in range(1,101):
	print '\nclear'
	print '#Beginning of cycle'
	print 'read_restart 	PM_zeta_'+zeta+'_equil_*'
	print 'fix		thermostat all langevin 1.0 1.0 25 3'
	print 'fix		timeintegration all nve'
	print 'dielectric	'+eps
	print 'fix		anode all wall/lj93 zlo -251.0000 1. 1.165 1.0 units box'
	print 'fix		cathode all wall/lj93 zhi 251.0000 1. 1.165 1.0 units box'
	print 'timestep 	0.001'
	print 'thermo_style	custom step temp etotal pe ecoul evdwl press cpu'
	print 'thermo		1000'
	print 'pair_style	lj/cut/coul/long 3.0000 15.000'
	print 'pair_coeff	* * 1. 2.67270'
	print 'kspace_style	pppm 1E-4'
	print 'kspace_modify	slab 3.0'
	print 'pair_modify     shift yes'
	if i==1:
		print 'dump		Movie_dump all atom 1000 tr_zeta'+zeta+'_Movie.all'
	print 'run  		10000'
	print 'dump		transient_dump all custom 1000 '+filename+' id type x y z xu yu zu vx vy vz fx fy fz'
	print 'dump_modify	transient_dump append yes'
	print 'print "_"'
	print 'print "_"'
	print 'print "Part 1: No field cycle '+ str(i)+'/100"'
	print 'print "_"'
	print 'print "_"'
	print 'run  		10000'	
	print 'write_restart   PM_zeta_'+zeta+'_equil_*'
	print 'print "_"'
	print 'print "_"'
	print 'print "Part 2: With field cycle '+ str(i)+'/100"'
	print 'print "_"'
	print 'print "_"'
	print 'fix		A_field A addforce 0.0 0.0 '+str(Fz)
	print 'fix		C_field C addforce 0.0 0.0 -'+str(Fz)
	print 'run  		1000000'
	print '#End of charge cycle\n'
    return

def DischargingCycles(filename):
    """The following commented out text must be printed 100X (or however many...)
    so that it can be appended to the run files.


Input
    filename: str() the name of the file with MD trajectory data
    Bjerrum:  float() the Bjerrum length that was set for the system
    r_ion:    float() the WCA radius of the ion OR
    		list() of WCA radii r_ion =[smaller,larger]
    M_per_slice: number of insertions per bin for muex_EV calculations
    ##CONSIDER ADDING VALENCY AS A FUCNTION INPUT!
Output:
    None, only printed statements that are written to a .txt file from the command line
"""


    zeta = filename[-9:-4]
#    print 'This ONLY works for CS like, slight modifications required to make this more general\n'
    ##To increase generality, include the following hard set variables as function inputs
    eps = str(10.0000)+'   #=1/Bjerrum'


    if zeta=='0.000':
	Fz=0.0
    if zeta=='0.500':
	Fz=0.03368
    if zeta=='1.000':
	Fz=0.06948
    if zeta=='1.500':
	Fz=0.10964 
    if zeta=='1.750':
	Fz=0.13213 
    if zeta=='2.000':
	Fz=0.15669 
    if zeta=='2.500':
	Fz=0.21359 
    if zeta=='2.750':
	Fz=0.24682 
    if zeta=='3.000':
	Fz=0.28390 
    if zeta=='3.500':
	Fz=0.37206 					
    if zeta=='4.000':
	Fz=0.48358 
    if zeta=='5.000':
	Fz=0.80669
    if zeta=='6.000':
	Fz=1.33572 
    if zeta=='7.000':
	Fz=2.20568 


#    for i in range(1,101):
    for i in range(1,101):
	print '\nclear'
	print '#Beginning of cycle'
	print 'read_restart 	PM_zeta_'+zeta+'_equil_*'
	print 'fix		thermostat all langevin 1.0 1.0 1 3'
	print 'fix		timeintegration all nve'
	print 'dielectric	'+eps
	print 'fix		anode all wall/lj93 zlo -251.0000 1. 1.165 1.0 units box'
	print 'fix		cathode all wall/lj93 zhi 251.0000 1. 1.165 1.0 units box'
	print 'fix		A_field A addforce 0.0 0.0 '+str(Fz)
	print 'fix		C_field C addforce 0.0 0.0 -'+str(Fz)	
	print 'timestep 	0.001'
	print 'thermo_style	custom step temp etotal pe ecoul evdwl press cpu'
	print 'thermo		1000'
	print 'pair_style	lj/cut/coul/long 0.3000 15.000'
	print 'pair_coeff	* * 1. 0.267270'
	print 'kspace_style	pppm 1E-4'
	print 'kspace_modify	slab 3.0'
	print 'pair_modify     shift yes'
	if i==1:
		print 'dump		Movie_dump all atom 1000 tr_zeta'+zeta+'_Movie.all'
	print 'run  		10000'
	print 'dump		transient_dump all custom 1000 '+str(i)+'_'+filename+' id type x y z xu yu zu vx vy vz fx fy fz'
	print 'dump_modify	transient_dump append yes'
	print 'print "_"'
	print 'print "_"'
	print 'print "Part 1: With field cycle '+ str(i)+'/100"'
	print 'print "_"'
	print 'print "_"'
	print 'run  		10000'	
	print 'write_restart   PM_zeta_'+zeta+'_equil_*'
	print 'print "_"'
	print 'print "_"'
	print 'print "Part 2: Without field cycle '+ str(i)+'/100"'
	print 'print "_"'
	print 'print "_"'
	print 'unfix		A_field'
	print 'unfix		C_field'
	print 'run  		2500000'
	print '#End of charge cycle\n'
    return
#DischargingCycles('1000_tr_zeta_3.500.txt')
#DischargingCycles('1000_tr_zeta_2.000.txt')


def DischargingCycles_ratios(filename,letter):
    """The following commented out text must be printed 100X (or however many...)
    so that it can be appended to the run files.


Input
    filename: str() the name of the file with MD trajectory data
    Bjerrum:  float() the Bjerrum length that was set for the system
    r_ion:    float() the WCA radius of the ion OR
    		list() of WCA radii r_ion =[smaller,larger]
    M_per_slice: number of insertions per bin for muex_EV calculations
    ##CONSIDER ADDING VALENCY AS A FUCNTION INPUT!
Output:
    None, only printed statements that are written to a .txt file from the command line
"""


    zeta = filename[-9:-4]
#    print 'This ONLY works for CS like, slight modifications required to make this more general\n'
    ##To increase generality, include the following hard set variables as function inputs
    eps = str(10.0000)+'   #=1/Bjerrum'


    if zeta=='0.000':
	Fz=0.0
    if zeta=='0.500':
	Fz=0.03368
    if zeta=='1.000':
	Fz=0.06948
    if zeta=='1.500':
	Fz=0.10964 
    if zeta=='1.750':
	Fz=0.13213 
    if zeta=='2.000':
	Fz=0.15669 
    if zeta=='2.500':
	Fz=0.21359 
    if zeta=='2.750':
	Fz=0.24682 
    if zeta=='3.000':
	Fz=0.28390 
    if zeta=='3.500':
	Fz=0.37206 					
    if zeta=='4.000':
	Fz=0.48358 
    if zeta=='5.000':
	Fz=0.80669
    if zeta=='6.000':
	Fz=1.33572 
    if zeta=='7.000':
	Fz=2.20568 

    if letter=='a':
    	rando = '3'
    elif letter=='b':
    	rando = '4'
    elif letter=='c':
    	rando = '5'
    elif letter=='d':
    	rando = '6'
    	
#    for i in range(1,101):
    for i in range(1,26):
	print '\nclear'
	print '#Beginning of cycle'
	print 'read_restart 	PM_zeta_'+zeta+'_equil'+letter+'_*'
	print 'fix		thermostat A langevin 1.0 1. 5 '+rando
	print 'fix		thermostat C langevin 1.0 1. 0.5 '+rando	
	print 'fix		timeintegration all nve'
	print 'dielectric	'+eps
	print 'fix		anode all wall/lj93 zlo -251.0000 1. 1.165 1.0 units box'
	print 'fix		cathode all wall/lj93 zhi 251.0000 1. 1.165 1.0 units box'
	print 'fix		A_field A addforce 0.0 0.0 '+str(Fz)
	print 'fix		C_field C addforce 0.0 0.0 -'+str(Fz)	
	print 'timestep 	0.001'
	print 'thermo_style	custom step temp etotal pe ecoul evdwl press cpu'
	print 'thermo		1000'
	print 'pair_style	lj/cut/coul/long 3.0000 15.000' #print 'pair_style	lj/cut/coul/long 0.3000 15.000'
	print 'pair_coeff	* * 1. 2.67270' #	print 'pair_coeff	* * 1. 0.267270'
	print 'kspace_style	pppm 1E-4'
	print 'kspace_modify	slab 3.0'
	print 'pair_modify     shift yes'
	if i==1:
		print 'dump		Movie_dump all atom 200 tr_zeta'+zeta+'_Movie.all'
	print 'run  		10000'
	print 'dump		transient_dump all custom 200 '+str(i)+'_'+filename[:-4]+letter+'.txt id type x y z vz'
	print 'print "_"'
	print 'print "_"'
	print 'print "Part 1: With field cycle '+ str(i)+'/100"'
	print 'print "_"'
	print 'print "_"'
	print 'run  		20000'	
	print 'write_restart   PM_zeta_'+zeta+'_equil'+letter+'_*'
	print 'print "_"'
	print 'print "_"'
	print 'print "Part 2: Without field cycle '+ str(i)+'/100"'
	print 'print "_"'
	print 'print "_"'
	print 'unfix		A_field'
	print 'unfix		C_field'
	print 'run  		5000000'
	print '#End of charge cycle\n'
    return
#DischargingCycles_ratios('200_tr_zeta_3.500.txt','d')


def  P2_Plot():
    """
	Plotting routines for publication 2.

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
#    import random	
#    random.shuffle(filenameS)

 
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
					print '\n\n',filename,sigWCA
					
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


    characteristic_lengthS,xlabel=LDs,r'$z/\lambda_{D}$'


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





		#This could surely be consolidated within the lines of code below...
	        z_density = [(x+y)/(L_z*2.) for (x,y) in zip(z_positions[0:len(z_positions)-1],z_positions[1:len(z_positions)])]
	        z_den_theory=ztheory#[z/L_z for z in z_positions]  
		z_density = np.array(z_density)*L_z-0.5*sigWCA #This accounts for all of the wall
		z_density=z_density/L_z

		special_marker='s'
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
#		ax1.errorbar(np.array(z_density[1:])*(L_z/characteristic_length),Nf[1:],yerr=None,marker=special_marker,ms=4.5,color=colors[i],ls=special_ls,label=special_label)


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
		ax1.errorbar(np.array(z_den_theory)*(L_z/characteristic_length),-np.array(rf_T),yerr=None,ls='-',color=color_scale)#,label=r'$\~ \Sigma$'+' = ' + str(round(Sigma_s/(dielectric*temperature/(valency*lam_D)),3)))
		ax1.errorbar(np.array(z_density[1:])*(L_z/characteristic_length),Nf[1:],yerr=None,marker=marker_special,mew=0.4,ms=4.5,color=color_scale,ls='None',label=labelS[i])


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
	        Nf_for_Eff=[]
	        for (y1,y2,z) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos):
	        	if z*L_z>z_wall:  #ONe more thing, try this without the equalds
	        		SigEff.append(Integrand) #There's something wrong with the first value assigned to this function
	        		Integrand=Integrand + 0.5*(y1+y2)*L_bin	
	        		if filename in GC_LB_LD_5_20:
	       				Nf_for_Eff.append(0.5*(y1+y2))	        		    		
	        		else:
	       				Nf_for_Eff.append(y1)
        	SigEffS.append(np.array(SigEff)/(dielectric*temperature/(valency*lam_D)))
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
    ax1.set_xlim(0,4)
    ax1.set_ylim(-0.5,35)
#    ax1.set_ylim(-0.5,15)
#    ax1.set_ylim(-0.5,50)
    xlabel=r'$\~ z = (z-\sigma/2)/\lambda_{D}$'
    ylabel=r'$-\rho/2 q n^{\infty]$'
    ax1.set_xlabel(xlabel,size='small')
    ax1.set_ylabel(ylabel,size='small')
    plt.setp(ax1.get_xticklabels(), fontsize='small')
    plt.setp(ax1.get_yticklabels(), fontsize='small')
    fig.set_size_inches(3.37,2)
#    plt.show()
#    plt.savefig('P1_PRL_F2_raw.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 
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


#	###LaTex output
#    print r"\begin{tabular}{|c|c|c|c|}"
#    print '\hline'
#    print r"$\Sigma_{app}$ & $\lambda_{D}$ & $N_{+}+N{-}$ \zeta_{+}^{measured}\times k_{B}T/q_{\pm} \\ \hline"
##    old_ratio = sigWCA_S[0]/BjerrumS[0]
#    for (Sig_s,lam_D,sig,lam_B,N_p,N_m,V) in zip(SIGMAS,LDs,sigWCA_S,BjerrumS,N_plusS,N_minusS,V_corS):
##      if old_ratio !=sig/lam_B:
##      		print '\n'
##      		old_ratio=old_ratio
#      print r"%1.4f & %1.3f & %i & %1.2f \\ \hline" % (Sig_s,lam_D,int(round(np.sum(N_p)+np.sum(N_m))),V)

    ratioS,legend_key=np.array(BjerrumS)/np.array(LDs),'ele'

###    ##This plots mu_ex^EV(z) for simulation and modified theories
    muexEV_bulkS=[]
    nonDim='yes'
#    Bikerman='yes'
    Bikerman='no'
#    CarnahanStarling='yes'
    CarnahanStarling='no'
    if nonDim=='yes':
      print 'Plotting NonDim mu_ex^EV(z)...'
    else:
      print 'Plotting Dim mu_ex^EV(z)...'	
    if Bikerman=='yes':
		print '\t\t...with Bikerman theory...'
		#Note to user: Bikerman theory must already be evaluated using MATLAB code, which is in ~/sims/PBik_Solver.m
    if CarnahanStarling=='yes':
	      print '\t...with Carnahan-Starling theory...'
    i=-1
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    fig.subplots_adjust(right=0.90)
    for (z_positions,characteristic_length,muex_EV,L_z) in zip(z_S,characteristic_lengthS,muexEV_S,L_zS):
    	i+=1

	if i==len(markers): #This resets the markers
		i=0

	if i==0:
		graphname = 'muex_EV_z' + '.pdf'
		ylabel=r'$\~\mu_{ex}^{EV}$'
		
	z_positions=[z/L_z for z in z_positions] #NonDim distances

        muexEV_bulk=np.mean(muex_EV[len(muex_EV)/2-2:len(muex_EV)/2+2])
        muexEV_bulkS.append(muexEV_bulk)

	##The code below includes bulk values:
	ax1.errorbar(np.array(z_positions)*(L_z/characteristic_length),np.array(muex_EV),yerr=None,marker=markers[i],ms=3.0,color=colors[i],ls='None')
	
#	#The code below subtracts off bulk values:
#	ax1.errorbar(np.array(z_positions)*(L_z/characteristic_length),np.array(muex_EV)-muexEV_bulk,yerr=None,marker=markers[i],ms=3.0,color=colors[i],ls='None')
	
	if Bikerman=='yes':
	  Bik_file='_Bik_zeta_' + ze + '.txt'
	  x_Bik=[[float(x) for x in line.split()] for line in file('x' + Bik_file,"r").readlines()]
	  x_Bik=np.array(x_Bik)*lam_D/characteristic_length
	  c_Bik=[[float(x) for x in line.split()] for line in file('counter' + Bik_file,"r").readlines()]
	  co_Bik=[[float(x) for x in line.split()] for line in file('co' + Bik_file,"r").readlines()]
	  ax1.plot(x_Bik[0],-0.5*np.log(np.array(co_Bik[0])*np.array(c_Bik[0])),color=colors[i],lw=1.5,ls='--') 
	if CarnahanStarling=='yes':
	  CS_file='_CS_zeta_' + ze + '.txt'
          x_CS=[[float(x) for x in line.split()] for line in file('x'+CS_file,"r").readlines()]
          x_CS=np.array(x_CS)*lam_D/characteristic_length
	  c_CS=[[float(x) for x in line.split()] for line in file('counter' + CS_file,"r").readlines()]
	  co_CS=[[float(x) for x in line.split()] for line in file('co' + CS_file,"r").readlines()]	
	  ax1.plot(x_CS[0],-0.5*np.log(np.array(co_CS[0])*np.array(c_CS[0])),color=colors[i],lw=1.5,ls='-.')      
  		##Keep this code for sometime - 02/22/12 14:43:58 
#	  if sig==1:
#	  	ax1.plot(x_CS[0],0.5*np.log(np.array(co_CS[0])*np.array(c_CS[0])),color=colors[i],lw=1.5,ls='-.')
#	  else:
#	  	ax1.plot(x_CS[0],-0.5*np.log(np.array(co_CS[0])*np.array(c_CS[0])),color=colors[i],lw=1.5,ls='-.')
    if xlabel==r'$z/L_{z}$':
      #These must be set manually! 
      test=0
#      ax1.set_xlim(0,0.2)       
#      ax1.set_ylim(0,10) 
    elif xlabel==r'$z/ \sigma_{WCA}$':
      #These must be set manually! 
      test=0
#      ax1.set_xlim(0,0.2)       
#      ax1.set_ylim(0,10)      
    ax1.set_xlabel(xlabel,size='x-large')
    ax1.set_ylabel(ylabel,size='x-large')
#    ax1.legend(loc=0) 
    plt.savefig(graphname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    plt.show()
    plt.close()
#    #Done plotting mu_excess^EV(z)

    print 'Plotting muexEV = muexEV(Phi_bulk)...'
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    i=-1
    j=-1
    PhiWCAs=[]
    Phi_shifted=[]
    Phi_HS=[]
    MAX = -1e9
    for (mexb,Xi,n0,sig_WCA,sigHS) in zip(muexEV_bulkS,XiS,n0s,sigWCA_S,sigHS_S):
    	i+=1
    	j+=1
	if i==len(markers): #This resets the markers
			i=0
    	Phi_WCA = 2.*n0*(np.pi/6)*sig_WCA**3
#    	print '%1.5f\t\t%1.5f' % (Phi_WCA,mexb)
    	MAX = np.max([MAX,Phi_WCA])
    	PhiWCAs.append(Phi_WCA)
    	Phi_shifted.append(Phi_WCA*0.953259114**3)  #Unweighted = 0.981989165  0.954028946 
    	Phi_HS.append(2.*n0*(np.pi/6)*sigHS**3)
    	if j==0:
    		qqq=1
    		ax1.errorbar(Phi_WCA,mexb,yerr=None,mfc='white',mec='blue',marker='o',ms=3.0,color='blue',ls='None',label=r'$\~ \mu_{ex}^{EV}(\sigma_{WCA})$')
#####    		ax1.errorbar(Phi_WCA,mexb,yerr=None,color='blue',ls='-',label=r'$\~ \mu_{ex}^{EV}(\sigma_{WCA})$')
	else:
		qqq=1
    		ax1.errorbar(Phi_WCA,mexb,yerr=None,mfc='white',mec='blue',marker='o',ms=3.0,color='blue',ls='None')
#####    		ax1.errorbar(Phi_WCA,mexb,yerr=None,color='blue',ls='-')
    Phi_theory = np.linspace(0,MAX,1000)
    ax1.errorbar(Phi_shifted,muexEV_bulkS,yerr=None,color='blue',ls='None',marker='*',ms=6.0,label=r'$\~ \mu_{ex}^{EV}(0.953*\sigma_{WCA})$') 

    ax1.errorbar(np.array(Phi_shifted)*(0.984796597/0.953259114)**3,muexEV_bulkS,yerr=None,color='orange',ls='None',marker='^',ms=3.0,label=r'$\~ \mu_{ex}^{EV}(0.955*\sigma_{WCA})$') 
    
    ax1.errorbar(Phi_HS,muexEV_bulkS,yerr=None,color='blue',marker='o',ms=3.0,ls='None',label=r'$\~ \mu_{ex}^{EV}(\sigma_{HS})$') 
    ax1.errorbar(Phi_theory,Phi_theory*(8-9*Phi_theory+3*Phi_theory**2)/(1-Phi_theory)**3,yerr=None,color='k',ls='-.',lw=1.5,label=r'$\mu_{ex}^{CS}(\Phi)$')#marker=markers[i],color=colors[i]    
    ax1.errorbar(Phi_theory,-np.log(1-Phi_theory/0.65),yerr=None,color='r',ls='--',label=r'$-\ln(1-\Phi/0.65)$')#marker=markers[i],color=colors[i]    # /0.65)
    ax1.errorbar(Phi_theory,-np.log(1-Phi_theory),yerr=None,color='k',ls='--',label=r'$-\ln(1-\Phi)$')#marker=markers[i],color=colors[i]    # /0.65)

    ax1.set_xlabel(r'$\Phi^{bulk}$',size='x-small')
    ax1.set_ylabel(r'$\~\mu_{ex}^{EV}(\Phi^{bulk})$',size='x-small')
	##Important fitting text below
#    coeffs = np.polyfit(PhiWCAs,np.exp(-np.array(muexEV_bulkS)),deg=3)
#    x_fit = np.linspace(0,MAX*1.1,100)
#    print coeffs
#    leg_info = r'$-\ln(%s \Phi_{WCA}^{3} + %s \Phi_{WCA}^{2}  %s\Phi_{WCA}+%s)$' % (str(round(coeffs[0],1)),str(round(coeffs[1],1)),str(round(coeffs[2],2)),str(round(coeffs[3],2)))
#    ax1.plot(x_fit,-np.log(np.polyval(np.poly1d(coeffs),x_fit)),color='b',lw=2.0,ls='-',label=r'$-\ln(1-f[\Phi])$')
#    ymean = np.mean(np.exp(-np.array(muexEV_bulkS)))
#    R_sq = 1. - (sum([(yi-fi)**2 for (yi,fi) in zip(np.exp(-np.array(muexEV_bulkS)),np.polyval(np.poly1d(coeffs),PhiWCAs))]) / sum([(yi-ymean)**2 for yi in np.exp(-np.array(muexEV_bulkS))]))
#    print 'R_sq = ',R_sq
#    ax1.text(1,12,'Fitting \~\mu_{ex}^{EV}(\sigma_{WCA}):\nR^2 = %s' % str(R_sq))

    plt.setp(ax1.get_xticklabels(), fontsize='x-small')
    plt.setp(ax1.get_yticklabels(), fontsize='x-small')

    ax1.set_xlim(0,MAX)       
    ax1.set_ylim(0,np.max(muexEV_bulkS))
    
    fig.set_size_inches(3.37,3.5)
    plt.savefig('muexEV_vs_Phi_B.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
##    plt.show()
#    ax1.legend(loc=0)#    fig.set_size_inches(11,7)

    ax1.set_xlim(0,0.05805)
    ax1.set_ylim(0,0.5) 
    plt.savefig('muexEV_vs_Phi_B_LT1.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 
    plt.close()



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




def  AsymReformat():
    """
	Plotting routines for publication 2, Figure 1 ONLY
	This is an attempt at adding modularity to my code.
Input: 
      XXXXfilenameS - A vector of all the filenames to be included in the plots

Output: Plots as meantioned above.
Notes: Would be nice to extende Bikerman theory, too.
"""
    import numpy as np
#    import matplotlib
#    import matplotlib.pyplot as plt
#    from matplotlib.colors import LinearSegmentedColormap
#    from mpl_toolkits.mplot3d import Axes3D
    import sys
#    from matplotlib.ticker import MaxNLocator
    from scipy import interpolate
    from scipy import optimize
    

    filenameS=[]

    for f in [[x for x in line.split()] for line in file(sys.argv[1],"r").readlines()]:
    	filenameS.append(f[0])

 
#    Volts=[]
#    zminS=[]
#    zmaxS=[]
#    Fields=[]
#    N_plusS=[]
#    N_minusS=[]
#    originalN_plusS,originalN_minusS = [],[]
#    LDs=[]
#    n0s=[]
#    SIGMAS=[]
#    BjerrumS=[]
#    sigWCA_S=[]
#    N_totS=[]
#    XiS=[]
#    L_zS=[]
#    areaS=[]
#    muexEV_S=[]
#    muexEVtot_S=[]
#    z_S,z_originalS=[],[]	#These are the z_pos for plotting is for voltage, field, and muex
#    L_binS=[]
#    PhiBulkS=[]
#    sigHS_S=[]
#    PhiWCAtot_S=[]
    
    Bikerman='n0'

   #These variables below are (reasonably) assumed to always equal these values
    z_wall=1.	
    temperature = 1.0
    valency = 1.0
    eps_wall=1.0

    #This loop build all the necessary data to superimpose plots. 
    print 'Loading Data for ',len(filenameS),' simulations...'
    i=0
    for filename in filenameS:  	
#    	print filename
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
#				N_totS.append(N_tot)
			if j==2:
				Xi=float(value)
#				XiS.append(Xi)
			if j==3:
				Bjerrum=float(value)
#				BjerrumS.append(Bjerrum)
			if j==4:
				sigWCA=float(value)
			
				if sigWCA==1.0:
#					print filename
					a=1 #?
#				sigWCA_S.append(sigWCA)

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
#					sig_HS = Noro_WCA(1,sigWCA, 10**-20, sigWCA, 10**5) #This should be used when needed!	
					sig_HS = 1E9 ##THis is turned off for faster code
					if i==1:
						print '\t**Hard sphere diameters are being calculated. This takes time.'	
#					print '(WCA,HS) = (%1.3f,%1.3f)' % (sigWCA,sig_HS)
#					print filename,sig_HS
#				sigHS_S.append(sig_HS)	
			if j==5:#This length will actually be determined from LAMMPS files
				test2=0 #The L_z that includes wall thickness is more accurately obtained from z_positions below	
				##Keep these uncommented and eventually delete - 03/07/12 14:34:06 
#				L_z=float(value) #This length DOES NOT include the walls... - 03/07/12 13:02:14 
#				L_zS.append(L_z) 
			if j==6:
				test2=0 ##This Lxy is already embeded in the files
			if j==7:
				sigBIG = float(value)
			value=''

	L_xy=float(value[:-4])
	area = L_xy*L_xy
    	z_50=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()])
    
	#Trim incoming data
	z_positions=z_50[:,0]
        SIGMA = z_positions[-1]
        L_bin=z_positions[1]-z_positions[0]
        L_binS = L_bin
        zpos = z_positions[:-1].tolist()
        z_originalS = zpos

        L_z = zpos[-1]
        
	Psi_tot=z_50[:,1]
	zmin = Psi_tot[-1]
	Volt = Psi_tot[:-1].tolist()
	
	E_tot=z_50[:,2]
	zmax = E_tot[-1]

	Field = E_tot[:-1].tolist()

	N_plus=z_50[:,3]
	originalN_plus = N_plus[:-2]

	N_minus=z_50[:,4]
	originalN_minus = N_minus[:-2]
	lam_D=N_minus[-2]
	n0 = (8.*np.pi*Bjerrum*lam_D**2)**-1

	PhiBulk = 2.*n0*(np.pi/6)*sigWCA**3
	

	muex_EV1=z_50[:,5]
	muexEV1tot = muex_EV1[-1]
	muex_EVz_small = muex_EV1[:-1].tolist()

	muex_EV2=z_50[:,6]
	muexEV2tot = muex_EV2[-1]
	muex_EVz_big = muex_EV2[:-1].tolist()


	PhiWCA = ((np.array(N_plus[:-2]) + np.array(N_minus[:-2]))/((L_xy*L_xy*L_bin))*(np.pi/6)*sigWCA**3)


	print 'HEre ',len(zpos),len(Volt),len(Field),len(originalN_plus),len(originalN_minus),len(muex_EVz_small),len(muex_EVz_big)
	print 'zpos LEFT = ',zpos[:len(zpos)/2+1]
	print len(zpos[:len(zpos)/2])

	print 'zpos RIGHT = ',zpos[len(zpos)/2:]
	print len(zpos[len(zpos)/2:])
	
##	print '\nvolt = ',Volt
#	print '\nnplus = ',originalN_plus
#	print '\nleft half nplus = ',originalN_plus[:len(originalN_plus)/2]
#	print len(originalN_plus[:len(originalN_plus)/2])
#	print '\nright half nplus = ',originalN_plus[len(originalN_plus)/2:]
#	print len(originalN_plus[len(originalN_plus)/2:])
##	print '\nexEV = ',muex_EVz_big

	##Assemble profiles from left side of the box
#	zposL = zpos[:len(zpos)/2+1]
#	VoltL = Volt[:len(Volt)/2+1]
#	FieldL = Field[:len(Field)/2+1]
#	counterL = originalN_plus[:len(originalN_plus)/2]
#	coL = originalN_minus[:len(originalN_minus)/2]
#	counterexL = muex_EVz_small[:len(muex_EVz_small)/2+1]
#	coexL = muex_EVz_big[:len(muex_EVz_big)/2+1]

	counterL,coR = [],[]
	coL,counterR = [],[]
	vL,vR = [],[]
	for (z,psi,efield,np,nm,musmall,mubig) in zip(zpos,Volt,Field,originalN_plus,originalN_minus,muex_EVz_small,muex_EVz_big):
		print z,'\t',psi,'\t',efield,'\t',np,'\t',nm,'\t',musmall,'\t',mubig
		if len(counterL)<len(zpos[:len(zpos)/2]):
			counterL.append(np)
			coL.append(nm)
		else:
			counterR.append(nm)
			coR.append(np)


		if len(vL)<=len(zpos[:len(zpos)/2]):
			vL.append(psi)
		else:
			vR.append(psi)
	vR.append(Volt[-1])

	vR = [vL[-1]] + vR
	
	print '\ncounterL = ',counterL
	print '\ncoR = ',coR
	print len(counterL),len(coR),len(counterL+coR)
	print '\nvolL = ',vL
	print '\nvolR = ',vR
	print len(vL),len(vR),len(vL+vR),len(Volt)

	for (x,y) in zip(vR,coR):
		print x,y
	print vR[-1]
	

#	print '\n\n\nzposL  = ',zposL
#	print '\ncounterL = ',counterL
#	print '\nleft values = ',len(zposL),len(VoltL),len(FieldL),len(counterL),len(coL),len(counterexL),len(coexL)
#	print '\nVoltL = ',VoltL



######	    AvInWallZ = r_wall - 0.5*(abs(np.mean(wallZ_L)+L_z*0.5)+(np.mean(wallZ_R)-L_z*0.5))
######	    Xi = Sigma_s*2.*np.pi*Bjerrum**2
######	    newname = 'Analyzed_GC_'+str(N_tot)+'_'+str(Xi)+'_'+str(Bjerrum)+'_'+str(2*r_ion[0])+'_'+str(L_z)+'_'+str(round(np.sqrt(area),2))+'_'+str(2*r_ion[1])
######	    total_output= file(newname+".txt","w") 
######	    for (z,V,E,Np,Nm,muexEV,muexHS) in zip(Shell_bin.tolist()+[Sigma_s],Psi_tot.tolist()+[z_min],E_tot.tolist()+[z_max],A_count.tolist()+[lambda_D,r_ion[0]/r_ion[1]],C_count.tolist()+[area,Bjerrum],muex_EV1.tolist()+[AvInWallZ],muex_EV2.tolist()+[AvInWallZ]):
######		total_output.write("%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\n" % (z,V,E,Nm,Np,muexEV,muexHS))
######	    total_output.close() 

	

	##This stuff below should happen anyway after new file is loaded in...




    return
#AsymReformat()


def  P2F1():
    """
	Plotting routines for publication 2, Figure 1 ONLY
	This is an attempt at adding modularity to my code.
Input: 
      XXXXfilenameS - A vector of all the filenames to be included in the plots

Output: Plots as meantioned above.
Notes: Would be nice to extende Bikerman theory, too.
"""
    import numpy as np
#    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.mlab import griddata
    from matplotlib.ticker import MaxNLocator,FormatStrFormatter
    from mpl_toolkits.mplot3d import Axes3D
    import sys
    from scipy import interpolate
    from scipy import optimize
#    from scipy.interpolate import UnivariateSpline
    from scipy.optimize import fsolve
#####    from scipy.optimize import root
    from scipy import stats
    
    
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
    oldF3_Data=['Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1312_0.0241789536991_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt']
    F3_Data=['Analyzed_GC_1312_0.0241789536991_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt','Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt']  
#    filenameS=GC_sig8
    filenameS=[]

    for f in [[x for x in line.split()] for line in file(sys.argv[1],"r").readlines()]:
    	filenameS.append(f[0])
#    	if '_P1_' in f[0]:
#    		print f[0]
#    	print 'submit "P1.py %s"\nsleep 2' % str(f[0])
#    print "THIS IS BEING REVERSED"
#    filenameS.reverse()
#    if len(filenameS)>100:
#    import random	
#    random.shuffle(filenameS)

 
    Volts=[]
    zminS=[]
    zmaxS=[]
    Fields=[]
    N_plusS=[]
    N_minusS=[]
    originalN_plusS,originalN_minusS = [],[]
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
    z_S,z_originalS=[],[]	#These are the z_pos for plotting is for voltage, field, and muex
    L_binS=[]
    PhiBulkS=[]
    sigHS_S=[]
    PhiWCAtot_S=[]
    
    Bikerman='n0'

   #These variables below are (reasonably) assumed to always equal these values
    z_wall=1.	
    temperature = 1.0
    valency = 1.0
    eps_wall=1.0

    #This loop build all the necessary data to superimpose plots. 
    print 'Loading Data for ',len(filenameS),' simulations...'
    i=0
    for filename in filenameS:  	
#    	print filename
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
#					print 'do nothing'
#					sig_HS = Noro_WCA(1,sigWCA, 10**-20, sigWCA, 10**5) #This should be used when needed!	
					sig_HS = 1E9 ##THis is turned off for faster code
#					if i==1:
#						print '\t**Hard sphere diameters are being calculated. This takes time.'	
#					print '(WCA,HS) = (%1.3f,%1.3f)' % (sigWCA,sig_HS)
#					print filename,sig_HS
				sigHS_S.append(sig_HS)	
			if j==5:#This length will actually be determined from LAMMPS files
				test2=0 #The L_z that includes wall thickness is more accurately obtained from z_positions below	
				##Keep these uncommented and eventually delete - 03/07/12 14:34:06 
#				L_z=float(value) #This length DOES NOT include the walls... - 03/07/12 13:02:14 
#				L_zS.append(L_z) 
			value=''
	L_xy=float(value[:-4])
#	print L_xy
	areaS.append(L_xy*L_xy)

    	z_50=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()])
    
	#Trim incoming data
	z_positions=z_50[:,0]
        SIGMAS.append(z_positions[-1])  
        L_bin=z_positions[1]-z_positions[0]
        L_binS.append(L_bin)
        zpos = z_positions[:-1].tolist()
        z_originalS.append(zpos)

        L_zS.append(zpos[-1])
        
	Psi_tot=z_50[:,1]
	zminS.append(Psi_tot[-1])
	Volts.append(Psi_tot[:-1].tolist())
	
	E_tot=z_50[:,2]
	zmaxS.append(E_tot[-1])
#	print filename,'\t',E_tot[-1],'\t',Psi_tot[-1],'\t',L_bin,'\t',zpos[-1]
	Fields.append(E_tot[:-1].tolist())

	N_plus=z_50[:,3]
	originalN_plusS.append(N_plus[:-2])	

	N_minus=z_50[:,4]
	originalN_minusS.append(N_minus[:-2])
	lam_D=N_minus[-2]
	LDs.append(N_minus[-2])
	n0 = (8.*np.pi*Bjerrum*lam_D**2)**-1
	n0s.append(n0)

	PhiBulkS.append(2.*n0*(np.pi/6)*sigWCA**3)
	

	muex_EV=z_50[:,5]
	muexEVtot_S.append(muex_EV[-1])
	muex_EVz = muex_EV[:-1].tolist()

	PhiWCA = ((np.array(N_plus[:-2]) + np.array(N_minus[:-2]))/((L_xy*L_xy*L_bin))*(np.pi/6)*sigWCA**3)


	if '_P1_' in filename:  #Could technically do this for all (and base it on exact zmin/zmax)...
		print filename
		ntot_wall = 0.
		np_wall,nm_wall = 0.,0.
		P1_PhiWCA,P1_muex_EVz,exEV_wallS,P1_z = [],[],[],[]
		P1_Nm,P1_Np=[],[]
		for (z,phi,ev,np1,nm1) in zip(zpos,PhiWCA,muex_EVz,N_plus[:-2],N_minus[:-2]):
			if z<=1.:
				ntot_wall = ntot_wall + (phi*(L_xy*L_xy*L_bin))*(np.pi/6)*sigWCA**3
				np_wall+=np1
				nm_wall+=nm1
				exEV_wallS.append(ev)
				zwall = [0.0,z]
			else:
				P1_PhiWCA.append(phi)
				P1_muex_EVz.append(ev)
				P1_z.append(z)
				P1_Nm.append(nm1)
				P1_Np.append(np1)

		P1_z = zwall + P1_z
				
		PhiWall = (ntot_wall/(L_xy*L_xy*(1.-0.299602925)))*(np.pi/6)*sigWCA**3 #Could technically use exact zmin/zmax...
		P1_PhiWCA = [PhiWall] + P1_PhiWCA

		NmWall = nm_wall / (L_xy*L_xy*(1.-0.299602925))
		P1_Nm = [NmWall] + P1_Nm
		
		NpWall = np_wall / (L_xy*L_xy*(1.-0.299602925))
		P1_Np = [NpWall] + P1_Np
		
		expEV_wallS = np.exp(-np.array(exEV_wallS))
		EV_wall = -np.log(sum(expEV_wallS)/len(expEV_wallS))
		P1_muex_EVz = [EV_wall] + P1_muex_EVz

		z_S.append(P1_z)
		PhiWCAtot_S.append(P1_PhiWCA)
		muexEV_S.append(P1_muex_EVz)
        	N_minusS.append(P1_Nm)
		N_plusS.append(P1_Np)
	else:
		PhiWCAtot_S.append(PhiWCA)
		muexEV_S.append(muex_EVz)
	        z_S.append(zpos)
        	N_minusS.append(N_minus[:-2])
		N_plusS.append(N_plus[:-2])


    ##Turns out Np and Nm notation is somehow switched... This is fixed below:
#    print "TEMPORARY FOR LAMMPS DEBUG"
    temp = N_plusS
    N_plusS = N_minusS
    N_minusS = temp

    temp = originalN_plusS
    originalN_plusS = originalN_minusS
    originalN_minusS = temp


#    print r"\begin{tabular}{|c|c|c|c|}"
#    print '\hline'
#    print r"$\Sigma_{app}$ & $\lambda_{D}$ & $N_{+}+N{-}$ \zeta_{+}^{measured}\times k_{B}T/q_{\pm} \\ \hline"
##    old_ratio = sigWCA_S[0]/BjerrumS[0]
#    for (Sig_s,lam_D,sig,lam_B,N_p,N_m,V) in zip(SIGMAS,LDs,sigWCA_S,BjerrumS,N_plusS,N_minusS,Volts):
##      if old_ratio !=sig/lam_B:
##      		print '\n'
##      		old_ratio=old_ratio
#      print lam_B,'\t',r"%1.4f & %1.3f & %i & %1.2f \\ \hline" % (Sig_s,lam_D,int(round(np.sum(N_p)+np.sum(N_m))),V[0])

#    print np.min(np.array(sigWCA_S)),np.max(np.array(sigWCA_S))

#    ratios = np.array(sigHS_S)/np.array(sigWCA_S)
#    print ratios
#    print np.mean(ratios)
#    print np.min(ratios)
#    print np.max(ratios)
#    print np.mean(ratios)/np.min(ratios),np.mean(ratios)/np.max(ratios)




    
    characteristic_lengthS,xlabel=LDs,r'$z/\lambda_{D}$'


#    print 'Phi Bulks = ',PhiBulkS
#    print 'Mean Phi = ',np.mean(PhiBulkS)

    markers=['v','o','>','s','^','d','<','*','3','D','4','p','h','1','H','2']*15
    colors=[]
    for step in range(0,len(filenameS)):
      colors.append(ROYGBIV_map(float(step),len(filenameS)))
#    colors.reverse()

#	##Turn this off for Figure 2 and 4, I think
#    colors=[]
#    test=[]
#    for (Xi,lam_D,Bjerrum,filename) in zip(XiS,LDs,BjerrumS,filenameS):
#    		dielectric=(4*np.pi*Bjerrum)**-1
##    		print filename,round(Xi/(2.*np.pi*Bjerrum**2)/(dielectric*1.0/(1.0*lam_D)),1),lam_D
##    		if Bjerrum==20.0:
##    	    		test.append([round(Xi/(2.*np.pi*Bjerrum**2)/(dielectric*1.0/(1.0*lam_D)),1),4,1])
##    		else:
#    	    	test.append([round(Xi/(2.*np.pi*Bjerrum**2)/(dielectric*1.0/(1.0*lam_D)),1),17,1])
#    for t in test:
#    	  print t
#          colors.append(ROYGBIV_map(t[0],t[1],t[2]))

#    ii=0
#    print 'PhiS = ',PhiBulkS[ii*17:(ii+1)*17],'\n'
#    ii+=1
#    print 'PhiS = ',PhiBulkS[ii*17:(ii+1)*17],'\n'

    Debugs = 'No'
    if Debugs == 'Yes':
		##Computing ND charge densities
	    print 'CHEMICAL POTENTIAL DEBUG, NOT COUL'
	    fig=plt.figure()
	    fig.set_size_inches(3.37,4.5)
	    ax1=fig.add_subplot(411)
	    ax2=fig.add_subplot(412)
	    ax3=fig.add_subplot(413)
	    ax4=fig.add_subplot(414)

	    NfS=[]
	    NDSigmaS=[]
	    SigEffS=[]
	#    V_corS=[]
	#    Max_effS=[]
	    GCs=[]
	    BikS=[]
	    S_CS=[]
	    i=-1
	    Nf_for_EffS=[]
	    VEffS=[]
	    NtS=[]
	    SCS0s=[]
	    zEffS=[]
	    one_time_iterator=0
	    z_forS=[]
	    rhoplusS,rhominusS=[],[]
	#    testlcor =   [0.47773365833300002, 0.36662254722199999, 0.25551143611100002, 0.25551143611100002, 0.25551143611100002, 0.144400325, 0.144400325, 0.144400325, 0.0332892138889, 0.0332892138889, 0.0332892138889]
	    for (Nm,Np,n0,Sigma_s,area,lam_D,L_bin,z_positions,L_z,Bjerrum,filename,sigWCA,V,muex,Npeff,Nmeff,Phi) in zip(originalN_plusS,originalN_minusS,n0s,SIGMAS,areaS,LDs,L_binS,z_originalS,L_zS,BjerrumS,filenameS,sigWCA_S,Volts,muexEV_S,N_plusS,N_minusS,PhiWCAtot_S):
			i+=1
			dielectric=(4*np.pi*Bjerrum)**-1
			
			Nf=[(npl-nm)/(area*L_bin*dielectric*temperature/(valency*lam_D**2)) for (npl,nm) in zip(Np,Nm)]  #This is consistent with other derivatoins

	#		print 'May need to reactiviate this code!!!'
	#		if '_GC_' in filename:
	#			muex=np.array(muex[:-1])
	#		else:
	#			muex=np.array(muex)

			Integrand = 0.
			correction=0.
			Nf=np.array([(npl-nm)/(area*L_bin) for (npl,nm) in zip(Np,Nm)])
			Nf=Nf[:len(Nf)/2]
			z_restore = np.array([(x+y)/(L_z*2.) for (x,y) in zip(z_positions[0:len(z_positions)-1],z_positions[1:len(z_positions)])]) #SO, this was the original z_density required as input to the code below
			z_pos = z_restore[:len(z_restore)/2]
			for (y1,y2,z) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos):
				Integrand=Integrand + 0.5*(y1+y2)*L_bin	  
				##The code below allows one to integrate out a specific region, i.e. l_corr could be used for this
	#	        	if z*L_z<z_wall:
	#	        	       	Integrand=Integrand + 0.5*(y1+y2)*L_bin
	#	        	       	correction=correction+Integrand
	#			else:
	#			    	Integrand=Integrand + 0.5*(y1+y2)*L_bin	  
			Sigma_meas = Integrand  #This must be my reportable Sigma, but NDd
			eff=(-Sigma_meas+correction)/(dielectric*temperature/(valency*lam_D))
			ND_Sigma = abs(Sigma_meas)/(dielectric*temperature/(valency*lam_D))
			NDSigmaS.append(ND_Sigma)

	#	        Integrand = eff*(dielectric*temperature/(valency*lam_D))
	#	        SigEff=[]
	#	        VEff=[]
	#	    	volt_correction = Sigma_s/(dielectric*temperature/(valency*lam_D))  - abs(eff)
	#	        Nf_for_Eff=[]
	#	        zEff = []
	#	        rhoplus,rhominus,EVexcess = [],[],[]
	#	        for (y1,y2,z,volt,p,m,x) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos,V[:len(Nf)],Npeff[:len(Nf)],Nmeff[:len(Nf)],muex[:len(Nf)]):  #This is the original line of code -- 07/30/13 16:51:36 
	#	        		zEff.append(z)
	#	        		SigEff.append(Integrand) #There's something wrong with the first value assigned to this function
	#	        		test = volt-volt_correction*z
	#	        		VEff.append(test)
	#	        		rhoplus.append(p)
	#	        		rhominus.append(m)
	#	        		EVexcess.append(x)
	#	        		Integrand=Integrand + 0.5*(y1+y2)*L_bin	
	#	        		if filename in GC_LB_LD_5_20:
	#	       				Nf_for_Eff.append(0.5*(y1+y2))	        		    		
	#	        		else:
	#	       				Nf_for_Eff.append(y1)
	#		zEffS.append(zEff)
	#        	SigEffS.append(np.array(SigEff)/(dielectric*temperature/(valency*lam_D)))
	#        	VEffS.append(VEff)
	#        	Nf_for_EffS.append(np.array(Nf_for_Eff)/(dielectric*temperature/(valency*lam_D**2)))
	#        	NfS.append(Nf/(dielectric*temperature/(valency*lam_D**2)))

	##		rhoplus = np.array(rhoplus)/(area*L_bin*n0)
	##		rhominus = np.array(rhominus)/(area*L_bin*n0)

	##		rhoplusS.append(rhoplus)
	##		rhominusS.append(rhominus)

			if 'MC' not in filename:
				muexEV_bulk=np.mean(muex[len(muex)/2-2:len(muex)/2+2])
			else:
				muexEV_bulk=np.mean(muex[-6:-2])


			ls_x,lw_x='',0.
			ms_x = 4.0
			label_str=''
			alabel=''
			special_ms='x'
			edge_width=0.1
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
				label_str=r'$\tilde \Sigma = 7.0$'
	#			evcorS.append(3.46615)
	#			lcorS.append(0.477733658333)
				col = ROYGBIV_map(0.,10)	
			elif filename == 'Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt':
				edge_width=0.1
				special_ms='o'
				label_str=r'$5.4$'
	#			evcorS.append(2.58685576833)
	#			lcorS.append(0.366622547222)
				special_ms='*'	
				col = ROYGBIV_map(1.,10)		
			elif filename == 'Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt':
				edge_width=0.4
				special_ms='o'	
				label_str=r'$4.2$'
	#			evcorS.append(2.13433)
	#			lcorS.append(0.255511436111)
				col = ROYGBIV_map(2.,10)	
			elif filename == 'Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt':
				edge_width=0.1
				special_ms='s'	
				label_str=r'$3.6$'
	#			evcorS.append(1.78295)
	#			lcorS.append(0.255511436111)
				special_ms='d'		
				col = ROYGBIV_map(3.,10)	
			elif filename == 'Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1104_0.0106795_0.1_3.0_500.0_23.73.txt':
				edge_width=0.4
				special_ms='s'
				ms_x=3.5
				label_str=r'$3.1$'
	#			evcorS.append(1.47685)
	#			lcorS.append(0.255511436111)
				col = ROYGBIV_map(4.,10)
				Force = 0.21359	
			elif filename == 'Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt':
				edge_width=0.1
				special_ms='p'	
				ms_x=3.5	
				label_str=r'$2.3$'
	#			evcorS.append(1.01055)
	#			lcorS.append(0.144400325)
				col = ROYGBIV_map(5.,10)	
			elif filename == 'Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt':
				edge_width=0.4
				special_ms='p'	
				ms_x=3.5
				label_str=r'$2.0$'
	#			evcorS.append(0.833384464765)
	#			lcorS.append(0.144400325)
				special_ms='v'	
				col = ROYGBIV_map(6.,10)	
			elif filename == 'Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt':
				edge_width=0.1
				special_ms='H'
				ms_x=3.5
				label_str=r'$1.6$'
	#			evcorS.append(0.696928982207)
	#			lcorS.append(0.144400325)
				col = ROYGBIV_map(7.,10)	
			elif filename == 'Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt':
				edge_width=0.4
				special_ms='H'	
				label_str=r'$1.0$'
	#			evcorS.append(0.395260565885)
	#			lcorS.append(0.0332892138889)
				special_ms='>'	
				col = ROYGBIV_map(8.,10)	
				Force = 0.06948			
			elif filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1004_0.001684_0.1_0.3_500.0_23.73.txt':
				edge_width=0.1
				special_ms='D'	
				label_str=r'$0.5$'
				Force = 0.03368
	#			evcorS.append(0.316165923552)
	#			lcorS.append(0.0332892138889)
				col = ROYGBIV_map(9.,10)
			elif filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt':
				edge_width=0.4
				special_ms='<'
				label_str=r'$0.0$'
				col = 'k'
			elif filename == 'Analyzed_GC_1004_0.000842_0.1_3.0_500.0_23.73.txt' or filename =='Analyzed_GC_1004_0.000842_0.1_0.3_500.0_23.73.txt':
				edge_width=0.1
				special_ms='D'	
				label_str=r'$0.5$'
				Force = 0.01684
				col = ROYGBIV_map(2.,10)
			elif filename == 'Analyzed_GC_1004_0.000421_0.1_3.0_500.0_23.73.txt':
				edge_width=0.1
				special_ms='D'	
				label_str=r'$0.5$'
				Force = 0.00842
	#			Force =0.01684
	#			print "Have I been assigning the wrong force?"
				col = ROYGBIV_map(5.,10)
			elif filename == 'Analyzed_GC_1004_0.0_0.1_3.0_500.0_23.73.txt':
				edge_width=0.1
				special_ms='D'	
				label_str=r'$0.5$'
				Force = 0.0
				col = ROYGBIV_map(1.,10)
			else:
				col = ROYGBIV_map(abs(Sigma_s),np.max(SIGMAS)*1.1,1)
			ms_x = 2.5

			if filename == 'Analyzed_MS_1004_0.03368_0.1_3.0_502.0_23.73.txt':
				print 'SHOULD BE HERE?'
				edge_width=0.1
				special_ms='D'	
				label_str=r'$0.5$'
				Force = 0.03368
	#			evcorS.append(0.316165923552)
	#			lcorS.append(0.0332892138889)
				col = ROYGBIV_map(9.,10)
			elif filename == 'Analyzed_MS_1004_0.01684_0.1_3.0_502.0_23.73.txt':
				edge_width=0.1
				special_ms='D'	
				label_str=r'$0.5$'
				Force = 0.01684
				col = ROYGBIV_map(2.,10)


			##TOp
			if 'MC' in filename:
				Force = ''
				for x in filename[17:]:
					if x is not '_':
						Force+=x
					else:
						break
				Force = -float(Force)
				col = ROYGBIV_map(float(i),7.)
				special_ms='D'	
				edge_width=0.1
			
			muex=np.array(muex[:-1])
	#		print muex
			rhominus = np.array(Nmeff)

			rhotot = np.array(Npeff) + rhominus
			Phitot = (rhotot*(np.pi/6)*sigWCA**3) / (area*L_bin)
			PhiCS = np.array(Phitot)
			if 'MC' not in filename:
				PhiCS = PhiCS*0.954028946**3
			CS_ex = PhiCS*(8-9*PhiCS+3*PhiCS**2)/(1-PhiCS)**3

			zEff = np.array(z_positions[:-1])
			#Now only take half the data
	#		if 'MC' not in filename:
	#			muex = muex[:len(muex)/2] - muexEV_bulk
	#			rhominus = rhominus[:len(rhominus)/2]
	#			zEff = zEff[:len(zEff)/2]
	#		else:
	#			print 'muex = ',muex
	#			print 'bulk ex = ',muexEV_bulk		
	#			muex = muex - muexEV_bulk
	#			rhominus = rhominus
	#			zEff = zEff

			if 'MC' not in filename:
				Force = -Force
				print 'Changed force sign'
			muideal = np.log(rhominus) #Is rhominus correct?
			mufield = -Force *zEff
			muEV = np.array(muex)
			mutot = muideal + mufield + muEV


			print filename,Force,np.max(muideal)

	#        	print 'zeff = ',zEff
	#        	print '\n\nrho minus = ',rhominus
	#        	print 'total = ',mutot
	#        	print '\n\nmuex  = ',muEV


		       	ax1.errorbar(zEff/sigWCA,mutot,yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x)	
			ax2.errorbar(zEff/sigWCA,muideal,yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x)
			ax3.errorbar(zEff/sigWCA,mufield,yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x)
			ax4.errorbar(zEff/sigWCA,muEV,yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x)
			ax4.errorbar(zEff/sigWCA,CS_ex,yerr=None,color=col,ls='-',lw=1.0)
			
	#    print 'May need to corrrect axes about subtracting BULK'
	    ax1.set_ylabel(r'$\ln N_- - Fz + \mu^{\rm EV} $',fontsize=10.)
	    ax2.set_ylabel(r'$\ln N_- $',fontsize=10.)
	    ax3.set_ylabel(r'$-Fz $',fontsize=10.)
	    ax4.set_ylabel(r'$\mu^{\rm EV}$',fontsize=10.)
	    ax4.set_xlabel(r'$z / \sigma$',fontsize=10.)
	    
	    plt.setp(ax1.get_xticklabels(), fontsize=8.)
	    plt.setp(ax1.get_yticklabels(), fontsize=8.)
	    plt.setp(ax2.get_xticklabels(), fontsize=8.)
	    plt.setp(ax2.get_yticklabels(), fontsize=8.)
	    plt.setp(ax3.get_xticklabels(), fontsize=8.)
	    plt.setp(ax3.get_yticklabels(), fontsize=8.)
	    plt.setp(ax4.get_xticklabels(), fontsize=8.)
	    plt.setp(ax4.get_yticklabels(), fontsize=8.)

	#    ax1.set_ylim(-0.5,1.) 
	    plt.savefig('DEBUG_totalchempotential_NotCoul_NoLimits.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      

	    ax1.set_xlim(-0.25,80)
	#    ax1.set_ylim(-0.5,1.) 
	    ax2.set_xlim(-0.25,80)
	#    ax2.set_ylim(2,8) 
	    ax3.set_xlim(-0.25,80)
	#    ax3.set_ylim(-10,0.) 
	    ax4.set_xlim(-0.25,80)
	    ax4.set_ylim(0,5) 

	    plt.savefig('DEBUG_totalchempotential_NotCoul_Limits.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      

	    ax1.set_xlim(-0.25,30)
	    ax1.set_ylim(0.,8) 
	    ax2.set_xlim(-0.25,30)
	    ax2.set_ylim(0,3.5) 
	    ax3.set_xlim(-0.25,30)
	    ax3.set_ylim(0,6.) 
	    ax4.set_xlim(-0.25,30)
	    ax4.set_ylim(0,5) 

	    plt.savefig('DEBUG_totalchempotential_NotCoul_Zoom.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      

	##    ax1.legend(loc='best',numpoints=1,prop={"size":10},columnspacing=0.07,borderpad=0.18,labelspacing=0.10,handletextpad=0.07) 
	#    fig.set_size_inches(3.37,3.5)
	#    plt.savefig('F?_totalchempotential_z_BikSucks.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
	##    plt.show()
	#    plt.close()

    if len(filenameS) == 10: #10:
	    print 'Plotting chemical potential components...'
	    CarnahanStarling='yes'
	    CarnahanStarling='no'
	    fig=plt.figure()
	    fig.set_size_inches(3.37,6.5)
	    ax1=fig.add_subplot(411)
	    ax2=fig.add_subplot(412)
	    ax3=fig.add_subplot(413)
	    ax4=fig.add_subplot(414)
	    NfS=[]
	    NDSigmaS=[]
	    SigEffS=[]
	    PhiWCABulkS=[]
	    GCs=[]
	    BikS=[]
	    S_CS=[]
	    i=-1
	    Nf_for_EffS=[]
	    VEffS=[]
	    NtS=[]
	    SCS0s=[]
	    zEffS=[]
	    one_time_iterator=0
	    z_forS=[]
	    rhoplusS,rhominusS=[],[]
	    ccc=-1
	    ax4.errorbar([100,100],[100,100],yerr=None,color='k',marker='None',ls='-',label=r'${\rm CS}$')
	    for (Nm,Np,n0,Sigma_s,area,lam_D,L_bin,z_positions,L_z,Bjerrum,filename,sigWCA,V,muex,Npeff,Nmeff,Phi) in zip(originalN_plusS,originalN_minusS,n0s,SIGMAS,areaS,LDs,L_binS,z_originalS,L_zS,BjerrumS,filenameS,sigWCA_S,Volts,muexEV_S,N_plusS,N_minusS,PhiWCAtot_S):
			i+=1
			dielectric=(4*np.pi*Bjerrum)**-1
	
			Integrand = 0.
			correction=0.
			Nf=np.array([(npl-nm)/(area*L_bin) for (npl,nm) in zip(Np,Nm)])
			Nf=Nf[:len(Nf)/2]
			z_restore = np.array([(x+y)/(L_z*2.) for (x,y) in zip(z_positions[0:len(z_positions)-1],z_positions[1:len(z_positions)])]) #SO, this was the original z_density required as input to the code below
			z_pos = z_restore[:len(z_restore)/2]
			for (y1,y2,z) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos):
				Integrand=Integrand + 0.5*(y1+y2)*L_bin	  
			Sigma_meas = Integrand  #This must be my reportable Sigma, but NDd
			eff=(-Sigma_meas+correction)/(dielectric*temperature/(valency*lam_D))
			ND_Sigma = abs(Sigma_meas)/(dielectric*temperature/(valency*lam_D))
			NDSigmaS.append(ND_Sigma)

			Integrand = eff*(dielectric*temperature/(valency*lam_D))
			SigEff=[]
			VEff=[]
		    	volt_correction = Sigma_s/(dielectric*temperature/(valency*lam_D))  - abs(eff)
			Nf_for_Eff=[]
			zEff = []
			rhoplus,rhominus,EVexcess = [],[],[]
			for (y1,y2,z,volt,p,m,x) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos,V[:len(Nf)],Npeff[:len(Nf)],Nmeff[:len(Nf)],muex[:len(Nf)]):  #This is the original line of code -- 07/30/13 16:51:36 
					zEff.append(z)
					SigEff.append(Integrand) #There's something wrong with the first value assigned to this function
					test = volt-volt_correction*z
					VEff.append(test)
					rhoplus.append(p)
					rhominus.append(m)
					EVexcess.append(x)
					Integrand=Integrand + 0.5*(y1+y2)*L_bin	
					if filename in GC_LB_LD_5_20:
		       				Nf_for_Eff.append(0.5*(y1+y2))	        		    		
					else:
		       				Nf_for_Eff.append(y1)

			rhoplus = np.array(rhoplus)/(area*L_bin*n0)
			rhominus = np.array(rhominus)/(area*L_bin*n0)
			EVexcess = np.array(EVexcess) - np.mean(EVexcess[-5:])
			zEff = np.array(zEff)*L_z

			ls_x,lw_x='',0.
			ms_x = 4.0
			label_str=''
			alabel=''
			special_ms='x'
			edge_width=0.1
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
			elif 'Analyzed_GC_1312' in filename:
				edge_width=0.4
				special_ms='h'
	#			label_str=r'$\tilde \Sigma = 7.0$'
				col = ROYGBIV_map(0.,10)
	#			ze = ['6.72','4.66']	
				if sigWCA==3.0:
					label_str=r'$\tilde \Sigma = 7.1$'
					ze ='6.67'
				elif sigWCA==4.5:
					ze ='10.8'
					label_str=r'$\tilde \Sigma = 6.8$'
				elif sigWCA==6.1:
					ze ='16.23'
					label_str=r'$\tilde \Sigma = 6.4$'
			elif 'Analyzed_GC_1232' in filename:
				edge_width=0.1
				special_ms='o'
				if sigWCA==3.0:
					label_str=r'$5.4$'
					ze ='5.14'
				elif sigWCA==4.5:
					ze ='7.6'
					label_str=r'$5.3$'
				elif sigWCA==6.1:
					ze ='10.93'
					label_str=r'$5.1$'
				special_ms='*'	
				col = ROYGBIV_map(1.,10)		
			elif 'Analyzed_GC_1162' in filename:
				edge_width=0.4
				special_ms='o'	
				col = ROYGBIV_map(2.,10)	
	#			ze = ['4.06','3.28']
	#			ze = ['4.02','3.28']
				if sigWCA == 3.0:
					label_str=r'$4.2$'
					ze = '4.02'				
				elif sigWCA==4.5:
					ze ='5.41'
					label_str=r'$4.1$'
				elif sigWCA==6.1:
					ze ='7.08'
					label_str=r'$3.9$'
			elif 'Analyzed_GC_1136' in filename:
				edge_width=0.1
				special_ms='s'	
	#			label_str=r'$3.6$'
				special_ms='d'		
				col = ROYGBIV_map(3.,10)
	#			ze = ['3.54','2.86']	
				if sigWCA == 3.0:
					label_str=r'$3.6$'
					ze = '3.51'				
				elif sigWCA==4.5:
					ze ='4.51'
					label_str=r'$3.6$'
				elif sigWCA==6.1:
					ze ='5.59'
					label_str=r'$3.5$'
			elif 'Analyzed_GC_1104' in filename:
				edge_width=0.4
				special_ms='s'
				ms_x=3.5
				label_str=r'$3.1$'
				col = ROYGBIV_map(4.,10)
				Force = 0.21359	
	#			ze = ['3.08','2.48']
	#			ze = ['3.11','2.48']
				if sigWCA==3.0:
					ze = '3.11'
					label_str=r'$3.1$'
				elif sigWCA==4.5:
					ze ='3.81'
				elif sigWCA==6.1:
					ze = '4.599'
			elif 'Analyzed_GC_1066' in filename:
				edge_width=0.1
				special_ms='p'	
				ms_x=3.5	
				label_str=r'$2.3$'
				col = ROYGBIV_map(5.,10)
	#			ze = ['2.32','2.01']	
				if sigWCA == 3.0:
					ze = '2.30'
				elif sigWCA==4.5:
					ze ='2.62'
				elif sigWCA==6.1:
					ze = '2.98'
			elif 'Analyzed_GC_1048' in filename:
				edge_width=0.4
				special_ms='p'	
				ms_x=3.5
				label_str=r'$2.0$'
				special_ms='v'	
				col = ROYGBIV_map(6.,10)
	#			ze = ['2.01','1.76']	
				if sigWCA == 3.0:
					ze = '1.97'
				elif sigWCA==4.5:
					ze ='2.13'
				elif sigWCA==6.1:
					ze = '2.31'
			elif 'Analyzed_GC_1034' in filename:
				edge_width=0.1
				special_ms='H'
				ms_x=3.5
				label_str=r'$1.6$'
				col = ROYGBIV_map(7.,10)
	#			ze = ['1.65','1.44']	
				if sigWCA==3.0:
					ze = '1.65'
				elif sigWCA==4.5:
					ze ='1.7'
				elif sigWCA==6.1:
					ze = '1.74'
			elif 'Analyzed_GC_1016' in filename:
				edge_width=0.4
				special_ms='H'	
				label_str=r'$1.0$'
				special_ms='>'	
				col = ROYGBIV_map(8.,10)	
				Force = 0.06948		
				if sigWCA == 3.0:
					ze = '1.07'
				elif sigWCA==4.5:
					ze ='1.05'
				elif sigWCA==6.1:
					ze = '1.09'	
			elif 'Analyzed_GC_1004' in filename:
				edge_width=0.1
				special_ms='D'	
				label_str=r'$0.5$'
				Force = 0.03368
				col = ROYGBIV_map(9.,10)
	#			ze = ['0.50','0.47']
				if sigWCA == 3.0:
					ze = '0.54'
				elif sigWCA==4.5:
					ze ='0.4701'
				elif sigWCA==6.1:
					ze ='0.49'
					label_str=r'$0.49$'
			elif filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1000_0.0_0.1_4.5_500.0_23.73.txt':
				edge_width=0.4
				special_ms='<'
				label_str=r'$0.0$'
				col = 'k'
				ze = ['0.00','0.00']
			else:
	#			print '\n\n',filename,Sigma_s,np.max(SIGMAS)
				col = ROYGBIV_map(abs(Sigma_s),np.max(SIGMAS)*1.1,1)
			ms_x = 2.5
			edge_width=0.1

			rhotot = np.array(rhominus) + np.array(rhominus)			
			Phitot = (rhotot*(np.pi/6)*sigWCA**3) / (area*L_bin)
			PhiCS = np.array(Phitot)

			PhiBulk = np.mean(Phitot[-5:])
			PhiWCABulkS.append(PhiBulk)

	#		print 'file,bulk vol frac,zeta,Sigma = ',filename,np.mean(Phitot[-5:])*0.736419046,VEff[0],eff
	#		print eff
			if sigWCA==3.0 or sigWCA==4.5 or sigWCA==6.1:
				CarnahanStarling='yes'
			else:
				CarnahanStarling='no'

			if sigWCA==3.0:
				special_ms = 'v'
			elif sigWCA == 3.77976:
				special_ms = 'h'
			elif sigWCA ==4.5:
				special_ms = 'd'
			elif sigWCA == 4.7622:
				special_ms = 's'
			elif sigWCA == 5.12993:
				special_ms = '<'
			elif sigWCA == 5.56991:
				special_ms = 'p'
			elif sigWCA == 5.84609:
				special_ms = 'o'			
			elif sigWCA==6.1:
				special_ms = '^'
				
			if CarnahanStarling=='yes' and ze!='0.00':	
				  if sigWCA==3.0:
				  	CS_file='CS_0.037_' + ze + '.txt'
			  	  elif sigWCA==4.5:
			  	  	CS_file='CS_0.125_' + ze + '.txt'
			  	  elif sigWCA==6.1:
			  	  	CS_file='CS_0.311_' + ze + '.txt'
				  MMA_CS=[[float(x) for x in line.split()] for line in file(CS_file,"r").readlines()]
				  MMAz,PhiWCA = [],[]
				  MMA_co,MMA_c,MMA_volt = [],[],[]
				  for line in MMA_CS:
				  	MMAz.append(line[0])
				  	Phi = (line[1]+line[2])*n0*(np.pi/6)*sigWCA**3
				  	PhiWCA.append(Phi)
				  	MMA_co.append(line[1])
				  	MMA_c.append(line[2])
				  	MMA_volt.append(line[4])
		  	  	  PhiCS = np.array(PhiWCA)*0.903042806**3
				  CS_ex = PhiCS*(8-9*PhiCS+3*PhiCS**2)/(1-PhiCS)**3
				  MMAz = np.array(MMAz)*lam_D/sigWCA-(1+0.299602925)/sigWCA
				  MMA_ex = CS_ex - CS_ex[len(CS_ex)-1] 

			zPM = zEff/sigWCA-(1+0.299602925)/sigWCA
		       	if CarnahanStarling=='yes':
		       		ax3.errorbar(MMAz,MMA_volt,yerr=None,color=col,ls='-',lw=1.0)	
		       		ax1.errorbar(MMAz,np.log(MMA_co),yerr=None,color=col,ls='-',lw=1.0)
		       		ax2.errorbar(MMAz,np.log(MMA_c),yerr=None,color=col,ls='-',lw=1.0)
		       		ax4.errorbar(MMAz,MMA_ex,yerr=None,color=col,ls='-',lw=1.0)	       		

		    		t1 = np.linspace(MMAz[0],MMAz[-1],20)[1:-1]
		    		interp_ExFromZ = interpolate.LSQUnivariateSpline(MMAz,MMA_ex,t1,k=5)
		    		fit_ExFromZ = lambda zpos: interp_ExFromZ(zpos)[0]  

			ax1.errorbar(zEff/sigWCA-(1+0.299602925)/sigWCA,np.log(rhoplus),yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)
			ax2.errorbar(zEff/sigWCA-(1+0.299602925)/sigWCA,np.log(rhominus),yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)
		       	ax3.errorbar(zEff/sigWCA-(1+0.299602925)/sigWCA,VEff,yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)
			ax4.errorbar(zEff/sigWCA-(1+0.299602925)/sigWCA,EVexcess,yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)


	#		print '\n\n',filename
  		
	    ax3.set_ylabel(r'$(\phi - \phi_{\rm B})/ \phi_{\rm T}$',fontsize=10.)
	    ax1.set_ylabel(r'$\ln (n^{+}/n_{\rm B})$',fontsize=10.)
	    ax2.set_ylabel(r'$\ln (n^{-}/n_{\rm B})$',fontsize=10.)
	    ax4.set_ylabel(r'$(\mu^{\rm EV}_{\rm ex}-\mu^{\rm EV}_{\rm ex,B}) / k_{\rm B}T$',fontsize=9.)
	    ax4.set_xlabel(r'$(z-\delta_{\rm w}) / \sigma$',fontsize=10.)
	    
	    plt.setp(ax1.get_xticklabels(), fontsize=8.)
	    plt.setp(ax1.get_yticklabels(), fontsize=8.)
	    plt.setp(ax2.get_xticklabels(), fontsize=8.)
	    plt.setp(ax2.get_yticklabels(), fontsize=8.)
	    plt.setp(ax3.get_xticklabels(), fontsize=8.)
	    plt.setp(ax3.get_yticklabels(), fontsize=8.)
	    plt.setp(ax4.get_xticklabels(), fontsize=8.)
	    plt.setp(ax4.get_yticklabels(), fontsize=8.)

	    ax1.set_xlim(-0.5,10)
	    ax2.set_xlim(-0.5,10)
    	    ax3.set_xlim(-0.5,10)	
    	    ax4.set_xlim(-0.5,10)

    	    if sigWCA==4.5:   
	    	ax1.set_ylim(-15,1.0) 
	    	ax2.set_ylim(-1.0,2.5) 
	    	ax3.set_ylim(-0.5,12.) 
	    	ax4.set_ylim(-1.25,12.) 
    	    elif sigWCA==6.1:   
	    	ax1.set_ylim(-15,1.0) 
	    	ax2.set_ylim(-1.0,2.5) 
	    	ax3.set_ylim(-0.5,17) 
	    	ax4.set_ylim(-1.25,50.) 

	    ax2.yaxis.set_major_locator(MaxNLocator(5)) 

	#    ax4.set_ylim(-0.25,2.5) 
	
	    if sigWCA==3.0 or sigWCA==4.5 or sigWCA==6.1:
		    ax4.legend(loc='upper right',numpoints=1,prop={"size":8},columnspacing=0.05,borderpad=0.17,labelspacing=0.10,handletextpad=0.15,handlelength = 1.0 ,ncol=3,markerscale=2.0) 
	    if len(filenameS)==10:   
	    	print 'Not saving  VoltRhosEV... plot'
	    	plt.savefig('VoltRhosEV_'+str(np.round(np.mean(PhiWCABulkS),3))+'.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
	    plt.close()

#	    print fail

    print 'Plotting f1inding oscillations...'
    CarnahanStarling='yes'
#    CarnahanStarling='no'
    fig=plt.figure()
    fig.set_size_inches(3.37,3.5)
    ax1=fig.add_subplot(111)
#    ax2=fig.add_subplot(212)
    NfS=[]
    NDSigmaS=[]
    SigEffS=[]
    PhiWCABulkS = []
    GCs=[]
    BikS=[]
    S_CS=[]
    i=-1
    Nf_for_EffS=[]
    VEffS=[]
    NtS=[]
#    SCS0s=[]
    zEffS=[]
    one_time_iterator=0
    z_forS=[]
    rhoplusS,rhominusS=[],[]
    ccc=-1
    offset = np.array(range(len(filenameS)))*0.75
    offset = offset.tolist()
    offset.reverse()
    lcorS = []
    markerS=[]
    PhiLocalCorS=[]
    VoltCorS,SigCorS,SigDifS=[],[],[]
    rhofCorS,rhototCorS = [],[]
    zetaS=[]
    SigRefS=[]
    EvexCorS=[]
    for (Nm,Np,n0,Sigma_s,area,lam_D,L_bin,z_positions,L_z,Bjerrum,filename,sigWCA,V,muex,Npeff,Nmeff,Phi) in zip(originalN_plusS,originalN_minusS,n0s,SIGMAS,areaS,LDs,L_binS,z_originalS,L_zS,BjerrumS,filenameS,sigWCA_S,Volts,muexEV_S,N_plusS,N_minusS,PhiWCAtot_S):
		i+=1
		dielectric=(4*np.pi*Bjerrum)**-1
		
		Integrand = 0.
		correction=0.
		Nf=np.array([(npl-nm)/(area*L_bin) for (npl,nm) in zip(Np,Nm)])
		Nf=Nf[:len(Nf)/2]
		z_restore = np.array([(x+y)/(L_z*2.) for (x,y) in zip(z_positions[0:len(z_positions)-1],z_positions[1:len(z_positions)])]) #SO, this was the original z_density required as input to the code below
		z_pos = z_restore[:len(z_restore)/2]
		for (y1,y2,z) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos):
			Integrand=Integrand + 0.5*(y1+y2)*L_bin	  
		Sigma_meas = Integrand  #This must be my reportable Sigma, but NDd
		eff=(-Sigma_meas+correction)/(dielectric*temperature/(valency*lam_D))
		ND_Sigma = abs(Sigma_meas)/(dielectric*temperature/(valency*lam_D))
		NDSigmaS.append(ND_Sigma)
		SigRefS.append((dielectric*temperature/(valency*lam_D)))

	        Integrand = eff*(dielectric*temperature/(valency*lam_D))
	        SigEff=[]
	        VEff=[]
	    	volt_correction = Sigma_s/(dielectric*temperature/(valency*lam_D))  - abs(eff)
	        Nf_for_Eff=[]
	        zEff = []
	        rhoplus,rhominus,EVexcess = [],[],[]
	        for (y1,y2,z,volt,p,m,x) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos,V[:len(Nf)],Npeff[:len(Nf)],Nmeff[:len(Nf)],muex[:len(Nf)]):  #This is the original line of code -- 07/30/13 16:51:36 
	        		zEff.append(z)
	        		SigEff.append(Integrand) #There's something wrong with the first value assigned to this function
	        		test = volt-volt_correction*z
	        		VEff.append(test)
	        		rhoplus.append(p)
	        		rhominus.append(m)
	        		EVexcess.append(x)
	        		Integrand=Integrand + 0.5*(y1+y2)*L_bin	
	        		if filename in GC_LB_LD_5_20:
	       				Nf_for_Eff.append(0.5*(y1+y2))	        		    		
	        		else:
	       				Nf_for_Eff.append(y1)

		rhoplus = np.array(rhoplus)/(area*L_bin*n0)
		rhominus = np.array(rhominus)/(area*L_bin*n0)
		zEff = np.array(zEff)*L_z
		zetaS.append(VEff[0])	
		SigEff = np.array(SigEff)/(dielectric*temperature/(valency*lam_D))
		
		ls_x,lw_x='',0.
		ms_x = 4.0
		label_str=''
		alabel=''
		special_ms='x'
		edge_width=0.1
#		print 'Zetas might need to change when final measurements come out'
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
		elif 'Analyzed_GC_1312' in filename:
			edge_width=0.4
			special_ms='h'
#			label_str=r'$\tilde \Sigma = 7.0$'
			col = ROYGBIV_map(0.,10)
#			ze = ['6.72','4.66']	
			if sigWCA==3.0:
				label_str=r'$\tilde \Sigma = 7.1$'
				ze ='6.67'
			elif sigWCA==4.5:
				ze ='10.8'
				label_str=r'$\tilde \Sigma = 6.8$'
			elif sigWCA==6.1:
				ze ='16.23'
				label_str=r'$\tilde \Sigma = 6.4$'
		elif 'Analyzed_GC_1232' in filename:
			edge_width=0.1
			special_ms='o'
			if sigWCA==3.0:
				label_str=r'$5.4$'
				ze ='5.14'
			elif sigWCA==4.5:
				ze ='7.6'
				label_str=r'$5.3$'
			elif sigWCA==6.1:
				ze ='10.93'
				label_str=r'$5.1$'
			special_ms='*'	
			col = ROYGBIV_map(1.,10)		
		elif 'Analyzed_GC_1162' in filename:
			edge_width=0.4
			special_ms='o'	
			col = ROYGBIV_map(2.,10)	
#			ze = ['4.06','3.28']
#			ze = ['4.02','3.28']
			if sigWCA == 3.0:
				label_str=r'$4.2$'
				ze = '4.02'				
			elif sigWCA==4.5:
				ze ='5.41'
				label_str=r'$4.1$'
			elif sigWCA==6.1:
				ze ='7.08'
				label_str=r'$3.9$'
		elif 'Analyzed_GC_1136' in filename:
			edge_width=0.1
			special_ms='s'	
#			label_str=r'$3.6$'
			special_ms='d'		
			col = ROYGBIV_map(3.,10)
#			ze = ['3.54','2.86']	
			if sigWCA == 3.0:
				label_str=r'$3.6$'
				ze = '3.51'				
			elif sigWCA==4.5:
				ze ='4.51'
				label_str=r'$3.6$'
			elif sigWCA==6.1:
				ze ='5.59'
				label_str=r'$3.5$'
		elif 'Analyzed_GC_1104' in filename:
			edge_width=0.4
			special_ms='s'
			ms_x=3.5
			label_str=r'$3.1$'
			col = ROYGBIV_map(4.,10)
			Force = 0.21359	
#			ze = ['3.08','2.48']
#			ze = ['3.11','2.48']
			if sigWCA==3.0:
				ze = '3.11'
				label_str=r'$3.1$'
			elif sigWCA==4.5:
				ze ='3.81'
			elif sigWCA==6.1:
				ze = '4.599'
		elif 'Analyzed_GC_1066' in filename:
			edge_width=0.1
			special_ms='p'	
			ms_x=3.5	
			label_str=r'$2.3$'
			col = ROYGBIV_map(5.,10)
#			ze = ['2.32','2.01']	
			if sigWCA == 3.0:
				ze = '2.30'
			elif sigWCA==4.5:
				ze ='2.62'
			elif sigWCA==6.1:
				ze = '2.98'
		elif 'Analyzed_GC_1048' in filename:
			edge_width=0.4
			special_ms='p'	
			ms_x=3.5
			label_str=r'$2.0$'
			special_ms='v'	
			col = ROYGBIV_map(6.,10)
#			ze = ['2.01','1.76']	
			if sigWCA == 3.0:
				ze = '1.97'
			elif sigWCA==4.5:
				ze ='2.13'
			elif sigWCA==6.1:
				ze = '2.31'
		elif 'Analyzed_GC_1034' in filename:
			edge_width=0.1
			special_ms='H'
			ms_x=3.5
			label_str=r'$1.6$'
			col = ROYGBIV_map(7.,10)
#			ze = ['1.65','1.44']	
			if sigWCA==3.0:
				ze = '1.65'
			elif sigWCA==4.5:
				ze ='1.7'
			elif sigWCA==6.1:
				ze = '1.74'
		elif 'Analyzed_GC_1016' in filename:
			edge_width=0.4
			special_ms='H'	
			label_str=r'$1.0$'
			special_ms='>'	
			col = ROYGBIV_map(8.,10)	
			Force = 0.06948		
			if sigWCA == 3.0:
				ze = '1.07'
			elif sigWCA==4.5:
				ze ='1.05'
			elif sigWCA==6.1:
				ze = '1.09'	
		elif 'Analyzed_GC_1004' in filename:
			edge_width=0.1
			special_ms='D'	
			label_str=r'$0.5$'
			Force = 0.03368
			col = ROYGBIV_map(9.,10)
#			ze = ['0.50','0.47']
			if sigWCA == 3.0:
				ze = '0.54'
			elif sigWCA==4.5:
				ze ='0.4701'
			elif sigWCA==6.1:
				ze ='0.49'
				label_str=r'$0.49$'
		elif 'Analyzed_GC_1000' in filename:
			edge_width=0.4
			special_ms='<'
			label_str=r'$0.0$'
			col = 'green'
			ze = ['0.00','0.00']
		else:
#			print '\n\n',filename,Sigma_s,np.max(SIGMAS)
			col = ROYGBIV_map(abs(Sigma_s),np.max(SIGMAS)*1.1,1)
		ms_x = 2.5
		edge_width=0.1


		rhotot = np.array(rhominus) + np.array(rhoplus)
		Phitot = (rhotot*(np.pi/6)*sigWCA**3) / (area*L_bin)
		PhiCS = np.array(Phitot)

#		print 'file,bulk vol frac,zeta,Sigma = ',filename,np.mean(Phitot[-5:])*0.736419046,VEff[0],eff
#		print eff
		PhiBulk = np.mean(Phitot[-5:])
		PhiWCABulkS.append(PhiBulk)

		NDrhofree = -(np.array(rhoplus) - np.array(rhominus))/2

#		if sigWCA==3.0 or sigWCA==4.5 or sigWCA==6.1:
#			CarnahanStarling='yes'
#		else:
#			CarnahanStarling='no'

#		CarnahanStarling='no'

		if sigWCA==3.0:
			special_ms = 'v'
		elif sigWCA == 3.77976:
			special_ms = 'h'
		elif sigWCA ==4.5:
			special_ms = 'd'
		elif sigWCA == 4.7622:
			special_ms = 's'
		elif sigWCA == 5.12993:
			special_ms = '<'
		elif sigWCA == 5.56991:
			special_ms = 'p'
		elif sigWCA == 5.84609:
			special_ms = 'o'			
		elif sigWCA==6.1:
			special_ms = '^'
		markerS.append(special_ms)		

		if i==0: 
		    ax1.plot([-1,-1],[-1,-1],color = 'k', ls = '-',lw = 0.5,marker = 'None',label =r'$\tilde \rho_{\rm spline}$')               
		    ax1.plot([-1,-1],[-1,-1],color = 'k', ls = 'None', marker = special_ms, mew = 0.5, ms=ms_x*1.5,mec = 'grey',label =r'$\tilde \rho_{\rm spline}^\prime( z^{\blacktriangle} ) = 0$')               
		    ax1.plot([-1,-1],[-1,-1],color = '0.75', marker = 'None',label =r'$\tilde \rho_{\rm fit}( z^{\blacktriangle}+\sigma)$')               
		    ax1.plot([-1,-1],[-1,-1],color='white',marker='*',mew=edge_width*4.5,ms=ms_x*3.0,mec = 'k',ls= 'None',label =r'$\tilde \rho (\ell_{\rm cor})$')	
    		    ax1.plot([-1,-1],[-1,-1],color='k',ls='-',lw=1.25,label =r'$\tilde \rho_{\rm CS}$')							       	
		
		if CarnahanStarling=='yes' and ze!='0.00':	
			  if sigWCA==3.0:
			  	CS_file='CS_0.037_' + ze + '.txt'
		  	  elif sigWCA==4.5:
		  	  	CS_file='CS_0.125_' + ze + '.txt'
		  	  elif sigWCA==6.1:
		  	  	CS_file='CS_0.311_' + ze + '.txt'
			  MMA_CS=[[float(x) for x in line.split()] for line in file(CS_file,"r").readlines()]
			  MMAz = []
			  NDMMA_free,NDMMA_tot = [],[]
			  for line in MMA_CS:
			  	MMAz.append(line[0])
			  	free = (-line[1] + line[2])/2
			  	NDMMA_free.append(free)
				tot = line[1] + line[2]
				NDMMA_tot.append(tot)
			  MMAz = np.array(MMAz)*lam_D/sigWCA-(1+0.299602925)/sigWCA
	       		  ax1.errorbar(MMAz,np.array(NDMMA_free)+offset[i],yerr=None,color='k',ls='-',lw=1.25)	

		zPM = np.array(zEff)/sigWCA-(1+0.299602925)/sigWCA
#	       	ax1.errorbar(zPM,NDrhofree+offset[i],yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)

		kk=4
		if sigWCA>=6.1 or sigWCA == 5.84609: #0.42 and 0.37
		  	t1 = np.linspace(zPM[0],1,6)[1:].tolist() + np.linspace(1,7,34)[1:].tolist() + np.linspace(7,zPM[-1],25)[1:-1].tolist()
		  	kk=5
  		elif  sigWCA == 5.56991 or sigWCA == 5.12993 or sigWCA == 4.7622 or sigWCA == 4.5: #0.32 and 0.25 and 0.20 and 0.17
		  	t1 = np.linspace(zPM[0],4,15)[1:].tolist() + np.linspace(4,zPM[-1],25)[1:-1].tolist() 	 
  		elif  sigWCA == 3.77976: #0.10
		  	t1 = np.linspace(zPM[0],1,4)[1:].tolist() + np.linspace(1,2,4)[1:].tolist() + np.linspace(2,3,4)[1:].tolist() + np.linspace(3,5,4)[1:].tolist() + np.linspace(5,zPM[-1],20)[1:-1].tolist()
  		elif sigWCA == 3.0: #0.05
		  	t1 = np.linspace(zPM[0],2,5)[1:].tolist() + np.linspace(2,zPM[-1],25)[1:-1].tolist() 
#		  	t1 = np.linspace(zPM[0],1,3)[1:].tolist() + np.linspace(1,2,3)[1:].tolist() + np.linspace(2,zPM[-1],20)[1:-1].tolist()
		  	kk=5

    		free_spline = interpolate.LSQUnivariateSpline(zPM,np.array(NDrhofree),t1,k=kk)
		spline_prime = lambda z : free_spline.__call__(z,1)[0]
		single_free_spline = lambda z : free_spline(z)[0]

		PhiLocal_spline = interpolate.LSQUnivariateSpline(zPM,np.array(Phitot),t1,k=kk)
		single_PhiLocal = lambda z : PhiLocal_spline(z)[0]		

		Voltage_spline = interpolate.LSQUnivariateSpline(zPM,np.array(VEff),t1,k=kk)
		single_Voltage = lambda z : Voltage_spline(z)[0]	

		SigEff_spline = interpolate.LSQUnivariateSpline(zPM,SigEff,t1,k=kk)
		single_SigEff = lambda z : SigEff_spline(z)[0]	

		Evex_spline = interpolate.LSQUnivariateSpline(zPM,np.array(EVexcess),t1,k=kk)
		single_Evex = lambda z : Evex_spline(z)[0]	


		zPM_plot = np.linspace(zPM[0],zPM[-1],len(zPM)*2) # 2X as much spacing than the bins
	       	ax1.errorbar(zPM_plot,free_spline(zPM_plot)+offset[i],yerr=None,color=col,ls ='-',lw = 0.5)
   
		x0s = np.linspace(0.0,8,30)
		if label_str==r'$\tilde \Sigma = 6.4$':
			x1s = np.linspace(0.0,4.,20).tolist() + np.linspace(4.,6,5).tolist() + x0s[x0s>=6].tolist()
			x0s = np.array(x1s)

#		print '\n',filename,special_ms
		primeroots=[]
		for x0 in x0s:
			SolverResults = fsolve(spline_prime, x0,full_output = 1)
			if 'solution' in SolverResults[3]:
				root1 = np.round(SolverResults[0],4)	
				if (root1 not in primeroots) and root1<=9 and root1>0:
					primeroots.append(root1)
		primeroots.sort()
		primeroots = primeroots[::-1]


		if len(primeroots)!=0:
#			print 'Possible lcors = ',primeroots
			ax1.errorbar(primeroots,free_spline(primeroots)+offset[i],color='k',marker=special_ms,mec = col,mew=0.5,ms=ms_x*1.5,ls='None')#,label=label_str)#,mec = col)

		if i==0:
			lcor_set = []
			sigOld = sigWCA
		elif sigOld!=sigWCA:
#			print lcorS[-1],sigWCA
			lcor_set = []
			sigOld = sigWCA			
		
		pretest = len(lcor_set)
		jj=-1
		for (far,near) in zip(primeroots[0:len(primeroots)-1],primeroots[1:len(primeroots)]):
			if far<=0:
				far = 0.
				lcorS.append(far)
				j+=1

			spacing = far - near
		
			if jj==-1 and (spacing<=0.505 and spacing>0.15):

				rhonear = single_free_spline(near)
				rhofar = single_free_spline(far)
				rhodif = np.abs(rhofar - rhonear)
				bulkdifS = []
				for z in zPM[-90:-15]:
					bulkdif = single_free_spline(z) - single_free_spline(z+spacing)
					bulkdifS.append(bulkdif)
				bulkRMS = np.sqrt(np.mean(np.array(bulkdifS)**2))

				if rhodif > 2.*bulkRMS: #This tests for a significant fluctuation
					z_bulkwards,rho_bulkwards = [],[]
					z_wallwards,rho_wallwards = [],[]
					for (x,y) in zip(zPM,NDrhofree):
						if x>=near-1.0 and x<=far and x>0:
							z_wallwards.append(x)
							rho_wallwards.append(y)								
						elif x>far and x<=far+1.0:								
							z_bulkwards.append(x)
							rho_bulkwards.append(y)
						elif x>far+1.0:
							z_wallwards = np.array(z_wallwards)
							rho_wallwards = np.array(rho_wallwards)
							z_bulkwards = np.array(z_bulkwards)
							rho_bulkwards = np.array(rho_bulkwards)	
							break
					z_cor = np.linspace(z_wallwards[0],far,20)
					ztest = np.array(z_cor.tolist() + np.linspace(far+z_cor[1]-z_cor[0],z_bulkwards[-1],40).tolist())

					expslope, expintercept, expr_value, expp_value, expstd_err = stats.linregress(z_bulkwards,np.log(rho_bulkwards))
					exprhotest = np.exp(ztest*expslope + expintercept)
					dif_from_exp = np.array(free_spline(ztest)) - exprhotest
					MaxDeviationFromExpFit = np.max(np.abs(dif_from_exp[20:]))
					fitRMS_exp = np.sqrt(((np.exp(z_bulkwards*expslope+expintercept) - rho_bulkwards)** 2).mean())
					corRMS_exp = np.sqrt(((np.exp(z_wallwards*expslope+expintercept) - rho_wallwards)** 2).mean())
#					print 'EXPONENTIAL (expR^2,expslope,expintercept) = (',expr_value**2,',',expslope,',',expintercept,')'
#					print 'EXPONENTIAL (max, min, fit extrema) = (',max(dif_from_exp[:20]),',',min(dif_from_exp[:20]),',',MaxDeviationFromExpFit,')'
#					print 'EXPONENTIAL (fitRMS, corRMS, cor/fit ) = ',fitRMS_exp,',',corRMS_exp,',',corRMS_exp/fitRMS_exp,')'

						###Old way of doing, worth preserving
#					dif_spline = interpolate.LSQUnivariateSpline(z_cor,np.array(free_spline(z_cor)) - np.exp(z_cor*expslope + expintercept),np.linspace(z_cor[0],z_cor[-1],5)[1:-1],k=5)
#					dif_spline_prime = lambda z : dif_spline.__call__(z,1)[0]
#					x0s = np.linspace(z_cor[0],z_cor[-1],20)

					dif_spline = interpolate.LSQUnivariateSpline(ztest,np.array(free_spline(ztest)) - np.exp(ztest*expslope + expintercept),np.linspace(ztest[0],ztest[-1],5)[1:-1],k=5)
					dif_spline_prime = lambda z : dif_spline.__call__(z,1)[0]

					difroots = []
					x0s = np.linspace(z_cor[0],z_cor[-1],20)
					for x0 in x0s:
						SolverResults =  fsolve(dif_spline_prime, x0,full_output = 1)
						if 'solution' in SolverResults[3]:
							root1 = np.round(SolverResults[0],3)	
							if (root1 not in difroots) and root1<=9:
								difroots.append(root1)
					difroots.sort()
					difroots = np.array(difroots[::-1])

#					print 'Dif roots and spacing = ',difroots,difroots[:-1]-difroots[1:]
#					print 'DifValues at roots = ',dif_spline(difroots)
#					print 'Delta DifValues = ',dif_spline(difroots)[:-1]-dif_spline(difroots)[1:]
					LDAdifroots = []
					x0s = np.linspace(ztest[20],ztest[-1],20)
					for x0 in x0s:
						SolverResults =  fsolve(dif_spline_prime, x0,full_output = 1)
						if 'solution' in SolverResults[3]:
							root1 = np.round(SolverResults[0],3)	
							if (root1 not in LDAdifroots) and root1<=9:
								LDAdifroots.append(root1)
					LDAdifroots.sort()
					LDAdifroots = np.array(LDAdifroots[::-1])

#					print 'LDA Dif roots and spacing = ',LDAdifroots,LDAdifroots[:-1]-LDAdifroots[1:]
#					print 'LDA DifValues at roots = ',dif_spline(LDAdifroots)
#					print 'LDA Delta DifValues = ',dif_spline(LDAdifroots)[:-1]-dif_spline(LDAdifroots)[1:]					
					if max(dif_from_exp[:20]) > MaxDeviationFromExpFit and min(dif_from_exp[:20]) < -MaxDeviationFromExpFit:
						if corRMS_exp > 5.*fitRMS_exp:
							if len(lcor_set)==0:
								jj+=1								
								ax1.errorbar(ztest,exprhotest+offset[i],yerr=None,color='0.75',ls ='-')
							       	ax1.errorbar(far,single_free_spline(far)+offset[i],yerr=None,color='white',marker='*',mew=edge_width*4.5,ms=ms_x*3.0,mec = col)							       	

							       	
								dif_single = single_free_spline(far) - np.exp(far*expslope + expintercept)
							       	lcorS.append(far)
							       	lcor_set.append(far)	
						       	elif far <= lcor_set[-1]+0.5*0.954028946**3:
								jj+=1
								ax1.errorbar(ztest,exprhotest+offset[i],yerr=None,color='0.75',ls ='-')
#							       	ax1.errorbar(far,single_free_spline(far)+offset[i],yerr=None,color='white',marker=special_ms,mew=edge_width*1.75,ms=ms_x*2.25,mec = col)							       	

							       	ax1.errorbar(far,single_free_spline(far)+offset[i],yerr=None,color='white',marker='*',mew=edge_width*4.5,ms=ms_x*3.0,mec = col)							       	

								dif_single = single_free_spline(far) - np.exp(far*expslope + expintercept)
							       	lcorS.append(far)
								lcor_set.append(far)	
		if len(lcor_set) == pretest:
			lcorS.append(0.)
		       	lcor_set.append(0.)
#			print "WARNING NO LCOR FOUND!"

		VoltageCor = single_Voltage(lcorS[-1])

		if lcorS[-1]>0:
			for zpm in zPM[zPM<=lcorS[-1]]:
				if zpm>0:
					print zpm,'\t',zetaS[-1]-single_Voltage(zpm)


		
		VoltCorS.append(VoltageCor)
		SigCor = ND_Sigma  + single_SigEff(lcorS[-1]) # = -Sigma - SigEff(lcor) [The vars are negatives of e/o]
		corindex = len(zPM[zPM<=(lcorS[-1]-0.5*L_bin/sigWCA)])
		binweight = sigWCA*abs(lcorS[-1] - zPM[corindex])/L_bin
		NpCor = np.sum(rhoplus[:corindex])+binweight*rhoplus[corindex] #= Np/(Axy * n0 * Lbin)
		NmCor = np.sum(rhominus[:corindex])+binweight*rhominus[corindex]
		rhofCor = -L_bin*0.5*(NpCor - NmCor)/((lcorS[-1]-zPM[0])*sigWCA)
		EvexCorS.append(single_Evex(lcorS[-1]))
		SigCorS.append(SigCor*(sigWCA/lam_D))
		SigDif = ND_Sigma - SigCor   #This is probably nonsense right now
		SigDifS.append(SigDif)
		NDrhofcor = (rhofCor*n0)/(valency/(4*np.pi*Bjerrum*sigWCA**2))
#		NDrhofcor = rhofCor*n0#/(valency/sigWCA**3)
		rhofCorS.append(NDrhofcor)		
		rhototCor = L_bin*(NpCor + NmCor)/((lcorS[-1]-zPM[0])*sigWCA) #This is rhoTOT
		rhototCorS.append(rhototCor)
		PhiCor = rhototCor*n0*(np.pi/6)*sigWCA**3
		PhiLocalCorS.append(PhiCor)
#		lcorS[-1] = lcorS[-1]*sigWCA
    i=-1
    for f in filenameS:
	i+=1
#    	ax2.errorbar(zPM_plot,np.zeros(len(zPM_plot))+offset[i],yerr=None,color='grey',ls ='-',lw = 0.5)
    ax1.set_ylabel(r'$ -\rho / 2 qe n^{\rm B} + {\rm offset} $',fontsize=10.) # = (n_+ - n_-) / 2 \rho^{\rm B}
#    ax2.set_ylabel(r'$ \tilde \rho - \tilde \rho_{\rm fit}(z>\ell_{\rm cor})  + {\rm offset}$',fontsize=10.) 
    ax1.set_xlabel(r'$(z-\delta_{\rm w}) / \sigma$',fontsize=10.)
    
    plt.setp(ax1.get_xticklabels(), fontsize=8.)  
    plt.setp(ax1.get_yticklabels(), fontsize=8.)
#    plt.setp(ax2.get_xticklabels(), fontsize=8.)
#    plt.setp(ax2.get_yticklabels(), fontsize=8.)

    ax1.set_xlim(-0.5,10)
    ax1.set_ylim(-0.25,12.)
    ax1.xaxis.set_major_locator(MaxNLocator(11))  
#    ax2.set_xlim(-0.5,10)
#    ax2.set_ylim(-0.5,np.max(offset)*1.05) 

    if len(filenameS)==10:
	ax1.legend(loc='upper right',numpoints=1,prop={"size":8},columnspacing=0.07,borderpad=0.25,labelspacing=0.07,handletextpad=0.10,handlelength = 1.25,ncol=1,markerscale=2.0)     
    	plt.savefig('JustRhoFindingOscillations_'+str(np.round(np.mean(PhiWCABulkS),3))+'.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
    plt.close()

    print fail

    if len(filenameS)<=88:
   
	    print 'Plotting finding oscillations...'
	    CarnahanStarling='yes'
	    CarnahanStarling='no'
	    fig=plt.figure()
	    fig.set_size_inches(3.37,4.5)
	    ax1=fig.add_subplot(211)
	    ax2=fig.add_subplot(212)
	    NfS=[]
	    NDSigmaS=[]
	    SigEffS=[]
	    PhiWCABulkS = []
	    GCs=[]
	    BikS=[]
	    S_CS=[]
	    i=-1
	    Nf_for_EffS=[]
	    VEffS=[]
	    NtS=[]
	#    SCS0s=[]
	    zEffS=[]
	    one_time_iterator=0
	    z_forS=[]
	    rhoplusS,rhominusS=[],[]
	    ccc=-1
	    offset = np.array(range(len(filenameS)))*0.75
	    offset = offset.tolist()
	    offset.reverse()
	    lcorS = []
	    markerS=[]
	    PhiLocalCorS=[]
	    VoltCorS,SigCorS,SigDifS=[],[],[]
	    rhofCorS,rhototCorS = [],[]
	    zetaS=[]
	    SigRefS=[]
	    EvexCorS=[]
	    for (Nm,Np,n0,Sigma_s,area,lam_D,L_bin,z_positions,L_z,Bjerrum,filename,sigWCA,V,muex,Npeff,Nmeff,Phi) in zip(originalN_plusS,originalN_minusS,n0s,SIGMAS,areaS,LDs,L_binS,z_originalS,L_zS,BjerrumS,filenameS,sigWCA_S,Volts,muexEV_S,N_plusS,N_minusS,PhiWCAtot_S):
			i+=1
			dielectric=(4*np.pi*Bjerrum)**-1
		
			Integrand = 0.
			correction=0.
			Nf=np.array([(npl-nm)/(area*L_bin) for (npl,nm) in zip(Np,Nm)])
			Nf=Nf[:len(Nf)/2]
			z_restore = np.array([(x+y)/(L_z*2.) for (x,y) in zip(z_positions[0:len(z_positions)-1],z_positions[1:len(z_positions)])]) #SO, this was the original z_density required as input to the code below
			z_pos = z_restore[:len(z_restore)/2]
			for (y1,y2,z) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos):
				Integrand=Integrand + 0.5*(y1+y2)*L_bin	  
			Sigma_meas = Integrand  #This must be my reportable Sigma, but NDd
			eff=(-Sigma_meas+correction)/(dielectric*temperature/(valency*lam_D))
			ND_Sigma = abs(Sigma_meas)/(dielectric*temperature/(valency*lam_D))
			NDSigmaS.append(ND_Sigma)
			SigRefS.append((dielectric*temperature/(valency*lam_D)))

			Integrand = eff*(dielectric*temperature/(valency*lam_D))
			SigEff=[]
			VEff=[]
		    	volt_correction = Sigma_s/(dielectric*temperature/(valency*lam_D))  - abs(eff)
			Nf_for_Eff=[]
			zEff = []
			rhoplus,rhominus,EVexcess = [],[],[]
			for (y1,y2,z,volt,p,m,x) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos,V[:len(Nf)],Npeff[:len(Nf)],Nmeff[:len(Nf)],muex[:len(Nf)]):  #This is the original line of code -- 07/30/13 16:51:36 
					zEff.append(z)
					SigEff.append(Integrand) #There's something wrong with the first value assigned to this function
					test = volt-volt_correction*z
					VEff.append(test)
					rhoplus.append(p)
					rhominus.append(m)
					EVexcess.append(x)
					Integrand=Integrand + 0.5*(y1+y2)*L_bin	
					if filename in GC_LB_LD_5_20:
		       				Nf_for_Eff.append(0.5*(y1+y2))	        		    		
					else:
		       				Nf_for_Eff.append(y1)

			rhoplus = np.array(rhoplus)/(area*L_bin*n0)
			rhominus = np.array(rhominus)/(area*L_bin*n0)
			zEff = np.array(zEff)*L_z
			zetaS.append(VEff[0])	
			SigEff = np.array(SigEff)/(dielectric*temperature/(valency*lam_D))
		
			ls_x,lw_x='',0.
			ms_x = 4.0
			label_str=''
			alabel=''
			special_ms='x'
			edge_width=0.1
	#		print 'Zetas might need to change when final measurements come out'
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
			elif 'Analyzed_GC_1312' in filename:
				edge_width=0.4
				special_ms='h'
	#			label_str=r'$\tilde \Sigma = 7.0$'
				col = ROYGBIV_map(0.,10)
	#			ze = ['6.72','4.66']	
				if sigWCA==3.0:
					label_str=r'$\tilde \Sigma = 7.1$'
					ze ='6.67'
				elif sigWCA==4.5:
					ze ='10.8'
					label_str=r'$\tilde \Sigma = 6.8$'
				elif sigWCA==6.1:
					ze ='16.23'
					label_str=r'$\tilde \Sigma = 6.4$'
			elif 'Analyzed_GC_1232' in filename:
				edge_width=0.1
				special_ms='o'
				if sigWCA==3.0:
					label_str=r'$5.4$'
					ze ='5.14'
				elif sigWCA==4.5:
					ze ='7.6'
					label_str=r'$5.3$'
				elif sigWCA==6.1:
					ze ='10.93'
					label_str=r'$5.1$'
				special_ms='*'	
				col = ROYGBIV_map(1.,10)		
			elif 'Analyzed_GC_1162' in filename:
				edge_width=0.4
				special_ms='o'	
				col = ROYGBIV_map(2.,10)	
	#			ze = ['4.06','3.28']
	#			ze = ['4.02','3.28']
				if sigWCA == 3.0:
					label_str=r'$4.2$'
					ze = '4.02'				
				elif sigWCA==4.5:
					ze ='5.41'
					label_str=r'$4.1$'
				elif sigWCA==6.1:
					ze ='7.08'
					label_str=r'$3.9$'
			elif 'Analyzed_GC_1136' in filename:
				edge_width=0.1
				special_ms='s'	
	#			label_str=r'$3.6$'
				special_ms='d'		
				col = ROYGBIV_map(3.,10)
	#			ze = ['3.54','2.86']	
				if sigWCA == 3.0:
					label_str=r'$3.6$'
					ze = '3.51'				
				elif sigWCA==4.5:
					ze ='4.51'
					label_str=r'$3.6$'
				elif sigWCA==6.1:
					ze ='5.59'
					label_str=r'$3.5$'
			elif 'Analyzed_GC_1104' in filename:
				edge_width=0.4
				special_ms='s'
				ms_x=3.5
				label_str=r'$3.1$'
				col = ROYGBIV_map(4.,10)
				Force = 0.21359	
	#			ze = ['3.08','2.48']
	#			ze = ['3.11','2.48']
				if sigWCA==3.0:
					ze = '3.11'
					label_str=r'$3.1$'
				elif sigWCA==4.5:
					ze ='3.81'
				elif sigWCA==6.1:
					ze = '4.599'
			elif 'Analyzed_GC_1066' in filename:
				edge_width=0.1
				special_ms='p'	
				ms_x=3.5	
				label_str=r'$2.3$'
				col = ROYGBIV_map(5.,10)
	#			ze = ['2.32','2.01']	
				if sigWCA == 3.0:
					ze = '2.30'
				elif sigWCA==4.5:
					ze ='2.62'
				elif sigWCA==6.1:
					ze = '2.98'
			elif 'Analyzed_GC_1048' in filename:
				edge_width=0.4
				special_ms='p'	
				ms_x=3.5
				label_str=r'$2.0$'
				special_ms='v'	
				col = ROYGBIV_map(6.,10)
	#			ze = ['2.01','1.76']	
				if sigWCA == 3.0:
					ze = '1.97'
				elif sigWCA==4.5:
					ze ='2.13'
				elif sigWCA==6.1:
					ze = '2.31'
			elif 'Analyzed_GC_1034' in filename:
				edge_width=0.1
				special_ms='H'
				ms_x=3.5
				label_str=r'$1.6$'
				col = ROYGBIV_map(7.,10)
	#			ze = ['1.65','1.44']	
				if sigWCA==3.0:
					ze = '1.65'
				elif sigWCA==4.5:
					ze ='1.7'
				elif sigWCA==6.1:
					ze = '1.74'
			elif 'Analyzed_GC_1016' in filename:
				edge_width=0.4
				special_ms='H'	
				label_str=r'$1.0$'
				special_ms='>'	
				col = ROYGBIV_map(8.,10)	
				Force = 0.06948		
				if sigWCA == 3.0:
					ze = '1.07'
				elif sigWCA==4.5:
					ze ='1.05'
				elif sigWCA==6.1:
					ze = '1.09'	
			elif 'Analyzed_GC_1004' in filename:
				edge_width=0.1
				special_ms='D'	
				label_str=r'$0.5$'
				Force = 0.03368
				col = ROYGBIV_map(9.,10)
	#			ze = ['0.50','0.47']
				if sigWCA == 3.0:
					ze = '0.54'
				elif sigWCA==4.5:
					ze ='0.4701'
				elif sigWCA==6.1:
					ze ='0.49'
					label_str=r'$0.49$'
			elif 'Analyzed_GC_1000' in filename:
				edge_width=0.4
				special_ms='<'
				label_str=r'$0.0$'
				col = 'green'
				ze = ['0.00','0.00']
			else:
	#			print '\n\n',filename,Sigma_s,np.max(SIGMAS)
				col = ROYGBIV_map(abs(Sigma_s),np.max(SIGMAS)*1.1,1)
			ms_x = 2.5
			edge_width=0.1


			rhotot = np.array(rhominus) + np.array(rhoplus)
			Phitot = (rhotot*(np.pi/6)*sigWCA**3) / (area*L_bin)
			PhiCS = np.array(Phitot)

	#		print 'file,bulk vol frac,zeta,Sigma = ',filename,np.mean(Phitot[-5:])*0.736419046,VEff[0],eff
	#		print eff
			PhiBulk = np.mean(Phitot[-5:])
			PhiWCABulkS.append(PhiBulk)

			NDrhofree = -(np.array(rhoplus) - np.array(rhominus))/2

			if sigWCA==3.0 or sigWCA==4.5 or sigWCA==6.1:
				CarnahanStarling='yes'
			else:
				CarnahanStarling='no'

			CarnahanStarling='no'

			if sigWCA==3.0:
				special_ms = 'v'
			elif sigWCA == 3.77976:
				special_ms = 'h'
			elif sigWCA ==4.5:
				special_ms = 'd'
			elif sigWCA == 4.7622:
				special_ms = 's'
			elif sigWCA == 5.12993:
				special_ms = '<'
			elif sigWCA == 5.56991:
				special_ms = 'p'
			elif sigWCA == 5.84609:
				special_ms = 'o'			
			elif sigWCA==6.1:
				special_ms = '^'
			markerS.append(special_ms)		
		
			if CarnahanStarling=='yes' and ze!='0.00':	
				  if sigWCA==3.0:
				  	CS_file='CS_0.037_' + ze + '.txt'
			  	  elif sigWCA==4.5:
			  	  	CS_file='CS_0.125_' + ze + '.txt'
			  	  elif sigWCA==6.1:
			  	  	CS_file='CS_0.311_' + ze + '.txt'
				  MMA_CS=[[float(x) for x in line.split()] for line in file(CS_file,"r").readlines()]
				  MMAz = []
				  NDMMA_free,NDMMA_tot = [],[]
				  for line in MMA_CS:
				  	MMAz.append(line[0])
				  	free = (-line[1] + line[2])/2
				  	NDMMA_free.append(free)
					tot = line[1] + line[2]
					NDMMA_tot.append(tot)
				  MMAz = np.array(MMAz)*lam_D/sigWCA-(1+0.299602925)/sigWCA
		       		  ax1.errorbar(MMAz,np.array(NDMMA_free)+offset[i],yerr=None,color=col,ls='-',lw=1.0)	

			zPM = np.array(zEff)/sigWCA-(1+0.299602925)/sigWCA
	#	       	ax1.errorbar(zPM,NDrhofree+offset[i],yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)

			kk=4
			if sigWCA>=6.1 or sigWCA == 5.84609: #0.42 and 0.37
			  	t1 = np.linspace(zPM[0],1,6)[1:].tolist() + np.linspace(1,7,34)[1:].tolist() + np.linspace(7,zPM[-1],25)[1:-1].tolist()
			  	kk=5
	  		elif  sigWCA == 5.56991 or sigWCA == 5.12993 or sigWCA == 4.7622 or sigWCA == 4.5: #0.32 and 0.25 and 0.20 and 0.17
			  	t1 = np.linspace(zPM[0],4,15)[1:].tolist() + np.linspace(4,zPM[-1],25)[1:-1].tolist() 	 
	  		elif  sigWCA == 3.77976: #0.10
			  	t1 = np.linspace(zPM[0],1,4)[1:].tolist() + np.linspace(1,2,4)[1:].tolist() + np.linspace(2,3,4)[1:].tolist() + np.linspace(3,5,4)[1:].tolist() + np.linspace(5,zPM[-1],20)[1:-1].tolist()
	  		elif sigWCA == 3.0: #0.05
			  	t1 = np.linspace(zPM[0],2,5)[1:].tolist() + np.linspace(2,zPM[-1],25)[1:-1].tolist() 
	#		  	t1 = np.linspace(zPM[0],1,3)[1:].tolist() + np.linspace(1,2,3)[1:].tolist() + np.linspace(2,zPM[-1],20)[1:-1].tolist()
			  	kk=5

	    		free_spline = interpolate.LSQUnivariateSpline(zPM,np.array(NDrhofree),t1,k=kk)
			spline_prime = lambda z : free_spline.__call__(z,1)[0]
			single_free_spline = lambda z : free_spline(z)[0]

			PhiLocal_spline = interpolate.LSQUnivariateSpline(zPM,np.array(Phitot),t1,k=kk)
			single_PhiLocal = lambda z : PhiLocal_spline(z)[0]		

			Voltage_spline = interpolate.LSQUnivariateSpline(zPM,np.array(VEff),t1,k=kk)
			single_Voltage = lambda z : Voltage_spline(z)[0]	

			SigEff_spline = interpolate.LSQUnivariateSpline(zPM,SigEff,t1,k=kk)
			single_SigEff = lambda z : SigEff_spline(z)[0]	

			Evex_spline = interpolate.LSQUnivariateSpline(zPM,np.array(EVexcess),t1,k=kk)
			single_Evex = lambda z : Evex_spline(z)[0]	


			zPM_plot = np.linspace(zPM[0],zPM[-1],len(zPM)*2) # 2X as much spacing than the bins
		       	ax1.errorbar(zPM_plot,free_spline(zPM_plot)+offset[i],yerr=None,color='k',ls ='-',lw = 0.5)
	   
			x0s = np.linspace(0.0,8,30)
	#		print '\n',filename,special_ms
			primeroots=[]
			for x0 in x0s:
				SolverResults = fsolve(spline_prime, x0,full_output = 1)
				if 'solution' in SolverResults[3]:
					root1 = np.round(SolverResults[0],4)	
					if (root1 not in primeroots) and root1<=9 and root1>0:
						primeroots.append(root1)
			primeroots.sort()
			primeroots = primeroots[::-1]

			if '1000' in filename:
				primeroots=[]
				print 'Skipping for no field...'


			if len(primeroots)!=0:
	#			print 'Possible lcors = ',primeroots
				ax1.errorbar(primeroots,free_spline(primeroots)+offset[i],color='k',marker=special_ms,mew=0.0,ms=ms_x*1.05,ls=ls_x,lw=0)#,label=label_str)#,mec = col)

			if i==0:
				lcor_set = []
				sigOld = sigWCA
			elif sigOld!=sigWCA:
	#			print lcorS[-1],sigWCA
				lcor_set = []
				sigOld = sigWCA			
		
			pretest = len(lcor_set)
			jj=-1
			for (far,near) in zip(primeroots[0:len(primeroots)-1],primeroots[1:len(primeroots)]):
				if far<=0:
					far = 0.
					lcorS.append(far)
					j+=1

				spacing = far - near
		
				if jj==-1 and (spacing<=0.505 and spacing>0.15):

					rhonear = single_free_spline(near)
					rhofar = single_free_spline(far)
					rhodif = np.abs(rhofar - rhonear)
					bulkdifS = []
					for z in zPM[-90:-15]:
						bulkdif = single_free_spline(z) - single_free_spline(z+spacing)
						bulkdifS.append(bulkdif)
					bulkRMS = np.sqrt(np.mean(np.array(bulkdifS)**2))

					if rhodif > 2.*bulkRMS: #This tests for a significant fluctuation
						z_bulkwards,rho_bulkwards = [],[]
						z_wallwards,rho_wallwards = [],[]
						for (x,y) in zip(zPM,NDrhofree):
							if x>=near-1.0 and x<=far and x>0:
								z_wallwards.append(x)
								rho_wallwards.append(y)								
							elif x>far and x<=far+1.0:								
								z_bulkwards.append(x)
								rho_bulkwards.append(y)
							elif x>far+1.0:
								z_wallwards = np.array(z_wallwards)
								rho_wallwards = np.array(rho_wallwards)
								z_bulkwards = np.array(z_bulkwards)
								rho_bulkwards = np.array(rho_bulkwards)	
								break
						z_cor = np.linspace(z_wallwards[0],far,20)
						ztest = np.array(z_cor.tolist() + np.linspace(far+z_cor[1]-z_cor[0],z_bulkwards[-1],40).tolist())

						expslope, expintercept, expr_value, expp_value, expstd_err = stats.linregress(z_bulkwards,np.log(rho_bulkwards))
						exprhotest = np.exp(ztest*expslope + expintercept)
						dif_from_exp = np.array(free_spline(ztest)) - exprhotest
						MaxDeviationFromExpFit = np.max(np.abs(dif_from_exp[20:]))
						fitRMS_exp = np.sqrt(((np.exp(z_bulkwards*expslope+expintercept) - rho_bulkwards)** 2).mean())
						corRMS_exp = np.sqrt(((np.exp(z_wallwards*expslope+expintercept) - rho_wallwards)** 2).mean())
	#					print 'EXPONENTIAL (expR^2,expslope,expintercept) = (',expr_value**2,',',expslope,',',expintercept,')'
	#					print 'EXPONENTIAL (max, min, fit extrema) = (',max(dif_from_exp[:20]),',',min(dif_from_exp[:20]),',',MaxDeviationFromExpFit,')'
	#					print 'EXPONENTIAL (fitRMS, corRMS, cor/fit ) = ',fitRMS_exp,',',corRMS_exp,',',corRMS_exp/fitRMS_exp,')'

							###Old way of doing, worth preserving
	#					dif_spline = interpolate.LSQUnivariateSpline(z_cor,np.array(free_spline(z_cor)) - np.exp(z_cor*expslope + expintercept),np.linspace(z_cor[0],z_cor[-1],5)[1:-1],k=5)
	#					dif_spline_prime = lambda z : dif_spline.__call__(z,1)[0]
	#					x0s = np.linspace(z_cor[0],z_cor[-1],20)

						dif_spline = interpolate.LSQUnivariateSpline(ztest,np.array(free_spline(ztest)) - np.exp(ztest*expslope + expintercept),np.linspace(ztest[0],ztest[-1],5)[1:-1],k=5)
						dif_spline_prime = lambda z : dif_spline.__call__(z,1)[0]

						difroots = []
						x0s = np.linspace(z_cor[0],z_cor[-1],20)
						for x0 in x0s:
							SolverResults =  fsolve(dif_spline_prime, x0,full_output = 1)
							if 'solution' in SolverResults[3]:
								root1 = np.round(SolverResults[0],3)	
								if (root1 not in difroots) and root1<=9:
									difroots.append(root1)
						difroots.sort()
						difroots = np.array(difroots[::-1])

	#					print 'Dif roots and spacing = ',difroots,difroots[:-1]-difroots[1:]
	#					print 'DifValues at roots = ',dif_spline(difroots)
	#					print 'Delta DifValues = ',dif_spline(difroots)[:-1]-dif_spline(difroots)[1:]
						LDAdifroots = []
						x0s = np.linspace(ztest[20],ztest[-1],20)
						for x0 in x0s:
							SolverResults =  fsolve(dif_spline_prime, x0,full_output = 1)
							if 'solution' in SolverResults[3]:
								root1 = np.round(SolverResults[0],3)	
								if (root1 not in LDAdifroots) and root1<=9:
									LDAdifroots.append(root1)
						LDAdifroots.sort()
						LDAdifroots = np.array(LDAdifroots[::-1])

	#					print 'LDA Dif roots and spacing = ',LDAdifroots,LDAdifroots[:-1]-LDAdifroots[1:]
	#					print 'LDA DifValues at roots = ',dif_spline(LDAdifroots)
	#					print 'LDA Delta DifValues = ',dif_spline(LDAdifroots)[:-1]-dif_spline(LDAdifroots)[1:]					
						if max(dif_from_exp[:20]) > MaxDeviationFromExpFit and min(dif_from_exp[:20]) < -MaxDeviationFromExpFit:
							if corRMS_exp > 5.*fitRMS_exp:
								if len(lcor_set)==0:
									jj+=1								
									ax1.errorbar(ztest,exprhotest+offset[i],yerr=None,color='grey',ls ='-')
									ax2.errorbar(ztest,dif_from_exp+offset[i],yerr=None,color=col,ls ='-')
								       	ax1.errorbar(far,single_free_spline(far)+offset[i],yerr=None,color='white',marker=special_ms,mew=edge_width,ms=ms_x*1.75,mec = col)							       	
									dif_single = single_free_spline(far) - np.exp(far*expslope + expintercept)
								       	ax2.errorbar(far,dif_single+offset[i],yerr=None,color='white',marker=special_ms,mew=edge_width,ms=ms_x*1.75,mec = col)
									ax2.errorbar(difroots,dif_spline(difroots)+offset[i],yerr=None,color='k',marker=special_ms,ms=ms_x*1.05,mew=0.0,lw=0)
								       	lcorS.append(far)
								       	lcor_set.append(far)	
	#								ax2.errorbar(LDAdifroots,dif_spline(LDAdifroots)+offset[i],yerr=None,color='purple',marker=special_ms,ms=ms_x*1.05,mew=0.0,lw=0)
	#								print 'Lcor was found to be ',far									       								       		
							       	elif far <= lcor_set[-1]+0.5*0.954028946**3:
									jj+=1
									ax1.errorbar(ztest,exprhotest+offset[i],yerr=None,color='grey',ls ='-')
									ax2.errorbar(ztest,dif_from_exp+offset[i],yerr=None,color=col,ls ='-')
								       	ax1.errorbar(far,single_free_spline(far)+offset[i],yerr=None,color='white',marker=special_ms,mew=edge_width,ms=ms_x*1.75,mec = col)
									dif_single = single_free_spline(far) - np.exp(far*expslope + expintercept)
								       	ax2.errorbar(far,dif_single+offset[i],yerr=None,color='white',marker=special_ms,mew=edge_width,ms=ms_x*1.75,mec = col)
									ax2.errorbar(difroots,dif_spline(difroots)+offset[i],yerr=None,color='k',marker=special_ms,ms=ms_x*1.05,mew=0.0,lw=0)
								       	lcorS.append(far)
									lcor_set.append(far)	
	#								ax2.errorbar(LDAdifroots,dif_spline(LDAdifroots)+offset[i],yerr=None,color='purple',marker=special_ms,ms=ms_x*1.05,mew=0.0,lw=0)
	#								print 'Lcor was found to be ',far		
			if len(lcor_set) == pretest:
				lcorS.append(0.)
			       	lcor_set.append(0.)
				print "WARNING NO LCOR FOUND!"
	#		lcorS[-1] = lcorS[-1] + 0.5*sigWCA
	#		if filename== 'Analyzed_GC_1312_0.024179_0.1_3.77976_500.0_23.73.txt':
	#			lcorS[-1] = 2.205 
	#			print 'Artificial forcing!'

			VoltageCor = single_Voltage(lcorS[-1])
			VoltCorS.append(VoltageCor)
			SigCor = ND_Sigma  + single_SigEff(lcorS[-1]) # = -Sigma - SigEff(lcor) [The vars are negatives of e/o]
			corindex = len(zPM[zPM<=(lcorS[-1]-0.5*L_bin/sigWCA)])
			binweight = sigWCA*abs(lcorS[-1] - zPM[corindex])/L_bin
			NpCor = np.sum(rhoplus[:corindex])+binweight*rhoplus[corindex] #= Np/(Axy * n0 * Lbin)
			NmCor = np.sum(rhominus[:corindex])+binweight*rhominus[corindex]
		
			rhofCor = -L_bin*0.5*(NpCor - NmCor)/((lcorS[-1]-zPM[0])*sigWCA)
			NDrhofcor = (rhofCor*n0)/(valency/(4*np.pi*Bjerrum*sigWCA**2))
			rhofCorS.append(NDrhofcor)	
		
			SigCorS.append(SigCor*(sigWCA/lam_D))
			SigDif = ND_Sigma - SigCor   #This is probably nonsense right now
			SigDifS.append(SigDif)
	
			rhototCor = L_bin*(NpCor + NmCor)/((lcorS[-1]-zPM[0])*sigWCA) #This is rhoTOT
			rhototCorS.append(rhototCor)
			PhiCor = rhototCor*n0*(np.pi/6)*sigWCA**3
			PhiLocalCorS.append(PhiCor)
			EvexCorS.append(single_Evex(lcorS[-1]))
	#		print 'elcor, (sigcor/rhocor)  = ',lcorS[-1],(SigCor*(sigWCA/lam_D) / NDrhofcor)
	    i=-1
	    for f in filenameS:
		i+=1
	    	ax2.errorbar(zPM_plot,np.zeros(len(zPM_plot))+offset[i],yerr=None,color='grey',ls ='-',lw = 0.5)
	    ax1.set_ylabel(r'$ -\rho / 2 e n^{\rm B} + {\rm offset} $',fontsize=10.) # = (n_+ - n_-) / 2 \rho^{\rm B}
	    ax2.set_ylabel(r'$ \tilde \rho - \tilde \rho_{\rm fit}(z>\ell_{\rm cor})  + {\rm offset}$',fontsize=10.) 
	    ax2.set_xlabel(r'$(z-\delta_{\rm w}) / \sigma$',fontsize=10.)
	    
	    plt.setp(ax1.get_xticklabels(), fontsize=8.)  
	    plt.setp(ax1.get_yticklabels(), fontsize=8.)
	    plt.setp(ax2.get_xticklabels(), fontsize=8.)
	    plt.setp(ax2.get_yticklabels(), fontsize=8.)

	    ax1.set_xlim(-0.5,10)
	    ax1.set_ylim(-0.25,14.) 
	    ax2.set_xlim(-0.5,10)
	    ax2.set_ylim(-0.5,np.max(offset)*1.05) 

	    if len(filenameS)==10:
	#	ax1.legend(loc='upper right',numpoints=1,prop={"size":8},columnspacing=0.07,borderpad=0.18,labelspacing=0.10,handletextpad=0.07,ncol=2,markerscale=2.0)     
	    	plt.savefig('FindingOscillations_'+str(np.round(np.mean(PhiWCABulkS),3))+'.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
	    plt.close()

	    NDSigmaS=np.array(NDSigmaS) #x
	    lcorS = np.array(lcorS)

#	    print '\n\nFitting SigCor versus DelVoltcor:'
#	    for (elWCA,elLD,Sigtot,rhocorWCA,sigcor,voltcor) in zip(lcorS,np.array(lcorS)*np.array(sigWCA_S)/np.array(LDs),NDSigmaS,rhofCorS,SigCorS,np.array(zetaS) - np.array(VoltCorS)):
#	    	print elWCA,'\t',Sigtot,'\t',rhocorWCA,'\t',voltcor,'\t',sigcor

#	    print '\n\n'
#	    print 'MAX ND Sigma =',np.max(NDSigmaS)
#	    print 'plotSigma =',NDSigmaS.tolist()
##	    print fail
#	    
#	    print '\nplotlcorS = ',lcorS.tolist()
#	    print '\nVoltCor = ',VoltCorS
#	    print '\nPhiBulks, min, max = ',PhiWCABulkS,np.mean(PhiWCABulkS),np.min(PhiWCABulkS),np.max(PhiWCABulkS)
#	    print '\nLDs, mean LD =',LDs,np.mean(LDs),'\t\t',np.min(LDs),np.max(LDs),'\t\t',np.max([np.abs(np.mean(LDs)-np.min(LDs))/np.mean(LDs),np.abs(np.mean(LDs)-np.max(LDs))])/np.mean(LDs)

#	    
#	    print '\nn0s =',n0s,np.mean(n0s),np.min(n0s),np.max(n0s)
#	    print '\n\n'
	    
	    if len(filenameS)>10 and len(filenameS)!=343:
		    groupColorS = ['red','orange','purple','green','blue','cyan','violet','grey']
		    groupLabelS = [r'$\bar{\Phi}_{\rm B}=0.42$',r'$0.37$',r'$0.32$',r'$0.25$',r'$0.20$',r'$0.17$',r'$0.10$',r'$0.05$']
		    TheoryNames = ['_0.311_17.txt','_0.274_17.txt','_0.237_17.txt','_0.185_17.txt','_0.148_17.txt','_0.125_17.txt','_0.074_17.txt','_0.037_17.txt']
	   	    markerS = ['^','o','p','<','s','d','h','v']
	    
		    print 'Plotting PhiBulk,SigApp veruss l_cor^CS contour'
		    import matplotlib.colors as colors
		    import matplotlib.cm as cmx
		    fig=plt.figure()
		    ax1=fig.add_subplot(111)
		    MMA_Phi = [0.03,0.05,0.07,0.10,0.13,0.15,0.17,0.20,0.23,0.25,0.27,0.30,0.33,0.35,0.37,0.40,0.43,0.45,0.47,0.50]
		    MMA_fitcor = [2.123775331,1.695857713,1.43592424,1.173463673,0.993542747,0.90197506,0.826621907,0.735990472,0.665033731,0.62568418,0.591325724,0.547337166,0.510309068,0.488864226,0.4696895,0.443728052,0.420911332,0.407282579,0.394633661,0.377469508]            
		    MMA_PMcor = [2.809350692,2.337064891,2.045791971,1.735251232,1.500128432,1.369855345,1.25664516,1.114103312,0.999057484,0.934707058,0.878514728,0.806857238,0.747127877,0.712766209,0.682137173,0.641266922,0.605758148,0.58467426,0.565249428,0.539006824]            

		    ##TO INCLUDE ZEROS
		    MMA_Phi = [0.0]+MMA_Phi
		    MMA_fitcor =[100.0] + MMA_fitcor
		    MMA_PMcor =[100.0] + MMA_PMcor

	##		##Just the DATA
	#            MMA_Phi = [0.05,0.07,0.10,0.13,0.15,0.17,0.20,0.23,0.25,0.27,0.30,0.33,0.35,0.37,0.40,0.43]
	#            MMA_fitcor = [1.695857713,1.43592424,1.173463673,0.993542747,0.90197506,0.826621907,0.735990472,0.665033731,0.62568418,0.591325724,0.547337166,0.510309068,0.488864226,0.4696895,0.443728052,0.420911332]            
	#            MMA_PMcor = [2.337064891,2.045791971,1.735251232,1.500128432,1.369855345,1.25664516,1.114103312,0.999057484,0.934707058,0.878514728,0.806857238,0.747127877,0.712766209,0.682137173,0.641266922,0.605758148]            
	#	
	   
		    SigRange = np.linspace(0,7.5,100)
		    fixLD,fixLB = 14.6,0.1
		    SigCorCS,rhocor = MMA_PMcor,0.0956
	#	    SigCorCS,rhocor = MMA_fitcor,0.069471312
		    
		    CS_Sigmas,CS_PhiS,CS_lcorS = [],[],[]
		    for (PhiB,SigDif) in zip(MMA_Phi,SigCorCS):  
		    	for SigTot in SigRange:
		    		if (SigTot-SigDif)<=0:
		    			CS_lcorS.append(0)
	    			else:
	    				fixsig = (24*PhiB*fixLB*fixLD**2)**(1/3.)
	    				NDrhocor = 2*rhocor*(fixLD/fixsig)**2
	    				lcorCS = ((SigTot-SigDif) / NDrhocor)
	    				CS_lcorS.append(lcorCS)
		    		CS_Sigmas.append(SigTot)
		    		CS_PhiS.append(PhiB)
		    CS_Sigmas=np.array(CS_Sigmas) #x
		    CS_PhiS=np.array(CS_PhiS) #y
		    CS_lcorS=np.array(CS_lcorS) #z
		    xi = np.linspace(0,np.max(CS_Sigmas),25)
		    yi = np.linspace(0,np.max(CS_PhiS),25)
		    zi = griddata(CS_Sigmas,CS_PhiS,CS_lcorS,xi,yi,interp='nn') #For other gridding options see: http://matplotlib.org/api/mlab_api.html?highlight=griddata#matplotlib.mlab.griddata
		    shelves = 6#*2
		    CS = plt.contour(xi,yi,zi,shelves,linewidths=0.5,colors='k')
		    colscheme = cm = plt.get_cmap('gist_rainbow')
		    CS = plt.contourf(xi,yi,zi,shelves,cmap=colscheme,vmax=abs(zi).max(), vmin=-abs(zi).max())    
		    cb = plt.colorbar(format='%i')
		    for t in cb.ax.get_yticklabels():
		    	t.set_fontsize(8.0)
		    cb.set_label(r'$\ell_{\rm cor}^{\rm CS} / \sigma$',fontsize=10.)  

		    if len(filenameS) == 35:
		    	plt.plot(NDSigmaS,PhiWCABulkS,c='white',ls=':')
	    	    else:
			    plt.plot([0,np.max(NDSigmaS)],[np.min(PhiWCABulkS),np.min(PhiWCABulkS)],c='k',ls=':')
			    plt.plot([np.max(NDSigmaS),np.max(NDSigmaS)],[np.min(PhiWCABulkS),np.max(PhiWCABulkS)],c='k',ls=':')
			    plt.plot([np.max(NDSigmaS),np.max(NDSigmaS)],[np.min(PhiWCABulkS),np.max(PhiWCABulkS)],c='k',ls=':')
			    plt.plot([0,np.max(NDSigmaS)],[np.max(PhiWCABulkS),np.max(PhiWCABulkS)],c='k',ls=':')
			    plt.plot([0,0],[np.max(PhiWCABulkS),np.min(PhiWCABulkS)],c='k',ls=':')

		    if len(MMA_Phi)==16:
	    	        plt.xlim(0,7.5)
		    	plt.ylim(0,0.45) 
	    	    else:
	    	        plt.xlim(0,7.5)
		    	plt.ylim(0,0.50)  
		    ax1.set_xlabel(r'$\Sigma / \Sigma_{\rm ref}$' ,fontsize=10.)
		    ax1.set_ylabel(r'$\Phi_{\rm B}$',fontsize=10.)
		    plt.setp(ax1.get_xticklabels(), fontsize=8.)
		    plt.setp(ax1.get_yticklabels(), fontsize=8.)
		    fig.set_size_inches(3.37,2.5)
		    plt.savefig('F?_lcorCSPM_contour.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
		##Biksucks
		    plt.close()

		    print 'Plotting PhiBulk,SigApp veruss l_cor contour'
		    import matplotlib.colors as colors
		    import matplotlib.cm as cmx
		    fig=plt.figure()
		    ax1=fig.add_subplot(111)
		    NDSigmaS=np.array(NDSigmaS) #x
		    PhiWCABulkS=np.array(PhiWCABulkS) #y
		    lcorS=np.array(lcorS) #z
		    xi = np.linspace(0,np.max(NDSigmaS)*1.0,25)
		    yi = np.linspace(0,np.max(PhiWCABulkS)*1.0,25)
		    zi = griddata(NDSigmaS,PhiWCABulkS,lcorS,xi,yi,interp='nn') #For other gridding options see: http://matplotlib.org/api/mlab_api.html?highlight=griddata#matplotlib.mlab.griddata
		    shelves = 6
		    CS = plt.contour(xi,yi,zi,shelves,linewidths=0.5,colors='k')
		    colscheme = cm = plt.get_cmap('gist_rainbow')
		    CS = plt.contourf(xi,yi,zi,shelves,cmap=colscheme,vmax=abs(zi).max(), vmin=-abs(zi).max())    
		    cb = plt.colorbar()
		    for t in cb.ax.get_yticklabels():
		    	t.set_fontsize(8.0)
		    cb.set_label(r'$\ell_{\rm cor} / \sigma$',fontsize=10.)  
		    for (i,col) in enumerate(groupColorS):
		    	plotSigS = NDSigmaS[i*10:10*(i+1)]
			plotPhiBulkS = PhiWCABulkS[i*10:10*(i+1)]
			if len(filenameS)==88:
			    	plotSigS = NDSigmaS[i*11:11*(i+1)]
				plotPhiBulkS = PhiWCABulkS[i*11:11*(i+1)]				
			mark = markerS[i]
			plt.scatter(plotSigS,plotPhiBulkS,marker=mark,linewidth='0',c='k',s=7,zorder=10) #Not sure what zorder is...
	#	    plt.xlim(0,np.max(NDSigmaS)*1.05)
	#	    plt.ylim(0,np.max(PhiWCABulkS)*1.05)    	
		    plt.xlim(0,7.5)
		    plt.ylim(0,0.45)  
		    ax1.set_xlabel(r'$\Sigma / \Sigma_{\rm ref}$' ,fontsize=10.)
		    ax1.set_ylabel(r'$\Phi_{\rm B}$',fontsize=10.)
		    plt.setp(ax1.get_xticklabels(), fontsize=8.)
		    plt.setp(ax1.get_yticklabels(), fontsize=8.)
		    fig.set_size_inches(3.37,2.5)
		    plt.savefig('F?_lCor_contour.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
		##Biksucks
		    plt.close()

	#	    print fail

		    print 'Plotting TOTAL Sigma versus Volt'
		    fig=plt.figure()
		    ax1=fig.add_subplot(111)
		    volt_theory = np.linspace(0,16)

		    ax1.errorbar(volt_theory,volt_theory,yerr=None,ls='--',lw=1.0,color='k',marker='None',label =r'${\rm DH}$')

		    GC_cap = 2*np.sinh(volt_theory/2.)
		    ax1.errorbar(volt_theory,GC_cap,yerr=None,ls=':',lw=1.0,color='k',marker='None',label =r'${\rm GC}$')
		    
		    QV_Bik=np.array([[float(x) for x in line.split()] for line in file('Bik_0.05_16.txt',"r").readlines()])
		    ax1.errorbar(QV_Bik[:,4],QV_Bik[:,5],yerr=None,ls='-.',lw=1.0,color='grey',marker='None')#,label =r'${\rm Bik}$')

		    QV_Bik=np.array([[float(x) for x in line.split()] for line in file('Bik_0.10_16.txt',"r").readlines()])
		    ax1.errorbar(QV_Bik[:,4],QV_Bik[:,5],yerr=None,ls='-.',lw=1.0,color='violet',marker='None')#,label =r'${\rm Bik}$')

	    	    ax1.plot([-1,-1],[-1,-1],color = 'k', ls = '-.',lw = 1.0,marker = 'None',label =r'${\rm Bik}$')               
	    	    ax1.plot([-1,-1],[-1,-1],color = 'k', ls = '-',lw = 1.0,marker = 'None',label =r'${\rm CS}$')
#	    	    ax1.plot([-1,-1],[-1,-1],color = 'k', ls = 'None',lw = 1.0,marker = 'None',label =r'$ $')
		    for (i,col) in enumerate(groupColorS):
		    	plotzetaS = zetaS[i*10:10*(i+1)]
		    	plotNDSigmaS = NDSigmaS[i*10:10*(i+1)]
		    	if len(filenameS)==88:
			    	plotzetaS = zetaS[i*11:11*(i+1)]
			    	plotNDSigmaS = NDSigmaS[i*11:11*(i+1)]		    		

			QV_CS=np.array([[float(x) for x in line.split()] for line in file('CS'+TheoryNames[i],"r").readlines()])
	   		ax1.errorbar(QV_CS[:,4],QV_CS[:,5],yerr=None,ls='-',lw=1.0,color=col,marker='None')#,label =groupLabelS[i])

				#The following adds Bikerman theory
	#		QV_CS=np.array([[float(x) for x in line.split()] for line in file('Bik'+TheoryNames[i],"r").readlines()])
	#   		ax1.errorbar(QV_CS[:,4],QV_CS[:,5],yerr=None,ls='-.',lw=1.0,color=col,marker='None')#,label =groupLabelS[i])

			mark = markerS[i]
		    	ax1.plot(plotzetaS,plotNDSigmaS,color=col,ls='',marker = mark)#,label =groupLabelS[i])
		    	ax1.plot([-1,-1],[-1,-1],color = col, ls = 'None',lw = 0.9,marker = mark,label =groupLabelS[i])

		    ax1.errorbar(volt_theory,GC_cap,yerr=None,ls=':',lw=1.0,color='k',marker='None')#,label =r'${\rm GC}$')
	#	    ax1.set_ylabel(r'$-\tilde \Sigma = ( \Sigma_{\rm cor} +  \Sigma_{\rm dif} ) / \Sigma_{\rm ref}$',fontsize=10.)
		    ax1.set_ylabel(r'$\Sigma / \Sigma_{\rm ref}$',fontsize=10.)
		    ax1.set_xlabel(r'$(\phi_0 - \phi_{\rm B}) / \phi_{\rm T}$',fontsize=10.)
		    plt.setp(ax1.get_xticklabels(), fontsize=8.)
		    plt.setp(ax1.get_yticklabels(), fontsize=8.)
		    ax1.legend(loc='best',numpoints=1,prop=dict(size=8.),columnspacing=0.11,borderpad=0.15,labelspacing=0.1,ncol=3,handletextpad=0.07,markerscale=0.75)     
		    ax1.set_xlim(0,16.5) 
#		    ax1.set_ylim(0,12.0)
		    ax1.set_ylim(0,8.0)
#		    fig.set_size_inches(3.37,3.0)
		    fig.set_size_inches(3.37,2.5)
		    plt.savefig('F?_ChargeVoltage_tot.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
	#	##Biksucks
		    plt.close()





		    print 'NEW Plotting V_cor veruss lcor'
		    fig=plt.figure()
		    ax1=fig.add_subplot(111)
		    for (i,col) in enumerate(groupColorS):
		    	plotlcorS = lcorS[i*10:10*(i+1)]
		    	plotDelVoltageCorS = np.array(zetaS[i*10:10*(i+1)]) - np.array(VoltCorS[i*10:10*(i+1)])
		    	if len(filenameS)==88:
		    		plotlcorS = lcorS[i*11:11*(i+1)]
			    	plotDelVoltageCorS = np.array(zetaS[i*11:11*(i+1)]) - np.array(VoltCorS[i*11:11*(i+1)])	
		    	
		    	ax1.plot(plotlcorS,plotDelVoltageCorS,color=col,ls='',lw=0.9,marker = markerS[i],label =groupLabelS[i])
		    ax1.set_xlabel(r'$\ell_{\rm cor} / \sigma$',fontsize=10.)
		    ax1.set_ylabel(r'$\phi_{\rm cor}/\phi_{\rm T}$',fontsize=10.)
		    plt.setp(ax1.get_xticklabels(), fontsize=8.)
		    plt.setp(ax1.get_yticklabels(), fontsize=8.)
		    ax1.legend(loc='best',numpoints=1,prop=dict(size=8.),columnspacing=0.10,borderpad=0.25,labelspacing=0.1,ncol=2,handletextpad=0.05,markerscale=0.75)     
		    ax1.set_xlim(0,7.5)       
		    ax1.set_ylim(0.,14.)
		    fig.set_size_inches(3.37,3.0)
#		    plt.savefig('F?_Vcor_vs_lCor.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
	#	##Biksucks
		    plt.close()

#		    print fail



		    print 'Plotting Sigma_cor versus Volt_cor'
		    fig=plt.figure()
		    ax1=fig.add_subplot(111)
		    for (i,col) in enumerate(groupColorS):
		    	plotSigCorS = SigCorS[i*10:10*(i+1)]
		    	plotDelVoltageCorS = np.array(zetaS[i*10:10*(i+1)]) - np.array(VoltCorS[i*10:10*(i+1)])
		    	if len(filenameS)==88:
			    	plotSigCorS = SigCorS[i*11:11*(i+1)]
			    	plotDelVoltageCorS = np.array(zetaS[i*11:11*(i+1)]) - np.array(VoltCorS[i*11:11*(i+1)])		    		
		    	mark = markerS[i]
		    	ax1.plot(plotDelVoltageCorS,plotSigCorS,color=col,ls='',marker=mark,label =groupLabelS[i])

		    fitVolt = np.linspace(0,15,100)

		    fitSig = np.sqrt(2*0.0956*fitVolt)
		    ax1.plot(fitVolt,fitSig,color='k',ls='-',lw=0.9,label =r'$\sqrt{2  \tilde \phi_{\rm cor}  \bar{\tilde \rho}_{\rm cor} }$')	    

		    fitSig = np.sqrt(2*0.069471312*fitVolt)
		    ax1.plot(fitVolt,fitSig,color='red',ls='-',lw=0.9,label =r'$\sqrt{2 \tilde \phi_{\rm cor}  \tilde \rho_{\rm cor}^{\rm fit}} $')	    
		    		    
		    ax1.set_ylabel(r'$\Sigma_{\rm cor} (\sigma/ \Sigma_{\rm ref} \lambda_{\rm D})$',fontsize=10.) # = \int_0^{\ell_{\rm cor} } \tilde \rho {\rm d} \tilde z 
		    ax1.set_xlabel(r'$(\phi_0 - \phi_{ \rm cor }) / \phi_{\rm T}$',fontsize=10.)
		    plt.setp(ax1.get_xticklabels(), fontsize=8.)
		    plt.setp(ax1.get_yticklabels(), fontsize=8.)
		    ax1.legend(loc='best',numpoints=1,prop=dict(size=8.),columnspacing=0.10,borderpad=0.25,labelspacing=0.1,ncol=2,handletextpad=0.05,markerscale=0.75)     
		    ax1.set_xlim(0,14)       
		    ax1.set_ylim(0,2)
	#	    fig.set_size_inches(3.37,3.00)
		    fig.set_size_inches(3.37,3.25)
#		    plt.savefig('F?_ChargeVoltage_cor.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
	#	##Biksucks
		    plt.close()








		    numpoints = 0
		    print 'NEW ONLY lcor>0: LOG LOG Plotting Sigma_cor versus Volt_cor'
		    fig=plt.figure()
#    		    fig.set_size_inches(3.37,2.50)
    		    fig.set_size_inches(1.15,0.95)
		    ax1=fig.add_subplot(111)
		    for (i,col) in enumerate(groupColorS):
		    	plotSigCorS = SigCorS[i*10:10*(i+1)]
		    	plotDelVoltageCorS = np.array(zetaS[i*10:10*(i+1)]) - np.array(VoltCorS[i*10:10*(i+1)])
		    	includelcorS = lcorS[i*10:10*(i+1)]
		    	if len(filenameS)==88:
			    	plotSigCorS = SigCorS[i*11:11*(i+1)]
			    	plotDelVoltageCorS = np.array(zetaS[i*11:11*(i+1)]) - np.array(VoltCorS[i*11:11*(i+1)])		
    			    	includelcorS = lcorS[i*11:11*(i+1)]    
    			tempX,tempY = [],[]		
			for (xx,yy,lcor) in zip(plotDelVoltageCorS,plotSigCorS,includelcorS):
				if lcor>0:
#					print xx,'\t',yy,'\t',lcor
					tempX.append(xx)
					tempY.append(yy)
				else:
					print lcor,yy
			plotDelVoltageCorS,plotSigCorS = tempX,tempY
			numpoints+=len(tempX)

		    	mark = markerS[i]		    	
		    	ax1.plot(plotDelVoltageCorS,plotSigCorS,color=col,ls='',marker=mark,label =groupLabelS[i])

		    print 'numpoints = ',numpoints
		    fitVolt = np.linspace(0,20,10000)
		    

		    fitSig = np.sqrt(2*0.0956*fitVolt)
		    ax1.plot(fitVolt,fitSig,color='k',ls='-',lw=0.9,label =r'$\sqrt{2  \tilde \phi_{\rm cor}  \bar{\tilde \rho}_{\rm cor} }$')	    

		    bbb=3./4
		    ax1.plot(fitVolt,(1.-bbb)*fitVolt**bbb,color='red',ls='-',lw=0.9,label =r'$  \frac{1}{4} \tilde \phi_{\rm cor}^{\frac{3}{4}}$')	

#		    fitSig = np.sqrt(2*0.069471312*fitVolt)
#		    ax1.plot(fitVolt,fitSig,color='red',ls='-',lw=0.9,label =r'$\sqrt{2 \tilde \phi_{\rm cor}  \tilde \rho_{\rm cor}^{\rm fit}} $')	    
		    		    
		    ax1.set_ylabel(r'$\log \left [ - \Sigma_{\rm cor} (\sigma/ \Sigma_{\rm ref} \lambda_{\rm D}) \right ]$',fontsize=10.) # = \int_0^{\ell_{\rm cor} } \tilde \rho {\rm d} \tilde z 
		    ax1.set_xlabel(r'$\log \left [(\phi_0 - \phi_{ \rm cor }) / \phi_{\rm T} \right ]$',fontsize=10.)
		    plt.setp(ax1.get_xticklabels(), fontsize=8.)
		    plt.setp(ax1.get_yticklabels(), fontsize=8.)
		    ax1.legend(loc='upper left',numpoints=1,prop=dict(size=8.),columnspacing=0.10,borderpad=0.25,labelspacing=0.1,ncol=2,handletextpad=0.10,markerscale=0.75,handlelength = 1.5)     

		    ax1.set_yscale('log')
		    ax1.set_xscale('log')
    		    fig.set_size_inches(3.37,3.00)
#		    plt.savefig('F?_nolimsLOGLOG_ChargeVoltage_cor_BIG.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
		    
#    		    ax1.set_xlim(0,14)       
#		    ax1.set_ylim(0,2)

    		    ax1.set_xlim(10**(-1.),15)       
		    ax1.set_ylim(0.05,2.5)
		    plt.savefig('F?_zoomLOGLOG_ChargeVoltage_cor_BIG.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=Falsetransparent=False,bbox_inches='tight')   
		    
#    		    fig.set_size_inches(1.15,0.95)
#    		    ax1.set_xlim(10**(-2.),20)       
#		    ax1.set_ylim(0,2)
#		    plt.savefig('F?_zoomlimsLOGLOG_ChargeVoltage_cor_small.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
#		    		    
#    		    fig.set_size_inches(1.15,0.95)
#    		    ##Save
#		    plt.savefig('F?_zoomlimsLOGLOG_ChargeVoltage_cor_small.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=True,bbox_inches='tight')   
#		    ax1.set_xlim(0,14)       
#		    ax1.set_ylim(0,2)
#		    plt.savefig('F?_LOGLOGb_ChargeVoltage_cor.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
	#	##Biksucks
		    plt.close()

#		    print fail

		    print 'Plotting rho(lcor) versus ~SigmaCor'
		    fig=plt.figure()
		    ax1=fig.add_subplot(111)
		    for (i,col) in enumerate(groupColorS):
		    	plotSigCorS = SigCorS[i*10:10*(i+1)]
			plotrhofCorS = rhofCorS[i*10:10*(i+1)]
			mark = markerS[i]
		    	ax1.plot(plotSigCorS,plotrhofCorS,color=col,ls='',lw=0.9,marker = mark,ms=3.5)#,label =groupLabelS[i])    	    	
		    x = np.linspace(0,2,100)
		    y = np.ones(len(x))*0.0956
		    ax1.plot(x,y,color='k',ls=':',lw=0.9,label =r'$ \bar{\tilde \rho}_{\rm cor}^{\rm PM} $')	
		    y = np.ones(len(x))*0.0695
		    ax1.plot(x,y,color='red',ls=':',lw=0.9,label =r'$\tilde \rho_{\rm cor}^{\rm fit} $')	    
		    ax1.set_xlabel(r'$ \tilde \Sigma_{\rm cor}$',fontsize=8.) # = \int_0^{\ell_{\rm cor} } \tilde \rho {\rm d} \tilde z 
		    ax1.set_ylabel(r'$  \rho_{\rm cor}  (\sigma^2  / 2 \Sigma_{\rm ref} \lambda_{\rm D})$',fontsize=8.)
		    plt.setp(ax1.get_xticklabels(), fontsize=7.)
		    plt.setp(ax1.get_yticklabels(), fontsize=7.)
	#	    ax1.legend(loc='best',numpoints=1,prop=dict(size=8.),columnspacing=0.10,borderpad=0.10,labelspacing=0.1,ncol=1,handletextpad=0.05,markerscale=0.75,handlelength=1.5)     
		#    ax1.set_xlim(0,MAX*1.1)       
	#	    ax1.set_ylim(0,0.15)
		    ax1.yaxis.set_major_locator(MaxNLocator(5)) 
		    fig.set_size_inches(1.15,0.95)
		    plt.savefig('F?_rhocor_vs_SigCor.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=True,bbox_inches='tight')   
	#	##Biksucks
		    plt.close()
		    


		    print 'Plotting Phi_cor veruss lcor'
		    fig=plt.figure()
		    ax1=fig.add_subplot(111)
		    for (i,col) in enumerate(groupColorS):
		    	plotlcorS = lcorS[i*10:10*(i+1)]
			plotPhiLocalS = np.array(PhiLocalCorS[i*10:10*(i+1)])
		    	ax1.plot(plotlcorS,plotPhiLocalS,color=col,ls='',lw=0.9,marker = markerS[i],label =groupLabelS[i])
		    ax1.set_xlabel(r'$\ell_{\rm cor} / \sigma$',fontsize=10.)
		    ax1.set_ylabel(r'$\Phi_{\rm cor}$',fontsize=10.)
		    plt.setp(ax1.get_xticklabels(), fontsize=8.)
		    plt.setp(ax1.get_yticklabels(), fontsize=8.)
		    ax1.legend(loc='best',numpoints=1,prop=dict(size=8.),columnspacing=0.10,borderpad=0.15,labelspacing=0.1,ncol=2,handletextpad=0.05,markerscale=0.75)     
		    ax1.set_xlim(0,7.5)       
		    ax1.set_ylim(0.,0.75)
		    fig.set_size_inches(3.37,2.5)
#		    plt.savefig('F?_PhiCor_vs_lCor.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
	#	##Biksucks
		    plt.close()

		    print 'Plotting Phi(lcor) veruss ~Sigma'
		    fig=plt.figure()
		    ax1=fig.add_subplot(111)
		    for (i,col) in enumerate(groupColorS):
		    	plotNDSigmaS = NDSigmaS[i*10:10*(i+1)]
			plotPhiLocalS = PhiLocalCorS[i*10:10*(i+1)]
			mark = markerS[i]
		    	ax1.plot(plotNDSigmaS,plotPhiLocalS,color=col,ls='',lw=0.9,marker = mark,label =groupLabelS[i])
		    ax1.set_xlabel(r'$-\Sigma / \Sigma_{\rm ref}$',fontsize=10.)
		    ax1.set_ylabel(r'$\Phi_{\rm cor}$',fontsize=10.)
		    plt.setp(ax1.get_xticklabels(), fontsize=8.)
		    plt.setp(ax1.get_yticklabels(), fontsize=8.)
		    ax1.legend(loc='best',numpoints=1,prop=dict(size=8.),columnspacing=0.10,borderpad=0.15,labelspacing=0.1,ncol=2,handletextpad=0.05,markerscale=0.75)     
		#    ax1.set_xlim(0,MAX*1.1)       
		    ax1.set_ylim(0.,0.75)
		    fig.set_size_inches(3.37,2.5)
#		    plt.savefig('F?_PhiCor_vs_Sigma.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
	#	##Biksucks
		    plt.close()


		    print 'Plotting rho(lcor) versus ~Sigma'
		    fig=plt.figure()
		    ax1=fig.add_subplot(111)
		    for (i,col) in enumerate(groupColorS):
		    	plotNDSigmaS = NDSigmaS[i*10:10*(i+1)]
			plotrhofCorS = rhofCorS[i*10:10*(i+1)]
			mark = markerS[i]
		    	ax1.plot(plotNDSigmaS,plotrhofCorS,color=col,ls='',lw=0.9,marker = mark,label =groupLabelS[i])	    	
		    ax1.set_xlabel(r'$-\Sigma / \Sigma_{\rm ref}$',fontsize=10.)
		    ax1.set_ylabel(r'$  \rho_{\rm cor} (\sigma^2 / 2 \Sigma_{\rm ref} \lambda_{\rm D})$',fontsize=10.)
		    plt.setp(ax1.get_xticklabels(), fontsize=8.)
		    plt.setp(ax1.get_yticklabels(), fontsize=8.)
		    ax1.legend(loc='best',numpoints=1,prop=dict(size=8.),columnspacing=0.10,borderpad=0.15,labelspacing=0.1,ncol=2,handletextpad=0.05,markerscale=0.75)     
		#    ax1.set_xlim(0,MAX*1.1)       
	#	    ax1.set_ylim(-0.25,10.5)
		    fig.set_size_inches(3.37,2.5)
#		    plt.savefig('F?_rhocor_vs_Sig.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
	#	##Biksucks
		    plt.close()




		    print 'Plotting rho(lcor) versus ~lcor'
		    fig=plt.figure()
		    ax1=fig.add_subplot(111)
		    for (i,col) in enumerate(groupColorS):
		    	plotlcorS = lcorS[i*10:10*(i+1)]
			plotrhofCorS = rhofCorS[i*10:10*(i+1)]
			mark = markerS[i]
		    	ax1.plot(plotlcorS,plotrhofCorS,color=col,ls='',lw=0.9,marker = mark,label =groupLabelS[i])
		    	for (x,y) in zip(plotlcorS,plotrhofCorS):
	 			ax1.plot(x,y,color=col,marker=mark)
		    x = np.linspace(0,7,100)
		    y = np.ones(len(x))*0.0695
		    ax1.plot(x,y,color='red',ls=':',lw=0.9,label =r'$\rho_{\rm cor}^{\rm fit} $')	    
		    x = np.linspace(0,7,100)
		    y = np.ones(len(x))*0.0956
		    ax1.plot(x,y,color='k',ls=':',lw=0.9,label =r'$\bar{\rho}_{\rm cor}^{\rm PM} $')	    
		    ax1.set_xlabel(r'$\ell_{\rm cor}/\sigma$',fontsize=10.)
		    ax1.set_ylabel(r'$  \rho_{\rm cor}  (\sigma^2  / 2 \Sigma_{\rm ref} \lambda_{\rm D})$',fontsize=10.)
		    plt.setp(ax1.get_xticklabels(), fontsize=8.)
		    plt.setp(ax1.get_yticklabels(), fontsize=8.)
		    ax1.legend(loc='best',numpoints=1,prop=dict(size=8.),columnspacing=0.10,borderpad=0.15,labelspacing=0.1,ncol=2,handletextpad=0.05,markerscale=0.75)     
	#	    ax1.set_xlim(0,7.5)       
	#	    ax1.set_ylim(0,0.15)
		    fig.set_size_inches(3.37,2.5)
#		    plt.savefig('F?_rhocor_vs_lcor.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
	#	##Biksucks
		    plt.close()
		    
	#	    print 'NEW: Plotting DIFFUSE Sigma versus Volt'
	#	    fig=plt.figure()
	#	    ax1=fig.add_subplot(111)
	#	    for (i,col) in enumerate(groupColorS):
	#	    	plotdiffusezetaS = np.array(VoltCorS[i*10:10*(i+1)]) 
	#	    	plotSigDif = SigDifS[i*10:10*(i+1)]

	#		QV_CS=np.array([[float(x) for x in line.split()] for line in file('CS'+TheoryNames[i],"r").readlines()])
	#   		ax1.errorbar(QV_CS[:,4],QV_CS[:,5],yerr=None,ls='-',lw=1.0,color=col,marker='None')#,label =groupLabelS[i])

	#			##To add Bikerman theory
	##		QV_CS=np.array([[float(x) for x in line.split()] for line in file('Bik'+TheoryNames[i],"r").readlines()])
	##   		ax1.errorbar(QV_CS[:,4],QV_CS[:,5],yerr=None,ls='-.',lw=1.0,color=col,marker='None')#,label =groupLabelS[i])

	#		jj=-1
	#		mark = markerS[i]
	#		ax1.plot(plotdiffusezetaS,plotSigDif,color=col,marker=mark,ls='')
	#	    	ax1.plot([-1,-1],[-1,-1],color = col, ls = '-',lw = 0.9,marker = mark,label =groupLabelS[i])

	#	    ax1.set_ylabel(r'$\Sigma_{\rm dif} / \Sigma_{\rm ref}$',fontsize=10.) # = \int_{\ell_{\rm cor}}^{\infty} \tilde \rho {\rm d} \tilde z
	#	    ax1.set_xlabel(r'$(\phi_{ \ell_{\rm cor} } - \phi_{\rm B}) qe / k_{\rm B}T$',fontsize=10.)
	#	    plt.setp(ax1.get_xticklabels(), fontsize=8.)
	#	    plt.setp(ax1.get_yticklabels(), fontsize=8.)
	#	    ax1.legend(loc='best',numpoints=1,prop=dict(size=8.),columnspacing=0.10,borderpad=0.15,labelspacing=0.1,ncol=2,handletextpad=0.05,markerscale=0.75)     
	#            ax1.set_xlim(0,8.0) 
	#            ax1.set_ylim(0,6.0)
	#	    fig.set_size_inches(3.37,2.5)
	#	    plt.savefig('F?_ChargeVoltage_dif.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
	##	##Biksucks
	#	    plt.close()

	#	    print 'Plotting rho(lcor) versus PhiB'
	#	    fig=plt.figure()
	#	    ax1=fig.add_subplot(111)
	#	    for (i,col) in enumerate(groupColorS):
	#	    	plotPhiB = PhiWCABulkS[i*10:10*(i+1)]
	#	        plotrhofCorS = rhofCorS[i*10:10*(i+1)]
	#		mark = markerS[i]
	#	    	ax1.plot(plotPhiB,plotrhofCorS,color=col,ls='',lw=0.9,marker = mark,label =groupLabelS[i])	    	
	#	    ax1.set_xlabel(r'$\Phi_{\rm B}$',fontsize=10.)
	#	    ax1.set_ylabel(r'$  \rho_{\rm cor} (\sigma^2 / 2 \Sigma_{\rm ref} \lambda_{\rm D})$',fontsize=10.)
	#	    plt.setp(ax1.get_xticklabels(), fontsize=8.)
	#	    plt.setp(ax1.get_yticklabels(), fontsize=8.)
	#	    ax1.legend(loc='best',numpoints=1,prop=dict(size=8.),columnspacing=0.10,borderpad=0.15,labelspacing=0.1,ncol=2,handletextpad=0.05,markerscale=0.75)     
	#	#    ax1.set_xlim(0,MAX*1.1)       
	##	    ax1.set_ylim(-0.25,10.5)
	#	    fig.set_size_inches(3.37,2.5)
	#	    plt.savefig('F?_rhocor_vs_PhiB.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
	##	##Biksucks
	#	    plt.close()

	#	    print 'NEW: Plotting lcor versus SigApp'
	#	    fig=plt.figure()
	#	    ax1=fig.add_subplot(111)
	#	    for (i,col) in enumerate(groupColorS):
	#	    	plotNDSigmaS = NDSigmaS[i*10:10*(i+1)]	    
	#	    	plotlcorS = lcorS[i*10:10*(i+1)]
	#		mark = markerS[i]
	#	    	ax1.plot(plotNDSigmaS,plotlcorS,color=col,ls='-',lw=0.9,marker = mark,label =groupLabelS[i])
	#	    ax1.set_xlabel(r'$-\Sigma / \Sigma_{\rm ref}$',fontsize=10.)
	#	    ax1.set_ylabel(r'$\ell_{\rm cor}/\sigma$',fontsize=10.)
	#	    plt.setp(ax1.get_xticklabels(), fontsize=8.)
	#	    plt.setp(ax1.get_yticklabels(), fontsize=8.)
	#	    ax1.legend(loc='best',numpoints=1,prop=dict(size=8.),columnspacing=0.10,borderpad=0.15,labelspacing=0.1,ncol=2,handletextpad=0.05,markerscale=0.75)     
	#	#    ax1.set_xlim(0,MAX*1.1)       
	#	    ax1.set_ylim(0,7.5)
	#	    fig.set_size_inches(3.37,2.5)
	#	    plt.savefig('F?_lCor_vs_Sigma.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
	##	##Biksucks
	#	    plt.close()

#	    print fail





#    print 'Plotting rho_free and rho_tot...'
#    CarnahanStarling='yes'
##    CarnahanStarling='no'
#    fig=plt.figure()
#    fig.set_size_inches(3.37,4.5)
#    ax1=fig.add_subplot(211)
#    ax2=fig.add_subplot(212)
#    NfS=[]
#    NDSigmaS=[]
#    SigEffS=[]
#    GCs=[]
#    BikS=[]
#    S_CS=[]
#    i=-1
#    Nf_for_EffS=[]
#    VEffS=[]
#    NtS=[]
#    SCS0s=[]
#    zEffS=[]
#    one_time_iterator=0
#    z_forS=[]
#    rhoplusS,rhominusS=[],[]
#    ccc=-1
#    for (Nm,Np,n0,Sigma_s,area,lam_D,L_bin,z_positions,L_z,Bjerrum,filename,sigWCA,V,muex,Npeff,Nmeff,Phi) in zip(originalN_plusS,originalN_minusS,n0s,SIGMAS,areaS,LDs,L_binS,z_originalS,L_zS,BjerrumS,filenameS,sigWCA_S,Volts,muexEV_S,N_plusS,N_minusS,PhiWCAtot_S):
#		i+=1
#		dielectric=(4*np.pi*Bjerrum)**-1
#		
#		Integrand = 0.
#		correction=0.
#		Nf=np.array([(npl-nm)/(area*L_bin) for (npl,nm) in zip(Np,Nm)])
#		Nf=Nf[:len(Nf)/2]
#		z_restore = np.array([(x+y)/(L_z*2.) for (x,y) in zip(z_positions[0:len(z_positions)-1],z_positions[1:len(z_positions)])]) #SO, this was the original z_density required as input to the code below
#		z_pos = z_restore[:len(z_restore)/2]
#		for (y1,y2,z) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos):
#			Integrand=Integrand + 0.5*(y1+y2)*L_bin	  
#			##The code below allows one to integrate out a specific region, i.e. l_corr could be used for this
##	        	if z*L_z<z_wall:
##	        	       	Integrand=Integrand + 0.5*(y1+y2)*L_bin
##	        	       	correction=correction+Integrand
##			else:
##			    	Integrand=Integrand + 0.5*(y1+y2)*L_bin	  
#		Sigma_meas = Integrand  #This must be my reportable Sigma, but NDd
#		eff=(-Sigma_meas+correction)/(dielectric*temperature/(valency*lam_D))
#		ND_Sigma = abs(Sigma_meas)/(dielectric*temperature/(valency*lam_D))
#		NDSigmaS.append(ND_Sigma)

#	        Integrand = eff*(dielectric*temperature/(valency*lam_D))
#	        SigEff=[]
#	        VEff=[]
#	    	volt_correction = Sigma_s/(dielectric*temperature/(valency*lam_D))  - abs(eff)
#	        Nf_for_Eff=[]
#	        zEff = []
#	        rhoplus,rhominus,EVexcess = [],[],[]
#	        for (y1,y2,z,volt,p,m,x) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos,V[:len(Nf)],Npeff[:len(Nf)],Nmeff[:len(Nf)],muex[:len(Nf)]):  #This is the original line of code -- 07/30/13 16:51:36 
#	        		zEff.append(z)
#	        		SigEff.append(Integrand) #There's something wrong with the first value assigned to this function
#	        		test = volt-volt_correction*z
#	        		VEff.append(test)
#	        		rhoplus.append(p)
#	        		rhominus.append(m)
#	        		EVexcess.append(x)
#	        		Integrand=Integrand + 0.5*(y1+y2)*L_bin	
#	        		if filename in GC_LB_LD_5_20:
#	       				Nf_for_Eff.append(0.5*(y1+y2))	        		    		
#	        		else:
#	       				Nf_for_Eff.append(y1)
##			zEffS.append(zEff)
##	        	SigEffS.append(np.array(SigEff)/(dielectric*temperature/(valency*lam_D)))
##	        	VEffS.append(VEff)
##	        	Nf_for_EffS.append(np.array(Nf_for_Eff)/(dielectric*temperature/(valency*lam_D**2)))
##	        	NfS.append(Nf/(dielectric*temperature/(valency*lam_D**2)))

#		rhoplus = np.array(rhoplus)/(area*L_bin*n0)
#		rhominus = np.array(rhominus)/(area*L_bin*n0)
##		EVexcess = np.array(EVexcess) - np.mean(EVexcess[-5:])
#		zEff = np.array(zEff)*L_z


#		ls_x,lw_x='',0.
#		ms_x = 4.0
#		label_str=''
#		alabel=''
#		special_ms='x'
#		edge_width=0.1
##		print 'Zetas might need to change when final measurements come out'
#		if filename == 'Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt':
#			edge_width=0.1
#			ls_x,lw_x=':',0.5
#			special_ms='*'
#			ms_x = 3.5
#			label_str=r'$11.6$'		
#		elif filename == 'Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt':
#			edge_width=0.1
#			ls_x,lw_x=':',0.5
#			special_ms='h'
#			label_str=r'$11.1$'
#		elif 'Analyzed_GC_1312' in filename:
#			edge_width=0.4
#			special_ms='h'
##			label_str=r'$\tilde \Sigma = 7.0$'
#			col = ROYGBIV_map(0.,10)
##			ze = ['6.72','4.66']	
#			if sigWCA==3.0:
#				label_str=r'$-\tilde \Sigma = 7.1$'
#				ze ='6.67'
#			elif sigWCA==4.5:
#				ze ='10.8'
#				label_str=r'$-\tilde \Sigma = 6.8$'
#			elif sigWCA==6.1:
#				ze ='16.23'
#				label_str=r'$-\tilde \Sigma = 6.4$'
#		elif 'Analyzed_GC_1232' in filename:
#			edge_width=0.1
#			special_ms='o'
#			if sigWCA==3.0:
#				label_str=r'$5.4$'
#				ze ='5.14'
#			elif sigWCA==4.5:
#				ze ='7.6'
#				label_str=r'$5.3$'
#			elif sigWCA==6.1:
#				ze ='10.93'
#				label_str=r'$5.1$'
#			special_ms='*'	
#			col = ROYGBIV_map(1.,10)		
#		elif 'Analyzed_GC_1162' in filename:
#			edge_width=0.4
#			special_ms='o'	
#			col = ROYGBIV_map(2.,10)	
##			ze = ['4.06','3.28']
##			ze = ['4.02','3.28']
#			if sigWCA == 3.0:
#				label_str=r'$4.2$'
#				ze = '4.02'				
#			elif sigWCA==4.5:
#				ze ='5.41'
#				label_str=r'$4.1$'
#			elif sigWCA==6.1:
#				ze ='7.08'
#				label_str=r'$3.9$'
#		elif 'Analyzed_GC_1136' in filename:
#			edge_width=0.1
#			special_ms='s'	
##			label_str=r'$3.6$'
#			special_ms='d'		
#			col = ROYGBIV_map(3.,10)
##			ze = ['3.54','2.86']	
#			if sigWCA == 3.0:
#				label_str=r'$3.6$'
#				ze = '3.51'				
#			elif sigWCA==4.5:
#				ze ='4.51'
#				label_str=r'$3.6$'
#			elif sigWCA==6.1:
#				ze ='5.59'
#				label_str=r'$3.5$'
#		elif 'Analyzed_GC_1104' in filename:
#			edge_width=0.4
#			special_ms='s'
#			ms_x=3.5
#			label_str=r'$3.1$'
#			col = ROYGBIV_map(4.,10)
#			Force = 0.21359	
##			ze = ['3.08','2.48']
##			ze = ['3.11','2.48']
#			if sigWCA==3.0:
#				ze = '3.11'
#				label_str=r'$3.1$'
#			elif sigWCA==4.5:
#				ze ='3.81'
#			elif sigWCA==6.1:
#				ze = '4.599'
#		elif 'Analyzed_GC_1066' in filename:
#			edge_width=0.1
#			special_ms='p'	
#			ms_x=3.5	
#			label_str=r'$2.3$'
#			col = ROYGBIV_map(5.,10)
##			ze = ['2.32','2.01']	
#			if sigWCA == 3.0:
#				ze = '2.30'
#			elif sigWCA==4.5:
#				ze ='2.62'
#			elif sigWCA==6.1:
#				ze = '2.98'
#		elif 'Analyzed_GC_1048' in filename:
#			edge_width=0.4
#			special_ms='p'	
#			ms_x=3.5
#			label_str=r'$2.0$'
#			special_ms='v'	
#			col = ROYGBIV_map(6.,10)
##			ze = ['2.01','1.76']	
#			if sigWCA == 3.0:
#				ze = '1.97'
#			elif sigWCA==4.5:
#				ze ='2.13'
#			elif sigWCA==6.1:
#				ze = '2.31'
#		elif 'Analyzed_GC_1034' in filename:
#			edge_width=0.1
#			special_ms='H'
#			ms_x=3.5
#			label_str=r'$1.6$'
#			col = ROYGBIV_map(7.,10)
##			ze = ['1.65','1.44']	
#			if sigWCA==3.0:
#				ze = '1.65'
#			elif sigWCA==4.5:
#				ze ='1.7'
#			elif sigWCA==6.1:
#				ze = '1.74'
#		elif 'Analyzed_GC_1016' in filename:
#			edge_width=0.4
#			special_ms='H'	
#			label_str=r'$1.0$'
#			special_ms='>'	
#			col = ROYGBIV_map(8.,10)	
#			Force = 0.06948		
#			if sigWCA == 3.0:
#				ze = '1.07'
#			elif sigWCA==4.5:
#				ze ='1.05'
#			elif sigWCA==6.1:
#				ze = '1.09'	
#		elif 'Analyzed_GC_1004' in filename:
#			edge_width=0.1
#			special_ms='D'	
#			label_str=r'$0.5$'
#			Force = 0.03368
#			col = ROYGBIV_map(9.,10)
##			ze = ['0.50','0.47']
#			if sigWCA == 3.0:
#				ze = '0.54'
#			elif sigWCA==4.5:
#				ze ='0.4701'
#			elif sigWCA==6.1:
#				ze ='0.49'
#				label_str=r'$0.49$'
#		elif filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1000_0.0_0.1_4.5_500.0_23.73.txt':
#			edge_width=0.4
#			special_ms='<'
#			label_str=r'$0.0$'
#			col = 'k'
#			ze = ['0.00','0.00']
#		else:
##			print '\n\n',filename,Sigma_s,np.max(SIGMAS)
#			col = ROYGBIV_map(abs(Sigma_s),np.max(SIGMAS)*1.1,1)
#		ms_x = 2.5

##		print 'file,bulk vol frac,zeta,Sigma = ',filename,np.mean(Phitot[-5:])*0.736419046,VEff[0],eff
##		print eff


#		rhotot = np.array(rhominus) + np.array(rhoplus)
#		Phitot = (rhotot*(np.pi/6)*sigWCA**3) / (area*L_bin)
#		PhiCS = np.array(Phitot)


#		NDrhofree = -(np.array(rhoplus) - np.array(rhominus))/2


#		if CarnahanStarling=='yes' and ze!='0.00':	
#			  if sigWCA==3.0:
#			  	CS_file='CS_0.037_' + ze + '.txt'
#		  	  elif sigWCA==4.5:
#		  	  	CS_file='CS_0.125_' + ze + '.txt'
#		  	  elif sigWCA==6.1:
#		  	  	CS_file='CS_0.311_' + ze + '.txt'
#			  MMA_CS=[[float(x) for x in line.split()] for line in file(CS_file,"r").readlines()]
#			  MMAz = []
#			  NDMMA_free,NDMMA_tot = [],[]
#			  for line in MMA_CS:
#			  	MMAz.append(line[0])
##			  	Phi = (line[1]+line[2])*n0*(np.pi/6)*sigWCA**3
##			  	PhiWCA.append(Phi)
#			  	free = (-line[1] + line[2])/2
#			  	NDMMA_free.append(free)
#				tot = line[1] + line[2]
#				NDMMA_tot.append(tot)
#				
##	  	  	  PhiCS = np.array(PhiWCA)*0.903042806**3
##			  CS_ex = PhiCS*(8-9*PhiCS+3*PhiCS**2)/(1-PhiCS)**3
#			  MMAz = np.array(MMAz)*lam_D/sigWCA-(1+0.299602925)/sigWCA
##			  MMA_ex = CS_ex - CS_ex[len(CS_ex)-1] 

#		offset = 3. - 0.5*i
#	       	ax1.errorbar(zEff/sigWCA-(1+0.299602925)/sigWCA,np.array(NDrhofree)+offset,yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)
#		ax2.errorbar(zEff/sigWCA-(1+0.299602925)/sigWCA,np.array(rhotot)+offset,yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x)
#		zPM = np.array(zEff)/sigWCA-(1+0.299602925)/sigWCA


#		if sigWCA>=6.1:
#		  	t1 = np.linspace(zPM[0],7,40)[1:].tolist() + np.linspace(7,zPM[-1],25)[1:-1].tolist()
#  		elif sigWCA == 4.5:
#  			print 'heree 4.5'
#		  	t1 = np.linspace(zPM[0],4,15)[1:].tolist() + np.linspace(4,zPM[-1],25)[1:-1].tolist() 
#  		elif sigWCA == 3.0:
#  			print 'heree'
#		  	t1 = np.linspace(zPM[0],3,10)[1:].tolist() + np.linspace(3,zPM[-1],25)[1:-1].tolist() 
#	  	
##		print '\n\n',filename
#	       	if CarnahanStarling=='yes':
#	       		ax1.errorbar(MMAz,np.array(NDMMA_free)+offset,yerr=None,color=col,ls='-',lw=1.0)	
#	       		ax2.errorbar(MMAz,np.array(NDMMA_tot)+offset,yerr=None,color=col,ls='-',lw=1.0)

#       			##This is the beginning of analysis to find lcor via CS -- which isn't the most general

##	    		t1 = np.linspace(MMAz[0],MMAz[-1],20)[1:-1]
##	    		interp_ExFromZ = interpolate.LSQUnivariateSpline(MMAz,MMA_ex,t1,k=5)
##	    		fit_ExFromZ = lambda zpos: interp_ExFromZ(zpos)[0]    		


##			print np.array(SigEff[:10])/(dielectric*temperature/(valency*lam_D))
##	    		for (z,ev,volt,Sig) in zip(zPM,EVexcess,VEff,SigEff):
##	    			if z>25:
##	    				ccc+=1
##	    				break
##				else:
##					dif = ev - fit_ExFromZ(z)
##					if ccc>=0:
##						print dif
##					else:
##						print z,'\t',dif

##					if z <=7.73:
##						print z
##						print 'ev, DEL_vcor, DEL_Sigcor = ',ev,-volt+VEff[0],-Sig/(dielectric*temperature/(valency*lam_D))+eff
##						print 'DIFFUSE Volt, Sigma = ',volt,Sig/(dielectric*temperature/(valency*lam_D)),'\n'

##    ax2.errorbar([-10,-10],[-10,-10],yerr=None,color='k',ls='-',lw=1.0,label = r'${\rm CS}$')
#    ax1.set_ylabel(r'$ \rho / 2 e n^{\rm B} + {\rm offset} $',fontsize=10.) # = (n_+ - n_-) / 2 \rho^{\rm B}
#    ax2.set_ylabel(r'$(n_+ + n_-) / n^{\rm B} + {\rm offset}$',fontsize=10.) 
#    ax2.set_xlabel(r'$(z-\delta_{\rm w}) / \sigma$',fontsize=10.)
#    
#    plt.setp(ax1.get_xticklabels(), fontsize=8.)
#    plt.setp(ax1.get_yticklabels(), fontsize=8.)
#    plt.setp(ax2.get_xticklabels(), fontsize=8.)
#    plt.setp(ax2.get_yticklabels(), fontsize=8.)

#    ax1.set_xlim(-0.5,15)
##    ax1.set_ylim(-0.25,6.) 
#    ax2.set_xlim(-0.5,15)
##	    ax2.set_ylim(0,3.5) 

#    ax1.legend(loc='best',numpoints=1,prop={"size":8},columnspacing=0.07,borderpad=0.18,labelspacing=0.10,handletextpad=0.07,ncol=2,markerscale=1.5) 
##    ax2.legend(loc='best',numpoints=1,prop={"size":10},columnspacing=0.07,borderpad=0.18,labelspacing=0.10,handletextpad=0.07) 
#    ax2.legend(loc='best',numpoints=2,prop={"size":8},columnspacing=0.07,borderpad=0.18,labelspacing=0.10,handletextpad=0.07,ncol=1,markerscale=1.5) 
#    plt.savefig('Rho_FreeAndTot_z.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
##    fig.set_size_inches(3.37,3.5)
##    plt.savefig('F?_totalchempotential_z_BikSucks.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
###    plt.show()
#    plt.close()






	##Computing ND charge densities
    print 'Plotting total chemical potential versus distance...'
    fig=plt.figure()
    fig.set_size_inches(3.37,3.5)
    ax1=fig.add_subplot(211)
    ax2=fig.add_subplot(212)
    NfS=[]
    NDSigmaS=[]
    SigEffS=[]
#    V_corS=[]
#    Max_effS=[]
    GCs=[]
    BikS=[]
    S_CS=[]
    i=-1
    Nf_for_EffS=[]
    VEffS=[]
    NtS=[]
    SCS0s=[]
    zEffS=[]
    one_time_iterator=0
    z_forS=[]
    rhoplusS,rhominusS=[],[]
#    testlcor =   [0.47773365833300002, 0.36662254722199999, 0.25551143611100002, 0.25551143611100002, 0.25551143611100002, 0.144400325, 0.144400325, 0.144400325, 0.0332892138889, 0.0332892138889, 0.0332892138889]
    for (Nm,Np,n0,Sigma_s,area,lam_D,L_bin,z_positions,L_z,Bjerrum,filename,sigWCA,V,muex,Npeff,Nmeff) in zip(originalN_plusS,originalN_minusS,n0s,SIGMAS,areaS,LDs,L_binS,z_originalS,L_zS,BjerrumS,filenameS,sigWCA_S,Volts,muexEV_S,N_plusS,N_minusS):
		i+=1
		dielectric=(4*np.pi*Bjerrum)**-1
		
		Nf=[(npl-nm)/(area*L_bin*dielectric*temperature/(valency*lam_D**2)) for (npl,nm) in zip(Np,Nm)]  #This is consistent with other derivatoins

		if '_GC_' in filename:
			muex=np.array(muex[:-1])
		else:
			muex=np.array(muex)


		Integrand = 0.
		correction=0.
		Nf=np.array([(npl-nm)/(area*L_bin) for (npl,nm) in zip(Np,Nm)])
		Nf=Nf[:len(Nf)/2]
		z_restore = np.array([(x+y)/(L_z*2.) for (x,y) in zip(z_positions[0:len(z_positions)-1],z_positions[1:len(z_positions)])]) #SO, this was the original z_density required as input to the code below
		z_pos = z_restore[:len(z_restore)/2]
		for (y1,y2,z) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos):
			Integrand=Integrand + 0.5*(y1+y2)*L_bin	  
			##The code below allows one to integrate out a specific region, i.e. l_corr could be used for this
#	        	if z*L_z<z_wall:
#	        	       	Integrand=Integrand + 0.5*(y1+y2)*L_bin
#	        	       	correction=correction+Integrand
#			else:
#			    	Integrand=Integrand + 0.5*(y1+y2)*L_bin	  
		Sigma_meas = Integrand  #This must be my reportable Sigma, but NDd
		eff=(-Sigma_meas+correction)/(dielectric*temperature/(valency*lam_D))
		ND_Sigma = abs(Sigma_meas)/(dielectric*temperature/(valency*lam_D))
		NDSigmaS.append(ND_Sigma)

		Integrand = eff*(dielectric*temperature/(valency*lam_D))
		SigEff=[]
		VEff=[]
	    	volt_correction = Sigma_s/(dielectric*temperature/(valency*lam_D))  - abs(eff)
		Nf_for_Eff=[]
		zEff = []
		rhoplus,rhominus,EVexcess = [],[],[]
		for (y1,y2,z,volt,p,m,x) in zip(Nf[0:len(Nf)-1],Nf[1:len(Nf)],z_pos,V[:len(Nf)],Npeff[:len(Nf)],Nmeff[:len(Nf)],muex[:len(Nf)]):  #This is the original line of code -- 07/30/13 16:51:36 
				zEff.append(z)
				SigEff.append(Integrand) #There's something wrong with the first value assigned to this function
				test = volt-volt_correction*z
				VEff.append(test)
				rhoplus.append(p)
				rhominus.append(m)
				EVexcess.append(x)
				Integrand=Integrand + 0.5*(y1+y2)*L_bin	
				if filename in GC_LB_LD_5_20:
	       				Nf_for_Eff.append(0.5*(y1+y2))	        		    		
				else:
	       				Nf_for_Eff.append(y1)
		zEffS.append(zEff)
		SigEffS.append(np.array(SigEff)/(dielectric*temperature/(valency*lam_D)))
		VEffS.append(VEff)
		Nf_for_EffS.append(np.array(Nf_for_Eff)/(dielectric*temperature/(valency*lam_D**2)))
		NfS.append(Nf/(dielectric*temperature/(valency*lam_D**2)))

		rhoplus = np.array(rhoplus)/(area*L_bin*n0)
		rhominus = np.array(rhominus)/(area*L_bin*n0)

		rhoplusS.append(rhoplus)
		rhominusS.append(rhominus)

		muexEV_bulk=np.mean(muex[len(muex)/2-2:len(muex)/2+2])
		muex = np.array(EVexcess) - muexEV_bulk

		ls_x,lw_x='',0.
		ms_x = 4.0
		label_str=''
		alabel=''
		special_ms='x'
		edge_width=0.1
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
		elif filename == 'Analyzed_GC_1312_0.0241789536991_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1312_0.024179_0.1_4.5_500.0_23.73.txt':
			edge_width=0.4
			special_ms='h'
			label_str=r'$\tilde \Sigma = 7.0$'
#			evcorS.append(3.46615)
#			lcorS.append(0.477733658333)
			col = ROYGBIV_map(0.,10)
			ze = ['6.72','4.66']	
		elif filename == 'Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1232_0.018603_0.1_4.5_500.0_23.73.txt':
			edge_width=0.1
			special_ms='o'
			label_str=r'$\tilde \Sigma = 5.4$'
			if sigWCA==4.5:
				label_str=r'$\tilde \Sigma = 5.3$'
#			evcorS.append(2.58685576833)
#			lcorS.append(0.366622547222)
			ze = ['5.16','3.85']
			special_ms='*'	
			col = ROYGBIV_map(1.,10)		
		elif filename == 'Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1162_0.014195_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1162_0.014195_0.1_4.5_500.0_23.73.txt':
			edge_width=0.4
			special_ms='o'	
			label_str=r'$4.2$'
			if sigWCA==4.5:
				label_str=r'$4.1$'
			col = ROYGBIV_map(2.,10)	
			ze = ['4.06','3.28']
		elif filename == 'Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1136_0.012341_0.1_4.5_500.0_23.73.txt':
			edge_width=0.1
			special_ms='s'	
			label_str=r'$3.6$'
#			evcorS.append(1.78295)
#			lcorS.append(0.255511436111)
			special_ms='d'		
			col = ROYGBIV_map(3.,10)
			ze = ['3.54','2.86']	
		elif filename == 'Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1104_0.0106795_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1104_0.0106795_0.1_4.5_500.0_23.73.txt':
			edge_width=0.4
			special_ms='s'
			ms_x=3.5
			label_str=r'$3.1$'
#			evcorS.append(1.47685)
#			lcorS.append(0.255511436111)
			col = ROYGBIV_map(4.,10)
			Force = 0.21359	
			ze = ['3.08','2.48']
		elif filename == 'Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1066_0.0078345_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1066_0.0078345_0.1_4.5_500.0_23.73.txt':
			edge_width=0.1
			special_ms='p'	
			ms_x=3.5	
			label_str=r'$2.3$'
#			evcorS.append(1.01055)
#			lcorS.append(0.144400325)
			col = ROYGBIV_map(5.,10)
			ze = ['2.32','2.01']	
		elif filename == 'Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1048_0.0066065_0.1_4.5_500.0_23.73.txt':
			edge_width=0.4
			special_ms='p'	
			ms_x=3.5
			label_str=r'$2.0$'
#			evcorS.append(0.833384464765)
#			lcorS.append(0.144400325)
			special_ms='v'	
			col = ROYGBIV_map(6.,10)
			ze = ['2.01','1.76']	
		elif filename == 'Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1034_0.005482_0.1_4.5_500.0_23.73.txt':
			edge_width=0.1
			special_ms='H'
			ms_x=3.5
			label_str=r'$1.6$'
#			evcorS.append(0.696928982207)
#			lcorS.append(0.144400325)
			col = ROYGBIV_map(7.,10)
			ze = ['1.65','1.44']	
		elif filename == 'Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1016_0.003474_0.1_4.5_500.0_23.73.txt':
			edge_width=0.4
			special_ms='H'	
			label_str=r'$1.0$'
#			evcorS.append(0.395260565885)
#			lcorS.append(0.0332892138889)
			special_ms='>'	
			col = ROYGBIV_map(8.,10)	
			Force = 0.06948		
			ze = ['1.06','0.99']	
		elif filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1004_0.001684_0.1_0.3_500.0_23.73.txt' or filename=='Analyzed_GC_1004_0.001684_0.1_4.5_500.0_23.73.txt':
			edge_width=0.1
			special_ms='D'	
			label_str=r'$0.5$'
			Force = 0.03368
#			evcorS.append(0.316165923552)
#			lcorS.append(0.0332892138889)
			col = ROYGBIV_map(9.,10)
			ze = ['0.50','0.47']
		elif filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1000_0.0_0.1_4.5_500.0_23.73.txt':
			edge_width=0.4
			special_ms='<'
			label_str=r'$0.0$'
			col = 'k'
			ze = ['0.00','0.00']
		else:
			col = ROYGBIV_map(abs(Sigma_s),np.max(SIGMAS)*1.1,1)
		ms_x = 3.5

		mup = np.log(rhoplus) + np.array(VEff) + np.array(muex)+i
		mum = np.log(rhominus) - np.array(VEff) + np.array(muex)+i

		zEff = np.array(zEff)*L_z		

#		print 'Trying out this for NoCoul'
#		print np.log(rhominus[:10])
#        	mup = np.log(rhoplus*n0) - Force *zEff  + np.array(muex)
		
#        	mum = np.log(rhominus*n0) #+ Force *zEff + np.array(muex)
#        	mup = np.log(rhominus*n0) + Force *zEff + np.array(muex)

####		zEff/sigWCA-(1+0.299602925)/sigWCA
#		print zEff[:5],'\n',zEff[:5]/sigWCA
#		print zEff/sigWCA-(1+0.299602925)/sigWCA

		ax1.errorbar(zEff/sigWCA-(1+0.299602925)/sigWCA,mup,yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)	
		ax2.errorbar(zEff/sigWCA-(1+0.299602925)/sigWCA,mum,yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)	        	        	
			        	

    ax1.set_ylabel(r'$(\mu_{\rm co} - \mu_{\rm Bulk}) / k_{\rm B} T$',fontsize=10.)
    ax2.set_ylabel(r'$(\mu_{\rm counter} - \mu_{\rm Bulk}) / k_{\rm B} T$',fontsize=10.)
    ax2.set_xlabel(r'$z / \sigma$',fontsize=10.)
    
    plt.setp(ax1.get_xticklabels(), fontsize=8.)
    plt.setp(ax1.get_yticklabels(), fontsize=8.)
    plt.setp(ax2.get_xticklabels(), fontsize=8.)
    plt.setp(ax2.get_yticklabels(), fontsize=8.)

#    plt.savefig('F?_totalchempotential_z_BikSucks_nolims.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      

    ax1.set_xlim(-0.25,15)
    ax1.set_ylim(-0.05,2.0) 
    ax2.set_xlim(-0.25,15)
    ax2.set_ylim(-0.05,2.0)     

#    ax1.legend(loc='best',numpoints=1,prop={"size":10},columnspacing=0.07,borderpad=0.18,labelspacing=0.10,handletextpad=0.07) 
    fig.set_size_inches(3.37,3.5)
#    plt.savefig('F?_totalchempotential_z_BikSucks.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
#    plt.show()
    plt.close()

    print 'here now'
    print NDSigmaS
    print np.min(NDSigmaS)
    print np.max(NDSigmaS)


    ratioS,legend_key=np.array(BjerrumS)/np.array(LDs),'ele'

###    ##This plots mu_ex^EV(z) for simulation and modified theories
    muexEV_bulkS=[]
    muexEV_wallS,PhiWallS=[],[]
    muexEV_midS,PhimidS=[],[]
    colmark_wallS,colmark_midS = [],[]
    theorycolS=[]
    nonDim='yes'
#    Bikerman='yes'
#    Bikerman='no'
    CarnahanStarling='yes'
    CarnahanStarling='no'
    fromlcor = 'yes'
    fromlcor = 'no'
    if nonDim=='yes':
      print 'Plotting NonDim mu_ex^EV(z)...'
    else:
      print 'Plotting Dim mu_ex^EV(z)...'	
    if Bikerman=='yes':
		print '\t\t...with Bikerman theory...'
		#Note to user: Bikerman theory must already be evaluated using MATLAB code, which is in ~/sims/PBik_Solver.m
    if CarnahanStarling=='yes':
	      print '\t...with Carnahan-Starling theory...'
    i=-1
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    #BikSucks
    cc=-1
    ax1.errorbar([-1,-2],[-1,-2],yerr=None,color='k',ls='-',lw=1.0,label='  '+ r'$\tilde \mu^{\rm CS}(0.95 \cdot \sigma)$') #
#    ax1.errorbar([-1,-2],[-1,-2],yerr=None,color='k',ls='-.',lw=1.0,label='  '+ r'$\tilde \mu^{\rm Bik}$')
    evcorS,lcorS=[],[]
    for (z_positions,characteristic_length,muex_EV,L_z,sigWCA,PhiWCA,lamD,lamB,filename,n0) in zip(z_originalS,characteristic_lengthS,muexEV_S,L_zS,sigWCA_S,PhiWCAtot_S,LDs,BjerrumS,filenameS,n0s):
    	i+=1
	cc+=1
	if i==len(markers): #This resets the markers
		i=0

	col = ROYGBIV_map((sigWCA)/lamD,3,1) #*0.954028946

	rat = (sigWCA/lamB)
	if rat <= 0.05:
		mark = 'o'
	elif rat>0.05 and rat<=0.2:
		mark = 's'
	elif rat>0.2 and rat<=1.:
		mark = 'h'
	elif rat>1. and rat<=10.:
		mark = 'p'
	elif rat>10. and rat<=120.:
		mark = 'D'
#	elif rat>1.5:
#		mark = 'D'
	msize = 5.0	

	for (zz,mu,phi) in zip(z_positions,muex_EV,PhiWCA):
		if zz>=0.5*sigWCA and zz<=3*sigWCA: #This used to be 0<value < 3sigWCA
			muexEV_wallS.append(mu)
			PhiWallS.append(phi)#*0.954028946**3)  #Leave as is, scaling happens in the theory
			col_mark = [col,mark]
			colmark_wallS.append(col_mark)
		elif zz>3*sigWCA and zz<=6*sigWCA:
			muexEV_midS.append(mu)
			PhimidS.append(phi)#*0.954028946**3)
			if phi>1:
				print 'Warning Phi>1 for:  ',filename
			col_mark = [col,mark]
			colmark_midS.append(col_mark)

	ls_x,lw_x='',0.
	ms_x = 4.0
	label_str=''
	alabel=''
	special_ms='x'
	edge_width=0.1
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
		label_str=r'$\tilde \Sigma = 7.0$'
		evcorS.append(3.46615)
		lcorS.append(0.477733658333)
		col = ROYGBIV_map(0.,10)	
	elif filename == 'Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='o'
		label_str=r'$5.4$'
		evcorS.append(2.58685576833)
		lcorS.append(0.366622547222)
		special_ms='*'	
		col = ROYGBIV_map(1.,10)		
	elif filename == 'Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1162_0.014195_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='o'	
		label_str=r'$4.2$'
		evcorS.append(2.13433)
		lcorS.append(0.255511436111)
		col = ROYGBIV_map(2.,10)	
	elif filename == 'Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='s'	
		label_str=r'$3.6$'
		evcorS.append(1.78295)
		lcorS.append(0.255511436111)
		special_ms='d'		
		col = ROYGBIV_map(3.,10)	
	elif filename == 'Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1104_0.0106795_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='s'
		ms_x=3.5
		label_str=r'$3.1$'
		evcorS.append(1.47685)
		lcorS.append(0.255511436111)
		col = ROYGBIV_map(4.,10)	
	elif filename == 'Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1066_0.0078345_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='p'	
		ms_x=3.5	
		label_str=r'$2.3$'
		evcorS.append(1.01055)
		lcorS.append(0.144400325)
		col = ROYGBIV_map(5.,10)	
	elif filename == 'Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='p'	
		ms_x=3.5
		label_str=r'$2.0$'
		evcorS.append(0.833384464765)
		lcorS.append(0.144400325)
		special_ms='v'	
		col = ROYGBIV_map(6.,10)	
	elif filename == 'Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='H'
		ms_x=3.5
		label_str=r'$1.6$'
		evcorS.append(0.696928982207)
		lcorS.append(0.144400325)
		col = ROYGBIV_map(7.,10)	
	elif filename == 'Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='H'	
		label_str=r'$1.0$'
		evcorS.append(0.395260565885)
		lcorS.append(0.0332892138889)
		special_ms='>'	
		col = ROYGBIV_map(8.,10)	
	elif filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='D'	
		label_str=r'$0.5$'
		evcorS.append(0.316165923552)
		lcorS.append(0.0332892138889)
		col = ROYGBIV_map(9.,10)	
	elif filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='<'
		label_str=r'$0.0$'
		evcorS.append(0.281650119324)
		lcorS.append(0.0332892138889)
	ms_x = 4.0
		
	z_positions=[z/L_z for z in z_positions] #NonDim distances

        muexEV_bulk=np.mean(muex_EV[len(muex_EV)/2-2:len(muex_EV)/2+2])
        muexEV_bulkS.append(muexEV_bulk)

	ax1.errorbar(np.array(z_positions[:len(muex_EV)/2])*(L_z/sigWCA)-(1+0.299602925)/sigWCA,np.array(muex_EV[:len(muex_EV)/2])-muexEV_bulk,yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)
	theorycolS.append(col)

#    if sys.argv[1]=='F6_ex.txt':
##	ax1.errorbar(lcorS,evcorS,yerr=None,color='magenta')
#    	evcorS = np.array(evcorS) + np.array(muexEV_bulkS)
##    	print 'post correction = ',evcorS
#    	ax1.errorbar(lcorS,evcorS,yerr=None,color='k')    	
    
    cc=-1
    i=-1
    for (z_positions,characteristic_length,muex_EV,L_z,sigWCA,PhiWCA,lamD,lamB,filename,n0,col) in zip(z_originalS,characteristic_lengthS,muexEV_S,L_zS,sigWCA_S,PhiWCAtot_S,LDs,BjerrumS,filenameS,n0s,theorycolS):
    	i+=1
	cc+=1

	ls_x,lw_x='',0.
	ms_x = 3.0
	label_str=''
	alabel=''
	special_ms='x'
	edge_width=0.1
	if filename == 'Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt':
		label_str=r'$11.6$'		
	elif filename == 'Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt':
		label_str=r'$11.1$'
	elif filename == 'Analyzed_GC_1312_0.0241789536991_0.1_3.0_500.0_23.73.txt':
		label_str=r'$\tilde \Sigma = 6.0$'
		ze = ['6.72','4.66']
	elif filename == 'Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt':
		label_str=r'$5.0$'
		ze = ['5.16','3.85']
	elif filename == 'Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1162_0.014195_0.1_3.0_500.0_23.73.txt':
		label_str=r'$3.9$'
		ze = ['4.06','3.28']
	elif filename == 'Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt':
		label_str=r'$3.4$'
		ze = ['3.54','2.86']
	elif filename == 'Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1104_0.0106795_0.1_3.0_500.0_23.73.txt':
		label_str=r'$3.0$'
		ze = ['3.08','2.48']
	elif filename == 'Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1066_0.0078345_0.1_3.0_500.0_23.73.txt':
		label_str=r'$2.2$'
		ze = ['2.32','2.01']
	elif filename == 'Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt':	
		label_str=r'$1.9$'
		ze = ['2.01','1.76']
	elif filename == 'Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt':
		label_str=r'$1.6$'
		ze = ['1.65','1.44']
	elif filename == 'Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt':
		label_str=r'$1.0$'
		ze = ['1.06','0.99']
	elif filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt':
		label_str=r'$0.5$'
		ze = ['0.50','0.47']
	elif filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt':
		label_str=r'$0.0$'
		ze = ['0.00','0.00']

	if Bikerman=='yes' and ze!='0.00':
		##See modifications to CS code below, if Bikerman is ever needed in these plots - 09/27/13 08:24:24 
	  Bik_file='Bik_0.0438_' + ze + '.txt'
          MMA_Bik=[[float(x) for x in line.split()] for line in file(Bik_file,"r").readlines()]
          MMAz,PhiBik=[],[]
	  for line in MMA_Bik:
	  	MMAz.append(line[0])
	  	Phi = (line[1]+line[2])*n0*(np.pi/6)*sigWCA**3
	  	PhiBik.append(Phi) 
	  PhiBik = np.array(PhiBik)
	  Bik_ex = -np.log(1 - PhiBik)
	  Bik_ex = Bik_ex - Bik_ex[len(Bik_ex)-1]   
	  ax1.plot(np.array(MMAz)*lamD/sigWCA,Bik_ex,color='k',lw=1.5,ls='-')    
	  ax1.plot(np.array(MMAz)*lamD/sigWCA,Bik_ex,color=col,lw=1.0,ls='-.')    
	if CarnahanStarling=='yes':
		if ze[0]!='0.00':
			  if fromlcor == 'yes': #This is plotting theory from lcor onwards
			    	ze = ze[1]
			    	xshift = lcorS[i]
			  else:
				ze = ze[0]	
				xshift = -(1+0.299602925)/sigWCA	
			  CS_file='CS_0.0438_' + ze + '.txt'
#			  print filename,CS_file
			  MMA_CS=[[float(x) for x in line.split()] for line in file(CS_file,"r").readlines()]
			  MMAz,PhiWCA = [],[]
			  for line in MMA_CS:
			  	MMAz.append(line[0])
			  	Phi = (line[1]+line[2])*n0*(np.pi/6)*sigWCA**3
			  	PhiWCA.append(Phi)
			  PhiCS = np.array(PhiWCA)*0.954028946**3
			  CS_ex = PhiCS*(8-9*PhiCS+3*PhiCS**2)/(1-PhiCS)**3
			  CS_ex = CS_ex - CS_ex[len(CS_ex)-1]   
			  xplot = np.array(MMAz)*lamD/sigWCA + xshift
			  ax1.plot(xplot,CS_ex,color=col,lw=1.0,ls='-') ##Old
#	  else:
#	  		##Save this code, but it doens't actually work... - 09/30/13 13:28:31 
#		  CS_file='CS_0.0438_6.72.txt'
#		  MMA_CS=[[float(x) for x in line.split()] for line in file(CS_file,"r").readlines()]
#		  MMAz,MMAnp,MMAnm,MMAfree,MMAvolt,MMASig = [],[],[],[],[],[]
#		  for line in MMA_CS:
#		  	MMAz.append(line[0])
#		  	MMAnp.append(line[1])
#		  	MMAnm.append(line[2])
#		  	MMAfree.append(line[3])
#		  	MMAvolt.append(line[4])
#		  	MMASig.append(line[5])
#		  	##There's other stuff buried in here like S and Capacitance
#		  interp_z,interp_ntot = CSInterpolation_ntot(MMAz,MMAnp,MMAnm,MMAfree,MMAvolt,MMASig,float(ze))
#		  PhiWCA = np.array(interp_ntot)*n0*(np.pi/6)*sigWCA**3
#		  PhiCS = np.array(PhiWCA)*0.954028946**3
#		  CS_ex = PhiCS*(8-9*PhiCS+3*PhiCS**2)/(1-PhiCS)**3
#		  CS_ex = CS_ex - CS_ex[len(CS_ex)-1]   
#		  xplot = np.array(interp_z)*lamD/sigWCA-(1+0.299602925)/sigWCA
#		  ax1.plot(xplot,CS_ex,color=col,lw=1.0,ls='--') 
    ax1.set_xlabel(r'$z / \sigma$',fontsize=10.)
    ax1.set_ylabel(r'$(\mu^{\rm EV}[\tilde z] - \mu^{\rm EV}[\tilde z^{\rm B}]) / k_{\rm B} T$',fontsize=10.)
    plt.setp(ax1.get_xticklabels(), fontsize=8.)
    plt.setp(ax1.get_yticklabels(), fontsize=8.)

#    plt.savefig('F4_exV_vs_z_BikSucks_nolims.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      

    ax1.set_xlim(-1/sigWCA*1.5,10)            
    ax1.set_ylim(-0.25,6)

#    ax1.legend(loc='best',numpoints=1,prop={"size":10},columnspacing=0.07,borderpad=0.18,labelspacing=0.10,handletextpad=0.07) 
    fig.set_size_inches(3.37,3.5)
#    if fromlcor == 'yes': #This is plotting theory from lcor onwards
#	    plt.savefig('F4_exV_vs_z_BikSucks_lcor.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
#    else:
#	    plt.savefig('F4_exV_vs_z_BikSucks_zmin.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
#    plt.show()
    plt.close()
#    #Done plotting mu_excess^EV(z)


	
    print 'Max sig/lamD = ',np.max(np.array(sigWCA_S)/np.array(LDs))
    print 'Min sig/lamD = ',np.min(np.array(sigWCA_S)/np.array(LDs))

    print 'Plotting muexEV = muexEV(Phi_bulk)...'
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    i=-1
    j=-1
    PhiWCAs=[]
    Phi_shifted=[]
    Phi_HS=[]
    MAX = -1e9
    colS,sigratS=[],[]
    markS = []
    LBrat=[]
    for (mexb,n0,sig_WCA,filename,lamD,lamB) in zip(muexEV_bulkS,n0s,sigWCA_S,filenameS,LDs,BjerrumS):
    	i+=1
    	j+=1
	if i==len(markers): #This resets the markers
			i=0
    	Phi_WCA = 2.*n0*(np.pi/6)*sig_WCA**3
    	MAX = np.max([MAX,Phi_WCA])
    	PhiWCAs.append(Phi_WCA)
    	Phi_shifted.append(Phi_WCA*0.954028946**3)  #Unweighted = 0.981989165  0.954028946 

	col = ROYGBIV_map(sig_WCA/lamD,3.0,1)
		#Used to be col = ROYGBIV_map((sig_WCA)/lamD,3.0,1) # for some reason
	colS.append(col)
	sigratS.append(sig_WCA/lamD)

	rat = (sigWCA/lamB)
	LBrat.append(rat)
	if rat <= 0.10:
		mark = 'o'
	elif rat>0.10 and rat<=0.2:
		mark = 's'
	elif rat>0.2 and rat<=1.:
		mark = 'h'
	elif rat>1. and rat<=2.:
		mark = 'p'
	elif rat>2.:
		mark = 'D'

	ax1.errorbar(Phi_WCA,mexb,yerr=None,marker=mark,mew=0.1,ms=msize,color=col,ls='None')	
	markS.append(mark)
    Phi_theory = np.linspace(0,0.9,2000)

    print 'MAX and MIN sig/lamD'
    print np.max(sigratS),np.min(sigratS)
    print 'MAX and MIN sig/lamB'
    print np.max(LBrat),np.min(LBrat)    
    LBrat.sort()
    print len(LBrat),LBrat
    HSPhi_theory = Phi_theory*0.903042806**3  
    print "NEEED TO WAIT FOR ALL RUNS BEFORE DOING THIS FOR REAL REAL"  
    ax1.errorbar(Phi_theory,Phi_theory*(8-9*Phi_theory+3*Phi_theory**2)/(1-Phi_theory)**3,yerr=None,color='red',ls='-',lw=1.3,label=r'$\mu^{\rm CS}_{\rm ex}(\sigma)$')    
    ax1.errorbar(Phi_theory,HSPhi_theory*(8-9*HSPhi_theory+3*HSPhi_theory**2)/(1-HSPhi_theory)**3,yerr=None,color='k',ls='-',lw=1.3,label=r'$\mu^{\rm CS}_{\rm ex}(\sigma_{\rm eff})$')


#    cbar_x = np.linspace(0.1,0.3,1500)
#    cbar_val = np.linspace(0.001,2.5,1500)
#    for (x,c) in zip(cbar_x,cbar_val):
#    	col = ROYGBIV_map(c,3.0,1)
#    	ax1.plot(x,7.0,marker='|',ms=5.5,color=col)
#    	ax1.plot(x,7.05,marker='|',ms=5.5,color=col)
#    	ax1.plot(x,7.07,marker='|',ms=5.5,color=col)
#    	ax1.plot(x,7.09,marker='|',ms=5.5,color=col)
#    	ax1.plot(x,7.11,marker='|',ms=5.5,color=col)
#    	ax1.plot(x,7.13,marker='|',ms=5.5,color=col)
#    	ax1.plot(x,7.15,marker='|',ms=5.5,color=col)

#    ax1.errorbar(Phi_theory,-np.log(1-Phi_theory),yerr=None,color='k',ls='-.',label=r'$\mu^{\rm Bik}=-\ln(1-\Phi) $') #= -\ln (1-\Phi/0.65 )

    ax1.errorbar(Phi_theory,-np.log(1-Phi_theory),yerr=None,color='red',ls='-.',label=r'$\mu^{\rm Bik}_{\rm ex}(\sigma) $') #= -\ln (1-\Phi/0.65 )
    ax1.errorbar(Phi_theory,-np.log(1-Phi_theory/0.484),yerr=None,color='k',ls='-.',label=r'$\mu^{\rm Bik}_{\rm ex}( \sigma_{\rm fit})$') #= -\ln (1-\Phi/0.65 )

    WCA_file='WCAex.txt'
    WCAex=[[float(x) for x in line.split()] for line in file(WCA_file,"r").readlines()]
    WCAPhi,exWCA = [],[]
    for line in WCAex:
  	WCAPhi.append(line[0])
  	exWCA.append(line[1])  
    ax1.plot(WCAPhi,exWCA,color='0.90',lw=1.3,ls='--',label=r'$\mu^{\rm WCA}_{\rm ex}$')  

#    ax1.errorbar(-1,-1,yerr=None,marker='',mew=0.1,ms=3.5,color='k',ls='None',label =  r'$   $')
      
#    labelS = [r'$ \sigma / \lambda_{\rm B} = 0.10$', r'$ 0.10 < \tilde \sigma \leq 0.2 $', r'$ 0.2 < \tilde \sigma \leq 1$', r'$ 1 < \tilde \sigma \leq 2$', r'$ 2 < \tilde \sigma \leq 10$', r'$   $', r'$   $']    	
    labelS = [r'$ \sigma / \lambda_{\rm B} = 0.10$', r'$ 0.10 < \check \sigma \leq 0.2 $', r'$ 0.2 < \check \sigma \leq 1$', r'$ 1 < \check \sigma \leq 2$', r'$ 2 < \check \sigma \leq 10$'] 
    for (mark,leg) in zip(['o','s','h','p','D'],labelS):
        ax1.errorbar(-1,-1,yerr=None,marker=mark,mew=0.1,ms=3.5,color='k',ls='None',label = leg)
    ax1.legend(loc='upper left',numpoints=1,prop=dict(size=8.),columnspacing=0.07,borderpad=0.30,labelspacing=0.08,ncol=2,handletextpad=0.05)#,handlelength=1.85)


#	    ax1.legend(loc='best',numpoints=1,prop=dict(size=8.),columnspacing=0.10,borderpad=0.15,labelspacing=0.1,ncol=2,handletextpad=0.05,markerscale=0.75)     

    ax1.set_xlabel(r'$\Phi_{\rm B} = \sigma^3 / \left ( 24 \lambda_{\rm B} \lambda_{\rm D}^2 \right ) $',fontsize=10.)
    ax1.set_ylabel(r'$\mu^{\rm EV}_{\rm ex,B} / k_{\rm B} T$',fontsize=10.)
    plt.setp(ax1.get_xticklabels(), fontsize=8.)
    plt.setp(ax1.get_yticklabels(), fontsize=8.)

#    ax1.set_xlim(0,MAX*1.1)       
#    ax1.set_ylim(0,np.max(muexEV_bulkS)*1.1)

    ax1.set_xlim(0,0.6)       
    ax1.set_ylim(0,9)

    fig.set_size_inches(3.37,2.5)
    plt.savefig('F1_BikSucks_Bulk.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=True,bbox_inches='tight')   
#    ax1.set_xlim(0,0.05805)
#    ax1.set_ylim(0,0.5) 
#    plt.savefig('muexEV_vs_Phi_B_LT1.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False) 
    plt.close()


#    print fail

	##ALL PEFECTLY FINE CODE, commented out for speed up
####    colmark_wallS,colmark_midS = [],[]

    print 'Plotting muexEV = muexEV(Phi_mid)...'
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    MAX = -1e9
    for (mexm,Phimid,colmark) in zip(muexEV_midS,PhimidS,colmark_midS):
    	MAX = np.max([MAX,Phimid])
	ax1.errorbar(Phimid,mexm,yerr=None,marker=colmark[1],mew=0.1,ms=3.5,color=colmark[0],ls='None')	
    HSPhi_theory = Phi_theory*0.903042806**3
    ax1.errorbar(Phi_theory,HSPhi_theory*(8-9*HSPhi_theory+3*HSPhi_theory**2)/(1-HSPhi_theory)**3,yerr=None,color='k',ls='-',lw=1.3,label=r'$\mu^{\rm CS}(\sigma_{\rm eff})$')
    Phi_theory = np.linspace(0,1.0,1000)
#    ax1.errorbar(Phi_theory,-np.log(1-Phi_theory),yerr=None,color='red',ls='-.',label=r'$\mu^{\rm Bik}(\sigma)$') #= -\ln (1-\Phi/0.65 )
#    ax1.errorbar(Phi_theory,Phi_theory*(8-9*Phi_theory+3*Phi_theory**2)/(1-Phi_theory)**3,yerr=None,color='red',ls='-',lw=1.3,label=r'$\mu^{\rm CS}(\sigma)$')    
#    ax1.errorbar(Phi_theory,-np.log(1-Phi_theory/0.484),yerr=None,color='green',ls='-.',label=r'$\mu^{\rm Bik}_{\rm fit}(1.27 \cdot \sigma)$') #= -\ln (1-\Phi/0.65 )
    ax1.set_xlabel(r'$\Phi \left [ 3  < z/\sigma < 6 \right ] $',fontsize=10.)
    ax1.set_ylabel(r'$\mu^{\rm EV}_{\rm ex} \left [ 3  < z/\sigma < 6 \right ] / k_{\rm B} T$',fontsize=10.)
    plt.setp(ax1.get_xticklabels(), fontsize=8.)
    plt.setp(ax1.get_yticklabels(), fontsize=8.)
#    ax1.legend(loc='best',numpoints=1,prop=dict(size=8.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
    ax1.set_xlim(0,1.0)       
    ax1.set_ylim(0,20)
    fig.set_size_inches(3.37,2.25)
    plt.savefig('F1_BikSucks_mid.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=True,bbox_inches='tight')   

##Biksucks
    plt.close()


    print 'Plotting muexEV = muexEV(Phi_wall)...'
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    MAX = -1e9
    for (mexw,Phiwall,colmark) in zip(muexEV_wallS,PhiWallS,colmark_wallS):
    	MAX = np.max([MAX,Phiwall])
    	if Phiwall<1.0:
    		ax1.errorbar(Phiwall,mexw,yerr=None,marker=colmark[1],mew=0.1,ms=3.5,color=colmark[0],ls='None')	
    Phi_theory = np.linspace(0,1.0,1000)
#    HSPhi_theory = Phi_theory*0.954028946**3
    HSPhi_theory = Phi_theory*0.903042806**3
    ax1.errorbar(Phi_theory,HSPhi_theory*(8-9*HSPhi_theory+3*HSPhi_theory**2)/(1-HSPhi_theory)**3,yerr=None,color='k',ls='-',lw=1.3,label='  '+r'$\mu^{\rm CS}(\sigma_{\rm eff})$')
#    ax1.errorbar(Phi_theory,-np.log(1-Phi_theory),yerr=None,color='red',ls='-.',label='  '+r'$\mu^{\rm Bik}  $') #= -\ln (1-\Phi/0.65 )

    ax1.set_xlabel(r'$\Phi \left [ z < 3 \sigma \right ] $',fontsize=10.)
    ax1.set_ylabel(r'$\mu^{\rm EV}_{\rm ex} \left [ z < 3 \sigma \right ] / k_{\rm B} T$',fontsize=10.)
    plt.setp(ax1.get_xticklabels(), fontsize=8.)
    plt.setp(ax1.get_yticklabels(), fontsize=8.)

    ax1.set_xlim(0,1.0)       
    ax1.set_ylim(0,20)
    fig.set_size_inches(3.37,2.25)
    plt.savefig('F1_BikSucks_wall.png', dpi=1000, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=True,bbox_inches='tight')   
##Biksucks
    plt.close()

#    print fail
#    print 'Plotting Excess minus Bulk...'
#    fig=plt.figure()
#    fig.set_size_inches(3.37,3.5)
#    ax1=fig.add_subplot(211)
#    ax2=fig.add_subplot(212)
#    PM_Phi_Ex_file='PM_Phi_Ex.txt'
#    PhiEx=[[float(x) for x in line.split()] for line in file(PM_Phi_Ex_file,"r").readlines()]
#    xPhi,yEx = [],[]
#    i=-1
#    for x in PhiEx[0]:
#    	i+=1
#	if (i % 2)==0:
#		xPhi.append(x)
#	else:
#  		yEx.append(x)		
#    t1 = np.linspace(xPhi[0],xPhi[-1],20)[1:-1]
#    interp_Ex_Phi = interpolate.LSQUnivariateSpline(xPhi,yEx,t1,k=5)
#    fit_Ex_Phi = lambda Phi: interp_Ex_Phi(Phi)[0]
##    print 'Max allowable Phi to interp function =',max(xPhi)

#    i=-1
#    bulknoise_fit = []
#    for (Evz,PhiWCAz,sigWCA,lamB,zpos,filename,L_z,ndSigma_s,ndSigeff,V,PhiB) in zip(muexEV_S,PhiWCAtot_S,sigWCA_S,BjerrumS,z_S,filenameS,L_zS,NDSigmaS,SigEffS,VEffS,PhiBulkS):
#	if '_GC_' in filename:
#		Evz=np.array(Evz[:-1])
#	else:
#		Evz=np.array(Evz)
#    	zpos=np.array(zpos[:-1])
#    	PhiWCAz = np.array(PhiWCAz)

##    	if PhiB<0.1:
##    		print filename

#	fitEV = np.array([fit_Ex_Phi(Phi) for Phi in PhiWCAz]) #Can't use this because of error...
#    	fit_difference = -Evz + fitEV

##	PhiCS = PhiWCAz*0.954028946**3
#	PhiCS = PhiWCAz*0.903042806**3	
#	CSEV = PhiCS*(8-9*PhiCS+3*PhiCS**2)/(1-PhiCS)**3
#    	CS_difference = -Evz + CSEV    	

#	CSdifnoise,fitdifnoise = [],[]
#	difference = fit_difference
#	for (z,ev,phi) in zip(zpos,Evz,PhiWCAz):
#		if z/L_z<=0.5 and z/L_z>=0.3:
#			fitEV = fit_Ex_Phi(phi)
##			PhiCS = phi*0.954028946**3
#			PhiCS = phi*0.903042806**3		
#			CSEV = PhiCS*(8-9*PhiCS+3*PhiCS**2)/(1-PhiCS)**3

#			CSdifnoise.append(-ev+CSEV)
#			fitdifnoise.append(-ev+fitEV) 

#	tol = np.max(abs(np.array(fitdifnoise)))
#	if tol<0.1:
#		bulknoise_fit.append(tol)	
#		

#	difference = fit_difference
#	dif_plot,z_plot,Phi_plot = [],[],[]
#	for (dif,z,phi) in zip(difference,zpos,PhiWCAz):		
#		if ((z/zpos[-1])<=0.5):
#			dif_plot.append(dif)
#			z_plot.append(z)
#			Phi_plot.append(phi)
#	rat = (sigWCA/lamB)
#	if rat <= 0.05:
#		mark = 'o'
#	elif rat>0.05 and rat<=0.2:
#		mark = 's'
#	elif rat>0.2 and rat<=1.:
#		mark = 'h'
#	elif rat>1. and rat<=10.:
#		mark = 'p'
#	elif rat>10. and rat<=120.:
#		mark = 'D'

#	if len(filenameS)==446:
#		col = ROYGBIV_map(abs(ndSigma_s),98.28,1) #98.28=np.max(corsig)*1.2
#	else:
#		col = ROYGBIV_map(abs(ndSigma_s),np.max(NDSigmaS)*1.1,1)
#	##Below is the code to use -- 09/14/13 12:01:51 
###    	ax1.errorbar(-np.array(dif_plot),np.array(z_plot)/sigWCA,yerr=None,marker=mark,mew=0.1,ms=3.5,color=col,ls='None')
##    	ax1.errorbar(-np.array(dif_plot),np.array(z_plot),yerr=None,marker=mark,mew=0.1,ms=3.5,color=col,ls='None')
##    	ax2.errorbar(-np.array(dif_plot),Phi_plot,yerr=None,marker=mark,mew=0.1,ms=3.5,color=col,ls='None')

##    ax1.set_ylabel(r'$ z / \sigma$',fontsize=10.)
#    ax1.set_ylabel(r'$ z $',fontsize=10.)

#    ax2.set_ylabel(r'$\Phi[\tilde z]$',fontsize=10.)
#    ax2.set_xlabel(r'$ (\mu^{\rm EV}[\Phi] - \mu^{\rm EV}_{\rm Bulk}[\Phi]) / k_{\rm B} T$',fontsize=10.)

#    plt.setp(ax1.get_xticklabels(), fontsize=8.)
#    plt.setp(ax1.get_yticklabels(), fontsize=8.)
#    plt.setp(ax2.get_xticklabels(), fontsize=8.)
#    plt.setp(ax2.get_yticklabels(), fontsize=8.)

#    ax1.set_xlim(-0.75,0.75)    
#    ax1.set_ylim(-1,30) 
#    ax2.set_xlim(-0.75,0.75)    
#    ax2.set_ylim(0,0.5) 
#    
##    ax1.legend(loc='best',numpoints=1,prop=dict(size=10.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
##    ax1.set_xlim(0,1)       
##    ax1.set_ylim(0,np.max(muexEV_midS)*1.1)
##    plt.show()
##    plt.savefig('newF_BikSucks_ex-bulk_zdimensional-30.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
##    ax1.set_ylim(-1,50) 
##    plt.ylim(ymin = -1) 
##    ax1.set_xlim(-0.15,0.15) 
##    plt.savefig('newF_BikSucks_ex-bulk_zdimensional-50.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
###Biksucks
#    plt.close()


##########Contour and correlation length data
		##ABSOLUTELY SAVE SAVE SAVE!!!
    PM_Phi_Ex_file='PM_Phi_Ex.txt'
    PhiEx=[[float(x) for x in line.split()] for line in file(PM_Phi_Ex_file,"r").readlines()]
    xPhi,yEx = [],[]
    i=-1
    for x in PhiEx[0]:
    	i+=1
	if (i % 2)==0:
		xPhi.append(x)
	else:
  		yEx.append(x)		
    t1 = np.linspace(xPhi[0],xPhi[-1],20)[1:-1]
    interp_Ex_Phi = interpolate.LSQUnivariateSpline(xPhi,yEx,t1,k=5)
    fit_Ex_Phi = lambda Phi: interp_Ex_Phi(Phi)[0]

		
    cont_z,cont_phi,cont_dif = [],[],[]
    plotdata_neg,plotdata_0,plotdata_positive=[],[],[]
#    print bulknoise_fit
#    bulknoise_fit = np.max(bulknoise_fit)

    bulknoise_fit = 0.0835704567412
#    print '\nNEED TO FIX V and SigEff!!!\n'
    print 'ENFORCED max noise value (factors in all runs) = ',bulknoise_fit
    i=-1
    for (Evz,PhiWCAz,sigWCA,lamB,zpos,filename,L_z,ndSigma_s,ndSigeff,V,PhiB,Lbin,mark) in zip(muexEV_S,PhiWCAtot_S,sigWCA_S,BjerrumS,z_S,filenameS,L_zS,NDSigmaS,SigEffS,VEffS,PhiBulkS,L_binS,markS):
	if '_GC_' in filename:
		Evz=np.array(Evz[:-1])
	else:
		Evz=np.array(Evz)
    	zpos=np.array(zpos[:-1])
    	PhiWCAz = np.array(PhiWCAz)

	cor_z,cor_phi,cor_ev,cor_CSdif,cor_fitdif = [],[],[],[],[]
	CSdifnoise,fitdifnoise,evnoise = [],[],[]
	cor_V,cor_Sig =[],[]
	veff = -1E6
	sigeff = veff
	for (z,phi,ev,veff,sigeff) in zip(zpos,PhiWCAz,Evz,V,ndSigeff):	
		if z/L_z<=0.5:
			zfromwall = (z-1+0.299602925)/sigWCA
	
			difEV = fit_Ex_Phi(phi) - ev

			##These are unused...
#	    		cont_z.append(zfromwall)
#	    		cont_phi.append(phi)
#	    		cont_dif.append(dif)
	    		
			if '_GC_' in filename and z<1:
				phi = phi*Lbin/0.299602925
				difEV = fit_Ex_Phi(phi) - ev

			fitEV = fit_Ex_Phi(phi)
#			PhiCS = phi*0.954028946**3
			PhiCS = phi*0.903042806**3	
			CSEV = PhiCS*(8-9*PhiCS+3*PhiCS**2)/(1-PhiCS)**3
			difCS = CSEV - ev
#			print 'Phi,EVdif,CSdif = ',phi,difEV,difCS

			##Might be able to get rid of this
#	    		cor_z.append((z-1.+0.299602925)/sigWCA)
#	    		cor_phi.append(phi)
#	    		cor_ev.append(ev)
#	    		cor_CSdif.append(-ev+CSEV)
#	    		cor_fitdif.append(-ev+fitEV)
#	    		cor_V.append(veff)
#	    		cor_Sig.append(sigeff)
	    		
	    		temp = [zfromwall,phi,difEV,difCS,mark]  #This is for color map
	    		if difEV<=0.3 and difEV>-0.1: #Used to be abs(dif)<=0.3
	    			plotdata_0.append(temp)
    			elif difEV<-0.1:
	    			plotdata_neg.append(temp)    
    			else:
    				plotdata_positive.append(temp)
		else:
			break

#	cor_ev = cor_ev[::-1]
#	cor_z = cor_z[::-1]
#	cor_phi = cor_phi[::-1]
#	cor_CSdif = cor_CSdif[::-1]
#	cor_fitdif = cor_fitdif[::-1]
#	cor_V = cor_V[::-1]
#	cor_Sig = cor_Sig[::-1]



    print 'Trying to plot a color map -3 to 3'
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    colscheme = cm = plt.get_cmap('RdBu')
    cNorm = colors.Normalize(vmin = -3,vmax = 3)  
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=colscheme)  
    ifEV = 1
    ifEV = 0
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    from operator import itemgetter
    if ifEV:
        plotdata_neg.sort(key = itemgetter(2),reverse=True)
   	plotdata_0.sort(key = itemgetter(2))
   	plotdata_positive.sort(key = itemgetter(2))
    else:
        plotdata_neg.sort(key = itemgetter(3),reverse=True)
   	plotdata_0.sort(key = itemgetter(3))
   	plotdata_positive.sort(key = itemgetter(3))    	
    plotdata = plotdata_0  + plotdata_positive + plotdata_neg
    count = 0
    min_dif,max_dif = 1E9,-1E9
    overallMaxEV,overallMaxCS = -1E9,-1E9
    for xyz in plotdata:
    	x = xyz[0]
    	y = xyz[1]
    	zEV = xyz[2]
    	zCS = xyz[3]
    	mark = xyz[4]
    	if abs(zEV)<=3 and abs(zCS)<=3:
	    	overallMaxEV = np.max([zEV,overallMaxEV])
	    	overallMaxCS = np.max([zCS,overallMaxCS])
    		if ifEV:
    			z = zEV
		else:
	    		z = zCS
    	    	col = scalarMap.to_rgba(z)
    	    	min_dif,max_dif = np.min([z,min_dif]),np.max([z,max_dif])
	    	ax1.plot(x,y,marker=mark,ms=2.5,color=col,mec = col,mew = 0.)
    	else:
    		count+=1
    cbar_x = np.linspace(-0.40,4.9,3500)
    cbar_col = np.linspace(min_dif,max_dif,3500)
    print min_dif,max_dif
    print overallMaxEV,overallMaxCS
#    for (x,c) in zip(cbar_x,cbar_col):
#    	col = scalarMap.to_rgba(c)
#    	ax1.plot(x,0.540+0.01,marker='|',ms=2.5,color=col)
#    	ax1.plot(x,0.545+0.01,marker='|',ms=2.5,color=col)
#    	ax1.plot(x,0.55+0.01,marker='|',ms=2.5,color=col)
#    	ax1.plot(x,0.555+0.01,marker='|',ms=2.5,color=col)
#    	ax1.plot(x,0.56+0.01,marker='|',ms=2.5,color=col)
#    	ax1.plot(x,0.565+0.01,marker='|',ms=2.5,color=col)
#    	ax1.plot(x,0.57+0.01,marker='|',ms=2.5,color=col)		
    ax1.set_xlabel(r'$ (z-\delta_{\rm w}) / \sigma$',fontsize=10.)
    ax1.set_ylabel(r'$\Phi(\tilde z)$',fontsize=10.)
    ax1.set_xlim(-0.5,5)    
    ax1.xaxis.set_major_locator(MaxNLocator(5)) 
    ax1.yaxis.set_major_locator(MaxNLocator(7)) 
    ax1.set_ylim(0,0.60) 
    plt.setp(ax1.get_xticklabels(), fontsize=8.)
    plt.setp(ax1.get_yticklabels(), fontsize=8.)
    fig.set_size_inches(3.37*0.5,3.5)
    #    plt.show()   
    if ifEV:
        plt.savefig('EV_BikSucks_ColorMap-RdBu_-3to3.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
    else:
    	plt.savefig('CS_BikSucks_ColorMap-RdBu_-3to3.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
    plt.close()
    print 'Number of points not added = ',count






    print 'Plotting CorCap vs l_cor  expressions'    
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    #For square inset
    i,cc=-1,-1  #Rando iteraters
    labels=[]
    maxX = -1E9
    maxY = maxX
    xy_data=[[],[]]   
###    P = '0.0438'  #AvERAGE
####    P = '0.0436891008773'
###    Tname = 'CS_0.0438_7.txt'
###    P=float(P)
###    QV_CS=np.array([[float(x) for x in line.split()] for line in file(Tname,"r").readlines()])
###    ax1.errorbar(QV_CS[:,4],np.array(QV_CS[:,5])/np.array(QV_CS[:,4]),yerr=None,ls='-',lw=0.75/(1.-P),color='blue',marker='None',label=r'${\rm CS}$')#

###    Tname = 'Bik_0.0438_5.txt'
###    P=float(P)
###    QV_CS=np.array([[float(x) for x in line.split()] for line in file(Tname,"r").readlines()])
###    ax1.errorbar(QV_CS[:,4],np.array(QV_CS[:,5])/np.array(QV_CS[:,4]),yerr=None,ls='--',lw=0.75/(1.-P),color='k',marker='None',label=r'${\rm Bik}$')#
###    
###    theory = np.linspace(0,5)
###    ax1.errorbar(theory,2*np.sinh(theory/2)/theory,yerr=None,ls='-',lw=0.3,color='k',marker='None',label=r'${\rm GC}$')#
    for (Sigmaz,Vz,Sigma,filename,zpos) in zip(SigEffS,VEffS,NDSigmaS,filenameS,zEffS):
#    	print zpos[:5]
    	
	ls_x,lw_x='',0.
	ms_x = 3.5
	label_str=''
	legS=[]
	alabel=''
	skip = 0
	col = 'k'
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
		label_str=r'$\tilde \Sigma = 7.0$'
		alabel=r'$0.92$'	
		## NEW
		col = ROYGBIV_map(0.,10)
		x,y  = 0.477733658,1.079078143
	elif filename == 'Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='o'
		label_str=r'$5.4$'
		alabel=r'$0.8$'
		## NEW P2 Markers
		special_ms='*'	
		col = ROYGBIV_map(1.,10)	
		x,y = 0.366622547,1.087171253	
	elif filename == 'Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='o'	
		label_str=r'$4.2$'
		alabel=r'$0.6$'
		col = ROYGBIV_map(2.,10)
		x,y = 0.255511436,1.087410851	
	elif filename == 'Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='s'	
		label_str=r'$3.6$'
		alabel=r'$0.5$'
		## NEW P2 Markers
		special_ms='d'	
		col = ROYGBIV_map(3.,10)	
		x,y = 0.255511436,1.066611769
	elif filename == 'Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='s'
		ms_x=3.5
		label_str=r'$3.1$'
		alabel=r'$0.4$'
		col = ROYGBIV_map(4.,10)	
		x,y = 0.255511436,1.043452831
	elif filename == 'Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='p'	
		ms_x=3.5	
		label_str=r'$2.3$'
		alabel=r'$0.3$'
		col = ROYGBIV_map(5.,10)
		x,y = 0.144400325,0.975297948	
	elif filename == 'Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='p'	
		ms_x=3.5
		label_str=r'$2.0$'
		alabel=r'$0.5$'
		## NEW P2 Markers
		special_ms='v'	
		col = ROYGBIV_map(6.,10)
		x,y = 0.144400325,0.932793621	
	elif filename == 'Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='H'
		ms_x=3.5
		label_str=r'$1.6$'
		alabel=r'$0.22$'
		col = ROYGBIV_map(7.,10)	
		x,y = 0.144400325,0.90746003
	elif filename == 'Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='H'	
		label_str=r'$1.0$'
		alabel=r'$0.22$'
		## NEW P2 Markers
		special_ms='>'	
		col = ROYGBIV_map(8.,10)
		x,y = 0.033289214,0.644541654	
	elif filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt':
		edge_width=0.1
		special_ms='D'	
		label_str=r'$0.5$'
		alabel=r'$0.13$'
		col = ROYGBIV_map(9.,10)
		x,y = 0.033289214,0.644541654	
	elif filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt':
		edge_width=0.4
		special_ms='<'
		label_str=r'$0.0$'
		alabel=r'$0.06$'
	else:
		skip = 1
		
	if not skip:
#	    	col = ROYGBIV_map(abs(Sigma),10)
		ax1.errorbar(x,y,yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)

    ax1.set_ylabel(r'$ -\~ \Sigma  / \~ \phi$',fontsize=10.)
    ax1.set_xlabel(r'$ (\phi - \phi^{\rm B})  qe / k_{\rm B}T $',fontsize=10.)

    plt.setp(ax1.get_xticklabels(), fontsize=8.)
    plt.setp(ax1.get_yticklabels(), fontsize=8.)
#    ax1.set_xlim(0,maxX*1.1) 
#    ax1.set_ylim(0,maxY*1.1)
#    ax1.set_xlim(0,7) 
#    ax1.set_ylim(0,7.5)
#    ax1.xaxis.set_major_locator(MaxNLocator(6)) 
#Biksucks
    fig.set_size_inches(3.37,3.5)
    fig.subplots_adjust(right=0.98) ##Lower puts more whitespace on the right
    fig.subplots_adjust(bottom=0.10)##Larger adds whitespace
    fig.subplots_adjust(left=0.13) #Larger puts more whitespace on the left
    fig.subplots_adjust(top=0.98) ##Smaller adds whitespace to tops
    ax1.legend(loc='best',numpoints=1,prop={"size":9},columnspacing=0.07,borderpad=0.18,labelspacing=0.10,handletextpad=0.07)

#    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.05,borderpad=0.25,labelspacing=0.05,ncol=1,handletextpad=0.05)    
#    plt.savefig('F??_BikSucks_CorCap_vs_lcor.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight') 
#    plt.show()
    plt.close() 





    print 'Plotting WCA minus HS and CS...'
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    MAX = -1e9
    exHS,exCS=[],[]
    test2=[]
    for (z_pos,WCA_B,PhiWCA_B,filename) in zip(z_S,muexEV_bulkS,PhiWCAs,filenameS):
#      	PhiCS_B = PhiWCA_B
#    	CS_B = PhiCS_B*(8-9*PhiCS_B+3*PhiCS_B**2)/(1-PhiCS_B)**3
#      	test2.append(CS_B)
      	  	
    	PhiCS_B = PhiWCA_B*0.954028946**3
    	CS_B = PhiCS_B*(8-9*PhiCS_B+3*PhiCS_B**2)/(1-PhiCS_B)**3    	
    	exCS.append(CS_B)

    z_HS=np.array([[float(x) for x in line.split()] for line in file('MeasuredHS_sigWCA.txt',"r").readlines()])
#    PhiS_sigWCA =z_HS[:,0]
    HS_sigWCA = z_HS[:,1]
    z_95HS=np.array([[float(x) for x in line.split()] for line in file('MeasuredHS_0.95sigWCA.txt',"r").readlines()])
#    PhiS_sigWCA =z_95HS[:,0]
    HS_95sigWCA = z_95HS[:,1]    


#    ax1.errorbar(PhiWCAs,np.array(muexEV_bulkS) - np.array(HS_sigWCA),yerr=None,color='red',ls='-',label=' '+r'$\mu^{\rm WCA} - \mu^{\rm HS}$')	
    ax1.errorbar(PhiWCAs,np.array(muexEV_bulkS) - np.array(HS_95sigWCA),yerr=None,color='k',ls='-',label=' '+r'$\mu^{\rm WCA} - \mu^{\rm HS}(0.95 \cdot \sigma)$')	
    ax1.errorbar(PhiWCAs,np.array(muexEV_bulkS) - np.array(exCS),yerr=None,color='blue',ls='-',label=' '+r'$\mu^{\rm WCA} - \mu^{\rm CS}(0.95 \cdot \sigma) $')	
    ax1.errorbar(PhiWCAs,np.array(HS_95sigWCA) - np.array(exCS),yerr=None,color='green',ls='-',label=' '+r'$\mu^{\rm HS} - \mu^{\rm CS}$')	

    ax1.set_xlabel(r'$\Phi^{\rm B}$',fontsize=10.)
    ax1.set_ylabel(r'$\Delta \mu^{\rm EV} / k_{\rm B} T$',fontsize=10.)
    plt.setp(ax1.get_xticklabels(), fontsize=8.)
    plt.setp(ax1.get_yticklabels(), fontsize=8.)
    ax1.legend(loc='best',numpoints=1,prop=dict(size=10.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
#    ax1.set_xlim(0,1)       
#    ax1.set_ylim(0,np.max(muexEV_midS)*1.1)
    fig.set_size_inches(3.37,2.5)
    plt.savefig('F4_BikSucks_HS_WCA_CS.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
##Biksucks
    plt.close()


	###Move this below

####    ##This plots rhoF(z) for simulation and modified theories
#    nonDim='yes'
##    Bikerman='yes'
##    Bikerman='no'
#    CarnahanStarling='yes'
#    CarnahanStarling='no'
##    fromlcor = 'yes'
##    fromlcor = 'no'
#    if nonDim=='yes':
#      print 'Plotting NonDim rhoF(z)...'
#    else:
#      print 'Plotting Dim rhoF(z)...'	
#    if Bikerman=='yes':
#		print '\t\t...with Bikerman theory...'
#		#Note to user: Bikerman theory must already be evaluated using MATLAB code, which is in ~/sims/PBik_Solver.m
#    if CarnahanStarling=='yes':
#	      print '\t...with Carnahan-Starling theory...'
#    i=-1
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    #BikSucks
#    cc=-1
#    ax1.errorbar([-1,-2],[-1,-2],yerr=None,color='k',ls='-',lw=1.0,label='  '+ r'$\tilde \mu^{\rm CS}(0.95 \cdot \sigma)$') #
##    ax1.errorbar([-1,-2],[-1,-2],yerr=None,color='k',ls='-.',lw=1.0,label='  '+ r'$\tilde \mu^{\rm Bik}$')
##    evcorS,lcorS=[],[]
#    for (z_positions,characteristic_length,muex_EV,L_z,sigWCA,PhiWCA,lamD,lamB,filename,n0,free) in zip(z_originalS,characteristic_lengthS,muexEV_S,L_zS,sigWCA_S,PhiWCAtot_S,LDs,BjerrumS,filenameS,n0s,NfS):
#    	i+=1
#	cc+=1

#	col = ROYGBIV_map((sigWCA)/lamD,3,1) #*0.954028946

#	rat = (sigWCA/lamB)
#	if rat <= 0.05:
#		mark = 'o'
#	elif rat>0.05 and rat<=0.2:
#		mark = 's'
#	elif rat>0.2 and rat<=1.:
#		mark = 'h'
#	elif rat>1. and rat<=10.:
#		mark = 'p'
#	elif rat>10. and rat<=120.:
#		mark = 'D'
##	elif rat>1.5:
##		mark = 'D'
#	msize = 5.0	

#	ls_x,lw_x='',0.
#	ms_x = 4.0
#	label_str=''
#	alabel=''
#	special_ms='x'
#	edge_width=0.1
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
#		label_str=r'$\tilde \Sigma = 7.0$'
#		evcorS.append(3.46615)
#		col = ROYGBIV_map(0.,10)	
#	elif filename == 'Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='o'
#		label_str=r'$5.4$'
#		special_ms='*'	
#		col = ROYGBIV_map(1.,10)		
#	elif filename == 'Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1162_0.014195_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='o'	
#		label_str=r'$4.2$'
#		col = ROYGBIV_map(2.,10)	
#	elif filename == 'Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='s'	
#		label_str=r'$3.6$'
#		special_ms='d'		
#		col = ROYGBIV_map(3.,10)	
#	elif filename == 'Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1104_0.0106795_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='s'
#		ms_x=3.5
#		label_str=r'$3.1$'
#		col = ROYGBIV_map(4.,10)	
#	elif filename == 'Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1066_0.0078345_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='p'	
#		ms_x=3.5	
#		label_str=r'$2.3$'
#		col = ROYGBIV_map(5.,10)	
#	elif filename == 'Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='p'	
#		ms_x=3.5
#		label_str=r'$2.0$'
#		special_ms='v'	
#		col = ROYGBIV_map(6.,10)	
#	elif filename == 'Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='H'
#		ms_x=3.5
#		label_str=r'$1.6$'
#		col = ROYGBIV_map(7.,10)	
#	elif filename == 'Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='H'	
#		label_str=r'$1.0$'
#		special_ms='>'	
#		col = ROYGBIV_map(8.,10)	
#	elif filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='D'	
#		label_str=r'$0.5$'
#		col = ROYGBIV_map(9.,10)	
#	elif filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='<'
#		label_str=r'$0.0$'
#		col = 'k'
#	ms_x = 4.0

#	z_positions=[z/L_z for z in z_positions] #NonDim distances

#	ax1.errorbar(np.array(z_positions[:len(free)/2])*(L_z/sigWCA)-(1+0.299602925)/sigWCA,free[:len(free)/2],yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)

##    if sys.argv[1]=='F6_ex.txt':
###	ax1.errorbar(lcorS,evcorS,yerr=None,color='magenta')
##    	evcorS = np.array(evcorS) + np.array(muexEV_bulkS)
###    	print 'post correction = ',evcorS
##    	ax1.errorbar(lcorS,evcorS,yerr=None,color='k')
#    
#    cc=-1
#    i=-1
#    for (z_positions,characteristic_length,muex_EV,L_z,sigWCA,PhiWCA,lamD,lamB,filename,n0,col) in zip(z_originalS,characteristic_lengthS,muexEV_S,L_zS,sigWCA_S,PhiWCAtot_S,LDs,BjerrumS,filenameS,n0s,theorycolS):
#    	i+=1
#	cc+=1

#	ls_x,lw_x='',0.
#	ms_x = 3.0
#	label_str=''
#	alabel=''
#	special_ms='x'
#	edge_width=0.1
#	if filename == 'Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$11.6$'		
#	elif filename == 'Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$11.1$'
#	elif filename == 'Analyzed_GC_1312_0.0241789536991_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$\tilde \Sigma = 6.0$'
#		ze = ['6.72','4.66']
#	elif filename == 'Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$5.0$'
#		ze = ['5.16','3.85']
#	elif filename == 'Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1162_0.014195_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$3.9$'
#		ze = ['4.06','3.28']
#	elif filename == 'Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$3.4$'
#		ze = ['3.54','2.86']
#	elif filename == 'Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1104_0.0106795_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$3.0$'
#		ze = ['3.08','2.48']
#	elif filename == 'Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1066_0.0078345_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$2.2$'
#		ze = ['2.32','2.01']
#	elif filename == 'Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt':	
#		label_str=r'$1.9$'
#		ze = ['2.01','1.76']
#	elif filename == 'Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$1.6$'
#		ze = ['1.65','1.44']
#	elif filename == 'Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$1.0$'
#		ze = ['1.06','0.99']
#	elif filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$0.5$'
#		ze = ['0.50','0.47']
#	elif filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$0.0$'
#		ze = ['0.00','0.00']

#	if Bikerman=='yes' and ze!='0.00':
#		##See modifications to CS code below, if Bikerman is ever needed in these plots - 09/27/13 08:24:24 
#	  Bik_file='Bik_0.0438_' + ze + '.txt'
#          MMA_Bik=[[float(x) for x in line.split()] for line in file(Bik_file,"r").readlines()]
#          MMAz,PhiBik=[],[]
#	  for line in MMA_Bik:
#	  	MMAz.append(line[0])
#	  	Phi = (line[1]+line[2])*n0*(np.pi/6)*sigWCA**3
#	  	PhiBik.append(Phi) 
#	  PhiBik = np.array(PhiBik)
#	  Bik_ex = -np.log(1 - PhiBik)
#	  Bik_ex = Bik_ex - Bik_ex[len(Bik_ex)-1]   
#	  ax1.plot(np.array(MMAz)*lamD/sigWCA,Bik_ex,color='k',lw=1.5,ls='-')    
#	  ax1.plot(np.array(MMAz)*lamD/sigWCA,Bik_ex,color=col,lw=1.0,ls='-.')    
#	if CarnahanStarling=='yes' and ze[0]!='0.00':
#		  if fromlcor == 'yes': #This is plotting theory from lcor onwards
#		    	ze = ze[1]
#		    	xshift = lcorS[i]
#		  else:
#			ze = ze[0]	
#			xshift = 0		
#		  CS_file='CS_0.0438_' + ze + '.txt'
#		  MMA_CS=[[float(x) for x in line.split()] for line in file(CS_file,"r").readlines()]
#		  MMAz,freeCS = [],[]
#		  for line in MMA_CS:
#		  	MMAz.append(line[0])
#		  	freeCS.append(-line[3])
#		  xplot = np.array(MMAz)*lamD/sigWCA-(1+0.299602925)/sigWCA + xshift
#		  ax1.plot(xplot,freeCS,color=col,lw=1.0,ls='-') ##Old
#        elif CarnahanStarling=='yes' and ze[0]=='0.00':
#        	xplot = np.array(z_positions)*lamD/sigWCA-(1+0.299602925)/sigWCA
#        	ax1.plot(xplot,np.zeros(len(xplot)),color='k',lw=1.0,ls='-')
#    ax1.set_xlabel(r'$z / \sigma$',fontsize=10.)
#    ax1.set_ylabel(r'$-\tilde \rho_{\rm f}$',fontsize=10.)
#    plt.setp(ax1.get_xticklabels(), fontsize=8.)
#    plt.setp(ax1.get_yticklabels(), fontsize=8.)

#    ax1.set_xlim(-1/sigWCA*1.5,10)            
##    ax1.set_ylim(-0.25,10)

##    ax1.legend(loc='best',numpoints=1,prop={"size":10},columnspacing=0.07,borderpad=0.18,labelspacing=0.10,handletextpad=0.07) 
#    fig.set_size_inches(3.37,3.5)
##    if fromlcor == 'yes': #This is plotting theory from lcor onwards
##	    plt.savefig('F4_rhoF_vs_z_BikSucks_lcor.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
##    else:
##	    plt.savefig('F4_rhoF_vs_z_BikSucks_zmin.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
##    plt.show()
#    plt.close()
##    #Done plotting rhoFree(z)


	###Move this below


####    ##This plots Voltage(z) for simulation and modified theories
##    Bikerman='yes'
##    Bikerman='no'
#    CarnahanStarling='yes'
#    CarnahanStarling='no'
##    fromlcor = 'yes'
#    fromlcor = 'no'
#    print 'Plotting Voltage(z)...'
#    if Bikerman=='yes':
#		print '\t\t...with Bikerman theory...'
#		#Note to user: Bikerman theory must already be evaluated using MATLAB code, which is in ~/sims/PBik_Solver.m
#    if CarnahanStarling=='yes':
#	      print '\t...with Carnahan-Starling theory...'
#    i=-1
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    #BikSucks
#    cc=-1
#    ax1.errorbar([-1,-2],[-1,-2],yerr=None,color='k',ls='-',lw=1.0,label='  '+ r'$\tilde \mu^{\rm CS}(0.95 \cdot \sigma)$') #
##    ax1.errorbar([-1,-2],[-1,-2],yerr=None,color='k',ls='-.',lw=1.0,label='  '+ r'$\tilde \mu^{\rm Bik}$')
##    evcorS,lcorS=[],[]
#    for (z_positions,characteristic_length,muex_EV,L_z,sigWCA,PhiWCA,lamD,lamB,filename,n0,volt) in zip(z_originalS,characteristic_lengthS,muexEV_S,L_zS,sigWCA_S,PhiWCAtot_S,LDs,BjerrumS,filenameS,n0s,VEffS):
#    	i+=1
#	cc+=1

#	col = ROYGBIV_map((sigWCA)/lamD,3,1) #*0.954028946

#	rat = (sigWCA/lamB)
#	if rat <= 0.05:
#		mark = 'o'
#	elif rat>0.05 and rat<=0.2:
#		mark = 's'
#	elif rat>0.2 and rat<=1.:
#		mark = 'h'
#	elif rat>1. and rat<=10.:
#		mark = 'p'
#	elif rat>10. and rat<=120.:
#		mark = 'D'
##	elif rat>1.5:
##		mark = 'D'
#	msize = 5.0	

#	ls_x,lw_x='',0.
#	ms_x = 4.0
#	label_str=''
#	alabel=''
#	special_ms='x'
#	edge_width=0.1
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
#		label_str=r'$\tilde \Sigma = 7.0$'
#		evcorS.append(3.46615)
#		col = ROYGBIV_map(0.,10)	
#	elif filename == 'Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='o'
#		label_str=r'$5.4$'
#		special_ms='*'	
#		col = ROYGBIV_map(1.,10)		
#	elif filename == 'Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1162_0.014195_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='o'	
#		label_str=r'$4.2$'
#		col = ROYGBIV_map(2.,10)	
#	elif filename == 'Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='s'	
#		label_str=r'$3.6$'
#		special_ms='d'		
#		col = ROYGBIV_map(3.,10)	
#	elif filename == 'Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1104_0.0106795_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='s'
#		ms_x=3.5
#		label_str=r'$3.1$'
#		col = ROYGBIV_map(4.,10)	
#	elif filename == 'Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1066_0.0078345_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='p'	
#		ms_x=3.5	
#		label_str=r'$2.3$'
#		col = ROYGBIV_map(5.,10)	
#	elif filename == 'Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='p'	
#		ms_x=3.5
#		label_str=r'$2.0$'
#		special_ms='v'	
#		col = ROYGBIV_map(6.,10)	
#	elif filename == 'Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='H'
#		ms_x=3.5
#		label_str=r'$1.6$'
#		col = ROYGBIV_map(7.,10)	
#	elif filename == 'Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='H'	
#		label_str=r'$1.0$'
#		special_ms='>'	
#		col = ROYGBIV_map(8.,10)	
#	elif filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='D'	
#		label_str=r'$0.5$'
#		col = ROYGBIV_map(9.,10)	
#	elif filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='<'
#		label_str=r'$0.0$'
#		col = 'k'
#	ms_x = 4.0

#	z_positions=[z/L_z for z in z_positions] #NonDim distances
#	ax1.errorbar(np.array(z_positions[:len(volt)/2])*(L_z/sigWCA)-(1+0.299602925)/sigWCA,volt[:len(volt)/2],yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)

##    if sys.argv[1]=='F6_ex.txt':
##	ax1.errorbar(lcorS,evcorS,yerr=None,color='magenta')
##    	evcorS = np.array(evcorS) + np.array(muexEV_bulkS)
###    	print 'post correction = ',evcorS
##    	ax1.errorbar(lcorS,evcorS,yerr=None,color='k')
#    cc=-1
#    i=-1
#    for (z_positions,characteristic_length,muex_EV,L_z,sigWCA,PhiWCA,lamD,lamB,filename,n0,col) in zip(z_originalS,characteristic_lengthS,muexEV_S,L_zS,sigWCA_S,PhiWCAtot_S,LDs,BjerrumS,filenameS,n0s,theorycolS):
#    	i+=1
#	cc+=1

#	ls_x,lw_x='',0.
#	ms_x = 3.0
#	label_str=''
#	alabel=''
#	special_ms='x'
#	edge_width=0.1
#	if filename == 'Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$11.6$'		
#	elif filename == 'Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$11.1$'
#	elif filename == 'Analyzed_GC_1312_0.0241789536991_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$\tilde \Sigma = 6.0$'
#		ze = ['6.72','4.66']
#	elif filename == 'Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$5.0$'
#		ze = ['5.16','3.85']
#	elif filename == 'Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1162_0.014195_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$3.9$'
#		ze = ['4.06','3.28']
#	elif filename == 'Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$3.4$'
#		ze = ['3.54','2.86']
#	elif filename == 'Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1104_0.0106795_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$3.0$'
#		ze = ['3.08','2.48']
#	elif filename == 'Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1066_0.0078345_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$2.2$'
#		ze = ['2.32','2.01']
#	elif filename == 'Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt':	
#		label_str=r'$1.9$'
#		ze = ['2.01','1.76']
#	elif filename == 'Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$1.6$'
#		ze = ['1.65','1.44']
#	elif filename == 'Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$1.0$'
#		ze = ['1.06','0.99']
#	elif filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$0.5$'
#		ze = ['0.50','0.47']
#	elif filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$0.0$'
#		ze = ['0.00','0.00']

#	if Bikerman=='yes' and ze!='0.00':
#		##See modifications to CS code below, if Bikerman is ever needed in these plots - 09/27/13 08:24:24 
#	  Bik_file='Bik_0.0438_' + ze + '.txt'
#          MMA_Bik=[[float(x) for x in line.split()] for line in file(Bik_file,"r").readlines()]
#          MMAz,PhiBik=[],[]
#	  for line in MMA_Bik:
#	  	MMAz.append(line[0])
#	  	Phi = (line[1]+line[2])*n0*(np.pi/6)*sigWCA**3
#	  	PhiBik.append(Phi) 
#	  PhiBik = np.array(PhiBik)
#	  Bik_ex = -np.log(1 - PhiBik)
#	  Bik_ex = Bik_ex - Bik_ex[len(Bik_ex)-1]   
#	  ax1.plot(np.array(MMAz)*lamD/sigWCA,Bik_ex,color='k',lw=1.5,ls='-')    
#	  ax1.plot(np.array(MMAz)*lamD/sigWCA,Bik_ex,color=col,lw=1.0,ls='-.')    
#	if CarnahanStarling=='yes' and ze[0]!='0.00':
#		  if fromlcor == 'yes': #This is plotting theory from lcor onwards
#		    	ze = ze[1]
#		    	xshift = lcorS[i]
#		  else:
#			ze = ze[0]	
#			xshift = 0		
#		  CS_file='CS_0.0438_' + ze + '.txt'
#		  MMA_CS=[[float(x) for x in line.split()] for line in file(CS_file,"r").readlines()]
#		  MMAz,voltCS = [],[]
#		  for line in MMA_CS:
#		  	MMAz.append(line[0])
#		  	voltCS.append(line[4])
#		  xplot = np.array(MMAz)*lamD/sigWCA-(1+0.299602925)/sigWCA + xshift
#		  ax1.plot(xplot,voltCS,color=col,lw=1.0,ls='-') ##Old
#        elif CarnahanStarling=='yes' and ze[0]=='0.00':
#        	ax1.plot(xplot,np.zeros(len(xplot)),color='k',lw=1.0,ls='-')
#    ax1.set_xlabel(r'$z / \sigma$',fontsize=10.)
#    ax1.set_ylabel(r'$ \tilde \phi$',fontsize=10.)
#    plt.setp(ax1.get_xticklabels(), fontsize=8.)
#    plt.setp(ax1.get_yticklabels(), fontsize=8.)
##    plt.savefig('F?_Volt_vs_z_BikSucks_nolims.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
#    ax1.set_xlim(-1,15)            
#    ax1.set_ylim(-0.25,7)
##    ax1.legend(loc='best',numpoints=1,prop={"size":10},columnspacing=0.07,borderpad=0.18,labelspacing=0.10,handletextpad=0.07) 
#    fig.set_size_inches(3.37,3.5)
##    if fromlcor == 'yes': #This is plotting theory from lcor onwards
##	    plt.savefig('F?_Volt_vs_z_BikSucks_lcor.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
##    else:
##	    plt.savefig('F?_Volt_vs_z_BikSucks_zmin.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
##    plt.show()
#    plt.close()



####    ##This plots RhoCounterCo(z) for simulation and modified theories
##    Bikerman='yes'
##    Bikerman='no'
#    CarnahanStarling='yes'
#    CarnahanStarling='no'
##    fromlcor = 'yes'
#    fromlcor = 'no'
#    print 'Plotting RhoCounterCo(z)...'
#    if Bikerman=='yes':
#		print '\t\t...with Bikerman theory...'
#		#Note to user: Bikerman theory must already be evaluated using MATLAB code, which is in ~/sims/PBik_Solver.m
#    if CarnahanStarling=='yes':
#	      print '\t...with Carnahan-Starling theory...'
#    i=-1
#    fig=plt.figure()
#    fig.set_size_inches(3.37,3.5)
#    ax1=fig.add_subplot(211)
#    ax2=fig.add_subplot(212)
#    #BikSucks
#    cc=-1
#    ax1.errorbar([-1,-2],[-1,-2],yerr=None,color='k',ls='-',lw=1.0,label='  '+ r'$\tilde \mu^{\rm CS}(0.95 \cdot \sigma)$') #
##    ax1.errorbar([-1,-2],[-1,-2],yerr=None,color='k',ls='-.',lw=1.0,label='  '+ r'$\tilde \mu^{\rm Bik}$')
##    evcorS,lcorS=[],[]
#    for (z_positions,characteristic_length,muex_EV,L_z,sigWCA,PhiWCA,lamD,lamB,filename,n0,rhoco,rhoc) in zip(z_originalS,characteristic_lengthS,muexEV_S,L_zS,sigWCA_S,PhiWCAtot_S,LDs,BjerrumS,filenameS,n0s,rhoplusS,rhominusS):
#    	i+=1
#	cc+=1

#	col = ROYGBIV_map((sigWCA)/lamD,3,1) #*0.954028946

#	rat = (sigWCA/lamB)
#	if rat <= 0.05:
#		mark = 'o'
#	elif rat>0.05 and rat<=0.2:
#		mark = 's'
#	elif rat>0.2 and rat<=1.:
#		mark = 'h'
#	elif rat>1. and rat<=10.:
#		mark = 'p'
#	elif rat>10. and rat<=120.:
#		mark = 'D'
##	elif rat>1.5:
##		mark = 'D'
#	msize = 5.0	

#	ls_x,lw_x='',0.
#	ms_x = 4.0
#	label_str=''
#	alabel=''
#	special_ms='x'
#	edge_width=0.1
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
#		label_str=r'$\tilde \Sigma = 7.0$'
#		evcorS.append(3.46615)
#		col = ROYGBIV_map(0.,10)	
#	elif filename == 'Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='o'
#		label_str=r'$5.4$'
#		special_ms='*'	
#		col = ROYGBIV_map(1.,10)		
#	elif filename == 'Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1162_0.014195_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='o'	
#		label_str=r'$4.2$'
#		col = ROYGBIV_map(2.,10)	
#	elif filename == 'Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='s'	
#		label_str=r'$3.6$'
#		special_ms='d'		
#		col = ROYGBIV_map(3.,10)	
#	elif filename == 'Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1104_0.0106795_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='s'
#		ms_x=3.5
#		label_str=r'$3.1$'
#		col = ROYGBIV_map(4.,10)	
#	elif filename == 'Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1066_0.0078345_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='p'	
#		ms_x=3.5	
#		label_str=r'$2.3$'
#		col = ROYGBIV_map(5.,10)	
#	elif filename == 'Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='p'	
#		ms_x=3.5
#		label_str=r'$2.0$'
#		special_ms='v'	
#		col = ROYGBIV_map(6.,10)	
#	elif filename == 'Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='H'
#		ms_x=3.5
#		label_str=r'$1.6$'
#		col = ROYGBIV_map(7.,10)	
#	elif filename == 'Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='H'	
#		label_str=r'$1.0$'
#		special_ms='>'	
#		col = ROYGBIV_map(8.,10)	
#	elif filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.1
#		special_ms='D'	
#		label_str=r'$0.5$'
#		col = ROYGBIV_map(9.,10)	
#	elif filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt':
#		edge_width=0.4
#		special_ms='<'
#		label_str=r'$0.0$'
#		col = 'k'
#	ms_x = 4.0

#	z_positions=[z/L_z for z in z_positions] #NonDim distances

#	ax1.errorbar(np.array(z_positions[:len(rhoco)])*(L_z/sigWCA)-(1+0.299602925)/sigWCA,rhoco,yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x,label=label_str)
#	ax2.errorbar(np.array(z_positions[:len(rhoc)])*(L_z/sigWCA)-(1+0.299602925)/sigWCA,rhoc,yerr=None,color=col,marker=special_ms,mew=edge_width,ms=ms_x,ls=ls_x,lw=lw_x)

##    if sys.argv[1]=='F6_ex.txt':
##	ax1.errorbar(lcorS,evcorS,yerr=None,color='magenta')
##    	evcorS = np.array(evcorS) + np.array(muexEV_bulkS)
###    	print 'post correction = ',evcorS
##    	ax1.errorbar(lcorS,evcorS,yerr=None,color='k')
#    
#    cc=-1
#    i=-1
#    for (z_positions,characteristic_length,muex_EV,L_z,sigWCA,PhiWCA,lamD,lamB,filename,n0,col) in zip(z_originalS,characteristic_lengthS,muexEV_S,L_zS,sigWCA_S,PhiWCAtot_S,LDs,BjerrumS,filenameS,n0s,theorycolS):
#    	i+=1
#	cc+=1

#	ls_x,lw_x='',0.
#	ms_x = 3.0
#	label_str=''
#	alabel=''
#	special_ms='x'
#	edge_width=0.1
#	if filename == 'Analyzed_GC_1352_0.0403345_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$11.6$'		
#	elif filename == 'Analyzed_GC_1412_0.0403345_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$11.1$'
#	elif filename == 'Analyzed_GC_1312_0.0241789536991_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$\tilde \Sigma = 6.0$'
#		ze = ['6.72','4.66']
#	elif filename == 'Analyzed_GC_1232_0.018603_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$5.0$'
#		ze = ['5.16','3.85']
#	elif filename == 'Analyzed_GC_1162_0.014194972246_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1162_0.014195_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$3.9$'
#		ze = ['4.06','3.28']
#	elif filename == 'Analyzed_GC_1136_0.0123408042618_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$3.4$'
#		ze = ['3.54','2.86']
#	elif filename == 'Analyzed_GC_1104_0.0106795300666_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1104_0.0106795_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$3.0$'
#		ze = ['3.08','2.48']
#	elif filename == 'Analyzed_GC_1066_0.00783450375952_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1066_0.0078345_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$2.2$'
#		ze = ['2.32','2.01']
#	elif filename == 'Analyzed_GC_1048_0.0066065_0.1_3.0_500.0_23.73.txt':	
#		label_str=r'$1.9$'
#		ze = ['2.01','1.76']
#	elif filename == 'Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt' or filename=='Analyzed_GC_1034_0.005482_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$1.6$'
#		ze = ['1.65','1.44']
#	elif filename == 'Analyzed_GC_1016_0.003474_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$1.0$'
#		ze = ['1.06','0.99']
#	elif filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt' or filename == 'Analyzed_GC_1004_0.001684_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$0.5$'
#		ze = ['0.50','0.47']
#	elif filename == 'Analyzed_GC_1000_0.0_0.1_3.0_500.0_23.73.txt':
#		label_str=r'$0.0$'
#		ze = ['0.00','0.00']

#	if Bikerman=='yes' and ze!='0.00':
#		##See modifications to CS code below, if Bikerman is ever needed in these plots - 09/27/13 08:24:24 
#	  Bik_file='Bik_0.0438_' + ze + '.txt'
#          MMA_Bik=[[float(x) for x in line.split()] for line in file(Bik_file,"r").readlines()]
#          MMAz,PhiBik=[],[]
#	  for line in MMA_Bik:
#	  	MMAz.append(line[0])
#	  	Phi = (line[1]+line[2])*n0*(np.pi/6)*sigWCA**3
#	  	PhiBik.append(Phi) 
#	  PhiBik = np.array(PhiBik)
#	  Bik_ex = -np.log(1 - PhiBik)
#	  Bik_ex = Bik_ex - Bik_ex[len(Bik_ex)-1]   
#	  ax1.plot(np.array(MMAz)*lamD/sigWCA,Bik_ex,color='k',lw=1.5,ls='-')    
#	  ax1.plot(np.array(MMAz)*lamD/sigWCA,Bik_ex,color=col,lw=1.0,ls='-.')    
#	if CarnahanStarling=='yes' and ze[0]!='0.00':
#		  if fromlcor == 'yes': #This is plotting theory from lcor onwards
#		    	ze = ze[1]
#		    	xshift = lcorS[i]
#		  else:
#			ze = ze[0]	
#			xshift = 0		
#		  CS_file='CS_0.0438_' + ze + '.txt'
#		  MMA_CS=[[float(x) for x in line.split()] for line in file(CS_file,"r").readlines()]
#		  MMAz,rhocCS,rhocoCS = [],[],[]
#		  for line in MMA_CS:
#		  	MMAz.append(line[0])
#		  	rhocCS.append(line[1])
#		  	rhocoCS.append(line[2])
#		  xplot = np.array(MMAz)*lamD/sigWCA-(1+0.299602925)/sigWCA + xshift
#		  ax2.plot(xplot,rhocoCS,color=col,lw=1.0,ls='-') ##Old
#		  ax1.plot(xplot,rhocCS,color=col,lw=1.0,ls='-') ##Old
#		  
#        elif CarnahanStarling=='yes' and ze[0]=='0.00':
#        	ax1.plot(xplot,np.ones(len(xplot)),color='k',lw=1.0,ls='-')
#        	ax2.plot(xplot,np.ones(len(xplot)),color='k',lw=1.0,ls='-')	

#    ax1.set_ylabel(r'$\tilde \rho_{\rm co}$',fontsize=10.)
#    ax2.set_ylabel(r'$\tilde \rho_{\rm c}$',fontsize=10.)
#    ax2.set_xlabel(r'$z / \sigma$',fontsize=10.)
#    
#    plt.setp(ax1.get_xticklabels(), fontsize=8.)
#    plt.setp(ax1.get_yticklabels(), fontsize=8.)
#    plt.setp(ax2.get_xticklabels(), fontsize=8.)
#    plt.setp(ax2.get_yticklabels(), fontsize=8.)
##    plt.savefig('F?_RhoCounterCo_vs_z_BikSucks_nolims.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
#    ax1.set_xlim(-0.5,15)
#    ax1.set_ylim(-0.25,1.25) 
#    ax2.set_xlim(-0.5,15)
##    ax2.set_ylim(0,15) 
##    ax1.legend(loc='best',numpoints=1,prop={"size":10},columnspacing=0.07,borderpad=0.18,labelspacing=0.10,handletextpad=0.07) 
#    fig.set_size_inches(3.37,3.5)
##    if fromlcor == 'yes': #This is plotting theory from lcor onwards
##	    plt.savefig('F?_RhoCounterCo_vs_z_BikSucks_lcor.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
##    else:
##	    plt.savefig('F?_RhoCounterCo_vs_z_BikSucks_zmin.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')      
##    plt.show()
#    plt.close()


#    print 'Plotting 3D plot'
#    fig=plt.figure()
#    ax = Axes3D(fig)
#    ax.scatter(cont_z,cont_phi,cont_dif,marker='o',color='r')
##    ax = fig.add_subplot(111, projection='3d')
##    ax1 = fig.add_subplot(111,projection='3d')
##    for (x,y,z) in zip(cont_z,cont_phi,cont_dif):
###    	col = ROYGBIV_map(z,np.max(cont_dif)*1.2,1)
##    	print x,y,z
##    	ax.scatter(x,y,z,marker='o',color='r')

#    ax.set_xlabel(r'$ z / \sigma$',fontsize=10.)
#    ax.set_ylabel(r'$\Phi[\tilde z]$',fontsize=10.)
#    ax.set_zlabel(r'$ (\mu^{\rm EV}[\Phi] - \mu^{\rm EV}_{\rm Bulk}[\Phi]) / k_{\rm B} T$',fontsize=10.)
#s#    plt.show()
#    plt.savefig('newF_BikSucks_3D.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    plt.close()
        
#    from matplotlib.mlab import griddata
#    plt.figure()
#    cont_phi = np.arange(40)
#    cont_z = np.arange(30)
#    cont_dif = np.linspace(-40,40,100)
#    xi = np.linspace(np.min(cont_phi),np.max(cont_phi),30)
#    yi = np.linspace(np.min(cont_z),np.max(cont_z),30)
#    zi = griddata(cont_phi,cont_z,cont_dif,xi,yi,interp='linear')
#    CS = plt.contour(xi,yi,zi,5,linewidths=0.5,colors='k')
#    CS = plt.contourf(xi,yi,zi,5,cmap=plt.cm.rainbow,vmax=abs(zi).max(), vmin=-abs(zi).max())
#    plt.colorbar()
    
#    matplotlib.rcParams['xtick.direction'] = 'out'
#    matplotlib.rcParams['ytick.direction'] = 'out'
#    CS = plt.contour(cont_phi, cont_z, cont_dif)
#    plt.clabel(CS, inline=1, fontsize=10)

#    plt.savefig('NewF_contour.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   












###    muexEV_wallS,PhiWallS=[],[]



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
    return  ##Above code for P2F1()
#P2F1()



#def  P2tr1_old():
#    """
#	Plotting routines for publication 2, Figure 1 ONLY
#	This is an attempt at adding modularity to my code.

#	This processes old output format and needs to be erased aftert he other is going

##Input: 
##      XXXXfilenameS - A vector of all the filenames to be included in the plots

##Output: Plots as meantioned above.
##Notes: Would be nice to extende Bikerman theory, too.
##"""
#    import numpy as np
#    import matplotlib.pyplot as plt
#    from matplotlib.colors import LinearSegmentedColormap
#    from mpl_toolkits.mplot3d import Axes3D
#    import sys

#    from scipy import optimize
#       
#    #Special groups of filenames
#    filenameS=[]

#    for f in [[x for x in line.split()] for line in file(sys.argv[1],"r").readlines()]:
#    	filenameS.append(f[0])
##    	print 'submit "P1.py %s"\nsleep 2' % str(f[0])
##    filenameS.reverse()
##    import random	
##    random.shuffle(filenameS)

#    GClike=['Analyzed_tr_1004_0.001684_0.1_0.1_502.0_23.73.txt','Analyzed_tr_1016_0.003474_0.1_0.1_502.0_23.73.txt','Analyzed_tr_1048_0.0066065_0.1_0.1_502.0_23.73.txt','Analyzed_tr_1162_0.014195_0.1_0.1_502.0_23.73.txt','Analyzed_tr_1232_0.018603_0.1_0.1_502.0_23.73.txt','Analyzed_tr_32_0.0106795_0.1_0.1_502.0_23.73.txt','Analyzed_tr_626_0.012341_0.1_0.1_502.0_23.73.txt']
#    CSlike=['Analyzed_tr_1016_0.003474_0.1_3.0_502.0_23.73.txt','Analyzed_tr_1066_0.0078345_0.1_3.0_502.0_23.73.txt','Analyzed_tr_1232_0.018603_0.1_3.0_502.0_23.73.txt','Analyzed_tr_1162_0.014195_0.1_3.0_502.0_23.73.txt']
#    print 'Add to these lists once python analysis is fully completed'
# 
##    V_tr=[]
##    E_tr=[]
##    Np_tr=[]
##    Nm_tr=[]
##    tcountS=[]
##    LD_tr=[]
##    EV_tr=[]
#    
#    Fields=[]
#    N_plusS=[]
#    N_minusS=[]
#    LDs=[]
#    BjerrumS=[]
#    sigWCA_S=[]
#    N_totS=[]
#    XiS=[]
#    areaS=[]
#    muexEV_S=[]
#    z_S=[]	#These are the z_pos for plotting is for voltage, field, and muex
#    L_binS=[]
#   # QEDLs=[]

#    print 'Capture all of the QEDL plots and put them onto the same figure'
#    ##This has so far been standard for all transient simulations...
#    Numbins = 502

#    
#   #These variables below are (reasonably) assumed to always equal these values
#    z_wall=1.	
#    temperature = 1.0
#    valency = 1.0
#    eps_wall=1.0
#    t_tot=2510 ##This is the same for everything

#    #This loop build all the necessary data to superimpose plots. 
#    print 'Loading Data for ',len(filenameS),' simulations...'
#    i=0
#    for filename in filenameS:  	
##    	print filename
#    	i+=1

#        #Extract relevant parameters from filename
#        j=0
#        value=''
#        for letter in filename[12:]:
#    		if letter!='_':
#    			value+=letter
#		else:
#			j+=1
#			if j==1:
##				if 'Vtr' in filename[12:]:
##					print 'here'
##				else:
#				N_tot=int(value)
#				N_totS.append(N_tot)
#			if j==2:
#				Xi=float(value)
#				XiS.append(Xi)
#			if j==3:
#				Bjerrum=float(value)
#				BjerrumS.append(Bjerrum)
#			if j==4:
#				sigWCA=float(value)
#				sigWCA_S.append(sigWCA)	
#			if j==5:#This length will actually be determined from LAMMPS files	
#				L_z=float(value) - 2*z_wall
#			value=''
#	L_xy=float(value[:-4])
#	area = L_xy*L_xy
#	areaS.append(area)

#	print 'Before file is loaded'
#    	z_50=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()]) ###This loads EVERYTHING
#    	print 'After file is loaded'
#    
#	#Trim incoming data
#	if i==1: 	##Z positions will always stay constant, so only need to do this once
#		z_positions=z_50[:,0]
#		zpos = z_positions[:(Numbins+1)].tolist()
#		z_density = np.array([(x+y)/2. for (x,y) in zip(zpos[0:len(zpos)-1],zpos[1:len(zpos)])])
#		L_bin=zpos[1]-zpos[0]
#		zpos = np.array(zpos)

#	Psi_tot=z_50[:,1]
#	ilow = 0
#	V_tr=[]
#	while (ilow+Numbins+1)<len(Psi_tot):
#		inst = Psi_tot[ilow:(ilow+Numbins+1)].tolist()
#		V_tr.append(inst)
#		ilow+=(Numbins+1)
#	
#	E_tot=z_50[:,2]
#	ilow = 0
#	E_tr=[]
#	while (ilow+Numbins+1)<len(E_tot):
#		inst = E_tot[ilow:(ilow+Numbins+1)].tolist()
#		E_tr.append(inst)
#		ilow+=(Numbins+1)

#	N_plus=z_50[:,3]
#	ilow = 0
#	Np_tr,LD_tr=[],[]
#	while (ilow+Numbins+1)<len(N_plus):
#		inst = N_plus[ilow:(ilow+Numbins+1)].tolist()
#		Np_tr.append(inst[:-1])
#		LD_tr.append(inst[-1])
#		ilow+=(Numbins+1)	

#	N_minus=z_50[:,4]
#	ilow = 0
#	Nm_tr,tcountS=[],[]
#	while (ilow+Numbins+1)<len(N_minus):
#		inst = N_minus[ilow:(ilow+Numbins+1)].tolist()
#		Nm_tr.append(inst[:-1])
#		tcountS.append(int(inst[-1]))
#		ilow+=(Numbins+1)

#	muex_EV=z_50[:,5]
#	ilow = 0
#	EV_tr=[]
#	while (ilow+Numbins+1)<len(muex_EV):
#		inst = muex_EV[ilow:(ilow+Numbins+1)].tolist()
#		EV_tr.append(inst)
#		ilow+=(Numbins+1)

#    	Nf_tr,Q_tr=[],[]
#    	Nt_tr=[]
#	for (Npt,Nmt) in zip(Np_tr,Nm_tr):
#	    Nf_tr.append(np.array(Npt)-np.array(Nmt))  
#	    Nt_tr.append(np.array(Npt)+np.array(Nmt))
#	    Q_tr.append(sum(Npt[:len(Npt)/2])-sum(Nmt[:len(Nmt)/2]))  	

#    ##Turns out Np and Nm notation is somehow switched... This is fixed below:
##    temp = N_plusS
##    N_plusS = N_minusS
##    N_minusS = temp

##       
#    	print 'Plotting voltage at key times'
#    	fig=plt.figure()
#    	ax1=fig.add_subplot(111)
#	fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    	fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
#    	fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
#    	fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
#    	i=-1
#    	j=-1
#    	maxy = -1e9
#    	colS=[]
#    	percentages = [0.,0.001,0.005,0.01,0.03,0.05,0.1,0.15,0.25,0.70,0.95]
#    	for per in percentages:
#		i+=1
#		if per>0:
#			tindex = int(per*(t_tot-10))+10
#		else:	
#			tindex = 0
#	    	col = ROYGBIV_map(float(i)/len(percentages),1.05,1)
#		colS.append(col)

#		volt_t = V_tr[tindex]
#		LD_t = LD_tr[tindex]

#		maxy = np.max([np.max(volt_t),maxy])
#	    	ax1.errorbar(zpos/LD_t,volt_t,yerr=None,color=col,ls='-',lw=1.10,label=str(per*100)+r'$\%$')

#    	ax1.set_xlabel(r'$z/\lambda_{\rm D}(t)$')
#    	ax1.set_ylabel(r'$\phi / \phi^{\rm T}$')
#    	plt.setp(ax1.get_xticklabels(), fontsize=6.)
#    	plt.setp(ax1.get_yticklabels(), fontsize=6.)

#    	ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
#    	ax1.set_xlim(0,5)       
#    	ax1.set_ylim(-0.5,maxy*1.1)
#    
#    	fig.set_size_inches(3.37,3.5)
#    	graphname = 'Tr_volt.png'
#    	if filename in GClike:
#    		graphname = 'GC_'+str(N_tot)+'_'+graphname
#	elif filename in CSlike:
#    		graphname = 'CS_'+str(N_tot)+'_'+graphname
#    	plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    	plt.close()

#    	print 'Plotting field at key times'
#    	fig=plt.figure()
#    	ax1=fig.add_subplot(111)
#	fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    	fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
#    	fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
#    	fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
#    	i=-1
#    	j=-1
#    	maxy = -1e9
#    	colS=[]
#    	percentages = [0.,0.001,0.005,0.01,0.03,0.05,0.1,0.15,0.25,0.70,0.95]
#    	for per in percentages:
#		i+=1
#		if per>0:
#			tindex = int(per*(t_tot-10))+10
#		else:	
#			tindex = 0
#	    	col = ROYGBIV_map(float(i)/len(percentages),1.05,1)
#		colS.append(col)

#		field_t = E_tr[tindex]
#		LD_t = LD_tr[tindex]
####-np.mean(field_t[Numbins/2-2:Numbins/2+2]
#		maxy = np.max([np.max(field_t),maxy])
#	    	ax1.errorbar(zpos/LD_t,field_t,yerr=None,color=col,ls='-',lw=1.10,label=str(per*100)+r'$\%$')

#    	ax1.set_xlabel(r'$z/\lambda_{\rm D}(t)$')
#    	ax1.set_ylabel(r'$E $')
#    	plt.setp(ax1.get_xticklabels(), fontsize=6.)
#    	plt.setp(ax1.get_yticklabels(), fontsize=6.)

#    	ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
#    	ax1.set_xlim(0,5)       
#    	ax1.set_ylim(0.,maxy*1.1)
#    
#    	fig.set_size_inches(3.37,3.5)
#    	graphname = 'Tr_field.png'
#    	if filename in GClike:
#    		graphname = 'GC_'+str(N_tot)+'_'+graphname
#	elif filename in CSlike:
#    		graphname = 'CS_'+str(N_tot)+'_'+graphname
#    	plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    	plt.close()
#    	
#    	print 'Plotting free charge density at key times'
#    	fig=plt.figure()
#    	ax1=fig.add_subplot(111)
#	fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    	fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
#    	fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
#    	fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
#    	i=-1
#    	j=-1
#    	maxy = -1e9
#    	colS=[]
#    	percentages = [0.,0.001,0.005,0.01,0.03,0.05,0.1,0.15,0.25,0.70,0.95]
#    	for per in percentages:
#		i+=1
#		if per>0:
#			tindex = int(per*(t_tot-10))+10
#		else:	
#			tindex = 0
#	    	col = ROYGBIV_map(float(i)/len(percentages),1.05,1)
#		colS.append(col)

#		FreeChargeDensity = Nf_tr[tindex] / (area*L_bin)
#		LD_t = LD_tr[tindex]
#		ND_free = FreeChargeDensity *4.*np.pi*Bjerrum*LD_t**2
#		maxy = np.max([np.max(ND_free),maxy])
#	    	ax1.errorbar(z_density[1:]/LD_t,ND_free[1:],yerr=None,color=col,ls='-',lw=1.10,label=str(per*100)+r'$\%$')

#    	ax1.set_xlabel(r'$z/\lambda_{\rm D}(t)$')
#    	ax1.set_ylabel(r'$\rho / 2 q en^{\rm B}(t)$')
#    	plt.setp(ax1.get_xticklabels(), fontsize=6.)
#    	plt.setp(ax1.get_yticklabels(), fontsize=6.)

#    	ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
#    	ax1.set_xlim(0,5)       
#    	ax1.set_ylim(-0.5,maxy*1.1)
#    
#    	fig.set_size_inches(3.37,3.5)
#    	graphname = 'Tr_free.png'
#    	if filename in GClike:
#    		graphname = 'GC_'+str(N_tot)+'_'+graphname
#	elif filename in CSlike:
#    		graphname = 'CS_'+str(N_tot)+'_'+graphname
#    	plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    	plt.close()

#    	
#    	print 'Plotting total charge density at key times'
#    	fig=plt.figure()
#    	ax1=fig.add_subplot(111)
#	fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    	fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
#    	fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
#    	fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
#    	i=-1
#    	j=-1
#    	maxy = -1e9
#    	colS=[]
#    	percentages = [0.,0.001,0.005,0.01,0.03,0.05,0.1,0.15,0.25,0.70,0.95]
#    	for per in percentages:
#		i+=1
#		if per>0:
#			tindex = int(per*(t_tot-10))+10
#		else:	
#			tindex = 0
#	    	col = ROYGBIV_map(float(i)/len(percentages),1.05,1)
#		colS.append(col)

#		TotChargeDensity = Nt_tr[tindex] / (area*L_bin)
#		LD_t = LD_tr[tindex]
#		ND_tot = np.array(TotChargeDensity) *4.*np.pi*Bjerrum*LD_t**2
#		maxy = np.max([np.max(ND_tot),maxy])
#	    	ax1.errorbar(z_density[1:]/LD_t,ND_tot[1:],yerr=None,color=col,ls='-',lw=1.10,label=str(per*100)+r'$\%$')

#    	ax1.set_xlabel(r'$z/\lambda_{\rm D}(t)$')
#    	ax1.set_ylabel(r'$\rho_{\rm tot} / 2 q en^{\rm B}(t)$')
#    	plt.setp(ax1.get_xticklabels(), fontsize=6.)
#    	plt.setp(ax1.get_yticklabels(), fontsize=6.)

#    	ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
#    	ax1.set_xlim(0,7)       
#    	ax1.set_ylim(-0.5,maxy*1.1)
#    
#    	fig.set_size_inches(3.37,3.5)
#    	graphname = 'Tr_tot.png'
#    	if filename in GClike:
#    		graphname = 'GC_'+str(N_tot)+'_'+graphname
#	elif filename in CSlike:
#    		graphname = 'CS_'+str(N_tot)+'_'+graphname
#    	plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    	plt.close()


#    	print 'Plotting free charge density at key times normalize by final time'
#    	fig=plt.figure()
#    	ax1=fig.add_subplot(111)
#    	fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    	fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
#    	fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
#    	fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
#    	i=-1
#    	j=-1

#    	colS=[]
#    	percentages = [0.,0.001,0.005,0.01,0.03,0.05,0.1,0.15,0.25,0.70,0.95]
#    	ND_free_final = (np.array(Nf_tr[t_tot-1]) / (area*L_bin))*4.*np.pi*Bjerrum*LD_tr[t_tot-1]**2
#    	for per in percentages:
#		i+=1
#		if per>0:
#			tindex = int(per*(t_tot-10))+10
#		else:	
#			tindex = 0
#	    	col = ROYGBIV_map(float(i)/len(percentages),1.05,1)
#		colS.append(col)
#		FreeChargeDensity = Nf_tr[tindex] / (area*L_bin)
#		LD_t = LD_tr[tindex]
#		ND_free = np.array(FreeChargeDensity) *4.*np.pi*Bjerrum*LD_t**2
#	    	ax1.errorbar(z_density[1:]/LD_t,ND_free[1:]/ND_free_final[1:],yerr=None,color=col,ls='-',lw=1.10,label=str(per*100)+r'$\%$')

#    	ax1.set_xlabel(r'$z/\lambda_{\rm D}(t)$')
#    	ax1.set_ylabel(r'$\rho / \rho(t^{\infty})$')
#    	plt.setp(ax1.get_xticklabels(), fontsize=6.)
#    	plt.setp(ax1.get_yticklabels(), fontsize=6.)

#    	ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
#    	ax1.set_xlim(0,3)       
#    	ax1.set_ylim(-1.75,1.75)
#    
#    	fig.set_size_inches(3.37,3.5)
#    	graphname = 'Tr_freeOverfinal.png'
#    	if filename in GClike:
#    		graphname = 'GC_'+str(N_tot)+'_'+graphname
#	elif filename in CSlike:
#    		graphname = 'CS_'+str(N_tot)+'_'+graphname
#    	plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    	plt.close()

#    	print 'Plotting Q_EDL length as a function of time'
#    	fig=plt.figure()
#    	ax1=fig.add_subplot(111)
#    	fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    	fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
#    	fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
#    	fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
#    	ax1.errorbar(np.array(tcountS)/float(t_tot),Q_tr,yerr=None,color='k',ls='-',lw=1.10)
##    	ax1.errorbar(np.array(tcountS)/float(t_tot),LD_tr,yerr=None,color='red',ls='-',lw=1.10)
#    	ax1.set_xlabel(r'$t_{\rm MD} / t^{\rm \infty} $')
#    	ax1.set_ylabel(r'$Q_{\rm EDL}$')
#    	plt.setp(ax1.get_xticklabels(), fontsize=6.)
#    	plt.setp(ax1.get_yticklabels(), fontsize=6.)
##    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
#    	ax1.set_xlim(0,1)    
#        plt.ylim(ymin=-5)   
##    ax1.set_ylim(-1.75,1.75)
#    	fig.set_size_inches(3.37,3.5)
#    	graphname = 'Tr_Q.png'
#    	if filename in GClike:
#    		graphname = 'GC_'+str(N_tot)+'_'+graphname
#	elif filename in CSlike:
#    		graphname = 'CS_'+str(N_tot)+'_'+graphname
#    	plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    	plt.close()

#    	print 'Plotting Q_EDL length as a function of time'
#    	fig=plt.figure()
#    	ax1=fig.add_subplot(111)
#    	fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    	fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
#    	fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
#    	fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
##    	ax1.errorbar(np.array(tcountS)/float(t_tot),Q_tr,yerr=None,color='k',ls='-',lw=1.10)
#    	ax1.errorbar(np.array(tcountS)/float(t_tot),LD_tr,yerr=None,color='red',ls='-',lw=1.10)
#    	ax1.set_xlabel(r'$t_{\rm MD} / t^{\rm \infty} $')
#    	ax1.set_ylabel(r'$\lambda_{\rm D}(\tilde t)$')
#    	plt.setp(ax1.get_xticklabels(), fontsize=6.)
#    	plt.setp(ax1.get_yticklabels(), fontsize=6.)
##    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
#    	ax1.set_xlim(0,1)     
##    ax1.set_ylim(-1.75,1.75)
#    	fig.set_size_inches(3.37,3.5)
#    	graphname = 'Tr_lamD.png'
#    	if filename in GClike:
#    		graphname = 'GC_'+str(N_tot)+'_'+graphname
#	elif filename in CSlike:
#    		graphname = 'CS_'+str(N_tot)+'_'+graphname
#    	plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    	plt.close()
#    
#    return  ##Above code for P2tr1()
##P2tr1_old()

def PlotSratio():
    import numpy as np
#    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.mlab import griddata
    from matplotlib.ticker import MaxNLocator,FormatStrFormatter
    from mpl_toolkits.mplot3d import Axes3D
    import sys
    from scipy import interpolate
    from scipy import optimize
#    from scipy.interpolate import UnivariateSpline
    from scipy.optimize import fsolve
#####    from scipy.optimize import root
    from scipy import stats

    print 'Plotting S ratios'
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
    lo,mi,hi=-1,-1,-1
#    PhiLabel = [r'$0.01$',r'$0.05$',r'$0.07$',r'$0.10$',r'$0.20$',r'$0.30$',r'$0.40$',r'$0.50$']
    PhiLabel = [r'$0.01$',r'$0.07$',r'$0.10$',r'$0.20$',r'$0.30$',r'$0.40$',r'$0.50$']
    for filename in filenameS:
    	i+=1
    	if i==len(PhiLabel):
    		i=0

    	print filename
#       	MMA=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()])

       	MMA=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()])



	print filename[4:8]
	
	if filename[4:8] == '0.05':
		ls_x = '--'
		lo+=1
	elif filename[4:8] == '0.1_':
		ls_x='-'
		mi+=1
	elif filename[4:8] == '0.15':
		ls_x='-.'
		hi+=1

	label_x = PhiLabel[i]
	
	if lo==0:
       		ax1.errorbar([6.1,6.1],[0.8,0.8],yerr=None,ls=ls_x,lw=1.0,color='white',marker='None',label=r'$\Phi_{\pm}=0.05$')
       		label_x = r'$\Phi_{\rm sol}=$'+	label_x	
	if mi==0:
       		ax1.errorbar([6.1,6.1],[0.8,0.8],yerr=None,ls=ls_x,lw=1.0,color='white',marker='None',label=r'$\Phi_{\pm}=0.10$')
       		label_x = r'$\Phi_{\rm sol}=$'+	label_x	
	if hi==0:
       		ax1.errorbar([6.1,6.1],[0.8,0.8],yerr=None,ls=ls_x,lw=1.0,color='white',marker='None',label=r'$\Phi_{\pm}=0.15$')       		
       		label_x = r'$\Phi_{\rm sol}=$'+	label_x	
	value = ''
	for letter in filename[9:]:
		if letter is 't':
			break
		value+=letter
	PhiSol = value[:-1]
       	col = ROYGBIV_map(float(PhiSol),0.75,1)
       	colS.append(col)

       	ax1.errorbar(MMA[:,2],MMA[:,3],yerr=None,ls=ls_x,lw=1.0,color=col,marker='None',label=label_x)


#    if hi!=0:
#    	ncols = 3
#    elif mi!=0:
#    	ncols = 2
#    else:
#    	ncols = 1
    	
    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=2,handletextpad=0.1) 
    ax1.set_xlabel(r'$S_{\rm CS}$',size='small')
    ax1.set_ylabel(r'$S_{\rm MFS}/S_{\rm CS}$',size='small')
    plt.setp(ax1.get_xticklabels(), fontsize='small')
    plt.setp(ax1.get_yticklabels(), fontsize='small')
#    ax1.xaxis.set_major_locator(MaxNLocator(6))
#    ax1.yaxis.set_major_locator(MaxNLocator(6))
    ax1.set_xlim(0,6)     
    plt.ylim(ymax=1.05)
#    ax1.set_ylim(-,1.75)
    fig.set_size_inches(3.37,3.5)
    graphname = 'Sratios.png'
    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    plt.close() 

    return
#PlotSratio()


def PlotMFSandCS():
    import numpy as np
#    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.mlab import griddata
    from matplotlib.ticker import MaxNLocator,FormatStrFormatter
    from mpl_toolkits.mplot3d import Axes3D
    import sys
    from scipy import interpolate
    from scipy import optimize
#    from scipy.interpolate import UnivariateSpline
    from scipy.optimize import fsolve
#####    from scipy.optimize import root
    from scipy import stats

    print 'Plotting MFS & CS ratios'
    filenameS=[]
    for f in [[x for x in line.split()] for line in file(sys.argv[1],"r").readlines()]:
    	filenameS.append(f[0])

#    filenameS=['MFS_0.05_0.2.txt']

    distanceS = []
    plusS,minusS = [],[]
    freeS = []
    voltageS = []
    SigEffS = []
    S_ofzS = []
    dCapS = []
    solventS=[]
    PhipmS=[]
    PhisolS=[]
    colS=[]
    ls_xS=[]
    labelS=[]

    i=-1
    lo,mi,hi=-1,-1,-1
    PhiLabel = [r'$0.01$',r'$0.07$',r'$0.10$',r'$0.20$',r'$0.30$',r'$0.40$',r'$0.50$']

    for filename in filenameS:
    	i+=1
    	if i==len(PhiLabel):
    		i=0

    	print filename
#       	MMA=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()])

       	MMA=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()])

       	distanceS.append(MMA[:,0])
       	plusS.append(MMA[:,1])
       	minusS.append(MMA[:,2])
       	freeS.append(MMA[:,3])
       	voltageS.append(MMA[:,4])
       	SigEffS.append(MMA[:,5])
       	S_ofzS.append(MMA[:,6])
       	dCapS.append(MMA[:,7])
 

	if filename[:2] == 'CS':
		ls_xS.append('--')
		#Extract relevant parameters from filename
		j=0
		value=''
		for letter in filename[8:]:
	    		if letter!='_':
	    			value+=letter
			else:
				j+=1
				if j==1:
#					Phiion = float(value)
#					PhipmS.append(Phiion)
					value=''
#		Phisolvent = float(value[:-4])
#		PhisolS.append(Phisolvent)
#		colS.append(ROYGBIV_map(float(0),0.75,1))	
#		labelS.append(str(Phiion))
	elif filename[:4] == 'MFS_':
		ls_xS.append('-')
	       	solventS.append(MMA[:,8])
	       	j=0
		value=''
		for letter in filename[4:]:
	    		if letter!='_':
	    			value+=letter
			else:
				j+=1
				if j==1:
#					Phiion = float(value)
#					PhipmS.append(Phiion)
					value=''
#		Phisolvent = float(value[:-4])
#		PhisolS.append(Phisolvent)
#		colS.append(ROYGBIV_map(float(Phisolvent),0.75,1))	
#		labelS.append(str(Phiion)+','+str(Phisolvent))
	elif filename[:4] == 'MFSS':
		ls_xS.append('-.')
	       	solventS.append(MMA[:,8])
	       	j=0
		value=''
		for letter in filename[5:]:
	    		if letter!='_':
	    			value+=letter
			else:
				j+=1
				if j==1:
#					Phiion = float(value)
#					PhipmS.append(Phiion)
					value=''
				if j==2:
#					Phisolvent = float(value)
#					PhisolS.append(Phisolvent)
					value=''
		w = float(value[1:-4])				
#		colS.append(ROYGBIV_map(float(Phisolvent),0.75,1))	
#		labelS.append(str(Phiion)+','+str(Phisolvent))
					
#	###rhoCS versus rhoMFS
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
#    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
#    i=-1
#    i=-1
#    for (free,S,filename,col,ls_x,label_x) in zip(freeS,S_ofzS,filenameS,colS,ls_xS,labelS):  
#	i+=1
#   	ax1.errorbar(freeS[0],free,yerr=None,ls=ls_x,lw=1.0,color=col,marker='None',label=label_x)	
#    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
#    ax1.set_ylabel(r'$\tilde \rho_{\rm MFS}$',size='small')
#    ax1.set_xlabel(r'$\tilde \rho_{\rm CS}(\Phi_{\pm}=0.15)$',size='small')
#    plt.setp(ax1.get_xticklabels(), fontsize='small')
#    plt.setp(ax1.get_yticklabels(), fontsize='small')
##    ax1.xaxis.set_major_locator(MaxNLocator(6))
##    ax1.yaxis.set_major_locator(MaxNLocator(6))
##    ax1.set_xlim(0,1.5)     
#    plt.xlim(xmax=0)
##    ax1.set_ylim(0,1.5)
#    fig.set_size_inches(3.37,3.5)
#    graphname = 'rho_vs_rho.png'
#    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    plt.close() 


	##all versus z
    print 'At this plot --- '
    fig=plt.figure()
    print 'now here'
    print len(distanceS)
    
    ax1=fig.add_subplot(211)
    ax2=fig.add_subplot(212)
    i=-1
    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='-',lw=1.0,color='k',marker='None',label=r'${\rm CS-Solvent}$')
    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls=':',lw=1.0,color='k',marker='None',label=r'${\rm CS}$')



    ax2.errorbar([-1,-1],[-1,-1],yerr=None,ls='-',lw=1.0,color='orange',marker='None',label=r'$n_{\rm Co-ion}$')
    ax2.errorbar([-1,-1],[-1,-1],yerr=None,ls='-',lw=1.0,color='blue',marker='None',label=r'$n_{\rm Counter-ion}$')
    ax2.errorbar([-1,-1],[-1,-1],yerr=None,ls='-',lw=1.0,color='grey',marker='None',label=r'$n_{\rm Solvent}$')


    print 'now here'
    print len(distanceS)
    
    for (distance,plus,filename,minus,voltage,free,sig) in zip(distanceS,plusS,filenameS,minusS,voltageS,freeS,SigEffS):  
	i+=1
	print len(distance)
	print filename

	if 'CS' in filename:
		ls_x = ':'
   		ax1.errorbar(voltage,sig,yerr=None,ls=ls_x,lw=1.0,color='k',marker='None')
   		

   		ax2.errorbar(distance,plus,yerr=None,ls=ls_x,lw=1.0,color='orange',marker='None')
   		ax2.errorbar(distance,minus,yerr=None,ls=ls_x,lw=1.0,color='blue',marker='None')
	else:
		ls_x = '-'
   		ax1.errorbar(voltage,sig,yerr=None,ls=ls_x,lw=1.0,color='k',marker='None')
   		
   		ax2.errorbar(distance,plus,yerr=None,ls=ls_x,lw=1.0,color='orange',marker='None')
   		ax2.errorbar(distance,minus,yerr=None,ls=ls_x,lw=1.0,color='blue',marker='None')
#   		ax1.errorbar(distance,-np.array(free),yerr=None,ls=ls_x,lw=1.0,color='grey',marker='None')#,label=r'${\rm |Free charge|}$')	
   		ax2.errorbar(distance,solventS[0],yerr=None,ls=ls_x,lw=1.0,color='grey',marker='None')	

    ax1.legend(loc='best',numpoints=1,prop=dict(size=8.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
    ax2.legend(loc='best',numpoints=1,prop=dict(size=8.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 

    ax1.set_xlabel(r'$\phi / \phi_{\rm T}$',fontsize=10.)
    ax1.set_ylabel(r'$\Sigma / \Sigma_{\rm ref}$',fontsize=10.)

    ax2.set_xlabel(r'$n_{i} / n_{\rm B}$',fontsize=10.)
    ax2.set_ylabel(r'$\Sigma / \Sigma_{\rm ref}$',fontsize=10.)
    
    plt.setp(ax1.get_xticklabels(), fontsize=8.)
    plt.setp(ax1.get_yticklabels(), fontsize=8.)
    plt.setp(ax2.get_xticklabels(), fontsize=8.)
    plt.setp(ax2.get_yticklabels(), fontsize=8.)
#    ax1.set_xlim(0,1.5)     
#    plt.xlim(xmin=0)
#    plt.ylim(ymin=-0.1)

    ax2.set_xlim(0,10)
    ax2.set_ylim(-0.1,10)
    fig.set_size_inches(5.0,4.0)
    graphname = 'MFS_Comparison.pdf'
    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight') 
    plt.close() 
    
    print fail


	###all versus z
    fig=plt.figure()
    ax1=fig.add_subplot(121)
    ax2=fig.add_subplot(122)
    i=-1
    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='--',lw=1.0,color='k',marker='None',label=r'$\Phi_{\pm}=0.15$')
    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='-',lw=1.0,color='k',marker='None',label=r'$(\Phi_{\pm},\Phi_{\rm sol})=(0.15,0.50)$')
    for (distance,plus,filename,col,ls_x,label_x,minus,voltage,free) in zip(distanceS,plusS,filenameS,colS,ls_xS,labelS,minusS,voltageS,freeS):  
	i+=1
	if i==0:	
	   	ax1.errorbar(distance,voltage,yerr=None,ls=ls_x,lw=1.0,color='blue',marker='None')#	,label=r'$Voltage$')	
#   		ax1.errorbar(distance,plus,yerr=None,ls=ls_x,lw=1.0,color='green',marker='None')
#   		ax1.errorbar(distance,minus,yerr=None,ls=ls_x,lw=1.0,color='red',marker='None')
   		ax1.errorbar(distance,-np.array(free),yerr=None,ls=ls_x,lw=1.0,color='grey',marker='None')
	else:
	   	ax1.errorbar(distance,voltage,yerr=None,ls=ls_x,lw=1.0,color='blue',marker='None',label=r'${\rm Voltage}$')	
#   		ax1.errorbar(distance,plus,yerr=None,ls=ls_x,lw=1.0,color='green',marker='None',label=r'${\rm Co-ion}$')	
#   		ax1.errorbar(distance,minus,yerr=None,ls=ls_x,lw=1.0,color='red',marker='None',label=r'${\rm Counter-ion}$')	
   		ax1.errorbar(distance,-np.array(free),yerr=None,ls=ls_x,lw=1.0,color='grey',marker='None',label=r'${\rm |Free Charge|}$')	
#   		ax1.errorbar(distance,solventS[0],yerr=None,ls=ls_x,lw=1.0,color='purple',marker='None',label=r'${\rm Solvent}$')	
    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
#    ax1.set_ylabel(r'$\tilde \psi,\tilde \rho_{\rm \pm}, \tilde \rho_{\rm sol}$',size='small')
    ax1.set_ylabel(r'$\phi / \phi_{\rm T}$',fontsize=10.)
    ax1.set_xlabel(r'$z / \lambda_{\rm D}$',fontsize=10.)
    
    plt.setp(ax1.get_xticklabels(), fontsize=8.)
    plt.setp(ax1.get_yticklabels(), fontsize=8.)
    plt.setp(ax2.get_xticklabels(), fontsize=8.)
    plt.setp(ax2.get_yticklabels(), fontsize=8.)


#    ax1.set_xlim(0,1.5)     
    plt.xlim(xmin=0)
    plt.ylim(ymin=-0.1)
#    ax1.set_ylim(0,1.5)
    fig.set_size_inches(3.37,4.0)
#    graphname = 'freeVolt_vs_distance.png'
    plt.savefig('FREEVOLT_profiles_MFS.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
    plt.close() 



	###rho vs S
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    i=-1
    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='--',lw=1.0,color='white',marker='None',label=r'$\Phi_{\pm},\Phi_{\rm sol}$')
    for (free,S,filename,col,ls_x,label_x) in zip(freeS,S_ofzS,filenameS,colS,ls_xS,labelS):  
	i+=1
   	ax1.errorbar(S,-np.array(free),yerr=None,ls=ls_x,lw=1.0,color=col,marker='None',label=label_x)	
    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
    ax1.set_ylabel(r'$-\tilde \rho$',size='small')
    ax1.set_xlabel(r'$S_{\rm MFS}$',size='small')
    plt.setp(ax1.get_xticklabels(), fontsize='small')
    plt.setp(ax1.get_yticklabels(), fontsize='small')
#    ax1.xaxis.set_major_locator(MaxNLocator(6))
#    ax1.yaxis.set_major_locator(MaxNLocator(6))
#    ax1.set_xlim(0,1.5)     
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
#    ax1.set_ylim(0,1.5)
    fig.set_size_inches(3.37,3.5)
    graphname = 'rho_vs_S.png'
    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    plt.close() 

	###dCap versus S
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    i=-1
#    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='--',lw=1.0,color='white',marker='None',label=r'$\Phi_{\pm}$')
    i=-1
    for (dC,S,filename,col,ls_x,label_x) in zip(dCapS,S_ofzS,filenameS,colS,ls_xS,labelS):  
	i+=1
#	if i==0:
#		label_x='0,1'
#	if i==1:
#		label_x='0.15,0'		
#	if i==9:
#	    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='-',lw=1.0,color='white',marker='None',label=r'$\Phi_{\pm},\Phi_{\rm sol}$')
   	ax1.errorbar(S,dC,yerr=None,ls=ls_x,lw=1.0,color=col,marker='None',label='0,1')	
   	
    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
    ax1.set_xlabel(r'$S$',size='small')
    ax1.set_ylabel(r'$C\prime (\lambda_{\rm D}/\epsilon)$',size='small')
    plt.setp(ax1.get_xticklabels(), fontsize='small')
    plt.setp(ax1.get_yticklabels(), fontsize='small')
#    ax1.xaxis.set_major_locator(MaxNLocator(6))
#    ax1.yaxis.set_major_locator(MaxNLocator(6))
#    ax1.set_,xlim(0,6)     
    plt.xlim(xmin=0)
    ax1.set_ylim(0,1.5)
    fig.set_size_inches(3.37,3.5)
    graphname = 'dCap_vs_S.png'
    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    plt.close() 


	###dCap versus voltage
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    i=-1
#    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='--',lw=1.0,color='white',marker='None',label=r'$\Phi_{\pm},\Phi_{\rm sol}$')
    ax1.errorbar(voltage,np.cosh(np.array(voltage)/2),yerr=None,ls='-',lw=1.0,color='k',marker='None',label=r'$\Phi_{\pm} \to 0,\Phi_{\rm sol \to 1$')
    i=-1
    for (dC,voltage,filename,col,ls_x,label_x) in zip(dCapS,voltageS,filenameS,colS,ls_xS,labelS):  
	i+=1
#	if i==9:
#	    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='-',lw=1.0,color='white',marker='None',label=r'$\Phi_{\pm},\Phi_{\rm sol}$')

	if i==0:
		label_x='0.15,0'

	print filename,i
   	ax1.errorbar(voltage,dC,yerr=None,ls=ls_x,lw=1.0,color=col,marker='None',label=label_x)	
   	
    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
    ax1.set_xlabel(r'$\tilde \psi$',size='small')
    ax1.set_ylabel(r'$C^{\prime} (\lambda_{\rm D}/\epsilon)$',size='small')
    plt.setp(ax1.get_xticklabels(), fontsize='small')
    plt.setp(ax1.get_yticklabels(), fontsize='small')
#    ax1.xaxis.set_major_locator(MaxNLocator(6))
#    ax1.yaxis.set_major_locator(MaxNLocator(6))
#    ax1.set_,xlim(0,6)     
    plt.xlim(xmin=0)
    ax1.set_ylim(0,1.5)
    fig.set_size_inches(3.37,3.5)
    graphname = 'dCap_vs_voltage.png'
    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    plt.close() 




    return
#PlotMFSandCS()

def  P2tr1():
    """
	Plotting routines for publication 2, Figure 1 ONLY
	This is an attempt at adding modularity to my code.
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
    filenameS=[]

    for f in [[x for x in line.split()] for line in file(sys.argv[1],"r").readlines()]:
    	filenameS.append(f[0])

#    filenameS.reverse()
#    import random	
#    random.shuffle(filenameS)
	
	
#    filenameS = ['Analyzed_tV_1048_0.0066065_0.1_0.3_502.0_23.73.txt']
    
    GClike=['Analyzed_tV_1016_0.003474_0.1_0.3_502.0_23.73.txt','Analyzed_tV_1034_0.005482_0.1_0.3_502.0_23.73.txt','Analyzed_tV_1048_0.0066065_0.1_0.3_502.0_23.73.txt','Analyzed_tV_1066_0.0078345_0.1_0.3_502.0_23.73.txt','Analyzed_tV_1162_0.014195_0.1_0.3_502.0_23.73.txt','Analyzed_tV_1232_0.018603_0.1_0.3_502.0_23.73.txt','Analyzed_tV_1312_0.024179_0.1_0.3_502.0_23.73.txt']
    CSlike=['Analyzed_tV_1004_0.001684_0.1_3.0_502.0_23.73.txt','Analyzed_tV_1016_0.003474_0.1_3.0_502.0_23.73.txt','Analyzed_tV_1034_0.005482_0.1_3.0_502.0_23.73.txt','Analyzed_tV_1048_0.0066065_0.1_3.0_502.0_23.73.txt','Analyzed_tV_1066_0.0078345_0.1_3.0_502.0_23.73.txt','Analyzed_tV_1104_0.0106795_0.1_3.0_502.0_23.73.txt','Analyzed_tV_1136_0.012341_0.1_3.0_502.0_23.73.txt','Analyzed_tV_1162_0.014195_0.1_3.0_502.0_23.73.txt','Analyzed_tV_1232_0.018603_0.1_3.0_502.0_23.73.txt','Analyzed_tV_1312_0.024179_0.1_3.0_502.0_23.73.txt']

    print 'Add to these lists once python analysis is fully completed'
 
#    V_tr=[]
#    E_tr=[]
#    Np_tr=[]
#    Nm_tr=[]
#    tcountS=[]
#    LD_tr=[]
#    EV_tr=[]
    
    Fields=[]
    N_plusS=[]
    N_minusS=[]
    LDs=[]
    BjerrumS=[]
    sigWCA_S=[]
    N_totS=[]
    XiS=[]
    areaS=[]
    muexEV_S=[]
    vzMS,vzPS=[],[]  
  
    z_S=[]	#These are the z_pos for plotting is for voltage, field, and muex
    L_binS=[]
   # QEDLs=[]



    print 'Capture all of the QEDL plots and put them onto the same figure'
    ##This has so far been standard for all transient simulations...
    Numbins = 502

#    Numbins = 250
    
   #These variables below are (reasonably) assumed to always equal these values
    z_wall=1.	
    temperature = 1.0
    valency = 1.0
    eps_wall=1.0
    t_tot=2510 ##This is the same for everything

    #This loop build all the necessary data to superimpose plots. 
    print 'Loading Data for ',len(filenameS),' simulations...'
    i=0
    for filename in filenameS:  
    	if filename in CSlike:
    		tauc = 4600.
	elif filename in GClike:
		tauc = 	4090.
	else:
		tauc = 6685.

    	i+=1

        #Extract relevant parameters from filename
        j=0
        value=''
	print 'here is a temp change'
        for letter in filename[12:]:	#This is the more general case
#        for letter in filename[14:]:  ##USe this for specific files,
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
				sigWCA_S.append(sigWCA)	
			if j==5:#This length will actually be determined from LAMMPS files	
				L_z=float(value) - 2*z_wall
			value=''
	L_xy=float(value[:-4])
	area = L_xy*L_xy
	areaS.append(area)

	print 'Before file is loaded'
    	z_50=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()]) ###This loads EVERYTHING
    	print 'After file is loaded'
    
	#Trim incoming data
	if i==1: 	##Z positions will always stay constant, so only need to do this once
		z_positions=z_50[:,0]
		zpos = z_positions[:(Numbins+1)].tolist()
		z_density = np.array([(x+y)/2. for (x,y) in zip(zpos[0:len(zpos)-1],zpos[1:len(zpos)])])
		L_bin=zpos[1]-zpos[0]
		zpos = np.array(zpos)

	Psi_tot=z_50[:,1]
	ilow = 0
	V_tr=[]
	while (ilow+Numbins+1)<len(Psi_tot):
		inst = Psi_tot[ilow:(ilow+Numbins+1)].tolist()
		V_tr.append(inst)
		ilow+=(Numbins+1)
	
	E_tot=z_50[:,2]
	ilow = 0
	E_tr=[]
	while (ilow+Numbins+1)<len(E_tot):
		inst = E_tot[ilow:(ilow+Numbins+1)].tolist()
		E_tr.append(inst)
		ilow+=(Numbins+1)

	N_plus=z_50[:,3]
	ilow = 0
	Np_tr,LD_tr=[],[]
	while (ilow+Numbins+1)<len(N_plus):
		inst = N_plus[ilow:(ilow+Numbins+1)].tolist()
		Np_tr.append(inst[:-1])
		LD_tr.append(inst[-1])
		ilow+=(Numbins+1)	

	N_minus=z_50[:,4]
	ilow = 0
	Nm_tr,tcountS=[],[]
	while (ilow+Numbins+1)<len(N_minus):
		inst = N_minus[ilow:(ilow+Numbins+1)].tolist()
		Nm_tr.append(inst[:-1])
		tcountS.append(int(inst[-1]))
		ilow+=(Numbins+1)

	muex_EV=z_50[:,5]
	ilow = 0
	EV_tr=[]
	while (ilow+Numbins+1)<len(muex_EV):
		inst = muex_EV[ilow:(ilow+Numbins+1)].tolist()
		EV_tr.append(inst)
		ilow+=(Numbins+1)

    	Nf_tr,Q_tr=[],[]
    	Nt_tr=[]
	for (Npt,Nmt) in zip(Np_tr,Nm_tr):
	    Nf_tr.append(np.array(Npt)-np.array(Nmt))  
	    Nt_tr.append(np.array(Npt)+np.array(Nmt))
	    Q_tr.append(sum(Npt[:len(Npt)/2])-sum(Nmt[:len(Nmt)/2]))  	


	VzP=z_50[:,6]
	ilow = 0
	vzP_tr=[]
	QL_tr=[]
	while (ilow+Numbins+1)<len(VzP):
		inst = VzP[ilow:(ilow+Numbins+1)].tolist()
		vzP_tr.append(inst[:-1])
		QL_tr.append(inst[-1])
		ilow+=(Numbins+1)

	VzM=z_50[:,7]
	ilow = 0
	vzM_tr=[]
	QR_tr=[]
	while (ilow+Numbins+1)<len(VzM):
		inst = VzM[ilow:(ilow+Numbins+1)].tolist()
		vzM_tr.append(inst[:-1])
		QR_tr.append(inst[-1])
		ilow+=(Numbins+1)


    ##Turns out Np and Nm notation is somehow switched... This is fixed below:
#    temp = N_plusS
#    N_plusS = N_minusS
#    N_minusS = temp

#       
#    	print 'Plotting voltage at key times'
#    	fig=plt.figure()
#    	ax1=fig.add_subplot(111)
#	fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    	fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
#    	fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
#    	fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
#    	i=-1
#    	j=-1
#    	maxy = -1e9
#    	colS=[]
#    	percentages = [0.,0.001,0.005,0.01,0.03,0.05,0.1,0.15,0.25,0.70,0.95]
#    	for per in percentages:
#		i+=1
#		if per>0:
#			tindex = int(per*(t_tot-10))+10
#		else:	
#			tindex = 0
#	    	col = ROYGBIV_map(float(i)/len(percentages),1.05,1)
#		colS.append(col)

#		volt_t = V_tr[tindex]
#		LD_t = LD_tr[tindex]

#		maxy = np.max([np.max(volt_t),maxy])
#	    	ax1.errorbar(zpos/LD_t,volt_t,yerr=None,color=col,ls='-',lw=1.10,label=str(per*100)+r'$\%$')

#    	ax1.set_xlabel(r'$z/\lambda_{\rm D}(t)$')
#    	ax1.set_ylabel(r'$\phi / \phi^{\rm T}$')
#    	plt.setp(ax1.get_xticklabels(), fontsize=6.)
#    	plt.setp(ax1.get_yticklabels(), fontsize=6.)

#    	ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
#    	ax1.set_xlim(0,5)       
#    	ax1.set_ylim(-0.5,maxy*1.1)
#    
#    	fig.set_size_inches(3.37,3.5)
#    	graphname = 'Tr_volt.png'
#    	if filename in GClike:
#    		graphname = 'GC_'+str(N_tot)+'_'+graphname
#	elif filename in CSlike:
#    		graphname = 'CS_'+str(N_tot)+'_'+graphname
#    	plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    	plt.close()


	##Name string
#	name_string = str(N_tot)
	name_string = str(filename[:3])

    	print 'Plotting field at key times'
    	fig=plt.figure()
    	ax1=fig.add_subplot(111)
	fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    	fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    	fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    	fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    	i=-1
    	j=-1
    	maxy = -1e9
    	colS=[]
    	percentages = [0.,0.001,0.005,0.01,0.03,0.05,0.1,0.15,0.25,0.70,0.95]
    	#percentages = [0.,0.001,0.005,0.01,0.03]#,0.05,0.1,0.15,0.25,0.70,0.95]

    	for per in percentages:
		i+=1
		if per>0:
			tindex = int(per*(t_tot-10))+10
		else:	
			tindex = 0
	    	col = ROYGBIV_map(float(i)/len(percentages),1.05,1)
		colS.append(col)

		field_t = E_tr[tindex]
		LD_t = LD_tr[tindex]
###-np.mean(field_t[Numbins/2-2:Numbins/2+2]
		maxy = np.max([np.max(field_t),maxy])
	    	ax1.errorbar(zpos,field_t,yerr=None,color=col,ls='-',lw=1.10,label=str(round(100*tindex/tauc,1)))


    	    	    	
#    	for (z,V,E,Np,Nm,muexEV,vA,vC) in zip(Shell_bin,Psi_tr[tcount],Efield_tr[tcount],A_count.tolist()+[tcount],C_count.tolist()+[lambda_D],muex_EV_tr[tcount],vzA.tolist()+[ql],vzC.tolist()+[qr]):
#		total_output.write("%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\t\t\t\t%-1.8f\n" % (z,V,E,Nm,Np,muexEV,vA,vC))

    	ax1.set_xlabel(r'$z$')
    	ax1.set_ylabel(r'$E $')
    	plt.setp(ax1.get_xticklabels(), fontsize=6.)
    	plt.setp(ax1.get_yticklabels(), fontsize=6.)

    	ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
#    	ax1.set_xlim(0,5)       
    	ax1.set_ylim(0.,maxy*1.1)
    
    	fig.set_size_inches(3.37,3.5)
    	graphname = 'Tr_field.png'
    	if filename in GClike:
    		graphname = 'GC_'+name_string+'_'+graphname
	elif filename in CSlike:
    		graphname = 'CS_'+name_string+'_'+graphname
	else:
    		graphname = 'temp_'+name_string+'_'+graphname
    	plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    	plt.close()


    	print 'Plotting overlaid field at key times'
    	fig=plt.figure()
    	ax1=fig.add_subplot(111)
	fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    	fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    	fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    	fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    	i=-1
    	j=-1
    	maxy = -1e9
    	colS=[]
    	percentages = [0.,0.001,0.005,0.01,0.03,0.05,0.1,0.15,0.25,0.70,0.95]
    	#percentages = [0.,0.001,0.005,0.01,0.03]#,0.05,0.1,0.15,0.25,0.70,0.95]

    	for per in percentages:
		i+=1
		if per>0:
			tindex = int(per*(t_tot-10))+10
		else:	
			tindex = 0
	    	col = ROYGBIV_map(float(i)/len(percentages),1.05,1)
		colS.append(col)

		field_t = E_tr[tindex]
		LD_t = LD_tr[tindex]
###-np.mean(field_t[Numbins/2-2:Numbins/2+2]
		maxy = np.max([np.max(field_t),maxy])

		zposL = zpos[:len(zpos)/2]
		field_tL = field_t[:len(field_t)/2]
		field_tR = field_t[len(field_t)/2+1:]
#		print len(field_tR),len(field_tL),len(zposL)
		field_tR.reverse()				
		
	    	ax1.errorbar(zposL,field_tL,yerr=None,color=col,ls='-',lw=1.0,label=str(round(100*tindex/tauc,1)))
	    	ax1.errorbar(zposL,field_tR,yerr=None,color=col,ls='-',lw=1.50)#,label=str(round(100*tindex/tauc,1)))

    	ax1.set_xlabel(r'$z$')
    	ax1.set_ylabel(r'$E $')
    	plt.setp(ax1.get_xticklabels(), fontsize=6.)
    	plt.setp(ax1.get_yticklabels(), fontsize=6.)

    	ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
#    	ax1.set_xlim(0,5)       
    	ax1.set_ylim(0.,maxy*1.1)
    
    	fig.set_size_inches(3.37,3.5)
    	graphname = 'Tr_fieldoverlay.png'
    	if filename in GClike:
    		graphname = 'GC_'+name_string+'_'+graphname
	elif filename in CSlike:
    		graphname = 'CS_'+name_string+'_'+graphname
	else:
    		graphname = 'temp_'+name_string+'_'+graphname
    	plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    	plt.close()


  	print 'Plotting field AVERAGE at key times'
    	fig=plt.figure()
    	ax1=fig.add_subplot(111)
	fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    	fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    	fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    	fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    	i=-1
    	j=-1
    	maxy = -1e9
    	colS=[]
    	percentages = [0.,0.001,0.005,0.01,0.03,0.05,0.1,0.15,0.25,0.70,0.95]
    	#percentages = [0.,0.001,0.005,0.01,0.03]#,0.05,0.1,0.15,0.25,0.70,0.95]

    	for per in percentages:
		i+=1
		if per>0:
			tindex = int(per*(t_tot-10))+10
		else:	
			tindex = 0
	    	col = ROYGBIV_map(float(i)/len(percentages),1.05,1)
		colS.append(col)

		field_t = E_tr[tindex]
		LD_t = LD_tr[tindex]
###-np.mean(field_t[Numbins/2-2:Numbins/2+2]
		maxy = np.max([np.max(field_t),maxy])

		zposL = zpos[:len(zpos)/2]
		field_tL = field_t[:len(field_t)/2]
		field_tR = field_t[len(field_t)/2+1:]
#		print len(field_tR),len(field_tL),len(zposL)
		field_tR.reverse()				
		
	    	ax1.errorbar(zposL,[(x+y)*0.5 for (x,y) in zip(field_tL,field_tR)],yerr=None,color=col,ls='-',lw=1.0,label=str(round(100*tindex/tauc,1)))

    	ax1.set_xlabel(r'$z$')
    	ax1.set_ylabel(r'$E $')
    	plt.setp(ax1.get_xticklabels(), fontsize=6.)
    	plt.setp(ax1.get_yticklabels(), fontsize=6.)

    	ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
#    	ax1.set_xlim(0,5)       
    	ax1.set_ylim(0.,maxy*1.1)
    
    	fig.set_size_inches(3.37,3.5)
    	graphname = 'Tr_fieldoverlayAVG.png'
    	if filename in GClike:
    		graphname = 'GC_'+name_string+'_'+graphname
	elif filename in CSlike:
    		graphname = 'CS_'+name_string+'_'+graphname
	else:
    		graphname = 'temp_'+name_string+'_'+graphname
    	plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    	plt.close()




    	print 'Plotting free charge density at key times'
    	fig=plt.figure()
    	ax1=fig.add_subplot(111)
	fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    	fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    	fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    	fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    	i=-1
    	j=-1
    	maxy = -1e9
    	colS=[]
    	percentages = [0.,0.001,0.005,0.01,0.03,0.05,0.1,0.15,0.25,0.70,0.95]
#    	percentages = [0.,0.001,0.005,0.01,0.03]#,0.05,0.1,0.15,0.25,0.70,0.95]
    	for per in percentages:
		i+=1
		if per>0:
			tindex = int(per*(t_tot-10))+10
		else:	
			tindex = 0
	    	col = ROYGBIV_map(float(i)/len(percentages),1.05,1)
		colS.append(col)

		FreeChargeDensity = Nf_tr[tindex] / (area*L_bin)
		LD_t = LD_tr[tindex]
		ND_free = FreeChargeDensity *4.*np.pi*Bjerrum*LD_t**2
		maxy = np.max([np.max(ND_free),maxy])
	    	ax1.errorbar(z_density[1:],ND_free[1:],yerr=None,color=col,ls='-',lw=1.10,label=str(per*100)+r'$\%$')

    	ax1.set_xlabel(r'$z $')
    	ax1.set_ylabel(r'$\rho / 2 q en^{\rm B}(t)$')
    	plt.setp(ax1.get_xticklabels(), fontsize=6.)
    	plt.setp(ax1.get_yticklabels(), fontsize=6.)

    	ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
#    	ax1.set_xlim(0,5)       
    	ax1.set_ylim(-maxy*1.1,maxy*1.1)
    
    	fig.set_size_inches(3.37,3.5)
    	graphname = 'Tr_free.png'
    	if filename in GClike:
    		graphname = 'GC_'+name_string+'_'+graphname
	elif filename in CSlike:
    		graphname = 'CS_'+name_string+'_'+graphname
	else:
    		graphname = 'temp_'+name_string+'_'+graphname
    	plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    	plt.close()





    	
    	print 'Plotting total charge density at key times'
    	fig=plt.figure()
    	ax1=fig.add_subplot(111)
	fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    	fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    	fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    	fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    	i=-1
    	j=-1
    	maxy = -1e9
    	colS=[]
    	percentages = [0.,0.001,0.005,0.01,0.03,0.05,0.1,0.15,0.25,0.70,0.95]
#    	percentages = [0.,0.001,0.005,0.01,0.03]#,0.05,0.1,0.15,0.25,0.70,0.95]
    	for per in percentages:
		i+=1
		if per>0:
			tindex = int(per*(t_tot-10))+10
		else:	
			tindex = 0
	    	col = ROYGBIV_map(float(i)/len(percentages),1.05,1)
		colS.append(col)

		TotChargeDensity = Nt_tr[tindex] / (area*L_bin)
		LD_t = LD_tr[tindex]
		ND_tot = np.array(TotChargeDensity) *4.*np.pi*Bjerrum*LD_t**2
		maxy = np.max([np.max(ND_tot),maxy])
	    	ax1.errorbar(z_density[1:],ND_tot[1:],yerr=None,color=col,ls='-',lw=1.10,label=str(per*100)+r'$\%$')

    	ax1.set_xlabel(r'$z $')
    	ax1.set_ylabel(r'$\rho_{\rm tot} / 2 q en^{\rm B}(t)$')
    	plt.setp(ax1.get_xticklabels(), fontsize=6.)
    	plt.setp(ax1.get_yticklabels(), fontsize=6.)

    	ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
#    	ax1.set_xlim(0,7)       
    	ax1.set_ylim(-0.5,maxy*1.1)
    
    	fig.set_size_inches(3.37,3.5)
    	graphname = 'Tr_tot.png'
    	if filename in GClike:
    		graphname = 'GC_'+name_string+'_'+graphname
	elif filename in CSlike:
    		graphname = 'CS_'+name_string+'_'+graphname
	else:
    		graphname = 'temp_'+name_string+'_'+graphname
    	plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    	plt.close()

    	print 'Plotting Q_EDL as a function of time'
    	fig=plt.figure()
    	ax1=fig.add_subplot(111)
    	fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    	fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    	fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    	fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    	ax1.errorbar(np.array(tcountS)/float(tauc),Q_tr,yerr=None,color='k',ls='-',lw=1.10)
#    	ax1.errorbar(np.array(tcountS)/float(t_tot),LD_tr,yerr=None,color='red',ls='-',lw=1.10)
    	ax1.set_xlabel(r'$t_{\rm MD} / \tau_{c} $')
    	ax1.set_ylabel(r'$Q_{\rm EDL}$')
    	plt.setp(ax1.get_xticklabels(), fontsize=6.)
    	plt.setp(ax1.get_yticklabels(), fontsize=6.)
#    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
    	ax1.set_xlim(0,1)    
        plt.ylim(ymin=-5)   
#    ax1.set_ylim(-1.75,1.75)
    	fig.set_size_inches(3.37,3.5)
    	graphname = 'Tr_Q.png'
    	if filename in GClike:
    		graphname = 'GC_'+name_string+'_'+graphname
	elif filename in CSlike:
    		graphname = 'CS_'+name_string+'_'+graphname
	else:
    		graphname = 'temp_'+name_string+'_'+graphname
    	plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    	plt.close()

    	print 'Plotting screening length as a function of time'
    	fig=plt.figure()
    	ax1=fig.add_subplot(111)
    	fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    	fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    	fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    	fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
#    	ax1.errorbar(np.array(tcountS)/float(t_tot),Q_tr,yerr=None,color='k',ls='-',lw=1.10)
    	ax1.errorbar(np.array(tcountS)/float(t_tot),LD_tr,yerr=None,color='red',ls='-',lw=1.10)
    	ax1.set_xlabel(r'$t_{\rm MD} / t^{\rm \infty} $')
    	ax1.set_ylabel(r'$\lambda_{\rm D}(\tilde t)$')
    	plt.setp(ax1.get_xticklabels(), fontsize=6.)
    	plt.setp(ax1.get_yticklabels(), fontsize=6.)
#    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
    	ax1.set_xlim(0,1)     
#    ax1.set_ylim(-1.75,1.75)
    	fig.set_size_inches(3.37,3.5)
    	graphname = 'Tr_lamD.png'
    	if filename in GClike:
    		graphname = 'GC_'+name_string+'_'+graphname
	elif filename in CSlike:
    		graphname = 'CS_'+name_string+'_'+graphname
	else:
    		graphname = 'temp_'+name_string+'_'+graphname
    	plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    	plt.close()



#    	print 'Plotting vz profiles at key times'
#    	fig=plt.figure()
#    	ax1=fig.add_subplot(111)
#	fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    	fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
#    	fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
#    	fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
#    	i=-1
#    	j=-1
#    	maxy = -1e9
#    	colS=[]
#    	percentages = [0.,0.001,0.005,0.01,0.03,0.05,0.1,0.15,0.25,0.70,0.95]
#    	for per in percentages:
#		i+=1
#		if per>0:
#			tindex = int(per*(t_tot-10))+10
#		else:	
#			tindex = 0
#	    	col = ROYGBIV_map(float(i)/len(percentages),1.05,1)
#		colS.append(col)

#		vzP_t = vzP_tr[tindex]
#		LD_t = LD_tr[tindex]
#		maxy = np.max([np.max(field_t),maxy])
#	    	ax1.errorbar(z_density[1:]/LD_t,vzP_t[1:],yerr=None,color=col,ls='-',lw=1.10,label=str(per*100)+r'$\%$')

#    	ax1.set_xlabel(r'$z/\lambda_{\rm D}(t)$')
#    	ax1.set_ylabel(r'$v_{z} $')
#    	plt.setp(ax1.get_xticklabels(), fontsize=6.)
#    	plt.setp(ax1.get_yticklabels(), fontsize=6.)

#    	ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.07,borderpad=0.25,labelspacing=0.1,ncol=1,handletextpad=0.05)
##    	ax1.set_xlim(0,5)       
##    	ax1.set_ylim(0.,maxy*1.1)
#    
#    	fig.set_size_inches(3.37,3.5)
#    	graphname = 'Tr_vz_plus.png'
#    	if filename in GClike:
#    		graphname = 'GC_'+str(N_tot)+'_'+graphname
#	elif filename in CSlike:
#    		graphname = 'CS_'+str(N_tot)+'_'+graphname
#    	plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    	plt.close()
    return  ##Above code for P2tr1()
#P2tr1()


def PlotMFS_MFSS():
    import matplotlib.pyplot as plt
#    from matplotlib.colors import LinearSegmentedColormap
#    from mpl_toolkits.mplot3d import Axes3D
    import sys

    print 'Plotting MFS & MFSS ratios'
    filenameS=[]
    for f in [[x for x in line.split()] for line in file(sys.argv[1],"r").readlines()]:
    	filenameS.append(f[0])

#    filenameS=['MFS_0.05_0.2.txt']

    distanceS = []
    plusS,minusS = [],[]
    freeS = []
    voltageS = []
    SigEffS = []
    S_ofzS = []
    dCapS = []
    solventS=[]
    PhipmS=[]
    PhisolS=[]
    colS=[]
    ls_xS=[]
    labelS=[]

    i=-1
    lo,mi,hi=-1,-1,-1
    PhiLabel = [r'$0.01$',r'$0.07$',r'$0.10$',r'$0.20$',r'$0.30$',r'$0.40$',r'$0.50$']

    for filename in filenameS:
    	i+=1
    	if i==len(PhiLabel):
    		i=0

    	print filename
#       	MMA=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()])

       	MMA=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()])

       	distanceS.append(MMA[:,0])
       	plusS.append(MMA[:,1])
       	minusS.append(MMA[:,2])
       	freeS.append(MMA[:,3])
       	voltageS.append(MMA[:,4])
       	SigEffS.append(MMA[:,5])
       	S_ofzS.append(MMA[:,6])
       	dCapS.append(MMA[:,7])
 

	if filename[:2] == 'CS':
		ls_xS.append('--')
		#Extract relevant parameters from filename
		j=0
		value=''
		for letter in filename[8:]:
	    		if letter!='_':
	    			value+=letter
			else:
				j+=1
				if j==1:
					Phiion = float(value)
					PhipmS.append(Phiion)
					value=''
		Phisolvent = float(value[:-4])
		PhisolS.append(Phisolvent)
		colS.append(ROYGBIV_map(float(0),0.75,1))	
		labelS.append(str(Phiion))
	elif filename[:4] == 'MFS_':
		ls_xS.append('-')
	       	solventS.append(MMA[:,8])
	       	j=0
		value=''
		for letter in filename[4:]:
	    		if letter!='_':
	    			value+=letter
			else:
				j+=1
				if j==1:
					Phiion = float(value)
					PhipmS.append(Phiion)
					value=''
		Phisolvent = float(value[:-4])
		PhisolS.append(Phisolvent)
		colS.append(ROYGBIV_map(float(Phisolvent),0.75,1))	
		labelS.append(str(Phiion)+','+str(Phisolvent))
	elif filename[:4] == 'MFSS':
		ls_xS.append('-.')
	       	solventS.append(MMA[:,8])
	       	j=0
		value=''
		for letter in filename[5:]:
	    		if letter!='_':
	    			value+=letter
			else:
				j+=1
				if j==1:
					Phiion = float(value)
					PhipmS.append(Phiion)
					value=''
				if j==2:
					Phisolvent = float(value)
					PhisolS.append(Phisolvent)
					value=''
		w = float(value[1:-4])				
		colS.append(ROYGBIV_map(float(Phisolvent),0.75,1))	
		labelS.append(str(Phiion)+','+str(Phisolvent))
					
#	###rhoCS versus rhoMFS
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
#    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
#    i=-1
#    i=-1
#    for (free,S,filename,col,ls_x,label_x) in zip(freeS,S_ofzS,filenameS,colS,ls_xS,labelS):  
#	i+=1
#   	ax1.errorbar(freeS[0],free,yerr=None,ls=ls_x,lw=1.0,color=col,marker='None',label=label_x)	
#    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
#    ax1.set_ylabel(r'$\tilde \rho_{\rm MFS}$',size='small')
#    ax1.set_xlabel(r'$\tilde \rho_{\rm CS}(\Phi_{\pm}=0.15)$',size='small')
#    plt.setp(ax1.get_xticklabels(), fontsize='small')
#    plt.setp(ax1.get_yticklabels(), fontsize='small')
##    ax1.xaxis.set_major_locator(MaxNLocator(6))
##    ax1.yaxis.set_major_locator(MaxNLocator(6))
##    ax1.set_xlim(0,1.5)     
#    plt.xlim(xmax=0)
##    ax1.set_ylim(0,1.5)
#    fig.set_size_inches(3.37,3.5)
#    graphname = 'rho_vs_rho.png'
#    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    plt.close() 



    print "rhoS_vs_distance"
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    i=-1
#    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='--',lw=1.0,color='k',marker='None',label=r'$\Phi_{\pm}=0.18$')
#    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='-',lw=1.0,color='k',marker='None',label=r'$(\Phi_{\pm},\Phi_{\rm sol})=(0.18,0.50)$')
    for (distance,plus,filename,col,ls_x,label_x,minus,voltage,free) in zip(distanceS,plusS,filenameS,colS,ls_xS,labelS,minusS,voltageS,freeS):  
	i+=1
	if i==0:	
#	   	ax1.errorbar(distance,voltage,yerr=None,ls=ls_x,lw=1.0,color='blue',marker='None')#	,label=r'$Voltage$')	
   		ax1.errorbar(distance,plus,yerr=None,ls=ls_x,lw=1.0,color='green',marker='None',label=r'${\rm Co-ion}$')	
   		ax1.errorbar(distance,minus,yerr=None,ls=ls_x,lw=1.0,color='red',marker='None',label=r'${\rm Counter-ion}$')	
   		ax1.errorbar(distance,solventS[0],yerr=None,ls=ls_x,lw=1.0,color='purple',marker='None',label=r'${\rm Solvent}$')
   		   		
#   		ax1.errorbar(distance,-np.array(free),yerr=None,ls=ls_x,lw=1.0,color='grey',marker='None')
	else:
#	   	ax1.errorbar(distance,voltage,yerr=None,ls=ls_x,lw=1.0,color='blue',marker='None',label=r'${\rm Voltage}$')	
   		ax1.errorbar(distance,plus,yerr=None,ls=ls_x,lw=1.0,color='green',marker='None')#,label=r'${\rm Co-ion}$')	
   		ax1.errorbar(distance,minus,yerr=None,ls=ls_x,lw=1.0,color='red',marker='None')#,label=r'${\rm Counter-ion}$')	
#   		ax1.errorbar(distance,-np.array(free),yerr=None,ls=ls_x,lw=1.0,color='grey',marker='None')#,label=r'${\rm |Free charge|}$')	
   		ax1.errorbar(distance,solventS[0],yerr=None,ls=ls_x,lw=1.5,color='purple',marker='None')#,label=r'${\rm Solvent}$')	
    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
    ax1.set_ylabel(r'$\tilde \rho_{\rm i}$',size='small')
    ax1.set_xlabel(r'$\tilde z$',size='small')
    plt.setp(ax1.get_xticklabels(), fontsize='small')
    plt.setp(ax1.get_yticklabels(), fontsize='small')
#    ax1.xaxis.set_major_locator(MaxNLocator(6))
#    ax1.yaxis.set_major_locator(MaxNLocator(6))
#    ax1.set_xlim(0,1.5)     
    plt.xlim(xmin=0)
    plt.ylim(ymin=-0.1)
#    ax1.set_ylim(0,1.5)
    fig.set_size_inches(3.37,3.5)
    graphname = 'rhoS_vs_distance.png'
    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    plt.close() 

    print 'freeVolt_vs_distance'
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    i=-1
#    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='--',lw=1.0,color='k',marker='None',label=r'$\Phi_{\pm}=0.18$')
#    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='-',lw=1.0,color='k',marker='None',label=r'$(\Phi_{\pm},\Phi_{\rm sol})=(0.18,0.50)$')
    for (distance,plus,filename,col,ls_x,label_x,minus,voltage,free) in zip(distanceS,plusS,filenameS,colS,ls_xS,labelS,minusS,voltageS,freeS):  
	i+=1
	if i==1:	
	   	ax1.errorbar(distance,voltage,yerr=None,ls=ls_x,lw=1.0,color='blue',marker='None')#	,label=r'$Voltage$')	
#   		ax1.errorbar(distance,plus,yerr=None,ls=ls_x,lw=1.0,color='green',marker='None')
#   		ax1.errorbar(distance,minus,yerr=None,ls=ls_x,lw=1.0,color='red',marker='None')
   		ax1.errorbar(distance,-np.array(free),yerr=None,ls=ls_x,lw=1.0,color='grey',marker='None')
	else:
	   	ax1.errorbar(distance,voltage,yerr=None,ls=ls_x,lw=1.0,color='blue',marker='None',label=r'${\rm Voltage}$')	
#   		ax1.errorbar(distance,plus,yerr=None,ls=ls_x,lw=1.0,color='green',marker='None',label=r'${\rm Co-ion}$')	
#   		ax1.errorbar(distance,minus,yerr=None,ls=ls_x,lw=1.0,color='red',marker='None',label=r'${\rm Counter-ion}$')	
   		ax1.errorbar(distance,-np.array(free),yerr=None,ls=ls_x,lw=1.0,color='grey',marker='None',label=r'${\rm |Free Charge|}$')	
#   		ax1.errorbar(distance,solventS[0],yerr=None,ls=ls_x,lw=1.0,color='purple',marker='None',label=r'${\rm Solvent}$')	
    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
#    ax1.set_ylabel(r'$\tilde \psi,\tilde \rho_{\rm \pm}, \tilde \rho_{\rm sol}$',size='small')
    ax1.set_ylabel(r'$\tilde \psi,-\tilde \rho_{\rm free}$',size='small')
    ax1.set_xlabel(r'$\tilde z$',size='small')
    plt.setp(ax1.get_xticklabels(), fontsize='small')
    plt.setp(ax1.get_yticklabels(), fontsize='small')
#    ax1.xaxis.set_major_locator(MaxNLocator(6))
#    ax1.yaxis.set_major_locator(MaxNLocator(6))
#    ax1.set_xlim(0,1.5)     
    plt.xlim(xmin=0)
    plt.ylim(ymin=-0.1)
#    ax1.set_ylim(0,1.5)
    fig.set_size_inches(3.37,3.5)
    graphname = 'freeVolt_vs_distance.png'
    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    plt.close() 



#	###rho vs S
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
#    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
#    i=-1
#    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='--',lw=1.0,color='white',marker='None',label=r'$\Phi_{\pm},\Phi_{\rm sol}$')
#    for (free,S,filename,col,ls_x,label_x) in zip(freeS,S_ofzS,filenameS,colS,ls_xS,labelS):  
#	i+=1
#   	ax1.errorbar(S,-np.array(free),yerr=None,ls=ls_x,lw=1.0,color=col,marker='None',label=label_x)	
#    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
#    ax1.set_ylabel(r'$-\tilde \rho$',size='small')
#    ax1.set_xlabel(r'$S_{\rm MFS}$',size='small')
#    plt.setp(ax1.get_xticklabels(), fontsize='small')
#    plt.setp(ax1.get_yticklabels(), fontsize='small')
##    ax1.xaxis.set_major_locator(MaxNLocator(6))
##    ax1.yaxis.set_major_locator(MaxNLocator(6))
##    ax1.set_xlim(0,1.5)     
#    plt.xlim(xmin=0)
#    plt.ylim(ymin=0)
##    ax1.set_ylim(0,1.5)
#    fig.set_size_inches(3.37,3.5)
#    graphname = 'rho_vs_S.png'
#    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    plt.close() 

#	###dCap versus S
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
#    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
#    i=-1
##    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='--',lw=1.0,color='white',marker='None',label=r'$\Phi_{\pm}$')
#    i=-1
#    for (dC,S,filename,col,ls_x,label_x) in zip(dCapS,S_ofzS,filenameS,colS,ls_xS,labelS):  
#	i+=1
##	if i==0:
##		label_x='0,1'
##	if i==1:
##		label_x='0.15,0'		
##	if i==9:
##	    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='-',lw=1.0,color='white',marker='None',label=r'$\Phi_{\pm},\Phi_{\rm sol}$')
#   	ax1.errorbar(S,dC,yerr=None,ls=ls_x,lw=1.0,color=col,marker='None',label='0,1')	
#   	
#    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
#    ax1.set_xlabel(r'$S$',size='small')
#    ax1.set_ylabel(r'$C\prime (\lambda_{\rm D}/\epsilon)$',size='small')
#    plt.setp(ax1.get_xticklabels(), fontsize='small')
#    plt.setp(ax1.get_yticklabels(), fontsize='small')
##    ax1.xaxis.set_major_locator(MaxNLocator(6))
##    ax1.yaxis.set_major_locator(MaxNLocator(6))
##    ax1.set_,xlim(0,6)     
#    plt.xlim(xmin=0)
#    ax1.set_ylim(0,1.5)
#    fig.set_size_inches(3.37,3.5)
#    graphname = 'dCap_vs_S.png'
#    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    plt.close() 


	###dCap versus voltage
    print 'dCap versus voltage'
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    i=-1
#    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='--',lw=1.0,color='white',marker='None',label=r'$\Phi_{\pm},\Phi_{\rm sol}$')
    ax1.errorbar(voltage,np.cosh(np.array(voltage)/2),yerr=None,ls='-',lw=1.0,color='k',marker='None',label=r'${\rm GC}$')
    i=-1
    for (dC,voltage,filename,col,ls_x,label_x) in zip(dCapS,voltageS,filenameS,colS,ls_xS,labelS):  
	i+=1

	if i==0:
		label_x=r'${\rm 2}$'+' '+r'${\rm Species}:$'+' '+r'HS$'
	if i==1:
		label_x=r'${\rm 3}$'+' '+r'${\rm Species}:$'+' '+r'HS$'
	if i==2:
		label_x=r'${\rm 3}$'+' '+r'${\rm Species}:$'+' '+r'Sticky Solvent$'		

	print filename,i,label_x,ls_x
   	ax1.errorbar(voltage,dC,yerr=None,ls=ls_x,lw=1.0,color=col,marker='None',label=label_x)	
   	
    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
    ax1.set_xlabel(r'$\tilde \psi$',size='small')
    ax1.set_ylabel(r'$C^{\prime} (\lambda_{\rm D}/\epsilon)$',size='small')
    plt.setp(ax1.get_xticklabels(), fontsize='small')
    plt.setp(ax1.get_yticklabels(), fontsize='small')
#    ax1.xaxis.set_major_locator(MaxNLocator(6))
#    ax1.yaxis.set_major_locator(MaxNLocator(6))
#    ax1.set_,xlim(0,6)     
    plt.xlim(xmin=0)
    ax1.set_ylim(0,1.5)
    fig.set_size_inches(3.37,3.5)
    graphname = 'dCap_vs_voltage.png'
    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    plt.close() 




    return
#PlotMFS_MFSS()


def PlotMFS_MFSS_cap():
    import matplotlib.pyplot as plt
#    from matplotlib.colors import LinearSegmentedColormap
#    from mpl_toolkits.mplot3d import Axes3D
    import sys

    print 'Plotting MFS & MFSS ratios'
    filenameS=[]
    for f in [[x for x in line.split()] for line in file(sys.argv[1],"r").readlines()]:
    	filenameS.append(f[0])

#    filenameS=['MFS_0.05_0.2.txt']

    distanceS = []
    plusS,minusS = [],[]
    freeS = []
    voltageS = []
    

    i=-1
    lo,mi,hi=-1,-1,-1
    PhiLabel = [r'$0.01$',r'$0.07$',r'$0.10$',r'$0.20$',r'$0.30$',r'$0.40$',r'$0.50$']

    for filename in filenameS:
    	i+=1
    	if i==len(PhiLabel):
    		i=0

    	print filename
#       	MMA=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()])

       	MMA=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()])

       	voltage = np.array(MMA[:,0])
       	CS =  MMA[:,1]
       	MFS = MMA[:,2]
       	MFSS= MMA[:,3]

       	if i==0:
       		num = '0.0'
	elif i==1:
		num = '-0.05'
	elif i==2:
		num = '-0.5'
	elif i==3:
		num = '-1.0'
	elif i==4:
		num = '-5.0'								


	###dCap versus voltage
    print 'dCap versus voltage'
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    i=-1
#    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='--',lw=1.0,color='white',marker='None',label=r'$\Phi_{\pm},\Phi_{\rm sol}$')
    ax1.errorbar(voltage,np.cosh(np.array(voltage)/2),yerr=None,ls='-',lw=1.0,color='k',marker='None',label=r'${\rm GC}$')
    ax1.errorbar(voltage,CS,yerr=None,ls='--',lw=1.0,color='red',marker='None',label=r'${\rm CS}$')
    ax1.errorbar(voltage,MFS,yerr=None,ls='-',lw=1.0,color='blue',marker='None',label=r'${\rm MFS}$')  
    ax1.errorbar(voltage,MFSS,yerr=None,ls='-.',lw=1.0,color='blue',marker='None',label=r'${\rm MFSS}$') 
    
    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
    ax1.set_xlabel(r'$\tilde \psi$',size='small')
    ax1.set_ylabel(r'$C^{\prime} (\lambda_{\rm D}/\epsilon)$',size='small')
    plt.setp(ax1.get_xticklabels(), fontsize='small')
    plt.setp(ax1.get_yticklabels(), fontsize='small')
#    ax1.xaxis.set_major_locator(MaxNLocator(6))
#    ax1.yaxis.set_major_locator(MaxNLocator(6))
#    ax1.set_,xlim(0,6)     
    plt.xlim(xmin=0)
    ax1.set_ylim(0,1.5)
    fig.set_size_inches(3.37,3.5)
    graphname = 'dCap_vs_voltage_all.png'
    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    plt.close() 


    return
#PlotMFS_MFSS_cap()

def PlotMFS_MFSS_cap2():
    import matplotlib.pyplot as plt
#    from matplotlib.colors import LinearSegmentedColormap
#    from mpl_toolkits.mplot3d import Axes3D
    import sys

    print 'Plotting MFS & MFSS ratios'
#    filenameS=[]
#    for f in [[x for x in line.split()] for line in file(sys.argv[1],"r").readlines()]:
#    	filenameS.append(f[0])

#    filenameS=['MFS_0.05_0.2.txt']

#    distanceS = []
#    plusS,minusS = [],[]
#    freeS = []
#    voltageS = []
#    

#    i=-1
#    lo,mi,hi=-1,-1,-1
#    PhiLabel = [r'$0.01$',r'$0.07$',r'$0.10$',r'$0.20$',r'$0.30$',r'$0.40$',r'$0.50$']

#    for filename in filenameS:
#    	i+=1
#    	if i==len(PhiLabel):
#    		i=0

#    	print filename
##       	MMA=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()])

#       	MMA=np.array([[float(x) for x in line.split()] for line in file(filename,"r").readlines()])

#       	voltage = np.array(MMA[:,0])
#       	CS =  MMA[:,1]
#       	MFS = MMA[:,2]
#       	MFSS= MMA[:,3]

#       	if i==0:
#       		num = '0.0'
#	elif i==1:
#		num = '-0.05'
#	elif i==2:
#		num = '-0.5'
#	elif i==3:
#		num = '-1.0'
#	elif i==4:
#		num = '-5.0'								


	###dCap versus voltage
    print 'dCap versus voltage'
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    i=-1
#    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='--',lw=1.0,color='white',marker='None',label=r'$\Phi_{\pm},\Phi_{\rm sol}$')

    MMA=np.array([[float(x) for x in line.split()] for line in file('dCap_0.2_0.5_w_0.0.txt',"r").readlines()])
    voltage = np.array(MMA[:,0])
    CS =  MMA[:,1]
    MFS = MMA[:,2]
    MFSS= MMA[:,3]
    ax1.errorbar(voltage,np.cosh(np.array(voltage)/2),yerr=None,ls='-',lw=1.0,color='k',marker='None',label=r'${\rm GC}$')
    ax1.errorbar(voltage,np.sinh(voltage) /((1. + 2*0.20*np.sinh(voltage/2)**2)*np.sqrt((2./0.20)*np.log(1+2.*0.20*np.sinh(voltage/2)**2))),yerr=None,ls='-.',lw=1.0,color='k',marker='None',label=r'${\rm Bik}$')
    ax1.errorbar(voltage,CS,yerr=None,ls='--',lw=1.0,color='red',marker='None',label=r'${\rm CS}$')

    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='--',lw=1.0,color='white',marker='None',label=r'${\rm MFSS}$')
    ax1.errorbar(voltage,MFSS,yerr=None,ls='-',lw=0.1,color='blue',marker='None',label=r'$w=0$')  


    MMA=np.array([[float(x) for x in line.split()] for line in file('dCap_0.2_0.5_w_-0.05.txt',"r").readlines()])
    voltage = np.array(MMA[:,0])
    CS =  MMA[:,1]
    MFS = MMA[:,2]
    MFSS= MMA[:,3]
#    ax1.errorbar(voltage,np.cosh(np.array(voltage)/2),yerr=None,ls='-',lw=1.0,color='k',marker='None',label=r'${\rm GC}$')
#    ax1.errorbar(voltage,CS,yerr=None,ls='--',lw=1.0,color='red',marker='None',label=r'${\rm CS}$')
    ax1.errorbar(voltage,MFSS,yerr=None,ls='-',lw=0.15,color='blue',marker='None',label=r'$w=-0.05$')  

    MMA=np.array([[float(x) for x in line.split()] for line in file('dCap_0.2_0.5_w_-0.5.txt',"r").readlines()])
    voltage = np.array(MMA[:,0])
    CS =  MMA[:,1]
    MFS = MMA[:,2]
    MFSS= MMA[:,3]
#    ax1.errorbar(voltage,np.cosh(np.array(voltage)/2),yerr=None,ls='-',lw=1.0,color='k',marker='None',label=r'${\rm GC}$')
#    ax1.errorbar(voltage,CS,yerr=None,ls='--',lw=1.0,color='red',marker='None',label=r'${\rm CS}$')
    ax1.errorbar(voltage,MFSS,yerr=None,ls='-',lw=0.5,color='blue',marker='None',label=r'$w=-0.5$') 


    MMA=np.array([[float(x) for x in line.split()] for line in file('dCap_0.2_0.5_w_-1.0.txt',"r").readlines()])
    voltage = np.array(MMA[:,0])
    CS =  MMA[:,1]
    MFS = MMA[:,2]
    MFSS= MMA[:,3]
#    ax1.errorbar(voltage,np.cosh(np.array(voltage)/2),yerr=None,ls='-',lw=1.0,color='k',marker='None',label=r'${\rm GC}$')
#    ax1.errorbar(voltage,CS,yerr=None,ls='--',lw=1.0,color='red',marker='None',label=r'${\rm CS}$')
    ax1.errorbar(voltage,MFSS,yerr=None,ls='-',lw=1.0,color='blue',marker='None',label=r'$w=-1.0$') 

    MMA=np.array([[float(x) for x in line.split()] for line in file('dCap_0.2_0.5_w_-5.0.txt',"r").readlines()])
    voltage = np.array(MMA[:,0])
    CS =  MMA[:,1]
    MFS = MMA[:,2]
    MFSS= MMA[:,3]
#    ax1.errorbar(voltage,np.cosh(np.array(voltage)/2),yerr=None,ls='-',lw=1.0,color='k',marker='None',label=r'${\rm GC}$')
#    ax1.errorbar(voltage,CS,yerr=None,ls='--',lw=1.0,color='red',marker='None',label=r'${\rm CS}$')
    ax1.errorbar(voltage,MFSS,yerr=None,ls='-',lw=1.5,color='blue',marker='None',label=r'$w=-5.0$') 
        
    
    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
    ax1.set_xlabel(r'$\tilde \psi$',size='small')
    ax1.set_ylabel(r'$C^{\prime} (\lambda_{\rm D}/\epsilon)$',size='small')
    plt.setp(ax1.get_xticklabels(), fontsize='small')
    plt.setp(ax1.get_yticklabels(), fontsize='small')
#    ax1.xaxis.set_major_locator(MaxNLocator(6))
#    ax1.yaxis.set_major_locator(MaxNLocator(6))
#    ax1.set_,xlim(0,6)     
    plt.xlim(xmin=0)
    ax1.set_ylim(0,1.5)
    fig.set_size_inches(3.37,3.5)
    graphname = 'dCap_vs_voltage_all2.png'
    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    plt.close() 


    return
#PlotMFS_MFSS_cap2()

def PlotMFS_MFSS_S():
    import matplotlib.pyplot as plt
#    from matplotlib.colors import LinearSegmentedColormap
#    from mpl_toolkits.mplot3d import Axes3D
    import sys

	###dCap versus voltage
    print 'rhoFREE versus S_MFSS'
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
    i=-1
    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='-',lw=1.0,color='white',marker='None',label=r'$(\Phi_{\rm sol},w)=(0.50,-0.5)$')

    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.15_0.5_w_-0.5.txt',"r").readlines()])
    free = -np.array(MMA[:,3])
    S =  MMA[:,6]
    ax1.errorbar(S,free,yerr=None,ls='-',lw=0.15,color='red',marker='None',label=r'$\Phi_{\pm}=0.15$')

#    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='--',lw=1.0,color='white',marker='None',label=r'${\rm MFSS}$')
#    ax1.errorbar(voltage,MFSS,yerr=None,ls='-',lw=0.1,color='blue',marker='None',label=r'$w=0$')  


    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.18_0.5_w_-0.5.txt',"r").readlines()])
    free = -np.array(MMA[:,3])
    S =  MMA[:,6]
    ax1.errorbar(S,free,yerr=None,ls='-',lw=0.2,color='red',marker='None',label=r'$\Phi_{\pm}=0.18$')


    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.2_0.5_w_-0.5.txt',"r").readlines()])
    free = -np.array(MMA[:,3])
    S =  MMA[:,6]
    ax1.errorbar(S,free,yerr=None,ls='-',lw=0.5,color='red',marker='None',label=r'$\Phi_{\pm}=0.20$')

    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.25_0.5_w_-0.5.txt',"r").readlines()])
    free = -np.array(MMA[:,3])
    S =  MMA[:,6]
    ax1.errorbar(S,free,yerr=None,ls='-',lw=1.0,color='red',marker='None',label=r'$\Phi_{\pm}=0.25$')

    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='--',lw=1.0,color='white',marker='None',label=r'$(\Phi_{\rm sol},\Phi_{\pm})=(0.50,0.20)$')

    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.2_0.5_w_0.txt',"r").readlines()])
    free = -np.array(MMA[:,3])
    S =  MMA[:,6]
    ax1.errorbar(S,free,yerr=None,ls='-.',lw=0.20,color='k',marker='None',label=r'$w = 0$')

    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.2_0.5_w_-0.05.txt',"r").readlines()])
    free = -np.array(MMA[:,3])
    S =  MMA[:,6]
    ax1.errorbar(S,free,yerr=None,ls='--',lw=0.25,color='k',marker='None',label=r'$w = -0.05$')

    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.2_0.5_w_-0.5.txt',"r").readlines()])
    free = -np.array(MMA[:,3])
    S =  MMA[:,6]
    ax1.errorbar(S,free,yerr=None,ls='--',lw=0.5,color='k',marker='None',label=r'$w = -0.50$')

    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.2_0.5_w_-1.txt',"r").readlines()])
    free = -np.array(MMA[:,3])
    S =  MMA[:,6]
    ax1.errorbar(S,free,yerr=None,ls='--',lw=1.0,color='k',marker='None',label=r'$w = -1.0$')   

    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.2_0.5_w_-5.txt',"r").readlines()])
    free = -np.array(MMA[:,3])
    S =  MMA[:,6]
    ax1.errorbar(S,free,yerr=None,ls='--',lw=1.5,color='k',marker='None',label=r'$w = -5.0$')       
        
    
    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
    ax1.set_xlabel(r'$S_{\rm MFSS}$',size='small')
    ax1.set_ylabel(r'$-\rho/2qen^{\rm B}$',size='small')
    plt.setp(ax1.get_xticklabels(), fontsize='small')
    plt.setp(ax1.get_yticklabels(), fontsize='small')
#    ax1.xaxis.set_major_locator(MaxNLocator(6))
#    ax1.yaxis.set_major_locator(MaxNLocator(6))
    ax1.set_xlim(0,10)     
    ax1.set_ylim(0,4)
    fig.set_size_inches(3.37,3.5)
    graphname = 'rhoF_vs_S_all2.png'
    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
    plt.close() 


    return
#PlotMFS_MFSS_S()

def LCorDemo():
    import matplotlib.pyplot as plt
#    from matplotlib.colors import LinearSegmentedColormap
#    from mpl_toolkits.mplot3d import Axes3D
    import sys




    print 'NEW: Lcor Demo'
    fig=plt.figure()
    ax1=fig.add_subplot(111)

    
    volt_theory = np.linspace(0,16)
    GC_cap = 2*np.sinh(volt_theory/2.)
    ax1.errorbar(volt_theory,GC_cap,yerr=None,ls=':',lw=1.0,color='k',marker='None',label =r'${\rm GC}$')
    
    QV_Bik=np.array([[float(x) for x in line.split()] for line in file('Bik_0.05_16.txt',"r").readlines()])
    ax1.errorbar(QV_Bik[:,4],QV_Bik[:,5],yerr=None,ls='-.',lw=1.0,color='grey',marker='None')#,label =r'${\rm Bik}$')

    QV_Bik=np.array([[float(x) for x in line.split()] for line in file('Bik_0.10_16.txt',"r").readlines()])
    ax1.errorbar(QV_Bik[:,4],QV_Bik[:,5],yerr=None,ls='-.',lw=1.0,color='violet',marker='None')#,label =r'${\rm Bik}$')

    ax1.plot([-1,-1],[-1,-1],color = 'k', ls = '-.',lw = 1.0,marker = 'None',label =r'${\rm Bik}$')               
    ax1.plot([-1,-1],[-1,-1],color = 'k', ls = '-',lw = 1.0,marker = 'None',label =r'${\rm CS}$')
    ax1.plot([-1,-1],[-1,-1],color = 'k', ls = 'None',lw = 1.0,marker = 'None',label =r'$ $')
    for (i,col) in enumerate(groupColorS):
    	plotzetaS = zetaS[i*10:10*(i+1)]
    	plotNDSigmaS = NDSigmaS[i*10:10*(i+1)]
    	if len(filenameS)==88:
	    	plotzetaS = zetaS[i*11:11*(i+1)]
	    	plotNDSigmaS = NDSigmaS[i*11:11*(i+1)]		    		

	QV_CS=np.array([[float(x) for x in line.split()] for line in file('CS'+TheoryNames[i],"r").readlines()])
	ax1.errorbar(QV_CS[:,4],QV_CS[:,5],yerr=None,ls='-',lw=1.0,color=col,marker='None')#,label =groupLabelS[i])

		#The following adds Bikerman theory
#		QV_CS=np.array([[float(x) for x in line.split()] for line in file('Bik'+TheoryNames[i],"r").readlines()])
#   		ax1.errorbar(QV_CS[:,4],QV_CS[:,5],yerr=None,ls='-.',lw=1.0,color=col,marker='None')#,label =groupLabelS[i])

	mark = markerS[i]
    	ax1.plot(plotzetaS,plotNDSigmaS,color=col,ls='',marker = mark)#,label =groupLabelS[i])
    	ax1.plot([-1,-1],[-1,-1],color = col, ls = 'None',lw = 0.9,marker = mark,label =groupLabelS[i])

    ax1.errorbar(volt_theory,GC_cap,yerr=None,ls=':',lw=1.0,color='k',marker='None')#,label =r'${\rm GC}$')
#	    ax1.set_ylabel(r'$-\tilde \Sigma = ( \Sigma_{\rm cor} +  \Sigma_{\rm dif} ) / \Sigma_{\rm ref}$',fontsize=10.)
    ax1.set_ylabel(r'$\Sigma / \Sigma_{\rm ref}$',fontsize=10.)
    ax1.set_xlabel(r'$(\phi_0 - \phi_{\rm B}) / \phi_{\rm T}$',fontsize=10.)
    plt.setp(ax1.get_xticklabels(), fontsize=8.)
    plt.setp(ax1.get_yticklabels(), fontsize=8.)
    ax1.legend(loc='best',numpoints=1,prop=dict(size=8.),columnspacing=0.11,borderpad=0.15,labelspacing=0.1,ncol=3,handletextpad=0.07,markerscale=0.75)     
    ax1.set_xlim(0,16.5) 
#		    ax1.set_ylim(0,12.0)
    ax1.set_ylim(0,8.0)
#		    fig.set_size_inches(3.37,3.0)
    fig.set_size_inches(3.37,2.5)
    plt.savefig('lCor_demo.png', dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False,bbox_inches='tight')   
#	##Biksucks
    plt.close()




#	###dCap versus voltage
#    print 'rhoFREE versus S_MFSS'
#    fig=plt.figure()
#    ax1=fig.add_subplot(111)
#    fig.subplots_adjust(right=0.95) ##Lower puts more whitespace on the right
#    fig.subplots_adjust(bottom=0.13)##Larger adds whitespace
#    fig.subplots_adjust(left=0.15) #Larger puts more whitespace on the left
#    fig.subplots_adjust(top=0.97) ##Smaller adds whitespace to top
#    i=-1
#    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='-',lw=1.0,color='white',marker='None',label=r'$(\Phi_{\rm sol},w)=(0.50,-0.5)$')

#    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.15_0.5_w_-0.5.txt',"r").readlines()])
#    free = -np.array(MMA[:,3])
#    S =  MMA[:,6]
#    ax1.errorbar(S,free,yerr=None,ls='-',lw=0.15,color='red',marker='None',label=r'$\Phi_{\pm}=0.15$')

##    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='--',lw=1.0,color='white',marker='None',label=r'${\rm MFSS}$')
##    ax1.errorbar(voltage,MFSS,yerr=None,ls='-',lw=0.1,color='blue',marker='None',label=r'$w=0$')  


#    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.18_0.5_w_-0.5.txt',"r").readlines()])
#    free = -np.array(MMA[:,3])
#    S =  MMA[:,6]
#    ax1.errorbar(S,free,yerr=None,ls='-',lw=0.2,color='red',marker='None',label=r'$\Phi_{\pm}=0.18$')


#    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.2_0.5_w_-0.5.txt',"r").readlines()])
#    free = -np.array(MMA[:,3])
#    S =  MMA[:,6]
#    ax1.errorbar(S,free,yerr=None,ls='-',lw=0.5,color='red',marker='None',label=r'$\Phi_{\pm}=0.20$')

#    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.25_0.5_w_-0.5.txt',"r").readlines()])
#    free = -np.array(MMA[:,3])
#    S =  MMA[:,6]
#    ax1.errorbar(S,free,yerr=None,ls='-',lw=1.0,color='red',marker='None',label=r'$\Phi_{\pm}=0.25$')

#    ax1.errorbar([-1,-1],[-1,-1],yerr=None,ls='--',lw=1.0,color='white',marker='None',label=r'$(\Phi_{\rm sol},\Phi_{\pm})=(0.50,0.20)$')

#    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.2_0.5_w_0.txt',"r").readlines()])
#    free = -np.array(MMA[:,3])
#    S =  MMA[:,6]
#    ax1.errorbar(S,free,yerr=None,ls='-.',lw=0.20,color='k',marker='None',label=r'$w = 0$')

#    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.2_0.5_w_-0.05.txt',"r").readlines()])
#    free = -np.array(MMA[:,3])
#    S =  MMA[:,6]
#    ax1.errorbar(S,free,yerr=None,ls='--',lw=0.25,color='k',marker='None',label=r'$w = -0.05$')

#    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.2_0.5_w_-0.5.txt',"r").readlines()])
#    free = -np.array(MMA[:,3])
#    S =  MMA[:,6]
#    ax1.errorbar(S,free,yerr=None,ls='--',lw=0.5,color='k',marker='None',label=r'$w = -0.50$')

#    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.2_0.5_w_-1.txt',"r").readlines()])
#    free = -np.array(MMA[:,3])
#    S =  MMA[:,6]
#    ax1.errorbar(S,free,yerr=None,ls='--',lw=1.0,color='k',marker='None',label=r'$w = -1.0$')   

#    MMA=np.array([[float(x) for x in line.split()] for line in file('MFSS_0.2_0.5_w_-5.txt',"r").readlines()])
#    free = -np.array(MMA[:,3])
#    S =  MMA[:,6]
#    ax1.errorbar(S,free,yerr=None,ls='--',lw=1.5,color='k',marker='None',label=r'$w = -5.0$')       
#        
#    
#    ax1.legend(loc='best',numpoints=1,prop=dict(size=6.),columnspacing=0.12,borderpad=0.3,labelspacing=0.1,ncol=1,handletextpad=0.1) 
#    ax1.set_xlabel(r'$S_{\rm MFSS}$',size='small')
#    ax1.set_ylabel(r'$-\rho/2qen^{\rm B}$',size='small')
#    plt.setp(ax1.get_xticklabels(), fontsize='small')
#    plt.setp(ax1.get_yticklabels(), fontsize='small')
##    ax1.xaxis.set_major_locator(MaxNLocator(6))
##    ax1.yaxis.set_major_locator(MaxNLocator(6))
#    ax1.set_xlim(0,10)     
#    ax1.set_ylim(0,4)
#    fig.set_size_inches(3.37,2.0)
#    graphname = 'LCor_Demo.png'
#    plt.savefig(graphname, dpi=1600, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format=None, transparent=False)   
#    plt.close() 


    return
#LCorDemo()
P2F1()



    
#import random
#Ntot = 1004
#Lx = 23.
#Ly = Lx
#Lz = 502.
#Numbins = 502
#bin=np.array([((-Lz/2.)+i*(Lz/float(Numbins))) for i in range(Numbins)]+[Lz*0.5])
#bin+=1
##print bin
#bin=bin[0:len(bin)/2] #But we only need bins from one wall to the MIDDLE of the box
#count = 0
#zp =[]
#for (zlo,zhi) in zip(bin[0:len(bin)-1],bin[1:len(bin)]):
#		x=random.uniform(-0.5*Lx,0.5*Lx)
#		y=random.uniform(-0.5*Ly,0.5*Ly)
#		z=random.uniform(zlo,zhi)
#		print "create_atoms	%i single %1.2f %1.2f %1.2f units box" % (2,x,y,z)
#		print "create_atoms	%i single %1.2f %1.2f %1.2f units box" % (1,x,y,-z)
#		x=random.uniform(-0.5*Lx,0.5*Lx)
#		y=random.uniform(-0.5*Ly,0.5*Ly)
#		z=random.uniform(zlo,zhi)
#		print "create_atoms	%i single %1.2f %1.2f %1.2f units box" % (2,x,y,-z)
#		print "create_atoms	%i single %1.2f %1.2f %1.2f units box" % (1,x,y,z)
#		zp.append(z)
#		zp.append(-z)
#		count+=4
#for i in range((Ntot - count)/2):
#		x=random.uniform(-0.5*Lx,0.5*Lx)
#		y=random.uniform(-0.5*Ly,0.5*Ly)
#		z=random.uniform(-250,250)
#		print "create_atoms	%i single %1.2f %1.2f %1.2f units box" % (2,x,y,z)
#		print "create_atoms	%i single %1.2f %1.2f %1.2f units box" % (1,x,y,-z)	
#		zp.append(-z)
#		count+=2	
#print "##%i atoms uniformally inserted" % count
#fuck = MyHist(zp,np.array([((-Lz/2.)+i*(Lz/float(Numbins))) for i in range(Numbins)]+[Lz*0.5]))
#for i in fuck:
#	print i
