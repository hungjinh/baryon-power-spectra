import illustris_python as il
import numpy as np


#import sys
#sys.path.insert(0,'/home/hungjinh/python_lib/CAMB/pycamb')
#import camb
#from camb import model, initialpower

temp_save_dir="/home/hungjinh/Research/baryon_proj/code/temp_storage/"


def build_den_cube_DMO(snap,basePath,boxsize,resol):

	# boxsize [Mpc]
	# Mptl : DM particle mass [Msun]

	fields=["Coordinates"]
	dm_pos = il.snapshot.loadSubset(basePath,snap,'dm',fields=fields);
	dm_pos=dm_pos/1000     # change unit to Mpc/h

	hist,edges = np.histogramdd(dm_pos, resol, range=[[0,boxsize],[0,boxsize],[0,boxsize]])

	avgN_per_pix=np.float(len(dm_pos))/resol**3
	den = (hist - avgN_per_pix)/avgN_per_pix

	## save a backup file...
	#	#den=np.load(temp_save_dir+"den_DMO.npy")
	#global temp_save_dir
	#np.save(temp_save_dir+"ill1_DMO_den.npy", den)

	return den


def cal_Amp_FFTden(den,resol):

	FFTden = np.fft.fftn(den)
	Amp_FFTden=np.absolute(FFTden)/resol**3

	## save a backup file...
	#global temp_save_dir
	#np.save(temp_save_dir+"Amp_FFTden.npy", Amp_FFTden)

	return Amp_FFTden


def cal_freq_FFT(boxsize,resol):
	
	delx=boxsize/np.double(resol)    
	fxyz  = np.fft.fftfreq(resol,d=delx)[np.mgrid[0:resol,0:resol,0:resol]]
	#kxyz  = fxyz*2*np.pi
	#kx=kxyz[0] ; ky=kxyz[1] ; kz=kxyz[2]
	Wk = np.sinc(fxyz[0]*delx)*np.sinc(fxyz[1]*delx)*np.sinc(fxyz[2]*delx)  
	k  = np.sqrt(fxyz[0]**2 + fxyz[1]**2 + fxyz[2]**2)*2*np.pi
	
	print "k_max:" ,np.max(k)
	print "kx_max:",np.max(fxyz[0])*2*np.pi
	print "kx_min:",np.min(fxyz[0])*2*np.pi

	return k,Wk



def cal_Pk(kbins,Amp_FFTden,k,Vbox):
	
	Nbins=len(kbins)-1
	avgPk = np.empty(Nbins,dtype=float)
	avgk  = np.empty(Nbins,dtype=float)
	Nk    = np.empty(Nbins,dtype=float)

	for j in range(Nbins):
		takeout_ID=np.where((k > kbins[j]) & (k <= kbins[j+1]))
		k_now = k[takeout_ID]
		Amp_FFTden_now=Amp_FFTden[takeout_ID]
		FFTden2_now  = Amp_FFTden_now * Amp_FFTden_now

		avgk[j]   = np.mean(k_now)
		avgPk[j]  = np.mean(FFTden2_now)*Vbox
		Nk[j]     = len(k_now)

	return avgk,avgPk,Nk


def gen_camb_Pk_illustris_cosmology(redshifts=[0]):

	h=0.704
	Omega_b=0.0456
	Omega_l=0.7274
	Omega_m=0.2726   # total matter density
	ns=0.963
	sigma8=0.809
	As=2.174e-9

	Omega_c=Omega_m-Omega_b

	ombh2=Omega_b*h**2
	omch2=Omega_c*h**2

	pars = camb.CAMBparams()
	pars.set_cosmology(H0=100*h, ombh2=ombh2, omch2=omch2, mnu=0., omk=0, tau=0.06)
	pars.InitPower.set_params(ns=ns, r=0, As=As)

	## matching sigma8 by tuning As in above
	pars.set_matter_power(redshifts=redshifts, kmax=2.0)
	results = camb.get_results(pars)
	#s8 = np.array(results.get_sigma8())
	#print s8/sigma8

	#Non-Linear spectra (Halofit)
	pars.NonLinear = model.NonLinear_both
	results.calc_power_spectra(pars)
	kh_nonlin, z_nonlin, pk_nonlin = results.get_matter_power_spectrum(minkh=0.0316, maxkh=39, npoints = 1000)

	print z_nonlin

	return kh_nonlin,pk_nonlin

def gen_camb_Pk_MBII_cosmology(redshifts=[0]):

	h=0.701
	Omega_b=0.046
	Omega_l=0.725
	Omega_m=0.275   # total matter density
	ns=0.968
	sigma8=0.816
	As=2.208e-9

	Omega_c=Omega_m-Omega_b

	ombh2=Omega_b*h**2
	omch2=Omega_c*h**2

	pars = camb.CAMBparams()
	pars.set_cosmology(H0=100*h, ombh2=ombh2, omch2=omch2, mnu=0., omk=0, tau=0.06)
	pars.InitPower.set_params(ns=ns, r=0, As=As)

	## matching sigma8 by tuning As in above
	pars.set_matter_power(redshifts=redshifts, kmax=2.0)
	results = camb.get_results(pars)
	#s8 = np.array(results.get_sigma8())
	#print s8/sigma8

	#Non-Linear spectra (Halofit)
	pars.NonLinear = model.NonLinear_both
	results.calc_power_spectra(pars)
	kh_nonlin, z_nonlin, pk_nonlin = results.get_matter_power_spectrum(minkh=0.0316, maxkh=30, npoints = 1000)

	return kh_nonlin,pk_nonlin



def build_mass_cube_hydro_DM(snap,basePath,Mptl,boxsize,resol):

	fields=["Coordinates"]
	dm_pos = il.snapshot.loadSubset(basePath,snap,'dm',fields=fields);
	dm_pos=dm_pos/1000     # change unit to Mpc/h

	Mptl=[Mptl]*len(dm_pos)

	mass_cube,edges = np.histogramdd(dm_pos, resol, range=[[0,boxsize],[0,boxsize],[0,boxsize]],weights=Mptl)

	#global temp_save_dir
	#np.save(temp_save_dir+"MassCube_DM.npy", mass_cube)

	return mass_cube


def build_mass_cube_hydro_baryon(snap,basePath,ptl_type,boxsize,resol):

	fields     =["Coordinates","Masses"]
	baryon_cat = il.snapshot.loadSubset(basePath,snap,ptl_type,fields=fields);
	baryon_pos = baryon_cat["Coordinates"]/1000.   # Mpc
	
	Mptl=baryon_cat["Masses"]*1e10 

	mass_cube,edges = np.histogramdd(baryon_pos, resol, range=[[0,boxsize],[0,boxsize],[0,boxsize]],weights=Mptl)


	#global temp_save_dir
	#np.save(temp_save_dir+"MassCube_"+ptl_type+".npy", mass_cube)

	return mass_cube



