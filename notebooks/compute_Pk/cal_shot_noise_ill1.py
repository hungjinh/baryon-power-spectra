import numpy as np
import sys
import time


snap=int(sys.argv[1])

data_dir="/home/hungjinh/Research/baryon_proj/data/hydro_box/illustris1/snap%03d/"%snap
save_dir="/home/hungjinh/Research/baryon_proj/data/Pk_data/ill1_hydro/shot_noise_hydro/"
print snap


Lbox = 75.  # Mpc/h
resol = 1024
Vbox = 75.**3  # Mpc/h
Nptl_dmo = 1820.**3

Tst=time.time()

print "----- Now dm particle -----"
dm_pos = np.load(data_dir+'dm_pos.npy')
Nptl_dm = len(dm_pos)
Mptl_dm = 4.40896524361e6   # Msun/h
Mptl_dm_list = np.array([Mptl_dm]*Nptl_dm,dtype='float64')

print "Nptl_dm:", Nptl_dm
print "Mptl_dm_list[0:20]:",Mptl_dm_list[0:20]



print "----- Now gas particle -----"
Mptl_gas_list = np.load(data_dir+'gas_Mptl.npy')    # Msun/h
Nptl_gas = len(Mptl_gas_list)

print "Nptl_gas:", Nptl_gas
print "Mptl_gas_list[0:20]:",Mptl_gas_list[0:20]



print "----- Now star particle -----"
Mptl_star_list = np.load(data_dir+'star_Mptl.npy')    # Msun/h
Nptl_star = len(Mptl_star_list)

print "Nptl_star:", Nptl_star
print "Mptl_star_list[0:20]:",Mptl_star_list[0:20]



print "----- Now bh particle -----"
Mptl_bh_list = np.load(data_dir+'bh_Mptl.npy')    # Msun/h
Nptl_bh = len(Mptl_bh_list)

print "Nptl_bh:", Nptl_bh
print "Mptl_bh_list[0:20]:",Mptl_bh_list[0:20]


##### shot noise calculation #####
Sigma_Mptl2 = np.sum(Mptl_dm_list**2)+np.sum(Mptl_gas_list**2)+np.sum(Mptl_star_list**2)+np.sum(Mptl_bh_list**2)

Sigma_Mptl = np.sum(Mptl_dm_list)+np.sum(Mptl_gas_list)+np.sum(Mptl_star_list)+np.sum(Mptl_bh_list)


print "Sigma_Mptl2:",Sigma_Mptl2
print "Sigma_Mptl:",Sigma_Mptl
print "Vbox:",Vbox

shot_noise_hydro = Vbox*Sigma_Mptl2/(Sigma_Mptl**2)

print "shot_noise_hydro:",shot_noise_hydro

np.save(save_dir+"shot_noise_hydro_%03d.npy"%snap,shot_noise_hydro)

print "shot_noise_dmo:",Vbox/Nptl_dmo

Ttot = (time.time()-Tst)/60.  ; print "total shot noise cal. time:",Ttot,"mins"
