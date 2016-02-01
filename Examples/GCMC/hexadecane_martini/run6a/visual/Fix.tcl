


#first Box
package require psfgen
resetpsf 
topology ./Top_Martini.inp
segment AC16 {
 pdb T_590_K_u_5240_r6a_l_BOX_0_restart.pdb
 auto angles dihedrals
 first none
 last none
 }



# first none, last none req'd to prevent auto termination of residues.

#Read coordinates
coordpdb T_590_K_u_5240_r6a_l_BOX_0_restart.pdb  AC16


writepsf Angle_START_AC16_sys_BOX_0.psf
writepdb Angle_START_AC16_sys_BOX_0.pdb

