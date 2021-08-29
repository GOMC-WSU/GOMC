# GOMC Example for the Gibbs Ensemble (GEMC) using MoSDeF [1, 2, 5-10, 13-17]

# Note: In this specific example, we will be using the GEMC_NVT ensemble.


# Import the required packages and specify the force field (FF) being used. 

# Note: For GOMC, the residue names are treated as molecules, so the residue names must be unique for each different molecule. [1, 2, 13-17]

# Note: Each residue can be set to a different FF, which is done by setting the residue name to a FF in a dictionary (FF_Dict).  The FF selection can be a FF name (set from foyer FF repositor) or a specified FF xml file. [1, 2, 13-17]

import mbuild as mb
from foyer import Forcefield
import mbuild.formats.charmm_writer as mf_charmm
import mbuild.formats.gomc_conf_writer as gomc_control


FF_file_water = '../common/spce.xml'
water = mb.load('O', smiles=True)
water.name = 'H2O'
water.energy_minimize(forcefield=FF_file_water, steps=10**5)

FF_file_ethanol = 'oplsaa'
ethanol = mb.load('CCO', smiles=True)
ethanol.name = 'ETO'
ethanol.energy_minimize(forcefield=FF_file_ethanol, steps=10**5)

FF_dict = {water.name: FF_file_water, ethanol.name: FF_file_ethanol}

residues_list = [ethanol.name, water.name]

fix_bonds_angles_residues = [water.name]


# Build the main simulation liquid box (box 0) and the vapor (box 1) for the simulation [1, 2, 13-17]


water_ethanol_box_liq = mb.fill_box(compound=[water, ethanol],
                                    density= 950,
                                    compound_ratio=[0.8, 0.2] ,
                                    box=[3.0, 3.0, 3.0])

water_ethanol_box_vap = mb.fill_box(compound=[water,ethanol],
                                    density= 100,
                                    compound_ratio=[0.8, 0.2],
                                    box=[8, 8, 8])


## Build the Charmm object, which is required to write the FF (.inp), psf, pdb, and GOMC control files [1, 2, 5-10, 13-17]

## The reorder_res_in_pdb_psf command reorders the psf and pdb to the order residues variable (i.e., the residues_list in this case) [1, 2, 13-17].  

## The fix_res_bonds_angles command fixes the angles and bonds for water in the Charmm FF file.  Note: This is specific to GOMC, as it sets the bond and angle k-values to 999999999999 [1, 2, 5-10, 13-17].


charmm = mf_charmm.Charmm(water_ethanol_box_liq,
                          'GEMC_NVT_water_ethanol_liq',
                          structure_box_1=water_ethanol_box_vap,
                          filename_box_1='GEMC_NVT_water_ethanol_vap',
                          ff_filename="GEMC_NVT_water_ethanol_FF",
                          forcefield_selection=FF_dict,
                          residues=residues_list,
                          bead_to_atom_name_dict=None,
                          fix_residue=None,
                          gomc_fix_bonds_angles=fix_bonds_angles_residues,
                          reorder_res_in_pdb_psf=True
                          )


# Write the write the FF (.inp), psf, pdb, and GOMC control files [1, 2, 5-10, 13-17]

charmm.write_inp()

charmm.write_psf()

charmm.write_pdb()


gomc_control.write_gomc_control_file(charmm, 'in_GCMC.conf',  'GCMC', 100, 300,
                                     input_variables_dict = {
                                         "ChemPot" : {"H2O" : -4000, "ETO" : -8000}
                                     }
                                     )



