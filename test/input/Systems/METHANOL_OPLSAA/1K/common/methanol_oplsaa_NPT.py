import mbuild as mb
from foyer import Forcefield
import mbuild.formats.charmm_writer as mf_charmm
import mbuild.formats.gomc_conf_writer as gomc_control
from mbuild.lib.recipes.polymer import Polymer

seed_number = 4


methanol_FF = 'oplsaa'

methanol = mb.load('CO', smiles=True)
methanol.name = 'MTO'

residues_list = [methanol.name]


# Build the main simulation liquid box (box 0) and the vapor (box 1) for the simulation [1, 2, 13-17]

methanol_box = mb.fill_box(compound=methanol,
                           density=825,
                           box=[3, 3, 3],
                           seed=seed_number
                           )
#left_box_liq.energy_minimize(forcefield=DTPA_methanol_FF, steps=10**5)
print("INFO: Finished Building the DTPA and methanol box ")


#empty_graphene_pore.save('test.xyz', overwrite=True)

print("INFO: Started Building the Charmm Object ")
# build the charmm object
charmm = mf_charmm.Charmm(methanol_box,
                          'methanol_box_0',
                          structure_box_1=methanol_box,
                          filename_box_1='methanol_box_1',
                          ff_filename="GOMC_oplsaa_methanol_FF",
                          forcefield_selection=methanol_FF,
                          residues=residues_list,
                          bead_to_atom_name_dict=None,
                          fix_residue=None,
                          gomc_fix_bonds_angles=None,
                          reorder_res_in_pdb_psf=False,
                          )
print("INFO: Finished Building the Charmm Object ")

# Write the write the FF (.inp), psf, pdb, and GOMC control files [1, 2, 5-10, 13-17]

print("INFO: Started Writing the PDB, PSF, and FF files ")
charmm.write_inp()

charmm.write_psf()

charmm.write_pdb()
print("INFO: Finished Writing the PDB, PSF, and FF files ")



print("INFO: Started Writing the GOMC Control File ")
# write the GOMC control file
gomc_output_data_every_X_steps = 1 * 10**4
output_true_list_input = [
    True,
    int(gomc_output_data_every_X_steps),
]
output_false_list_input = [
    False,
    int(gomc_output_data_every_X_steps),
]

gomc_output_coord_every_X_steps = 1 * 10**5
output_coord_true_list_input = [
    True,
    int(gomc_output_data_every_X_steps),
]

gomc_control.write_gomc_control_file(charmm, 'in_GCMC.conf', 'GCMC', 60*10**6, 298,
                                     input_variables_dict={"ChemPot": {methanol.name: -3000},
                                                           "OutputName": "GCMC_output",
                                                           "VDWGeometricSigma": True,
                                                           "Potential": "SWITCH",
                                                           "LRC": False,
                                                           "RcutLow": 0.0,
                                                           "Rcut": 16,
                                                           "Rswitch": 14,
                                                           "RcutCoulomb_box_0": 14,
                                                           "PRNG": seed_number,
                                                           "VolFreq": 0.00,
                                                           "MultiParticleFreq": 0.01,
                                                           "DisFreq": 0.19,
                                                           "RotFreq": 0.2,
                                                           "IntraSwapFreq": 0.1,
                                                           "SwapFreq": 0.20,
                                                           "RegrowthFreq": 0.15,
                                                           "CrankShaftFreq": 0.15,
                                                           "EqSteps": 1000,
                                                           "PressureCalc": output_true_list_input,
                                                           "RestartFreq": output_true_list_input,
                                                           "ConsoleFreq": output_true_list_input,
                                                           "BlockAverageFreq": output_true_list_input,
                                                           "HistogramFreq": output_true_list_input,
                                                           "CheckpointFreq": output_coord_true_list_input,
                                                           "CoordinatesFreq": output_coord_true_list_input,
                                                           "DCDFreq": output_coord_true_list_input ,
                                                           "CBMC_First": 12,
                                                           "CBMC_Nth": 10,
                                                           "CBMC_Ang": 50,
                                                           "CBMC_Dih": 50,
                                                           }
                                    )
print("INFO: Finished Writing the GOMC Control File ")


