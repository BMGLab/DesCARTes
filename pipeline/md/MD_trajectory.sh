# step 1: PBC correction and centering
gmx_mpi trjconv -s step7_production_500ns.tpr \
                -f step7_production_500ns.xtc \
                -o step7_pbc_fixed.xtc \
                -pbc mol -center -ur compact


#Center: Protein (enter)
#Output : System (enter)

# step 2: Fitting (rotation and translation correction)
gmx_mpi trjconv -s step7_production_500ns.tpr \
                -f step7_pbc_fixed.xtc \
                -o step7_production_fixed.xtc \
                -fit rot+trans


#Fit: Backbone (enter)
#Output : System (enter)


#RMSD:
bashgmx_mpi rms -s step7_production_500ns.tpr \
            -f step7_production_fixed.xtc \
            -o rmsd_backbone_fixed.xvg \
            -tu ns

#first group: Backbone
#second group: Backbone

#RMSF:
bashgmx_mpi rmsf -s step7_production_500ns.tpr \
             -f step7_production_fixed.xtc \
             -o rmsf_fixed.xvg \
             -res

#Backbone or C-alpha

#Radius of Gyration:
bashgmx_mpi gyrate -s step7_production_500ns.tpr \
               -f step7_production_fixed.xtc \
               -o rg_fixed.xvg
# Protein



#PyMOL'de Görüntüleme
# python PyMOL'de:

load step7_production_fixed.xtc
load step5_input.pdb