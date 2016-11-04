"""
Created on Sat Oct 08 16:42:17 2016

@author: Madeline and Pearl
"""
import matplotlib.pyplot as plt
import pandas as pd
import tellurium as te

# Question 2
r2 = te.loada(
'''
   # Transcription
R1: DNA_lacI => DNA_lacI + RNA_lacI; k1 * DNA_lacI;
R2: DNA_glk => RNA_glk + DNA_glk; k2 * DNA_glk;
R3: DNA_tetR => RNA_tetR + DNA_tetR; k3 * DNA_tetR

k1 = 10; k2 = 10; k3 = 10;

# Translation
R4: RNA_lacI => lacI; k4 * RNA_lacI;
R5: RNA_glk => glk; k5 * RNA_glk;
# Want to measure glk; plot of y axis

k4 = 10; k5 = 10;

# Inhibition
R6: lacI + DNA_glk => lacI_DNA_glk_complex; k6 * lacI * DNA_glk;
R7: lacI_DNA_glk_complex => DNA_glk; k7 * lacI_DNA_glk_complex;
R8: RNA_tetR + DNA_lacI => RNA_tetR_DNA_lacI_complex; k8 * RNA_tetR * DNA_lacI;
R9: RNA_tetR_DNA_lacI_complex => DNA_lacI; k9 * RNA_tetR_DNA_lacI_complex;
R10: aTc + RNA_tetR => aTc_RNA_tetR_complex; k10 * aTc * RNA_tetR;
R11: aTc_RNA_tetR_complex => RNA_tetR; k11 * aTc;

R12: IPTG + lacI => IPTG_lacI_complex; k12 * IPTG * lacI;
R13: IPTG_lacI_complex => lacI; k13 * IPTG_lacI_complex;

# Degradation, dilution and cell division
R14: RNA_tetR => ; k_cell_div
R15: RNA_lacI => ; k_cell_div
R16: RNA_tetR => ; k_cell_div
R17: lacI => ; k_cell_div
R18: glk => ; k_cell_div

#variables that get set by python code
DNA_lacI = 1; DNA_glk = 1; DNA_tetR = 1;
RNA_tetR = 0; RNA_lacI = 0; RNA_glk = 0;
tetR = 0; lacI = 0; glk = 0; 

# Dilution occurs as the cell divides, therefore k_cell_div
# is calculated from standard cell doubling time
k_cell_div = 1.28 *(10 ^(-4)); 
''')

# get the stoichiometry matrix
print r2.getFullStoichiometryMatrix()
matrix = pd.DataFrame((r2.getFullStoichiometryMatrix()))
matrix.to_csv('./HW4_Stoich_Mat.csv')

# Question 3
modelstring3 = (
'''
# Transcription
R1: DNA_lacI => DNA_lacI + RNA_lacI; k1 * DNA_lacI;
R2: DNA_glk => RNA_glk + DNA_glk; k2 * DNA_glk;
R3: DNA_tetR => RNA_tetR + DNA_tetR; k3 * DNA_tetR

k1 = 10; k2 = 10; k3 = 10;

# Translation
R4: RNA_lacI => lacI; k4 * RNA_lacI;
R5: RNA_glk => glk; k5 * RNA_glk;
# Want to measure glk; plot of y axis

k4 = 10; k5 = 10;

# Inhibition
R6: lacI + DNA_glk => lacI_DNA_glk_complex; k6 * lacI * DNA_glk;
R7: lacI_DNA_glk_complex => DNA_glk; k7 * lacI_DNA_glk_complex;
R8: RNA_tetR + DNA_lacI => RNA_tetR_DNA_lacI_complex; k8 * RNA_tetR * DNA_lacI;
R9: RNA_tetR_DNA_lacI_complex => DNA_lacI; k9 * RNA_tetR_DNA_lacI_complex;
R10: aTc + RNA_tetR => aTc_RNA_tetR_complex; k10 * aTc * RNA_tetR;
R11: aTc_RNA_tetR_complex => RNA_tetR; k11 * aTc;

R12: IPTG + lacI => IPTG_lacI_complex; k12 * IPTG * lacI;
R13: IPTG_lacI_complex => lacI; k13 * IPTG_lacI_complex;

# Degradation, dilution and cell division
R14: RNA_tetR => ; k_cell_div
R15: RNA_lacI => ; k_cell_div
R16: RNA_tetR => ; k_cell_div
R17: lacI => ; k_cell_div
R18: glk => ; k_cell_div

#variables that get set by python code
DNA_lacI = 1; DNA_glk = 1; DNA_tetR = 1;
RNA_tetR = 0; RNA_lacI = 0; RNA_glk = 0;
tetR = 0; lacI = 0; glk = 0; 

 # Dilution occurs as the cell divides, therefore k_cell_div
 # is calculated from standard cell doubling time
k_cell_div = 1.28 *(10 ^(-4)); 
''')

r3 = te.loada(modelstring3)
r3.timeCourseSelections = ['time', '[glk]']
model3 = r3.simulate(0, 400, steps=1000)
r3.plot(model3, title="glk trend")


# Question 4
modelstring4 = (
'''
# Transcription
R1: DNA_lacI => DNA_lacI + RNA_lacI; k1 * DNA_lacI;
R2: DNA_glk => RNA_glk + DNA_glk; k2 * DNA_glk;
R3: DNA_tetR => RNA_tetR + DNA_tetR; k3 * DNA_tetR

k1 = 10; k2 = 10; k3 = 10;

# Translation
R4: RNA_lacI => lacI; k4 * RNA_lacI;
R5: RNA_glk => glk; k5 * RNA_glk;
# Want to measure glk; plot of y axis

k4 = 10; k5 = 10;

# Inhibition
R6: lacI + DNA_glk => lacI_DNA_glk_complex; k6 * lacI * DNA_glk;
R7: lacI_DNA_glk_complex => DNA_glk; k7 * lacI_DNA_glk_complex;
R8: RNA_tetR + DNA_lacI => RNA_tetR_DNA_lacI_complex; k8 * RNA_tetR * DNA_lacI;
R9: RNA_tetR_DNA_lacI_complex => DNA_lacI; k9 * RNA_tetR_DNA_lacI_complex;
R10: aTc + RNA_tetR => aTc_RNA_tetR_complex; k10 * aTc * RNA_tetR;
R11: aTc_RNA_tetR_complex => RNA_tetR; k11 * aTc;
# k7 is the degredation rate of lacI; we vary this and plot on x axis

R12: IPTG + lacI => IPTG_lacI_complex; k12 * IPTG * lacI;
R13: IPTG_lacI_complex => lacI; k13 * IPTG_lacI_complex;

# Degradation, dilution and cell division
R14: RNA_tetR =>; k_cell_div
R15: RNA_lacI => ; k_cell_div
R16: RNA_glk => ; k_cell_div
R17: lacI => ; k_cell_div
R18: glk => ; k_cell_div

#variables that get set by python code
DNA_lacI = 1; DNA_glk = 1; DNA_tetR = 1;
RNA_tetR = 0; RNA_lacI = 0; RNA_glk = 0;
tetR = 0; lacI = 0; glk = 0; 

# Dilution occurs as the cell divides, therefore k_cell_div
# is calculated from standard cell doubling time
k_cell_div = 1.28 *(10 ^(-4)); 
''')

r4 = te.loada(modelstring4)
r4.timeCourseSelections = ['time', '[glk]']
model4 = r4.simulate(0, 400, steps=1000)
r4.plot(model4, title="Level of glk enzyme as a function of lacI degradation rate")
