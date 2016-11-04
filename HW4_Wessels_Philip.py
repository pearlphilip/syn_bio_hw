"""
Created on Sat Oct 08 16:42:17 2016

@author: Madeline and Pearl
"""
import matplotlib.pyplot as plt
import numpy
import pandas as pd
import tellurium as te

# Question 2
r2 = te.loada(
'''
# Transcription
R1: DNA_tetR => RNA_tetR + DNA_tetR; k1 * DNA_tetR;
R2: DNA_lacI => DNA_lacI + RNA_lacI; k2 * DNA_lacI;
R3: DNA_glk => RNA_glk + DNA_glk; k3 * DNA_glk;

# Translation
R4: RNA_lacI => lacI; k4 * RNA_lacI;
R5: RNA_glk => glk; k5 * RNA_glk;

# Inhibition
R6: lacI + DNA_glk => lacI_DNA_glk_complex; k6 * lacI * DNA_glk;
R7: lacI_DNA_glk_complex => DNA_glk + lacI; k7 * lacI_DNA_glk_complex;
R8: RNA_tetR + DNA_lacI => RNA_tetR_DNA_lacI_complex; k8 * RNA_tetR * DNA_lacI;
R9: RNA_tetR_DNA_lacI_complex => DNA_lacI + RNA_tetR; k9 * RNA_tetR_DNA_lacI_complex;
R10: aTc + RNA_tetR => aTc_RNA_tetR_complex; k10 * aTc * RNA_tetR;
R11: aTc_RNA_tetR_complex => RNA_tetR + aTc; k11 * aTc_RNA_tetR_complex;
R12: IPTG + lacI => IPTG_lacI_complex; k12 * IPTG * lacI;
R13: IPTG_lacI_complex => lacI + IPTG; k13 * IPTG_lacI_complex;

# Dilution due to cell division
R14: RNA_tetR => ; k_cell_div;
R15: RNA_lacI => ; k_cell_div;
R16: RNA_tetR => ; k_cell_div;
R17: lacI => ; k_cell_div;
R18: glk => ; k_cell_div;

# Variables that get set by python code
DNA_lacI = 1; DNA_glk = 1; DNA_tetR = 1;
RNA_tetR = 0; RNA_lacI = 0; RNA_glk = 0;
lacI = 0; glk = 0; 
lacI_DNA_glk_complex = 0; RNA_tetR_DNA_lacI_complex = 0; 
aTc_RNA_tetR_complex = 0; IPTG_lacI_complex = 0;

# Dilution occurs as the cell divides, therefore k_cell_div
# is calculated from standard cell doubling time
k_cell_div = 1.28 *(10 ^(-4)); 

k1 = 10; k2 = 10; k3 = 10;
k4 = 10; k5 = 10;
k6 = 10; k7 = 10; k8 = 10; k9 = 10; k10 = 10; k11 = 10;
k12 = 10; k13 = 10;
k14 = 10; k15 = 10; k16 = 10; k17 = 10; k18 = 10;

''')
aTc = 0; IPTG = 0;

# get the stoichiometry matrix
print r2.getFullStoichiometryMatrix()
matrix = pd.DataFrame((r2.getFullStoichiometryMatrix()))
matrix.columns = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'R10', 
                  'R11', 'R12', 'R13', 'R14', 'R15', 'R16', 'R17', 'R18']
'''matrix = matrix.reindex(['DNA_tetR', 'RNA_tetR', 'DNA_lacI', 'RNA_lacI', 
                           'DNA_glk', 'RNA_glk', 'lacI', 'glk', 
                           'lacI_DNA_glk_complex', 'RNA_tetR_DNA_lacI_complex',
                           'aTc', 'aTc_RNA_tetR_complex', 'IPTG', 
                           'IPTG_lacI_complex'])'''
matrix.to_csv('./HW4_Stoich_Mat.csv')

# Question 3
modelstring3 = (
'''
# Transcription
R1: DNA_tetR => RNA_tetR + DNA_tetR; k1 * DNA_tetR;
R2: DNA_lacI => DNA_lacI + RNA_lacI; k2 * DNA_lacI;
R3: DNA_glk => RNA_glk + DNA_glk; k3 * DNA_glk;

# Translation
R4: RNA_lacI => lacI; k4 * RNA_lacI;
R5: RNA_glk => glk; k5 * RNA_glk;

# Inhibition
R6: lacI + DNA_glk => lacI_DNA_glk_complex; k6 * lacI * DNA_glk;
R7: lacI_DNA_glk_complex => DNA_glk + lacI; k7 * lacI_DNA_glk_complex;
R8: RNA_tetR + DNA_lacI => RNA_tetR_DNA_lacI_complex; k8 * RNA_tetR * DNA_lacI;
R9: RNA_tetR_DNA_lacI_complex => DNA_lacI + RNA_tetR; k9 * RNA_tetR_DNA_lacI_complex;
R10: aTc + RNA_tetR => aTc_RNA_tetR_complex; k10 * aTc * RNA_tetR;
R11: aTc_RNA_tetR_complex => RNA_tetR + aTc; k11 * aTc_RNA_tetR_complex;
R12: IPTG + lacI => IPTG_lacI_complex; k12 * IPTG * lacI;
R13: IPTG_lacI_complex => lacI + IPTG; k13 * IPTG_lacI_complex;

# Dilution due to cell division
R14: RNA_tetR => ; k_cell_div;
R15: RNA_lacI => ; k_cell_div;
R16: RNA_tetR => ; k_cell_div;
R17: lacI => ; k_cell_div;
R18: glk => ; k_cell_div;

# Variables that get set by python code
DNA_lacI = 1; DNA_glk = 1; DNA_tetR = 1;
RNA_tetR = 0; RNA_lacI = 0; RNA_glk = 0;
lacI = 0; glk = 0; 
lacI_DNA_glk_complex = 0; RNA_tetR_DNA_lacI_complex = 0; 
aTc_RNA_tetR_complex = 0; IPTG_lacI_complex = 0;

# Dilution occurs as the cell divides, therefore k_cell_div
# is calculated from standard cell doubling time
k_cell_div = 1.28 *(10 ^(-4)); 

# We need to vary k1 to k18 till we get a sigmoid dependancy of glk expression 
# levels on increasing aTc concentrations

k1 = 10;
k2 = 90;
k3 = 80;
k4 = 70;
k5 = 60;
k6 = 50;
k7 = 45;
k8 = 40;
k9 = 35;
k10 = 35; 
k11 = 25;
k12 = 20;
k13 = 15;
k14 = 10;
k15 = 8;
k16 = 6;
k17 = 4;
k18 = 2;
''')
modelstring_aTc_IPTG = modelstring3 + "aTc = 0; IPTG = 0;"

r3 = te.loada(modelstring_aTc_IPTG)
model3 = r3.simulate(0, 1000, steps=1000)
plt.plot(model3[:,0], model3[:,1], label='DNA_tetR')
plt.plot(model3[:,0], model3[:,2], label='RNA_tetR')
plt.plot(model3[:,0], model3[:,3], label='RNA_lacI')
plt.plot(model3[:,0], model3[:,4], label='DNA_glk')
plt.plot(model3[:,0], model3[:,5], label='RNA_glk')
plt.plot(model3[:,0], model3[:,6], label='lacI')
plt.plot(model3[:,0], model3[:,7], label='glk')
plt.plot(model3[:,0], model3[:,8], label='lacI_DNA_glk_complex')
plt.plot(model3[:,0], model3[:,9], label='RNA_tetR_DNA_lacI_complex')
plt.plot(model3[:,0], model3[:,10], label='aTc')
plt.plot(model3[:,0], model3[:,11], label='aTc_RNA_tetR_complex')
plt.plot(model3[:,0], model3[:,12], label='IPTG')
plt.plot(model3[:,0], model3[:,13], label='IPTG_lacI_complex')
plt.title("All species outputs from a simple model of the inverter under min \
aTc and IPTG")
plt.xlabel('Time in seconds')
plt.ylabel('Species density')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()

plt.plot(model3[:,0], model3[:,7], label='glk')
plt.title("Glk output from a simple model of the inverter under min \
aTc and IPTG")
plt.xlabel('Time in seconds')
plt.ylabel('Species density')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()
    
aTc_matrix = numpy.linspace(0, 2.16*(10**(-7)), num=6); 
IPTG_matrix = numpy.linspace(0, 10**(-4), num=6);
glk_aTc = []
for i in aTc_matrix:
    aTc_value = i
    IPTG_value = 0
    modelstring_aTc_IPTG = modelstring3 + "aTc = %d; IPTG = %d;" % (aTc_value, IPTG_value)
    r3 = te.loada(modelstring_aTc_IPTG)
    r3.timeCourseSelections = ['time', '[glk]']
    model3 = r3.simulate(0, 27*60*60, steps=1000)
    glk_aTc.append(model3[-1, 1])
plt.plot(aTc_matrix, glk_aTc)
plt.title('Variation of glk activity levels with aTc concentrations')
plt.xlabel('aTc concentrations')
plt.ylabel('glk species density')
plt.show()

glk_IPTG = []
for j in IPTG_matrix:
    IPTG_value = j
    aTc_value = 0
    modelstring_aTc_IPTG = modelstring3 + "aTc = %d; IPTG = %d;" % (aTc_value, IPTG_value)
    r3 = te.loada(modelstring_aTc_IPTG)
    r3.timeCourseSelections = ['time', '[glk]']
    model3 = r3.simulate(0, 27*60*60, steps=1000)
    glk_IPTG.append(model3[-1, 1])
plt.plot(IPTG_matrix, glk_IPTG)
plt.title('Variation of glk activity levels with IPTG concentrations')
plt.xlabel('IPTG concentrations')
plt.ylabel('glk species density')
plt.show()

# Question 4
modelstring4 = (
'''
# Transcription
R1: DNA_tetR => RNA_tetR + DNA_tetR; k1 * DNA_tetR;
R2: DNA_lacI => DNA_lacI + RNA_lacI; k2 * DNA_lacI;
R3: DNA_glk => RNA_glk + DNA_glk; k3 * DNA_glk;

# Translation
R4: RNA_lacI => lacI; k4 * RNA_lacI;
R5: RNA_glk => glk; k5 * RNA_glk;
# Want to measure glk; plot of y axis

# Inhibition
R6: lacI + DNA_glk => lacI_DNA_glk_complex; k6 * lacI * DNA_glk;
R7: lacI_DNA_glk_complex => DNA_glk + lacI; k7 * lacI_DNA_glk_complex;
R8: RNA_tetR + DNA_lacI => RNA_tetR_DNA_lacI_complex; k8 * RNA_tetR * DNA_lacI;
R9: RNA_tetR_DNA_lacI_complex => DNA_lacI + RNA_tetR; k9 * RNA_tetR_DNA_lacI_complex;
R10: aTc + RNA_tetR => aTc_RNA_tetR_complex; k10 * aTc * RNA_tetR;
R11: aTc_RNA_tetR_complex => RNA_tetR + aTc; k11 * aTc_RNA_tetR_complex;
R12: IPTG + lacI => IPTG_lacI_complex; k12 * IPTG * lacI;
R13: IPTG_lacI_complex => lacI + IPTG; k13 * IPTG_lacI_complex;

# Dilution due to cell division
R14: RNA_tetR =>; k_cell_div;
R15: RNA_lacI => ; k_cell_div;
R16: RNA_glk => ; k_cell_div;
R17: lacI => ; k_cell_div;
R18: glk => ; k_cell_div;

# Degradation 
R19: lacI => ; kdeg * lacI;

# Variables that get set by python code
DNA_lacI = 1; DNA_glk = 1; DNA_tetR = 1;
RNA_tetR = 0; RNA_lacI = 0; RNA_glk = 0;
lacI = 0; glk = 0;
lacI_DNA_glk_complex = 0; RNA_tetR_DNA_lacI_complex = 0; 
aTc_RNA_tetR_complex = 0; IPTG_lacI_complex = 0; 
aTc = 0; IPTG = 0;

# Dilution occurs as the cell divides, therefore k_cell_div
# is calculated from standard cell doubling time
k_cell_div = 1.28 *(10 ^(-4)); 
k1 = 10; k2 = 10; k3 = 10;
k4 = 10; k5 = 10;
k6 = 10; k7 = 10; k8 = 10; k9 = 10; k10 = 10; k11 = 10;
k12 = 10; k13 = 10;
k14 = 10; k15 = 10; k16 = 10; k17 = 10; k18 = 10;

''')
kdeg_matrix = numpy.linspace(10**(-6), 10, num=10);
glk_kdeg = []
for i in kdeg_matrix:
    kdeg_value = i
    modelstring_kdeg = modelstring4 + "kdeg = %d;" % (kdeg_value)
    r4 = te.loada(modelstring_kdeg)
    r4.timeCourseSelections = ['time', '[glk]']
    model4 = r4.simulate(0, 100, steps=100)
    glk_kdeg.append(model4[-1, 1])
plt.plot(kdeg_matrix, glk_kdeg)
plt.title('Variation of glk activity levels with lacI degradation rates')
plt.xlabel('lacI deg rate constant')
plt.ylabel('glk species density')
plt.show()