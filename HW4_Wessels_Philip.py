"""
Created on Sat Oct 08 16:42:17 2016

@author: Madeline and Pearl
"""
import matplotlib.pyplot as plt
import pandas as pd
import tellurium as te

# Question 2
r2 = te.loada('''
                # Transcription
                R1: => RNA_tetR; k1;
                R2: => RNA_lacI; k2;
                R3: => RNA_glk; k3;
                
                k1 = 10; k2 = 10; k3 = 10;
                
                # Translation                
                R4: RNA_lacI => lacI; k4 * RNA_lacI;
                R5: RNA_glk => glk; k5 * RNA_glk;
                
                k4 = 10; k5 = 10;
            
                # Degradation, dilution and cell division
                R6: RNA_tetR => ; k_cell_div
                R7: RNA_lacI => ; k_cell_div
                R8: RNA_tetR => ; k_cell_div
                R9: lacI => ; k_cell_div
                R10: glk => ; k_cell_div 
 
                #variables that get set by python code
                RNA_tetR = 0; RNA_lacI = 0; RNA_glk = 0;
                tetR = 0; lacI = 0; glk = 0; 
                
                # Dilution occurs as the cell divides, therefore k_cell_div = 
                # k_dilution and is calculated from standard cell doubling
                # time
                k_cell_div = 1.28 *(10 ^(-4)); k_dilution = 1.28 *(10 ^(-4));
                ''')

# get the stoichiometry matrix
print r2.getFullStoichiometryMatrix()
matrix = pd.DataFrame((r2.getFullStoichiometryMatrix()))
matrix.to_csv('./HW4_Stoich_Mat.csv')

# Question 3
modelstring3 = ('''
                # Transcription
                R1: => RNA_tetR; k1;
                R2: => RNA_lacI; k2;
                R3: => RNA_glk; k3;
                
                k1 = 10; k2 = 10; k3 = 10;
                
                # Translation                
                R4: RNA_lacI => lacI; k4 * RNA_lacI;
                R5: RNA_glk => glk; k5 * RNA_glk;
                
                k4 = 10; k5 = 10;
            
                # Degradation, dilution and cell division
                R6: RNA_tetR => ; k_cell_div
                R7: RNA_lacI => ; k_cell_div
                R8: RNA_tetR => ; k_cell_div
                R9: lacI => ; k_cell_div
                R10: glk => ; k_cell_div 
 
                #variables that get set by python code
                RNA_tetR = 0; RNA_lacI = 0; RNA_glk = 0;
                tetR = 0; lacI = 0; glk = 0; 
                
                # Dilution occurs as the cell divides, therefore k_cell_div = 
                # k_dilution and is calculated from standard cell doubling
                # time
                k_cell_div = 1.28 *(10 ^(-4)); k_dilution = 1.28 *(10 ^(-4));
                ''')

r3 = te.loada(modelstring3)
r3.timeCourseSelections = ['time', '[glk]']
model3 = r3.simulate(0, 400, steps=1000)
r3.plot(model3, title="glk trend in a system without aTc")