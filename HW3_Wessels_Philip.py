"""
Created on Sat Oct 08 16:42:17 2016

@author: Madeline and Pearl
"""
import tellurium as te

# Question 2
r2 = te.loada('''
R10: FadR + FA => FadR_FA; k6 * FadR * FA 
R11: FadR + acyl_CoA => FadR_acyl_CoA; k7 * FadR * acyl_CoA
R12: FadR_FA => FadR + FA; k8_FA * FadR_FA
R13: FadR_acyl_CoA => FadR + acyl_CoA; k8_acyl_CoA * FadR_acyl_CoA
R16: FadR + P3 => FadR_P3; k11 * FadR * P3
R17: FadR_P3 => FadR + P3; k12 * FadR_P3
R18: FadR + P4 => FadR_P4; k13 * FadR * P4
R19: FadR_P4 => FadR + P4; k14 * FadR_P4

k6 = 1; k7 = 2; k8_FA = 3; k8_acyl_CoA = 4; 
k11 = 5; k12 = 6; k13 = 7; k14 = 8;

FA = 5;
acyl_CoA = 8;
FadR = 4;
FadR_FA = 3;
FadR_acyl_CoA = 7;
P3 = 2;
FadR_P3 = 1;
P4 = 9;
FadR_P4 = 6;
''')

# reaction rate vector
print(r2.getReactionRates())

# get concentrations of X & Y
print(r2.getFloatingSpeciesConcentrations())

# get the stoichiometry matrix
print(r2.getFullStoichiometryMatrix())

# Question 3
modelstring3 =  '''
                model question3()

                // Gene expression
                R1: => FadR; K5
                R2: => Pdc_AdhB; K15
                R3: => AtfA; K16
                R4: => TesA; K17
                R5: => FadD; K18
                
                // Note K5, K15 to K18 are first order rate constants for gene 
                // expression from DNA species P1 to P5, calculated in the HW 
                // document
                K5 = 3.48 * (10 ** (-2))
                K15 = 9.31 * (10 ** (-3))
                K16 = 1.43 * (10 ** (-2))
                K17 = 3.99 * (10 ** (-2))
                K18 = 1.48 * (10 ** (-2))

                // FAEE biosynthesis
                R6: => FA; K1
                    // TesA => TesA + FA; TesA is a catalyst
                R7: FA => acyl_CoA; FA * K2
                    // FadD + FA => FadF + acyl_CoA; FA is a catalyst
                R8: ETOH + acyl_CoA => FAEE; ETOH * acyl_CoA * K3
                    // AtfA + ETOH + acyl_CoA => AtfA + FAEE; AtfA is a catalyst
                R9: => ETOH; K4
                    // Pdc_adhB => Pdc_AdhB + ETOH; Pdc_adhB is a catalyst
                    
                K1 = 4.9 * (10 ** (1))
                K2 = 5.5 * (10 ** (-3))
                K3 = 3.4 * (10 ** (-3))
                K4 = 3.8 * (10 ** (2))

                // FadR_ligand binding
                R10: FadR + FA => FadR_FA; FadR * FA * K6
                R11: FadR + acyl_CoA => FadR_acyl_CoA; FadR * acyl_CoA * K7
                R12: FadR_FA => FadR + FA; FadR_DA * K8_FA
                R13: FadR_acyl_CoA => FadR + acyl_CoA; FadR_acyl_CoA * K8_acyl_CoA
                
                // K6 and K7 are constants for association
                // K8_FA and K8_acyl_CoA are constants for dissociation, 
                // source shown in the HW document
                K6 = 1.9 * (10 ** (-2))
                K7 = 1.9 * (10 ** (-2))
                K8_FA = 9.3 * (10 ** (-8))
                K8_acyl_CoA = 9.3 * (10 ** (-11))
                
                // FadR_Promoter binding
                R14: FadR + P2 => FadR_P2; FadR * P2 * K9
                R15: FadR_P2 => FadR + P2; FadR_P2 * K10
                R16: FadR + P3 => FadR_P3; FadR * P3 * K11
                R17: FadR_P3 => FadR + P3; FadR_P3 * K12
                R18: FadR + P4 => FadR_P4; FadR * P4 * K13
                R19: FadR_P4 => FadR + P4; FadR_P4 * K14
                
                // K9,11,13 are constants for association
                // K10,12,14 are constants for dissociation
                // source shown in the HW document
                K9 = 1.9 * (10 ** (-2))
                K10 = 1.3 * (10 ** (-1))
                K11 = 1.9 * (10 ** (-2))
                K12 = 1.3 * (10 ** (-1))
                K13 = 1.9 * (10 ** (-2))
                K14 = 1.3 * (10 ** (-1))

                // Reactions R20 to R37 account for loss of species due to
                // cell growth associated dilution, degradation, and secretion
                // from the cell

                // Cell growth/Species dilution
                R20: TesA => ; K19
                R21: Pdc_AdhB => ; K19
                R22: AtfA => ; K19
                R23: FadD => ; K19
                R24: FadR => ; K19
                R25: FadR_FA => ; K19
                R26: FadR_acyl_CoA => ; K19
                R27: P1 => ; K20
                R28: P2 => ; K20
                R29: P3 => ; K20
                R30: P4 => ; K20
                R31: P5 => ; K20
                R32: FadR_P2 => ; K20
                R33: FadR_P3 => ; K20
                R34: FadR_P4 => ; K20
                R35: FA => ; K21
                R36: acyl_CoA => ; K22
                R37: ETOH => ; K23
                
                K19 = 6.3 * (10 ** (-5))
                K20 = 6.3 * (10 ** (-5))
                K21 = 6.3 * (10 ** (-5))
                K22 = 6.3 * (10 ** (-5))
                K23 = 6.3 * (10 ** (-5))

                // DNA replication
                
                R38: P1 => P1 + P1; K24
                R39: P2 => P2 + P2; K24
                R40: P3 => P3 + P3; K24
                R41: P4 => P4 + P4; K24
                R42: P5 => P5 + P5; K24
                R43: FadR_P2 => FadR_P2 + FadR_P2; K24
                R44: FadR_P3 => FadR_P3 + FadR_P3; K24
                R45: FadR_P4 => FadR_P4 + FadR_P4; K24
                
                P1 = 1.2 * (10 ** (-4)); 
                P2 = 5.2 * (10 ** (-3)); 
                P3 = 5.4 * (10 ** (-3));
                P4 = 6.0 * (10 ** (-3)); 
                P5 = 3.7 * (10 ** (-2));
                
                // K24 is the cell division rate constant, from calculations of  
                // cell doubling time
                K24 = 1.28 * (10 ** (-4))

                end
                '''

r3 = te.loada(modelstring3)
model3 = r3.simulate(0, 4500, 1000)
r3.plot(model3, title=" ")


# Question 4
modelstring4 = '''
                model question4()

                '''

r4 = te.loada(modelstring4)
model4 = r4.simulate(0, 4500, 1000)
r4.plot(model4, title=" ")
