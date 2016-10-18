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
                R1: => FadR; k5
                R2: => Pdc_AdhB; k15
                R3: => AtfA; k16
                R4: => TesA; k17
                R5: => Fad; k18

                // FAEE biosynthesis
                R6: => FA; k1
                    // TesA => TesA + FA; TesA is a catalyst
                R7: FA => acyl-CoA; FA * k2
                    // FadD + FA => FadF + acyl-CoA; FA is a catalyst
                R8: ETOH + acyl-CoA => FAEE; ETOH * acyl-CoA * k3
                    // AtfA + ETOH + acyl-CoA => AtfA + FAEE; AtfA is a catalyst
                R9: => ETOH; k4
                    // Pdc_adhB => Pdc_AdhB + ETOH; Pdc_adhB is a catalyst

                // Note: FADR-ligand binding and FadR-promoter binding constants
                // can be found by the following equation
                // kx * (Ki / (Ki + [IC]))
                // where [IC] is inhibitor concentration
                // Ki is enzyme inhibition constant
                // See description of Supplementary figure 6 for more info

                // FadR-ligand binding
                R10: FadR + FA => FadR-FA; FadR * FA * k6
                R11: FadR + acyl-CoA => FadR-acyl-CoA; FadR * acyl-CoA * k7
                R12: FadR-FA => FadR + FA; FadR-DA * k8_FA
                R13: FadR-acyl-CoA => FadR + acyl-CoA; FadR-acyl-CoA * k8_acyl-CoA

                // FadR-Promoter binding
                R14: FadR + P2 => FadR-P2; FadR * P2 * k9
                R15: FadR-P2 => FadR + P2; FadR-P2 * k10
                R16: FadR + P3 => FadR-P3; FadR * P3 * k11
                R17: FadR-P3 => FadR + P3; FadR-P3 * k12
                R18: FadR + P4 => FadR-P4; FadR * P4 * k13
                R19: FadR-P4 => FadR + P4; FadR-P4 * k14

                // Note: Reactions R19-R37 account for loss of species due to
                // cell-growth associated dilution, degradation, and secretion
                // from the cell

                // Cell growth/Species dilution
                R20: TesA => ; k19
                R21: Pdc_AdhB => ; k19
                R22: AtfA => ; k19
                R23: FaD => ; k19
                R24: FadR => ; k19
                R25: FadR-FA => ; k19
                R26: FadR-acyl-CoA => ; k19
                R27: P1 => ; k20
                R28: P2 => ; k20
                R29: P3 => ; k20
                R30: P4 => ; k20
                R31: P5 => ; k20
                R32: FadR-P2 => ; k20
                R33: FadR-P3 => ; k20
                R34: FadR-P4 => ; k20
                R35: FA => ; k21
                R36: acyl-CoA => ; k22
                R37: ETOH => ; k23

                // DNA replication
                R38: P1 => P1 + P1; k24
                R39: P2 => P2 + P2; k24
                R40: P3 => P3 + P3; k24
                R41: P4 => P4 + P4; k24
                R42: P5 => P5 + P5; k24
                R43: FadR-P2 = > FadR-P2 + FadR-P2; k24
                R44: FadR-P3 = > FadR-P3 + FadR-P3; k24
                R45: FadR-P4 = > FadR-P4 + FadR-P4; k24
                
                P1 = 1.2 * (10 ^ (-4)); 
                P2 = 5.2 * (10 ^ (-3)); 
                P3 = 5.4 * (10 ^ (-3));
                P4 = 6.0 * (10 ^ (-3)); 
                P5 = 3.7 * (10 ^ (-2));

                // Set rate constants

                // These can be calculated using the following equation
                // 1 / (Px_init + (ORF / kpol))
                //  where Px_init is transcription initiation rate
                // and kpol is a rate for coupled transcription and traslation
                // and ORF is open reading frame
                // See supplementary figure 6 description for more info

                // Note k5...k19 are first order for gene expression from DNA species P1 to P5
                k5 =
                k15 =
                k16 =
                k17 =
                k18 =
                k19 =

                k1 = 4.9 * (10 ^ (1))
                k2 = 5.5 * (10 ^ (-3))
                k3 = 3.4 * (10 ^ (-3))
                k4 = 3.8 * (10 ^ (2))

                // k6 and k7 are constants for association
                // k8_Fa and k8_acyl-CoA are constants for dissociation
                k6 =
                k7 =
                k8_FA =
                k8_acyl-CoA =

                // k9,11,13 are constants for association
                // k10,12,14 are constants for dissociation
                k9 =
                k10 =
                k11 =
                k12 =
                k13 =
                k14 =

                k19 =
                k20 =
                k21 =
                k22 =
                k23 =

                k24 =

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
