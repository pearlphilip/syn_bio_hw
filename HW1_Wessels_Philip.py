# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 13:14:07 2016

@author: Madeline and Pearl
"""
import tellurium as te

### Question 2
modelstring1 = '''
                model question2()
                
                // The two processes have been collapsed into one expression
                // DNA -> protein = ( rateOfRNA/DNA * DNA ) * rateOfProtein/RNA
                
                k1 = 1/10; k2 = 1/20; k3 = 1/30;
                k4 = 1/10; k5 = 1/100; k6 = 1/1000;
                J1: -> S3; (k4 * S1) * k1;
                J2: -> S4; (k4 * S1) * k2;
                J3: -> S5; (k4 * S1) * k3;
                J4: -> S6; (k5 * S1) * k1;
                J5: -> S7; (k5 * S1) * k2;
                J6: -> S8; (k5 * S1) * k3;
                J7: -> S9; (k6 * S1) * k1;
                J8: -> S10; (k6 * S1) * k2;
                J9: -> S11; (k6 * S1) * k3;
                
                // S1 is DNA, S2 is RNA, S3-S11 are protein counts
                // Initial amounts of each species:
                S1 = 1; S2 = 0; S3 = 0; S4 = 0; S5 = 0; S6 = 0;
                S7 = 0; S8 = 0; S9 = 0; S10 = 0; S11 = 0;
                
                end
                '''

r1 = te.loada(modelstring1)

### Question 3
model1 = r1.simulate(0,18 * 60 * 60,10000) 
# End time: 18 hr * 60 min/hr * 60 s/min
r1.plot(model1)

### Question 4
# Re-use model code above, removing the head line "model question2()" and replacing with our new model
modelstring1_body_info = '\n'.join(modelstring1.split('\n')[2:])

modelstring2 = '''
                model question4()                
                
                at (time % (20 * 60) == 0): S1 = 2*S1 
                // DNA doubles every 20 minutes (in seconds)
                
                // Protein degrades every hour (in seconds) 
                // We use a rate constant to denote this

                k7 = 1/(1*60*60);
                J10: S3 -> ; k7 * S3
                J11: S4 -> ; k7 * S4
                J12: S5 -> ; k7 * S5
                J13: S6 -> ; k7 * S6
                J14: S7 -> ; k7 * S7
                J15: S8 -> ; k7 * S8
                J16: S9 -> ; k7 * S9
                J17: S10 -> ; k7 * S10
                J18: S11 -> ; k7 * S11
             ''' 

modelstring2 += modelstring1_body_info

r2 = te.loada(modelstring2)
model2 = r2.simulate(0,18 * 60 * 60,10000) 
# End time: 18 hr * 60 min/hr * 60 s/min
r2.plot(model2)
