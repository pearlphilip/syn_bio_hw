# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 13:14:07 2016

@author: Madeline and Pearl
"""
import tellurium as te

# Question 2
modelstring1 = '''
                model question2()
                
                // DNA -> protein = ( rateOfRNA/DNA * DNA ) * rateOfProtein/RNA
                J1: -> S3; (1/10 * S1) * 1/10;
                J2: -> S3; (1/10 * S1) * 1/20;
                J3: -> S3; (1/10 * S1) * 1/30;
                J4: -> S3; (1/100 * S1) * 1/10;
                J5: -> S3; (1/100 * S1) * 1/20;
                J6: -> S3; (1/100 * S1) * 1/30;
                J7: -> S3; (1/1000 * S1) * 1/10;
                J8: -> S3; (1/1000 * S1) * 1/20;
                J9: -> S3; (1/1000 * S1) * 1/30;
                
                // S1 is DNA, S2 is RNA, S3 is protein
                S1 = 1; S2 = 0; S3 = 0;
                
                end
                '''

r1 = te.loada(modelstring1)

# Question 3
model = r1.simulate(0,18 * 60 * 60,10000) 
# End time: 18 hr * 60 min/hr * 60 s/min
r1.plot(model)

# Question 4
modelstring2 = '''
                model question4()                
                
                at (time % (60 * 20) == 0): S1 = 2*S1 
                // doubles DNA every 20 minutes
                at (time % (60 * 60) == 0): S3 = S3 - 1 
                // degrades protein every hour
                
                
                // DNA -> protein = ( rateOfRNA/DNA * DNA ) * rateOfProtein/RNA
                J1: -> S3; (1/10 * S1) * 1/10;
                J2: -> S3; (1/10 * S1) * 1/20;
                J3: -> S3; (1/10 * S1) * 1/30;
                J4: -> S3; (1/100 * S1) * 1/10;
                J5: -> S3; (1/100 * S1) * 1/20;
                J6: -> S3; (1/100 * S1) * 1/30;
                J7: -> S3; (1/1000 * S1) * 1/10;
                J8: -> S3; (1/1000 * S1) * 1/20;
                J9: -> S3; (1/1000 * S1) * 1/30;
                
                // S1 is DNA, S2 is RNA, S3 is protein
                S1 = 1; S2 = 0; S3 = 0;
                
                protein = S3;
                end
                '''

r2 = te.loada(modelstring2)
model = r2.simulate(0,18 * 60 * 60,10000) 
# End time: 18 hr * 60 min/hr * 60 s/min
r2.plot(model)