# -*- coding: utf-8 -*-
"""
Created on Sat Oct 08 16:42:17 2016

@author: pearl
"""
import tellurium as te


### Question 1b
modelstring1b =  '''
                model question1b()
                
                k1 = 1; k2 = 1; k3 = 1;
                J1: A + B -> 2B; k1 * A * B;
                J2: B + C -> 2C; k2 * B * C;
                J3: C + A -> 2A; k3 * C * A;

                A = 1; B = 1; C = 1;

                end
                '''

r1 = te.loada(modelstring1b);
model1 = r1.simulate(0, 10, 1000);
r1.plot(model1);

### Question 1c, Trial 1
modelstring1c1 =  '''
                model question1c1()
                
                k1 = 2; k2 = 2; k3 = 2;
                J1: A + B -> 2B; k1 * A * B;
                J2: B + C -> 2C; k2 * B * C;
                J3: C + A -> 2A; k3 * C * A;

                A = 1; B = 1; C = 1;

                end
                '''

r2 = te.loada(modelstring1c1);
model2 = r2.simulate(0, 10, 1000);
r2.plot(model2);

### Question 1c, Trial 2
modelstring1c2 =  '''
                model question1c2()
                
                k1 = 2; k2 = 1; k3 = 1;
                J1: A + B -> 2B; k1 * A * B;
                J2: B + C -> 2C; k2 * B * C;
                J3: C + A -> 2A; k3 * C * A;

                A = 1; B = 1; C = 1;

                end
                '''

r3 = te.loada(modelstring1c2);
model3 = r3.simulate(0, 10, 1000);
r3.plot(model3);

### Question 1c, Trial 3
modelstring1c3 =  '''
                model question1c3()
                
                k1 = 2; k2 = 3; k3 = 5;
                J1: A + B -> 2B; k1 * A * B;
                J2: B + C -> 2C; k2 * B * C;
                J3: C + A -> 2A; k3 * C * A;

                A = 1; B = 1; C = 1;

                end
                '''

r4 = te.loada(modelstring1c3);
model4 = r4.simulate(0, 10, 1000);
r4.plot(model4);
