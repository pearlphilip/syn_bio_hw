# -*- coding: utf-8 -*-
"""
Created on Sat Oct 08 16:42:17 2016

@author: pearl
"""
import numpy as np
import os
import pandas as pd
import tellurium as te
import matplotlib.pyplot as plt
from scipy.integrate import odeint


### Question 1b
modelstring1b =  '''
                model question1b()
                
                // When k1 = k2 = k3 = 1
                
                k1 = 1; k2 = 1; k3 = 1;
                J1: A + B -> 2B; k1 * A * B;
                J2: B + C -> 2C; k2 * B * C;
                J3: C + A -> 2A; k3 * C * A;

                A = 1; B = 1; C = 1;

                end
                '''

r1 = te.loada(modelstring1b);
model1 = r1.simulate(0, 10, 1000);
r1.plot(model1, title="Dynamics of A, B & C when k1 = 1, k2 = 1, k3 = 1");

### Question 1c, Trial 1
modelstring1c1 =  '''
                model question1c1()
                
                // When k1 = k2 = k3 != 1
                
                k1 = 2; k2 = 2; k3 = 2;
                J1: A + B -> 2B; k1 * A * B;
                J2: B + C -> 2C; k2 * B * C;
                J3: C + A -> 2A; k3 * C * A;

                A = 1; B = 1; C = 1;

                end
                '''

r2 = te.loada(modelstring1c1);
model2 = r2.simulate(0, 10, 1000);
r2.plot(model2, title="Dynamics of A, B & C when k1 = 2, k2 = 2, k3 = 2");

### Question 1c, Trial 2
modelstring1c2 =  '''
                model question1c2()
                
                // When k1 != k2 = k3
                
                k1 = 2; k2 = 1; k3 = 1;
                J1: A + B -> 2B; k1 * A * B;
                J2: B + C -> 2C; k2 * B * C;
                J3: C + A -> 2A; k3 * C * A;

                A = 1; B = 1; C = 1;

                end
                '''

r3 = te.loada(modelstring1c2);
model3 = r3.simulate(0, 10, 1000);
r3.plot(model3, title="Dynamics of A, B & C when k1 = 2, k2 = 1, k3 = 1");

### Question 1c, Trial 3
modelstring1c3 =  '''
                model question1c3()
                
                // When k1 != k2 != k3
                
                k1 = 2; k2 = 3; k3 = 5;
                J1: A + B -> 2B; k1 * A * B;
                J2: B + C -> 2C; k2 * B * C;
                J3: C + A -> 2A; k3 * C * A;

                A = 1; B = 1; C = 1;

                end
                '''

r4 = te.loada(modelstring1c3);
model4 = r4.simulate(0, 10, 1000);
r4.plot(model4, title="Dynamics of A, B & C when k1 = 2, k2 = 3, k3 = 5");

### Question 2d
def ode(X, t):
    x, y = X
    dXdt = [-0.5 * (x ** 2) * y,  0.5 * (x ** 2) * y]
    return dXdt

X0 = [1, 2]
t = np.linspace(0, 10, 101)
sol = odeint(ode, X0, t) 
plt.plot(t, sol[:, 0], 'b', label='x(t)')
plt.plot(t, sol[:, 1], 'g', label='y(t)')
plt.legend(loc='best')
plt.xlabel('Time t')
plt.ylabel('x and y concentrations')
plt.title('Dynamics of reaction system with [x] = 1 and [y] = 2 at t = 0')
plt.grid()
plt.show()

### Question 3
modelstring3 = '''
                # Protein1 is Insulin, Protein2 is GFP
                R1: -> RNA1; k1;
                R2: -> Protein1; k2 * RNA1;
                R3: -> RNA2; k3;
                R4: -> Protein2; k4 * RNA2
                R5: RNA1 -> ; k_deg1 * RNA1 
                R6: RNA2 -> ; k_deg2 * RNA1 
                R7: -> Vol; k_cell_div * Vol
 
                #variables that get set by python code
                RNA1 = 0; Protein1 = 0; RNA2 = 0; Protein2 = 0; Vol = 10
                k1 = 0; k2 = 0; k3 = 0; k4 = 0; 
                k_deg1 = 1/3600; k_deg2 = 2/3600; k_cell_div = 3/3600;
                
                '''
r3 = te.loada(modelstring3)
# Reciprocals of the rate constants
reciprocal_k1 = [10, 100, 1000]
reciprocal_k2 = [10, 20, 30]
reciprocal_k3 = [200, 300, 500]
reciprocal_k4 = [50, 100, 150]

for i in reciprocal_k1:
    for j in reciprocal_k2:
        for k in reciprocal_k3:
            for l in reciprocal_k4:
                r3.RNA1 = 0
                r3.Protein1 = 0
                r3.RNA2 = 0
                r3.Protein2 = 0
                r3.k1 = 1.0/i
                r3.k2 = 1.0/j
                r3.k3 = 1.0/k
                r3.k4 = 1.0/l
                
        
                model3 = r3.simulate(0, 100, 1000) 
                # End time: 18 hr * 60 min/hr * 60 s/min
                r3.plot(model3, title="k1 = every %is k2 = every %is k3 = every %is k4 = every %is"%(i, j, k, l))
                
                # Saving the model into a csv file
                df_model3 = pd.DataFrame(model3)
                df_model3.columns = ['Time(s)', 'RNA1', 'Protein1', 'RNA2', 'Protein2', 'Vol']
                output_filename = 'k1_%is_k2_%is_k3_%is_k4_%is.csv'%(i, j, k, l)
                output_path = os.path.join('./csv_files', output_filename)
                df_model3.to_csv(output_path)
"""
# simulate returns a numpy array
# request columns 'time','RNA','Protein','vol'
m = r.simulate(0,100,100,['time','RNA','Protein','vol'])

# plotting code
import matplotlib.pyplot as plt
# want to plot RNA and protein per cell
plt.plot(m[:,0], m[:,1]/m[:,3], label='RNA') # divide 2nd col by 4th col
plt.plot(m[:,0], m[:,2]/m[:,3], label='Protein') # divide 3rd col by 4th col
plt.legend()
plt.show()"""