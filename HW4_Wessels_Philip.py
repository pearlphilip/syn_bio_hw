"""
Created on Sat Oct 08 16:42:17 2016

@author: Madeline and Pearl
"""
import matplotlib.pyplot as plt
import os
import pandas as pd
import pprint
import tellurium as te

# Question 2
r2 = te.loada('''
                
                ''')

# get the stoichiometry matrix
print(r2.getFullStoichiometryMatrix())

# Question 3
modelstring3 =  '''
                
                '''

r3 = te.loada(modelstring3)
r3.timeCourseSelections = ['time', '[FAEE]']
model3 = r3.simulate(0, 4500, steps=1000)
r3.plot(model3, title="FAEE trend in DSRS control system")