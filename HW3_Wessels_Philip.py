# -*- coding: utf-8 -*-
"""
Created on Sat Oct 08 16:42:17 2016

@author: pearl
"""
import tellurium as te


# Question 3
modelstring3 = '''
                model question3()

                '''

r3 = te.loada(modelstring3)
model3 = r3.simulate(0, 4500, 1000)
r3.plot(model3, title=" ")
