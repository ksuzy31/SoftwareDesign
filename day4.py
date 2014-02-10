# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 14:53:47 2014

@author: ksuzy31
"""

def sum_1_to_n(n):
 s=0   
 for i in range (n+1):
     s=s+i
 return s

def factorial(n):
    s=1
    for i in range (1,n+1):
        s= s*i
    return s
    
def sum_integers(x,y):
    s=0    
    for i in range (x,y):
        s=x+i+y+i
    return s        

def sum_1_to_n_whileloop(n):
    s=0
    i=0
    while i<=n:
        s=s+i
        i=i+1


        