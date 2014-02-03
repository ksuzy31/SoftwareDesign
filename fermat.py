# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 00:08:21 2014

@author: ksuzy31
"""
def box():
    print ('+ ' + '- '*4 +'+ '+ '- '*4 + '+')+((('\n'+'|'+' '*9 +'|'+ ' '*9+'|')*4+('\n'+'+ ' + '- '*4 +'+ '+ '- '*4 + '+'))*2)   #this will print a 2*2 box
   
box ()

def box2 (n):
    print (('+ ' +'- '*4)*n+'+\n'+(('|'+' '*9)*n+'|\n')*4)*n+('+ ' +'- '*4)*n+'+'  # this will print an n by n matrix 
box2(n) #input a number for n in the terminal 
