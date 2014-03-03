# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 11:34:57 2014

@author: kevin suzuki
"""


import random
from math import *
import Image


def build_random_function(min_depth, max_depth):
    
    """This function takes min_depth and max_depth as upper and lower bound inputs and creates a function that returns a random number of lists based upon the inputs"""
    
    p=[['x'],['y']]
    if max_depth<=1:
        return p[random.randint(0,1)]
    a=build_random_function(min_depth-1,max_depth-1)
    b=build_random_function(min_depth-1,max_depth-1)
    product=['prod',a,b]    
    sin_pi=['sin_pi',a]
    cos_pi=['cos_pi',a]
    squared=['^2',a]
    cubed=['^3',a]
    flist=[product,sin_pi,cos_pi,squared,cubed,['x'],['y']]
    if min_depth>1:
        return flist[random.randint(0,4)]
    else:
        return flist[random.randint(0,6)]

def evaluate_random_function(f, x, y):
    
   
"""This function evaluates a list generated by build_random_function(min_depth,max_depth),x and y and outputs a random value based upon the three inputs"""
    
    d=f[0]
    if d=='x':
        return x
    if d=='y':
        return y
    if d=='prod':
        return evaluate_random_function(f[1],x,y)*evaluate_random_function(f[2],x,y)
    if d=='sin_pi':
        return sin(pi*evaluate_random_function(f[1],x,y))
    if d=='cos_pi':
        return cos(pi*evaluate_random_function(f[1],x,y))
    if d=='^2':
        return (evaluate_random_function(f[1],x,y))**2
    if d=='^3':
        return (evaluate_random_function(f[1],x,y))**3


        
def remap_interval(val, input_interval_start, input_interval_end, output_interval_start, output_interval_end):
    """ Maps the input value that is in the interval [input_interval_start, input_interval_end]
        to the output interval [output_interval_start, output_interval_end].  The mapping
        is an affine one (i.e. output = input*c + b).
    
        TODO: please fill out the rest of this docstring
    """
    
    """this code allows you to remap the input interval and maps it to an output interval. It will find the change in the interval. 
    """
    
    val=float(val)
    input_interval_start=float(input_interval_start)
    input_interval_end=float(input_interval_end)
    output_interval_start=float(output_interval_start)
    output_interval_end=float(output_interval_end)
    x=val-input_interval_start
    y=input_interval_end-input_interval_start
    z=output_interval_end-output_interval_start
    c=output_interval_start
    f=(x/y)*(z)+c
    return f


def draw(xsize,ysize,min_depth,max_depth):
    
    """This code will take 4 inputs: xsize,ysize,min_depth and max_depth as inputs and will return a visual image based off the combination of the 4 inputs."""
    xsize=int(xsize)
    ysize=int(ysize)
    red=build_random_function(min_depth, max_depth)
    blue=build_random_function(min_depth, max_depth)
    green=build_random_function(min_depth, max_depth)
    im = Image.new("RGB",(int(xsize),int(ysize)))
    for i in range (int(xsize)):
        x= remap_interval(i,0.0,xsize,-1.0,1.0)
        for j in range (int(ysize)):
            y=remap_interval(i,0.0,ysize,-1.0,1.0)

            rd=evaluate_random_function(red,x,y)
            bl=evaluate_random_function(blue,x,y)
            gn=evaluate_random_function(green,x,y)
            redmap=remap_interval(rd,-1,1,0,255)
            bluemap=remap_interval(bl,-1,1,0,255)
            greenmap=remap_interval(gn,-1,1,0,255)
            im.putpixel((i,j),(int(redmap),int(bluemap),int(greenmap)))
    im.save("image"+str(number)+".bmp")
