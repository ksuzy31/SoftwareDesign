# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 13:51:32 2014

@author: ksuzy31
"""

#print "aaaaaaaaaabaaaaaaaaaabaaaaaaaaaabaaaaaaaaaabaaaaaaaaaabdaaaaaaaaaacdaaaaaaaaaacdaaaaaaaaaacdaaaaaaaaaacdaaaaaaaaaac"

def nuts (x):
    if int(x)<=100 and int(x)>=0:
        print "hello"
    elif int(x)>100 and int(x)<500:
        print "goodbye"
    elif int(x)>=600 and int(x)<=1000:
        print "ciao"
    else:
        print "You doing this wrong"


def table():
    x=int(raw_input("What is x?\n"))
    y=raw_input('What is y?\n')
    if (int(x)>=0 and int(x)<=100) or y=='01052':
        print "Yes"


def hypotenuse ():
    a = int(raw_input("What is a?\n"))    
    b = int(raw_input("What is b?\n")) 
    c=math.sqrt(a**2+b**2)
    return c
    
def right_justify (s):
    x= " "*(70-len(s))+s
    return x

def eat_pie():
     return 'I like eating pie'
    
def do_twice(f):
    x=f()
    y=f()
    return(x,y)
    