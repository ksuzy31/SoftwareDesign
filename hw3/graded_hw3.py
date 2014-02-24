# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 11:24:42 2014

@author: YOUR NAME HERE
kevin suzuki 
"""
import random
aAcids = ['F','L','I','M','V','S','P','T','A','Y',
      '|','H','Q','N','K','D','E','C','W','R',
      'G']

codon = [['TTT', 'TTC'],
          ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
          ['ATT', 'ATC', 'ATA'],
          ['ATG'],
          ['GTT', 'GTC', 'GTA', 'GTG'],
          ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
          ['CCT', 'CCC', 'CCA', 'CCG'],
          ['ACT', 'ACC', 'ACA', 'ACG'],
          ['GCT', 'GCC', 'GCA', 'GCG'],
          ['TAT', 'TAC'],
          ['TAA', 'TAG', 'TGA'],
          ['CAT', 'CAC'],
          ['CAA', 'CAG'],
          ['AAT', 'AAC'],
          ['AAA', 'AAG'],
          ['GAT', 'GAC'],
          ['GAA', 'GAG'],
          ['TGT', 'TGC'],
          ['TGG'],
          ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
          ['GGT', 'GGC', 'GGA', 'GGG']]
              

# you may find it useful to import these variables (although you are not required to use them)
#from amino_acids import aa, codons

def collapse(L):
    """ Converts a list of strings to a string by concatenating all elements of the list """
    output = ""
    for s in L:
        output = output + s
    return output


def coding_strand_to_AA(rawin):
    if len(rawin) %3 == 0:
        print rawin
    else:
        print 'bad'
    m=[]
    amino=''
    for i in range (len(rawin)/3):
        v= rawin[i*3:i*3+3]
        m.append(v)  
    for a in range (len(m)):
        for x in range (len(codon)):
            temp=codon[x]
            for n in range (len(temp)):
                if temp[n] == m[a]:
                    amino+=aAcids[x]
    return amino
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).
        
        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
    """
    from amino_acids import codons


def coding_strand_to_AA_unit_tests():
    """ Unit tests for the coding_strand_to_AA function """
        
    a = "AAA"
    b = "CAGCGTTGGATGCAA"
    c = "ATAATGTGTAATCA"
    aeo = "K"
    beo = "QRWMQ"
    ceo = "IMCN"
    print "input: " + a + ", expected output: " + aeo + ", actual output: " + coding_strand_to_AA(a)
    print "input: " + b + ", expected output: " + beo + ", actual output: " + coding_strand_to_AA(b)
    print "input: " + c + ", expected output: " + ceo + ", actual output: " + coding_strand_to_AA(c)

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence
    
        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    """
    v=""
    j=""    
    for p in range(len(dna)):
        if dna[p]=='A':
            v+='T'
        if dna[p]=='T':
            v+='A'
        if dna[p]=='G':
            v+='C'
        if dna[p]=='C':
            v+='G'
    for q in range (len(v)):
        j+=v[-(q+1)]
    return j
    
    
    
def get_reverse_complement_unit_tests():
    """ Unit tests for the get_complement function """
    
    a = "AAA"
    b = "CAGCGTTGGATGCAA"
    c = "ATAATGTGTAATCA"
    aeo = "TTT"
    beo = "TTGCATCCAACGCTG"
    ceo = "TGATTACACATTAT"
    print "input: " + a + ", expected output: " + aeo + ", actual output: " + get_reverse_complement(a)
    print "input: " + b + ", expected output: " + beo + ", actual output: " + get_reverse_complement(b)
    print "input: " + c + ", expected output: " + ceo + ", actual output: " + get_reverse_complement(c)    
    
def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start codon and returns
        the sequence up to but not including the first in frame stop codon.  If there
        is no in frame stop codon, returns the whole string.
        
        dna: a DNA sequence
        returns: the open reading frame represented as a string
    """
    r=[]
    m=""
    stop=codons[10]
    for i in range (len(dna)/3):
        v= dna[i*3:i*3+3]
        if v == stop[0] or v== stop[1] or v==stop[2]:        
            break        
        m+=(v)
    return m

def rest_of_ORF_unit_tests():
    """ Unit tests for the rest_of_ORF function """
        
    a = "ATG"
    b = "ATGTTTTGATGAA"
    c = "ATGATGTGTAGATCA"
    aeo = "ATG"
    beo = "ATGTTT"
    ceo = "ATGATGTGTAGATCA"
    print "input: " + a + ", expected output: " + aeo + ", actual output: " + rest_of_ORF(a)
    print "input: " + b + ", expected output: " + beo + ", actual output: " + rest_of_ORF(b)
    print "input: " + c + ", expected output: " + ceo + ", actual output: " + rest_of_ORF(c)    
       
        
def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence and returns
        them as a list.  This function should only find ORFs that are in the default
        frame of the sequence (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """

    m=[]
    n=0
    while n<=len(dna):
        v= dna[n*3:n*3+3]
        if v == 'ATG':
            x=dna[n*3:len(dna)]
            dnastring=rest_of_ORF(x)
            m+=[dnastring]
            n+=len(dnastring)/3
        n+=1
    return m
    
            
     
def find_all_ORFs_oneframe_unit_tests():
    """ Unit tests for the find_all_ORFs_oneframe function """

    a = "ATGGGGGGGTGA"
    b = "ATGTTTTAAATGAAATAG"
    c = "ATGGGGTAGGGGATGCCCTAA"
    aeo = "[ATGGGGGGG]"
    beo = "[ATGTTT,ATGAAA]"
    ceo = "[ATGGGG,ATGCCC]"
    print "input: " + a + ", expected output: " + aeo + ", actual output: " + str(find_all_ORFs_oneframe(a))
    print "input: " + b + ", expected output: " + beo + ", actual output: " + str(find_all_ORFs_oneframe(b))
    print "input: " + c + ", expected output: " + ceo + ", actual output: " + str(find_all_ORFs_oneframe(c))
    

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in all 3
        possible frames and returns them as a list.  By non-nested we mean that if an
        ORF occurs entirely within another ORF and they are both in the same frame,
        it should not be included in the returned list of ORFs.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
    p=dna[1:len(dna)]
    o=dna[2:len(dna)]     
    firstframe=find_all_ORFs_oneframe(dna)
    secondframe=find_all_ORFs_oneframe(p)
    thirdframe=find_all_ORFs_oneframe(o)
    return(firstframe,secondframe,thirdframe)
        
     

def find_all_ORFs_unit_tests():
    """ Unit tests for the find_all_ORFs function """
        
    a = "ATGGATGGTGAATGAGTGA"
    b = "ATGTTTTAATGAAATAG"
    c = "ATGGGGTAGTGAGGTATGAGCTACCTATGAA"
    aeo = "[ATGGATGGTGA,ATGGTGAATGAG,ATGAGTGA]"
    beo = "[ATGTTT,ATGAAA]"
    ceo = "[ATGGGG,ATGAGCTACCTA,ATGAA]"
    print "input: " + a + ", expected output: " + aeo + ", actual output: " + str(find_all_ORFs(a))
    print "input: " + b + ", expected output: " + beo + ", actual output: " + str(find_all_ORFs(b))
    print "input: " + c + ", expected output: " + ceo + ", actual output: " + str(find_all_ORFs(c))  



def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.
        
        dna: a DNA sequence
        returns: a list of non-nested ORFs
    """
    p=get_reverse_complement(dna)
    reverse=find_all_ORFs(p)
    strand=find_all_ORFs(dna)
    return(strand, reverse)

def find_all_ORFs_both_strands_unit_tests():
    """ Unit tests for the find_all_ORFs_both_strands function """

a = "ATGGATGGTGAATGAGTGA"
b = "ATGTTTTAATGAAATAG"
c = "ATGGGGTAGTGAGGTATGAGCTACCTATGAA"
aeo = "[ATGGATGGTGA,ATGGTGAATGAG,ATGAGTGA]"
beo = "[ATGTTT, ATGAAA]"
ceo = "[ATGGGG,ATGAGCTACCTA,ATGAA]"
print "input: " + a + ", expected output: " + aeo + ", actual output: " + str(find_all_ORFs_both_strands(a))
print "input: " + b + ", expected output: " + beo + ", actual output: " + str(find_all_ORFs_both_strands(b))
print "input: " + c + ", expected output: " + ceo + ", actual output: " + str(find_all_ORFs_both_strands(c))

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string"""    
    t=find_all_ORFs_both_strands(dna)
    a=t[0]
    q=a[0]
    if not q:
        longest=""
    else:
        longest=q[0]
    for i in range (len(t)):
        temp=t[i]        
        for b in range (len(temp)):
            temp2=temp[b]            
            for x in range (len(temp2)):        
                if (len(temp2[x]))>(len(longest)):
                    longest=temp2[x]
    return longest
    

def longest_ORF_unit_tests():
    """ Unit tests for the longest_ORF function """

    a = "ATGGATGGTGAATGAGTGA"
    b = "ATGTTTTAATGAAATAG"
    c = "ATGGGGTAGTGAGGTATGAGCTACCTATGAA"
    aeo = "ATGGATGGTGAA"
    beo = "ATGTTT"
    ceo = "ATGAGCTACCTA"
    print "input: " + a + ", expected output: " + aeo + ", actual output: " + longest_ORF(a)
    print "input: " + b + ", expected output: " + beo + ", actual output: " + longest_ORF(b)
    print "input: " + c + ", expected output: " + ceo + ", actual output: " + longest_ORF(c)


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence
        
        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    k=[]
    o=0
    v=0
    while v < num_trials:
        for x in range (len(dna)):
            k+=[dna[x]]
        shuffle (k)
        h=collapse(k)
        k=[]
        c=longest_ORF(h)
        v+=1
        if o < len(c):
            o=len(c)
    return (o)

def gene_finder(dna, threshold):
    """ Returns the amino acid sequences coded by all genes that have an ORF
        larger than the specified threshold.
        
        dna: a DNA sequence
        threshold: the minimum length of the ORF for it to be considered a valid
                   gene.
        returns: a list of all amino acid sequences whose ORFs meet the minimum
                 length specified.
    """
    s=[]
    t=find_all_ORFs_both_strands(dna)
    for i in range (len(t)):
        temp=t[i]        
        for b in range (len(temp)):
            temp2=temp[b]            
            for x in range (len(temp2)):        
                if (len(temp2[x]))>(threshold):
                    s+=[coding_strand_to_AA(temp2[x])]
    return s