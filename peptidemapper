# -*- coding: utf-8 -*-
"""
@author: Karl Barber
"""

import csv

OLS_file='FILE.csv' ##Supplementary Data 1 with 'K' added to end of <=31aa sequence column
DNAList=[]
def OLS_reader(x):
    txtReader=open(x)
    DNAReader=csv.reader(txtReader, delimiter=',')

    for element in DNAReader:
        DNAList.append(element[17])

OLS_reader(OLS_file)


file1='MAXQUANT_1.csv' ##MaxQuant search results file (all .raw files uploaded on Proteomexchange ID PXD008707)
file2='MAXQUANT_2.csv' ##MaxQuant search results file (all .raw files uploaded on Proteomexchange ID PXD008707)

choice1=''
choice1=input('Do you want to look at results obtained from samples made with a) SepOTS b) supD c)both?\nResponse: ')

choice2=''
choice2=input('Do you want to look at a) ERLIC results b) TiO2 results or c) all results?\nResponse: ')
if choice2=='a':
    choice_dec='ERLIC'
if choice2=='b':
    choice_dec='TiO2'
if choice2=='c':
    choice_dec=''

protList=[]
def compiled_file_start(x):
    txtReader=open(x)
    protReader=csv.reader(txtReader, delimiter=',')
    
    if choice1=='a':  
        for element in protReader:
            if 'Sep' in element[21] and choice_dec in element[21] or 'Sequence' in element[0]:
                if 'MULTI-MATCH' not in element[18]:
                    protList.append(element)
    if choice1=='b':
        for element in protReader:
            if 'Ser' in element[21] and choice_dec in element[21] or 'Sequence' in element[0]:
                if 'MULTI-MATCH' not in element[18]:
                    protList.append(element)
    if choice1=='c':
        for element in protReader:
            if choice_dec in element[21] or 'Sequence' in element[0]:
                if 'MULTI-MATCH' not in element[18]:
                    protList.append(element)
    
print('Initiating evidence file compilation 1/2...')
compiled_file_start(file1)
print('Initiating evidence file compilation 2/2...')
compiled_file_start(file2)


protListGlobal=protList


class MS_hits(object):
    
    def __init__(self):
        self.protList=protListGlobal        
        self.protListTemp=self.protList
        self.hits=[]
        self.counter=0
        self.breakcount=0
        
    def seqImport(self): #imports desired column as 'hits' list
        print('Importing '+str(self.seq_type)+' column from compiled evidence file...')        
        for i in self.protList[0]:
            if str(i)==str(self.seq_type):
                self.breakcount=1
                break
            self.counter+=1
        if self.breakcount==0:
            raise ValueError('Invalid column name')
        for i in self.protListTemp:
            self.hits.append(i[self.counter])
            
    def removeContam(self): #run before seqImport to remove contaminants
        print('Removing all E. coli proteins, contaminants, and reverse hits...')
        self.protListTemp=[]
        for i in self.protList[1:]:
            if i[16][0]=='s':
                self.protListTemp.append(i)

class ModSeq(MS_hits): #to trim Modified Sequences for analysis
    
    def __init__(self):
        MS_hits.__init__(self)
        self.modHits=[]
        self.seq_type='Modified sequence'
        
    def ModTrim(self): #removes underscores from modified sequences
        print('Trimming modified sequences for analysis...')        
        self.modHits.append(self.hits[0]) #imports unaltered header
        for i in self.hits:
            self.modHits.append(i[1:-1])
            
    def SepFilter(self): #returns only unique pSer hits
        self.SepList=[]
        self.tempStr=''
        print('Calculating the number of unique Ser phosphopeptides...')
        for i in self.modHits:
            self.tempStr=i.replace('(ac)','')
            self.tempStr=self.tempStr.replace('(de)','')
            self.tempStr=self.tempStr.replace('(ox)','')
            self.tempStr=self.tempStr.replace('T(ph)','T')
            self.tempStr=self.tempStr.replace('Y(ph)','Y')
            if 'S(ph)' in i:
                self.SepList.append(self.tempStr.replace('S(ph)','s')) #pSer reported as s
        self.SepList=set(self.SepList)
        self.SepList=list(self.SepList)
        self.SepList.sort()
        print('The number of unique Ser phosphopeptides is: '+str(len(self.SepList)))
        
    def SepMap(self):
        self.uniqueHits=[]
        self.unmappedSepHits=[]
        self.mappedSepHits=[]
        self.mappedSepHits_PEPTIDES=[]
        print('Initiating Sep mapping process. This may take a while...')
        for i in self.SepList:
            self.IDcounter=0
            self.breaker=0
            for j in DNAList:
                self.IDcounter+=1
                ##if i[:-2] in j: #discounts last residues in case artificial tryptic present
                if i in j:
                    self.mappedSepHits.append(j)
                    self.mappedSepHits_PEPTIDES.append(i)
                    self.uniqueHits.append(self.IDcounter)
                    self.breaker=1
                    break
            if self.breaker==0:
                self.unmappedSepHits.append(i)
        self.uniqueHits=set(self.uniqueHits)
        self.uniqueHits=list(self.uniqueHits)
        
                
        
class ProtCountDiv(MS_hits): 
    """ 
    Input data consists of one or multiple protein IDs per line separated by ';'
    makeListDiv attribute splits multiple protins IDs into single entries
        appended to same list
    It then removes duplicates, removes contaminants, and sorts by ID number
    It returns the number of unique OLS proteins identified
    """

    def __init__(self):
        MS_hits.__init__(self)
        self.divCount=1
        self.seq_type='Leading Razor Protein' #This is used for MOST STRINGENT analysis, and altered for other less stringent forms

    def makeListDiv(self): #to divide entries such as 'proteins' delineated by ';'
        print('Calculating number of unique proteins by '+str(self.seq_type)+' column...')        
        for i in self.hits[1:]:
            if ';' in i: #splits multiple entries at ';', adds them to hits list
                self.hits[self.divCount]=self.hits[self.divCount].split(';')
                while len(self.hits[self.divCount])>=2:
                    self.hits.append(self.hits[self.divCount][-1])
                    self.hits[self.divCount].pop()
            self.divCount+=1
        self.divCount=1
        for i in self.hits[1:]: #converts list type entries to strings
            if type(i)==list:
                self.hits[self.divCount]=str(i).strip('[]')[1:-1] #remove brackets and excess quotation marks
            self.divCount+=1
        self.divList=set(self.hits) #turns into set (remove duplicates)
        self.divList=list(self.divList) #turns back into list (scriptable)
        if '' in self.divList: #remove '' empty strings
            self.divList.remove('')
        self.divList.sort() #puts in numberical order
        self.divList_noContam=[]
        for i in self.divList:
            if 'sp|' in i:
                self.divList_noContam.append(i)
        print('The number of unique proteins in the list is: '+str(len(self.divList_noContam)))


a=''
a=input('What do you want to see? \na = pSer hits\nb = total hits\nEnter your choice: ')

if a=='a':
    b=ModSeq()
    b.removeContam()
    b.seqImport()
    b.ModTrim()
    b.SepFilter()
    
    determiner=''
    determiner=input('Do you want to see unique number of pSer hits that MAP to the encoded phosphosites (pSer at TAG)?\nType "Yes" or "No"\n')
    if determiner=='Yes':
        b.SepMap()
        print('The number of unique library-mapped pSer mode #1 phosphopeptides is: '+str(len(b.uniqueHits)))
    if determiner=='No':
        print('Ok, well bye then.')

if a=='b':
    b=ProtCountDiv()
    b.removeContam()
    b.seqImport()
    b.makeListDiv()
