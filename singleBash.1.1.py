#!/usr/bin/env python3
#IMPORTant stuff :p
#native imports that don't need installed
import os # operation system functions
import sys # system functions
import subprocess #subproccess systm package. similar to os functions
import time #native import capable of keeping track of time for me
import random as rand
import io

#start tracking the time the program takes to run
timerStart = time.perf_counter()

#try to import all packages, install if not
isPass = False
while isPass == False:
    try:
        #import bunch of programs needed
        import multiprocessing #multiprocess package
        import pandas as pd # gives excel like functions
        import seaborn as sns # pretty plot package
        from matplotlib import pyplot as plt # basic plotting package
        import glob
        import re
        import csv
        import gc
        import math
        import argparse

        # configure where plots are displayed
        #%matplotlib inline 
        import numpy as np # numerical functions
        import scipy.stats as stat #statistical engine
        import matplotlib.lines as mlines # for midlines
        sys.setrecursionlimit(3000)
        
        # exit the loop
        isPass = True
    #split up the error data and grab just the module name, then install that module
    #with mamba... hopefully
    except ModuleNotFoundError:
        exception_type, exception_object, exception_traceback = sys.exc_info()
        splitObjects = str(exception_object).split("'")
        
        #now try to install it
        os.system("mamba install "+splitObjects[1])
                  
#trouble shooting optional exit function
def stopProgram():
    repeat = True
    while repeat == True:
        contFunc = input("Continue running program? (y/n):\n")
        if contFunc.casefold() == "n":
            sys.exit("\nExiting...\n")
        elif contFunc.casefold() == "y":
            print("Continuing...\n")

parser=argparse.ArgumentParser()

parser.add_argument("--inputSeq", "-i", help="the FULL sequence you are looking at")
parser.add_argument("--bashSeq", "-b",help="sequence you want insert into your input DNA sequence and mutate your DNA with.")
parser.add_argument("--targetSeq", "-t",help="region you wish to bash in the form of sequence. ex: AACTTGGTA")
parser.add_argument("--zeroMotifs", "-z", help="Loops through null sequence until no motifs found in mutant region.",action="count",default=0)
parser.add_argument("--verbose","-v",help="increase the verbosity of outputs to get running updates and the offending motifs causing atempts to fail, along with location ID the motif start and genomic sequence found.", action="store_true")

args=parser.parse_args()
    
# pull from parser
inputSeq = args.inputSeq
bashSeqFull = args.bashSeq
targetSeq = args.targetSeq

#split the bashseq into two equal parts
bashSeq1, bashSeq2 = bashSeqFull[:len(bashSeqFull)//2], bashSeqFull[len(bashSeqFull)//2:]

#print running status for first MT
print('\nTrying for bash seq of first mutant type')
counter = 0
runStatus = True
while runStatus == True:
    if args.zeroMotifs > 0:
        ##use the null sequence
        nullSeq = bashSeq1
        bashSeq = nullSeq[counter:counter+len(targetSeq)]
        counter += 1
        
    #error out if the counter + len(targetSeq) is bigger than the nullSeq
    if counter+len(targetSeq) > len(nullSeq):
        print(f'first_MT_sequence:\nnotFound')
        newSeq1 = f'not found after {counter} iterations'
        print(newSeq1)
        runStatus = False
        continue


    if len(bashSeq) != len(targetSeq):
        # stop here because and error will ensue later
        print(f'len of {bashSeq} = {len(bashSeq)}\nlen of {targetSeq} = {len(targetSeq)}')
        print(f'length of bash sequence and target region un-equal')
        exit()

    bashReg = inputSeq.find(targetSeq)
    bashReg = bashReg - len(inputSeq)
    bashRegEnd = int(bashReg) + len(bashSeq)
    #print(f'{bashReg} to {bashRegEnd}')
    
    #make an output_Files directory to work into:
    pwd = os.getcwd()
    if not os.path.exists(f'{pwd}/output_Files'):
        os.system('mkdir output_Files')
    
    #store the user defined regions to a fasta file for analyzing with the other program
    wDir = os.getcwd()
    with open(wDir+'/seq_'+str(len(inputSeq))+'.fa', 'w') as targetFa:
        nameLine = '>seq_'+str(len(inputSeq))
        #print(f'legnth of original sequence: {len(inputSeq)}')
        targetFa.write(f'{nameLine}\n{inputSeq}')
        targetFa.close()
        
    #use the info above to find motifs in selected region
    inputFasta = wDir+'/seq_'+str(len(inputSeq))+'.fa'
    allMotifsFile = '/data/lab/duttke/software/HOMER/motifs/combined/homer_all.motifs'
    
    def findAllMotifs(targetDNA, motifGroup):
        #define a name for the output file
        targetName = targetDNA.split('.')[0]
        targetName = targetName.split('/')[-1]
        motifName = motifGroup.split('.')[0]
        motifName = motifName.split('/')[-1]
        name = targetName +'_'+ motifName
        os.system(f'findMotifs.pl {targetDNA} fasta output_Files/tmpFiles/ -find {motifGroup} > output_Files/{name}.out.txt 2>/dev/null')
        return(f'output_Files/{name}.out.txt')
    
    outNameWT = findAllMotifs(inputFasta, allMotifsFile)
                 
    # now that motifs have been found for the base sequnece
    # bash in the in the bashSeq and run motifs again
    # ref back through input seq at the point of bashReg to the point of end of bash seq
    newSeq = inputSeq.replace(targetSeq, bashSeq)
    #print(f'\nsequence removed:\n{targetSeq}')
    #print(f'\nsequence to insert:\n{bashSeq}')
    #print(f'\ninput sequence:\n{inputSeq}')
    #print(f'\nmutated sequence:\n{newSeq}\n')
    
    ## make a new fasta with it
    with open(wDir+'/seq_'+str(len(inputSeq))+'_mod.fa', 'w') as targetFa:
        linesToWrite = f'>{str(len(inputSeq))}_mod\n{newSeq}'
        #print(f'new sequence length of {len(newSeq)}')
        targetFa.write(linesToWrite)
        targetFa.close()
    
    # run for the new fasta file
    inputFasta = wDir+'/seq_'+str(len(inputSeq))+'_mod.fa'
    outNameMT = findAllMotifs(inputFasta, allMotifsFile)
    
    #now for the complex part. open the modified output file of motifs and find all motifs present based
    #on the area designated by user to find motifs present.
    with open(outNameMT, 'r') as outMotifs:
        #do the math for the "middle index" area
        lenInput = len(inputSeq)
        seqIndx = lenInput/2
        
        #start and end regions to search in
        miPosStart = (lenInput + bashReg) - seqIndx
        miPosEnd = (lenInput + bashRegEnd) - seqIndx
        
        #print(f'({lenInput} + {bashReg}) - {seqIndx} = miPosStart = {miPosStart}')
        #print(f'({lenInput} + {bashRegEnd}) - {seqIndx} = miPosEnd = {miPosEnd}')
            
        lineCount = 0
        motifsGained = []
        print(f'\nRep: {counter}') 
        for eachLine in outMotifs.readlines():
            if lineCount > 0:
                motifStrt = eachLine.split('\t')[1]
                motifStrt = int(motifStrt)
                
                # make a check to two types of index here. Forward strand (+) are index in a positive manner while 
                # reverse strand (-) are indexed backwards
                if eachLine.split('\t')[4] == '-':
                    motifEnd = motifStrt - len(eachLine.split('\t')[2])
                ## motif start is positive
                elif eachLine.split('\t')[4] == '+':
                    motifEnd = motifStrt + len(eachLine.split('\t')[2])
                ## motif start is 0 (and everthing the first two dont catch)
                else:
                    motifEnd = motifStrt + len(eachLine.split('\t')[2])
                
                qualCheck = eachLine.split('\t')[3]
                qualCheck = qualCheck.split('/')[0]
                #print(f'Checking motif {qualCheck}')
                
                #output the offending motif name and sequnece found in the desired region
                if ((motifStrt >= miPosStart) and (motifStrt <= miPosEnd)) or ((motifEnd >= miPosStart) and (motifEnd <= miPosEnd)):
                    offMotif = eachLine.split('\t')[1:4]
                    motifsGained.append(offMotif)
                    #print(f'start motif region index: {motifStrt}')
                    #print(f'endmotif region index: {motifEnd}')
                    print(f'Found offending motif of:\n{offMotif[-1].split("/")[0]} at {offMotif[0]}')
            lineCount += 1
        outMotifs.close()  
        
        
    # with both both output files present, open each and compare in pandas DF
    tarWT = pd.read_csv(outNameWT,sep='\t')
    tarMT = pd.read_csv(outNameMT,sep='\t')

    tarWTMotifs = []
    for index, row in tarWT.iterrows():
        smolLst = [row['Motif Name'],row['Sequence'],row['Offset'],row['Strand']]
        #print(smolLst)
        tarWTMotifs.append(smolLst)

    tarMTMotifs = []  
    for index, row in tarMT.iterrows():
        smolLst = [row['Motif Name'],row['Sequence'],row['Offset'],row['Strand']]
        #print(smolLst)
        tarMTMotifs.append(smolLst)
    
    ## simplified
    ######################################################### FINAL OUTPUTS #######################################################################################
    motifsLost = []
    for eachMotif in tarWTMotifs:
        if eachMotif not in tarMTMotifs:
            motifsLost.append(eachMotif)
    
    newMotifs = []
    for eachMotif in tarMTMotifs:
        if eachMotif not in tarWTMotifs:
            newMotifs.append(eachMotif)
   
    if args.zeroMotifs >= 1:
        ## check to see if no motifs show in region.
        ## if 0, continue on the thing below, if not, pass to top of loop ref 1 to the right
        if len(motifsGained) >= 1:
            continue
        
    print(f'first_MT_sequence:\n{newSeq}\nIteraions:{counter}')
    
    if len(motifsLost) > 0:
        print(f'\nMotifs lost from WT:')
        for eachMotif in motifsLost:
            print(str(eachMotif[2])+'('+eachMotif[1]+')'+eachMotif[0].split('(')[0])

    if len(newMotifs) > 0:
        print(f'\nMotifs gained in MT:')
        for eachMotif in newMotifs:
            print(str(eachMotif[2])+'('+eachMotif[1]+')'+eachMotif[0].split('(')[0])
    newSeq1 = newSeq
    #end the running loop
    runStatus = False

print('\nTrying for bash seq of second mutant type')
if args.zeroMotifs >= 1:
    #run the the second half null sequence to find a second unique 


    inputSeq = args.inputSeq
    bashSeqFull = args.bashSeq
    targetSeq = args.targetSeq

    counter = 0
    runStatus = True
    while runStatus == True:

        if args.zeroMotifs > 0:
            ##use the null sequence
            nullSeq = bashSeq2
            bashSeq = nullSeq[counter:counter+len(targetSeq)]
            counter += 1

        #error out if the counter + len(targetSeq) is bigger than the nullSeq
        if counter+len(targetSeq) > len(nullSeq):
            print(f'second_MT_sequence:\nnotFound')
            newSeq2 = f'not found after {counter} iterations'
            print(newSeq2)
            runStatus = False
            continue


        if len(bashSeq) != len(targetSeq):
            # stop here because and error will ensue later
            print(f'len of {bashSeq} = {len(bashSeq)}\nlen of {targetSeq} = {len(targetSeq)}')
            print(f'length of bash sequence and target region un-equal')
            exit()

        bashReg = inputSeq.find(targetSeq)
        bashReg = bashReg - len(inputSeq)
        bashRegEnd = int(bashReg) + len(bashSeq)
        #print(f'{bashReg} to {bashRegEnd}')

        #make an output_Files directory to work into:
        if not os.path.exists(f'{pwd}/output_Files'):
            os.system('mkdir output_Files')

        #store the user defined regions to a fasta file for analyzing with the other program
        wDir = os.getcwd()
        with open(wDir+'/seq_'+str(len(inputSeq))+'.fa', 'w') as targetFa:
            nameLine = '>seq_'+str(len(inputSeq))
            #print(f'legnth of original sequence: {len(inputSeq)}')
            targetFa.write(f'{nameLine}\n{inputSeq}')
            targetFa.close()

        #use the info above to find motifs in selected region
        inputFasta = wDir+'/seq_'+str(len(inputSeq))+'.fa'
        allMotifsFile = '/data/lab/duttke/software/HOMER/motifs/combined/homer_all.motifs'

        def findAllMotifs(targetDNA, motifGroup):
            #define a name for the output file
            targetName = targetDNA.split('.')[0]
            targetName = targetName.split('/')[-1]
            motifName = motifGroup.split('.')[0]
            motifName = motifName.split('/')[-1]
            name = targetName +'_'+ motifName
            os.system(f'findMotifs.pl {targetDNA} fasta output_Files/tmpFiles/ -find {motifGroup} > output_Files/{name}.out.txt 2>/dev/null' )
            return(f'output_Files/{name}.out.txt')
        
        #eat the stdout of this function so IT QUITS PRINTING JUNK
        text_trap = io.StringIO()
        sys.stdout = text_trap
        outNameWT = findAllMotifs(inputFasta, allMotifsFile)

        # now that motifs have been found for the base sequnece
        # bash in the in the bashSeq and run motifs again
        # ref back through input seq at the point of bashReg to the point of end of bash seq
        newSeq = inputSeq.replace(targetSeq, bashSeq)
        #print(f'\nsequence removed:\n{targetSeq}')
        #print(f'\nsequence to insert:\n{bashSeq}')
        #print(f'\ninput sequence:\n{inputSeq}')
        #print(f'\nmutated sequence:\n{newSeq}\n')

        ## make a new fasta with it
        with open(wDir+'/seq_'+str(len(inputSeq))+'_mod.fa', 'w') as targetFa:
            linesToWrite = f'>{str(len(inputSeq))}_mod\n{newSeq}'
            #print(f'new sequence length of {len(newSeq)}')
            targetFa.write(linesToWrite)
            targetFa.close()

        # run for the new fasta file
        inputFasta = wDir+'/seq_'+str(len(inputSeq))+'_mod.fa'
        outNameMT = findAllMotifs(inputFasta, allMotifsFile)
        
        #restore the stdout
        sys.stdout = sys.__stdout__
        
        #now for the complex part. open the modified output file of motifs and find all motifs present based
        #on the area designated by user to find motifs present.
        with open(outNameMT, 'r') as outMotifs:
            #do the math for the "middle index" area
            lenInput = len(inputSeq)
            seqIndx = lenInput/2

            #start and end regions to search in
            miPosStart = (lenInput + bashReg) - seqIndx
            miPosEnd = (lenInput + bashRegEnd) - seqIndx

            #print(f'({lenInput} + {bashReg}) - {seqIndx} = miPosStart = {miPosStart}')
            #print(f'({lenInput} + {bashRegEnd}) - {seqIndx} = miPosEnd = {miPosEnd}')

            lineCount = 0
            motifsGained = []
            print(f'\nRep: {counter}') 
            for eachLine in outMotifs.readlines():
                if lineCount > 0:
                    motifStrt = eachLine.split('\t')[1]
                    motifStrt = int(motifStrt)
                    
                    ## motif start is negative strand 
                    if eachLine.split('\t')[4] == '-':
                        motifEnd = motifStrt - len(eachLine.split('\t')[2])
                    ## motif start is positive strand
                    elif eachLine.split('\t')[4] == '+':
                        motifEnd = motifStrt + len(eachLine.split('\t')[2])
                    ## motif start is 0 (and everthing the first two dont catch)
                    else:
                        motifEnd = motifStrt + len(eachLine.split('\t')[2])

                    qualCheck = eachLine.split('\t')[3]
                    qualCheck = qualCheck.split('/')[0]
                    #print(f'Checking motif {qualCheck}')

                    #output the offending motif name and sequnece found in the desired region
                    if ((motifStrt >= miPosStart) and (motifStrt <= miPosEnd)) or ((motifEnd >= miPosStart) and (motifEnd <= miPosEnd)):
                        offMotif = eachLine.split('\t')[1:4]
                        motifsGained.append(offMotif)
                        #print(f'start motif region index: {motifStrt}')
                        #print(f'endmotif region index: {motifEnd}')
                        print(f'Found offending motif of:\n{offMotif[-1].split("/")[0]}  at {offMotif[0]}')                
                lineCount += 1
            outMotifs.close()


        # with both both output files present, open each and compare in pandas DF
        tarWT = pd.read_csv(outNameWT,sep='\t')
        tarMT = pd.read_csv(outNameMT,sep='\t')

        tarWTMotifs = []
        for index, row in tarWT.iterrows():
            smolLst = [row['Motif Name'],row['Sequence'],row['Offset'],row['Strand']]
            #print(smolLst)
            tarWTMotifs.append(smolLst)

        tarMTMotifs = []
        for index, row in tarMT.iterrows():
            smolLst = [row['Motif Name'],row['Sequence'],row['Offset'],row['Strand']]
            #print(smolLst)
            tarMTMotifs.append(smolLst)

        ## simplified
        ######################################################### FINAL OUTPUTS #######################################################################################
        motifsLost = []
        for eachMotif in tarWTMotifs:
            if eachMotif not in tarMTMotifs:
                motifsLost.append(eachMotif)

        newMotifs = []
        for eachMotif in tarMTMotifs:
            if eachMotif not in tarWTMotifs:
                newMotifs.append(eachMotif)

        if args.zeroMotifs >= 1:
            ## check to see if no motifs show in region.
            ## if 0, continue on the thing below, if not, pass to top of loop ref 1 to the right
            if len(motifsGained) >= 1:
                continue
        

        print(f'second_MT_sequence:\n{newSeq}\nIterations: {counter}')
        if len(motifsLost) > 0:
            print(f'\nMotifs lost from WT:')
            for eachMotif in motifsLost:
                print(str(eachMotif[2])+'('+eachMotif[1]+')'+eachMotif[0].split('(')[0])

        if len(newMotifs) > 0:
            print(f'\nMotifs gained in MT:')
            for eachMotif in newMotifs:
                print(str(eachMotif[2])+'('+eachMotif[1]+')'+eachMotif[0].split('(')[0])
        newSeq2 = newSeq
        #end the running loop
        runStatus = False

##save both the newSeq outputs to a file seperated by a new line
with open('mtSeq.txt', 'w') as mtSeq:
    mtSeq.write(f'WT_Sequence:{inputSeq}\nMT_Sequence_1:{newSeq1}\nMT_Sequence_2:{newSeq2}')
    mtSeq.close()



