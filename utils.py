#!/usr/bin/env python3.9
# -*- coding: utf-8 -*-
''' This module contains utilities for working with fastq files containing sequences and for cleaning purposes '''

#Necessary modules
from datetime import datetime
import statistics, time, re
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def hora():
    ''' This function returns the corresponding date and time '''
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    return("[" + dt_string + "] ")

def count_seq(datos):
    ''' This function counts the number of sequences '''
    n_seq = 0
    for record in datos:
        n_seq += 1
    return n_seq

def edit_seq(datos,limit_N,limit_Q,n_seq):
    ''' This function removes the sequences with number of Ns above a quantity and with quality below a value'''
    #We create the lists where the good sequences are stored
    valid_N = []
    valid = []
    
    for record in datos:
        N = 0
        media = statistics.mean(record.letter_annotations["phred_quality"])
        seq = str(record.seq)
        for s in seq:
            if s == "N":
                N += 1
            else:
                continue
        if N <= limit_N:
            valid_N.append(record)
            if media >= limit_Q:
                valid.append(record)
            
    #We show the number of remaining sequences
    rm_N = n_seq - len(valid_N)   
    rm_Q = len(valid_N) - len(valid)
    n_seq = n_seq - rm_N - rm_Q      
    print(hora() + str(rm_N) + " sequences with more than " + str(limit_N) + " Ns have been removed\n")
    time.sleep(4)    
    print(hora() + str(rm_Q) + " sequences with less than " + str(limit_Q) + " of quality have been removed\n")
    time.sleep(2)
    print(hora() + str(n_seq) + " sequences are kept\n")
    return valid

def analisis_organismos(valid,genomes_f,n_dir):
    ''' This function analyzes the resulting sequences and shows how many belong to each organism '''
    
    #List to save Human sequences
    Humano = []

    #We create a dictionary to store the number of sequences of each organism
    organismo = {
        "Human": 0,
        "Escherichia_coli": 0,
        "Helicobacter_pylori": 0,
        "Pseudomonas_aeruginosa": 0,
        "Sinorhizobium_meliloti": 0,
        "Staphylococcus_aureus": 0
        }

    #For each valid sequence, the Ns are changed to . to be able to search for them in each of the files
    #If it goes through all the files and is not found, it belongs to Human    
    for vsequence in valid:
        sequence = str(vsequence.seq)
        sequence = re.sub("N","(A|C|G|T)",sequence)
        n_file = 0
        found = False
        for file in genomes_f:
            n_file += 1
            for record_seq in SeqIO.parse(n_dir + "/" + file,"fasta"):
                if n_file < len(genomes_f) and found == False:
                    if re.search(sequence,str(record_seq.seq)):
                        organismo[record_seq.id] += 1
                        found = True   
                            
                elif n_file == len(genomes_f) and found == False:
                            organismo["Human"] += 1
                            Humano.append(vsequence)
          
    #We show the results on the screen
    for o in organismo:
        print("\n Found " + str(organismo[o]) + " records belonging to " + o)                   
        
    print("\n")
    
    #To calculate barplot ranges
    maximo1 = 0
    maximo2 = 0
    
    for o in organismo:
        if organismo[o] > maximo1:
            maximo1 = organismo[o]
        elif organismo[o] > maximo2:
            maximo2 = organismo[o]
    
    #We transform the dictionary into a pandas DataFrame
    for o in organismo:
        organismo[o] = [organismo[o]]
  
    df = pd.DataFrame.from_dict(organismo, orient = "index")
    df.columns=["Reads"]
    
    #We create the barplot 
    fig, (ax1,ax2) = plt.subplots(2,1,sharex=True,
                         figsize=(5,6))
    ax1.spines['bottom'].set_visible(False)
    ax1.tick_params(axis='x',which='both',bottom=False)
    ax2.spines['top'].set_visible(False)
    
    #Barplot ranges
    rango11 = maximo1 - 500
    rango12 = maximo1 +  500
    rango2 = maximo2 + 100
    ax2.set_ylim(0,rango2)
    ax1.set_ylim(rango11,rango12)
    ax1.set_yticks(np.arange(rango11,rango12,rango2))
    
    df.plot(ax=ax1,kind='bar')
    df.plot(ax=ax2,kind='bar', legend = False)
    for tick in ax2.get_xticklabels():
        tick.set_rotation(15)
    d = .015  
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (-d, +d), **kwargs)      
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    kwargs.update(transform=ax2.transAxes)  
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
    
    #We save it
    plt.savefig("organism_plot.png")

    #We export it to excel
    writer = pd.ExcelWriter("organism.xlsx")
    df.to_excel(writer)
    writer.save()
    
    print(hora() + "organism_plot.png and organism.xlsx files have been created\n")
    time.sleep(2)
    
    #We print the Human sequences in a fasta file
    SeqIO.write(Humano,"Human_seqs.fq","fastq")
    
    print(hora() + "The file Human_seqs.fq with " + str(len(Humano)) + " cleaned human sequences has been created\n")
    time.sleep(2)
    print(hora() + "Analysis finished\n")    
    return
