#!/usr/bin/env python3.9
# -*- coding: utf-8 -*-

#necessary modules
from utils import hora, count_seq, edit_seq, analisis_organismos
import argparse, sys, os, time
import colorama as cl
from Bio import SeqIO

#First we ask for the parameters we need
parser = argparse.ArgumentParser(
    description = "This program accepts a fastq file to work based on a directory with genome files to compare and some quality parameters ")

#We ask for the file
parser.add_argument(
                    '-file',
                    action="store",
                    dest="file_in",
                    type=str,
                    help="Input file")
#We ask for the directory
parser.add_argument(
                    '-dir',
                    action="store",
                    dest="dir",
                    type=str,
                    help="Genome Directory")
#We ask for the Ns cut number
parser.add_argument(
                    '-N',
                    action="store",
                    dest="N_lim",
                    type=int,
                    help="Ns limit")
#We ask for the quality limit
parser.add_argument(
                    '-Q',
                    action="store",
                    dest="Q_value",
                    type=int,
                    help="Minimum quality of the sequences")
#We save the arguments
results = parser.parse_args()

#If the arguments are not inserted, an error will be raised for each argument and, in addition, it will print the help
if results.file_in == None:
    print(cl.Fore.GREEN + "\nPlease enter the fastq file to analyse\n" + cl.Fore.RESET)
    parser.print_help()
elif results.dir == None:
    print(cl.Fore.GREEN + "\nPlease enter the directory with the genomes to compare\n" + cl.Fore.RESET)
    parser.print_help()
elif results.N_lim == None or results.N_lim < 0:
    print(cl.Fore.GREEN + "\nPlease enter the limit of Ns of the sequences\n" + cl.Fore.RESET)
    parser.print_help()
elif results.Q_value == None or results.Q_value < 0 or results.Q_value > 40:
    print(cl.Fore.GREEN + "\nPlease enter the quality limit of the sequences (0-40)\n" + cl.Fore.RESET)
    parser.print_help()
else:
    #We enter here if the parameters are entered correctly
    
    #We separate all the console arguments with a space and display it on the screen
    print("\n" + hora() + "Entered data: " + " ".join(sys.argv) + "\n")
       
    #We check the input file
    (filename, extension) = os.path.splitext(results.file_in)
    if extension != ".fq":
        print(cl.Fore.GREEN + "\nPlease enter a valid fastq file\n" + cl.Fore.RESET)
    else:
        #We check that the files in the genome directory are valid
        genomes = os.listdir(results.dir)
        genomes_ok = True
        for file in genomes:
            (filename_q, extension_q) = os.path.splitext(file)
            if extension_q != ".fna":
                genomes_ok = False
        #If there is any wrong file, we show an error        
        if genomes_ok != True:
            print(cl.Fore.GREEN + "\nPlease check that the genome files are valid\n" + cl.Fore.RESET)
        else:      
            print(hora() + "Reading the file: " + results.file_in + "\n")
            time.sleep(2)

            #We show the number of sequences
            secuencias = count_seq(SeqIO.parse(results.file_in,"fastq"))
            print(hora() + "Found " + str(secuencias) + " sequences to be analysed\n")

            time.sleep(2)
            print(hora() + "Deleting sequences with more than " + str(results.N_lim) + " Ns and less of" + str(results.Q_value) + " of quality\n")
            time.sleep(4)
            #We create the list of good sequences
            valid_seq = []           
            #We filter by number of Ns and quality and obtain the final result
            valid_seq = edit_seq(SeqIO.parse(results.file_in,"fastq"),results.N_lim,results.Q_value,secuencias)
     
            print(hora() + "Analyzing contaminants\n")   

            #We analyze to which genome each sequence corresponds
            analisis_organismos(valid_seq, genomes, results.dir)

            
