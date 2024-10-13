#!/usr/bin/env python


###
### VCF Filtering and Fasta Creation Version 0.1
###

"""
Version:        0.1.1
Author:         Toni de Dios Martinez
Contact:        tonidedios94@gmail.com
Date:           07/10/2024
Citation:       Biancamaria Bonucci, Toni de-Dios, Rémi Barbieri, Jess Emma Thompson, Sofia Panella, Francesca Radina, Sandra Savilli, Helja Kabral, Anu Solnik, Mary Anne Tafuri, John Robb, Christiana Lyn Scheib
                Ancient DNA and protein preservation and transfer in concretions covering human remains.
                iScience

Usage:          python vcfFFCv0.1.py -i <Sorted Multi Sample VCF> -o <Output file> -r <Heterozyogus Ratio> 
                                     -p <Min SNP distance in BP> -q <Min Calling Quality> -d <Min Depth> 
                                     -m <Min % of covered SNPs>

Example:        

Output:         
"""



import os, sys, argparse
import math
import re

### 
### Usage ADDMA_VCFfilter.py -i -r 
### Will asume that the VCF is properly sorted.
parser = argparse.ArgumentParser()

parser.add_argument( "-i", "--input",   dest = "vcf",     required=True,    help="Input VCF file"                 )
parser.add_argument( "-o", "--output",  dest = "out",     default="FFout",  help="Output files"                   )
parser.add_argument( "-r", "--ratio",   dest = "het",     default = 0.9,    help="Heterozygous ratio"             )
parser.add_argument( "-p", "--dist",    dest = "dst",     default = 3,      help="Minumum Distance between snps"  )
parser.add_argument( "-q", "--qual",    dest = "qual",    default = 1,      help="Calling quality"                )
parser.add_argument( "-d", "--depth",   dest = "depth",   default = 1,      help="Depth"                          )
parser.add_argument( "-m", "--missing", dest = "miss",    default = 0,      help="Missingness"                    )


args = parser.parse_args()

### Arguments

try:    ### Open VCF
    with open(args.vcf,"r") as vcf:
        vcf = open(args.vcf,"r")
except FileNotFoundError:
    print('The VCF file does not exist')
    quit()

try:    ### Het Ratio
    float(args.het)
    if 0.5 < float(args.het) <= 1:
        het=float(args.het)
    else:
        raise ValueError
except ValueError:
    print('Het ratio must be between 0.5 and 1')
    quit()
    
try:    ### Depth
    int(args.depth)
    if int(args.depth) >= 1:
        depth=int(args.depth)
    else:
        raise ValueError
except ValueError:
    print('Depth must be a positive integer bigger than 1')
    quit()

try:    ### Qual
    int(args.qual)
    if int(args.qual) >= 1:
        qual=int(args.qual)
    else:
        raise ValueError
except ValueError:
    print('Quality must be a positive integer bigger than 1')
    quit()

try:    ### SNP distance
    int(args.dst)
    if int(args.dst) >= 0:
        dst=int(args.dst)
    else:
        raise ValueError
except ValueError:
    print('SNP distance must be a positive integer')
    quit()

try:    ### Missingness
    float(args.miss)
    if 0 <= float(args.miss) <= 1:
        miss=float(args.miss)
    else:
        raise ValueError
except ValueError:
    print('Missing value should be between 0 (all missing) or 1 (no missing)')
    quit()

try:    ### Output
    str(args.out)
    output=str(args.out)

    fasta_out=".fasta"
    fasta_out=output+fasta_out

    vcf_out=".vcf"
    vcf_out=output+vcf_out
    
    stats_out=".stats.txt"
    stats_out=output+stats_out

    params_out=".params.txt"
    params_out=output+params_out
except ValueError:
    print('Something wrong happened lol.')
    quit()    

### Results Variables

newVCF=[]
sample=[]
fasta=[]
fastaPos=[]
heterozyosityPos=[]
heterozyosity=[]
missingPos=[]

### Functions

def check_ratio(x):
    y=[ int(j) for j in x ]
    return(min(y)/(sum(y)))

def check_ref(x):
    y=[ int(j) for j in x ]
    return(y[0]/sum(y))

def check_alt(x):
    y=[ int(j) for j in x ]
    return(y[1]/sum(y))

def add_2fasta(x):
    fasta.append(x)

def add_pos2fa(x):
    fastaPos.append(x)

def mean_het(x):
    try:
        return(sum(x)/len(x))
    except:
        return("NA")

def merge_result(x):
    x=x

def merge_geno(x):
    x=x

def trans_matrix(M):
    return [[M[x][y] for x in range(len(M))] for y in range(len(M[0]))]

def all_equal(iterator):
    iterator = iter(iterator)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == x for x in iterator)    

SNPposition=["NO_POT_SER_QUE_EL_CHROMOSOMA",int(-1-dst)]
SNPpositionNew=[0,0]
record=0

# num_lines = sum(1 for _ in vcf)

# print("Reading",vcf,"\n")

print("Opening File\n")

for line in vcf: ### Reading VCF
    li=line.strip()


    if li.startswith("#"):  
        newVCF.append(line.strip("\n")) ### Storing Header
        
        if "#CHROM" in line:
            header=(line.strip("\n")).split("\t")
            to_append=header[9:(len(header))]
            sample.append(to_append)    ### Storing Samples For Fasta Use

    if not li.startswith("#"):

        record=(record+1)
        print(record,' positions analysed ...', end='\r')


        field=line.split("\t")

        new_genotype=[] ### Where new filtered genotypes are gonna be stored.

        chr=field[0]
        pos=field[1]
        sID=field[2]
        ref=field[3]
        alt=field[4]
        cQL=field[5]
        fil=field[6]
        INF=field[7]
        fmt=field[8]
        fastaPos=[]
        heterozyosityPos=[]

        if (len(alt.split(",")))>1:
            fil="Multiallelic"
            
            for sampleIDX in range(9,len(field)):
                geno = field[sampleIDX].split(":")
                GT=geno[0]
                AD=geno[1]
                try:
                    DP=int(geno[2])
                except:
                    continue    
                GT="./."
                NucleotideAll="N"
            continue
        if (len(ref) > 1) or (len(alt) > 1):
            fil="Indel"
            
            for sampleIDX in range(9,len(field)):
                geno = field[sampleIDX].split(":")
                GT=geno[0]
                AD=geno[1]
                try:
                    DP=int(geno[2])
                except:
                    continue    
                GT="./."
                NucleotideAll="N"
            continue
        if float(cQL) < float(qual):
            fil="Low_Quality"
            
            for sampleIDX in range(9,len(field)):
                geno = field[sampleIDX].split(":")
                GT=geno[0]
                AD=geno[1]
                try:
                    DP=int(geno[2])
                except:
                    continue    
                GT="./."
                NucleotideAll="N"
            continue
        if (alt == '.'):
            fil="Monoallelic"
           
            for sampleIDX in range(9,len(field)):
                geno = field[sampleIDX].split(":")
                GT=geno[0]
                AD=geno[1]
                try:
                    DP=int(geno[2])
                except:
                    GT="./."
                    NucleotideAll="N"
                    continue
                if (DP < depth):
                    GT="./."
                    NucleotideAll="N"
                    continue
                else:
                    NucleotideAll=str(ref)
                    continue
            continue
        else:

            for sampleIDX in range(9,len(field)):

                geno = field[sampleIDX].split(":")

                GT=geno[0]
                AD=geno[1]

                try:
                    DP=int(geno[2])
                except:
                    heterozyosityPos.append("NA")
                    Nucleotide="N"
                    NucleotideAll=Nucleotide
                    fastaPos.append(Nucleotide)
                    continue 

                                                         

                alleles=AD.split(",") ### 0/1:3,10:13 > 3 and 10
                

                if (len(alleles) < 2) or (DP < depth) or (GT == "./."):
                    GT="./."
                    Nucleotide="N"
                    NucleotideAll=Nucleotide
                    fastaPos.append(Nucleotide)
                    heterozyosityPos.append("NA")
                    continue

                if alleles[1] == ".":
                    GT="0/0"
                    Nucleotide=str(ref)
                    NucleotideAll=Nucleotide
                    fastaPos.append(Nucleotide)
                   
                    heterozyosityPos.append(0)
                    continue

                else:
                    try:
                        ratioHere=check_ratio(alleles)
                        heterozyosityPos.append(ratioHere)
                    except:
                        heterozyosityPos.append(0)

                    if check_ref(alleles) >= het:
                        GT="0/0"
                        Nucleotide=str(ref)
                        NucleotideAll=Nucleotide
                        fastaPos.append(Nucleotide)
                        continue
                    if check_alt(alleles) >= het:
                        GT="1/1"
                        Nucleotide=str(alt)
                        NucleotideAll=Nucleotide
                        fastaPos.append(Nucleotide)
                        continue
                    if check_alt(alleles) < het and check_ref(alleles) < het:
                        GT="./."
                        Nucleotide="N"
                        NucleotideAll=Nucleotide
                        fastaPos.append(Nucleotide)
                        continue
                    else:
                        GT="./."
                        Nucleotide="N"
                        NucleotideAll=Nucleotide
                        fastaPos.append(Nucleotide)
                        continue

            
            if (((len(fastaPos))-(fastaPos.count("N")))/len(fastaPos)) >= miss:
                check=[x for x in fastaPos if x != 'N']
                if all_equal(check) != True:
                    SNPpositionNew=[chr,pos]
                    if (SNPpositionNew[0]) != (SNPposition[0]):
                        SNPposition=[chr,pos]
                        add_2fasta(fastaPos)
                        heterozyosity.append(heterozyosityPos)
                        #print(chr,pos,fastaPos,len(fastaPos))
                        continue
                    else:
                        if (int(SNPpositionNew[1])) > ((int(SNPposition[1]))+int(dst)):
                            SNPposition=[chr,pos]
                            add_2fasta(fastaPos)
                            heterozyosity.append(heterozyosityPos)
                            continue

                            #print(chr,pos,fastaPos,len(fastaPos))
                        else:
                            continue
                else:
                    continue
            else:
                continue

print(record,' positions analysed ...\n')        
print("Filtering Finished!\n")


##### Het

print("Calculating Heterozygous ratio ...\n")


fasta=trans_matrix(fasta)
fasta=[ "".join(x) for x in fasta ]


heterozyosity=trans_matrix(heterozyosity)
heterozyosity=[[x for x in a if x != 'NA'] for a in heterozyosity]
heterozyosity=[mean_het(x) for x in heterozyosity]

sample=[x for y in sample for x in y]

#print(len(heterozyosity),len(fasta),len(sample))

final_stats=[]
final_flvcf=[]
final_fasta=[]

print("Formating Output ...\n")


### Output

for i in range(0,len(heterozyosity)):
    
    result_stats=[str(sample[i]),"\t",str((heterozyosity[i])*100),"\t",str(100-((str(fasta[i])).count("N")/len(str(fasta[i]))*100)),"\t",str(len(str(fasta[i])))]
    result_fasta=[">",str(sample[i]),"\n",str(fasta[i])]
    
    result_stats="".join(result_stats)
    final_stats.append(result_stats)

    result_fasta="".join(result_fasta)
    final_fasta.append(result_fasta)

    result_flvcf=""

    #print("".join(result))
    #head=[">", head1, "|",head2 ]
    #head="".join(head)
    #print("\n".join(head,"".join(fasta[i])))

header="Sample\tHetRatio\tBreadthCoverage\tPositions"

final_stats.insert(0,header)

final_stats.append("")

final_stats="\n".join(final_stats)

final_fasta.append("")
final_fasta="\n".join(final_fasta)


print("Printing Output ...\n")

f1=open(fasta_out, "w")
f1.write(final_fasta)
f1.close()

f2=open(stats_out, "w")
f2.write(final_stats)
f2.close()


print("Finished!\n")

#heterozyosity=trans_matrix(heterozyosity)


