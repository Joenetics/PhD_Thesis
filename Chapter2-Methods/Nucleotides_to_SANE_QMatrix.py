## make clustering based on % split. (e.g; 1% = 100 individual splits  and belonging)
#SANE: Simulating Ancestry through Nucleotide distance Equations
import sys
import os
import math
import datetime
import random
import numpy as np
import multiprocessing
from multiprocessing import process
import json
import statistics

# various options listed below...


#Psiko_BinaryInputFile_Plates_1-2-3-4-5-6-7SacchOnlyNoMaskedDoubleReduced_CIGAR.geno
#full_data = True  # put StDev, etc, into output Q file.
full_data = False  # do NOT put StDev etc into output Q file.
#IncludingMasked = "NoMasked"  # dont put masked SNPs in final output
IncludingMasked = ""  # put masked in final output
#OnlyS28C = "S28C_only"  # only use S28C genes. Build Q-matrix with this!
OnlyS28C = ""  # for everything else.
CoreOnly = "CoreOnly"
CoreOnly = ""
cigar_only = "_CIGAR"
#cigar_only = ""
SC = ""
#SC = "SacchOnly"
if SC == "SacchOnly":
    #SCset = "Full"
    #SCset = "Reduced"
    SCset = "DoubleReduced"
else:
    SCset = ""
NumberOfCores = 3
#branch_number = input("How many branches/columns you want in final q-matrix?")
branch_number = 3
global which_calculation
#which_calculation = input("Which calculation do you want to do? Percentage/TND/JCD/K2P/TamD")
#which_calculation = "Percentage"
#which_calculation = "TND"
#which_calculation = "JCD"
which_calculation = "K2P"
#which_calculation = "TamD"


print(f"{datetime.datetime.now()}: Started... ")
percent_accuracy = 5  # this can be 5,...
#percent_accuracy = input("Choose percent accuracy in %. e.g, 1,2,5,10...")  # this can be 5,...
# percent_accuracy = int(percent_accuracy)  # uncomment with previous

chosenfile= "Joined_MAF_files\\" + CoreOnly + "ATGC_File_AllStrains_Plates_1-2-3-4-5-6-7" + SC + IncludingMasked + OnlyS28C + SCset + cigar_only + ".csv"  # sacch strains and only S28C genes.
#chosenfile= "PlateDEVCF_CIGAR_AfterR\\" + CoreOnly + "ATGC_File_AllStrains_PlateDE" + SC + IncludingMasked + OnlyS28C + SCset + cigar_only + ".csv"  # sacch strains and only S28C genes.
#DE_yes = "DE"
DE_yes = ""
data_file = open(chosenfile, "r")
#output_for_tree = open(chosenfile.split(".")[0] + "_ATGCTree" + str(percent_accuracy) + ".q", "w")
output_for_tree = open(chosenfile.split(".")[0] + "_ATGCTree_" + which_calculation + str(percent_accuracy) + ".q", "w")
#file_for_Q_clusters = open(chosenfile.split(".")[0] + "_ATGCClusters" + str(percent_accuracy) + ".q", "w")
file_for_Q_clusters = open(chosenfile.split(".")[0] + "_ATGCClusters_" + which_calculation + str(percent_accuracy) + ".q", "w")



def estimate_nucleotide_frequencies(seq):
    seq = seq.replace('-','').upper()
    A = seq.count('A')
    C = seq.count('C')
    G = seq.count('G')
    T = seq.count('T')
    length = float(len(seq))
    return [ x/length for x in [A,C,G,T] ]
from math import log

def pdistance(seq1, seq2):
    p = 0
    pairs = []
    for x in zip(seq1,seq2):
        if '-' not in x:
            pairs.append(x)
    #for (x,y) in zip(seq1,seq2):
    for (x, y) in pairs:
        if x != y:
            p += 1
    #length = (len(seq1) + len(seq2)) / 2
    length = len(pairs)
    if (float(p) / length) >= 0.75:  # if more than 75% dissimilar, then maths fails. Do this to avoid errors. And big score will fail match anyways
        return 0.749
    else:
        return float(p) / length


def Jukes_Cantor_distance(seq1, seq2):
    """
    distance = -b log(1 - p / b)
    where:
    b = 3/4
    and p = p-distance, i.e. uncorrected distance between seq1 and seq2
    """
    from math import log
    b = 0.75
    p = pdistance(seq1, seq2)
    try:
        d = -b * log(1-p/b)
    except ValueError:
        print("Tried to take log of a negative number")

        return None
    return d


def TNdistance(seq1, seq2):
    """
    Tajima-Nei distance = -b log(1 - p / b)
    where:
    b = 0.5 * [ 1 - Sum i from A to T(Gi^2+p^2/h) ]
    h = Sum i from A to G( Sum j from C to T (Xij^2/2*Gi*Gj))
    p = p-distance, i.e. uncorrected distance between seq1 and seq2
    Xij = frequency of pair (i,j) in seq1 and seq2, with gaps removed
    Gi = frequency of base i over seq1 and seq2 """
    from math import log

    ns = ['A','C','G','T']
    G = estimate_nucleotide_frequencies(seq1 + seq2)
    p = pdistance(seq1, seq2)
    pairs = []
    h = 0

    #collect ungapped pairs
    for x in zip(seq1,seq2):
        if '-' not in x:
            pairs.append(x)

    #pair frequencies are calculated for AC, AG, AT, CG, CT, GT (and reverse order)
    for i in range(len(ns)-1):
        for j in range(i+1,len(ns)):
            if i != j:
                paircount = pairs.count( (ns[i], ns[j]) ) + pairs.count( (ns[j], ns[i]) )
                Xij_sq = (float(paircount)/len(pairs))**2
                GiGj = G[i]*G[j]
                h += 0.5*(Xij_sq/GiGj)  #h value used to calculate b
    b = 0.5*(1-sum([x**2 for x in G])+p**2/h)
    try:
        d = -b * log(1 - p/b)
    except ValueError:
        return -5
        print("Tried to take log of a negative number")
        # return None  # old
    return d

def K2Pdistance(seq1,seq2):
    """
    Kimura 2-Parameter distance = -0.5 log( (1 - 2p -q) * sqrt( 1 - 2q ) )
    where:
    p = transition frequency
    q = transversion frequency
    """
    from math import log, sqrt
    pairs = []

    #collect ungapped pairs
    for x in zip(seq1,seq2):
        if '-' not in x:
            pairs.append(x)

    ts_count=0
    tv_count=0
    length = len(pairs)

    transitions = [ "AG", "GA", "CT", "TC"]
    transversions = [ "AC", "CA", "AT", "TA",
                      "GC", "CG", "GT", "TG" ]

    for (x,y) in pairs:
        if x+y in transitions:
            ts_count += 1
        elif x+y in transversions:
            tv_count += 1

    p = float(ts_count) / length
    q = float(tv_count) / length
    try:
        d = -0.5 * log( (1 - 2*p - q) * sqrt( 1 - 2*q ) )
    except ValueError:
        return -5
        print("Tried to take log of a negative number")
        return None
    return d

def Tamuradistance(seq1,seq2):
    """
    Tamura distance = -C log( 1 - P/C - Q ) - 0.5( 1 - C )log( 1 - 2Q )
    where:
    P = transition frequency
    Q = transversion frequency
    C = GC1 + GC2 - 2 * GC1 * GC2
    GC1 = GC-content of sequence 1
    GC2 = GC-coontent of sequence 2
    """
    from math import log
    pairs = []

    #collect ungapped pairs
    for x in zip(seq1,seq2):
        if '-' not in x:
            pairs.append(x)

    ts_count=0
    tv_count=0
    length = len(pairs)

    transitions = [ "AG", "GA", "CT", "TC"]
    transversions = [ "AC", "CA", "AT", "TA",
                      "GC", "CG", "GT", "TG" ]

    for (x,y) in pairs:
        if x+y in transitions:
            ts_count += 1
        elif x+y in transversions:
            tv_count += 1

    p = float(ts_count) / length
    q = float(tv_count) / length
    gc1 = sum(estimate_nucleotide_frequencies(seq1)[1:3])
    gc2 = sum(estimate_nucleotide_frequencies(seq2)[1:3])
    c = gc1 + gc2 - 2 * gc1 * gc2

    try:
        d = -c * log( 1 - p/c - q) - 0.5 * ( 1 - c ) * log ( 1 - 2*q )
    except ValueError:
        return -5
        print("Tried to take log of a negative number")
    return d


def average_of_strains(all_genomes):
    blank_genome_counterA = [0] * len(all_genomes[0])
    blank_genome_counterT = [0] * len(all_genomes[0])
    blank_genome_counterG = [0] * len(all_genomes[0])
    blank_genome_counterC = [0] * len(all_genomes[0])
    averaged_genome = [0] * len(all_genomes[0])
    for genome in range(0, len(all_genomes)):
        snp_list = all_genomes[genome]
        for snp in range(0, len(snp_list)):
            if snp_list[snp] == "A":
                blank_genome_counterA[snp] += 1
            elif snp_list[snp] == "T":
                blank_genome_counterT[snp] += 1
            elif snp_list[snp] == "G":
                blank_genome_counterG[snp] += 1
            elif snp_list[snp] == "C":
                blank_genome_counterC[snp] += 1
    for snp in range(0, len(blank_genome_counterA)):
        ## assign probability to nucleotide based on frequency in population.
        ## if 10% frequency in population, 10% of being ancestral here.
        threshold_for_A = float(blank_genome_counterA[snp]/len(all_genomes))
        threshold_for_T = float(blank_genome_counterT[snp]/len(all_genomes)) + threshold_for_A
        threshold_for_G = float(blank_genome_counterG[snp]/len(all_genomes)) + threshold_for_T
        threshold_for_C = float(blank_genome_counterC[snp]/len(all_genomes))
        current_score = random.random()
        if current_score <= threshold_for_A:
            averaged_genome[snp] = "A"
        elif threshold_for_A < current_score <= threshold_for_T:
            averaged_genome[snp] = "T"
        elif threshold_for_T < current_score <= threshold_for_G:
            averaged_genome[snp] = "G"
        elif current_score > threshold_for_G:
            averaged_genome[snp] = "C"




    averaged_genome = " ".join(averaged_genome)
    return averaged_genome


def next_branch(current_branches_dictionary, all_data_dictionary):
    current_branches = []
    for eachbranch in current_branches_dictionary:
        current_branches.append(eachbranch)
    adding_scores_dictionary = {}      # This dictionary will store sum of distances of each potential branch from all current branches
    for branch in current_branches:
        for potential_branch in all_data_dictionary:
            if potential_branch in current_branches:
                pass  # dont test branch against itself!
            else:
                if which_calculation == "Percentage":
                    if (branch + "_" + potential_branch) in every_strain_matched_dictionary:
                        if potential_branch in adding_scores_dictionary:
                            adding_scores_dictionary[potential_branch] = float(adding_scores_dictionary[potential_branch] + every_strain_matched_dictionary[(branch + "_" + potential_branch)])
                        else:
                            adding_scores_dictionary[potential_branch] = every_strain_matched_dictionary[(branch + "_" + potential_branch)]

                    else:
                        print("DD")
                        if potential_branch in adding_scores_dictionary:
                            adding_scores_dictionary[potential_branch] = float(adding_scores_dictionary[potential_branch] + pdistance(all_data_dictionary[branch], all_data_dictionary[potential_branch]))
                        else:
                            adding_scores_dictionary[potential_branch] = pdistance(all_data_dictionary[branch], all_data_dictionary[potential_branch])

                elif which_calculation == "JCD":
                    if (branch + "_" + potential_branch) in every_strain_matched_dictionary:
                        if potential_branch in adding_scores_dictionary:
                            adding_scores_dictionary[potential_branch] = float(adding_scores_dictionary[potential_branch] + every_strain_matched_dictionary[(branch + "_" + potential_branch)])
                        else:
                            adding_scores_dictionary[potential_branch] = every_strain_matched_dictionary[(branch + "_" + potential_branch)]

                    else:
                        print("DD")
                        if potential_branch in adding_scores_dictionary:
                            adding_scores_dictionary[potential_branch] = float(adding_scores_dictionary[potential_branch] + Jukes_Cantor_distance(all_data_dictionary[branch], all_data_dictionary[potential_branch]))
                        else:
                            adding_scores_dictionary[potential_branch] = Jukes_Cantor_distance(all_data_dictionary[branch], all_data_dictionary[potential_branch])
                elif which_calculation == "TND":
                    if (branch + "_" + potential_branch) in every_strain_matched_dictionary:
                        if potential_branch in adding_scores_dictionary:
                            adding_scores_dictionary[potential_branch] = float(adding_scores_dictionary[potential_branch] + every_strain_matched_dictionary[(branch + "_" + potential_branch)])
                        else:
                            adding_scores_dictionary[potential_branch] = every_strain_matched_dictionary[(branch + "_" + potential_branch)]

                    else:
                        print("DD")
                        if potential_branch in adding_scores_dictionary:
                            adding_scores_dictionary[potential_branch] = float(adding_scores_dictionary[potential_branch] + TNdistance(all_data_dictionary[branch], all_data_dictionary[potential_branch]))
                        else:
                            adding_scores_dictionary[potential_branch] = TNdistance(all_data_dictionary[branch], all_data_dictionary[potential_branch])
                elif which_calculation == "K2P":
                    if (branch + "_" + potential_branch) in every_strain_matched_dictionary:
                        if potential_branch in adding_scores_dictionary:
                            adding_scores_dictionary[potential_branch] = float(adding_scores_dictionary[potential_branch] + every_strain_matched_dictionary[(branch + "_" + potential_branch)])
                        else:
                            adding_scores_dictionary[potential_branch] = every_strain_matched_dictionary[(branch + "_" + potential_branch)]

                    else:
                        print("DD")
                        if potential_branch in adding_scores_dictionary:
                            adding_scores_dictionary[potential_branch] = float(adding_scores_dictionary[potential_branch] + K2Pdistance(all_data_dictionary[branch], all_data_dictionary[potential_branch]))
                        else:
                            adding_scores_dictionary[potential_branch] = K2Pdistance(all_data_dictionary[branch], all_data_dictionary[potential_branch])
                elif which_calculation == "TamD":
                    if (branch + "_" + potential_branch) in every_strain_matched_dictionary:
                        if potential_branch in adding_scores_dictionary:
                            adding_scores_dictionary[potential_branch] = float(adding_scores_dictionary[potential_branch] + every_strain_matched_dictionary[(branch + "_" + potential_branch)])
                        else:
                            adding_scores_dictionary[potential_branch] = every_strain_matched_dictionary[(branch + "_" + potential_branch)]

                    else:
                        print("DD")
                        if potential_branch in adding_scores_dictionary:
                            adding_scores_dictionary[potential_branch] = float(adding_scores_dictionary[potential_branch] + Tamuradistance(all_data_dictionary[branch], all_data_dictionary[potential_branch]))
                        else:
                            adding_scores_dictionary[potential_branch] = Tamuradistance(all_data_dictionary[branch], all_data_dictionary[potential_branch])


    best_branch = max(adding_scores_dictionary, key=adding_scores_dictionary.get)  # highest score = highest nucleotide difference. Find diverse strains to determine cluster origins
    return best_branch


def adding_sub_strain(current_branches_dictionary, all_data_dictionary):
    current_branches = []
    for eachbranch in current_branches_dictionary:
        current_branches.append(eachbranch)
    adding_scores_dictionary = {}      # This dictionary will store sum of distances of each potential branch from all current branches
    eachbranch_dictioary = {}
    for branch in current_branches:
        for potential_branch in all_data_dictionary:
            if potential_branch in current_branches:
                pass  # dont test branch against itself!
            else:
                if which_calculation == "Percentage":
                    #if (branch + "_" + potential_branch) in every_strain_matched_dictionary:
                   adding_scores_dictionary[potential_branch] = every_strain_matched_dictionary[branch + "_" + potential_branch]

                elif which_calculation == "JCD":
                    #if (branch + "_" + potential_branch) in every_strain_matched_dictionary:
                    adding_scores_dictionary[potential_branch] = every_strain_matched_dictionary[branch + "_" + potential_branch]
                elif which_calculation == "TND":
                    #if (branch + "_" + potential_branch) in every_strain_matched_dictionary:
                    adding_scores_dictionary[potential_branch] = every_strain_matched_dictionary[branch + "_" + potential_branch]
                elif which_calculation == "K2P":
                    #if (branch + "_" + potential_branch) in every_strain_matched_dictionary:
                    adding_scores_dictionary[potential_branch] = every_strain_matched_dictionary[branch + "_" + potential_branch]
                elif which_calculation == "TamD":
                    #if (branch + "_" + potential_branch) in every_strain_matched_dictionary:
                    adding_scores_dictionary[potential_branch] = every_strain_matched_dictionary[branch + "_" + potential_branch]
        # adding_scores_dictionary has a branch of main tree matched to find what branch to add next. lowest score ==most similarity. store.
        eachbranch_dictioary[branch + "-" + min(adding_scores_dictionary, key=adding_scores_dictionary.get)] = adding_scores_dictionary[min(adding_scores_dictionary, key=adding_scores_dictionary.get)]  # parentbranch-childbranch = score
        adding_scores_dictionary = {}
    best_branch = str(min(eachbranch_dictioary, key=eachbranch_dictioary.get)) + "-" + str(eachbranch_dictioary[min(eachbranch_dictioary, key=eachbranch_dictioary.get)])  # parent_branch-subbranch-score  # only add new strain where it fits best on present tree (parent == on current tree, child == potential strain addition)
    return best_branch


def first_branch(all_data_dictionary):  # Best starting branch/trunk will be strain with lowest similarity to another. Meaning, diverse strains already identified
    adding_scores_dictionary = {}      # This dictionary will store sum of distances of each potential branch from all current branches
    best_strain_score = 0  # initialise it. every strain will match more than 0%, so should be fine.
    for branch1 in all_data_dictionary:
        for branch2 in all_data_dictionary:
            if branch1 == branch2:
                pass
            else:
                if which_calculation == "Percentage":
                    if (branch1 + "_" + branch2) in every_strain_matched_dictionary:
                        adding_scores_dictionary[branch2] = every_strain_matched_dictionary[(branch1 + "_" + branch2)]
                    else:
                        print("Wrong JSON was read in")
                        adding_scores_dictionary[branch2] = pdistance(all_data_dictionary[branch1], all_data_dictionary[branch2])

                elif which_calculation == "JCD":
                    if (branch1 + "_" + branch2) in every_strain_matched_dictionary:
                        adding_scores_dictionary[branch2] = every_strain_matched_dictionary[(branch1 + "_" + branch2)]
                    else:
                        print("Wrong JSON read in")
                        adding_scores_dictionary[branch2] = Jukes_Cantor_distance(all_data_dictionary[branch1], all_data_dictionary[branch2])
                elif which_calculation == "TND":
                    if (branch1 + "_" + branch2) in every_strain_matched_dictionary:
                        adding_scores_dictionary[branch2] = every_strain_matched_dictionary[(branch1 + "_" + branch2)]
                    else:
                        print("Wrong JSON read in")
                        adding_scores_dictionary[branch2] = TNdistance(all_data_dictionary[branch1], all_data_dictionary[branch2])
                elif which_calculation == "K2P":
                    if (branch1 + "_" + branch2) in every_strain_matched_dictionary:
                        adding_scores_dictionary[branch2] = every_strain_matched_dictionary[(branch1 + "_" + branch2)]
                    else:
                        print("Wrong JSON read in")
                        adding_scores_dictionary[branch2] = K2Pdistance(all_data_dictionary[branch1], all_data_dictionary[branch2])
                elif which_calculation == "TamD":
                    if (branch1 + "_" + branch2) in every_strain_matched_dictionary:
                        adding_scores_dictionary[branch2] = every_strain_matched_dictionary[(branch1 + "_" + branch2)]
                    else:
                        print("Wrong JSON read in")
                        adding_scores_dictionary[branch2] = Tamuradistance(all_data_dictionary[branch1], all_data_dictionary[branch2])

        current_strain_high_score_key = max(adding_scores_dictionary, key=adding_scores_dictionary.get)
        current_strain_high_score = adding_scores_dictionary[current_strain_high_score_key]
        if float(current_strain_high_score) >= float(best_strain_score):
            best_strain = max(adding_scores_dictionary, key=adding_scores_dictionary.get)
            best_strain_score = adding_scores_dictionary[best_strain]
            best_strain = branch1

        adding_scores_dictionary = {}
    return best_strain

print(f"{datetime.datetime.now()}: All functions loaded... ")
strains_dictionary = {}
empty_score_dictionary = {}
print(f"{datetime.datetime.now()}: Reading in strains from file... ")
for line in data_file:
    data = line.strip().split(",")
    strains_dictionary[data[0]] = ''.join(data[1:])
    empty_score_dictionary[data[0]] = float(0)

print(f"{datetime.datetime.now()}: Calculating genetic distances between strains... ")

# Save
#np.save(chosenfile.split("\\")[1].split(".")[0]+ ".npy", every_strain_matched_dictionary)

# Load
global every_strain_matched_dictionary  # massive dictionary of matching all strains to eachother
every_strain_matched_dictionary = {}

try:
    f = open("JSON_Files\\"+ CoreOnly + SC + IncludingMasked + OnlyS28C + SCset + cigar_only + which_calculation + DE_yes + '.json')
    every_strain_matched_dictionary = json.load(f)
    f.close()
    print(f"{datetime.datetime.now()}: Main dictionary loaded from previous save... ")
except:
    for parent_branch in strains_dictionary:
            for child_branch in strains_dictionary:
                if parent_branch == child_branch:
                    pass
                else:
                    if which_calculation == "Percentage":
                        every_strain_matched_dictionary[parent_branch + "_" + child_branch] = pdistance(strains_dictionary[parent_branch], strains_dictionary[child_branch])
                    elif which_calculation == "JCD":
                        try:
                            every_strain_matched_dictionary[parent_branch + "_" + child_branch] = Jukes_Cantor_distance(strains_dictionary[parent_branch], strains_dictionary[child_branch])
                        except:
                            print("somethign failed. " +  strains_dictionary[parent_branch] + " " + strains_dictionary[child_branch])
                            sys.exit()
                            every_strain_matched_dictionary[parent_branch + "_" + child_branch] = 0.00001
                    elif which_calculation == "TND":
                        try:
                            every_strain_matched_dictionary[parent_branch + "_" + child_branch] = TNdistance(strains_dictionary[parent_branch], strains_dictionary[child_branch])
                        except:
                            every_strain_matched_dictionary[parent_branch + "_" + child_branch] = 0.00001
                    elif which_calculation == "K2P":
                        try:
                            every_strain_matched_dictionary[parent_branch + "_" + child_branch] = K2Pdistance(strains_dictionary[parent_branch], strains_dictionary[child_branch])
                        except:
                            every_strain_matched_dictionary[parent_branch + "_" + child_branch] = 0.00001
                    elif which_calculation == "TamD":
                        try:
                            every_strain_matched_dictionary[parent_branch + "_" + child_branch] = Tamuradistance(strains_dictionary[parent_branch], strains_dictionary[child_branch])
                        except:
                            every_strain_matched_dictionary[parent_branch + "_" + child_branch] = 0.00001
    #open(CoreOnly + SC + IncludingMasked + OnlyS28C + '.json')

    print(f"{datetime.datetime.now()}: Distances dictionary ready to save... ")
    with open("JSON_Files\\" + CoreOnly + SC + IncludingMasked + OnlyS28C + SCset + cigar_only + which_calculation + DE_yes+ '.json', 'w') as f:
        json.dump(every_strain_matched_dictionary, f)
    f.close()

    print(f"{datetime.datetime.now()}: Main distances dictionary saved and loaded... ")


starting_branches = {first_branch(strains_dictionary): ""}  # add first branch
while len(starting_branches) < branch_number:  # keep adding braches, which are strains most dissimilar to current strains
    starting_branches[next_branch(starting_branches, strains_dictionary)] = ""


all_branches = []
parent_tree_dictionary = {}
for eachbranch in starting_branches:
    all_branches.append(eachbranch)
    parent_tree_dictionary[eachbranch] = ""

print(f"{datetime.datetime.now()}: Main tree roots constructed... ")

child_to_parent_dictionary = {}
while len(all_branches) < len(strains_dictionary):
    sub_tree_placing = adding_sub_strain(all_branches, strains_dictionary)
    data = sub_tree_placing.split("-")
    parent_branch = data[0]
    child_branch = data[1]
    matching_score = data[2]
    child_to_parent_dictionary[data[1]] = data[0]
    while True:
        child_branch = child_branch
        parent_branch = child_to_parent_dictionary[child_branch]
        if parent_branch in parent_tree_dictionary:
            break
        child_branch = parent_branch
    if parent_branch in parent_tree_dictionary:
        if parent_tree_dictionary[parent_branch] == "":
            parent_tree_dictionary[parent_branch] = sub_tree_placing
        else:
            parent_tree_dictionary[parent_branch] = parent_tree_dictionary[parent_branch] + ";" + sub_tree_placing
    else:
        parent_tree_dictionary[data[0]] = sub_tree_placing

    all_branches.append(data[1])  # start filling up used branches

print(f"{datetime.datetime.now()}: Tree built... ")
## now, trees built. Cluster percentage time!

Number_of_AncestralGenomes = 10  # change this number the number of average genomes.
Eachstrain_Master_Cluster_Dictionary = {}

for AncestralGenome in range(0, Number_of_AncestralGenomes):
    main_branch_average_genomes = {}
    cluster_names = []
    length_of_genome = 0
    for mainbranch in parent_tree_dictionary:
        branch_genomes = []
        branch_genomes.append(strains_dictionary[mainbranch]) # this one first to add main one!
        if parent_tree_dictionary[mainbranch] == "":
            genome = strains_dictionary[mainbranch]  # parentbranch-childbranch-score
            branch_genomes.append(genome)
            cluster_names.append(mainbranch)
            cluster_genome = average_of_strains(branch_genomes)
            main_branch_average_genomes[mainbranch] = cluster_genome.replace(" ", "")
        else:

            for minibranch in parent_tree_dictionary[mainbranch].split(";"):
                genome = strains_dictionary[minibranch.split("-")[1]]  # parentbranch-childbranch-score
                branch_genomes.append(genome)
            cluster_names.append(mainbranch)
            cluster_genome = average_of_strains(branch_genomes)
            main_branch_average_genomes[mainbranch] = (cluster_genome).replace(" ", "")
            length_of_genome = len(cluster_genome.replace(" ", ""))
    ## try Number 2
    number_of_segments = float(100/percent_accuracy)
    segment_stepsize = int(str(float(length_of_genome/number_of_segments)).split(".")[0])

    eachstrain_matched_to_clusters = {}

    counter_for_test = 1
    for eachstrain in strains_dictionary:
        strain_counter_list = [0] * len(main_branch_average_genomes)
        matchup_to_previous = 0  # let default be 1... will bring minor error if first segment has multiple same-score matches.
        for step in range(0, length_of_genome - 1-segment_stepsize, segment_stepsize):
            endstep = step + segment_stepsize
            if endstep + segment_stepsize >= length_of_genome:  # if final segment before remainder, add remainder to final segment
                endstep = length_of_genome - 1
            strain_clustering = []
            for mainbranch in main_branch_average_genomes:
                if main_branch_average_genomes[mainbranch][step:endstep] == "":
                    print("******\nSubset data:")
                    print(mainbranch)
                    print(main_branch_average_genomes[mainbranch][step:endstep])
                    print(eachstrain)
                    print(strains_dictionary[eachstrain][step:endstep])
                    print("*****")

                if (main_branch_average_genomes[mainbranch][step:endstep] or strains_dictionary[eachstrain][step:endstep]) == "":
                    pass
                else:
                    if which_calculation == "Percentage":
                        try:
                            strain_clustering.append(pdistance(main_branch_average_genomes[mainbranch][step:endstep], strains_dictionary[eachstrain][step:endstep]))  # match averaged genome of cluster to that of strain. get score. append score. now clustered by %!
                        except:
                            strain_clustering.append(1000)
                    elif which_calculation == "JCD":
                        try:
                            strain_clustering.append(Jukes_Cantor_distance(main_branch_average_genomes[mainbranch][step:endstep], strains_dictionary[eachstrain][step:endstep]))
                            #print(main_branch_average_genomes[mainbranch][step:endstep] + "___" + str(Jukes_Cantor_distance(main_branch_average_genomes[mainbranch][step:endstep], strains_dictionary[eachstrain][step:endstep])))
                        except:
                            strain_clustering.append(1000)
                    elif which_calculation == "TND":
                        try:
                            strain_clustering.append(TNdistance(main_branch_average_genomes[mainbranch][step:endstep], strains_dictionary[eachstrain][step:endstep]))
                        except:
                            strain_clustering.append(1000)
                    elif which_calculation == "K2P":
                        try:
                            strain_clustering.append(K2Pdistance(main_branch_average_genomes[mainbranch][step:endstep], strains_dictionary[eachstrain][step:endstep]))
                        except:
                            strain_clustering.append(1000)
                    elif which_calculation == "TamD":
                        try:
                            strain_clustering.append(Tamuradistance(main_branch_average_genomes[mainbranch][step:endstep], strains_dictionary[eachstrain][step:endstep]))
                        except:
                            strain_clustering.append(1000)
            for i in range(0, len(strain_clustering)):
                if strain_clustering[i] == None:
                    strain_clustering[i] = 1000
            smallestvalue = sorted(strain_clustering)[0]  # sort values, get smallest (small value = similar DNA!!)

            if strain_clustering.count(smallestvalue) == 1:
                for eachhit in range(0, len(strain_clustering)):
                    if strain_clustering[eachhit] == smallestvalue:
                        strain_counter_list[eachhit] += 1
                        matchup_to_previous = eachhit
            else:
                strain_counter_list[matchup_to_previous] += 1
            if endstep == (length_of_genome - 1):
                break
            counter_for_test += 1
        sum_of_scores = sum(strain_counter_list)
        for eachcluster in range(0, len(strain_counter_list)):
            strain_counter_list[eachcluster] = str(float(strain_counter_list[eachcluster]/sum_of_scores))
        eachstrain_matched_to_clusters[eachstrain] = "\t".join(strain_counter_list)
        counter_for_test = 0
        if eachstrain in Eachstrain_Master_Cluster_Dictionary:
            Eachstrain_Master_Cluster_Dictionary[eachstrain] = Eachstrain_Master_Cluster_Dictionary[eachstrain] + ["\t".join(strain_counter_list)]
        else:
            Eachstrain_Master_Cluster_Dictionary[eachstrain] = ["\t".join(strain_counter_list)]
    print(f"{datetime.datetime.now()}: Average genomes constructed... {str(float(int(AncestralGenome + 1)/Number_of_AncestralGenomes) * 100)}%")


print(f"{datetime.datetime.now()}: Length of main branch average genomes = {len(main_branch_average_genomes)}")
eachstrain_matched_to_clusters = {}
for eachstrain in Eachstrain_Master_Cluster_Dictionary:
    specific_clusters = [[]] * branch_number
    mean_cluster = [""] * branch_number
    stdev_cluster = [""] * branch_number
    for eachcluster in range(0, len(Eachstrain_Master_Cluster_Dictionary[eachstrain])):
        cluster_group = Eachstrain_Master_Cluster_Dictionary[eachstrain][eachcluster].split("\t")
        for cluster in range(0, len(cluster_group)):

            specific_clusters[cluster] = specific_clusters[cluster] + [float(cluster_group[cluster])]
    for cluster in range(0, len(specific_clusters)):
        mean_cluster[cluster] = str(float(sum(specific_clusters[cluster])/len(specific_clusters[cluster])))
        stdev_cluster[cluster] = str(statistics.stdev(specific_clusters[cluster]))
    eachstrain_matched_to_clusters[eachstrain] = "\t".join(mean_cluster) + "\t" + "\t".join(stdev_cluster)


##to here

if full_data:
    print("Strain", "\t".join(cluster_names), "Stdev_cluster1", "Stdev_cluster1", "Stdev_cluster1", sep="\t", file=file_for_Q_clusters)
    for eachstrain in eachstrain_matched_to_clusters:
        print(eachstrain, eachstrain_matched_to_clusters[eachstrain], sep="\t", file=file_for_Q_clusters)
else:
    print("Strain", "\t".join(cluster_names), sep="\t", file=file_for_Q_clusters)
    for eachstrain in eachstrain_matched_to_clusters:
        print(eachstrain, "\t".join(eachstrain_matched_to_clusters[eachstrain].split("\t")[0:3]), sep="\t", file=file_for_Q_clusters)