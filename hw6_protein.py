"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    file=open(filename,"r")
    lines=file.read()
    file.close()
    nlines=lines.splitlines()
    empty_string=""
    for line in nlines:
        empty_string+=line
    return empty_string

'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    x=dna.replace("T", "U")
    empty_lst=[]
    for i in range(startIndex, len(x), 3):
        empty_lst.append(x[i:i+3])
        if x[i:i+3]=="UAG" or x[i:i+3]=="UAA" or x[i:i+3]=="UGA":
            break
    return empty_lst

'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):
    import json
    dictionary={}
    file=open(filename)
    new_file=json.load(file)
    for i,j in new_file.items():
        for k in j:
            dictionary[k.replace("T", "U")]=i
    return dictionary


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    protein_list=[]
    if codons[0]=="AUG":
        protein_list.append("Start")
    for i in range(1, len(codons)):
        if codons[i] in codonD.keys():
            protein_list.append(codonD[codons[i]])
    return protein_list


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    file1=readFile(dnaFilename)
    file2=makeCodonDictionary(codonFilename)
    count=0
    new_list=[]
    i=0
    while i<len(file1):
        if file1[i:i+3]=="ATG":
            file3=dnaToRna(file1, i)
            proteins=generateProtein(file3, file2)
            new_list.append(proteins)
            i=i+3*len(file3)
        else:
            i=i+1
            count+=1
    return new_list


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    unique_list=[]
    for i in proteinList1:
        for j in proteinList2:
            if i==j and i not in unique_list:
                unique_list.append(i)
    return unique_list


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    combine_list=[]
    for i in proteinList:
        for j in i:
            if i not in combine_list:
                combine_list.append(j)
    return combine_list


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    amino_list={}
    for i in aaList:
        if i not in amino_list:
            amino_list[i]=1
        else:
            amino_list[i]+=1
    return amino_list


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    combine_list1=combineProteins(proteinList1)
    combine_list2=combineProteins(proteinList2)
    amino_dict1=aminoAcidDictionary(combine_list1)
    amino_dict2=aminoAcidDictionary(combine_list2)
    frequency_dict1={}
    frequency_dict2={}
    temporary_list=[]
    frequency_difference=[]
    for i in amino_dict1:
        frequency_dict1[i]=amino_dict1[i]/len(combine_list1)
        if i not in temporary_list and i!="Start" and i!="Stop":
            temporary_list.append(i)
    for j in amino_dict2:
        frequency_dict2[j]=amino_dict2[j]/len(combine_list2)
        if j not in temporary_list and j!="Start" and j!="Stop":
            temporary_list.append(j)
    for k in temporary_list:
        frequency1=0
        frequency2=0
        if k in frequency_dict1:
            frequency1=frequency_dict1[k]
        if k in frequency_dict2:
            frequency2=frequency_dict2[k]
        difference=frequency2-frequency1
        if difference < -cutoff or difference > cutoff:
            frequency_difference.append([k, frequency1, frequency2])
    return frequency_difference


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    print("The following proteins occurred in both DNA sequences")
    for i in commonalities:
        proteins_count=""
        list1=i[1:(len(i)-1)]
        count=0
        for j in list1:
            proteins_count+=j
            count+=1
            if count!=len(list1):
                proteins_count+="-"
        if len(proteins_count)!=0:
            print(proteins_count)
    print("The following amino acids occurred at very different ratesin two DNA sequences")
    for i in differences:
        x=i[0]
        frequency1=round(i[1]*100, 2)
        frequency2=round(i[2]*100, 2)
        print(str(x)+" "+str(frequency1)+" % in Seq1"+","+str(frequency2)+"% in Seq2")       
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    return


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    return


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    import numpy
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    test.week1Tests()
    # test.testSynthesizeProteins()
    print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    runWeek1()

    ## Uncomment these for Week 2 ##
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    # test.testFindAminoAcidDifferences()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()


    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """
