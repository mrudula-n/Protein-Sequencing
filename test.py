# data="data\human_p53.txt"
# f=open(data)
# file=f.read()
# f.close()
# nlines=file.splitlines()
# print(nlines)

# strng=""
# for line in nlines:
#     strng+=line
# print(strng)

# print(file)
dna="ATGGATGGACTCTAACTCATGCCCTTTTAG"
startIndex=0
def dnaToRna(dna, startIndex):
    x=dna.replace("T", "U")
    empty_lst=[]
    for i in range(startIndex, len(x), 3):
        empty_lst.append(x[i:i+3])
        if x[i:i+3]=="UAG" or x[i:i+3]=="UAA" or x[i:i+3]=="UGA":
            break
    return empty_lst

print(dnaToRna(dna, startIndex))