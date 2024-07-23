'''1. Column 1 (Gene): lists the names of genes that are found in both the reference sequence and the  variant. The list should be sorted alphabetically using Pythons sorted() function. 
The portion above shows 5 genes, but the complete table should have 11 genes. For the gene name, I used the first 10  letters after the bar in the comment line. 

2. Column 2 (NC length): the length of the gene in the reference sequence. Remember to remove the  newline characters. 

3. Column 3 (OR length): similarly, the length of the variant gene. 

4. Column 4 (Differences): A list of the first four differences between the reference gene and the variant gene. For example, in Row 2, for ORF10 prot, there were no differences. 
In Row 3, ORF1a poly, the  first difference shown is A139G_M, which is supposed to mean that at position 139 the A in the  reference has become G in the variant. The ‘M’ indicates that this would be a “missense” variation,  meaning that the amino acid changes. We are assuming that the codons start with the first letter, so  ATG is the first codon for each. The codon containing position 139 would be at positions 138, 139,  140, since 138%3 is 0. To know if the change is missence (M), silent (S) or nonsense (N), you have to  compare the amino acid produced by the codon at that position in the referencesequence with the  corresponding amino acid produced by the codon in the variant seqence at that spot. Silent (S), means  no change in amino acid; Missense (M) means the amino acids are different and Nonsense (N) means  the amino acid changed to a Stop, i.e., the codon changed to one of TAG, TGA, or TAA. 5. Column 5 (Diff count): is the total number of differences. 
'''
f = open("sequencesfasta.txt")
seq = f.readline()
seq = f.read()
seq = seq.replace('\n', '')
n = len(seq)

aa = {'GGG' : 'G', 'GGA' : 'G', 'GGC' : 'G', 'GGT' : 'G', 
      'AGG' : 'R', 'AGA' : 'R', 'AGC' : 'S', 'AGT' : 'S', 
      'CGG' : 'R', 'CGA' : 'R', 'CGC' : 'R', 'CGT' : 'R', 
      'TGG' : 'W', 'TGA' : 'X', 'TGC' : 'C', 'TGT' : 'C',
      
      'GAG' : 'E', 'GAA' : 'E', 'GAC' : 'D', 'GAT' : 'D',
      'AAG' : 'K', 'AAA' : 'K', 'AAC' : 'N', 'AAT' : 'N', 
      'CAG' : 'Q', 'CAA' : 'Q', 'CAC' : 'H', 'CAT' : 'H', 
      'TAG' : 'X', 'TAA' : 'X', 'TAC' : 'Y', 'TAT' : 'Y', 
      
      'GCG' : 'A', 'GCA' : 'A', 'GCC' : 'A', 'GCT' : 'A',
      'ACG' : 'T', 'ACA' : 'T', 'ACC' : 'T', 'ACT' : 'T', 
      'CCG' : 'P', 'CCA' : 'P', 'CCC' : 'P', 'CCT' : 'P', 
      'TCG' : 'S', 'TCA' : 'S', 'TCC' : 'S', 'TCT' : 'S', 
      
      'GTG' : 'V', 'GTA' : 'V', 'GTC' : 'V', 'GTT' : 'V',
      'ATG' : 'M', 'ATA' : 'I', 'ATC' : 'I', 'ATT' : 'I', 
      'CTG' : 'L', 'CTA' : 'L', 'CTC' : 'L', 'CTT' : 'L', 
      'TTG' : 'L', 'TTA' : 'L', 'TTC' : 'F', 'TTT' : 'F', 
    }  # Amino Acids Dictionary
nc = [] # NC titles
orf = [] #orf titles
nc_seq = [] #sequences
orf_seq = [] 
nucleo = ['A', 'T', 'C', 'G']
current_seq = ''
current_title = ''

for line in seq:
    if line.startswith('>'): #if the line is a title
        if current_seq: #if there's a current sequence, save it by the last title
            if current_title.startswith('NC') or current_title.startswith('join(NC'):
                nc_seq.append(current_seq)
            elif current_title.startswith('ORF') or current_title.startswith('join(ORF'):
                orf_seq.append(current_seq)
            current_seq = ''
        current_title = line[1:].strip() #removes '>'
        if current_title.startswith('NC') or current_title.startswith('join(NC'):
            title = current_title.split('|', 1)[1]
            title = title.split('[', 1)[0].strip()
            orf.append(title)
    
    else: #if the line is a sequence
        for char in line:
            if char.upper() in nucleo:
                current_seq += char

if current_seq: # Operating on the last sequence
    if current_title.startswith('NC') or current_title.startswith('join(NC'):
        nc_seq.append(current_seq)
    elif current_title.startswith('ORF') or current_title.startswith('join(ORF'):
        orf_seq.append(current_seq)

#Comparing NC and ORF sequences and finding differences
def mutation(cod1, cod2):
    a1 = aa.get(cod1)
    a2 = aa.get(cod2)

    if a1 == a2:
        return 'S' # silent
    elif len(a1) != len(a2):
        return 'N' # nonsense
    else:
        return 'M' # missense

def get_diff(seq1, seq2): #function to determine difference between the 2 sequences
    diff = []
    diff_count = 0
    for i in range(0, min(len(seq1), len(seq2)), 3): # 3-steps for codons
        cod1 = seq1[i:i+3]
        cod2 = seq2[i:i+3]

        if cod1 != cod2:
            stat = mutation(cod1, cod2)
            for j in range(3):
                if i+j < min(len(seq1), len(seq2)) and seq1[i+j] != seq2[i+j]: #make sure it's within boundaries
                    diff.append(f"{seq1[i+j]}{i+j}{seq2[i+j]}_{stat}")
                    diff_count += 1
    return diff, diff_count

results = [] #store results in list

for i, orf_title in enumerate(orf):
    orf_seq_arr = orf_seq[i]
    nc_seq_arr = None
    for j, nc_title in enumerate(nc):
        if orf_title in nc_title:
            nc_seq_arr = nc_seq[j]
            break

    if nc_seq_arr:
        diff, diff_count = get_diff(nc_seq_arr, orf_seq_arr)
        results.append({
            "Gene": orf_title,
            "NC length": len(nc_seq_arr),
            "ORF length": len(orf_seq_arr),
            "Differences": diff,
            "Difference count": diff_count })
formatted_op = []
header = ["Gene", "NC length","Differences", "Differences count"]
formatted_op.append(header)

for result in results:
    row = [
        result["Gene"],
        result["NC length"],
        result["ORF length"],
        result["Differences"] if len(result["Differences"]) <= 4 else result["Differences"][:4],
        result["Differences count"]
    ]
    formatted_op.append(row)

formatted_op = sorted(formatted_op)

#print headers
for i in range(0, len(formatted_op [0])):
    if i == 3:
        print(f"{formatted_op[0][i]:<40}", end='\t')
    else:
        print(f"{formatted_op[0][i]:<40}", end= '\t')
print('\n')

#print gene data
for i in range(1, len(formatted_op)):
    for j in range(0, len(formatted_op[i])):
        if j ==3:
            print(f"{str(formatted_op[i][j]):<40}", end = '\t')
        else:
            print(f"{str(formatted_op[i][j]):<25}", end ='\t')
    print('\n')