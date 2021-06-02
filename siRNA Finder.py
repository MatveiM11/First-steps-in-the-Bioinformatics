DNA_seq = input("")

RNA_seq  = DNA_seq.replace("T", "U")

N = RNA_seq
iRNA = (N.translate(str.maketrans({"A": "U", "G": "C", "U": "A", "C": "G"})))
iRNA_str = iRNA 

K = 22
iRNA_comb = [iRNA[i: j] for i in range(len(iRNA_str)) for j in range(i + 1, len(iRNA_str) + 1) if len(iRNA_str[i:j]) == K]
print("All iRNA combinations", iRNA_comb)

siRNAs = iRNA_comb

for x in siRNAs:
    print('siRNA: ', x, ' percentage: ',((x.count("C") + x.count("G"))) / len(x) * 100, '%')

output = [x for x in siRNAs if 30 <= ((x.count("C") + x.count("G"))) / len(x) * 100 <=50]

print('output: ', output)
