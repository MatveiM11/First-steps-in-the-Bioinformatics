#Сделали - Матвей Мочалов(https://github.com/MatveiM11) и Митя Каштанов(https://github.com/duoapassivist)

# sequence.fasta 
get = input("Введите полное имя файла: ") # В файле должна быть только одна последовательность мРНК 
gene = open(get, 'r') 
 
lenth = len(gene.readline()) 
line = str(gene.read())[lenth:] 
 
for i in line: 
    if i == "\n": 
        line = line.replace("\n","") 
 
line = line.replace("T", "U") #производим замену нуклеотидов (тимина на урацил) 
 
gene.close() 
 
# code region 
small = [] 
for n in range(len(line)): 
    if line.find("AUG") == n: 
        line = line[n:] 
        break 
for i in range(101, len(line)): 
    if line.find('UAG') == i or line.find('UAA') == i or line.find('UGA') == i: 
        codes_RNA = line[:i+3] 
        UTR = line[i+3:] 
        break 
 
senseRNA = (codes_RNA.translate(str.maketrans({ 
    "A": "U", 
    "G": "C", 
    "U": "A", 
    "C": "G" 
})))[::-1]  #создаем комплементарную рнк 
 
K = 21 # количество нуклеотидов в siRNA 
res = [senseRNA[i: j] for i in range(len(senseRNA)) for j in range(i + 1, len(senseRNA) + 1) if len(senseRNA[i:j]) == K] 
#код, основная идея которого, в том что мы просто смещаем рамку в 21 нуклеотид на один нуклеотид каждый раз 
 
# seed region 
for x in res: 
    seed = x 
    if UTR.translate(str.maketrans({"A":"U", "G":"C", "U":"A", "C":"G"}))[::-1].count(seed[2:9]) != 0: 
        res.remove(x) 
 
top = dict() 
 
# трансформация в shRNA 
antisense_siRNAs = [] 
 
for s in range(len(res)): 
    antisense_siRNAs.append(str(res[s].translate(str.maketrans({"A": "T", "G": "C", "U": "A", "C": "G"})))) 
 
res_DNA = [] 
for s in range(len(res)): 
    res_DNA.append(str(res[s].replace("U", "T"))) 
 
loop = 'TCAAGAG' 
shRNAs = [] 
 
for i in range(len(res)): 
    shRNAs.append(str(res_DNA[i]) + str(loop) + str(antisense_siRNAs[i])) 
 
for n in range (len(res)): 
    a = str(res[n]) 
    score = 0 # подсчет баллов 
    rep = 0 # подсчет повторений 
 
    if lambda x: (.3 <= (x.count('G') + x.count('C')) / len(x) <= .52) and (len(x) == 21): 
        score += 1 
 
    # Ui-Tei rule 
    if a[-1] == "A" or a[-1] != "U": # A/U на конце 
        score += 1 
    if a[18] == "C" or a[18] == "G": # C/G на 19 
        score += 1 
    if a[:6].count("A") + a[:6].count("U") >= 4: # A/U >= 4 в местах 1-7 
        score += 1 
    if a[9] != "C" or a[9] != "G": # C/G не 10 
        score += 1 
 
    # Reynolds rule 
    for z in range(3, 11): # нет повторений (больше чем 3 нуклеотида) 
        for i in range(len(a) - z): 
            if a[i:i + z] == a[i + z:i + z * 2]: 
                rep += 1 
    if rep == 0: 
        score += 1 
    if a[0] != "C" or a[0] != "G": # C/G не 1 
        score += 1 
    if [16] == "A": # A на 17 
        score += 1 
    if a[9] == "U": # U на 10 
        score += 1 
    if a[:5].count("A") + a[:5].count("U") >= 3: # на 1-5 местах A/U больше или равно 3 
        score += 1 
 
    # Amarzguioui rule 
    if a[6] != "G": # G не 7 
        score += 1 
    if a[0] == "A" or a[0] == "U": # A/U на 1 
        score += 1 
    if a[13] == "A" or a[13] == "U": # A/U на 14 
        score += 1 
    if a[18] != "U": # U не 19 
        score += 1 
 
    # Параметры эффективности 
    if a.count("AAC") >= 1 or a.count("UC") >= 1 or a.count("UG") >= 1 or a.count("AAG") >= 1 or a.count("AGC") >= 1 or a.count("UCU") >= 1 or a.count("UCCG") >= 1 or a.count("CUU") >= 1 or a.count("CU") >= 1 or a.count("GUU") >= 1 or a.count("UCC") >= 1 or a.count("CG") >= 1 or a.count("AUC") >= 1 or a.count("GCG") >= 1 or a.count("UUU") >= 1 or a.count("ACA") >= 1 or a.count("UUC") >= 1 or a.count("CAA") >= 1: 
        score += 1 
    if a.count("CUU") == 0 or a.count("CUA") == 0 or a.count("GUU") == 0 or a.count("GU") == 0 or a.count("GAU") == 0 or a.count("ACGA") == 0 or a.count("GCC") == 0 or a.count("GUGC") == 0 or a.count("GGC") == 0 or a.count("CCG") == 0 or a.count("CAG") == 0 or a.count("GAG") == 0 or a.count("GCA") == 0 or  a.count("AUA") == 0 or a.count("CUG") == 0 or a.count("AG") == 0 or a.count("GG") == 0 or a.count("GGA") == 0:  # Avoid motifs 
        score += 1 
 
    rep = 0 
 
    siRNA = [] 
    gins = [] 
    # энергия Гиббса для пар нуклеотидов -

GIBS 
    for x in res: 
        g = 0 
        if 'AA' or 'UU' in x: 
            g += x.count('AA') * (-0.9) + x.count('UU') * (-0.9) 
        if 'AU' in x: 
            g += x.count('AU') * (-0.9) 
        if 'UA' in x: 
            g += x.count('UA') * (-1.1) 
        if 'CA' or 'UG' in x: 
            g += x.count('CA') * (-1.8) + x.count('UG') * (-1.8) 
        if 'CU' or 'AG' in x: 
            g += x.count('CU') * (-1.7) + x.count('AG') * (-1.7) 
        if 'GA' or 'UC' in x: 
            g += x.count('GA') * (-2.3) + x.count('UC') * (-2.3) 
        if 'GU' or 'AC' in x: 
            g += x.count('GU') * (-2.1) + x.count('AC') * (-2.1) 
        if 'CG' in x: 
            g += x.count('CG') * (-2.0) 
        if 'GC' in x: 
            g += x.count('GC') * (-3.4) 
        if 'GG' or 'CC' in x: 
            g += x.count('GG' or 'CC') * (-2.9) 
        gins.append(round(g, 1)) 
 
    d = str(res[n]) 
    score = round((score / 17)*100, 1) 
    if score >= 50: 
        top[d] = (score, gins[n], shRNAs[n]) 
 
top_sorted_descending = dict(sorted(top.items(), key=lambda item: (item[1][0], item[1][1], item[1][2]), reverse=True)) 
 
out = open("gene_with_shRNA.txt", "w") 
out.write("       Sequence siRNA       :  Score  :  delta G  :      shRNA\n                                 %      kcal/mol\n") 
for element in top_sorted_descending: 
    out.write('as 5\' ' + element + ' 3\'  ' + str(top_sorted_descending[element][0]) + '     ' + str(top_sorted_descending[element][1]) + '    : 5\' ' + 
                str(top_sorted_descending[element][2]) + ' 3\'\n' + 'ss 3\' ' + element.translate(str.maketrans({"A": "U", "G": "C", "U": "A", "C": "G"})) + 
              ' 5\'\n\n') 
out.close()
