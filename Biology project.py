while True:

    print("What type of the data you want to process?", "RNA? DNA?", sep = "\n")
    A = input("")
    if A == "RNA":
          print("What you want to do with RNA?", "(a1)RNA to Protein? (a2)RNA to DNA? (a3)RNA Length? (a4)RNA Splicing? (a5)RNA nucleotides (a6)ratio? RNAi? ", sep = "\n")
          break

    if A == "DNA":
          print("What you want to do with DNA", "(b1)DNA to RNA? (b2)DNA length? (b3)DNA nucleotides ratio? (b4)Complimentary DNA chain?", sep = "\n")
          break



    else:
         print("Nope") #TF2 references
         continue


while A == "RNA":
    X = input("")
    if X == "a1":
     print("Input RNA")
     import ptranslation #Written by myself module for RNA translation
     RNA_translation = input("")
     print (RNA_translation)
     break
    if X == "a2":
     print("Input RNA")
     RNA_seq = input("")
     DNA_seq  = RNA_seq.replace("U", "T")
     print (DNA_seq)
     break
    if X == "a3":
     print ("Input RNA")
     RNA = input("")
     RNAo = int (len(RNA)//2)
     print(RNAo)
     break
    if X == "a4":
     print ("Input RNA")
     import re
     RNA = input("")
     regex = r"GU(?:\w{0,}?)AG" #regex black magic Uga-Buga, still need to undestand it better
     exons = re.sub(regex, '', RNA)
     print(exons)
     break
    if X == ("a5"):
     print("Input RNA")
     RNAn = input ("")
     X = (RNAn.count("A")) + (RNAn.count("U")) + (RNAn.count("C")) + (RNAn.count ("G"))
     A = (RNAn.count("A")/X)
     U = (RNAn.count("U")/X)
     G = (RNAn.count("G")/X)
     C = (RNAn.count("C")/X)
     print ("A =", A, "U =", U, "G =", G, "C =", C)
     break
    if X == ("a6"):
        print("Input RNA")
        N = input("")
        print(N.translate(str.maketrans({"A": "U", "G": "C", "U": "A", "C": "G"})))
        break
    else:
      print("NOPE! Input correct option!")
      continue
      break

while A == "DNA":
    Y = input("")
    if Y == "b1":
     print ("Input DNA ")
     DNA_seq = input("")
     RNA_seq  = DNA_seq.replace("T", "U")
     print (RNA_seq)
     break
    if Y == ("b2"):
        print("Input DNA")
        DNA = input("")
        DNAo = int (len(DNA)//2)
        print(DNAo)
        break
    if Y == ("b3"):
        print("Input DNA")
        DNAn = input ("")
        X = (DNAn.count("A")) + (DNAn.count("T")) + (DNAn.count("C")) + (DNAn.count ("G"))
        A = (DNAn.count("A")/X)
        T = (DNAn.count("T")/X)
        G = (DNAn.count("G")/X)
        C = (DNAn.count("C")/X)
        print ("A =", A, "T =", T, "G =", G, "C =", C)
        break
    if Y == ("b4"):
        print("Input DNA")
        Nd = input("")
        print(Nd.translate(str.maketrans({"A": "T", "G": "C", "T": "A", "C": "G"})))
        break
    else:
      print("NOPE! Input correct option!")
      continue
      break
#psyhological support in process of  the coding - M.Kashtanov
