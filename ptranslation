import re

RNA_translation = input("")
RNA_dict = {
  "GCA":"A", "GCC":"A", "GCG":"A", "GCU":"A",
  "UGC":"C", "UGU":"C", "GAC":"D", "GAU":"D",
  "GAA":"E", "GAG":"E", "UUC":"F", "UUU":"F",
  "GGA":"G", "GGC":"G", "GGG":"G", "GGU":"G",
  "CAC":"H", "CAU":"H", "AUA":"I", "AUC":"I",
  "AUU":"I", "AAA":"K", "AAG":"K", "UUA":"L",
  "UUG":"L", "CUA":"L", "CUC":"L", "CUG":"L",
  "CUU":"L", "AUG":"START", "AAC":"N", "AAU":"N",
  "CCA":"P", "CCC":"P", "CCG":"P", "CCU":"P",
  "CAA":"Q", "CAG":"Q", "AGA":"R", "AGG":"R",
  "CGA":"R", "CGC":"R", "CGU":"R", "CGG":"R",
  "AGC":"S", "AGU":"S", "UCA":"S", "UCC":"S",
  "UCG":"S", "UCU":"S", "ACA":"T", "ACC":"T",
  "ACG":"T", "ACU":"T", "GUA":"V", "GUC":"V",
  "GUG":"V", "GUU":"V", "UGG":"W", "UAC":"Y",
  "UAU":"Y", "UAG":"STOP", "UAA":"STOP", "UGA":"STOP"
}
# inv_dict = {v: k for k, v in ini_dict.items()} - Useful regex Uga-Buga, if you need to invert the dictionary

regex = "|".join(RNA_dict.keys())
translation = re.sub(regex, lambda m: RNA_dict[m.group()], RNA_translation)

print(translation)
