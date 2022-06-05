#!/usr/bin/python3

# Integrantes:
# Paloma Abigail Ramírez Pizaña
# Daniel Linares Gil

import re

CODONES_TRADUCCION = {
    # Alanina
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    # Cisteina
    "UGU":"C", "UGC":"C",
    # Acido aspartico
    "GAU":"D", "GAC":"D",
    # Acido glutamico
    "GAA":"E", "GAG":"E",
    # Fenilalanina
    "UUU":"F", "UUC":"F",
    # Glicina
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    # Histidina
    "CAU":"H", "CAC":"H",
    # Isoleucina
    "AUA":"I", "AUU":"I", "AUC":"I",
    # Lisina
    "AAA":"K", "AAG":"K",
    # Leucina
    "UUA":"L", "UUG":"L", "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    # Metionina
    "AUG":"M", 
    # Aspargina
    "AAU":"N", "AAC":"N",
    # Prolina
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    # Glutamina
    "CAA":"Q", "CAG":"Q",
    # Arginina
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R", "AGA":"R", "AGG":"R",
    # Serina
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S", "AGU":"S", "AGC":"S",
    # Treonina
    "ACU":"U", "ACC":"U", "ACA":"U", "ACG":"U",
    # Valina
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    # Triptofano
    "UGG":"W",
    # Tirosina
    "UAU":"Y", "UAC":"Y",
    # Stop
    "UAA":"_", "UAG":"_", "UGA":"_"
}

CODON_INICIO = "ATG"

CODONES_TERMINO = {"TAA", "TAG", "TGA"}

NUCLEOTIDOS = {"A", "T", "C", "G"}

class Gen:
    def __init__(self, gene_str):
        self.gene_str = gene_str

    def __getitem__(self, key):
        start_pos = (key % len(self)) * 3
        return self.gene_str[start_pos:start_pos + 3]

    def __len__(self):
        return len(self.gene_str) // 3

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

patron_adn = re.compile('ATG([ATCG]{3})+(TAG|TAA|TGA)$', re.IGNORECASE)

def validar_secuencia(dna):
    if len(dna) % 3 != 0:
        return "Longitud invalida"
    if dna[0:3] != CODON_INICIO:
        return "Codón de inicio invalido"
    if not patron_adn.match(dna):
        return "Secuencia invalida"
    return None

def complemento(dna):
    dnacomplemento = []
    for i in (dna):
        if i in ["A","a"]:
            dnacomplemento.append("T")
        elif i in ["T", "t"]:
            dnacomplemento.append("A")
        elif i in ["G", "g"]:
            dnacomplemento.append("C")
        elif i in ["C","c"]:
            dnacomplemento.append("G")

    return ''.join(dnacomplemento)

def transcrito(dnacomplemento):
    rnatranscrito = []
    for i in (dnacomplemento):
        if i == "A":
            rnatranscrito.append("U")
        elif i == "T":
            rnatranscrito.append("A")
        elif i == "G":
            rnatranscrito.append("C")
        elif i == "C":
            rnatranscrito.append("G")

    return ''.join(rnatranscrito)

def cadena_aa(rna):
    return ''.join(CODONES_TRADUCCION.get(c, '-') for c in Gen(rna))

if __name__ == '__main__':
    adn = input("Introduzca la sequencia de ADN: ").upper()
    e = validar_secuencia(adn)
    if not e:
        com = complemento(adn)
        arn = transcrito(com)
        aa = cadena_aa(arn)

        print("Complemento: ", com)
        print("ARN transcrito: ", arn)
        print("Aminoácidos: ", aa)
    else:
        print(e)

