"""
Protein Structural Feature Predictor
------------------------------------
This script loads DNA sequences from a FASTA file, translates them into
peptides, and performs several lightweight structural analyses including:

- Motif scanning using regex patterns
- Disorder prediction based on residue type composition
- Secondary structure prediction using residue propensities
- Detection of hydrophobic core residues
- Domain boundary inference using disordered regions
- Approximate plDDT-like confidence scoring

This is a simplified, self-contained imitation of features often seen
in protein structure prediction tools.
"""

import re

aaMass = {
'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C': 121.158,
'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
}

aaHydrophobic = "AILMFWVPG"
aaPolar = "STNQY"
aaCharged = "KRDEH"

aaHelix = "ALMQEKH"
aaSheet = "VIYFWT"
aaCoil = "PGSND"

stopCodons = ["TAA", "TAG", "TGA"]

rnaCodonTable = {
# RNA codon table
# U
'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C', # UxU
'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C', # UxC
'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-', # UxA
'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W', # UxG
# C
'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R', # CxU
'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', # CxC
'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', # CxA
'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', # CxG
# A
'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S', # AxU
'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', # AxC
'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', # AxA
'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', # AxG
# G
'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G', # GxU
'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', # GxC
'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', # GxA
'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G' # GxG
}

motifDict = {
    "Cys-rich zinc-binding": "CXXC",
    "N-glycosylation": "NXS",
    "Phosphorylation": "[ST]P",
    "Acidic patch": "E[DE]E",
    "Basic cluster": "KKXK",
    "Tyrosine phosphorylation": "Y[ED]",
    "Cys-zic cluster": "CXXXC",
    "Basic NLS": "K[KR]{2.}",
    "RNA-binding": "RXXR",
    "Protease cleavage sight": "QG",
}

class sequenceFetcher:
    """
    Class responsible for loading FASTA entries.
    Attributes:
        filename (str): Path to the FASTA file.
    Methods:
        getEntries() -> list:
            Returns a list of (header, sequence) tuples parsed by fastaReader.
    """
    def __init__(self, filename):
        self.filename = filename

    def getEntries(self):
        reader = fastaReader(self.filename)
        return list(reader.readFasta())

class fastaReader:
    """
    Basic FASTA file reader.
    This generator function streams each FASTA entry one at a time and yields
    (header, sequence) pairs.
    """
    def __init__(self, filename):
        self.filename = filename

    def readFasta(self):
        """
        Reads the FASTA file line by line and yields (header, sequence)
        tuples. Works with multi-line sequences and ignores blank lines.
        """
        header = ''
        sequence = ''

        try:
            file = open(self.filename)
        except:
            return

        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header:
                    yield header, sequence
                header = line[1:]
                sequence = ''
            else:
                sequence += line.upper()

        if header:
            yield header, sequence

        file.close()



class Protein:
    """
    Protein object that stores DNA/peptide data and performs structural analysis.
    Key analyses performed:
    - Translation of DNA → peptide
    - Motif scanning
    - Disorder prediction from amino-acid composition
    - Secondary structure prediction via residue propensity rules
    - Hydrophobic core detection
    - Domain prediction using disorder boundaries
    - Approximate plDDT confidence scoring
    """
    def __init__(self, header, dna):
        self.header = header
        self.dna = dna
        self.peptide = ""
        self.motifs = []
        self.disordered = []
        self.hydrophobicCore = []
        
    def findDisorder(self):
        """
        Predicts disordered regions based on long stretches (≥8) of polar or charged
        residues. This is a simple heuristic meant to approximate real disorder
        predictors like IUPred or AlphaFold pLDDT values.
        """
        pep = self.peptide
        region = ""
        start = 0

        for i, aa in enumerate(pep):
            if aa in aaPolar or aa in aaCharged:
                if region == "":
                    start = i
                region += aa
            else:
                if len(region) >= 8:
                    self.disordered.append((start, i-1))
                region = ""

        if len(region) >= 8:
            self.disordered.append((start, len(pep)-1))

    def predictSecondaryStructure(self):
        """
        Assigns a secondary structure state (H, E, C) to each residue purely based
        on amino-acid preference groups. This is not a true prediction model, but
        works as a heuristic summary.
        """
        ss = []
        for aa in self.peptide:
            if aa in aaHelix:
                ss.append("H")
            elif aa in aaSheet:
                ss.append("E")
            elif aa in aaCoil:
                ss.append("C")
            else:
                ss.append("C")   # default to coil
        self.secondaryStructure = "".join(ss)

    def secondaryStructureSummary(self):
        """
        Returns a summary of secondary structure counts:
        H = helix, E = sheet, C = coil
        """
        if not self.secondaryStructure:
            return {"H": 0, "E": 0, "C": 0}
        
        counts = {"H": 0, "E": 0, "C": 0}
        for ss in self.secondaryStructure:
            if ss in counts:
                counts[ss] += 1
            else:
                counts["C"] += 1  # default to coil if unknown
        return counts

    
    def plddtSummary(self):
        """
        Approximate high and low confidence residues based on disorder.
        - Low confidence: residues in disordered regions
        - High confidence: residues outside disordered regions
        Returns (high_count, low_count)
        """
        low_indices = set()
        for start, end in self.disordered:
            low_indices.update(range(start, end+1))

        total_residues = len(self.peptide)
        low = len(low_indices)
        high = total_residues - low
        return high, low

    def avgConfidence(self, start=None, end=None):
        """
        Computes an approximate confidence score over a region by assigning
        residue-type weights. Charged regions are treated as lower confidence
        and hydrophobic residues as higher.
        Intended as a rough heuristic for reporting.
        """
        if not self.peptide:
            return 0
    
        scores = []
        if start is None and end is not None:
            for i in range(start, end+1):
                aa = self.peptide[i]
                if aa in aaCharged:
                    score = 30
                elif aa in aaPolar:
                    score = 50
                else:
                    score = 60
                scores.append(score)
        else:
            for s, e in self.disordered:
                for i in range(s, e+1):
                    aa = self.peptide[i]
                    if aa in aaCharged:
                        score = 30
                    elif aa in aaPolar:
                        score = 50
                    else:
                        score = 60
                    scores.append(score)
    
        return sum(scores) / len(scores) if scores else None
    
    def translate(self):
        """
        Translates the DNA sequence into a peptide using the RNA codon table.
        Translation stops at the first stop codon encountered.
        """
        pep = ""
        for i in range(0, len(self.dna)-2, 3):
            codon = self.dna[i:i+3]
            aa = rnaCodonTable.get(codon, 'X')
            if aa == '*':
                break
            pep += aa
        self.peptide = pep

    def analyze(self):
        """
        Runs the main analysis steps:
        - Motif detection
        - Disorder prediction
        - Hydrophobic core identification
        """
        self.findMotifs()
        self.findDisorder()
        self.findHydrophobicCore()

    def findMotifs(self):
        """
        Uses regex scanning to detect known motifs from `motifDict`. The 'X'
        wildcard in motif definitions is converted into a '.' regex wildcard.
        Reports motif name and residue index.
        """
        pep = self.peptide
        self.motifs = []
        
        for motifName, pattern in motifDict.items():
                regex = pattern.replace("X", ".")
                match = re.finditer(regex, pep)
                for m in match:
                    self.motifs.append((motifName, m.start()))

    def findHydrophobicCore(self):
        """
        Identifies hydrophobic residues (as defined in aaHydrophobic).
        This is a minimal approximation of core-forming positions.
        """
        pep = self.peptide
        for i, aa in enumerate(pep):
            if aa in aaHydrophobic:
                self.hydrophobicCore.append(i)

    def domainFinder(self):
        """
        Predicts domain boundaries by splitting the protein at disordered regions.
        Structured regions between disordered segments are treated as domains.
        Returns:
            List of (start, end) domain coordinate pairs (1-indexed).
        """
        if not self.disordered:
            # Whole protein is one domain
            return [(1, len(self.peptide))]
    
        domains = []
        disorder = sorted(self.disordered)  # ensure ordered
        length = len(self.peptide)
    
        # 1. Start → first disordered region
        first_start = disorder[0][0]
        if first_start > 0:
            domains.append((1, first_start))
    
        # 2. Middle structured regions
        for i in range(len(disorder) - 1):
            end_current = disorder[i][1]
            start_next = disorder[i+1][0]
    
            if start_next - end_current > 1:
                domains.append((end_current+1, start_next))
    
        # 3. Last disordered → end
        last_end = disorder[-1][1]
        if last_end < length:
            domains.append((last_end+1, length))
    
        return domains

class outputFormatter:
    """
    Formats all computed protein information into a user-readable text report.
    """
    def __init__(self, protein):
        self.protein = protein

    def buildReport(self):
        """
        Constructs a structured report summarizing:
        - Motifs
        - Secondary structure counts
        - Disordered regions + average confidence
        - Overall high/low confidence residue counts
        - Hydrophobic residues
        - Domain predictions
        """
        out = []
        out.append("--- Structural Summary ---")
        out.append("File: " + self.protein.header)
        out.append("Total Residues: " + str(len(self.protein.peptide)))
        out.append("Peptide Preview: " + self.protein.peptide[:60] + "...")

        out.append("\nMotifs:")
        if not self.protein.motifs:
            out.append("None")
        else:
            for motifName, idx in self.protein.motifs:
                out.append(f"{motifName} at position {idx}")

        out.append("\nSecondary Structure Summary:")
        summary = self.protein.secondaryStructureSummary()
        out.append(f"-- Alpha Helices: {summary['H']} residues")
        out.append(f"-- Beta-Sheets: {summary['E']} residues")
        out.append(f"-- Coils/Loops: {summary['C']} residues")
        
        out.append("\nDisordered Regions:")
        if not self.protein.disordered:
            out.append("None")
        else:
            for i, (start, end) in enumerate(self.protein.disordered, 1):
                score = self.protein.avgConfidence(start, end)
                out.append(f"-- Disordered region {i}: {start}-{end} (Score: {score:.2f})")

        high, low = self.protein.plddtSummary()
        out.append("\nOverall plDDT Confidence")
        out.append(f"-- High Confidence (>90) residues: {high}")
        out.append(f"-- Low Confidence (<50) residues: {low}")
        
        out.append("\nHydrophobic Core Residues (indices):")
        if not self.protein.hydrophobicCore:
            out.append("None")
        else:
            out.append("" + ", ".join(str(i) for i in self.protein.hydrophobicCore))

        out.append("\n--- Domain Analysis ---")
        out.append("Domains:")
        domains = self.protein.domainFinder()

        if not domains:
            out.append("None")
        else:
            for i, (start, end) in enumerate(domains, 1):
                out.append(f"Domain {i} (Residues {start}-{end})")

        return "\n".join(out)

def main():
    """
    Entry point for the program.
    Loads sequences from 'sequences.fasta', translates each into a protein,
    performs structural analysis, and prints a formatted report.
    """
    fetcher = sequenceFetcher("sequences.fasta")   # your dictionary file
    entries = fetcher.getEntries()

    for title, dna in entries:
        protein = Protein(title, dna)
        protein.translate()
        protein.predictSecondaryStructure()
        protein.analyze()
        high, low = protein.plddtSummary()

        print("-" * 100)
        
        formatter = outputFormatter(protein)
        print(formatter.buildReport())

if __name__ == "__main__":
    main()
