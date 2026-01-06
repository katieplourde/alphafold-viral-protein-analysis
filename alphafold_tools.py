from __future__ import annotations #allows for edits on future classes without completion of previous classes
from dataclasses import dataclass
from pathlib import Path
import subprocess 
from typing import Dict, Optional, Tuple 

import json

GENETIC_CODE = {
    "ATA":"I", "ATC":"I", "ATT":"I", "ATG":"M",  #Reference to take the codon reading frame [i:i+3] to the proteins needed for a submittable FASTA file 
    "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
    "AAC":"N", "AAT":"N", "AAA":"K", "AAG":"K",
    "AGC":"S", "AGT":"S", "AGA":"R", "AGG":"R",
    "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L",
    "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
    "CAC":"H", "CAT":"H", "CAA":"Q", "CAG":"Q",
    "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R",
    "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V",
    "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A",
    "GAC":"D", "GAT":"D", "GAA":"E", "GAG":"E",
    "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G",
    "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S",
    "TTC":"F", "TTT":"F", "TTA":"L", "TTG":"L",
    "TAC":"Y", "TAT":"Y", "TAA":"STOP", "TAG":"STOP",
    "TGC":"C", "TGT":"C", "TGA":"STOP", "TGG":"W",
}
def translate_nuc_to_aa(nuc_seq: str, frame: int = 0):
    nuc_seq = nuc_seq.upper().replace("U", "T") #If the found sequence is RNA puts it into DNA for the reference table , also uppercases it. small edits
    aas: list[str] = [] #begins with an empty frame
    for i in range(frame, len(nuc_seq) - 2,3):
        codon = nuc_seq [i:i+3] #previously referenced codon reading frame
        aa = GENETIC_CODE.get(codon, "X")
        if aa == "STOP":
            break
        aas.append(aa) #append the amino acids
    return "".join(aas) #join together the retrieved amino acids 



class SequenceRecord:
    """ 
    Stores the information about a given accession number that is pulled from the FASTA file in my dictionary. 

    Handles one stored sequence at a time, parsing for its header and associated sequence, recording the associate accession ID, virus name, segment labels, and sequences. 

    accession: string value that is taken from the user and parsed through the dictionary until found, assocaited with ">"

    virus_name: string value that is found [:-2] in the header line because it omits "segment __", if not found it will default to UNKNOWN

    segment: string value that is the very last value (S, M, or L) indicating the size of the sequence that is associated with that accession ID

    sequence: string value that is below the header line and is the nucleotides that will be stored to translate and form the analysis needed

    __repr__: reutrns the record of the string, reads the accession number and the segment 

    length: the count of how many nucleotides are found in the retrieved FASTA sequence

    from_header_and_seq: parses the header and sequence into the components needed to be nice and pretty and readable. Handles the dropping of ">", joining the lines together, and deriving and the important pieces of information from the dictionary 
    """

    def __init__(
        self, 
        accession: str,
        virus_name: str, 
        segment: str,
        gene_name: Optional[str],
        sequence: str,
    ):
        self.accession = accession 
        self.virus_name = virus_name
        self.segment = segment
        self.gene_name = gene_name 
        self.sequence = sequence 

    def __repr__(self):
        return f"SequenceRecord(accession={self.accession!r}, segment={self.segment!r})"

    def length(self):
        return len(self.sequence) #counts the number of bases in the retrieved sequence 
        
    
    def from_header_and_seq(header: str, sequence: str) -> SequenceRecord:
        parts = header.split()

        accession = parts[0]
        virus_name = " ".join(parts[1:-2]) if len(parts) > 3 else "UNKNOWN_VIRUS"
        segment = parts[-1] if len(parts) >= 2 else "UNKNOWN"
        clean_seq = sequence.replace(" ","").replace("\n","") #Again, clean up the sequence so it has minimal errors with the input and reference 

        return SequenceRecord(
            accession = accession, 
            virus_name = virus_name,
            segment = segment, 
            gene_name = None,
            sequence = clean_seq,
        )
        

class SequenceDatabase:
    """ 
    This class will read the FASTA file and provide lookup methods and iterations.

    fasta_path: library path ; maps the accession ID but defualts to data/sequences.fasta

    records_by_accession: map the IDs to corresponding sequences according to the chosen library 

    get_by_accession: returns the sequence retrieved by the accession 

    Key components of the library that are required for the class to work: 
        - Headers begin with ">"
        -Sequence lines are stripped of any accidental spaces and are joined together
        """
    
    def __init__(self, fasta_path: str | Path):
        self.fasta_path = Path(fasta_path)
        self.records_by_accession: Dict[str, SequenceRecord] = {}

    def load(self):
        if not self.fasta_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {self.fasta_path}") #If the FASTA doesn't exist then it will return the f"
        header: Optional[str] = None
        seq_lines: list[str] = []

        with self.fasta_path.open() as fh:
            for raw_line in fh:
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(">"): #I made the dictionary based off a file of data given and cross referencing to NCBI database and each accession ID begins with '>'
                    if header is not None:
                        sequence = "\n".join(seq_lines)
                        record = SequenceRecord.from_header_and_seq(header, sequence)
                        self.records_by_accession[record.accession] = record

                    header = line[1:] #Skips the >
                    seq_lines = [] #Empty list 
                else:
                    seq_lines.append(line)

        if header is not None:
            sequence = "\n".join(seq_lines) #Stiches the lines together 
            record = SequenceRecord.from_header_and_seq(header, sequence)
            self.records_by_accession[record.accession] = record

    def get_by_accession(self, accession: str):
        return self.records_by_accession[accession]

    def __len__(self):
        return len(self.records_by_accession) #the count of the sequence based of the accession record 

    def accessions(self):
        return self.records_by_accession.keys() #the  key in the dictionary 

    def records(self):
        return self.records_by_accession.values() #the value associated with the key in the dictionary 
                


class AlphaFoldManager: 
    """ Class will determine whether the sequence is under or over the maximum limit set by ColabFold, that being 1400 aa, and either:
    a.) output a FASTA or JSON file that can be copy and pasted directly into AlphaFold3
    b.) submit a translated FASTA file from my premade data file that will output downloadable files 
    Features:
    - Writes FASTA file for a SequenceRecord
    - Handles long sequences without submitting an error and be helpful
    -Calls 'colabfold_batch' using subprocess 
    -Has a set dictinary for every outcome
    """
    
    def __init__(self, max_auto_len: int = 1400, results_dir: str | Path= "colabfold_results", colabfold_exe: str = "colabfold_batch"): #The max amino acid sequence for ColabFold is 1400 amino acids, could work with 1500 but becomes inaccurate and struggles
        self.max_auto_len = max_auto_len #length threshold for ColabFold
        self.results_dir = Path(results_dir) #output directory 
        self.colabfold_exe = colabfold_exe #colabfold_batch executable; ubuntu
        self.results_dir.mkdir(parents = True, exist_ok = True) #ensures that the results actually exist

    def make_fasta(self, record: "SequenceRecord", region: Optional[Tuple[int, int]] = None, description: str = "", seq_override: Optional[str] = None):
        base_seq = seq_override if seq_override is not None else record.sequence
        
        if region is None:
            seq = base_seq 
            id_part = record.accession
            region_info = f"full_length = {record.length()}"
        else: 
            start, end = region 
            seq = base_seq[start:end]
            id_part = f"{record.accession}_{start+1}_{end}"
            region_info = f"residues={start+1}-{end}"

        desc = [f"virus: {record.virus_name.replace(' ','_')}", f"segment = {record.segment}", region_info] #Building the header components ie. the virus name and the segment (S, M, L) and the length of the segment

        if description: 
            desc.append(description)

        header = ">" + id_part + " " + " ".join(desc) #In the dictionary I created the accession number is begun with '>' and then followed with the virus name and segment and then in a separate line the FASTA sequence, so this outputs somewhat the same setup 

        lines = [header] 
        for i in range(0, len(seq), 70): #The number of nucleotides across was 70, so I am keeping the same structure from the NCIB database to keep a more constant output and input 
            lines.append(seq[i:i+70])

        return "\n".join(lines) #Joins the lines together

    def run_or_prepare(self, record: "SequenceRecord") -> Dict[str, object]:
        aa_full = translate_nuc_to_aa(record.sequence)
        
        if len(aa_full) > self.max_auto_len:  #Handles if the sequence cannot be registered into ColabFold and instead needs to output a JSON or copy and pasteable file for AlphaFold3 
        
            aa_region = aa_full[:self.max_auto_len]
            region = (0, len(aa_region))
            fasta_str = self.make_fasta(record, region=region, description="AlphaFold3_manual_input_protein", seq_override=aa_full,
            )

            json_path = self.results_dir / f"{record.accession}_af3_manual.json" #make a JSON file for input into AlphaFold 3, makes the copy and pastable process faster 
            json_data = {                           #The data needed in the JSON file 
                "mode": "manual_af3",
                "accession": record.accession,
                "virus_name": record.virus_name,
                "aa_length_full": len(aa_full),
                "aa_length_used": len(aa_region), 
                "max_auto_len": self.max_auto_len, 
                "fasta": fasta_str,
            }
            json_path.write_text(json.dumps(json_data, indent=2))

            manual_fasta_path = self.results_dir / f"{record.accession}_af3_manual.fasta"
            manual_fasta_path.write_text(fasta_str)
            
            msg_lines = [ 
                f"Sequence {record.accession} translated to protein is {len(aa_full)} aa,", 
                f"which exceeds the automated limit of {self.max_auto_len} aa.",
                "",
                "Copy and Paste the FASTA below into AlphaFold3 web:",
                "",
                fasta_str,
                "",
                f"(Metadata + FASTA also saved to {json_path})", 
            ]
            
            return {
                "mode": "manual_af3",
                "message": "\n".join(msg_lines),
                "fasta": fasta_str,
                "fasta_path": None,
                "results_dir": None,
                "json_path": json_path,
                "aa_sequence": aa_full,
            }
#The second case in which it is acceptable for ColabFold and is directly run into it
        fasta_path = self.results_dir / f"{record.accession}.fasta" #Creates the path 
        fasta_str = self.make_fasta(record, seq_override=aa_full) #Generates and writes the FASTA file
        fasta_path.write_text(fasta_str)
#This is the subprocess command, the call for ColabFold

        accession_outdir = self.results_dir / record.accession
        accession_outdir.mkdir(parents=True, exist_ok=True)
        
        cmd = [
            self.colabfold_exe, 
            str(fasta_path), 
            str(accession_outdir), #output directory 
        ]

        print("[AlphaFoldManager] Running:", " ".join(cmd))

        try: 
            completed = subprocess.run(cmd, check = False)

            if completed.returncode != 0:
                msg = (
                    "WARNING: colabfold_batch was found but returned a non-zero exit code"
                    f"({completed.returncode}).\n"
                    "Falling back to manual AlphaFold3 FASTA/JSON output.\n"
                )
            msg = (
                f"Ran ColabFold on {record.accession} ({len(aa_full)} aa).\n"
                f"FASTA file: {fasta_path}\n"
                f"Results directory: {self.results_dir}"
            )
            return {
                "mode": "colabfold_run",
                "message": msg,
                "fasta": fasta_str,
                "fasta_path": fasta_path, 
                "results_dir": self.results_dir, 
            }
        except FileNotFoundError: #Arises when colabfold_batch is not installed, so will not work on a laptop or system unless the program is downloaded 
            msg = (
                "ERROR: 'colabfold_batch' executable was not found.\n"
            )
            print(msg)
            return manual_af3_output()
    
def run_tools( #Runs all the tools created above
    fasta_path: str, 
    accession: Optional[str] = None, 
    max_auto_len: int = 1400, 
    colabfold_exe: str = "colabfold_batch",  #Again the sequence maximum is 1400 
):
    db = SequenceDatabase("data/sequences.fasta") #uploads the data file 
    db.load()

    if accession is not None:
        if accession not in db.accessions():
            print(f"Accession {accession} was not found in sequences.")
            return 

        rec = db.get_by_accession(accession)
        af = AlphaFoldManager( 
            max_auto_len=max_auto_len, 
            results_dir="colabfold_results",
            colabfold_exe=colabfold_exe,
        )
        result = af.run_or_prepare(rec)
        print("Mode:", result["mode"])
        print(result["message"])
        return

    while True: #little loop 
        acc = input("Enter an accession ID or QUIT: ").strip() #prompted input 
        if acc.upper() == "QUIT":
            print("Bye!") #an out to quit the loop 
            break

        if not acc:
            print("Enter valid accession or QUIT.")
            continue

        if acc not in db.accessions():
            print(f"Accession{acc} not found. Please check for errors or NCBI database.")
            continue

        rec = db.get_by_accession(acc)
        af = AlphaFoldManager( 
            max_auto_len=max_auto_len, 
            results_dir="colabfold_results",
            colabfold_exe=colabfold_exe,
        )
        result = af.run_or_prepare(rec) #registers the run or preapre option, whether to output a JSON or a colabfold file 


        print("\nMode:", result["mode"])
        print(result["message"])
        print("-" * 70) #little breaker in between results 

    

def main(): 
    import argparse
    parser = argparse.ArgumentParser(description="AlphaTools / ColabFold helper") #Parsing for the various pieces of the argument 
    parser.add_argument(
        "--fasta",
        "-f",
        default = "data/sequences.fasta",
        help = "Path to nuclotide FASTA file",
    )
    parser.add_argument(
        "--accession",
        "-a",
        help = "Accession ID to process",
    )
    parser.add_argument(
        "--colabfold-exe",
        default = "colabfold_batch",
        help = "Name / path of the colabfold_batch executable", 
    )
    parser.add_argument(
        "--max-auto-len",
        type = int,
        default = 1400, 
        help = "Maximum amino acid length for ColabFold submission",
    )
    args = parser.parse_args()

    run_tools(
        fasta_path = args.fasta, 
        accession = args.accession,
        max_auto_len = args.max_auto_len, 
        colabfold_exe = args.colabfold_exe
    )

if __name__ == "__main__":
    main()

  
