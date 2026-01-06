from __future__ import annotations #allows me to reference a class that has yet to be defined 

import argparse #python's standard library for argument parsing, so the code can read lines like how i run it in ubuntu; allows us to use it for the add argument section for the fasta and accession, it allows us to parse the input Accession ID

from alphafold_tools import SequenceDatabase, AlphaFoldManager #Importing the classes from my section that output the ideal analysis, ie. the methods on finding the accession ID, stripping it of '>' and finding the sequence lines; the alpha fold manager that decides whether to run or prepare, ie. the colabfold subprocess
from entireReader import Protein, outputFormatter #Imports from Oktawian's code that formats the readily available analysis on motif, disordered reagions, etc.

def analyze_record(rec, af_manager: AlphaFoldManager) -> None: #Rec is the record that is storing the information about the ID (accession, virus name, segment, and nucleotide sequence from alphafold_tools
    print("=" * 120)
    print(f"Accession: {rec.accession}") #outputs the recorded accession
    print(f"Virus: {rec.virus_name}") #output the recorded virus_name
    print(f"Segment: {rec.segment}") #outputs the recorded segment 
    print("-" * 120) #a border for clarity and pretty-ness

    af_result = af_manager.run_or_prepare(rec) or {} #the AlphaFold result calls upon the af_manager from alphafold_tools and it returns the dictionary for the mode, message, etc. that is stored in the af_manager module; the empty dictionary just ensures that if an accession number is not found then it won't completely crash in an external code runner 

    print("[AlphaFoldTools]")
    print(f"Mode: {af_result.get('mode', 'UNKNOWN')}") #returns whether the af_mananer resulted in a colab fold run or an UNKNOWN ID, yet to be placed into the dictionary
    message = af_result.get("message", "").strip() #references the message portrion for either result from the af_manager module
    if message:
        print(message)
    else:
        print("(No message returned from AlphaFoldManager)") #handles the case if no message or case was found, so again it doesn't crash or produce an error

    print("-" * 120) #a border



    protein = Protein(rec.accession, rec.sequence) #Protein was written in entireReader.py, sotring all the info like the translated aa, secondary structure pLDDT scores, etc.

    if hasattr(protein, "translate"): #"has attribute" translate, then it will translate the protein, loophole around the attribute errors that were arising 
        protein.translate()

    if hasattr(protein, "predictSecondaryStructure"): #same logic carries down for all hasattr
        protein.predictSecondaryStructure()

    if hasattr(protein, "analyze"):
        protein.analyze()

    if hasattr(protein, "plddtSummary"):
        try:
            high, low = protein.plddtSummary()
            print(f"  High pLDDT: {high}")
            print(f"  Low pLDDT:  {low}")
        except Exception as e:  
            print(f"(Could not compute pLDDT summary: {e})")

    print("-" * 120)
    print("[Sequence / Structure Analysis Report]") #the sequence name and then the whole analysis report the entireReader produce with a border
    
    try:
        formatter = outputFormatter(protein)
        report_text = formatter.buildReport()
    except Exception as e:  
        report_text = f"(Error while building formatted report: {e})"

    print(report_text) #this takes the build report and outputs it 
    print("=" * 120)
    print()

def run_tools(
    fasta_path: str, 
    accession: str | None = None, 
    max_auto_len: int = 1400, 
    colabfold_exe: str = "colabfold_batch",
) -> None: 
    db = SequenceDatabase(fasta_path) #i was previously hard coding the file to be my sequences.fasta data file, but fasta_path keeps it ambigious and allows for alteration of accession ID dictionaries without failing 
    db.load() #loads the dictionary 

    af_manager = AlphaFoldManager(
        max_auto_len = max_auto_len, 
        results_dir = "colabfold_results",
        colabfold_exe = colabfold_exe, 
    )

    if accession is not None: #If there is an input 
        if accession not in db.accessions(): #but the accession was not found in the database
            print(f"[ERROR] Accession {accession} was not found in {fasta_path}.") #print and error but allow for retyping incase of a typo or not realizing it wasn't in the dictionary yet; kind of a limiting factor is it can't connect to the NCBI database or any other premade infinite dictionary, so it's me typing and looking up the nucleotide sequences 
            return

        rec = db.get_by_accession(accession) #retrieves the accession and then analyzes it using the alphafold manager
        analyze_record(rec, af_manager)
        return

    while True: #sets the loop so multiple inputs can happen until user decides to quit
        acc = input("Enter an accession ID or QUIT: ").strip() # prompts and then strips it to look it up 

        if acc.upper() == "QUIT": #if the user wants to break, uppercases it so if they input quit it still works; prevents pickiness 
            print("Bye!") 
            break #breaks the loop 

        if not acc: #if it wasn't a valid ID or QUIT, then it prompts a change without automatically breaking 
            print("Please enter a valid Accession ID or QUIT.")
            continue
        
        if acc not in db.accessions():
            print(f"[ERROR] Accession {accession} was not found in {fasta_path}.")
            continue
        

        rec = db.get_by_accession(acc)
        analyze_record(rec, af_manager)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run AlphaFold/ColabFold + custom analysis for one or more accessions."
    )

    parser.add_argument(
        "--fasta",
        "-f",
        default="data/sequences.fasta",
        help="Path to nucleotide FASTA file (same format as used by alphafold_tools.py).",
    ) #can overraide the data sequence, while there was no hardcoded data/sequences.fasta and instead was fasta_path, so the user can input their out data file if necessary, if not it will default to that file

    parser.add_argument(
        "--accession",
        "-a",
        help=(
            "Accession ID to process. Can be given multiple times. "
            "If omitted, all accessions in the FASTA will be processed."
        ),
    ) #again allows the user to override with their accession ID

    parser.add_argument(
        "--colabfold-exe",
        default="colabfold_batch",
        help="Name/path of the colabfold_batch executable.",
    ) #when i run it in ubuntu i need to locate the .py in my user folders so this allows the for the path ie. "home/katie/localcolabfold/colabfold-conda...

    parser.add_argument(
        "--max-auto-len",
        type=int,
        default=1400,
        help="Maximum amino acid length allowed for automatic ColabFold submission.",
    ) #can change the max amino acids if need be, but defaults to colabfold's 1400 range

    args = parser.parse_args() #reads what the user input and begins to fill in all the args. below

    run_tools(
        fasta_path = args.fasta,
        accession = args.accession, 
        max_auto_len = args.max_auto_len, 
        colabfold_exe = args.colabfold_exe,
    )

if __name__ == "__main__":
    main()
