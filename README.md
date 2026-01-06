# alphafold-viral-protein-analysis
A python based pipeline that automates sequence retrieval from a given library, outputting a translated sequence, structural analysis, pLDDT confidence charts, and formatted files for 3D imaging through AlphaFold/ColabFold. 

## Example Output

###Predicted Protein Structure 
[OQ207703 Antu virus isolate YB_tick_2021_24 segment L] <img width="1064" height="865" alt="Screenshot 2025-12-11 081559" src="https://github.com/user-attachments/assets/b005bbd9-d2a4-41b5-8ae7-4e78e6531488" />
<img width="845" height="731" alt="Screenshot 2025-12-11 081032" src="https://github.com/user-attachments/assets/daa4fd35-9b0d-4493-918c-7ca629e7bc56" />
 ###pLDDT Confidence Plot
 High pLDDT values are indicated in blue correlatting to the confidence of certain regions of the predicted protein structure for the five given ranks. 
 <img width="2417" height="438" alt="OQ207703_virus__Antu_virus_isolate_YB_tick_2021_24_segment___L_full_length___12001_pae" src="https://github.com/user-attachments/assets/880c0961-7d97-41ee-8094-5e391e4f4a5e" />
<img width="1391" height="939" alt="OQ207703_virus__Antu_virus_isolate_YB_tick_2021_24_segment___L_full_length___12001_plddt" src="https://github.com/user-attachments/assets/08f9eaef-88ba-44b3-8a32-b72e69b513af" />

###Manual Strucutral Analysis 
--- Structural Summary ---
File: OQ207703
Total Residues: 4000
Peptide Preview: XXXXAHPNNTXNXXQRXAXTXXXXEQXXXQXQXSXXXRTXXTXHXDXXXRXXTXPXXXEQ...

Motifs:
N-glycosylation at position 3181
Phosphorylation at position 106
Phosphorylation at position 1751
Phosphorylation at position 1833
Phosphorylation at position 2557
Phosphorylation at position 2906
Phosphorylation at position 3178
Phosphorylation at position 3505
RNA-binding at position 246
RNA-binding at position 430
RNA-binding at position 546
RNA-binding at position 603
RNA-binding at position 747
RNA-binding at position 898
RNA-binding at position 1004
RNA-binding at position 1163
RNA-binding at position 1217
RNA-binding at position 1760
RNA-binding at position 1913
RNA-binding at position 1978
RNA-binding at position 2089
RNA-binding at position 2214
RNA-binding at position 2652
RNA-binding at position 3643
RNA-binding at position 3839
Protease cleavage sight at position 329
Protease cleavage sight at position 452
Protease cleavage sight at position 650
Protease cleavage sight at position 818
Protease cleavage sight at position 1131
Protease cleavage sight at position 2564
Protease cleavage sight at position 2974
Protease cleavage sight at position 3118
Protease cleavage sight at position 3139
Protease cleavage sight at position 3518
Protease cleavage sight at position 3523
Protease cleavage sight at position 3685

Secondary Structure Summary:
-- Alpha Helices: 740 residues
-- Beta-Sheets: 244 residues
-- Coils/Loops: 3016 residues

Disordered Regions:
None

Overall plDDT Confidence
-- High Confidence (>90) residues: 4000
-- Low Confidence (<50) residues: 0

Hydrophobic Core Residues (indices):
4, 6, 17, 54, 75, 76, 78, 91, 107, 111, 118, 122, 123, 129, 137, 139, 151, 157, 159, 171, 173, 178, 186, 198, 209, 238, 245, 248, 250, 251, 253, 254, 284, 306, 315, 318, 330, 345, 356, 359, 363, 375, 377, 412, 413, 420, 432, 437, 439, 453, 456, 457, 469, 473, 504, 505, 520, 529, 541, 542, 548, 575, 580, 591, 596, 611, 617, 618, 625, 629, 633, 647, 651, 658, 670, 680, 684, 697, 700, 717, 728, 734, 737, 764, 769, 770, 810, 819, 820, 860, 871, 875, 884, 905, 925, 928, 944, 945, 947, 953, 955, 965, 986, 993, 1011, 1027, 1035, 1051, 1055, 1059, 1071, 1082, 1085, 1099, 1106, 1117, 1119, 1126, 1132, 1138, 1149, 1159, 1165, 1187, 1190, 1198, 1213, 1222, 1228, 1244, 1246, 1267, 1268, 1298, 1300, 1305, 1315, 1337, 1372, 1378, 1383, 1385, 1387, 1397, 1405, 1416, 1425, 1460, 1472, 1479, 1483, 1491, 1492, 1499, 1500, 1502, 1522, 1526, 1532, 1550, 1569, 1572, 1581, 1592, 1602, 1608, 1650, 1657, 1660, 1664, 1671, 1680, 1681, 1689, 1691, 1703, 1711, 1715, 1717, 1726, 1728, 1730, 1735, 1744, 1748, 1752, 1772, 1777, 1826, 1830, 1832, 1834, 1837, 1846, 1847, 1855, 1856, 1867, 1869, 1881, 1888, 1904, 1914, 1920, 1925, 1932, 1933, 1955, 1958, 1983, 2014, 2019, 2031, 2048, 2055, 2063, 2090, 2101, 2117, 2120, 2132, 2141, 2154, 2156, 2162, 2171, 2210, 2211, 2239, 2254, 2263, 2265, 2267, 2270, 2282, 2283, 2286, 2317, 2338, 2342, 2344, 2374, 2380, 2388, 2392, 2394, 2398, 2400, 2401, 2402, 2405, 2410, 2421, 2423, 2447, 2454, 2456, 2472, 2482, 2500, 2502, 2512, 2549, 2551, 2552, 2558, 2565, 2580, 2585, 2623, 2646, 2680, 2691, 2702, 2707, 2710, 2715, 2724, 2726, 2736, 2737, 2742, 2755, 2765, 2775, 2783, 2785, 2794, 2804, 2807, 2820, 2821, 2835, 2840, 2842, 2844, 2845, 2859, 2866, 2885, 2887, 2894, 2898, 2902, 2907, 2918, 2921, 2936, 2967, 2975, 3014, 3016, 3019, 3020, 3023, 3025, 3032, 3041, 3044, 3048, 3054, 3062, 3080, 3091, 3101, 3107, 3109, 3119, 3126, 3132, 3140, 3146, 3150, 3158, 3168, 3179, 3190, 3200, 3201, 3210, 3217, 3226, 3232, 3234, 3235, 3246, 3249, 3251, 3265, 3277, 3278, 3284, 3302, 3313, 3318, 3324, 3344, 3351, 3359, 3364, 3373, 3377, 3383, 3392, 3395, 3396, 3397, 3404, 3405, 3417, 3427, 3428, 3436, 3441, 3451, 3463, 3467, 3481, 3491, 3493, 3506, 3519, 3524, 3527, 3533, 3543, 3548, 3551, 3553, 3581, 3591, 3624, 3627, 3628, 3634, 3645, 3650, 3686, 3689, 3698, 3705, 3707, 3713, 3719, 3743, 3744, 3747, 3748, 3749, 3754, 3773, 3777, 3786, 3795, 3806, 3814, 3819, 3824, 3831, 3835, 3840, 3853, 3874, 3882, 3895, 3898, 3905, 3909, 3936, 3944, 3988, 3990

--- Domain Analysis ---
Domains:
Domain 1 (Residues 1-4000)
========================================================================================================================

## Features
- Prompts the input of accession ID or "QUIT" initiating a loop
- Retrieves nucleotide sequence from library outputting virus name and segment
- Translates nucleotide sequence into correct amino acid sequence and count
- Runs manual ColabFold or AlphaFold3
- Generates JSON and PDB file for 3D imaging
- Computes and summarizes pLDDT confidence charts for the predicted ranks
- Constructs subfolder of given PNGs, JSONs, PDBs, and FASTAs titled with accession ID
- Uses biological principles to output structure analysis based on amino acid composition

## Tech Stack 
- Programming: Python 
- AlphaFold, ColabFold
- Jupyter, Jupyter Notebook, Linux(Ubuntu)
- FASTA, PDB, JSON, PNG

## How to Run
1. Copy and download pipeline
2. Pip install ColabFold
3. Run pipeline
   - cd ~/localcolabfold
     source localcolabfold/conda/etc/profile.d/conda.sh
     conda activate ~/localcolabfold/localcolabfold.colabfold-conda
     export JAX_PLATFORMS=cpu
     cd/mnt/c/Users/"user"/alphafold_project

   - cd/mnt/c/Users/"user"/alphafold_project
     rm -rf __pycache__

   - python combined_tools.py \
      --fasta data/sequences.fasta \
      --colabfold-exe/home/"user"/localcolabfold/localcolabfold/colabfold-conda/bin/colabfold_batch
     
