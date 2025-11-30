
import pandas as pd

# Input and Output paths
input_csv = r"D:\2026YJ\My_MR_Project\CMap_Results\Real_CMap_Top_Reversers_MOA.csv"
output_csv = r"D:\2026YJ\论文写作\spring模板\TableS5_CMap_Mechanisms.csv"

try:
    # Read the file
    df = pd.read_csv(input_csv)
    
    # Rename columns to match Table S5 style
    # Original: Mechanism, Connectivity Score, Example Drug
    df.columns = ['Mechanism', 'Connectivity_Score', 'Example_Drug']
    
    # Save
    df.to_csv(output_csv, index=False)
    print(f"Successfully generated Table S5 CSV at: {output_csv}")
    
except Exception as e:
    print(f"Error generating Table S5 CSV: {e}")
