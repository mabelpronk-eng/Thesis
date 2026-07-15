##############################################
# Same function as extract_signature.py script
##############################################
# Save signature for later use  
# Concatenate along columns (genes aligned by index) 
#
# Script can be used if snakemake is not used 
##############################################
import pandas as pd
import os
import sys


THIS_DIR = os.path.dirname(__file__)
SRC_DIR  = os.path.abspath(os.path.join(THIS_DIR, '..','..', 'src'))
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

#import Statescope.Statescope as scope
from Statescope.Statescope import Initialize_Statescope, Statescope

Statescope_model = Statescope.load('/net/beegfs/users/P086608/StatescopePro/tutorial/pseudobulk_level2/model_initialized.pkl')
signature = pd.concat([Statescope_model.scExp, Statescope_model.scVar], axis=1) 
 
# # # Add gene column from index  
signature["Gene"] = signature.index 

# # # Add logical marker column  
signature["IsMarker"] = signature["Gene"].isin(Statescope_model.Markers) 

# # # Reorder columns  
cols = ["Gene", "IsMarker"] + [c for c in signature.columns if c not in ["Gene", "IsMarker"]] 

signature = signature[cols] 

print(signature)
# # ## Write  
output_path = "/net/beegfs/users/P086608/StatescopePro/tutorial/pseudobulk_level2/Signature_16celltypes.txt" 

#  
signature.to_csv(output_path, sep="\t", index=False )
