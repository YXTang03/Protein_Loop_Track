from .utils import *

def main(
        work_dir:str, 
        rcsb_pdb_ids:str, 
        output_dir:str, 
        missing_structure_log_name:str, 
        full_log_name:str, 
        missing_structure_log_dir:str, 
        full_log_dir:str, 
        fix_threshold:int=10      
):
    analyzer = ProteinAnalyzer(work_dir)
    analyzer.batch_analyze(
        rcsb_pdb_ids=rcsb_pdb_ids,
        output_dir=output_dir,
        missing_structure_log_name= missing_structure_log_name,
        full_log_name= full_log_name,
        missing_structure_log_dir=missing_structure_log_dir,
        full_log_dir=full_log_dir,
        fix_threshold=fix_threshold
    )


