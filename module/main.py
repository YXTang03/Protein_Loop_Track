from pymol import cmd
from .utils import track_and_log
from .utils import setup_full_logger
from tqdm import tqdm

def main(
        work_dir:str,
        rcsb_pdb_ids:str ,
        full_log_dir:str ,
        output_dir:str , 
        output_residual_index:bool = True
        
):
    
    if full_log_dir:
        full_logger = setup_full_logger(full_log_dir)
    cmd.cd(work_dir)
    with open(rcsb_pdb_ids, 'r') as f:
        content = f.read()
        content = content.split(',')
        pbar = tqdm(
            content, 
            colour='CYAN',
            smoothing= 1.0
        )
        for id in pbar: 
            pbar.set_description(f"Processing PDB: {id}" )
            track_and_log(
                id, 
                output_dir,
                full_log_dir=full_logger, 
                output_residual_index = output_residual_index
            )
