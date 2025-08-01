from pymol import cmd
from .utils import track_and_log
from .utils import setup_full_logger

def main(
        work_dir:str,
        full_log_dir:str = 'E:/Coding/jupyter_root/ynby/pj3/repo/test/full_log.txt',
        rcsb_pdb_ids:str =  'E:/Coding/jupyter_root/ynby/pj3/repo/test/test.txt', 
        output_dir:str = 'E:/Coding/jupyter_root/ynby/pj3/repo/test/log.txt', 
        
):
    if full_log_dir:
        full_logger = setup_full_logger(full_log_dir)
    cmd.cd(work_dir)
    with open(rcsb_pdb_ids, 'r') as f:
        content = f.read()
        content = content.split(',')
        for id in content: 
                track_and_log(id, output_dir,full_log_dir=full_logger)