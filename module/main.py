from pymol import cmd

from .utils import pep_pretreatment
from .utils import pdb_to_mtx
from .utils import adj_matrix_to_list
from .utils import pdb_to_adj_list
from .utils import selector
from .utils import is_loop
from .utils import track_loop
from .utils import transform_ids
from .utils import track_and_log

def main(
        work_dir:str,
        rcsb_pdb_ids:str =  'E:/Coding/jupyter_root/ynby/pj3/repo/test/test.txt', 
        output_dir:str = 'E:/Coding/jupyter_root/ynby/pj3/repo/test/log.txt'
):
    cmd.cd(work_dir)
    with open(rcsb_pdb_ids, 'r') as f:
        content = f.read()
        content = content.split(',')
        for id in content: 
                track_and_log(id, output_dir)