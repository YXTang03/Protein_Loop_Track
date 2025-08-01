import numpy as np
from pymol import cmd
import logging
from contextlib import contextmanager

def pep_pretreatment(chain_label):
    
    '''
    Peptide pretreatment:
    1. Remove atoms composing rings on residuals.
    2. Remove solvent and ligand molecules in current file.
    '''

    aa_with_residual_ring = ['HIS', 'PRO', 'PHE', 'TYR', 'TRP']
    for aa in aa_with_residual_ring:
        cmd.select(f'{aa.lower()}', f'chain {chain_label} and resn {aa} and not backbone')
        cmd.remove(f'{aa.lower()}')
        cmd.remove('solvent')
        cmd.remove('organic')

    global model, index_map, atom_indices, index_map_r
    model = cmd.get_model(f'chain {chain_label}')
    atom_indices = [atom.index for atom in model.atom]
    index_map = {atom_idx: i  for i, atom_idx in enumerate(atom_indices)}
    index_map_r = {i: atom_idx  for i, atom_idx in enumerate(atom_indices)}


def pdb_to_mtx():

    '''
    Convert the protein structural from the pdb file to the undirectional graph represented by the adjacency matrix.
    '''

    n = len(atom_indices)
    adj_matrix = np.full((n, n), float('inf'))
    
    for i in range(len(model.atom)):
        for bond in model.bond:
            bond_atom0 = bond.index[0]
            bond_atom1 = bond.index[1]
            if bond_atom0 == i:
                adj_matrix[bond_atom0, bond_atom1] = 1
                adj_matrix[bond_atom1, bond_atom0] = 1
    
    return adj_matrix, model.atom


def adj_matrix_to_list(matrix):

    '''
    Convert the adjacency matrix to the adjacency list.
    '''

    n = matrix.shape[0]
    adj_list = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(n):
            if matrix[i][j] == 1:
                adj_list[i].append(j)
    return adj_list


def pdb_to_adj_list():

    '''
    Alternatively, and more concisely, the adjacency list can be created from the protein structure in one step.
    '''

    adj_list = {}
    for i in range(len(model.atom)):
        if i not in adj_list:
            adj_list[i] = []
        for j in range(len(model.bond)):
            bond_atom0 = model.bond[j].index[0]
            bond_atom1 = model.bond[j].index[1]
            if bond_atom0 == i and bond_atom1 not in adj_list[i]:
                adj_list[i].append(bond_atom1)
            if bond_atom1 == i and bond_atom0 not in adj_list[i]:
                adj_list[i].append(bond_atom0)
    return adj_list


def selector(
        chain_label, 
        full_logger:logging.Logger,
        residual_name = 'CYS',
        element_name = 'sg'
):
    
    '''
    Return the indexes of sulphur atoms in disulfide bonds

    In PyMol, the outcomes of this function **equivalent** to that of operation: 
    
    *L* -> *atom identifiers* -> *index* 
    
    after the selection of sulphur atoms on disulfide bonds by 

    *cmd.select("chain A and resn CYS and not backbone and name sg")*

    '''

    
    start_cys = cmd.get_model(f'chain {chain_label} and resn {residual_name} and not backbone and name {element_name}')
    dslf = [atom.index for atom in start_cys.atom]
    start = []
    for i, atom_1 in enumerate(dslf):
        for j, atom_2 in enumerate(dslf):
            bond_exists = any((i == a and j == b) for bond in start_cys.bond for a, b in [bond.index])
            if bond_exists:
                start.append(index_map[atom_2])
                if full_logger:
                    full_logger.info(f'Atom index: {atom_1} bonds with Atom index: {atom_2}')
    return start


def is_loop(adj_list, start):

    '''
    Discriminate the existence of rings containing disulfide bonds.

    Start searching with the sulphur atoms in disulfide bonds. 

        -> Parameter *start* = indexes of sulphur atoms in disulfide bonds
    '''

    print(f'start with {model.atom[start].name}')
    
    visited = set()
    found_cycle = [False]

    def dfs(node, parent):
        visited.add(node)
        for neighbor in adj_list[node]:
            if neighbor == parent:
                continue
            if neighbor == start and parent != start:
                found_cycle[0] = True
                return
            if neighbor not in visited:
                dfs(neighbor, node)

    dfs(start, None)
    
    return found_cycle[0]


def track_loop(adj_list:dict, start:int,full_logger:logging.Logger):
    '''
    Tracking the indexes in the adjacency matrix or list, given that there is a loop starting with a sulphur atom in the disulphide bond
    '''


    global cycle_path
    visited = set()
    parent = {}
    cycle_path = []

    
    if full_logger:
        full_logger.info(f'start with {model.atom[start].name}')


    def dfs(node, par):
        visited.add(node)
        parent[node] = par

        for neighbor in adj_list[node]:
            if neighbor == par:
                continue
            if neighbor in visited:
                path1 = []
                path2 = []
                x = node
                while x is not None:
                    path1.append(x)
                    if x == neighbor:
                        break
                    x = parent[x]

                y = neighbor
                while y not in path1:
                    path2.append(y)
                    y = parent[y]

                #Change in 01/08/2025 to fix the problem that miss selction in the last second residual
                path = path1[:path1.index(y)+1] + path2[::-1]
                #path = path1 + path2

                if start in path:
                    cycle_path.extend(path + [start]) 
                    return True
            else:
                if dfs(neighbor, node):
                    return True
        return False

    dfs(start, None)
    
    return cycle_path if cycle_path else None


def transform_ids(
        ids:list, output_residual_index = True
):
    '''
    Transform ids to 

    **either** 

    PyMOL residual indexes (prefered option)

    **or**

    PyMOL atom indexes (not recommond)
    '''
    if output_residual_index:
        loop_idx = list(dict.fromkeys([model.atom[i].resi for i in ids]))
    else:
        loop_idx = list(dict.fromkeys([index_map_r[i] for i in ids]))
    output = ''
    for id in loop_idx[:-1]:
        output += str(id) + '+'
    output += str(loop_idx[-1])
    return output


@contextmanager
def open_log_files(output_log_path, full_log_path):
    try: 
        output_log = open(output_log_path, 'a')
        if full_log:
            full_log = open(full_log_path, 'a')
        yield output_log, full_log

    finally:
        if output_log:
            output_log.close()
        if full_log:
            full_log.close()


def setup_full_logger(full_log):
    logger = logging.getLogger(f'full_log')
    logger.setLevel(logging.INFO)

    if not logger.handlers:
        handler = logging.FileHandler(full_log, mode='a', encoding='utf-8')
        
        logger.addHandler(handler)

    return logger


def track_and_log(
        pdb_id:str, 
        output_dir,
        full_log_dir,
        chain_label = None, 
        need_adj_matrix: bool = False, 
        output_residual_index = True
):
    


    if full_log_dir:
        full_logger = setup_full_logger(full_log_dir)

    try:
        try:
            cmd.fetch(pdb_id)
        except:
            pass        
        chain_labels = cmd.get_chains()



        for chain_label in chain_labels:
            if full_logger:
                full_logger.info(f"{'*'*20} Processing {pdb_id}, chain {chain_label} {'*'*20}\n")

            pep_pretreatment(chain_label)
            adj_list = pdb_to_adj_list()


            if need_adj_matrix and full_logger:
                full_logger.info(f"{'-'*10}Adjacency Matrix Info:{'-'*10}")
                adj_matrix, atoms = pdb_to_mtx(chain_label)
                full_logger.info(f"Adjacency Matrix Shape: {adj_matrix.shape}\nFirst 5 Atom in Current Peptide: ")
                for i in range(5):
                    full_logger.info(f"{i}: {atoms[i].resn} {atoms[i].name} ({atoms[i].resi})\n")
                    
            if full_logger:
                full_logger.info(f"{'-'*10}Disulphide Bonds Form between Atoms:{'-'*10}\n")
                
            dslf = selector(chain_label, full_logger)
            
            for atom_sg in dslf:
                if full_logger:
                    full_logger.info(f"{'-'*10}Following Info is Compatible with PyMOL cmd.selector:{'-'*10}\n")
                    full_logger.info(f"{'-'*10}Tracking Atom Indexes on the Loop:{'-'*10}\n")
                try:
                    ids = track_loop(adj_list, atom_sg, full_logger)
                    ids.sort()
                    select_index = transform_ids(ids, output_residual_index)
                    
                    if full_logger:
                        full_logger.info(f"{select_index}\n")


                    with open(output_dir, 'a') as log:
                        log.write(f'"{pdb_id}", {chain_label}, {select_index}\n')
                except:
                    if full_logger:
                        full_logger.info('Current Disulphide Bond Does Not in a Loop\n')
                    
    finally:
        cmd.delete(f'{pdb_id}')