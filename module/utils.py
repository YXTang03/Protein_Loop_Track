import os
import numpy as np
from pymol import cmd
import logging
import tqdm
from pdbfixer import PDBFixer
from openmm.app.pdbfile import PDBFile

class PepFix:
    """
    # Fix missing residuals in the protein structure analyzed currently.
    """
    
    def __init__(self, work_dir: str):
        self.work_dir = work_dir
    
    def fix_pdb(
            self, 
            pdb_id: str, 
            missing_structure_log_dir:str,
            missing_structure_log_name:str,
            fix_residual:bool = False,
            fix_atom:bool = True,
            fix_threshold: int = 10, 
    ) -> None:
        
        missing_structure_logger = Logger(
            missing_structure_log_dir, 
            missing_structure_log_name
        )
        missing_structure_log = missing_structure_logger.get_logger()
        
        file_path = os.path.join(self.work_dir, f"{pdb_id}.cif")
        #f'{self.work_dir}/{pdb_id}.cif'
        
        fixer = PDBFixer(file_path)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        
        all_missing_atom = fixer.missingAtoms
        all_missing_res = fixer.missingResidues
        if missing_structure_log:
            if all_missing_atom:
                missing_structure_log.info(f'{pdb_id}: \nMissing Atom{all_missing_atom}\n')
            if all_missing_res:
                missing_structure_log.info(f'{pdb_id}: \nMissing Residual{all_missing_res}\n')
        if fix_residual is False:
            all_missing_res = {}
        if fix_atom is False:
            all_missing_atom = {}

        for missing_res in all_missing_res.values():
            if len(missing_res) <= fix_threshold:
                fixer.addMissingAtoms()
                output_file_name = os.path.join(self.work_dir, f"{pdb_id}.pdb")
                with open(output_file_name, 'w') as f:
                    PDBFile.writeFile(
                        fixer.topology, 
                        fixer.positions, 
                        f
                    )


class ProteinStructure:
    """蛋白质结构数据的封装类"""
    
    def __init__(self, chain_label: str):
        self.chain_label = chain_label
        self.model = None
        self.atom_indices = []
        self.index_map = {}
        self.index_map_r = {}
        self.adj_list = {}
        
    def preprocess_peptide(self) -> None:
        """
        Peptide pretreatment:
        1. Remove atoms composing rings on residuals.
        2. Remove solvent and ligand molecules in current file.
        """
        aa_with_residual_ring = ['HIS', 'PRO', 'PHE', 'TYR', 'TRP']
        
        for aa in aa_with_residual_ring:
            cmd.select(f'{aa.lower()}', f'chain {self.chain_label} and resn {aa} and not backbone')
            cmd.remove(f'{aa.lower()}')
        
        cmd.remove('solvent')
        cmd.remove('organic')
        
        # 获取模型数据
        self.model = cmd.get_model(f'chain {self.chain_label}')
        self.atom_indices = [atom.index for atom in self.model.atom]
        self.index_map = {atom_idx: i for i, atom_idx in enumerate(self.atom_indices)}
        self.index_map_r = {i: atom_idx for i, atom_idx in enumerate(self.atom_indices)}
    
    def build_adjacency_matrix(self) -> tuple[np.ndarray, list[any]]:
        """Convert the protein structural from the pdb file to the undirectional graph 
        represented by the adjacency matrix."""
        n = len(self.atom_indices)
        adj_matrix = np.full((n, n), float('inf'))
        
        for i in range(len(self.model.atom)):
            for bond in self.model.bond:
                bond_atom0 = bond.index[0]
                bond_atom1 = bond.index[1]
                if bond_atom0 == i:
                    adj_matrix[bond_atom0, bond_atom1] = 1
                    adj_matrix[bond_atom1, bond_atom0] = 1
        
        return adj_matrix, self.model.atom
    
    def build_adjacency_list(self) -> dict[int, list[int]]:
        """
        Convert the adjacency matrix to the adjacency list.
        """
        self.adj_list = {}
        
        for i in range(len(self.model.atom)):
            if i not in self.adj_list:
                self.adj_list[i] = []
            
            for j in range(len(self.model.bond)):
                bond_atom0 = self.model.bond[j].index[0]
                bond_atom1 = self.model.bond[j].index[1]
                
                if bond_atom0 == i and bond_atom1 not in self.adj_list[i]:
                    self.adj_list[i].append(bond_atom1)
                if bond_atom1 == i and bond_atom0 not in self.adj_list[i]:
                    self.adj_list[i].append(bond_atom0)
        
        return self.adj_list
    
    @staticmethod
    def matrix_to_list(matrix: np.ndarray) -> dict[int, list[int]]:
        """将邻接矩阵转换为邻接表"""
        n = matrix.shape[0]
        adj_list = {i: [] for i in range(n)}
        
        for i in range(n):
            for j in range(n):
                if matrix[i][j] == 1:
                    adj_list[i].append(j)
        
        return adj_list


class DisulfideBondAnalyzer:
    """二硫键分析器"""
    
    def __init__(
            self, 
            protein_structure: ProteinStructure, 
            logger:logging.Logger = None
    ):
        self.protein = protein_structure
        self.logger = logger
        
    def find_disulfide_atoms(self, residual_name: str = 'CYS', element_name: str = 'sg') -> list[int]:
        """
        返回二硫键中硫原子的索引
        等效于PyMol操作: cmd.select("chain A and resn CYS and not backbone and name sg")
        """
        start_cys = cmd.get_model(
            f'chain {self.protein.chain_label} and resn {residual_name} and not backbone and name {element_name}'
        )
        dslf = [atom.index for atom in start_cys.atom]
        start = []
        
        for i, atom_1 in enumerate(dslf):
            for j, atom_2 in enumerate(dslf):
                bond_exists = any(
                    (i == a and j == b) for bond in start_cys.bond for a, b in [bond.index]
                )
                if bond_exists:
                    start.append(self.protein.index_map[atom_2])
                    if self.logger:
                        self.logger.info(f'Atom index: {atom_1} bonds with Atom index: {atom_2}')
        
        return start
    
    def has_loop(self, start: int) -> bool:
        """判断是否存在包含二硫键的环结构"""
        if self.logger:
            self.logger.info(f'start with {self.protein.model.atom[start].name}')
        
        visited = set()
        found_cycle = [False]
        
        def dfs(node: int, parent: int) -> None:
            visited.add(node)
            for neighbor in self.protein.adj_list[node]:
                if neighbor == parent:
                    continue
                if neighbor == start and parent != start:
                    found_cycle[0] = True
                    return
                if neighbor not in visited:
                    dfs(neighbor, node)
        
        dfs(start, None)
        return found_cycle[0]
    
    def track_loop_path(self, start: int) :# -> None | list:
        """跟踪环路径，返回环中所有原子的索引"""
        visited = set()
        parent = {}
        cycle_path = []
        
        if self.logger:
            self.logger.info(f'start with {self.protein.model.atom[start].name}')
        
        def dfs(node: int, par: int) -> bool:
            visited.add(node)
            parent[node] = par
            
            for neighbor in self.protein.adj_list[node]:
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
                    
                    path = path1[:path1.index(y)+1] + path2[::-1]
                    
                    if start in path:
                        cycle_path.extend(path + [start])
                        return True
                else:
                    if dfs(neighbor, node):
                        return True
            return False
        
        dfs(start, None)
        return cycle_path if cycle_path else None


class IndexTransformer:
    """索引转换器"""
    
    def __init__(self, protein_structure: ProteinStructure):
        self.protein = protein_structure
    
    def transform_ids(self, ids: list[int], output_residual_index: bool = True) -> str:
        """
        转换索引为PyMOL残基索引（推荐）或原子索引
        """
        if output_residual_index:
            loop_idx = list(dict.fromkeys([self.protein.model.atom[i].resi for i in ids]))
        else:
            loop_idx = list(dict.fromkeys([self.protein.index_map_r[i] for i in ids]))
        
        if not loop_idx:
            return ""
        
        output = ''
        for id in loop_idx[:-1]:
            output += str(id) + '+'
        output += str(loop_idx[-1])
        
        return output


class Logger:
    """Create Log. To initialize a log, parameter path and name are needed."""
    
    def __init__(
            self, 
            log_path: str = None, 
            log_name: str = None
    ):
        self.full_log_path = log_path
        self.log_name = log_name
        self.logger = None
        
        if log_path:
            self.logger = self._setup_logger()
    
    def _setup_logger(self) -> logging.Logger:
        """设置完整日志记录器"""
        logger = logging.getLogger(self.log_name)
        logger.setLevel(logging.INFO)
        
        if not logger.handlers:
            handler = logging.FileHandler(self.full_log_path, mode='a', encoding='utf-8')
            logger.addHandler(handler)
        
        return logger
    
    def get_logger(self) -> logging.Logger:
        """获取日志记录器"""
        return self.logger


class ProteinAnalyzer:
    """主分析器类，协调各个组件完成蛋白质分析"""
    
    def __init__(self, work_dir: str):
        self.work_dir = work_dir
        self.pdb_fixer = PepFix(work_dir)
        
    def analyze_protein(
        self, 
        pdb_id: str, 
        output_dir: str,
        full_log_name: str,
        full_log_dir: str = None,
        need_adj_matrix: bool = False,
        output_residual_index: bool = True
    ) -> None:
        """分析单个蛋白质的二硫键环路"""
        
        # 设置日志
        logger_manager = Logger(
            full_log_dir, 
            full_log_name
        )
        full_logger = logger_manager.get_logger()
        
        try:
            
            try:
                cmd.load(os.path.join(self.work_dir, f"{pdb_id}.pdb"))
            except:
                cmd.load(os.path.join(self.work_dir, f"{pdb_id}.cif"))
            
            chain_labels = cmd.get_chains()
            
            for chain_label in chain_labels:
                if full_logger:
                    full_logger.info(f"{'*'*20} Processing {pdb_id}, chain {chain_label} {'*'*20}\n")
                
                
                protein = ProteinStructure(chain_label)
                protein.preprocess_peptide()
                protein.build_adjacency_list()
                
                
                if need_adj_matrix and full_logger:
                    full_logger.info(f"{'-'*10}Adjacency Matrix Info:{'-'*10}")
                    adj_matrix, atoms = protein.build_adjacency_matrix()
                    full_logger.info(f"Adjacency Matrix Shape: {adj_matrix.shape}\nFirst 5 Atom in Current Peptide: ")
                    
                    for i in range(min(5, len(atoms))):
                        full_logger.info(f"{i}: {atoms[i].resn} {atoms[i].name} ({atoms[i].resi})\n")
                
                # 分析二硫键
                analyzer = DisulfideBondAnalyzer(protein, full_logger)
                transformer = IndexTransformer(protein)
                
                if full_logger:
                    full_logger.info(f"{'-'*10}Disulphide Bonds Form between Atoms:{'-'*10}\n")
                
                dslf = analyzer.find_disulfide_atoms()
                
                for atom_sg in dslf:
                    if full_logger:
                        full_logger.info(f"{'-'*10}Following Info is Compatible with PyMOL cmd.selector:{'-'*10}\n")
                        full_logger.info(f"{'-'*10}Tracking Atom Indexes on the Loop:{'-'*10}\n")
                    
                    try:
                        ids = analyzer.track_loop_path(atom_sg)
                        if ids:
                            ids.sort()
                            select_index = transformer.transform_ids(ids, output_residual_index)
                            
                            if full_logger:
                                full_logger.info(f"{select_index}\n")
                            
                            with open(output_dir, 'a') as log:
                                log.write(f'"{pdb_id}", {chain_label}, {select_index}\n')
                    except Exception as e:
                        if full_logger:
                            full_logger.info(f'Current Disulphide Bond Does Not in a Loop: {str(e)}\n')
        
        finally:
            cmd.delete(f'{pdb_id}')
    
    def batch_analyze(
        self,
        rcsb_pdb_ids: str,
        output_dir: str,
        missing_structure_log_name,
        full_log_name,
        missing_structure_log_dir: str = None,
        full_log_dir: str = None,
        fix_threshold: int = 10
    ) -> None:
        """批量分析蛋白质"""
        
        cmd.cd(self.work_dir)
        
        with open(rcsb_pdb_ids, 'r') as f:
            content = f.read()
            pdb_ids = content.split(',')

            pbar = tqdm(
            pdb_ids, 
            colour='CYAN',
            smoothing= 1.0
        )
            
            for pdb_id in pbar:
                pdb_id = pdb_id.strip()
                if pdb_id:
                    # 修复PDB文件
                    self.pdb_fixer.fix_pdb(
                        pdb_id, 
                        missing_structure_log_dir,
                        missing_structure_log_name,
                        fix_residual=False, 
                        fix_atom=True,
                        fix_threshold=fix_threshold                  
                    )
                    
                    # 分析蛋白质
                    self.analyze_protein(
                        pdb_id, 
                        output_dir, 
                        full_log_name,
                        full_log_dir=full_log_dir
                    )