"""
Polypeptide PDB cleaning script for OpenFF compatibility.

Cleans VMD-generated polypeptide PDB files by:
- Fixing non-standard atom names (ASN: N0→ND2, HN21→HD21, HN22→HD22)
- Fixing C-terminal atom names (ALA: OD1→O, OD2→OXT)
- Converting NH2 capping group to standard COO- terminus (OXT)
- Renaming N-terminus HN→H1 and adding missing H2, H3
- Removing TB/lanthanide ions (handled separately by project.py)
- Sorting residues and renumbering atoms

Usage:
    python unscrew_polypeptide.py LBT3-.pdb
    python unscrew_polypeptide.py LBT5-.pdb
    python unscrew_polypeptide.py --all  # Process all .pdb files in directory

Output:
    {input_stem}_cleanedPDB.pdb
"""

import argparse
import numpy as np
from pathlib import Path


CLEANED_PDB_SUFFIX = "_cleanedPDB"

# Known lanthanide/metal residue names to skip
SKIP_RESIDUES = {'TB', 'ND', 'LA', 'CE', 'PR', 'SM', 'EU', 'GD', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU'}


def generate_nterm_hydrogens(n_pos, h1_pos, ca_pos):
    """Generate H2, H3 positions for NH3+ terminus using tetrahedral geometry."""
    n_to_h1 = h1_pos - n_pos
    n_to_h1 = n_to_h1 / np.linalg.norm(n_to_h1)

    n_to_ca = ca_pos - n_pos
    n_to_ca = n_to_ca / np.linalg.norm(n_to_ca)

    perp1 = np.cross(n_to_h1, n_to_ca)
    perp1 = perp1 / np.linalg.norm(perp1)
    perp2 = np.cross(n_to_h1, perp1)
    perp2 = perp2 / np.linalg.norm(perp2)

    bond_len = 1.01  # N-H bond length in Angstroms

    h2_dir = -0.33 * n_to_h1 + 0.94 * perp1
    h3_dir = -0.33 * n_to_h1 - 0.47 * perp1 + 0.81 * perp2

    h2_pos = n_pos + bond_len * h2_dir / np.linalg.norm(h2_dir)
    h3_pos = n_pos + bond_len * h3_dir / np.linalg.norm(h3_dir)

    return h2_pos, h3_pos


def clean_polypeptide_pdb(input_path: Path) -> Path:
    """
    Clean a polypeptide PDB file for OpenFF compatibility.

    Args:
        input_path: Path to input PDB file

    Returns:
        Path to output cleaned PDB file
    """
    output_path = input_path.parent / f"{input_path.stem}{CLEANED_PDB_SUFFIX}.pdb"

    residue_atoms = {}  # (chain, resnum, resname) -> list of atom lines
    other_lines = []
    nh2_n_coords = None  # Store NH2 N position to create OXT
    nterm_info = {}  # Store N-terminal atom coords for adding H2, H3

    with open(input_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21:22].strip() or 'A'
                res_num = int(line[22:26])

                # Skip lanthanide/metal ions entirely
                if res_name.upper() in SKIP_RESIDUES or atom_name.upper() in SKIP_RESIDUES:
                    continue

                # Capture NH2 N position for OXT placement, then skip NH2 entirely
                if res_name == 'NH2':
                    if atom_name == 'N':
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        nh2_n_coords = (x, y, z)
                    continue

                # Fix atom names for ASN residues
                if res_name == 'ASN':
                    if atom_name == 'N0':
                        line = line[:12] + ' ND2' + line[16:]
                        atom_name = 'ND2'
                    elif atom_name == 'HN21':
                        line = line[:12] + 'HD21' + line[16:]
                        atom_name = 'HD21'
                    elif atom_name == 'HN22':
                        line = line[:12] + 'HD22' + line[16:]
                        atom_name = 'HD22'

                # Fix C-terminal ALA carboxylate atom names (OD1/OD2 -> O/OXT)
                if res_name == 'ALA':
                    if atom_name == 'OD1':
                        line = line[:12] + '   O' + line[16:]
                        atom_name = 'O'
                    elif atom_name == 'OD2':
                        line = line[:12] + ' OXT' + line[16:]
                        atom_name = 'OXT'

                # Rename HN to H1 for N-terminus compatibility
                if atom_name == 'HN':
                    line = line[:12] + '  H1' + line[16:]
                    atom_name = 'H1'

                key = (chain, res_num, res_name)
                if key not in residue_atoms:
                    residue_atoms[key] = []
                residue_atoms[key].append(line)

                # Store terminal atom info for adding missing hydrogens/OXT
                if atom_name in ('N', 'H1', 'H2', 'H3', 'CA', 'OXT'):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    if key not in nterm_info:
                        nterm_info[key] = {}
                    nterm_info[key][atom_name] = np.array([x, y, z])

            elif line.startswith('CRYST1'):
                other_lines.append(line)

    # Sort residues by chain and residue number
    sorted_keys = sorted(residue_atoms.keys(), key=lambda x: (x[0], x[1]))

    if not sorted_keys:
        raise ValueError(f"No valid residues found in {input_path}")

    n_term_key = sorted_keys[0]
    c_term_key = sorted_keys[-1]

    with open(output_path, 'w') as f:
        # Write CRYST1 if present
        for line in other_lines:
            f.write(line)

        atom_serial = 1
        for key in sorted_keys:
            for line in residue_atoms[key]:
                new_line = f'{line[:6]}{atom_serial:5d}{line[11:]}'
                f.write(new_line)
                atom_serial += 1

            # Add H2, H3 for N-terminal NH3+ (only if not already present)
            if key == n_term_key and key in nterm_info:
                info = nterm_info[key]
                if 'N' in info and 'H1' in info and 'CA' in info:
                    if 'H2' not in info or 'H3' not in info:
                        h2_pos, h3_pos = generate_nterm_hydrogens(info['N'], info['H1'], info['CA'])
                        chain, res_num, res_name = key
                        if 'H2' not in info:
                            h2_line = f'ATOM  {atom_serial:5d}  H2  {res_name:3s} {chain:1s}{res_num:4d}    {h2_pos[0]:8.3f}{h2_pos[1]:8.3f}{h2_pos[2]:8.3f}  1.00  0.00           H\n'
                            f.write(h2_line)
                            atom_serial += 1
                        if 'H3' not in info:
                            h3_line = f'ATOM  {atom_serial:5d}  H3  {res_name:3s} {chain:1s}{res_num:4d}    {h3_pos[0]:8.3f}{h3_pos[1]:8.3f}{h3_pos[2]:8.3f}  1.00  0.00           H\n'
                            f.write(h3_line)
                            atom_serial += 1

            # Add OXT atom for C-terminal carboxylate (only if not already present)
            if key == c_term_key:
                cterm_info = nterm_info.get(key, {})
                if 'OXT' not in cterm_info and nh2_n_coords is not None:
                    chain, res_num, res_name = key
                    oxt_line = f'ATOM  {atom_serial:5d}  OXT {res_name:3s} {chain:1s}{res_num:4d}    {nh2_n_coords[0]:8.3f}{nh2_n_coords[1]:8.3f}{nh2_n_coords[2]:8.3f}  1.00  0.00           O\n'
                    f.write(oxt_line)
                    atom_serial += 1

        f.write('TER\n')
        f.write('END\n')

    print(f"Cleaned: {input_path.name} -> {output_path.name}")
    return output_path


def main():
    parser = argparse.ArgumentParser(
        description='Clean polypeptide PDB files for OpenFF compatibility'
    )
    parser.add_argument(
        'input_files',
        nargs='*',
        help='Input PDB file(s) to clean'
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='Process all .pdb files in current directory (excluding already cleaned files)'
    )

    args = parser.parse_args()

    if args.all:
        current_dir = Path(__file__).parent
        pdb_files = [
            p for p in current_dir.glob('*.pdb')
            if CLEANED_PDB_SUFFIX not in p.stem
        ]
        if not pdb_files:
            print("No PDB files found to process")
            return
        for pdb_file in pdb_files:
            clean_polypeptide_pdb(pdb_file)
    elif args.input_files:
        for input_file in args.input_files:
            input_path = Path(input_file)
            if not input_path.exists():
                print(f"Error: {input_file} not found")
                continue
            clean_polypeptide_pdb(input_path)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
