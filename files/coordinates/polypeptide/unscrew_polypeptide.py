"""
Polypeptide PDB cleaning script for OpenFF compatibility.

Uses PDBFixer (OpenMM) to automatically:
- Add missing heavy atoms
- Add missing hydrogens at specified pH
- Handle terminal residues (NH3+, COO-)

Manual preprocessing handles non-standard elements:
- Fixing VMD-generated atom names (N0→ND2, HN21→HD21, etc.)
- Removing TB/lanthanide ions (handled separately by project.py)
- Removing NH2 capping groups (PDBFixer adds proper C-terminus)

Usage:
    python unscrew_polypeptide.py LBT3-.pdb
    python unscrew_polypeptide.py LBT5-.pdb
    python unscrew_polypeptide.py --all  # Process all .pdb files in directory

Output:
    {input_stem}_cleanedPDB.pdb
"""

import argparse
import tempfile
from pathlib import Path

# PDBFixer for automatic structure repair
from pdbfixer import PDBFixer
from openmm.app import PDBFile


CLEANED_PDB_SUFFIX = "_cleanedPDB"
DEFAULT_PH = 7.0

# Known lanthanide/metal residue names to skip
SKIP_RESIDUES = {'TB', 'ND', 'LA', 'CE', 'PR', 'SM', 'EU', 'GD', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU'}

# Atom name fixes for VMD-generated PDBs
# Format: (residue_name, old_atom_name) -> new_atom_name
ATOM_NAME_FIXES = {
    ('ASN', 'N0'): 'ND2',
    ('ASN', 'HN21'): 'HD21',
    ('ASN', 'HN22'): 'HD22',
    ('ALA', 'OD1'): 'O',
    ('ALA', 'OD2'): 'OXT',
}


def preprocess_pdb(input_path: Path) -> str:
    """
    Preprocess PDB to fix non-standard atom names and remove non-standard residues.

    Groups atoms by residue to handle scattered atoms (e.g., ASN atoms appearing
    after C-terminus in VMD-generated files).

    Returns path to temporary preprocessed PDB file.
    """
    # Group atoms by (chain, resnum, resname) to handle scattered atoms
    residue_atoms = {}  # (chain, resnum, resname) -> list of lines
    header_lines = []

    with open(input_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21:22].strip() or 'A'
                res_num = int(line[22:26])

                # Skip lanthanide/metal ions
                if res_name.upper() in SKIP_RESIDUES or atom_name.upper() in SKIP_RESIDUES:
                    continue

                # Skip NH2 capping group (PDBFixer will add proper C-terminus)
                if res_name == 'NH2':
                    continue

                # Fix non-standard atom names
                fix_key = (res_name, atom_name)
                if fix_key in ATOM_NAME_FIXES:
                    new_name = ATOM_NAME_FIXES[fix_key]
                    # Pad atom name to 4 characters
                    if len(new_name) < 4:
                        new_name = f'{new_name:>4}'
                    line = line[:12] + new_name + line[16:]

                # Rename HN to H for backbone amide (PDBFixer expects H, not HN)
                if atom_name == 'HN':
                    line = line[:12] + ' H  ' + line[16:]

                # Group by residue
                key = (chain, res_num, res_name)
                if key not in residue_atoms:
                    residue_atoms[key] = []
                residue_atoms[key].append(line)

            elif line.startswith('CRYST1'):
                header_lines.append(line)

    # Sort residues by chain and residue number
    sorted_keys = sorted(residue_atoms.keys(), key=lambda x: (x[0], x[1]))

    # Build output with atoms grouped by residue
    lines_out = header_lines.copy()
    atom_serial = 1
    for key in sorted_keys:
        for line in residue_atoms[key]:
            # Renumber atoms sequentially
            new_line = f'{line[:6]}{atom_serial:5d}{line[11:]}'
            lines_out.append(new_line)
            atom_serial += 1

    lines_out.append('TER\n')
    lines_out.append('END\n')

    # Write to temporary file
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
    tmp.writelines(lines_out)
    tmp.close()

    return tmp.name


def clean_polypeptide_pdb(input_path: Path, ph: float = DEFAULT_PH) -> Path:
    """
    Clean a polypeptide PDB file using PDBFixer.

    Args:
        input_path: Path to input PDB file
        ph: pH for protonation states (default: 7.0)

    Returns:
        Path to output cleaned PDB file
    """
    output_path = input_path.parent / f"{input_path.stem}{CLEANED_PDB_SUFFIX}.pdb"

    # Step 1: Preprocess to fix atom names and remove non-standard residues
    preprocessed_path = preprocess_pdb(input_path)

    try:
        # Step 2: Use PDBFixer to repair the structure
        fixer = PDBFixer(filename=preprocessed_path)

        # Find and add missing heavy atoms
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

        # Add hydrogens at physiological pH
        # This automatically handles:
        # - N-terminus: NH3+ (adds H1, H2, H3)
        # - C-terminus: COO- (no extra H on carboxylate)
        # - Standard residue protonation states
        fixer.addMissingHydrogens(pH=ph)

        # Step 3: Write the cleaned PDB
        with open(output_path, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)

        print(f"Cleaned: {input_path.name} -> {output_path.name} (pH={ph})")

    finally:
        # Clean up temporary file
        Path(preprocessed_path).unlink(missing_ok=True)

    return output_path


def main():
    parser = argparse.ArgumentParser(
        description='Clean polypeptide PDB files for OpenFF compatibility using PDBFixer'
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
    parser.add_argument(
        '--ph',
        type=float,
        default=DEFAULT_PH,
        help=f'pH for protonation states (default: {DEFAULT_PH})'
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
            try:
                clean_polypeptide_pdb(pdb_file, ph=args.ph)
            except Exception as e:
                print(f"Error processing {pdb_file.name}: {e}")
    elif args.input_files:
        for input_file in args.input_files:
            input_path = Path(input_file)
            if not input_path.exists():
                print(f"Error: {input_file} not found")
                continue
            try:
                clean_polypeptide_pdb(input_path, ph=args.ph)
            except Exception as e:
                print(f"Error processing {input_file}: {e}")
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
