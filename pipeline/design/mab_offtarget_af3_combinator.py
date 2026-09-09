#!/usr/bin/env python3
"""
Generate AlphaFold3 input JSON files for mAb-antigen off-target analysis.

This script:
1. Parses mAb FASTA file (heavy and light chains)
2. Parses off-target antigen FASTA files
3. Generates AF3 JSON files for each mAb-off_target combination
4. Organizes outputs in structured directories
"""

import json
import os
from pathlib import Path
from typing import Dict, List, Tuple


class FASTAParser:
    """Parser for FASTA format files."""
    
    @staticmethod
    def parse_fasta(filepath: str) -> Dict[str, str]:
        """
        Parse a FASTA file and return a dictionary of sequences.
        
        Args:
            filepath: Path to the FASTA file
            
        Returns:
            Dictionary mapping sequence IDs to sequences
        """
        sequences = {}
        current_id = None
        current_seq = []
        
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                if line.startswith('>'):
                    # Save previous sequence
                    if current_id is not None:
                        sequences[current_id] = ''.join(current_seq)
                    
                    # Start new sequence
                    current_id = line[1:].strip()
                    current_seq = []
                else:
                    current_seq.append(line)
            
            # Save last sequence
            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)
        
        return sequences


class mAbParser:
    """Parser for monoclonal antibody sequences."""
    
    def __init__(self, fasta_filepath: str):
        """
        Initialize mAb parser with FASTA file.
        
        Args:
            fasta_filepath: Path to mAb FASTA file
        """
        self.sequences = FASTAParser.parse_fasta(fasta_filepath)
        self.mabs = self._parse_mabs()
    
    def _parse_mabs(self) -> Dict[str, Dict[str, str]]:
        """
        Parse mAb sequences into heavy and light chains.
        
        Returns:
            Dictionary mapping mAb names to their heavy and light chains
        """
        mabs = {}
        
        for seq_id, sequence in self.sequences.items():
            # Parse mAb name and chain type
            if '_H' in seq_id:
                mab_name = seq_id.replace('_H', '')
                chain_type = 'heavy'
            elif '_L' in seq_id:
                mab_name = seq_id.replace('_L', '')
                chain_type = 'light'
            else:
                continue
            
            # Initialize mAb entry if needed
            if mab_name not in mabs:
                mabs[mab_name] = {}
            
            mabs[mab_name][chain_type] = sequence
        
        return mabs
    
    def get_mab_names(self) -> List[str]:
        """Get list of all mAb names."""
        return list(self.mabs.keys())
    
    def get_mab_chains(self, mab_name: str) -> Tuple[str, str]:
        """
        Get heavy and light chain sequences for a mAb.
        
        Args:
            mab_name: Name of the mAb
            
        Returns:
            Tuple of (heavy_chain, light_chain) sequences
        """
        mab = self.mabs.get(mab_name, {})
        return mab.get('heavy', ''), mab.get('light', '')


class AF3InputGenerator:
    """Generator for AlphaFold3 input JSON files."""
    
    def __init__(self, model_seed: int = 42):
        """
        Initialize AF3 input generator.
        
        Args:
            model_seed: Random seed for AlphaFold3 modeling
        """
        self.model_seed = model_seed
    
    def create_af3_json(
        self,
        name: str,
        heavy_chain: str,
        light_chain: str,
        antigen_sequence: str
    ) -> dict:
        """
        Create AlphaFold3 input JSON structure.
        
        Args:
            name: Name for this structure prediction
            heavy_chain: Heavy chain sequence
            light_chain: Light chain sequence
            antigen_sequence: Antigen sequence
            
        Returns:
            Dictionary representing AF3 input JSON
        """
        return {
            "dialect": "alphafold3",
            "version": 1,
            "name": name,
            "modelSeeds": [self.model_seed],
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": heavy_chain,
                        "modifications": []
                    }
                },
                {
                    "protein": {
                        "id": "B",
                        "sequence": light_chain,
                        "modifications": []
                    }
                },
                {
                    "protein": {
                        "id": "C",
                        "sequence": antigen_sequence,
                        "modifications": []
                    }
                }
            ]
        }
    
    def save_json(self, data: dict, filepath: str):
        """
        Save JSON data to file.
        
        Args:
            data: Dictionary to save as JSON
            filepath: Output file path
        """
        with open(filepath, 'w') as f:
            json.dump(data, f, indent=2)


class OffTargetAnalysis:
    """Main class for off-target analysis setup."""
    
    # Mapping of on-targets to off-target FASTA files
    ON_TARGET_TO_FILE = {
        'CLDN1': 'CLDN_off_CLDN1.fasta',
        'CLDN6': 'CLDN_off_CLDN6.fasta',
        'CLDN18_2': 'CLDN_off_CLDN182.fasta'
    }
    
    # Mapping of therapeutic mAbs to their on-targets
    MAB_ON_TARGETS = {
        'lixudebart': 'CLDN1',
        'eclutatug': 'CLDN1',
        'ixotatug': 'CLDN6',
        'garetatug': 'CLDN18_2',
        'ciletatug': 'CLDN18_2',
        'omectatug': 'CLDN18_2',
        'osemitamab': 'CLDN18_2',
        'sonesitamab': 'CLDN18_2',
        'tecotabart': 'CLDN18_2',
        'gresonitamab': 'CLDN18_2',
        'zolbetuximab': 'CLDN18_2'
    }
    
    def __init__(
        self,
        mab_fasta: str,
        offtarget_dir: str,
        output_base_dir: str
    ):
        """
        Initialize off-target analysis.
        
        Args:
            mab_fasta: Path to mAb FASTA file
            offtarget_dir: Directory containing off-target FASTA files
            output_base_dir: Base directory for outputs
        """
        self.mab_parser = mAbParser(mab_fasta)
        self.offtarget_dir = Path(offtarget_dir)
        self.output_base_dir = Path(output_base_dir)
        self.af3_generator = AF3InputGenerator()
        
        # Create output directories
        self.af3_inputs_dir = self.output_base_dir / 'af3_inputs'
        self.af3_outputs_dir = self.output_base_dir / 'af3_outputs'
        
    def generate_all_inputs(self):
        """Generate all AF3 input JSON files."""
        print("Starting AF3 input generation...")
        print(f"Output directory: {self.output_base_dir}")
        print("-" * 80)
        
        total_files = 0
        
        for mab_name in self.mab_parser.get_mab_names():
            # Get mAb on-target
            on_target = self.MAB_ON_TARGETS.get(mab_name)
            
            if not on_target:
                print(f"Warning: No on-target defined for {mab_name}, skipping...")
                continue
            
            # Get corresponding off-target file
            offtarget_file = self.ON_TARGET_TO_FILE.get(on_target)
            
            if not offtarget_file:
                print(f"Warning: No off-target file for {on_target}, skipping {mab_name}...")
                continue
            
            offtarget_path = self.offtarget_dir / offtarget_file
            
            if not offtarget_path.exists():
                print(f"Warning: Off-target file not found: {offtarget_path}, skipping...")
                continue
            
            # Parse off-target sequences
            offtarget_sequences = FASTAParser.parse_fasta(str(offtarget_path))
            
            # Get mAb chains
            heavy_chain, light_chain = self.mab_parser.get_mab_chains(mab_name)
            
            if not heavy_chain or not light_chain:
                print(f"Warning: Missing chain sequences for {mab_name}, skipping...")
                continue
            
            # Create mAb-specific directory
            mab_input_dir = self.af3_inputs_dir / f"{mab_name}_offtarget_af3_input"
            mab_output_dir = self.af3_outputs_dir / f"{mab_name}_offtarget_af3_output"
            
            mab_input_dir.mkdir(parents=True, exist_ok=True)
            mab_output_dir.mkdir(parents=True, exist_ok=True)
            
            print(f"\nProcessing {mab_name} (on-target: {on_target})")
            print(f"  Input folder:  {mab_input_dir}")
            print(f"  Output folder: {mab_output_dir}")
            
            # Generate JSON for each off-target
            files_generated = 0
            for antigen_id, antigen_seq in offtarget_sequences.items():
                # Clean antigen ID for filename
                antigen_name = antigen_id.replace('_HUMAN', '').replace('>', '')
                
                # Create AF3 input JSON
                json_name = f"{mab_name}_vs_{antigen_name}"
                json_data = self.af3_generator.create_af3_json(
                    name=json_name,
                    heavy_chain=heavy_chain,
                    light_chain=light_chain,
                    antigen_sequence=antigen_seq
                )
                
                # Save JSON file
                json_filepath = mab_input_dir / f"{json_name}.json"
                self.af3_generator.save_json(json_data, str(json_filepath))
                
                files_generated += 1
                total_files += 1
            
            print(f"  Generated {files_generated} JSON files")
        
        print("\n" + "=" * 80)
        print(f"COMPLETE! Generated {total_files} AF3 input JSON files")
        print(f"Input directory:  {self.af3_inputs_dir}")
        print(f"Output directory: {self.af3_outputs_dir}")
        print("=" * 80)
        
    def create_summary_file(self):
        """Create a summary file of all generated inputs."""
        summary_path = self.output_base_dir / 'analysis_summary.txt'
        
        with open(summary_path, 'w') as f:
            f.write("AlphaFold3 Off-Target Analysis Summary\n")
            f.write("=" * 80 + "\n\n")
            
            f.write("Therapeutic mAbs and their targets:\n")
            f.write("-" * 80 + "\n")
            for mab, target in sorted(self.MAB_ON_TARGETS.items()):
                f.write(f"  {mab:<20} -> {target}\n")
            
            f.write("\n" + "=" * 80 + "\n")
            f.write("\nDirectory Structure:\n")
            f.write("-" * 80 + "\n")
            f.write(f"  af3_inputs/   - Contains mAb-specific input folders with JSON files\n")
            f.write(f"  af3_outputs/  - Empty folders ready for AF3 prediction results\n")
            
            f.write("\n" + "=" * 80 + "\n")
            f.write("\nHow to run AlphaFold3:\n")
            f.write("-" * 80 + "\n")
            f.write("For each mAb, run AF3 predictions:\n\n")
            
            for mab_name in sorted(self.MAB_ON_TARGETS.keys()):
                mab_input_dir = self.af3_inputs_dir / f"{mab_name}_offtarget_af3_input"
                mab_output_dir = self.af3_outputs_dir / f"{mab_name}_offtarget_af3_output"
                
                if mab_input_dir.exists():
                    f.write(f"{mab_name}:\n")
                    f.write(f"  Input:  {mab_input_dir}\n")
                    f.write(f"  Output: {mab_output_dir}\n")
                    f.write(f"  Command example:\n")
                    f.write(f"    alphafold3 predict --input-dir {mab_input_dir} --output-dir {mab_output_dir}\n\n")
        
        print(f"\nSummary file created: {summary_path}")


def main():
    """Main execution function."""
    # File paths - adjust these to your setup
    MAB_FASTA = '/home/ozangocmen/RFantibody/mab_OTR/mab_fasta_files/mab_list.fasta'
    OFFTARGET_DIR = '/home/ozangocmen/RFantibody/mab_OTR/fasta_files_offtarget/off_targets_by_mabs'
    OUTPUT_BASE_DIR = '/home/ozangocmen/RFantibody/mab_OTR/mab_offtarget_folders'
    
    # Create analysis instance
    analysis = OffTargetAnalysis(
        mab_fasta=MAB_FASTA,
        offtarget_dir=OFFTARGET_DIR,
        output_base_dir=OUTPUT_BASE_DIR
    )
    
    # Generate all inputs
    analysis.generate_all_inputs()
    
    # Create summary
    analysis.create_summary_file()
    
    print("\n✓ All done! Your AF3 input files are ready.")


if __name__ == '__main__':
    main()