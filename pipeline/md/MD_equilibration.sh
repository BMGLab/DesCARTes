#!/bin/bash
#
#SBATCH -p a100q
#SBATCH -t 0-01:00:00
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --gres=gpu:1
#SBATCH -A tvybyi
#SBATCH -J dengeleme
#SBATCH --output=dengeleme-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=91240000028@ogrenci.ege.edu.tr

echo "Islemci Sayisi: $SLURM_CPUS_PER_TASK"
echo "Kullanilan GPU: $CUDA_VISIBLE_DEVICES"

# Modulu yukle
module load GROMACS/gromacs-2024.3-openmpi-5.0.5-gcc-11.5.0-cuda-nvidia-hpc-sdk-24.7-cuda-11.8-a100q

# GROMACS'in tam yolu (Hatasiz calismasi icin)
GMX_BIN="/ari/progs/GROMACS/gromacs-2024.3-openmpi-5.0.5-gcc-11.5.0-cuda-nvidia-hpc-sdk-24.7-cuda-11.8-a100q/bin"

echo "Dengeleme Basliyor..."

# --- Dengeleme Adimlari ---

# Adım 6.0: Minimization
# DUZELTME: -r step5_input.gro eklendi
echo "Adim 6.0: Minimization..."
$GMX_BIN/gmx_mpi grompp -f step6.0_minimization.mdp -c step5_input.gro -r step5_input.gro -p topol.top -n index.ndx -o step6.0_minimization.tpr -maxwarn 2
$GMX_BIN/gmx_mpi mdrun -deffnm step6.0_minimization

# Adım 6.1: Equilibration (NVT)
echo "Adim 6.1: NVT..."
$GMX_BIN/gmx_mpi grompp -f step6.1_equilibration.mdp -c step6.0_minimization.gro -r step6.0_minimization.gro -p topol.top -n index.ndx -o step6.1_equilibration.tpr -maxwarn 2
$GMX_BIN/gmx_mpi mdrun -deffnm step6.1_equilibration

# Adım 6.2: NPT 1
echo "Adim 6.2: NPT 1..."
$GMX_BIN/gmx_mpi grompp -f step6.2_equilibration.mdp -c step6.1_equilibration.gro -r step6.1_equilibration.gro -p topol.top -n index.ndx -o step6.2_equilibration.tpr -maxwarn 2
$GMX_BIN/gmx_mpi mdrun -deffnm step6.2_equilibration

# Adım 6.3: NPT 2
echo "Adim 6.3: NPT 2..."
$GMX_BIN/gmx_mpi grompp -f step6.3_equilibration.mdp -c step6.2_equilibration.gro -r step6.2_equilibration.gro -p topol.top -n index.ndx -o step6.3_equilibration.tpr -maxwarn 2
$GMX_BIN/gmx_mpi mdrun -deffnm step6.3_equilibration

# Adım 6.4: NPT 3
echo "Adim 6.4: NPT 3..."
$GMX_BIN/gmx_mpi grompp -f step6.4_equilibration.mdp -c step6.3_equilibration.gro -r step6.3_equilibration.gro -p topol.top -n index.ndx -o step6.4_equilibration.tpr -maxwarn 2
$GMX_BIN/gmx_mpi mdrun -deffnm step6.4_equilibration

# Adım 6.5: NPT 4
echo "Adim 6.5: NPT 4..."
$GMX_BIN/gmx_mpi grompp -f step6.5_equilibration.mdp -c step6.4_equilibration.gro -r step6.4_equilibration.gro -p topol.top -n index.ndx -o step6.5_equilibration.tpr -maxwarn 2
$GMX_BIN/gmx_mpi mdrun -deffnm step6.5_equilibration

# Adım 6.6: NPT 5 (Son Dengeleme)
echo "Adim 6.6: NPT 5..."
$GMX_BIN/gmx_mpi grompp -f step6.6_equilibration.mdp -c step6.5_equilibration.gro -r step6.5_equilibration.gro -p topol.top -n index.ndx -o step6.6_equilibration.tpr -maxwarn 2
$GMX_BIN/gmx_mpi mdrun -deffnm step6.6_equilibration

echo "Dengeleme Basariyla Bitti."
