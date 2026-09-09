#!/bin/bash
#
#SBATCH -p a100q              # A100 GPU node
#SBATCH -t 3-00:00:00         
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64    # 64 CPU
#SBATCH --gres=gpu:1          # A100 GPU
#SBATCH -A project_name
#SBATCH -J rfab_CLDN4
#SBATCH --output=rfab_CLDN4-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mail

set -euo pipefail

############################################
# LOG + TIME
############################################
JOB_START=$(date +%s)
echo "=========================================================="
echo "JOB START : $(date)"
echo "NODE      : $(hostname)"
echo "JOB ID    : $SLURM_JOB_ID"
echo "=========================================================="

############################################
# PATHS
############################################
HOME_DIR="/ari/users/user/RFantibody"
SCRATCH_DIR="/YEREL/rf_${SLURM_JOB_ID}"

TARGET_PDB_SRC="${HOME_DIR}/input_files/9cmiA.pdb"
FRAMEWORK_PDB_SRC="${HOME_DIR}/scripts/examples/example_inputs/hu-4D5-8_Fv.pdb"

############################################
# SCRATCH PREP (KRİTİK ADIM)
############################################
echo "[1/6] SCRATCH hazırlanıyor..."
mkdir -p $SCRATCH_DIR

echo ">>> Dosyalar kopyalanıyor..."
# Temel dosyalar
cp $HOME_DIR/RFantibody_base.sif $SCRATCH_DIR/
cp $HOME_DIR/pyproject.toml $SCRATCH_DIR/
cp $HOME_DIR/poetry.lock $SCRATCH_DIR/

# Klasörler (ÖNCEDEN DERLENMİŞ HALİYLE)
# venv ve include klasörlerini kopyaladığımız için tekrar setup.sh çalıştırmaya gerek kalmayacak.
echo ">>> venv ve kütüphaneler taşınıyor (Bu işlem 1-2 dk sürebilir)..."
cp -r $HOME_DIR/scripts $SCRATCH_DIR/
cp -r $HOME_DIR/src $SCRATCH_DIR/
cp -r $HOME_DIR/bin $SCRATCH_DIR/         # Quiver tools (qvscorefile vb.)
cp -r $HOME_DIR/weights $SCRATCH_DIR/
cp -r $HOME_DIR/include $SCRATCH_DIR/     # USalign ve DGL burada
cp -r $HOME_DIR/venv $SCRATCH_DIR/        # Hazır Python ortamı

# Çalıştırma izinlerini garantiye al
chmod +x $SCRATCH_DIR/bin/*
# USalign binary'sinin iznini kontrol et (varsa)
if [ -f "$SCRATCH_DIR/include/USalign/USalign" ]; then
    chmod +x $SCRATCH_DIR/include/USalign/USalign
fi

echo ">>> Input dosyaları hazırlanıyor..."
mkdir -p $SCRATCH_DIR/inputs
cp $TARGET_PDB_SRC $SCRATCH_DIR/inputs/target.pdb
cp $FRAMEWORK_PDB_SRC $SCRATCH_DIR/inputs/framework.pdb

echo ">>> Output dizinleri oluşturuluyor..."
mkdir -p $SCRATCH_DIR/outputs/CLDN4

cd $SCRATCH_DIR

############################################
# PYTHON ENV CHECK
############################################
echo ""
echo "[1.5/6] Environment kontrolü..."
# setup.sh ÇALIŞTIRMIYORUZ. Kopyalanan venv'i kullanacağız.
# Poetry'nin venv'i gördüğünden emin olalım.

apptainer exec --nv --no-home \
  -B $SCRATCH_DIR:/opt/RFantibody \
  RFantibody_base.sif \
  bash -c "export VIRTUAL_ENV=/opt/RFantibody/venv; export PATH=\$VIRTUAL_ENV/bin:\$PATH; echo 'Python path:'; which python; echo 'Poetry check:'; poetry env info"

echo ">>> ✓ Hazır ortam doğrulandı."

############################################
# PARAMETERS
############################################
HOTSPOTS="[A45,A46,A146,A147,A47,A149,A150,A151,A152,A153,A158]"
LOOPS="[L1:8-13,L2:7,L3:9-11,H1:7,H2:6,H3:5-13]"
NUM_DESIGNS=1000
SEQS_PER_STRUCT=3

############################################
# STEP 1: RFdiffusion → Quiver
############################################
echo ""
echo "[2/6] RFdiffusion başlıyor..."

apptainer run --nv --no-home \
  -B $SCRATCH_DIR:/opt/RFantibody \
  RFantibody_base.sif \
  poetry run python /opt/RFantibody/scripts/rfdiffusion_inference.py \
    --config-name antibody \
    antibody.target_pdb=/opt/RFantibody/inputs/target.pdb \
    antibody.framework_pdb=/opt/RFantibody/inputs/framework.pdb \
    inference.ckpt_override_path=/opt/RFantibody/weights/RFdiffusion_Ab.pt \
    "ppi.hotspot_res=$HOTSPOTS" \
    "antibody.design_loops=$LOOPS" \
    inference.num_designs=$NUM_DESIGNS \
    inference.quiver=/opt/RFantibody/outputs/CLDN4/backbones.qv

############################################
# STEP 2: ProteinMPNN (Quiver → Quiver)
############################################
echo ""
echo "[3/6] ProteinMPNN başlıyor..."

apptainer run --nv --no-home \
  -B $SCRATCH_DIR:/opt/RFantibody \
  RFantibody_base.sif \
  poetry run python /opt/RFantibody/scripts/proteinmpnn_interface_design.py \
    -quiver /opt/RFantibody/outputs/CLDN4/backbones.qv \
    -outquiver /opt/RFantibody/outputs/CLDN4/sequences.qv \
    -seqs_per_struct $SEQS_PER_STRUCT

############################################
# STEP 3: RF2 (Quiver → Quiver)
############################################
echo ""
echo "[4/6] RF2 başlıyor..."

apptainer run --nv --no-home \
  -B $SCRATCH_DIR:/opt/RFantibody \
  RFantibody_base.sif \
  poetry run python /opt/RFantibody/scripts/rf2_predict.py \
    input.quiver=/opt/RFantibody/outputs/CLDN4/sequences.qv \
    output.quiver=/opt/RFantibody/outputs/CLDN4/predictions.qv

############################################
# STEP 4: FILTERING (Quiver → Top PDFs)
############################################
echo ""
echo "[5/6] Top designs extract ediliyor..."

# Scorefile oluştur
apptainer run --nv --no-home \
  -B $SCRATCH_DIR:/opt/RFantibody \
  RFantibody_base.sif \
  /opt/RFantibody/bin/qvscorefile /opt/RFantibody/outputs/CLDN4/predictions.qv > $SCRATCH_DIR/outputs/CLDN4/scores.txt

# Filter (pAE < 10, RMSD < 2.0) ve top 100 seç
cd $SCRATCH_DIR/outputs/CLDN4

awk 'NR>1 && $2!="" && $3!="" && $2<10 && $3<2.0 {print $1, $2, $3}' scores.txt | \
  sort -k2,2g -k3,3g | \
  head -n 100 > top100_list.txt

NUM_PASSED=$(wc -l < top100_list.txt)
echo ">>> $NUM_PASSED designs filtreyi geçti"

# Extract top 100 PDFs
mkdir -p top100_pdbs

if [ -s top100_list.txt ]; then
  apptainer run --nv --no-home \
    -B $SCRATCH_DIR:/opt/RFantibody \
    RFantibody_base.sif \
    bash -c "cut -d' ' -f1 /opt/RFantibody/outputs/CLDN4/top100_list.txt | /opt/RFantibody/bin/qvextractspecific /opt/RFantibody/outputs/CLDN4/predictions.qv"
   
  mv *.pdb top100_pdbs/ 2>/dev/null || echo ">>> Uyarı: PDB dosyaları bulunamadı"
else
  echo ">>> UYARI: Hiçbir design filtreyi geçemedi!"
fi

############################################
# STEP 5: PACKAGING (DISK-FRIENDLY)
############################################
echo ""
echo "[6/6] Sonuçlar paketleniyor..."

# SADECE top 100 PDFs + scores
tar -czf top100_CLDN4_${SLURM_JOB_ID}.tar.gz \
  top100_pdbs/ \
  scores.txt \
  top100_list.txt

# Tar dosyasını HOME'a kopyala
cp top100_CLDN4_${SLURM_JOB_ID}.tar.gz $HOME_DIR/

############################################
# CLEANUP
############################################
echo ">>> Temizlik yapılıyor..."
rm -f backbones.qv sequences.qv predictions.qv
cd /
rm -rf $SCRATCH_DIR

echo "=========================================================="
echo "JOB END   : $(date)"
echo "DURATION  : $((($(date +%s) - JOB_START) / 60)) minutes"
echo "=========================================================="
