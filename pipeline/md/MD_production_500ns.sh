#SBATCH -p a100q
#SBATCH -t 3-00:00:00           
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --gres=gpu:1
#SBATCH -A project_name
#SBATCH -J prod_500ns
#SBATCH --output=prod_500ns-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mail

echo "=========================================="
echo "500 ns Üretim Simülasyonu Başlatılıyor"
echo "GPU Node'una bağlanıldı: $(date)"
echo "=========================================="

# --- 1. MODÜL VE ORTAM AYARLARI ---
module purge
module load GROMACS/gromacs-2024.3-openmpi-5.0.5-gcc-11.5.0-cuda-nvidia-hpc-sdk-24.7-cuda-11.8-a100q

# GMXRC Source
echo "GROMACS ortamı yükleniyor..."
source /ari/progs/GROMACS/gromacs-2024.3-openmpi-5.0.5-gcc-11.5.0-cuda-nvidia-hpc-sdk-24.7-cuda-11.8-a100q/bin/GMXRC

# --- 2. SÜRE AYARLARI: 500 ns = 250,000,000 adım (2 fs timestep) ---
echo "=========================================="
echo "Simülasyon parametreleri ayarlanıyor..."
echo "Hedef süre: 500 ns (250 milyon adım)"
echo "=========================================="

# step7_production.mdp dosyasını 500 ns için düzenle
sed -i 's/^nsteps.*/nsteps = 250000000/' step7_production.mdp

# Checkpoint ve output frekanslarını optimize et (RMSD, RMSF, Rg, SASA, MM/PBSA için yeterli)
sed -i 's/^nstxout-compressed.*/nstxout-compressed = 500000/' step7_production.mdp  # Her 1 ns (500 frame toplam)
sed -i 's/^nstlog.*/nstlog = 50000/' step7_production.mdp                           # Log her 100 ps
sed -i 's/^nstenergy.*/nstenergy = 50000/' step7_production.mdp                     # Enerji her 100 ps

# --- 3. TPR DOSYASI OLUŞTURMA ---
echo "=========================================="
echo "TPR dosyası oluşturuluyor..."
echo "=========================================="
gmx_mpi grompp -f step7_production.mdp \
               -c step6.6_equilibration.gro \
               -p topol.top \
               -n index.ndx \
               -o step7_production_500ns.tpr \
               -maxwarn 2

# Kontrol: TPR oluşturuldu mu?
if [ $? -ne 0 ]; then
    echo "=========================================="
    echo "HATA: TPR dosyası oluşturulamadı!"
    echo "Lütfen mdp, gro, top ve index dosyalarını kontrol edin."
    echo "=========================================="
    exit 1
fi

echo "TPR dosyası başarıyla oluşturuldu."

# --- 4. SİMÜLASYON BAŞLATMA ---
echo "=========================================="
echo "500 ns Üretim Simülasyonu Başlıyor..."
echo "Başlangıç zamanı: $(date)"
echo "=========================================="

# GPU doğrudan iletişim özelliğini etkinleştir (isteğe bağlı performans iyileştirmesi)
export GMX_ENABLE_DIRECT_GPU_COMM=1

# Simülasyonu çalıştır
gmx_mpi mdrun -s step7_production_500ns.tpr \
              -deffnm step7_production_500ns \
              -v

# --- 5. TAMAMLANMA KONTROLÜ ---
if [ $? -eq 0 ]; then
    echo "=========================================="
    echo "SİMÜLASYON BAŞARIYLA TAMAMLANDI!"
    echo "Bitiş zamanı: $(date)"
    echo "=========================================="

    # Çıktı dosyalarını listele
    echo "Oluşturulan dosyalar:"
    ls -lh step7_production_500ns.*

else
    echo "=========================================="
    echo "UYARI: Simülasyon beklenmedik şekilde sonlandı!"
    echo "Checkpoint dosyalarından devam edilebilir."
    echo "=========================================="
fi

echo "=========================================="
echo "İşlem sona erdi: $(date)"
echo "=========================================="
