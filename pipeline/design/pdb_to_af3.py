#!/usr/bin/env python3
import os
import glob
import json
from Bio import SeqIO
from Bio.PDB import PDBParser, PPBuilder

# ================= AYARLAR =================
# 1. Antikor PDB dosyalarinin oldugu klasor (H, L ve T zincirleri iceren PDB'ler)
ANTIBODY_PDB_DIR = "/home/ozangocmen/RFantibody/CLDN4/output_CLDN4/selected_pdb"

# 2. Hedef (Target) proteinlerin oldugu FASTA dosyasi
TARGET_FASTA_FILE = "/home/biolab/Projects/DesCARTes_wd/binder/RFantibody/CLDN_family/all_results/CLDN_family.fasta"

# 3. JSON dosyalarinin kaydedilecegi klasor
OUTPUT_JSON_DIR = "/home/ozangocmen/RFantibody/CLDN_family_combinations/json_inputs"

# 4. AF3 Model Seeds
MODEL_SEEDS = [42] 
# ===========================================

def get_hl_sequences_only(pdb_path):
    """
    PDB dosyasindan SADECE 'H' ve 'L' zincirlerinin sekanslarini ceker.
    T zincirini veya baska zincirleri almaz.
    Geriye { 'H': 'sekans...', 'L': 'sekans...' } sozlugu dondurur.
    """
    parser = PDBParser(QUIET=True)
    ppb = PPBuilder()
    
    try:
        structure = parser.get_structure('temp', pdb_path)
    except Exception as e:
        print(f"HATA: PDB okunamadi {pdb_path}: {e}")
        return None

    extracted_chains = {}
    
    # Model icindeki zincirleri kontrol et
    for model in structure:
        for chain in model:
            chain_id = chain.id
            
            # Sadece H ve L zincirlerini istiyoruz
            if chain_id in ['H', 'L']:
                peptides = ppb.build_peptides(chain)
                if peptides:
                    # Polipeptit parcaciklarini birlestir
                    seq_str = "".join(str(pp.get_sequence()) for pp in peptides)
                    extracted_chains[chain_id] = seq_str
            
            # T zinciri veya baska zincirler (A, B, C vb.) gelirse buraya girmez, atlanir.

    # Kontrol: H veya L eksik mi?
    if 'H' not in extracted_chains or 'L' not in extracted_chains:
        # Eger H/L ID'leri yoksa (bazen A/B olarak gelir), ilk iki zinciri H ve L varsaymak icin 
        # asagidaki blogu acabilirsiniz. Ancak simdilik kati kural uyguluyoruz.
        return None
    
    return extracted_chains

def clean_filename(name):
    """Dosya isimleri icin guvenli olmayan karakterleri temizler."""
    return name.replace("|", "_").replace(" ", "_").replace("/", "-")

def main():
    os.makedirs(OUTPUT_JSON_DIR, exist_ok=True)

    # 1. Targetlari Yukle
    print(f"Target dosyasi okunuyor: {TARGET_FASTA_FILE}")
    targets = list(SeqIO.parse(TARGET_FASTA_FILE, "fasta"))
    print(f"Toplam {len(targets)} hedef protein bulundu.\n")

    # 2. PDB Dosyalarini Bul
    pdb_files = sorted(glob.glob(os.path.join(ANTIBODY_PDB_DIR, "*.pdb")))
    if not pdb_files:
        print(f"UYARI: {ANTIBODY_PDB_DIR} dizininde PDB bulunamadi!")
        return

    print(f"Toplam {len(pdb_files)} antikor PDB dosyasi islenecek.")
    
    count = 0

    # 3. Kombinasyon Dongusu
    for pdb_path in pdb_files:
        pdb_name = os.path.basename(pdb_path).replace(".pdb", "")
        
        # PDB'den SADECE H ve L sekanslarini al (T'yi at)
        hl_seqs = get_hl_sequences_only(pdb_path)
        
        if not hl_seqs:
            print(f"ATLANDI (H veya L zinciri bulunamadi): {pdb_name}")
            continue

        # Zincir sirasi: Once H, Sonra L
        seq_h = hl_seqs['H']
        seq_l = hl_seqs['L']

        for target in targets:
            target_name = target.id
            target_seq = str(target.seq)

            # Dosya ismi olustur
            safe_ab_name = clean_filename(pdb_name)
            safe_target_name = clean_filename(target_name)
            
            # Ornek: memb_CLDN4_..._vs_CLDN4_Human
            job_name = f"{safe_ab_name}_vs_{safe_target_name}"
            
            # --- AF3 JSON Yapisi ---
            af3_data = {
                "dialect": "alphafold3",
                "version": 1,
                "name": job_name,
                "modelSeeds": MODEL_SEEDS,
                "sequences": [
                    # 1. Zincir: Heavy Chain (PDB'den gelen H) -> ID: A
                    {
                        "protein": {
                            "id": "A",
                            "sequence": seq_h,
                            "modifications": []
                        }
                    },
                    # 2. Zincir: Light Chain (PDB'den gelen L) -> ID: B
                    {
                        "protein": {
                            "id": "B",
                            "sequence": seq_l,
                            "modifications": []
                        }
                    },
                    # 3. Zincir: Target Chain (FASTA'dan gelen) -> ID: C
                    {
                        "protein": {
                            "id": "C",
                            "sequence": target_seq,
                            "modifications": []
                        }
                    }
                ]
            }

            # JSON kaydet
            out_name = f"{job_name}.json"
            out_path = os.path.join(OUTPUT_JSON_DIR, out_name)

            try:
                with open(out_path, "w") as f:
                    json.dump(af3_data, f, indent=2)
                count += 1
            except Exception as e:
                print(f"Yazma hatasi: {e}")

    print(f"\n✅ İŞLEM TAMAMLANDI.")
    print(f"📂 Çıktılar: {OUTPUT_JSON_DIR}")
    print(f"📄 Toplam oluşturulan JSON: {count}")

if __name__ == "__main__":
    main()