# 1. TSV başlığını oluştur (TAB kullanmak için echo -e "... \t ...")
echo -e "Structure\tAffinity_kcal-mol" > /home/ozangocmen/RFantibody/CLDN_family_combinations/binding_affinity/prodigy_affinity.tsv

# 2. prodigy'nin ham çıktısını (zaten TAB ile ayrılmış) doğrudan dosyaya ekle
prodigy -q -np 48 /home/ozangocmen/RFantibody/CLDN_family_combinations/cif_outputs --selection A,B C >> /home/ozangocmen/RFantibody/CLDN_family_combinations/binding_affinity/prodigy_affinity.tsv