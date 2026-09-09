nextflow run main.nf \
  --inputs '/home/biolab/Projects/DesCARTes_wd/data/phase2_OTRtestPDBbind/*_json/*.json' \
  --db_dir /mnt/ssd0/alphafold3_db \
  --model_dir /home/biolab/alphafold3/models \
  -with-docker 'alphafold3:latest' \
  -resume