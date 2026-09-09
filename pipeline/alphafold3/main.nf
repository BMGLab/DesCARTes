nextflow.enable.dsl=2

// Varsayılan parametreler (Komut satırından ezilebilir)
params.inputs    = null // Komut satırından verilecek
params.out_base  = "/home/biolab/Projects/DesCARTes_wd/results/OTRtest/af3_outputs_99" // Çıktı ana dizini
params.db_dir    = null
params.model_dir = null
params.c_db_dir    = "/root/public_databases"
params.c_model_dir = "/root/models"

workflow {
    // 1. Inputları topla ve Klasör İsimlerini Analiz Et
    Channel
        .fromPath(params.inputs, checkIfExists: true)
        .map { file -> 
            def id = file.baseName // dosya adı (uzantısız): eclutatug_vs_CLD2
            def parent_dir = file.parent.name // üst klasör: eclutatug_offtarget_af3_input
            
            // "input" kelimesini "output" ile değiştirerek hedef klasör adını oluştur
            def target_subdir = parent_dir.replace('input', 'output') 
            
            // Tuple yapısı: [id, hedef_klasör_adı, dosya_yolu]
            return tuple(id, target_subdir, file)
        }
        .set { input_ch }

    // Stage 1: Data Pipeline (CPU)
    // sub_dir bilgisini taşıyoruz
    data_ch = AF3_DATA_PIPELINE(input_ch)

    // Stage 2: Inference (GPU)
    AF3_INFERENCE(data_ch)
}

process AF3_DATA_PIPELINE {
    tag { id }
    // Kaynak ayarları nextflow.config dosyasından alınır

    input:
        tuple val(id), val(sub_dir), path(json_file)

    output:
        // DEĞİŞİKLİK BURADA:
        // Klasör yapısı veya harf büyüklüğü ne olursa olsun "_data.json" ile biten dosyayı yakala
        tuple val(id), val(sub_dir), path("${id}/**/*_data.json")

    script:
    """
    set -euo pipefail
    mkdir -p ${id}

    python /app/alphafold/run_alphafold.py \
        --norun_inference \
        --json_path ${json_file} \
        --output_dir ${id} \
        --db_dir ${params.c_db_dir} \
        --model_dir ${params.c_model_dir}
    """
}

process AF3_INFERENCE {
    tag { id }
    
    // DİKKAT: Çıktı dizini burada dinamik olarak belirleniyor
    // Örn: .../af3_outputs/eclutatug_offtarget_af3_output/eclutatug_vs_CLD2/
    publishDir "${params.out_base}/${sub_dir}/${id}", mode: 'copy'

    input:
        tuple val(id), val(sub_dir), path(data_json)

    output:
        path("*") // Tüm sonuçları kopyala

    script:
    """
    set -euo pipefail
    export XLA_PYTHON_CLIENT_PREALLOCATE=false

    python /app/alphafold/run_alphafold.py \
        --norun_data_pipeline \
        --json_path ${data_json} \
        --output_dir . \
        --db_dir ${params.c_db_dir} \
        --model_dir ${params.c_model_dir} \
        --num_diffusion_samples 1
    """
}