CODE_DIR='/data/proj/GCB_MB/spatial_cut-tag/Yale_Apr_2022_nano_Tn5/code/'
RAW_DIR='/data/proj/GCB_MB/spatial_cut-tag/Yale_Apr_2022_nano_Tn5/raw_data/'

python3 $CODE_DIR/run_analysis.py --fastq $RAW_DIR/KI_1_1/ \
                                  --modalities H3K27ac H3K27me3 ATAC \
                                  --barcodes TATAGCCT ATAGAGGC AGATGTGT \
                                  --threads 8 \
                                  --condaenv bcdCT \
                                  --snakeargs='--profile htcondor --rerun-incomplete' &
# sleep 5
python3 $CODE_DIR/run_analysis.py --fastq $RAW_DIR/KI_1_2/ \
                                  --modalities H3K27ac H3K27me3 ATAC\
                                  --barcodes TATAGCCT ATAGAGGC AGATGTGT \
                                  --threads 8 \
                                  --condaenv bcdCT \
                                  --snakeargs='--profile htcondor --rerun-incomplete' &
# sleep 5
python3 $CODE_DIR/run_analysis.py --fastq $RAW_DIR/KI_1_3/ \
                                  --modalities H3K27ac H3K27me3 ATAC\
                                  --barcodes TATAGCCT ATAGAGGC AGATGTGT \
                                  --threads 8 \
                                  --condaenv bcdCT \
                                  --snakeargs='--profile htcondor --rerun-incomplete' &
# sleep 5
python3 $CODE_DIR/run_analysis.py --fastq $RAW_DIR/KI_1_4/ \
                                  --modalities H3K27ac H3K27me3 ATAC\
                                  --barcodes TATAGCCT ATAGAGGC AGATGTGT \
                                  --threads 8 \
                                  --condaenv bcdCT \
                                  --snakeargs='--profile htcondor --rerun-incomplete' &
# sleep 5
python3 $CODE_DIR/run_analysis.py --fastq $RAW_DIR/KI_1_5/ \
                                  --modalities H3K27ac H3K27me3 ATAC\
                                  --barcodes TATAGCCT ATAGAGGC AGATGTGT \
                                  --threads 8 \
                                  --condaenv bcdCT \
                                  --snakeargs='--profile htcondor --rerun-incomplete' &


python3 $CODE_DIR/run_analysis.py --fastq $RAW_DIR/KI_1_5/ \
                                  --modalities H3K27ac H3K27me3 ATAC\
                                  --barcodes TATAGCCT ATAGAGGC AGATGTGT \
                                  --threads 8 \
                                  --condaenv bcdCT \
                                  --snakeargs='--dag' | dot -Tpdf > workflow.pdf


python3 $CODE_DIR/run_analysis.py --fastq $RAW_DIR/KI_5/ \
                                  --modalities H3K27ac H3K27me3\
                                  --barcodes TATAGCCT ATAGAGGC \
                                  --threads 8 \
                                  --condaenv bcdCT \
                                  --snakeargs='--profile htcondor --rerun-incomplete' &


python3 $CODE_DIR/run_analysis.py --fastq $RAW_DIR/KI_6/ \
                                  --modalities H3K27ac H3K27me3\
                                  --barcodes TATAGCCT ATAGAGGC \
                                  --threads 8 \
                                  --condaenv bcdCT \
                                  --snakeargs='--profile htcondor --rerun-incomplete' &

python3 $CODE_DIR/run_analysis.py --fastq $RAW_DIR/KI_7/ \
                                  --modalities H3K27ac H3K27me3 ATAC\
                                  --barcodes TATAGCCT ATAGAGGC CCTATCCT \
                                  --threads 8 \
                                  --condaenv bcdCT \
                                  --snakeargs='--profile htcondor --rerun-incomplete' &

python3 $CODE_DIR/run_analysis.py --fastq $RAW_DIR/KI_8/ \
                                  --modalities H3K27ac H3K27me3 ATAC\
                                  --barcodes TATAGCCT ATAGAGGC CCTATCCT \
                                  --threads 8 \
                                  --condaenv bcdCT \
                                  --snakeargs='--profile htcondor --rerun-incomplete' &





#rsync -razv --include "*/" --include="*.broadPeak" --exclude="*"  --progress marek@monod.mbb.ki.se:/data/proj/GCB_MB/spatial_cut-tag/Yale_Apr_2022_nano_Tn5/results2 ./
#rsync -razv --include "*/" --include="*.narrowPeak" --exclude="*"  --progress marek@monod.mbb.ki.se:/data/proj/GCB_MB/spatial_cut-tag/Yale_Apr_2022_nano_Tn5/results2 ./
#rsync -razv --include "*/" --include="*.bed" --exclude="*"  --progress marek@monod.mbb.ki.se:/data/proj/GCB_MB/spatial_cut-tag/Yale_Apr_2022_nano_Tn5/results2 ./
#rsync -razv --include "*/" --include="*.bw" --exclude="*"  --progress marek@monod.mbb.ki.se:/data/proj/GCB_MB/spatial_cut-tag/Yale_Apr_2022_nano_Tn5/results2 ./