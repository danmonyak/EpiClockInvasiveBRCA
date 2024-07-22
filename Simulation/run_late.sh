i=$1
python -u late_metastasis.py $i > late_met/output_files/out_${i}.txt 2>&1 &
