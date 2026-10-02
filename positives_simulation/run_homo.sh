
export PYTHONPATH=/home/avraham/MaruvkaLab/msmutect_postprocessing
for v in 0.05 0.1 0.15 0.2 0.25; do
    python homo_homo_test.py "$v" &
done
wait