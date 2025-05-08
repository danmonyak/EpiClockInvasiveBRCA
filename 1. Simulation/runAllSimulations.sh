echo "\n######################################################################################"
echo "########### Running single_cell_metastasis.py, might take up to 20 minutes ###########"
for i in $(seq 1 4);
do
	echo "######################################################################################"
done
echo "\n\n\n"
python single_cell_metastasis.py

echo "\n######################################################################################"
echo "########### Running synchronous_met.py, might take up to 20 minutes ###########"
for i in $(seq 1 4);
do
	echo "######################################################################################"
done
echo "\n\n"
python synchronous_met.py


echo "\n######################################################################################"
echo "########### Running the 3 simulations set up in run_ensemble.py, might take up to an hour minutes ###########"
for i in $(seq 1 4);
do
	echo "######################################################################################"
done
echo "\n\n\n\n"

for i in $(seq 0 2);
do
	echo "######################################################################################"
	echo "Running simulation $i"
	echo "######################################################################################\n\n\n\n"
	python run_ensemble.py $i
done
