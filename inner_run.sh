g++ -o marginals marginals.c
./marginals
python ../collate_marginals.py marginals.csv computed_shapley_value.py
g++ -O3 -o run_SEBM run_SEBM.c
python ./run_methods.py
./run_SEBM 10
./run_SEBM 50
./run_SEBM 100
./run_SEBM 500
./run_SEBM 1000
echo "maleki"
python compute_error.py ./data_out_maleki_2250.csv
python compute_error.py ./data_out_maleki_11250.csv
python compute_error.py ./data_out_maleki_22500.csv
python compute_error.py ./data_out_maleki_112500.csv
python compute_error.py ./data_out_maleki_225000.csv
echo "simple"
python compute_error.py ./data_out_simple_2250.csv
python compute_error.py ./data_out_simple_11250.csv
python compute_error.py ./data_out_simple_22500.csv
python compute_error.py ./data_out_simple_112500.csv
python compute_error.py ./data_out_simple_225000.csv
echo "castro"
python compute_error.py ./data_out_castro_2250.csv
python compute_error.py ./data_out_castro_11250.csv
python compute_error.py ./data_out_castro_22500.csv
python compute_error.py ./data_out_castro_112500.csv
python compute_error.py ./data_out_castro_225000.csv
echo "sebm"
python compute_error.py ./data_out_Burgess_2250.csv
python compute_error.py ./data_out_Burgess_11250.csv
python compute_error.py ./data_out_Burgess_22500.csv
python compute_error.py ./data_out_Burgess_112500.csv
python compute_error.py ./data_out_Burgess_225000.csv
