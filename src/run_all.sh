
#!/bin/sh
python model_pro_batch1.py
python viz_spikes.py
python model_spiked_abiotic_batch_1.py
python model_abiotic_batch_1.py
python model_spiked_pro_batch1.py
python model_spiked_pro_batch2.py
python model_detoxers_batch.py 57 
python model_detoxers_batch.py 58
python model_detoxers_batch.py 59
python model_detoxers_batch.py 60
python model_syn_batch.py 53
python model_syn_batch.py 52
python model_syn_batch.py 28
python model_spiked_syn_batch.py 53
python model_spiked_syn_batch.py 52
python model_spiked_syn_batch.py 28
python model_spiked_detoxers_batch.py 57 
python model_spiked_detoxers_batch.py 58
python model_spiked_detoxers_batch.py 59
python model_spiked_detoxers_batch.py 60
python phi_vs_kdam.py
#python model_spiked_EZ55.py

