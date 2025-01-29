
#!/bin/sh
python model_spiked_EZ55.py
python model_spiked_detoxers_batch.py 57 
python model_spiked_detoxers_batch.py 58
python model_spiked_detoxers_batch.py 59
python model_spiked_detoxers_batch.py 60
python model_abiotic_batch_1.py # (Figure 3 & Figure S2) 
python model_spiked_abiotic_batch_1.py # (figure 4) 
python model_pro_batch1.py # (figure 5)	
python model_spiked_pro_batch1.py # (figure 6) 
python model_syn_batch.py 53 # (figure 7 & figure S3)
python model_spiked_syn_batch.py 53 # (figure 8) 
python phi_vs_kdam.py # (figure 9) 
python viz_spikes.py # (figure S1)
python model_spiked_pro_batch2.py # (Figure S4)
python model_syn_batch.py 52 # (figure S5)
python model_spiked_syn_batch.py 52 # (figure S6)
python model_syn_batch.py 28 # (figure S7)
python model_spiked_syn_batch.py 28 # (figure S8)

