## Learning the language of TCR-Epitope-MHC interactions with minimal finetuning (Interface prediction)

Download the `stcrdab_pdb.csv` file and the `contact_map` folder from this [repo](https://github.com/pengxingang/TEIM/tree/main/data). If you want to do sequence-level pretraining, you should also download the `binding_data` folder. Run `train_res.py` to train our model on this task. Altough we did not use sequence-level pretraining in our final experiments, you can run it using `train_seq.py`. 
