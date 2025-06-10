from baselines import ProtT5, ProGen, ESM

model = ProtT5(model_name='Rostlab/prot_t5_xl_bfd', layers='last', devices=[0], batch_size=2)
model = ProtT5(model_name='Rostlab/prot_t5_xl_uniref50', layers='last', devices=[0], batch_size=2)
model = ProGen(model_name='hugohrban/progen2-xlarge', layers='last', devices=[0], batch_size=2)
model = ProGen(model_name='hugohrban/progen2-large', layers='last', devices=[0], batch_size=2)
model = ESM(model_name='facebook/esm2_t36_3B_UR50D', layers='last', devices=[0], batch_size=2)
model = ESM(model_name='facebook/esm2_t33_650M_UR50D', layers='last', devices=[0], batch_size=2)
model = ESM(model_name='facebook/esm2_t30_150M_UR50D', layers='last', devices=[0], batch_size=2)
model = ESM(model_name='facebook/esm1b_t33_650M_UR50S', layers='last', devices=[0], batch_size=2)