# RPath
<p align="center">
  <img src="data/img/logo.jpg" alt="logo" width="200"/>
</p>

Source code and data presented in "Causal reasoning over knowledge graphs leveraging multimodal transcriptomic signatures for drug discovery".

## Structure

- `data`: contains the four transcriptomic datasets used in the paper as well as the two KGs. Furthermore, it contains the drug-disease pairs studied in clinical trials which are used as positive pairs for the validation of the RPath algorithm.
- `src`: contains the code to run the algorithm and its validation as well as to preprocess the datasets. It is divided between `data_preprocess` and `notebooks`. Notebooks contains the core of the analysis and validation.
