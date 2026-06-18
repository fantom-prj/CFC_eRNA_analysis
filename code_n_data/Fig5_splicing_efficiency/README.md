# Splicing efficiency

* splicing_junction.R -> main code describing all the analyses
* junction_run.py -> associated code for running splice AI
* Results are presented in Fig.5 and ext_Fig.7
* Data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_CFC

# Related Methods
* [Splicing efficiency and spliceAI](#SJ)


# <a name="SJ"></a>Splicing efficiency and spliceAI

The splicing efficiency was calculated for all the observed donor sites and acceptor sites from the aligned reads, number of splice / (number of splice + number of span). The splicing efficiency of all the observed splice junctions (SJs) was obtained by number of splice / (number of splice + number of span), where the “number of span” represents the average spanning event from the corresponding donor site and acceptor site (Fig. 5b). The calculated splicing efficiencies were grouped at gene levels, where the SJs belonging to the reads were linked to transcript models and to gene models, and the mean splicing efficiency normalized by intron count was used. <br>
All identified donor and acceptor splice sites were extended by ±40 nt and converted into FASTA format. These sequences containing both intron and exon were then analysed using spliceAI (v1.3.1) to predict splicing potential based on the surrounding sequence composition. The prediction was done using the approach for scoring custom sequences described on the official GitHub page. Briefly, it involves padding a custom input sequence to 10,000 nucleotides. The padded sequence is converted into a one-hot encoded format suitable for model input. Five pre-trained SpliceAI models are applied to the encoded sequence, and the predictions are averaged to obtain acceptor and donor site probabilities. <br>


