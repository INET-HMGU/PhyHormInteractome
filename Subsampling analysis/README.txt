# Pathway contact point subsampling analysis

A subsampling approach was used to compare the number of pathway contact points between the PhI<subscript>MAIN</subscript> network and literature curated proteins interaction networks from IntAct and BioGRID database. Therefore 1,000 times from each network 100 interactions were sampled and the number of pathway contact points were counted. The resulting distributions are compared to indentify different connectivity of phytohormone pathways in the networks. 

How to run: The supplementary table 2 must be put in the folder and renamed to "Extended Table 2 - InteractionData.xlsx". Then first run "G_plot_crosstalk_results.r" to count pathway contact points and do subsampling analysis. Afterwards run "J_monte_carlo_evaluation.r" to plot the results of the subsampling evaluation and calculate P-values. 
