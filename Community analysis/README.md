# Community analysis

Determine communities in the systematic network PhI<sub>MAIN</sub>. For community detection the Edge Betweenness algorithm described in *Girvan and Newman, PNAS, 2002* and implemented in R package igraph was used. 

The scripts in this folder 
- identify communities using edge betweenness or other algorithms in igraph package

- calculate degree and clustering coefficient distribution

- calculate shortest paths between and within hormone signaling pathways

- calculate hormone enrichment of communities with fisher's exact test

How to run:
The supplementary tables 1 and 2 must be put in the folder and renamed to "Extended Table 1 - Search space.xlsx" and "Extended Table 2 - InteractionData.xlsx", respectively. 
Then run "A_community_analysis.r" to do community analysis and several other analyses on the communities and hormone annotations.


