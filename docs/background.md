# Motivation 

When we perform a single-cell RNA-seq experiment, we often want to compare individual cells’ transcriptomes from one condition (e.g. mutant or disease) against controls (e.g. wild-type or healthy). In principle, doing so is straightforward: first, cluster the cells, then classify the clusters according to type based on the expression of discriminative marker genes. Next, count the number of cells of each type in each sample and compare these counts across sample groups to assess which cell types are more abundant in one group than another. However, this immediately raises new questions: if a particular type of cells disappear in diseased samples, where do they go? Do they transform into one or more pathological cell states? 


# Background

![overview](assets/how_it_works.png)

Hooke is implemented using the [PLNmodels package](https://pln-team.github.io/PLNmodels/index.html). PLN models are a multivariate mixed generalized linear model with a Poisson distribution, allowing them to overcome the computational challenges posed by count data. They provide a convenient framework to perform multivariate statistical regression to describe how environmental effects or perturbations alter the relative abundances of each species.


<!--# Simply computing raw correlations between their abundances can lead to extreme biases and false interactions between species.
# PLN Networks are an extension of the PLN modeling approach that describes how all pairs of species co-vary as a parsimonious network of partial correlations between species that directly interact.
# 
# 
# PLN network models will accurately quantify shifts in the distribution of cells over molecular states following genetic perturbations.
# 
# Hooke aims to directly quantify and clearly report how perturbations redirect cells to new molecular states.
# 
# Modeling the correlation structure between the abundances of different cell states could capture how cells transition between them. However, there are two challenges to doing so, both of which PLN networks will overcome. A first challenge is that cell abundances will co-vary simply because of sampling. PLN networks will distinguish “interesting” flows of cells from one state to the next that occur over time and or in response to perturbations from mundane correlations that occur within replicate samples. A second challenge is that cell counts for some states may be small and sparse, which could imply spurious changes and false correlations across samples.  For example, if we sampled 1,000 cells from a wild-type embryo and a mutant one that lacks an abundant cell type like muscle, rare cell types would appear abundant in the mutant, even though their proportions have only increased because another cell type was lost. PLN networks will address the second challenge because they are explicitly designed to model proportions of multiple species from count data and normalize for technical differences in depth of profiling across samples.-->
