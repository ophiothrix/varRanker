## So we've seen that we can lose non-core histone marks, including the limiting of H3K27ac and H3K9ac to binary marks. Which means we can probably gather enough positive training point for exact mark. So the new plan is:
## Gather all possible tissues from AS DHS with the perfect match to ROADMAP tissue.
## Annotate them with core marks. Where available use only binary H3K27ac and H3K9ac, where absent - use imputed information.
## Train a model on that perfect match, look at the cross-validation rates
## Use an imperfect match as a test set.
