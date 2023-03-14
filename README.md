# mars-biofoundry
Metabolic engineering towards a biofoundry on Mars



### Number of pathways: 113
### Number of pathways with modules: 113 - 72
### Number of modules without reactions: 2
### Number of modules processed: 169 - 2



# How to get all rules for a given organism

## Focusing on ECs:

1. By EC from genes: parse https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+gn:T04747

2. By iterating over each reaction: parse https://www.kegg.jp/entry/R01015

3. By getting EC from pathway: bioservices or https://www.kegg.jp/dbget-bin/www_bget?tez00450


## Focusing on reaction numbers:

1. Use SPARQL to retrieve IDs for RR1 in MetaNetX: https://rdf.metanetx.org/

2. Convert to EC (see point 2. above)
