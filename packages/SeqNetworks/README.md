## Seq Networks

Python class for building networks from similar strings

##installation
```
pip install --upgrade git+http://ghe-rss.roche.com/plsRED-Bioinformatics/SeqNetworks.git
```

## Utilizing the SeqNetwork class

1) *Instatiate the class and add dimensions*

```
aln_df = pandas.read_csv("/sc1/groups/pls-redbfx/pipeline_runs/iPETE/production/2019-07-15-Expt47_cell_lines_dedup_params/analysis/SUP_T1/filtered_reads/SUP_T1_on_target.tsv", sep="\t")

graph = SeqNetwork(1, 2, 2)

for read_name, uid, cdr3, uid_qual, read_qual in zip(aln_df["name"], aln_df["uid"],
    	       	    	  	    	      	     aln_df["cdr3"], aln_df["uid_qual"],
						     aln_df["avg_read_qual"]):
	quals = [int(x) for x in uid_qual.split(",")]
	if sum(quals)/len(quals) < 30:
		continue
	if read_qual < 30:
		continue
	if isinstance(cdr3, str) and isinstance(uid, str): 
	    graph.add_dimension("UMI", read_name, uid, 1)
	    graph.add_dimension("CDR3", read_name, cdr3, 1)

```

First, we instantiate the class defining the parameters used for deduping

```
graph = SeqNetwork(processors, bidRatio, maxSteps)
```

In the example above, two dimensions have been added to the graph `UMI` and `CDR3`.
To add a dimension to each read, four arguments are needed

```
graph.add_dimension(dimension_name, read_name, string, string_dist)
```

`dimension_name`: provides a label for the dimension
`read_name`: links all strings across all dimensions
`string`: is the string used for similarity searching (e.g. UMI)
`string_dist`: is the distance value used in the BKtree similarity search for each dimension.


2) *build the graph*
After all dimensions have been added, BKtrees are utilized to identify all similar strings. Here, nodes represent all reads that share identical `strings` between all dimensions. Likewise, edges connect reads which share the same `string_dist` constraints across all dimensions.  

```
graph.build_graph()
```

3) *partition the graph*
once the graph has built, it can be partitioned by comparing the weights between adjecent nodes. 
```
graph.remove_edges()
```

4) *report similarity families, or UMI families*
After the graph has been partitioned, all joined nodes are merged together. There reads are retured with their UMI family labels.

```
families = graph.collapse_graph()
family_labels = graph.retreive_read_labels(families)	

```

5) *summarize graph characteristics*
A high level summary of the graph characteristics are also available
```
graph.summarize_graph()
```

Example of output
```
number of reads = 1019889
number of dimensions = 2
	UMI unique seqs = 9812
	CDR3 unique seqs = 1086
number of nodes = 19930
total centroids = 1021

```












