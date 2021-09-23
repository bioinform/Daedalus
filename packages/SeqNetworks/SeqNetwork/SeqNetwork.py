import collections
import Levenshtein
from SeqNetwork.bktree import BKtree
import math
import pandas
import multiprocessing
import time

"""
A class for building similarity networks between one or more strings within a sequencing read.
Based on sequence similarity, reads are clustered together into 'UMI' families. While the main purpose of this class is for deduplication of Read data, the core functions can also be used to define sequence clusters in other contexts.
"""

class SeqNetwork:
    """
    Construct a similarity network for a set of strings.       
    """
    def __init__(self):
        self.graph = collections.OrderedDict()
        self.nodes = collections.OrderedDict()
        self.index = {}
        self.node_sizes = {}
        self.index_reads = {}
        self.dimensions = collections.OrderedDict()
        self.dimension_distances = collections.OrderedDict()
        self.edges = {}
        self.BKtrees = {}        
        self.total_dimensions = 0        
        self.total_centroids = 0

    def summarize_graph(self):
        """
        Summarize information about the graph
        Returns 
        -------
        Report : list
        A list of strings, reporting information about the graph.
        """
        report = []
        dimensions = list(self.dimensions.keys())
        report.append("number of reads = {}".format(len(self.dimensions[dimensions[0]])))
        report.append("number of dimensions = {}".format(self.total_dimensions))
        for dimName in dimensions:
            dimTotal = len(set(self.dimensions[dimName].values()))
            report.append("{} unique seqs = {}".format(dimName, dimTotal))
        report.append("number of nodes = {}".format(len(self.graph)))
        report.append("total centroids = {}".format(self.total_centroids))        
        return report
        
    def add_dimension(self, dimName, readName, seqString, stringDist):
        """
        Add dimensions (sequences) to the graph
        
        Parameters
        ----------
        dimName : string
        The name of the dimension to add (e.g. 'UMI', 'seqId', 'string1')
        readName : string
        A name for the sequence, which must be used across all dimensions (readNames)
        seqString : string
        The sequence that you want to store in this dimension
        stringDist : int
        The distance searched on this dimension
        """
        if not dimName in self.dimensions:
            self.dimensions[dimName] = {}
            self.dimension_distances[dimName] = stringDist
            self.total_dimensions += 1            
        if not readName in self.dimensions[dimName]:
            self.dimensions[dimName][readName] = seqString
            
    def define_nodes(self):
        """
        Collects strings across all dimensions, linking exact copies by name
        stores (seq1, seq2 ... seqN) as keys in a dictionary with the complete set of names as values
        """
        dimNames = list(self.dimensions.keys())
        readNames = list(self.dimensions[dimNames[0]].keys())                
        ##track read names in index
        for readName in readNames:            
            seqs = []
            for dimName in dimNames:
                seqs.append((dimName, self.dimensions[dimName][readName]))
            seqs = tuple(seqs)
            ## nodes contain all reads for exact matching dimensions
            if not seqs in self.nodes:
                self.nodes[seqs] = set()
            self.nodes[seqs].add(readName)    
        ##sort nodes by size
        self.nodes = collections.OrderedDict(
            sorted(self.nodes.items(), key = lambda x: len(x[1]), reverse=True)
        )
            
    def build_trees(self):
        """
        Build a BK tree for each dimension
        """
        for dimName in self.dimensions:
            dimensionSeqs = set(self.dimensions[dimName].values())
            self.BKtrees[dimName] = BKtree(iter(dimensionSeqs), Levenshtein.distance, 0)

    def build_index(self):
        """
        Create an index to cross-reference sequences with their names
        Additionally, store the sizes of each node in a dictionary
        """
        idx = 0
        for node in self.nodes:            
            self.node_sizes[idx] = len(self.nodes[node])
            self.index[node] = idx
            ## point from index value back to reads
            self.index_reads[idx] = self.nodes[node]
            for seqTup in node:
                dimName, seq = seqTup
                if not dimName in self.index:
                    self.index[dimName] = {}                
                if not seq in self.index[dimName]:
                    self.index[dimName][seq] = set()
                self.index[dimName][seq].add(idx)                
            idx += 1
        return 0

    #################################
    ## finding edges with single core
    #################################
    def find_edges(self):
        """
        Identify all edges within a specified search distance for each dimension

        Parameters
        ----------
        search_dist : int
        The Levenshtein edit distance allowed during the BKtree search over each dimension
        """
        for dimName in self.dimensions:
            self.edges[dimName] = {}
            tree = self.BKtrees[dimName]
            dimensionSeqs = set(self.dimensions[dimName].values())
            searchDist = self.dimension_distances[dimName]
            for seq in dimensionSeqs:
                similarSeqs = tree.find(seq, searchDist)
                self.edges[dimName][seq] = set.union(*[self.index[dimName][x] for x in similarSeqs]) 

    ## parallel implementation
    # def find_edges(self):
    #     """
    #     Identify all edges within a specified search distance for each dimension

    #     Parameters
    #     ----------
    #     search_dist : int
    #     The Levenshtein edit distance allowed during the BKtree search over each dimension
    #     """
    #     for dimName in self.dimensions:
    #         self.edges[dimName] = {}
    #         tree = self.BKtrees[dimName]
    #         dimensionSeqs = set(self.dimensions[dimName].values())
    #         searchDist = self.dimension_distances[dimName]
    #         ##break set into chunks for threads
    #         n = int(len(dimensionSeqs) / self.processors) + 2        
    #         seqGroups = [list(dimensionSeqs)[i * n:(i + 1) * n] for i in range((len(dimensionSeqs) + n - 1) // n )]
    #         with multiprocessing.Pool(processes=self.processors) as pool:
    #             res = pool.imap_unordered( self.parallel_search_tree, ( (seqs, tree, searchDist) for seqs in seqGroups) )
    #             edges = {}
    #             ##store index values as edges
    #             for item in res:
    #                 for seq in item:
    #                     self.edges[dimName][seq] = set.union(*[self.index[dimName][x] for x in item[seq]])
          
    # def parallel_search_tree(self, args):
    #     """
    #     Search for similar sequences in BKtree, in parallel
    #     Parameters
    #     ----------
    #     args : tuple
    #     seqs chunks, BKtree, search distance parameter, and the name of 
        
    #     """
    #     (seqs, tree, searchDist) = args
    #     res = {}
    #     for seq in seqs:                        
    #         similarSeqs = tree.find(seq, searchDist)
    #         res[seq] = similarSeqs
    #     return res
        
    def build_graph(self):
        """
        Construct the graph, build a plane in multi-dimensional space, connecting reads with similarities across all dimensions
        """
        self.define_nodes()
        self.build_index()
        self.define_centroids()
        self.build_trees()
        self.find_edges()
        ## join edges across all dimensions
        dimNames = list(self.dimensions.keys())        
        for node in self.nodes:
            nodeIdx = self.index[node]            
            ##intersect all edges from each dimension
            all_dimensions = []
            for i in range(len(dimNames)):
                all_dimensions.append(self.edges[dimNames[i]][node[i][1]])            
            self.graph[nodeIdx] = set.intersection(*all_dimensions) 
            ## the node should not have an edge to itself
            self.graph[nodeIdx] = self.graph[nodeIdx].difference(set([nodeIdx]))
              
    def walk_graph(self, node, max_steps):
        """
        Traverse all edges in decending order, based on node size.

        Parameters
        ----------
        node : int
        a node on the graph
        
        Returns
        -------
        walked : All the nodes found and the distance (number of edges) between them
        """
        step = 0
        walked = {}
        #walked[node] = 0
        edges = self.graph[node]
        while True:
            ##limit the walking
            if step >= max_steps:
                break
            step +=1
            ##avoid redundancy as we traverse the graph
            edges = set([x for x in edges if not x in walked])
            # return the edge found and the distance to each edge
            new_edges = set()
            for edge in edges:
                if edge == node:
                    continue
                ##only walk in descending order
                if not edge in walked:
                    walked[edge] = step
                ##do not continue if next node is larger
                if self.node_sizes[edge] > self.node_sizes[node]:
                    continue
                for new_edge in set(self.graph[edge]).difference(edges):
                    new_edges.add(new_edge)
            if len(new_edges) == 0:
                break
            else:
                edges = new_edges
        return walked

    def define_centroids(self):
        """
        determine centroids from node sizes by calculating an efficiency of data representation
        Efficiency is calculated as Work_output/Work_input. In terms of the graph, Work_output is defined as the 
        Cummulative fraction of reads contained in N nodes and Work_input is defined as the square root of N nodes 
        considered. 
        The most efficient representation of data is then found by maximizing this calculation. 
        The total number of N nodes at the maximum are centroids and are considered read events.

        """
        node_weights = [self.node_sizes[x] for x in self.node_sizes]
        node_weights = sorted(node_weights, reverse=True)
        node_stats = []
        cummulative_total = 0
        total = sum(node_weights)
        count = 0
        for weight in node_weights:
            count += 1
            cummulative_total += weight
            fraction = cummulative_total/total
            efficiency = fraction / math.sqrt(count)
            stats = {"weight":weight,
                     "cummulative_total":cummulative_total,
                     "fraction":fraction,
                     "efficiency": efficiency,
                     "count":count}
            node_stats.append(stats)
        node_stats = sorted(node_stats, key = lambda x: (-x["efficiency"], x["count"]))
        self.total_centroids = node_stats[0]["count"]
        ##if not definition found, set to zero
        if self.total_centroids >= 0.5*len(self.nodes):
            self.total_centroids = 0 
        
            
    def remove_edges(self, max_steps, minBiddingRatio):
        """
        Network partitioning steps.
        - Tranverses a graph comparing the weight between connected nodes (1-edit distance edges.)
        - Define Centroids as True nodes
        - remove edges between centroids
        - bid for smaller nodes (sequencing error)
        """

        ###############################
        ## calculate weights for nodes
        ###############################
        bids = {}    
        centroids = set()
        ## nodes are sorted in order
        node_count=0
        for node in self.graph:
            node_count += 1
            node_size = self.node_sizes[node]
            if node_count < self.total_centroids and node_size > 1:
                centroids.add(node)
            walked = self.walk_graph(node, max_steps)
            for edge in walked:
                edge_size = self.node_sizes[edge]                
                steps = walked[edge]                
                weight = (node_size**(1/steps))/edge_size
                bid = {"bidding_node":node,
                       "edge":edge,
                       "node_size":node_size,
                       "steps":steps,
                       "weight":weight}
                if not edge in bids:
                    bids[edge] = []
                bids[edge].append(bid)

        #################################
        ## accept offers, breaking edges
        #################################
        accepted = {}
        nodes = list(self.graph.keys())        
        for node in nodes:            
            ## centroids are considered real, they cannot accept offers
            if node in centroids:
                continue
            ## does the node have any bids for it
            if not node in bids:
                continue
            ##if node has already merged, continue
            if not node in self.graph:
                continue
            ##point to all edges
            allEdges = self.graph[node]
            ##get all offers
            offers = bids[node]
            offers = sorted(offers, key = lambda x : -x["weight"])
            for bid in offers:
                ## is the bid high enough to be accepted            
                if bid["weight"] >= minBiddingRatio:                                
                    accepted[node] = bid["bidding_node"]
                    ##if the bidding node has accepted an offer, look at other offers
                    if bid["bidding_node"] in accepted:
                        continue
                    ##commit node to bidding node
                    self.graph[node] = set([bid["bidding_node"]])
                    self.graph[bid["bidding_node"]].add(node)
                    ##remove node from all other edges
                    for edge in allEdges:                    
                        if not edge == bid["bidding_node"]:
                            self.graph[edge] = self.graph[edge].difference(set([node]))
                    break
                
        ##########################################                   
        ## remove edges between any two centroids
        ##########################################
        for node in self.graph:
            if node in centroids: 
                rmCents = self.graph[node].difference(centroids)                
                rmCents.add(node)
                self.graph[node] = rmCents                        
                                 
    def collapse_graph(self):
        """
        Traverse all edges in decending order, based on node size.
        walk all edges connected to each node
        
        Returns
        -------
        families : dict
        the collapsed graph, a dictionary with the largest node as a key and all connected smaller nodes as values
        """
        families = {}
        collapsed = {}
        for node in self.graph:
            if node in collapsed:
                continue
            else:
                collapsed[node] = 0
            step = 0            
            walked = set()
            #walked[node] = 0
            edges = self.graph[node]
            while True:                
                ##avoid redundancy as we traverse the graph, don't walk back
                edges = edges.difference(walked)
                # return the edge found and the distance to each edge
                new_edges = set()
                for edge in edges:
                    if edge == node:
                        continue
                    walked.add(edge)
                    for new_edge in set(self.graph[edge]).difference(edges):
                        new_edges.add(new_edge)
                if len(new_edges) == 0:
                    break
                else:
                    edges = new_edges
            walked.add(node)
            families[node] = walked
            for edge in walked:
                collapsed[edge] = 0
        return families

    def retreive_read_labels(self, families):
        """
        For a given set of families, return a dictionary of reads and read labels
        
        Parameters
        ----------
        families : dict
        A collapsed graph, where all values are the clusters against the key.
        """
        read_labels = {}
        family_idx = 0
        for node in families:
            family_idx += 1
            read_group = set()
            for idx in families[node]:
                read_group = read_group.union(self.index_reads[idx])
            for read in read_group:
                read_labels[read] = family_idx
        return read_labels
    
def main():
    """
    The following code serves as an example only, this package is not meant to be run on the command line, 
    but instead to be interfaced at the object level only.
    """
    
    print("read file")
    
    aln_df = pandas.read_csv("/sc1/groups/pls-redbfx/pipeline_runs/iPETE/production/2019-09-10-Ipete_redesign_Expt2/analysis/HUT78_100ng_10UMI_i7_nobc/filtered_reads/HUT78_100ng_10UMI_i7_nobc_on_target.tsv", sep="\t")

    
    ## minimum required columns
    ## feed name, seq, qual for each dimension
    graph = SeqNetwork(2, 2, 2)

    count = 0
    print("add dimensions")
    for name, umi1, cdr3, read_qual in zip(aln_df["name"], aln_df['uid'],
                                                 aln_df["cdr3"], aln_df["avg_read_qual"]):
        if read_qual < 30:
            continue
        if isinstance(cdr3, str) and isinstance(umi1, str):
            if len(cdr3) < 100:
                count += 1
                graph.add_dimension("UMI1", name, umi1, 1)
                ##can also add UMI2
                #graph.add_dimension("UMI2", name, umi2, 1)
                graph.add_dimension("CDR3", name, cdr3, 1)

                
    print("build graph")
    ##specify edit distances
    start = time.time()
    graph.build_graph()
    end = time.time()
    print("time", end-start)

    print("remove edges")
    start = time.time()
    graph.remove_edges()
    end = time.time()
    print("time", end-start)

    print("collapse_graph")
    families = graph.collapse_graph()
    
    for line in graph.summarize_graph():
        print(line)
    print("families", len(families))

    
if __name__ == '__main__':
    #arguments = make_arg_parser().parse_args()
    #logger = getLogger(
    #    arguments.alignment_summary.replace(".tsv", ""), debug=False)
    main()
