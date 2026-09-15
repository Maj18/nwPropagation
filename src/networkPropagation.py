import pandas as pd
import numpy as np
from tqdm.notebook import tqdm
from rwr_functions import *
from constants import *
import json
import os
from sknetwork.data import from_edge_list


class networkPropagation():
    """
    NetworkProgation
    """
    
    """
    Initiate class.
    Parameters:
        feature_pval_file: a processed data file, e.g. from differential analysis, it needs to have at least 3 columns: ´Pvalue´, ´NCBI_id´, ´Error´
    """
    def __init__(self, feature_pval_file):
        PACKAGE_ROOT = Path(__file__).parent
        network_dir = f"{PACKAGE_ROOT}/data/networks/"
        self.feature_pval_file = feature_pval_file
        with open(f"{network_dir}/networks/networks_n_edges.json", "r") as f:
            networks_n_edges = json.load(f)
            
        with open(f"{network_dir}/networks/networks_n_nodes.json", "r") as f:
            networks_n_nodes = json.load(f)
        
        network_names = networks_n_edges.keys()
        self.network_files = [
            f"{network_dir}/{network}/edges_list_ncbi.csv"
            for network in network_names]
        self.networks_n_nodes = networks_n_edges
        self.networks_n_edges = networks_n_edges
        self.NETWORKS = network_names
        
        
    """
    Load graph
    Parameters:
        network: a network name
    """
    def load_graph_nx(self, network):
        print("Loading {} graph - 7901_NetworkPropagation.py:46".format(network))
        df = pd.read_csv(network_files[network], dtype={'node1': str, 'node2': str})[["node1", "node2"]] ######
        graph = from_pandas_edgelist(df, source="node1", target="node2")
        # graph["names"] = graph["names"].astype(str)
        return graph

    """
    Load your processed data
    """
    def load_pegasus_results(self):
        data = pd.read_csv(self.feature_pval_file)
        data = data[~data["NCBI_id"].isna()]
        data["NCBI_id"] = data["NCBI_id"].astype(str)
        # data[data["NCBI_id"].isin(graph["names"])].sort_values(by="Pvalue")
        # min_pval = data["Pvalue"][data["Pvalue"]>0].min()
        # print(min_pval)
        pegasus_scores = {}
        for i, row in data.iterrows():
            # pv = np.maximum(min_pval, np.minimum(1, row["Pvalue"]))
            pv = row["Pvalue"] if row["Pvalue"]>0.0 else row["Error"]
            pegasus_scores[row["NCBI_id"]] = np.maximum(1e-16, -np.log10(pv))
        # pegasus_ncbi_genes = set(pegasus_scores.keys())
        data["Score"] = data["NCBI_id"].map(pegasus_scores)
        return data

    """
    Initiate graph node scores: -log10(Pvalue) if node in your data, otherwise 0
    """
    def init_rwr_scores_nx(self, network):
        graph = load_graph_nx(self, network=network)
        data = load_pegasus_results(self)
        node2idx = {str(n): i for i, n in enumerate(graph.nodes)}
        # idx2node = {v: k for k, v in node2idx.items()}
        ncbi2gene = dict(zip(data.NCBI_id, data.Gene))
        ncbi_genes = set(data.NCBI_id)
        pegasus_scores = dict(zip(data.NCBI_id, data.Score))
        pagerank_seeds = {}
        for node in graph.nodes: #["names"].astype(str):
            if node in ncbi_genes:
                pagerank_seeds[node] = pegasus_scores[node]
            else:
                pagerank_seeds[node] = 0
        return graph, data, pagerank_seeds


    """
    Perform network propagation
    Parameters:
        alpha: For a directed graph the degree shortcut no longer holds, so we find that keeps the walk ergodic. 
            The alpha is the teleportation probability in the PageRank random walk. The teleport prevents the walk from getting trapped in cycles or disconnected parts of the graph, 
            helping ensure a unique stationary distribution.
    """
    def perform_rwr_nx(self, alpha, network):
        graph, data, seeds = init_rwr_scores_nx(self, network=network)
        rwr_scores = pagerank(graph, alpha=alpha, personalization=seeds)
        return rwr_scores, graph, data, seeds, alpha

    """
    Process the network propagation results
    """
    def process_rwr_results_nx(self, alpha, network):
        scores, graph, data, seeds, alpha = perform_rwr_nx(self, alpha, network)
        node2idx = {str(n): i for i, n in enumerate(graph.nodes)}
        idx2node = {v: k for k, v in node2idx.items()}
        ncbi2gene = dict(zip(data.NCBI_id, data.Gene))
        seeds_vals = np.fromiter(seeds.values(), dtype="float")
        max_val = np.max(seeds_vals[~np.isinf(seeds_vals)])
        rwr_results = []
        for i, node in enumerate(graph.nodes):
            row = {}
            row["Idx"] = node2idx[node]
            row["Gene NCBI ID"] = node
            row["Symbol"] = ncbi2gene[node] if node in ncbi2gene.keys() else "-" #####????
            init_score = seeds[node]
            if np.isinf(init_score):
                init_score = max_val
            row["Initial Score"] = init_score
            row["Final Score"] = scores[node]
            rwr_results.append(row)
            
        rwr_results = pd.DataFrame(rwr_results).sort_values(by="Final Score", ascending=False)
        return scores, graph, data, seeds, alpha, rwr_results
    
    """
    Calculate metrics
    """  
    def calculate_metrics(rwr_res, K, gene_seeds, targets, net_name, alpha, disease, scoring="Score"):
        results = []
        genes = [str(s) for s in rwr_res["Gene NCBI ID"] if str(s) not in gene_seeds] #data[col].astype(str)
        # pak = precision_at_k(targets, genes, K)
        # apk = average_precision_at_k(targets, genes, K)
        apk = calc_apk(targets, genes, k_max=K)
        # results.append({"Network": net_name, "Alpha": alpha, "Metric": "Precision", "K": K, "Value": pak, "Method": scoring, "Disease": disease})
        results.append({"Network": net_name, "Alpha": alpha, "Metric": "Average Precision", "K": K, "Value": apk, "Method": scoring, "Disease": disease})
        results = pd.DataFrame(results)
        return results
    
    """
    Rank the nodes based on the network propagation results
    """   
    def rankNode(self, alpha, disease="Disease", Ks=[10, 20, 50]):
        genes_ranks = {}
        genes_scores = {}
        for netname in self.NETWORKS:
            scores, graph, data, seeds, alpha, rwr_results = process_rwr_results_nx(self, alpha=alpha, network=netname) 
            rank = 0  
            score_norm = np.sum(list(seeds.values())) * 0.01
            n_nodes = self.networks_n_nodes[netname]
            for i, row in rwr_results.iterrows():
                gn = str(row["Gene NCBI ID"])
                rank += 1
                if gn not in genes_ranks:
                    genes_ranks[gn] = [rank/n_nodes]
                else:
                    genes_ranks[gn].append(rank/n_nodes)

                if gn not in genes_scores:
                    genes_scores[gn] = [row["Final Score"]/score_norm]
                else:
                    genes_scores[gn].append(row["Final Score"]/score_norm)
                    
        genes_ranks_df = []
        for k, v in genes_ranks.items():
            if len(v)<2:
                continue
            genes_ranks_df.append({"Gene NCBI ID": k, "Avg. Rank": np.mean(v), "Method": "Score", "Disease": disease, "Alpha": alpha})
            
        genes_ranks_df = pd.DataFrame(genes_ranks_df).sort_values(by="Avg. Rank", ascending=True)
        genes_scores_df = []
        for k, v in genes_scores.items():
            if len(v)<2:
                continue
            genes_scores_df.append({"Gene NCBI ID": k, "Avg. Score": np.mean(v), "Method": "Score", "Disease": disease, "Alpha": alpha})
            
        genes_scores_df = pd.DataFrame(genes_scores_df).sort_values(by="Avg. Score", ascending=False)
        for K in Ks:
            metrics = calculate_metrics(genes_ranks_df, K, self.gene_seeds_ncbi, self.ncbi_targets, "Multilayer", alpha, disease, scoring="Avg. Rank")
            results.append(metrics)
            metrics = calculate_metrics(genes_scores_df, K, self.gene_seeds_ncbi, self.ncbi_targets, "Multilayer", alpha, disease, scoring="Avg. Score")
            results.append(metrics)
            
        return results

    """
    Run network propagation on mulitple alpha values
    """  
    def NetworkPropagationWmultiAlpha(self, ALPHAS, outdir, disease="Disease", Ks=[10, 20, 50]):
        for alpha in tqdm(ALPHAS):
            results = rankNode(self, alpha=alpha, disease=disease, Ks=Ks)
            
        results = pd.concat(results)
        os.makedirs(outdir, exist_ok=True) 
        results.to_csv(f"{outdir}/AVG_ensembles_metrics.csv", index=False)
        return results

