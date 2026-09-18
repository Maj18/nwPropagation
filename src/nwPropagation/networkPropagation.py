import pandas as pd
import numpy as np
# from tqdm.notebook import tqdm
from tqdm import tqdm
import json
import os
# from sknetwork.data import from_edge_list
from pathlib import Path
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sknetwork.ranking import PageRank
from networkx import pagerank, from_pandas_edgelist
import networkx as nx
import pickle
import gseapy as gp
import matplotlib.lines as mlines
import seaborn as sns
import mygene


class networkPropagation():
    """
    Differential-expression-weighted network propagation analysis.
    """
    
    """
    Initiate class.
    Parameters:
        feature_pval_file: a processed data file, e.g. from differential analysis, it needs to have at least 2 columns: ´Pvalue´, ´NCBI_id´, ´Uniprot´, ´Gene´, ´logFC´
    """
    def __init__(self, feature_pval_file, logFCcutoff=0.3, adjPcutoff=0.05):
        PACKAGE_ROOT = Path(__file__).parent
        print(PACKAGE_ROOT)
        network_dir = f"{PACKAGE_ROOT}/../../data/networks/"
        self.feature_pval_file = feature_pval_file
        with open(f"{network_dir}/networks_n_edges.json", "r") as f:
            networks_n_edges = json.load(f)
            
        with open(f"{network_dir}/networks_n_nodes.json", "r") as f:
            networks_n_nodes = json.load(f)
        
        network_names = networks_n_edges.keys()
        self.network_files = {
            network: f"{network_dir}/{network}/edges_list_ncbi.csv"
            for network in network_names}
        self.networks_n_nodes = networks_n_nodes
        self.networks_n_edges = networks_n_edges
        self.NETWORKS = network_names
        # self.known_disease_genes_ncbi = known_disease_genes_ncbi
        self.logFCcutoff = logFCcutoff
        self.adjPcutoff = adjPcutoff
        
        
    """
    Load graph
    Parameters:
        network: a network name
    """
    def load_graph_nx(self, network):
        print("Loading {} graph  7901_NetworkPropagation.py:51 - networkPropagation.py:63".format(network))
        print(self.network_files)
        df = pd.read_csv(self.network_files[network], dtype={'node1': str, 'node2': str})[["node1", "node2"]]
        graph = from_pandas_edgelist(df, source="node1", target="node2")
        # graph["names"] = graph["names"].astype(str)
        return graph

    """
    Load your processed data
    """
    def load_pegasus_results(self, FILTER=True):
        data = pd.read_csv(self.feature_pval_file, delimiter="\t", 
                dtype={'NCBI_id': str, 'Gene': str, 'Uniprot': str})
        data = data[~data["NCBI_id"].isna()]
        data["NCBI_id"] = data["NCBI_id"].astype(str)
        data = data.sort_values(by="Pvalue")
        data = data[
            data["NCBI_id"].notna()
        ].drop_duplicates(
            subset="NCBI_id"
        ).copy()
        if FILTER==True:
            data = data[
                (abs(data["logFC"]) > self.logFCcutoff) &
                (data["adj.P.Val"] < self.adjPcutoff)].copy()
            
        # data[data["NCBI_id"].isin(graph["names"])].sort_values(by="Pvalue")
        # min_pval = data["Pvalue"][data["Pvalue"]>0].min()
        # print(min_pval)
        pegasus_scores = {}
        for i, row in data.iterrows():
            # pv = np.maximum(min_pval, np.minimum(1, row["Pvalue"]))
            pv = row["Pvalue"] if row["Pvalue"]>0.0 else pd.NA
            pegasus_scores[row["NCBI_id"]] = np.maximum(1e-16, -np.log10(pv))
        # pegasus_ncbi_genes = set(pegasus_scores.keys())
        data["Score"] = data["NCBI_id"].map(pegasus_scores)
        return data

    """
    Initiate graph node scores: -log10(Pvalue) if node in your data, otherwise 0
    """
    def init_rwr_scores_nx(self, network):
        graph = self.load_graph_nx(network=network)
        data = self.load_pegasus_results(FILTER=True)
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
        graph, data, seeds = self.init_rwr_scores_nx(network=network)
        rwr_scores = nx.pagerank(graph, alpha=alpha, personalization=seeds)
        return rwr_scores, graph, data, seeds, alpha

    """
    Process the network propagation results
    """
    def process_rwr_results_nx(self, alpha, network):
        scores, graph, data, seeds, alpha = self.perform_rwr_nx(alpha, network)
        node2idx = {str(n): i for i, n in enumerate(graph.nodes)}
        idx2node = {v: k for k, v in node2idx.items()}
        ncbi2uniprot = dict(zip(data.NCBI_id, data.Uniprot))
        ncbi2gene = dict(zip(data.NCBI_id, data.Gene))
        seeds_vals = np.fromiter(seeds.values(), dtype="float")
        max_val = np.max(seeds_vals[~np.isinf(seeds_vals)])
        rwr_results = []
        for i, node in enumerate(graph.nodes):
            row = {}
            row["Idx"] = node2idx[node]
            row["Gene NCBI ID"] = node
            row["Uniprot"] = ncbi2uniprot[node] if node in ncbi2uniprot.keys() else "-"
            row["Symbol"] = ncbi2gene[node] if node in ncbi2gene.keys() else "-"
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
    def calc_apk(self, y_true, y_pred, k_max=0):
        # Check if all elements in lists are unique
        if len(set(y_true)) != len(y_true):
            raise ValueError("Values in y_true are not unique")
        if len(set(y_pred)) != len(y_pred):
            raise ValueError("Values in y_pred are not unique")
        if k_max != 0:
            y_pred = y_pred[:k_max]

        correct_predictions = 0
        running_sum = 0
        for i, yp_item in enumerate(y_pred):
            k = i+1 # our rank starts at 1       
            if yp_item in y_true:
                correct_predictions += 1
                running_sum += correct_predictions/k

        if correct_predictions==0:
            return 0.0
        
        return running_sum/min(len(y_true), k_max)
        #return running_sum/correct_predictions

    def calculate_metrics(self, rwr_res, K, gene_seeds, targets, net_name, alpha, disease, scoring="Score"):
        results = []
        genes = [str(s) for s in rwr_res["Gene NCBI ID"] if str(s) not in gene_seeds] #data[col].astype(str)
        # pak = precision_at_k(targets, genes, K)
        # apk = average_precision_at_k(targets, genes, K)
        apk = self.calc_apk(targets, genes, k_max=K)
        # results.append({"Network": net_name, "Alpha": alpha, "Metric": "Precision", "K": K, "Value": pak, "Method": scoring, "Disease": disease})
        results.append({"Network": net_name, "Alpha": alpha, 
                "Metric": "Average Precision", "K": K, "Value": apk, "Method": scoring, "Disease": disease})
        results = pd.DataFrame(results)
        return results
    
    """
    Rank the nodes based on the network propagation results
    """   
    def rankNode(self, alpha, disease="Disease", outDIR="./"):
        genes_ranks = {}
        genes_scores = {}
        for netname in self.NETWORKS:
            scores, graph, data, seeds, alpha, rwr_results = self.process_rwr_results_nx(alpha=alpha, network=netname) 
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
            genes_ranks_df.append({"Gene NCBI ID": k, 
                "Avg. Rank": np.mean(v), "Method": "Score", "Disease": disease, "Alpha": alpha})
            
        genes_ranks_df = pd.DataFrame(genes_ranks_df).sort_values(by="Avg. Rank", ascending=True)
        with open(f"{outDIR}/genes_ranks_df.pkl", "wb") as file:
            pickle.dump(genes_ranks_df, file)
        
        genes_ranks_df2 = genes_ranks_df.copy()
        genes_ranks_df2.columns.values[0] = "NCBI_id"
        data_orig = self.load_pegasus_results(FILTER=False)
        print(data_orig.shape)
        for df in [data_orig, genes_ranks_df2]:
            df.columns = df.columns.astype(str).str.strip()
            
        print(data_orig.columns)
        print(genes_ranks_df2.columns)
        combined = (
            genes_ranks_df2
            .merge(data_orig, on="NCBI_id", how="outer")
        ).sort_values(by="Avg. Rank", ascending=True)        
        combined.to_csv(f"{outDIR}/genes_ranks_df.csv", index=False)
                  
        genes_scores_df = []
        for k, v in genes_scores.items():
            print("k: Gene - networkPropagation.py:251")
            print("v: Propagation scores from each network, here we require the edges need to be supported by at least 2 networks - networkPropagation.py:252")
            if len(v)<2:
                continue
            
            genes_scores_df.append({"Gene NCBI ID": k, 
                "Avg. Score": np.mean(v), "Method": "Score", "Disease": disease, "Alpha": alpha})
            
        genes_scores_df = pd.DataFrame(genes_scores_df).sort_values(by="Avg. Score", ascending=False)
        with open(f"{outDIR}/genes_scores_df.pkl", "wb") as file:
            pickle.dump(genes_scores_df, file)
        genes_scores_df2 = genes_scores_df.copy()
        genes_scores_df2.columns.values[0] = "NCBI_id"
        for df in [data_orig, genes_scores_df2]:
            df.columns = df.columns.astype(str).str.strip()
            
        print(data_orig.columns)
        print(genes_scores_df2.columns)
        combined2 = (
            genes_scores_df2
            .merge(data_orig, on="NCBI_id", how="outer")
        ).sort_values(by="Avg. Score", ascending=False)             
        combined2.to_csv(f"{outDIR}/genes_scores_df.csv", index=False)
        print("Evaluation metrics at different top node nr cutoff K ... - networkPropagation.py:274")
        
        return data, genes_ranks_df, genes_scores_df
        
    
    def metrics_all(self, alpha, disease="Disease", Ks=[10, 20, 50], 
        outDIR="./", target_genes=[], data=None, 
        genes_ranks_df=None, genes_scores_df=None):
        if data is None:
            data, genes_ranks_df, genes_scores_df =  self.rankNode(
                alpha=alpha,
                outDIR=outDIR)
            
        results = []
        for K in Ks:
            ncbi_seeds = data["NCBI_id"].dropna().astype(str).tolist()
            if len(target_genes)==0:
                data_orig = self.load_pegasus_results(FILTER=False)
                filtered = data_orig[(abs(data_orig["logFC"]) > 0.3)&(data_orig["adj.P.Val"] < 0.05)].copy()
                ncbi_targets = (
                    filtered["NCBI_id"]
                    .dropna()
                    .astype(str)
                    .tolist())
            else:
                ncbi_targets = target_genes
                
            # ncbi_targets = target_genes + ncbi_targets0
            print(ncbi_targets)
            # print(ncbi_targets.head())
            # rwr_res, K, gene_seeds, targets, net_name, alpha, disease, scoring="Score"
            metrics = self.calculate_metrics(genes_ranks_df, 
                    K, ncbi_seeds, targets=ncbi_targets, 
                    net_name="Multilayer", alpha=alpha, disease=disease, scoring="Avg. Rank")
            results.append(metrics)
            metrics = self.calculate_metrics(genes_scores_df, 
                    K, ncbi_seeds, targets=ncbi_targets, 
                    net_name="Multilayer", alpha=alpha, disease=disease, scoring="Avg. Score")
            results.append(metrics)
            
        return results

    """
    Run network propagation on mulitple alpha values
    """  
    def NetworkPropagationWmultiAlpha(
        self,
        outdir,
        ALPHAS=None,
        disease="Disease",
        Ks=[10, 20, 50, 100]
    ):
        if ALPHAS is None:
            ALPHAS = np.concatenate(
                (
                    [0.0, 0.005, 0.1],
                    [0.15, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
                    # np.linspace(0.001, 0.01, 3),
                    # np.linspace(0.02, 0.9, 8)
                )
            )

        all_results = []
        for alpha in tqdm(ALPHAS):
            os.makedirs(f"{outdir}/Alpha_{alpha}/", exist_ok=True)
            result =  self.metrics_all(alpha, 
                disease=disease, Ks=Ks, 
                outDIR=f"{outdir}/Alpha_{alpha}/")
            all_results.append(result)

        
        with open(f"{outdir}/all_results.pkl", "wb") as file:
            pickle.dump(all_results, file)
        
        # with open(f"{outdir}/all_results.pkl", "rb") as file:
        #     all_results = pickle.load(file)
           
        flat_results = [
            df
            for result_group in all_results
            for df in result_group
        ]
        results = pd.concat(flat_results, ignore_index=True) 
        output_file = os.path.join(
            outdir,
            "AVG_ensembles_metrics.csv"
        )
        results.to_csv(output_file, index=False)

        return results
    


"""
Utility functions   
""" 

def convert_identifiers(
    df_file,
    delimiter="\t",
    FeatureCol="Features",
    Col2drop=[],
    from_type="uniprot",
    to_fields="entrezgene,symbol",
    species="human"
    ):
    # 1. Load your original data
    df = pd.read_csv(df_file, delimiter=delimiter)
    if len(Col2drop)>0:
        df = df.drop(columns=Col2drop)
    # 2. Parse 'to_fields' into a clean list of target keys
    if isinstance(to_fields, str):
        target_fields = [f.strip() for f in to_fields.split(",")]
    else:
        target_fields = list(to_fields)
    # 3. Query MyGeneInfo
    mg = mygene.MyGeneInfo()
    result = mg.querymany(
        df[FeatureCol].tolist(),
        scopes=from_type,
        fields=target_fields,
        species=species,
    )
    # 4. Handle the mapping part dynamically
    mapping_data = []
    for item in result:
        # Base row linking back to the original 'Features' query
        row = {FeatureCol: item.get("query")}
        # Dynamically pull whatever keys the user requested
        for field in target_fields:
            # Safely handle missing hits or nested sub-fields
            row[field] = item.get(field, None)
        mapping_data.append(row)
    # 5. Build dataframe and drop duplicate query rows if any exist
    mapping_df = pd.DataFrame(mapping_data).drop_duplicates(subset=[FeatureCol])
    # 6. Merge back seamlessly
    df = df.merge(mapping_df, on=FeatureCol, how="left")
    df.rename(columns={"entrezgene": "NCBI_id", "symbol":"Gene"}, inplace=True)
    return df


def plot_metrics(df, outdir):
    # Force K to a string/category so matplotlib treats it as discrete shapes
    df["K"] = df["K"].astype(str)
    # 2. Set style and initialize a single figure
    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(9, 6))  # Strictly ONE plot window
    # 3. Plot everything onto the same single set of axes
    # 'hue' splits by Method color, 'style' splits by K shape
    sns.lineplot(
        data=df,
        x="Alpha",
        y="Value",
        hue="Method",
        style="K",
        markers=True,  # Gives each K its own unique geometric shape point
        dashes=False,  # Keeps lines solid so it is easy to track colors
        linewidth=2.5,
        markersize=9,
        palette="Set1",
        ax=ax,  # Forces it into this exact single plot
    )
    # 4. Visual formatting
    plt.title(
        "Evaluation Summary: Impact of Alpha, Method, and K",
        fontsize=14,
        weight="bold",
        pad=15,
    )
    plt.xlabel("Alpha (Hyperparameter)", fontsize=11, labelpad=10)
    plt.ylabel("Metric Value", fontsize=11, labelpad=10)
    # Cleanly position the legend outside the graph area so lines never hide behind it
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)
    g.savefig(f"{outdir}/evaluation_summary.png", dpi=300, bbox_inches="tight")


def runORA(gene_list, 
    gene_sets="GO_Biological_Process_2023", 
    organism="human", outDIR=f"../../test/Alpha_0.3/",
    padj_cutoff=0.05, top_term=20,
    title="Enriched pathways"):
    os.makedirs(outDIR, exist_ok=True)
    # 2. Run Over-Representation Analysis using Enrichr
    # 'GO_Biological_Process_2023' is the standard, up-to-date library name
    ora_result = gp.enrichr(
        gene_list=gene_list,
        gene_sets=gene_sets,
        organism=organism,  # Change to 'mouse' or 'fly' if applicable
        outdir=outDIR,  # Set a path string here if you want text files auto-saved
    )
    # 3. Extract the results into a Pandas DataFrame
    ora_df = ora_result.res2d
    # Sort by Adjusted P-value to see the most significant pathways first
    ora_df = ora_df.sort_values("Adjusted P-value")
    # Display the top 5 enriched terms
    print(ora_df[["Term - networkPropagation.py:469", "Overlap", "P-value", "Adjusted P-value", "Genes"]].head())
    # 4. Save the full results table to a CSV file
    ora_df.to_csv(f"{outDIR}/ora_results.csv", index=False)
    ora_df_filtered = ora_df[
        ora_df["Adjusted P-value"] < padj_cutoff].copy()
    sig_path_nr = ora_df_filtered.shape[0]
    max_path_name_length = ora_df_filtered['Term'].astype(str).str.len().max()
    # 5. Generate a publication-ready Dot/Bubble Plot
    # This automatically plots the top enriched pathways
    # print(2+(sig_path_nr*10))
    # plt.figure(figsize=(2+(max_path_name_length)/3, 2+(sig_path_nr*3)))
    gp.plot.dotplot(
        ora_result.res2d,
        title=title,
        cmap="viridis",  # Color gradient for significance
        cutoff=padj_cutoff,  # Only show terms with Adjusted P-value < 0.05
        top_term=top_term,  # Limit plot to top 10 terms
    )
    # Save the plot securely
    print(max_path_name_length)
    plt.gcf().set_size_inches(2+(max_path_name_length)*0.005, 1.75+(sig_path_nr*0.23)) 
    plt.savefig(f"{outDIR}/dotplot.png", 
        dpi=300, bbox_inches="tight")
    # plt.show()
    return ora_df, ora_result


def merge_graphs_with_provenance(G1, G2, net1, net2):
    # Create a new graph of the same type (Graph, DiGraph, etc.)
    G_merged = type(G1)()
    # 1. Process Nodes
    nodes_G1 = set(G1.nodes())
    nodes_G2 = set(G2.nodes())
    all_nodes = nodes_G1.union(nodes_G2)
    for node in all_nodes:
        # Transfer existing attributes from original graphs
        attrs = {}
        if node in G1:
            attrs.update(G1.nodes[node])
        if node in G2:
            attrs.update(G2.nodes[node])       
        # Add provenance/source metadata
        in_g1 = node in G1
        in_g2 = node in G2     
        if in_g1 and in_g2:
            source = net1 + net2
        elif in_g1:
            source = net1
        else:
            source = net2     
        attrs['source'] = source
        G_merged.add_node(node, **attrs)
    # 2. Process Edges
    edges_G1 = set(G1.edges())
    edges_G2 = set(G2.edges())
    all_edges = edges_G1.union(edges_G2)    
    for u, v in all_edges:
        attrs = {}
        if G1.has_edge(u, v):
            attrs.update(G1.get_edge_data(u, v))
        if G2.has_edge(u, v):
            attrs.update(G2.get_edge_data(u, v))            
        in_g1 = G1.has_edge(u, v)
        in_g2 = G2.has_edge(u, v)       
        if in_g1 and in_g2:
            source = [net1, net2]
        elif in_g1:
            source = net1
        else:
            source = net2          
        attrs['source'] = source
        G_merged.add_edge(u, v, **attrs)        
    return G_merged


def DecorateNode(G, propagation_file, delimiter=",", N=100):
    df = pd.read_csv(propagation_file, delimiter=delimiter)
    df = df.sort_values(by="Avg. Rank", ascending=True)
    up_sig = df[
        (df["logFC"] > 0.3) &
        (df["adj.P.Val"] < 0.05)].copy()
    up = list(set(
        up_sig["NCBI_id"]
        .dropna()
        .astype(str)
        .tolist()))
    down_sig = df[
        (df["logFC"] < -0.3) &
        (df["adj.P.Val"] < 0.05)].copy()
    down = list(set(
        down_sig["NCBI_id"]
        .dropna()
        .astype(str)
        .tolist()))
    df_filtered = df[
        (abs(df["logFC"]) > 1) &
        (df["adj.P.Val"] < 0.001)].copy()
    seeds = list(set(
        df_filtered["NCBI_id"]
        .dropna()
        .astype(str)
        .tolist()))
    propagated = list(set(
        df["NCBI_id"].iloc[0:+N].dropna().astype(str).tolist()
    ))
    subG = G.subgraph(propagated)
    nodes_G = set(subG.nodes())
    for node in nodes_G:
        # Transfer existing attributes from original graphs
        subG.nodes[node]['status'] = "na"
        if node in propagated:
            subG.nodes[node]['status'] = "propagated"
        if node in seeds:
            subG.nodes[node]['status'] = "seed"
        subG.nodes[node]['regulation'] = "na"
        if node in up:
            subG.nodes[node]['regulation'] = "up"
        if node in down:
            subG.nodes[node]['regulation'] = "down" 
    return subG     


def plotGraph(G, outDIR, figsize=(30, 30), legend_loc="upper left"):
    # 1. Color mapping based on regulation
    color_map = {
        "up": "pink",
        "down": "skyblue",
        "na": "snow",
    }
    # 2. Compute Layout Positions
    pos = nx.spring_layout(G, seed=42)
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    # 3. Draw Edges in 'gray' color
    nx.draw_networkx_edges(
        G, pos, edge_color="gray", alpha=0.8, width=1.5, ax=ax
    )
    # 4. Draw Nodes Grouped by Shape (Status)
    # --- Seed Nodes -> Square ('s') ---
    seed_nodes = [
        n for n, attr in G.nodes(data=True) if attr.get("status") == "seed"
    ]
    seed_colors = [
        color_map.get(G.nodes[n].get("regulation"), color_map["na"])
        for n in seed_nodes
    ]
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=seed_nodes,
        node_shape="s",
        node_color=seed_colors,
        edgecolors="orange",
        linewidths=2.0,
        node_size=2000,
        ax=ax,
    )
    # --- Propagated / Other Nodes -> Circle ('o') ---
    propagated_nodes = [
        n for n, attr in G.nodes(data=True) if attr.get("status") != "seed"
    ]
    propagated_colors = [
        color_map.get(G.nodes[n].get("regulation"), color_map["na"])
        for n in propagated_nodes
    ]
    nx.draw_networkx_nodes(
        G,
        pos,
        nodelist=propagated_nodes,
        node_shape="o",
        node_color=propagated_colors,
        edgecolors="orange",
        linewidths=2.0,
        node_size=2000,
        ax=ax,
    )
    # 5. Draw Labels in Black
    nx.draw_networkx_labels(
        G, pos, font_color="black", font_weight="bold", font_size=10, ax=ax
    )
    # 6. Custom Legend for Shapes and Colors
    legend_elements = [
        mlines.Line2D(
            [],
            [],
            color="none",
            marker="s",
            markerfacecolor="gray",
            markeredgecolor="orange",
            markersize=12,
            label="Seed Node (Square)",
        ),
        mlines.Line2D(
            [],
            [],
            color="none",
            marker="o",
            markerfacecolor="gray",
            markeredgecolor="orange",
            markersize=12,
            label="Propagated Node (Circle)",
        ),
        mlines.Line2D(
            [],
            [],
            color="none",
            marker="o",
            markerfacecolor="pink",
            markeredgecolor="orange",
            markersize=12,
            label="Up-regulated (Red)",
        ),
        mlines.Line2D(
            [],
            [],
            color="none",
            marker="o",
            markerfacecolor="skyblue",
            markeredgecolor="orange",
            markersize=12,
            label="Down-regulated (Blue)",
        ),
        mlines.Line2D(
            [],
            [],
            color="none",
            marker="o",
            markerfacecolor="snow",
            markeredgecolor="orange",
            markersize=12,
            label="NA / Other (Snow)",
        ),
    ]
    legend = ax.legend(
        handles=legend_elements, loc=legend_loc, title="Legend", frameon=True
    )
    legend.get_frame().set_facecolor("#3B3B3B")
    legend.get_frame().set_edgecolor("grey")
    plt.setp(legend.get_texts(), color="white")
    plt.setp(legend.get_title(), color="white")
    plt.title(
        "Gene Propagation Network", fontsize=16, fontweight="bold", color="black"
    )
    plt.axis("off")
    # 7. Save and Show
    plt.savefig(
        f"{outDIR}/network_propagation.png",
        bbox_inches="tight",
        dpi=300,
        facecolor=fig.get_facecolor(),
    )
    plt.show()
    plt.close()
        
