import os
import gseapy as gp
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def perform_deg_enrichment(adata, annotation_column='annotation', 
                            gene_sets='MSigDB_Hallmark_2020', 
                            organism='Human', pvals_adj_thr = 0.05,
                            logFC_thr = 0.5, num_input_genes = 10,
                            compare_group = 'IPF UL',
                            reference_group = 'the rest'):
    print(f"Performing differential expression analysis for {annotation_column}...")
    key_added = f"{annotation_column}_{compare_group}_vs_{reference_group}"
    print(f"{compare_group} vs {reference_group}")
    if reference_group == 'the rest':
        sc.tl.rank_genes_groups(adata, groupby=annotation_column, method='wilcoxon',  groups=[compare_group], key_added=key_added)
    else:
        sc.tl.rank_genes_groups(adata, groupby=annotation_column, method='wilcoxon',  groups=[compare_group],reference=reference_group, key_added=key_added)

    deg_df = sc.get.rank_genes_groups_df(adata, group=compare_group, key=key_added)
    

    print(f"Performing Enrichr for {gene_sets}...")
    print(f"Enrichr gene selection criteria: pvals_adj < {pvals_adj_thr} & |logfoldchanges| > {logFC_thr}")

    pos_df = deg_df.query(f'pvals_adj < {pvals_adj_thr} & logfoldchanges > {logFC_thr}') # Filtering criteria
    pos_gene_list = pos_df['names'].tolist()

    pos_enr_df, neg_enr_df = None, None
    if len(pos_gene_list) > num_input_genes:
        pos_enr = gp.enrichr(
            gene_list=pos_gene_list,
            gene_sets=gene_sets,
            organism=organism,
            cutoff=0.05
        )
        pos_enr_df = pos_enr.results

    neg_df = deg_df.query(f'pvals_adj < {pvals_adj_thr} & logfoldchanges < -{logFC_thr}') # Filtering criteria
    neg_gene_list = neg_df['names'].tolist()
    if len(neg_gene_list) > num_input_genes:
        neg_enr = gp.enrichr(
            gene_list=neg_gene_list,
            gene_sets=gene_sets,
            organism=organism,
            cutoff=0.05
        )
        neg_enr_df = neg_enr.results

    
    return deg_df, pos_df, neg_df, pos_enr_df, neg_enr_df


def plot_volcano(
    df,
    lfc_thr=1.0,
    padj_thr=0.05,
    top_n=10,
    figsize=(10, 6),
    title=None,
    save_path=None,
    show=False
):
    """
    df: DEG dataframe with columns ['names', 'pvals_adj', 'logfoldchanges']
    """
    if df is None or df.empty:
        return

    df = df.dropna(subset=['pvals_adj', 'logfoldchanges']).copy()
    if df.empty:
        return

    df['-log10(padj)'] = -np.log10(df['pvals_adj'].astype(float) + 1e-300)

    df['sig'] = 'Not Sig'
    df.loc[(df['logfoldchanges'] >= lfc_thr) & (df['pvals_adj'] < padj_thr), 'sig'] = 'Up'
    df.loc[(df['logfoldchanges'] <= -lfc_thr) & (df['pvals_adj'] < padj_thr), 'sig'] = 'Down'

    top_up = df.query('sig == "Up"').sort_values('pvals_adj').head(top_n)
    top_down = df.query('sig == "Down"').sort_values('pvals_adj').head(top_n)
    top_genes = pd.concat([top_up, top_down], axis=0)

    plt.figure(figsize=figsize)
    palette = {'Up': '#D43F3AFF', 'Down': '#4477AFFF', 'Not Sig': '#CCCCCC80'}

    sns.scatterplot(
        data=df,
        x='logfoldchanges', y='-log10(padj)',
        hue='sig', palette=palette, s=15, linewidth=0
    )

    # threshold lines
    plt.axhline(-np.log10(padj_thr), color='black', linestyle='--', lw=0.8)
    plt.axvline(lfc_thr, color='black', linestyle='--', lw=0.8)
    plt.axvline(-lfc_thr, color='black', linestyle='--', lw=0.8)

    # label top genes (optional adjustText)
    texts = []
    for _, row in top_genes.iterrows():
        texts.append(plt.text(row['logfoldchanges'], row['-log10(padj)'], row['names'], fontsize=8))
    try:
        from adjustText import adjust_text
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
    except Exception:
        pass

    plt.xlabel('log$_2$ Fold Change', fontsize=12)
    plt.ylabel('-log$_{10}$ Adjusted P-value', fontsize=12)
    if title:
        plt.title(title, fontsize=13, pad=10)
    plt.legend(frameon=False, title='', loc='upper right')
    sns.despine()
    plt.tight_layout()

    if save_path is not None:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    plt.close()


def plot_enrichment_barplot(pos_enr_df, neg_enr_df=None, top_n=10, save_path=None, title="Enrichr", combined_score=100, gene_sets='KEGG_2021_Human', figsize=(10, 6), show=False):
    def prep(df, direction):
        if df is None or df.empty:
            return pd.DataFrame(columns=["clean_term", "log10_adjusted_p_value", "direction"])
        d = df.nsmallest(top_n, "Adjusted P-value").copy()
        # Clean term based on gene set type
        if "GO" in gene_sets:
            # Remove GO IDs like (GO:0008150)
            d["clean_term"] = d["Term"].astype(str).str.replace(r"\s*\(GO:\d+\)", "", regex=True)
        elif "Reactome" in gene_sets or "REACTOME" in gene_sets:
            # Remove REACTOME pathway IDs like (R-HSA-xxxxx)
            d["clean_term"] = d["Term"].astype(str).str.replace(r"\s*\(R-[A-Z]+-\d+\)", "", regex=True)
        else:
            # For other gene sets (Hallmark, etc.), use term as-is
            d["clean_term"] = d["Term"].astype(str)
        d["log10_adjusted_p_value"] = -np.log10(d["Adjusted P-value"].astype(float) + 1e-300)
        d["direction"] = direction
        if direction == "Down":
            d["log10_adjusted_p_value"] = -d["log10_adjusted_p_value"]  # Down to the left
        return d[["clean_term", "log10_adjusted_p_value", "direction"]]

    # Filter by combined score
    pos_enr_df = pos_enr_df[pos_enr_df["Combined Score"] > combined_score]
    if neg_enr_df is not None:
        neg_enr_df = neg_enr_df[neg_enr_df["Combined Score"] > combined_score]
    
    up = prep(pos_enr_df, "Up")
    down = prep(neg_enr_df, "Down") if neg_enr_df is not None else pd.DataFrame(columns=["clean_term", "log10_adjusted_p_value", "direction"])
    
    # Fix FutureWarning: only concatenate non-empty dataframes
    dfs_to_concat = [df for df in [up, down] if not df.empty]
    if not dfs_to_concat:
        return
    df = pd.concat(dfs_to_concat, ignore_index=True)

    fig, ax = plt.subplots(figsize=figsize)
    colors = df["direction"].map({"Up": "#c44b41", "Down": "#86a6c7"}).values

    ax.barh(range(len(df)), df["log10_adjusted_p_value"], color=colors, alpha=0.7)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["clean_term"])
    ax.set_xlabel("-log10(Adjusted P-value)")
    ax.set_title(title)
    ax.invert_yaxis()
    ax.axvline(0, color="black", lw=0.8)
    # separate x-limits based on up/down magnitudes
    up_max = up["log10_adjusted_p_value"].max() if not up.empty else 0
    down_max = (-down["log10_adjusted_p_value"]).max() if not down.empty else 0  # convert back to positive magnitude
    ax.set_xlim(-(down_max * 1.1 if down_max > 0 else 1), (up_max * 1.1 if up_max > 0 else 1))

    # Fix UserWarning: set ticks before labels
    ticks = ax.get_xticks()
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{abs(t):g}" for t in ticks])

    plt.subplots_adjust(left=0.4)

    if save_path is not None:
        dir_path = os.path.dirname(save_path)
        if dir_path:
            os.makedirs(dir_path, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)