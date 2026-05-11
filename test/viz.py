"""
test/viz.py - Step 5: Visualize synteny network (from cluster output)

Publication-quality visualization inspired by R ggraph / visNetwork / Cytoscape.

Features:
  - ggraph-style static network with curved edges (FancyArrowPatch),
    node sizing by degree, publication-ready fonts and legend
  - Cytoscape-style cluster tiling: each cluster rendered as a separate
    subgraph arranged in a grid layout with borders and backgrounds
  - Interactive HTML with adjustable physics, visual controls, and PNG export
  - Cluster overview: size distribution, species composition barplot

Usage:
    python viz.py -c Filtered.clusters.tsv -n Final_Network.tsv -s species.lst -d seqs
    python viz.py -c Filtered.clusters.tsv -n Final_Network.tsv -s species.lst --palette palette.tsv
    python viz.py -c Filtered.clusters.tsv -n Final_Network.tsv -s species.lst --plot-type static --layout circular
    python viz.py -c Filtered.clusters.tsv -n Final_Network.tsv -s species.lst --plot-type interactive
    python viz.py -c Filtered.clusters.tsv -n Final_Network.tsv -s species.lst --plot-type tiled --top-clusters 20
    python viz.py -c Filtered.clusters.tsv -n Final_Network.tsv -s species.lst --top-k 100
"""

import sys
import json
import argparse
from pathlib import Path
from typing import List, Dict, Set, Optional, Tuple
from collections import defaultdict
import math

try:
    import networkx as nx
    _HAS_NX = True
except ImportError:
    _HAS_NX = False

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.patches import FancyArrowPatch
    from matplotlib.lines import Line2D
    import matplotlib.colors as mcolors
    import matplotlib.path as mpath
    import numpy as np
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False

try:
    from pyvis.network import Network as PyVisNetwork
    _HAS_PYVIS = True
except ImportError:
    _HAS_PYVIS = False


SPECIES_PALETTE = [
    '#e6194b', '#3cb44b', '#4363d8', '#f58231', '#911eb4',
    '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990',
    '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3',
    '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#000000',
]

PLOT_STYLE = {
    'font_family': 'sans-serif',
    'font_family_serif': 'serif',
    'title_size': 14,
    'label_size': 10,
    'tick_size': 8,
    'legend_size': 9,
    'edge_color_light': '#b0b0b0',
    'edge_color_dark': '#505050',
    'bg_color': '#ffffff',
    'grid_color': '#f0f0f0',
    'cluster_bg': '#fafafa',
    'cluster_border': '#d0d0d0',
}


def log_info(msg):
    print(f"[INFO] {msg}", file=sys.stderr)


def log_warn(msg):
    print(f"[WARN] {msg}", file=sys.stderr)


def log_error(msg):
    print(f"[ERROR] {msg}", file=sys.stderr)


def load_species_list(list_file: str) -> List[str]:
    with open(list_file, 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]


def build_gene_species_map(species_list: List[str], work_dir: str = ".") -> Dict[str, str]:
    gene_map = {}
    wd = Path(work_dir)
    for sp in species_list:
        bed_file = wd / f"{sp}.bed"
        if not bed_file.exists():
            continue
        with open(bed_file, 'r') as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 4:
                    gene_map[parts[3]] = sp
    return gene_map


def infer_species(gene_id: str, gene_species_map: Dict[str, str]) -> Optional[str]:
    if gene_id in gene_species_map:
        return gene_species_map[gene_id]
    base = gene_id.split('.')[0]
    if base in gene_species_map:
        return gene_species_map[base]
    return None


def load_palette(palette_file: str, species_list: List[str]) -> Dict[str, str]:
    color_map = {}
    with open(palette_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                color_map[parts[0].strip()] = parts[1].strip()
    for sp in species_list:
        if sp not in color_map:
            log_warn(f"Species '{sp}' not in palette file, using default")
    return color_map


def build_species_color_map(species_list: List[str], palette_file: Optional[str] = None) -> Dict[str, str]:
    if palette_file:
        custom = load_palette(palette_file, species_list)
        color_map = {}
        for i, sp in enumerate(species_list):
            color_map[sp] = custom.get(sp, SPECIES_PALETTE[i % len(SPECIES_PALETTE)])
        log_info(f"Using custom palette from {palette_file} ({len(custom)} species defined)")
        return color_map
    return {sp: SPECIES_PALETTE[i % len(SPECIES_PALETTE)] for i, sp in enumerate(species_list)}


def load_cluster_tsv(cluster_tsv: str) -> Tuple[Dict[str, str], Dict[str, int], Dict[str, Set[str]]]:
    gene_to_cluster = {}
    cluster_sizes = {}
    cluster_genes = defaultdict(set)
    with open(cluster_tsv, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            gene_id, cluster_id = parts[0], parts[1]
            gene_to_cluster[gene_id] = cluster_id
            cluster_genes[cluster_id].add(gene_id)
    for cid, genes in cluster_genes.items():
        cluster_sizes[cid] = len(genes)
    return gene_to_cluster, cluster_sizes, dict(cluster_genes)


def load_network_tsv(network_tsv: str) -> Tuple[List[Tuple[str, str, float]], Set[str]]:
    edges = []
    nodes = set()
    with open(network_tsv, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            src, tgt = parts[0], parts[1]
            try:
                score = float(parts[2])
            except ValueError:
                score = 0.0
            edges.append((src, tgt, score))
            nodes.add(src)
            nodes.add(tgt)
    return edges, nodes


def build_nx_graph(edges, nodes, gene_species_map=None, gene_cluster_map=None):
    G = nx.Graph()
    for node in nodes:
        attrs = {}
        if gene_species_map:
            attrs['species'] = infer_species(node, gene_species_map) or "unknown"
        if gene_cluster_map and node in gene_cluster_map:
            attrs['cluster'] = gene_cluster_map[node]
        G.add_node(node, **attrs)
    for src, tgt, score in edges:
        G.add_edge(src, tgt, weight=score)
    return G


def _setup_pub_style():
    plt.rcParams.update({
        'font.family': PLOT_STYLE['font_family'],
        'font.size': PLOT_STYLE['label_size'],
        'axes.titlesize': PLOT_STYLE['title_size'],
        'axes.labelsize': PLOT_STYLE['label_size'],
        'xtick.labelsize': PLOT_STYLE['tick_size'],
        'ytick.labelsize': PLOT_STYLE['tick_size'],
        'legend.fontsize': PLOT_STYLE['legend_size'],
        'figure.facecolor': PLOT_STYLE['bg_color'],
        'axes.facecolor': PLOT_STYLE['bg_color'],
        'savefig.facecolor': PLOT_STYLE['bg_color'],
        'savefig.dpi': 300,
        'axes.spines.top': False,
        'axes.spines.right': False,
    })


def _compute_layout(G, layout="force", seed=42):
    n_nodes = G.number_of_nodes()
    if layout == "force":
        k_val = 2.0 / (n_nodes ** 0.5 + 1)
        pos = nx.spring_layout(G, k=k_val, iterations=120, seed=seed)
    elif layout == "circular":
        pos = nx.circular_layout(G)
    elif layout == "kamada_kawai":
        try:
            pos = nx.kamada_kawai_layout(G)
        except Exception:
            pos = nx.spring_layout(G, seed=seed)
    else:
        pos = nx.spring_layout(G, seed=seed)
    return pos


def _draw_curved_edges(ax, G, pos, edge_widths, edge_alphas, edge_color=None):
    if edge_color is None:
        edge_color = PLOT_STYLE['edge_color_light']

    edge_rad_map = {}
    for (src, tgt), width, alpha in zip(G.edges(), edge_widths, edge_alphas):
        key = tuple(sorted([src, tgt]))
        edge_rad_map[key] = edge_rad_map.get(key, 0) + 1

    edge_count = {}
    for (src, tgt), width, alpha in zip(G.edges(), edge_widths, edge_alphas):
        x0, y0 = pos[src]
        x1, y1 = pos[tgt]
        key = tuple(sorted([src, tgt]))
        idx = edge_count.get(key, 0)
        edge_count[key] = idx + 1

        total = edge_rad_map[key]
        if total <= 1:
            rad = 0.0
        else:
            rad = 0.15 * (idx - (total - 1) / 2.0)

        if abs(rad) < 0.01:
            ax.plot([x0, x1], [y0, y1],
                    color=edge_color, linewidth=width, alpha=alpha, zorder=1)
        else:
            arrow = FancyArrowPatch(
                (x0, y0), (x1, y1),
                connectionstyle=f"arc3,rad={rad}",
                arrowstyle="-",
                color=edge_color,
                linewidth=width,
                alpha=alpha,
                zorder=1,
            )
            ax.add_patch(arrow)


def plot_ggraph_network(
        edges, nodes, output_file,
        gene_species_map=None, species_list=None, species_color_map=None,
        layout="force", top_k=None,
        figsize=(14, 10), dpi=300,
        title=None,
):
    if not _HAS_MPL or not _HAS_NX:
        raise ImportError("matplotlib + networkx required")

    _setup_pub_style()

    G = build_nx_graph(edges, nodes, gene_species_map)

    if top_k and top_k > 0:
        degree_dict = dict(G.degree())
        top_nodes = sorted(degree_dict, key=degree_dict.get, reverse=True)[:top_k]
        G = G.subgraph(top_nodes).copy()
        log_info(f"Showing top {top_k} nodes by degree")

    n_nodes = G.number_of_nodes()
    pos = _compute_layout(G, layout)

    if species_color_map is None:
        species_color_map = {}
        if species_list:
            for i, sp in enumerate(species_list):
                species_color_map[sp] = SPECIES_PALETTE[i % len(SPECIES_PALETTE)]

    node_colors = []
    for node in G.nodes():
        sp = G.nodes[node].get('species', 'unknown')
        node_colors.append(species_color_map.get(sp, '#cccccc'))

    degrees = dict(G.degree())
    max_deg = max(degrees.values()) if degrees else 1
    min_deg = min(degrees.values()) if degrees else 1
    deg_range = max(max_deg - min_deg, 1)
    node_sizes = [30 + 270 * ((degrees[n] - min_deg) / deg_range) for n in G.nodes()]

    edge_weights = [G.edges[e].get('weight', 1.0) for e in G.edges()]
    max_w = max(edge_weights) if edge_weights else 1
    min_w = min(edge_weights) if edge_weights else 0
    w_range = max(max_w - min_w, 1)
    edge_widths = [0.3 + 2.0 * ((w - min_w) / w_range) for w in edge_weights]
    edge_alphas = [0.08 + 0.25 * ((w - min_w) / w_range) for w in edge_weights]

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    _draw_curved_edges(ax, G, pos, edge_widths, edge_alphas)

    xs = [pos[n][0] for n in G.nodes()]
    ys = [pos[n][1] for n in G.nodes()]
    ax.scatter(xs, ys, c=node_colors, s=node_sizes,
               alpha=0.85, edgecolors='white', linewidths=0.5, zorder=2)

    if n_nodes <= 80:
        for node in G.nodes():
            x, y = pos[node]
            ax.annotate(node, (x, y), fontsize=4, alpha=0.6,
                        ha='center', va='bottom', xytext=(0, 3),
                        textcoords='offset points', zorder=3)

    if species_list and species_color_map:
        legend_handles = [
            mpatches.Patch(facecolor=species_color_map[sp], edgecolor='white', label=sp)
            for sp in species_list if sp in species_color_map
        ]
        deg_sizes = [min_deg, (min_deg + max_deg) // 2, max_deg] if max_deg > min_deg else [min_deg]
        for d in deg_sizes:
            sz = 30 + 270 * ((d - min_deg) / deg_range)
            legend_handles.append(
                Line2D([0], [0], marker='o', color='w', markerfacecolor='#cccccc',
                       markersize=max(3, min(sz ** 0.5, 12)), label=f'deg={d}')
            )
        ax.legend(handles=legend_handles, loc='upper left', frameon=True,
                  framealpha=0.9, edgecolor='#cccccc', fancybox=True)

    if title:
        ax.set_title(title, fontweight='bold', pad=12)
    else:
        ax.set_title(f"Synteny Network ({n_nodes} nodes, {G.number_of_edges()} edges)",
                     fontweight='bold', pad=12)

    ax.text(0.99, 0.01,
            f"Nodes: {n_nodes}  Edges: {G.number_of_edges()}\n"
            f"Layout: {layout}  Degree: {min_deg}-{max_deg}",
            transform=ax.transAxes, fontsize=7, color='#888888',
            ha='right', va='bottom', family='monospace')

    ax.axis('off')
    plt.tight_layout()

    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_path), dpi=dpi, bbox_inches='tight')
    plt.close()
    log_info(f"ggraph-style plot saved: {out_path}")
    return str(out_path)


def plot_tiled_clusters(
        edges, nodes, gene_cluster_map, cluster_genes, cluster_sizes,
        output_file,
        gene_species_map=None, species_list=None, species_color_map=None,
        top_clusters=20, layout="force",
        figsize=None, dpi=300,
        title=None,
):
    if not _HAS_MPL or not _HAS_NX:
        raise ImportError("matplotlib + networkx required")

    _setup_pub_style()

    if species_color_map is None:
        species_color_map = {}
        if species_list:
            for i, sp in enumerate(species_list):
                species_color_map[sp] = SPECIES_PALETTE[i % len(SPECIES_PALETTE)]

    sorted_clusters = sorted(cluster_genes.items(), key=lambda x: len(x[1]), reverse=True)
    selected = sorted_clusters[:top_clusters]

    n_clusters = len(selected)
    n_cols = min(5, n_clusters)
    n_rows = math.ceil(n_clusters / n_cols)

    if figsize is None:
        fw = max(4 * n_cols, 12)
        fh = max(4 * n_rows, 8)
        figsize = (fw, fh)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_clusters == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)

    G_full = build_nx_graph(edges, nodes, gene_species_map, gene_cluster_map)

    edge_set = set()
    for src, tgt, _ in edges:
        edge_set.add((src, tgt))
        edge_set.add((tgt, src))

    for idx, (cid, genes) in enumerate(selected):
        row, col = divmod(idx, n_cols)
        ax = axes[row][col]

        sub_nodes = set(genes)
        sub_edges = [(s, t) for s, t in edge_set if s in sub_nodes and t in sub_nodes]

        G_sub = nx.Graph()
        G_sub.add_nodes_from(sub_nodes)
        for s, t in sub_edges:
            G_sub.add_edge(s, t)

        if G_sub.number_of_nodes() == 0:
            ax.axis('off')
            continue

        if layout == "force":
            k_val = 1.5 / (G_sub.number_of_nodes() ** 0.5 + 1)
            pos = nx.spring_layout(G_sub, k=k_val, iterations=50, seed=42)
        elif layout == "circular":
            pos = nx.circular_layout(G_sub)
        else:
            pos = nx.spring_layout(G_sub, seed=42)

        node_colors = []
        for node in G_sub.nodes():
            sp = G_full.nodes[node].get('species', 'unknown') if node in G_full else 'unknown'
            node_colors.append(species_color_map.get(sp, '#cccccc'))

        ax.set_facecolor(PLOT_STYLE['cluster_bg'])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color(PLOT_STYLE['cluster_border'])
            spine.set_linewidth(0.8)

        sub_edge_weights = []
        for s, t in G_sub.edges():
            w = G_full.edges[s, t].get('weight', 1.0) if G_full.has_edge(s, t) else 1.0
            sub_edge_weights.append(w)

        if sub_edge_weights:
            max_sw = max(sub_edge_weights)
            min_sw = min(sub_edge_weights)
            sw_range = max(max_sw - min_sw, 1)
            sub_widths = [0.3 + 1.5 * ((w - min_sw) / sw_range) for w in sub_edge_weights]
        else:
            sub_widths = []

        for (s, t), w in zip(G_sub.edges(), sub_widths):
            x0, y0 = pos[s]
            x1, y1 = pos[t]
            rad = 0.1 if G_sub.number_of_edges() < 30 else 0.05
            if abs(rad) < 0.01:
                ax.plot([x0, x1], [y0, y1], alpha=0.3,
                        color=PLOT_STYLE['edge_color_light'], linewidth=w, zorder=1)
            else:
                arrow = FancyArrowPatch(
                    (x0, y0), (x1, y1),
                    connectionstyle=f"arc3,rad={rad}",
                    arrowstyle="-",
                    color=PLOT_STYLE['edge_color_light'],
                    linewidth=w, alpha=0.3, zorder=1,
                )
                ax.add_patch(arrow)

        xs = [pos[n][0] for n in G_sub.nodes()]
        ys = [pos[n][1] for n in G_sub.nodes()]
        ax.scatter(xs, ys, c=node_colors, s=20,
                   edgecolors='white', linewidths=0.2, alpha=0.85, zorder=2)

        species_in_cluster = set()
        for gene in genes:
            sp = infer_species(gene, gene_species_map) if gene_species_map else None
            if sp:
                species_in_cluster.add(sp)

        ax.set_title(f"{cid} (n={len(genes)}, sp={len(species_in_cluster)})",
                     fontsize=8, fontweight='bold', pad=3)
        ax.set_xticks([])
        ax.set_yticks([])

    for idx in range(n_clusters, n_rows * n_cols):
        row, col = divmod(idx, n_cols)
        axes[row][col].axis('off')

    if species_list and species_color_map:
        legend_handles = [
            mpatches.Patch(facecolor=species_color_map[sp], edgecolor='white', label=sp)
            for sp in species_list if sp in species_color_map
        ]
        fig.legend(handles=legend_handles, loc='lower center', ncol=len(species_list),
                   frameon=True, framealpha=0.9, edgecolor='#cccccc', fancybox=True,
                   fontsize=PLOT_STYLE['legend_size'])

    if title:
        fig.suptitle(title, fontweight='bold', fontsize=PLOT_STYLE['title_size'], y=0.98)
    else:
        fig.suptitle(f"Cluster Overview (top {n_clusters} of {len(cluster_genes)})",
                     fontweight='bold', fontsize=PLOT_STYLE['title_size'], y=0.98)

    plt.tight_layout(rect=[0, 0.04, 1, 0.96])

    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_path), dpi=dpi, bbox_inches='tight')
    plt.close()
    log_info(f"Tiled cluster plot saved: {out_path}")
    return str(out_path)


_INTERACTIVE_HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>__TITLE__</title>
<script type="text/javascript" src="https://unpkg.com/vis-network@9.1.6/standalone/umd/vis-network.min.js"></script>
<style>
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:'Segoe UI',system-ui,sans-serif;display:flex;height:100vh;background:#f8f8f8}
#sidebar{width:280px;background:#fff;padding:16px;overflow-y:auto;border-right:1px solid #e0e0e0;flex-shrink:0}
#main{flex:1;position:relative;display:flex;flex-direction:column}
#toolbar{height:40px;background:#fff;border-bottom:1px solid #e0e0e0;display:flex;align-items:center;padding:0 12px;gap:8px}
#toolbar span{font-size:12px;color:#666}
#network-wrap{flex:1;position:relative;background:#fff}
#network-canvas{width:100%;height:100%}
.cg{margin-bottom:14px;padding-bottom:10px;border-bottom:1px solid #f0f0f0}
.cg h3{font-size:12px;margin-bottom:6px;color:#333;text-transform:uppercase;letter-spacing:0.5px}
.cg label{display:block;font-size:11px;margin-bottom:2px;color:#555}
.cg input[type=range]{width:100%;height:4px;-webkit-appearance:none;background:#e0e0e0;border-radius:2px;outline:none}
.cg input[type=range]::-webkit-slider-thumb{-webkit-appearance:none;width:14px;height:14px;border-radius:50%;background:#4363d8;cursor:pointer}
.cg select{width:100%;padding:4px 6px;font-size:11px;border:1px solid #ddd;border-radius:3px}
.cg input[type=checkbox]{margin-right:4px}
.btn{padding:6px 14px;border:none;border-radius:4px;cursor:pointer;font-size:11px;font-weight:500}
.btn-p{background:#4363d8;color:#fff}.btn-p:hover{background:#3453c0}
.btn-s{background:#e8e8e8;color:#333}.btn-s:hover{background:#d0d0d0}
.stats{font-size:11px;color:#666;padding:8px;background:#f5f5f5;border-radius:4px;margin-bottom:10px;line-height:1.6}
.legend{margin-top:10px}
.legend-item{display:flex;align-items:center;margin-bottom:3px;font-size:11px}
.legend-dot{width:10px;height:10px;border-radius:50%;margin-right:6px;flex-shrink:0}
#toast{position:fixed;bottom:20px;right:20px;background:#333;color:#fff;padding:10px 20px;border-radius:6px;font-size:12px;display:none;z-index:999;transition:opacity 0.3s}
.val{color:#4363d8;font-weight:600}
</style>
</head>
<body>
<div id="sidebar">
<h2 style="font-size:15px;margin-bottom:4px">SynNet Viz</h2>
<p style="font-size:10px;color:#999;margin-bottom:12px">Interactive Network Explorer</p>
<div class="stats" id="stats">Loading...</div>
<div class="cg">
<h3>Physics</h3>
<label>Gravity: <span class="val" id="v-grav">-8000</span></label>
<input type="range" id="r-grav" min="-20000" max="-500" value="-8000" step="500">
<label>Spring Length: <span class="val" id="v-slen">120</span></label>
<input type="range" id="r-slen" min="20" max="400" value="120" step="10">
<label>Spring Constant: <span class="val" id="v-scon">0.04</span></label>
<input type="range" id="r-scon" min="0.005" max="0.2" value="0.04" step="0.005">
<label>Damping: <span class="val" id="v-damp">0.09</span></label>
<input type="range" id="r-damp" min="0.01" max="0.5" value="0.09" step="0.01">
</div>
<div class="cg">
<h3>Visual</h3>
<label>Node Scale: <span class="val" id="v-nscale">1.0</span>x</label>
<input type="range" id="r-nscale" min="0.3" max="3.0" value="1.0" step="0.1">
<label>Edge Scale: <span class="val" id="v-escale">1.0</span>x</label>
<input type="range" id="r-escale" min="0.2" max="3.0" value="1.0" step="0.1">
<label><input type="checkbox" id="cb-label"> Show Labels</label>
<label><input type="checkbox" id="cb-curve" checked> Curved Edges</label>
</div>
<div class="cg">
<h3>Filter</h3>
<label>Min Degree: <span class="val" id="v-mindeg">0</span></label>
<input type="range" id="r-mindeg" min="0" max="__MAX_DEG__" value="0" step="1">
</div>
<div class="cg">
<h3>Export</h3>
<button class="btn btn-p" onclick="exportPNG()">&#128247; Save PNG</button>
<button class="btn btn-s" onclick="network.stabilize()">&#9654; Stabilize</button>
<button class="btn btn-s" onclick="network.stopSimulation()">&#9724; Stop</button>
<button class="btn btn-s" onclick="resetPhysics()">&#8635; Reset</button>
</div>
<div class="legend" id="legend"></div>
</div>
<div id="main">
<div id="toolbar">
<span id="tb-info"></span>
</div>
<div id="network-wrap">
<div id="network-canvas"></div>
</div>
</div>
<div id="toast"></div>
<script>
var nodesData=new vis.DataSet(__NODES_JSON__);
var edgesData=new vis.DataSet(__EDGES_JSON__);
var allNodes=nodesData.get();
var allEdges=edgesData.get();
var container=document.getElementById('network-canvas');
var data={nodes:nodesData,edges:edgesData};
var options={
  physics:{
    barnesHut:{gravitationalConstant:-8000,centralGravity:0.3,springLength:120,springConstant:0.04,damping:0.09},
    stabilization:{iterations:200}
  },
  nodes:{font:{size:0,color:'#333'},borderWidth:0.5,borderWidthSelected:2,shadow:{enabled:true,color:'rgba(0,0,0,0.1)',size:3}},
  edges:{smooth:{type:'continuous'},color:{color:'#c0c0c0',highlight:'#4363d8',hover:'#888'},shadow:false},
  interaction:{hover:true,tooltipDelay:200,navigationButtons:true,keyboard:true,multiselect:true}
};
var network=new vis.Network(container,data,options);
var baseNodeSizes={};
nodesData.forEach(function(n){baseNodeSizes[n.id]=n.size||10;});
var baseEdgeWidths={};
edgesData.forEach(function(e){baseEdgeWidths[e.id]=e.width||1;});

function updateStats(){
  var nc=nodesData.length,ec=edgesData.length;
  document.getElementById('stats').innerHTML='<b>Nodes:</b> '+nc+' &nbsp; <b>Edges:</b> '+ec;
  document.getElementById('tb-info').textContent='Nodes: '+nc+' | Edges: '+ec;
}
network.on('stabilized',updateStats);
updateStats();

function bindSlider(rid,vid,cb){
  var r=document.getElementById(rid),v=document.getElementById(vid);
  r.oninput=function(){v.textContent=this.value;cb(this.value);};
}
bindSlider('r-grav','v-grav',function(v){network.setOptions({physics:{barnesHut:{gravitationalConstant:+v}}});});
bindSlider('r-slen','v-slen',function(v){network.setOptions({physics:{barnesHut:{springLength:+v}}});});
bindSlider('r-scon','v-scon',function(v){network.setOptions({physics:{barnesHut:{springConstant:+v}}});});
bindSlider('r-damp','v-damp',function(v){network.setOptions({physics:{barnesHut:{damping:+v}}});});

bindSlider('r-nscale','v-nscale',function(v){
  var s=parseFloat(v);
  nodesData.forEach(function(n){nodesData.update({id:n.id,size:baseNodeSizes[n.id]*s});});
});
bindSlider('r-escale','v-escale',function(v){
  var s=parseFloat(v);
  edgesData.forEach(function(e){edgesData.update({id:e.id,width:baseEdgeWidths[e.id]*s});});
});

document.getElementById('cb-label').onchange=function(){
  network.setOptions({nodes:{font:{size:this.checked?10:0}}});
};
document.getElementById('cb-curve').onchange=function(){
  network.setOptions({edges:{smooth:this.checked?{type:'continuous'}:false}});
};

var degreeMap={};
allEdges.forEach(function(e){
  degreeMap[e.from]=(degreeMap[e.from]||0)+1;
  degreeMap[e.to]=(degreeMap[e.to]||0)+1;
});
bindSlider('r-mindeg','v-mindeg',function(v){
  var minD=parseInt(v);
  var keep=new Set();
  allNodes.forEach(function(n){if((degreeMap[n.id]||0)>=minD)keep.add(n.id);});
  var filteredNodes=allNodes.filter(function(n){return keep.has(n.id);});
  var filteredEdges=allEdges.filter(function(e){return keep.has(e.from)&&keep.has(e.to);});
  nodesData.clear();edgesData.clear();
  nodesData.add(filteredNodes);edgesData.add(filteredEdges);
  updateStats();
});

function exportPNG(){
  var canvas=container.querySelector('canvas');
  if(canvas){
    var link=document.createElement('a');
    link.download='synnet_network.png';
    link.href=canvas.toDataURL('image/png');
    link.click();
    toast('PNG saved!');
  }else{
    toast('Canvas not ready');
  }
}
function resetPhysics(){
  document.getElementById('r-grav').value=-8000;document.getElementById('v-grav').textContent='-8000';
  document.getElementById('r-slen').value=120;document.getElementById('v-slen').textContent='120';
  document.getElementById('r-scon').value=0.04;document.getElementById('v-scon').textContent='0.04';
  document.getElementById('r-damp').value=0.09;document.getElementById('v-damp').textContent='0.09';
  network.setOptions({physics:{barnesHut:{gravitationalConstant:-8000,springLength:120,springConstant:0.04,damping:0.09}}});
  network.stabilize();
  toast('Physics reset');
}
function toast(msg){
  var t=document.getElementById('toast');t.textContent=msg;t.style.display='block';
  setTimeout(function(){t.style.display='none';},2000);
}

var spList=__SPECIES_LIST__;
var spColors=__SPECIES_COLORS__;
var ld=document.getElementById('legend');
var lh='<h3 style="font-size:11px;margin-bottom:6px">Species</h3>';
spList.forEach(function(sp){
  lh+='<div class="legend-item"><div class="legend-dot" style="background:'+(spColors[sp]||'#ccc')+'"></div>'+sp+'</div>';
});
ld.innerHTML=lh;
</script>
</body>
</html>"""


def plot_interactive_app(
        edges, nodes, output_file,
        gene_species_map=None, gene_cluster_map=None,
        species_list=None, species_color_map=None,
        cluster_genes=None, cluster_sizes=None,
        top_k=None, title=None,
):
    if not _HAS_NX:
        raise ImportError("networkx required for interactive plots")

    G = build_nx_graph(edges, nodes, gene_species_map, gene_cluster_map)

    if top_k and top_k > 0:
        degree_dict = dict(G.degree())
        top_nodes = sorted(degree_dict, key=degree_dict.get, reverse=True)[:top_k]
        G = G.subgraph(top_nodes).copy()
        log_info(f"Showing top {top_k} nodes by degree")

    if species_color_map is None:
        species_color_map = {}
        if species_list:
            for i, sp in enumerate(species_list):
                species_color_map[sp] = SPECIES_PALETTE[i % len(SPECIES_PALETTE)]

    degrees = dict(G.degree())
    max_deg = max(degrees.values()) if degrees else 1

    vis_nodes = []
    for node in G.nodes():
        sp = G.nodes[node].get('species', 'unknown')
        color = species_color_map.get(sp, '#cccccc')
        cluster = G.nodes[node].get('cluster', '')
        deg = degrees.get(node, 0)
        size = 6 + 18 * (deg / max(max_deg, 1))
        title_str = f"{node}\nSpecies: {sp}\nDegree: {deg}"
        if cluster:
            title_str += f"\nCluster: {cluster}"
        vis_nodes.append({
            "id": node, "label": node, "color": color,
            "title": title_str, "size": round(size, 1),
            "borderWidth": 0.5, "borderWidthSelected": 2,
        })

    edge_weights = [G.edges[e].get('weight', 1.0) for e in G.edges()]
    max_w = max(edge_weights) if edge_weights else 1
    min_w = min(edge_weights) if edge_weights else 0
    w_range = max(max_w - min_w, 1)

    vis_edges = []
    for i, (src, tgt) in enumerate(G.edges()):
        w = G.edges[src, tgt].get('weight', 1.0)
        width = 0.5 + 2.5 * ((w - min_w) / w_range)
        vis_edges.append({
            "from": src, "to": tgt,
            "value": round(w, 1), "width": round(width, 2),
            "title": f"score: {w:.0f}",
        })

    html = _INTERACTIVE_HTML_TEMPLATE
    html = html.replace('__TITLE__', title or 'SynNet Interactive Visualization')
    html = html.replace('__NODES_JSON__', json.dumps(vis_nodes))
    html = html.replace('__EDGES_JSON__', json.dumps(vis_edges))
    html = html.replace('__SPECIES_LIST__', json.dumps(species_list or []))
    html = html.replace('__SPECIES_COLORS__', json.dumps(species_color_map or {}))
    html = html.replace('__MAX_DEG__', str(max(degrees.values()) if degrees else 10))

    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(html)

    log_info(f"Interactive app saved: {out_path}")
    log_info(f"Open in browser to adjust layout and export PNG")
    return str(out_path)


def plot_cluster_size_dist(cluster_sizes, output_file, figsize=(8, 5), dpi=300, bins=50, title=None):
    if not _HAS_MPL:
        raise ImportError("matplotlib required")

    _setup_pub_style()
    sizes = list(cluster_sizes.values())

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.hist(sizes, bins=bins, color='#4363d8', edgecolor='white', alpha=0.85, linewidth=0.5)
    ax.set_xlabel("Cluster size (number of genes)")
    ax.set_ylabel("Count")
    ax.set_title(title or f"Cluster size distribution (n={len(sizes)})", fontweight='bold')

    if max(sizes) > 10 * (sum(sizes) / len(sizes)):
        ax.set_yscale('log')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_path), dpi=dpi, bbox_inches='tight')
    plt.close()
    log_info(f"Cluster distribution plot saved: {out_path}")
    return str(out_path)


def plot_species_composition(
        cluster_genes, species_list, gene_species_map, output_file,
        species_color_map=None, figsize=(12, 6), dpi=300, top_n=30, title=None,
):
    if not _HAS_MPL:
        raise ImportError("matplotlib required")

    _setup_pub_style()

    if species_color_map is None:
        species_color_map = {sp: SPECIES_PALETTE[i % len(SPECIES_PALETTE)]
                             for i, sp in enumerate(species_list)}

    sorted_clusters = sorted(cluster_genes.items(), key=lambda x: len(x[1]), reverse=True)[:top_n]

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    cluster_labels = []
    bottom = np.zeros(len(sorted_clusters))

    for sp in species_list:
        counts = []
        for cid, genes in sorted_clusters:
            sp_count = sum(1 for gene in genes if infer_species(gene, gene_species_map) == sp)
            counts.append(sp_count)
        cluster_labels = [cid for cid, _ in sorted_clusters]

        ax.bar(cluster_labels, counts, bottom=bottom,
               color=species_color_map[sp], label=sp, edgecolor='white', linewidth=0.3)
        bottom = bottom + np.array(counts)

    ax.set_xlabel("Cluster ID")
    ax.set_ylabel("Number of genes")
    ax.set_title(title or f"Species composition (top {top_n} clusters)", fontweight='bold')
    ax.legend(frameon=True, framealpha=0.9, edgecolor='#cccccc')
    plt.xticks(rotation=45, ha='right', fontsize=7)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    out_path = Path(output_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(out_path), dpi=dpi, bbox_inches='tight')
    plt.close()
    log_info(f"Species composition plot saved: {out_path}")
    return str(out_path)


def run_viz(
        cluster_tsv, network_tsv, species_list_file,
        bed_dir=".", palette_file=None,
        output_dir=".", plot_type="all",
        layout="force", top_k=None, top_clusters=20,
        interactive=False,
        figsize=None, dpi=300,
):
    log_info("Step 5: SynNet Visualization")

    species_list = load_species_list(species_list_file)
    log_info(f"Loaded {len(species_list)} species: {' -> '.join(species_list)}")

    gene_species_map = build_gene_species_map(species_list, bed_dir)
    if gene_species_map:
        log_info(f"Built gene-species map: {len(gene_species_map)} genes")
    else:
        log_warn("No .bed files found, species coloring will be limited")

    species_color_map = build_species_color_map(species_list, palette_file)

    gene_cluster_map, cluster_sizes, cluster_genes = load_cluster_tsv(cluster_tsv)
    log_info(f"Loaded {len(cluster_sizes)} clusters from {cluster_tsv}")

    edges, nodes = load_network_tsv(network_tsv)
    log_info(f"Loaded {len(edges)} edges, {len(nodes)} nodes from {network_tsv}")

    filtered_nodes = {n for n in nodes if n in gene_cluster_map}
    filtered_edges = [(s, t, w) for s, t, w in edges if s in filtered_nodes and t in filtered_nodes]
    log_info(f"Filtered to {len(filtered_edges)} edges within clustered nodes")

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / Path(cluster_tsv).stem.replace(".clusters", "")

    results = {}

    if plot_type in ("static", "all"):
        if _HAS_MPL and _HAS_NX:
            path = plot_ggraph_network(
                filtered_edges, filtered_nodes,
                str(out_dir / f"{prefix.name}_network.png"),
                gene_species_map=gene_species_map,
                species_list=species_list,
                species_color_map=species_color_map,
                layout=layout, top_k=top_k,
                figsize=figsize or (14, 10), dpi=dpi,
            )
            results["static"] = path
        else:
            log_warn("matplotlib/networkx not available, skip static plot")

    if plot_type in ("tiled", "all"):
        if _HAS_MPL and _HAS_NX:
            path = plot_tiled_clusters(
                filtered_edges, filtered_nodes,
                gene_cluster_map, cluster_genes, cluster_sizes,
                str(out_dir / f"{prefix.name}_clusters_tiled.png"),
                gene_species_map=gene_species_map,
                species_list=species_list,
                species_color_map=species_color_map,
                top_clusters=top_clusters, layout=layout,
                figsize=figsize, dpi=dpi,
            )
            results["tiled"] = path
        else:
            log_warn("matplotlib/networkx not available, skip tiled plot")

    if plot_type in ("interactive", "all") or interactive:
        if _HAS_NX:
            path = plot_interactive_app(
                filtered_edges, filtered_nodes,
                str(out_dir / f"{prefix.name}_interactive.html"),
                gene_species_map=gene_species_map,
                gene_cluster_map=gene_cluster_map,
                species_list=species_list,
                species_color_map=species_color_map,
                cluster_genes=cluster_genes,
                cluster_sizes=cluster_sizes,
                top_k=top_k,
            )
            results["interactive"] = path
        else:
            log_warn("networkx not available, skip interactive plot")

    if _HAS_MPL:
        path = plot_cluster_size_dist(
            cluster_sizes,
            str(out_dir / f"{prefix.name}_cluster_dist.png"),
            dpi=dpi,
        )
        results["cluster_dist"] = path

        if gene_species_map:
            path = plot_species_composition(
                cluster_genes, species_list, gene_species_map,
                str(out_dir / f"{prefix.name}_species_composition.png"),
                species_color_map=species_color_map,
                dpi=dpi,
            )
            results["species_composition"] = path

    log_info("Visualization completed!")
    return {"success": True, "outputs": results}


def main():
    parser = argparse.ArgumentParser(
        description="Step 5: Visualize synteny network (from cluster output)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python viz.py -c Filtered.clusters.tsv -n Final_Network.tsv -s species.lst -d seqs
  python viz.py -c Filtered.clusters.tsv -n Final_Network.tsv -s species.lst --palette palette.tsv
  python viz.py -c Filtered.clusters.tsv -n Final_Network.tsv -s species.lst --plot-type static
  python viz.py -c Filtered.clusters.tsv -n Final_Network.tsv -s species.lst --plot-type tiled --top-clusters 20
  python viz.py -c Filtered.clusters.tsv -n Final_Network.tsv -s species.lst --plot-type interactive
  python viz.py -c Filtered.clusters.tsv -n Final_Network.tsv -s species.lst --layout circular
  python viz.py -c Filtered.clusters.tsv -n Final_Network.tsv -s species.lst --top-k 100
        """,
    )

    parser.add_argument("-c", "--cluster-tsv", required=True,
                        help="Cluster TSV file (from 'cluster' command, *.clusters.tsv)")
    parser.add_argument("-n", "--network-tsv", required=True,
                        help="Network TSV file (from 'network' command)")
    parser.add_argument("-s", "--species-list", required=True,
                        help="Species list file (for color coding and species mapping)")
    parser.add_argument("-d", "--bed-dir", default=".",
                        help="Directory containing .bed files (default: current dir)")
    parser.add_argument("--palette", default=None,
                        help="Custom species color palette file (TSV: species<tab>color)")
    parser.add_argument("--output-dir", default=".",
                        help="Output directory for plots (default: current dir)")
    parser.add_argument("--plot-type", type=str, default="all",
                        choices=["static", "tiled", "interactive", "all"],
                        help="Plot type: static (ggraph-style), tiled (Cytoscape-style), "
                             "interactive (HTML with controls), all (default: all)")
    parser.add_argument("--layout", type=str, default="force",
                        choices=["force", "circular", "kamada_kawai"],
                        help="Layout algorithm (default: force)")
    parser.add_argument("--top-k", type=int, default=None,
                        help="Show only top K nodes by degree in network view (default: show all)")
    parser.add_argument("--top-clusters", type=int, default=20,
                        help="Number of top clusters to show in tiled view (default: 20)")
    parser.add_argument("--interactive", action="store_true",
                        help="Force interactive HTML output")
    parser.add_argument("--dpi", type=int, default=300,
                        help="DPI for static plots (default: 300)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

    try:
        result = run_viz(
            cluster_tsv=args.cluster_tsv,
            network_tsv=args.network_tsv,
            species_list_file=args.species_list,
            bed_dir=args.bed_dir,
            palette_file=args.palette,
            output_dir=args.output_dir,
            plot_type=args.plot_type,
            layout=args.layout,
            top_k=args.top_k,
            top_clusters=args.top_clusters,
            interactive=args.interactive,
            dpi=args.dpi,
        )
        if result["success"]:
            for fmt, path in result.get("outputs", {}).items():
                log_info(f"  {fmt}: {path}")
            log_info("Done!")
        else:
            log_error("Visualization failed")
            sys.exit(1)
    except ImportError as e:
        log_error(f"Missing dependency: {e}")
        log_info("Install required packages: pip install matplotlib networkx numpy")
        sys.exit(1)


if __name__ == "__main__":
    main()
