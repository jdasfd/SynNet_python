"""
test/viz.py - Step 5: Visualize synteny network (interactive)

Input: Clusters.synnet.tsv from cluster.py
Output: Clusters.synnet.html (interactive network visualization)

Features:
  - Interactive network with species-colored nodes
  - Adjustable species colors in the interface (color picker)
  - Export to PNG/PDF directly from the browser
  - Physics simulation controls

Usage:
    python viz.py -i network_output/Clusters.synnet.tsv -s species.lst -b seqs
"""

import sys
import json
import argparse
import pathlib
from typing import List, Dict, Set, Tuple
from collections import defaultdict

try:
    import networkx as nx # type: ignore
    _HAS_NX = True
except ImportError:
    _HAS_NX = False


DEFAULT_PALETTE = [
    '#e6194b', '#3cb44b', '#4363d8', '#f58231', '#911eb4',
    '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990',
    '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3',
    '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#000000',
]


def log_info(msg):
    print(f"[INFO] {msg}", file=sys.stderr)


def log_warn(msg):
    print(f"[WARN] {msg}", file=sys.stderr)


def log_error(msg):
    print(f"[ERROR] {msg}", file=sys.stderr)


def load_species_list(list_file: str) -> List[str]:
    with open(list_file, 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]


def build_gene_species_map(species_list: List[str], bed_dir: pathlib.Path) -> Dict[str, str]:
    gene_map = {}
    for sp in species_list:
        bed_file = bed_dir / f"{sp}.bed"
        if not bed_file.exists():
            log_warn(f"BED file not found: {bed_file}")
            continue
        with open(bed_file, 'r') as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 4:
                    gene_map[parts[3]] = sp
    return gene_map


def load_synnet_tsv(synnet_file: pathlib.Path) -> Tuple[List[Tuple[str, int, str, str]], Set[str], Dict[str, Set[str]]]:
    edges = []
    nodes = set()
    cluster_genes = defaultdict(set)

    with open(synnet_file, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 4:
                continue
            cluster_id, score_str, node1, node2 = parts[0], parts[1], parts[2], parts[3]
            try:
                score = int(score_str)
            except ValueError:
                score = 0
            edges.append((cluster_id, score, node1, node2))
            nodes.add(node1)
            nodes.add(node2)
            cluster_genes[cluster_id].add(node1)
            cluster_genes[cluster_id].add(node2)

    return edges, nodes, dict(cluster_genes)


def build_species_color_map(species_list: List[str]) -> Dict[str, str]:
    return {sp: DEFAULT_PALETTE[i % len(DEFAULT_PALETTE)] for i, sp in enumerate(species_list)}


_INTERACTIVE_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>__TITLE__</title>
<script type="text/javascript" src="https://unpkg.com/vis-network@9.1.6/standalone/umd/vis-network.min.js"></script>
<style>
*{margin:0;padding:0;box-sizing:border-box}
body{font-family:'Segoe UI',system-ui,sans-serif;background:#fafafa}
#network{width:100%;height:900px;background:#fff}
#panel{position:absolute;top:10px;right:10px;background:#fff;padding:15px;border:1px solid #ccc;border-radius:8px;font-size:13px;z-index:1000;max-width:260px;box-shadow:0 2px 8px rgba(0,0,0,0.1)}
#panel h3{font-size:14px;margin-bottom:8px;color:#333}
#panel small{color:#666}
.color-row{margin:8px 0;display:flex;align-items:center}
.color-row input{width:35px;height:25px;border:1px solid #ccc;border-radius:3px;cursor:pointer}
.color-row span{margin-left:8px;font-weight:500}
.btn-row{margin:10px 0}
.btn{padding:8px 15px;border:none;border-radius:4px;cursor:pointer;font-size:12px;margin-right:5px}
.btn-primary{background:#4363d8;color:#fff}
.btn-secondary{background:#e8e8e8;color:#333}
.stat{margin:5px 0;font-size:12px}
.stat b{color:#4363d8}
hr{margin:12px 0;border:none;border-top:1px solid #eee}
</style>
</head>
<body>
<div id="network"></div>
<div id="panel">
<h3>Species Colors</h3>
<small>Click color box to change</small>
<div id="color-controls"></div>
<hr>
<h3>Export</h3>
<div class="btn-row">
<button class="btn btn-primary" onclick="exportPNG()">PNG</button>
<button class="btn btn-secondary" onclick="exportPDF()">PDF</button>
</div>
<hr>
<h3>Network Stats</h3>
<div class="stat">Nodes: <b id="stat-nodes">0</b></div>
<div class="stat">Edges: <b id="stat-edges">0</b></div>
<div class="stat">Clusters: <b id="stat-clusters">0</b></div>
</div>
<script>
var nodesData=new vis.DataSet(__NODES_JSON__);
var edgesData=new vis.DataSet(__EDGES_JSON__);
var container=document.getElementById('network');
var data={nodes:nodesData,edges:edgesData};
var options={
  physics:{barnesHut:{gravitationalConstant:-3000,centralGravity:0.3,springLength:50},minVelocity:0.75},
  nodes:{borderWidth:1,borderWidthSelected:2,opacity:0.9,font:{size:10,color:'#333'}},
  edges:{smooth:{type:'continuous'},color:{color:'#cccccc',opacity:0.6}},
  interaction:{hover:true,tooltipDelay:200,navigationButtons:true}
};
var network=new vis.Network(container,data,options);

document.getElementById('stat-nodes').textContent=nodesData.length;
document.getElementById('stat-edges').textContent=edgesData.length;
document.getElementById('stat-clusters').textContent='__N_CLUSTERS__';

var spColors=__SPECIES_COLORS__;
var spList=__SPECIES_LIST__;
var colorHtml='';
spList.forEach(function(sp){
  var c=spColors[sp]||'#cccccc';
  colorHtml+='<div class="color-row"><input type="color" id="color_'+sp+'" value="'+c+'" onchange="updateColor(\''+sp+'\',this.value)"><span>'+sp+'</span></div>';
});
document.getElementById('color-controls').innerHTML=colorHtml;

function updateColor(species,color){
  var updates=[];
  nodesData.forEach(function(n){
    if(n.group===species){updates.push({id:n.id,color:{background:color,border:color}});}
  });
  nodesData.update(updates);
}

function exportPNG(){
  var canvas=document.querySelector('#network canvas');
  if(canvas){
    var link=document.createElement('a');
    link.download='Clusters.synnet.png';
    link.href=canvas.toDataURL('image/png',1.0);
    link.click();
  }else{alert('Canvas not found');}
}

function exportPDF(){
  var canvas=document.querySelector('#network canvas');
  if(canvas){
    var imgData=canvas.toDataURL('image/png',1.0);
    var win=window.open('','_blank');
    win.document.write('<html><head><title>Clusters.synnet</title></head><body style="margin:0"><img src="'+imgData+'" style="max-width:100%"></body></html>');
    win.document.close();
    setTimeout(function(){win.print();},500);
  }else{alert('Canvas not found');}
}
</script>
</body>
</html>"""


def plot_interactive_network(
        edges: List[Tuple[str, int, str, str]],
        nodes: Set[str],
        cluster_genes: Dict[str, Set[str]],
        gene_species_map: Dict[str, str],
        species_list: List[str],
        species_color_map: Dict[str, str],
        output_file: pathlib.Path,
) -> None:
    if not _HAS_NX:
        raise ImportError("networkx required for interactive plot")

    G = nx.Graph()
    for node in nodes:
        sp = gene_species_map.get(node, "unknown")
        G.add_node(node, species=sp, group=sp)

    for cluster_id, score, node1, node2 in edges:
        if G.has_edge(node1, node2):
            continue
        G.add_edge(node1, node2, cluster=cluster_id)

    n_nodes = G.number_of_nodes()
    n_edges = G.number_of_edges()
    n_clusters = len(cluster_genes)

    log_info(f"Network: {n_nodes} nodes, {n_edges} edges, {n_clusters} clusters")

    degrees = dict(G.degree())
    max_deg = max(degrees.values()) if degrees else 1
    min_deg = min(degrees.values()) if degrees else 1
    deg_range = max(max_deg - min_deg, 1)

    vis_nodes = []
    for node in G.nodes():
        sp = G.nodes[node].get('species', 'unknown')
        color = species_color_map.get(sp, '#cccccc')
        deg = degrees.get(node, 1)
        size = 25 + 35 * ((deg - min_deg) / deg_range)

        cluster_ids = set()
        for neighbor in G.neighbors(node):
            if G.edges[node, neighbor].get('cluster'):
                cluster_ids.add(G.edges[node, neighbor]['cluster'])
        cluster_str = ', '.join(sorted(cluster_ids)[:3])
        if len(cluster_ids) > 3:
            cluster_str += f' (+{len(cluster_ids)-3} more)'

        vis_nodes.append({
            "id": node, "label": node, "group": sp,
            "color": {"background": color, "border": color},
            "size": round(size, 1),
            "title": f"<b>{node}</b><br>Species: {sp}<br>Degree: {deg}<br>Clusters: {cluster_str}",
            "font": {"size": 10, "color": "#333"},
        })

    vis_edges = []
    for src, tgt in G.edges():
        cluster = G.edges[src, tgt].get('cluster', '')
        vis_edges.append({
            "from": src, "to": tgt,
            "width": 1,
            "color": {"color": "#cccccc", "opacity": 0.6},
            "title": f"Cluster: {cluster}",
        })

    html = _INTERACTIVE_TEMPLATE
    html = html.replace('__TITLE__', 'SynNet Visualization')
    html = html.replace('__NODES_JSON__', json.dumps(vis_nodes))
    html = html.replace('__EDGES_JSON__', json.dumps(vis_edges))
    html = html.replace('__SPECIES_LIST__', json.dumps(species_list))
    html = html.replace('__SPECIES_COLORS__', json.dumps(species_color_map))
    html = html.replace('__N_CLUSTERS__', str(n_clusters))

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html)

    log_info(f"Interactive plot saved: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Step 5: Visualize synteny network (interactive)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Features:
  - Interactive network with species-colored nodes
  - Color picker to change species colors
  - Export to PNG/PDF from the browser
  - Physics controls in the interface

Examples:
  python viz.py -i network_output/Clusters.synnet.tsv -s species.lst -b seqs
        """,
    )

    parser.add_argument("-i", "--input", required=True,
                        help="Input Clusters.synnet.tsv file from cluster.py")
    parser.add_argument("-s", "--species-list", required=True,
                        help="Species list file")
    parser.add_argument("-b", "--bed-dir", required=True,
                        help="Directory containing .bed files")
    parser.add_argument("-o", "--output-dir", default=None,
                        help="Output directory (default: same as input file)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")

    args = parser.parse_args()

    if not _HAS_NX:
        log_error("networkx is required")
        log_info("Install: pip install networkx")
        sys.exit(1)

    log_info("Step 5: Visualize Synteny Network")

    input_file = pathlib.Path(args.input)
    if not input_file.exists():
        log_error(f"Input file not found: {input_file}")
        sys.exit(1)

    species_list = load_species_list(args.species_list)
    log_info(f"Species: {', '.join(species_list)}")

    bed_dir = pathlib.Path(args.bed_dir)
    gene_species_map = build_gene_species_map(species_list, bed_dir)
    log_info(f"Loaded {len(gene_species_map)} gene-species mappings")

    species_color_map = build_species_color_map(species_list)

    log_info(f"Loading: {input_file}")
    edges, nodes, cluster_genes = load_synnet_tsv(input_file)
    log_info(f"Loaded {len(edges)} edges, {len(nodes)} nodes, {len(cluster_genes)} clusters")

    if args.output_dir:
        output_dir = pathlib.Path(args.output_dir)
    else:
        output_dir = input_file.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    output_html = output_dir / "Clusters.synnet.html"
    plot_interactive_network(edges, nodes, cluster_genes, gene_species_map, species_list, species_color_map, output_html)

    log_info("Done!")


if __name__ == "__main__":
    main()
