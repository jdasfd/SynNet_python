import sys
import json
from pathlib import Path
from typing import List, Dict, Set, Optional, Tuple
from collections import defaultdict

from synnet.utils.logger import get_logger, info, warning, error, debug
from synnet.utils.io import (
    read_species_list,
    read_synnet_tsv,
    build_gene_species_map,
    ensure_dir,
)

logger = get_logger(__name__)

try:
    import networkx as nx
    _HAS_NX = True
except ImportError:
    _HAS_NX = False


DEFAULT_PALETTE = [
    '#e6194b', '#3cb44b', '#4363d8', '#f58231', '#911eb4',
    '#42d4f4', '#f032e6', '#bfef45', '#fabed4', '#469990',
    '#dcbeff', '#9A6324', '#fffac8', '#800000', '#aaffc3',
    '#808000', '#ffd8b1', '#000075', '#a9a9a9', '#000000',
]


def build_species_color_map(species_list: List[str], palette_file: Optional[str] = None) -> Dict[str, str]:
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
        output_file: Path,
        title: Optional[str] = None,
) -> str:
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

    info(f"Network: {n_nodes} nodes, {n_edges} edges, {n_clusters} clusters")

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
    html = html.replace('__TITLE__', title or 'SynNet Visualization')
    html = html.replace('__NODES_JSON__', json.dumps(vis_nodes))
    html = html.replace('__EDGES_JSON__', json.dumps(vis_edges))
    html = html.replace('__SPECIES_LIST__', json.dumps(species_list))
    html = html.replace('__SPECIES_COLORS__', json.dumps(species_color_map))
    html = html.replace('__N_CLUSTERS__', str(n_clusters))

    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html)

    info(f"Interactive plot saved: {output_file}")
    return str(output_file)


def visualize_synnet(
        synnet_file: str,
        species_list_file: str,
        bed_dir: str,
        output_dir: Optional[str] = None,
) -> Dict[str, any]:
    species_list = read_species_list(species_list_file)
    info(f"Species: {', '.join(species_list)}")

    gene_species_map = build_gene_species_map(species_list, Path(bed_dir))
    info(f"Loaded {len(gene_species_map)} gene-species mappings")

    species_color_map = build_species_color_map(species_list, None)

    synnet_path = Path(synnet_file)
    info(f"Loading: {synnet_path}")
    edges, nodes, cluster_genes = read_synnet_tsv(synnet_path)
    info(f"Loaded {len(edges)} edges, {len(nodes)} nodes, {len(cluster_genes)} clusters")

    if output_dir:
        out_dir = Path(output_dir)
    else:
        out_dir = synnet_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    output_html = out_dir / "Clusters.synnet.html"

    if not _HAS_NX:
        error("networkx is required for visualization")
        info("Install: pip install networkx")
        return {"success": False, "error": "Missing dependencies"}

    plot_interactive_network(
        edges, nodes, cluster_genes,
        gene_species_map, species_list, species_color_map,
        output_html,
    )

    info(f"Done! Output: {output_html}")
    return {"success": True, "output": str(output_html)}
