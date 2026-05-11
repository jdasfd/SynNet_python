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
import argparse
import pathlib
from typing import List, Dict, Set, Tuple
import collections

try:
    import networkx as nx # type: ignore
    _HAS_NX = True
except ImportError:
    _HAS_NX = False

try:
    from pyvis.network import Network # type: ignore
    _HAS_PYVIS = True
except ImportError:
    _HAS_PYVIS = False


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


def load_synnet_tsv(synnet_file: pathlib.Path) -> Tuple[List[Tuple[str, str, int, str]], Set[str], Dict[str, Set[str]]]:
    edges = []
    nodes = set()
    cluster_genes = collections.defaultdict(set)
    
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
            edges.append((cluster_id, node1, node2, score))
            nodes.add(node1)
            nodes.add(node2)
            cluster_genes[cluster_id].add(node1)
            cluster_genes[cluster_id].add(node2)
    
    return edges, nodes, dict(cluster_genes)


def build_species_color_map(species_list: List[str]) -> Dict[str, str]:
    return {sp: DEFAULT_PALETTE[i % len(DEFAULT_PALETTE)] for i, sp in enumerate(species_list)}


def generate_control_panel_html(species_list: List[str], species_color_map: Dict[str, str], n_nodes: int, n_edges: int, n_clusters: int) -> str:
    color_init_js = ", ".join([f'"{sp}": "{color}"' for sp, color in species_color_map.items()])
    
    html = """
<div id="control-panel" style="position:absolute; top:10px; right:10px; background:white; padding:15px; border:1px solid #ccc; border-radius:8px; font-size:13px; z-index:1000; max-width:250px;">
    <b style="font-size:14px;">Species Colors</b><br>
    <small style="color:#666;">Click color box to change</small><br>
"""
    for sp in species_list:
        color = species_color_map.get(sp, '#cccccc')
        html += f"""
    <div style="margin:8px 0; display:flex; align-items:center;">
        <input type="color" id="color_{sp}" value="{color}" onchange="updateColor('{sp}', this.value)" style="width:35px; height:25px; border:1px solid #ccc; cursor:pointer; border-radius:3px;">
        <span style="margin-left:8px; font-weight:500;">{sp}</span>
    </div>
"""
    
    html += """
    <hr style="margin:12px 0; border-color:#eee;">
    <b style="font-size:14px;">Export</b><br>
    <div style="margin:8px 0;">
        <button onclick="exportPNG()" style="padding:8px 15px; cursor:pointer; border-radius:4px; border:1px solid #ccc; background:#f5f5f5;">PNG</button>
        <button onclick="exportPDF()" style="padding:8px 15px; cursor:pointer; border-radius:4px; border:1px solid #ccc; background:#f5f5f5; margin-left:5px;">PDF</button>
    </div>
    <hr style="margin:12px 0; border-color:#eee;">
    <b style="font-size:14px;">Network Stats</b><br>
"""
    html += f"""
    <div style="margin:5px 0;">Nodes: <b>{n_nodes}</b></div>
    <div style="margin:5px 0;">Edges: <b>{n_edges}</b></div>
    <div style="margin:5px 0;">Clusters: <b>{n_clusters}</b></div>
</div>

<script>
var speciesColors = {{{color_init_js}}};

function initColors() {{
    if (typeof network === 'undefined' || !network.body) return;
    var nodes = network.body.data.nodes;
    var updates = [];
    nodes.forEach(function(node) {{
        var sp = node.group;
        if (sp && speciesColors[sp]) {{
            updates.push({{id: node.id, color: {{background: speciesColors[sp], border: speciesColors[sp]}}}});
        }}
    }});
    nodes.update(updates);
}}

if (typeof network !== 'undefined') {{
    network.on("stabilized", initColors);
    setTimeout(initColors, 1000);
}}

function updateColor(species, color) {{
    speciesColors[species] = color;
    if (typeof network === 'undefined') return;
    var nodes = network.body.data.nodes;
    var updates = [];
    nodes.forEach(function(node) {{
        if (node.group === species) {{
            updates.push({{id: node.id, color: {{background: color, border: color}}}});
        }}
    }});
    nodes.update(updates);
}}

function findCanvas() {{
    var canvases = document.getElementsByTagName('canvas');
    for (var i = 0; i < canvases.length; i++) {{
        if (canvases[i].width > 0 && canvases[i].height > 0) {{
            return canvases[i];
        }}
    }}
    return null;
}}

function exportPNG() {{
    var canvas = findCanvas();
    if (canvas) {{
        var link = document.createElement('a');
        link.download = 'Clusters.synnet.png';
        link.href = canvas.toDataURL('image/png', 1.0);
        link.click();
    }} else {{
        alert('Canvas not ready. Please wait for network to stabilize and try again.\\n\\nTip: You can also right-click on the network and select "Save image as..."');
    }}
}}

function exportPDF() {{
    var canvas = findCanvas();
    if (canvas) {{
        var imgData = canvas.toDataURL('image/png', 1.0);
        var win = window.open('', '_blank');
        win.document.write('<html><head><title>Clusters.synnet</title></head><body style="margin:0;text-align:center;">');
        win.document.write('<img src="' + imgData + '" style="max-width:100%;">');
        win.document.write('</body></html>');
        win.document.close();
        setTimeout(function() {{ win.print(); }}, 500);
    }} else {{
        alert('Canvas not ready. Please wait for network to stabilize and try again.\\n\\nTip: You can also use browser print (Ctrl+P) to save as PDF.');
    }}
}}
</script>
"""
    return html


def plot_interactive_network(
        edges: List[Tuple[str, str, int, str]],
        nodes: Set[str],
        cluster_genes: Dict[str, Set[str]],
        gene_species_map: Dict[str, str],
        species_list: List[str],
        species_color_map: Dict[str, str],
        output_file: pathlib.Path,
) -> None:
    if not _HAS_NX or not _HAS_PYVIS:
        raise ImportError("networkx and pyvis required for interactive plot")

    G = nx.Graph()
    for node in nodes:
        sp = gene_species_map.get(node, "unknown")
        G.add_node(node, species=sp)
    
    for cluster_id, node1, node2, score in edges:
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

    for node in G.nodes():
        sp = G.nodes[node].get('species', 'unknown')
        color = species_color_map.get(sp, '#cccccc')
        deg = degrees.get(node, 1)
        size = 25 + 35 * ((deg - min_deg) / deg_range)
        
        G.nodes[node]['color'] = color
        G.nodes[node]['size'] = size
        G.nodes[node]['group'] = sp
        G.nodes[node]['label'] = node
        G.nodes[node]['font'] = {'size': 10, 'color': '#333333'}
        G.nodes[node]['shape'] = 'dot'
        
        cluster_ids = set()
        for neighbor in G.neighbors(node):
            if G.edges[node, neighbor].get('cluster'):
                cluster_ids.add(G.edges[node, neighbor]['cluster'])
        cluster_str = ', '.join(sorted(cluster_ids)[:3])
        if len(cluster_ids) > 3:
            cluster_str += f' (+{len(cluster_ids)-3} more)'
        G.nodes[node]['title'] = f"<b>{node}</b><br>Species: {sp}<br>Degree: {deg}<br>Clusters: {cluster_str}"

    for src, tgt in G.edges():
        cluster = G.edges[src, tgt].get('cluster', '')
        G.edges[src, tgt]['width'] = 1
        G.edges[src, tgt]['color'] = '#cccccc'
        G.edges[src, tgt]['title'] = f"Cluster: {cluster}"

    net = Network(height='900px', width='100%', bgcolor='#fafafa', font_color='#333333', directed=False)
    
    for node, data in G.nodes(data=True):
        color = data.get('color', '#cccccc')
        net.add_node(
            node,
            label=data.get('label', node),
            title=data.get('title', ''),
            size=data.get('size', 25),
            group=data.get('group', ''),
            color=color,
            font=data.get('font', {'size': 10, 'color': '#333333'}),
            shape=data.get('shape', 'dot'),
        )
    
    for src, tgt, data in G.edges(data=True):
        net.add_edge(src, tgt, color='#cccccc', width=1, title=data.get('title', ''))
    
    net.show_buttons(filter_=['physics'])
    
    net.set_options("""
{
  "nodes": {
    "borderWidth": 1,
    "borderWidthSelected": 2,
    "opacity": 0.9
  },
  "edges": {
    "smooth": {
      "type": "continuous"
    }
  },
  "physics": {
    "barnesHut": {
      "gravitationalConstant": -3000,
      "centralGravity": 0.3,
      "springLength": 50
    },
    "minVelocity": 0.75
  }
}
""")

    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    net.save_graph(str(output_file))
    
    with open(output_file, 'r', encoding='utf-8') as f:
        html_content = f.read()
    
    html_content = html_content.replace(
        '<script src="lib/bindings/utils.js"></script>',
        ''
    )
    
    control_panel = generate_control_panel_html(species_list, species_color_map, n_nodes, n_edges, n_clusters)
    html_content = html_content.replace('</body>', control_panel + '</body>')
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
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

    if not _HAS_NX or not _HAS_PYVIS:
        log_error("networkx and pyvis are required")
        log_info("Install: pip install networkx pyvis")
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