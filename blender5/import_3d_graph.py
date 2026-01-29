"""
Import a 3D graph layout (from compute_3d_layout.py) into Blender.
Reads layout_3d.json, builds DNA strand tubes (Bezier curves with bevel)
colored by GC content, plus junction spheres at branch points.
"""

import bpy
import bmesh
import json
import os
import math
import colorsys
import numpy as np

# ── SETTINGS ──────────────────────────────────────────────────────────────────

SCRIPT_DIR = "/projects"
LAYOUT_PATH = os.path.join(SCRIPT_DIR, "layout_3d.json")

TUBE_RADIUS_MIN = 0.06
TUBE_RADIUS_MAX = 0.20
JUNCTION_RADIUS = 0.06
JUNCTION_RADIUS_MAX = 0.14
COLLECTION_NAME = "Assembly Graph 3D"

# Double-helix geometry
HELIX_RADIUS = 0.10       # distance from centerline to each backbone
BACKBONE_RADIUS = 0.025   # tube thickness of each backbone strand
RUNG_RADIUS = 0.012       # thickness of base-pair rungs
HELIX_PITCH = 1.0         # arc-length per full helical turn
RUNG_SPACING = 0.10       # arc-length between base-pair rungs
MIN_HELIX_PTS = 6         # minimum centerline points to build a helix

# ── GC COLOR RAMP ────────────────────────────────────────────────────────────

def gc_color(t_normalized):
    """Map a normalized value [0,1] to color.
    Low values (AT-rich) = warm amber/gold
    High values (GC-rich) = cool blue/violet
    Input should already be percentile-normalized.
    """
    stops = [
        (0.00, (1.00, 0.75, 0.15)),  # gold (very AT-rich)
        (0.20, (1.00, 0.50, 0.10)),  # amber
        (0.35, (0.85, 0.30, 0.30)),  # warm red
        (0.50, (0.60, 0.20, 0.60)),  # purple (balanced)
        (0.65, (0.30, 0.20, 0.80)),  # blue-violet
        (0.80, (0.15, 0.35, 0.90)),  # blue
        (1.00, (0.10, 0.50, 1.00)),  # bright blue (very GC-rich)
    ]

    t = max(0.0, min(1.0, t_normalized))
    if t <= 0:
        return stops[0][1] + (1.0,)
    if t >= 1:
        return stops[-1][1] + (1.0,)

    for i in range(len(stops) - 1):
        t0, c0 = stops[i]
        t1, c1 = stops[i + 1]
        if t0 <= t <= t1:
            frac = (t - t0) / (t1 - t0)
            r = c0[0] + (c1[0] - c0[0]) * frac
            g = c0[1] + (c1[1] - c0[1]) * frac
            b = c0[2] + (c1[2] - c0[2]) * frac
            return (r, g, b, 1.0)

    return stops[-1][1] + (1.0,)


def compute_gc_normalization(data):
    """Compute percentile-based normalization for GC values across all paths.
    Returns (gc_p5, gc_range) so that per-path GC can be mapped to [0,1]."""
    paths = data.get('paths', [])
    if not paths:
        return 0.0, 1.0
    all_avg_gc = [sum(p['gc']) / len(p['gc']) for p in paths if p['gc']]
    arr = np.array(all_avg_gc)
    p5 = np.percentile(arr, 5)
    p95 = np.percentile(arr, 95)
    rng = p95 - p5 if p95 > p5 else 0.01
    print(f"  GC normalization: p5={p5:.3f}, p95={p95:.3f}, range={rng:.3f}")
    return p5, rng


# ── BLENDER UTILITIES ────────────────────────────────────────────────────────

def clear_collection(name):
    if name in bpy.data.collections:
        col = bpy.data.collections[name]
        for obj in list(col.objects):
            bpy.data.objects.remove(obj, do_unlink=True)
        bpy.data.collections.remove(col)
    # Clean orphan data
    for mat in list(bpy.data.materials):
        if mat.name.startswith(("StrandMat", "JunctionMat", "EdgeMat", "NodeMat",
                                "KmerMat", "BackboneMat", "RungMat")):
            bpy.data.materials.remove(mat)
    for curve in list(bpy.data.curves):
        if curve.name.startswith("Strand_"):
            bpy.data.curves.remove(curve)
    for img in list(bpy.data.images):
        if img.name.startswith("SeqTex_"):
            bpy.data.images.remove(img)


def create_collection(name):
    col = bpy.data.collections.new(name)
    bpy.context.scene.collection.children.link(col)
    return col


def _bsdf_input(bsdf, *names):
    """Find a BSDF input by trying multiple names (handles API changes across versions)."""
    for name in names:
        if name in bsdf.inputs:
            return bsdf.inputs[name]
    return None


def _configure_principled(bsdf):
    """Configure a Principled BSDF for a wet biological/organic look."""
    bsdf.inputs['Roughness'].default_value = 0.25

    inp = _bsdf_input(bsdf, 'Specular IOR Level', 'Specular')
    if inp:
        inp.default_value = 0.6

    inp = _bsdf_input(bsdf, 'Subsurface Weight', 'Subsurface')
    if inp:
        inp.default_value = 0.15

    inp = _bsdf_input(bsdf, 'Subsurface Radius')
    if inp:
        inp.default_value = (0.1, 0.05, 0.03)

    inp = _bsdf_input(bsdf, 'Coat Weight', 'Clearcoat')
    if inp:
        inp.default_value = 0.3

    inp = _bsdf_input(bsdf, 'Coat Roughness', 'Clearcoat Roughness')
    if inp:
        inp.default_value = 0.1


def make_vcol_material(name):
    """Principled BSDF material driven by vertex color attribute."""
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    tree = mat.node_tree
    tree.nodes.clear()

    output = tree.nodes.new('ShaderNodeOutputMaterial')
    output.location = (400, 0)
    bsdf = tree.nodes.new('ShaderNodeBsdfPrincipled')
    bsdf.location = (100, 0)
    _configure_principled(bsdf)

    attr = tree.nodes.new('ShaderNodeAttribute')
    attr.location = (-200, 0)
    attr.attribute_name = "Col"
    attr.attribute_type = 'GEOMETRY'

    tree.links.new(attr.outputs['Color'], bsdf.inputs['Base Color'])
    ssc = _bsdf_input(bsdf, 'Subsurface Color')
    if ssc:
        tree.links.new(attr.outputs['Color'], ssc)
    tree.links.new(bsdf.outputs[0], output.inputs[0])
    return mat


def make_flat_material(name, color):
    """Principled BSDF material with a flat base color."""
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    tree = mat.node_tree
    tree.nodes.clear()

    output = tree.nodes.new('ShaderNodeOutputMaterial')
    output.location = (400, 0)
    bsdf = tree.nodes.new('ShaderNodeBsdfPrincipled')
    bsdf.location = (100, 0)
    _configure_principled(bsdf)
    bsdf.inputs['Base Color'].default_value = color
    ssc = _bsdf_input(bsdf, 'Subsurface Color')
    if ssc:
        ssc.default_value = color

    tree.links.new(bsdf.outputs[0], output.inputs[0])
    return mat


# ── NUCLEOTIDE COLORS ────────────────────────────────────────────────────────

NUCLEOTIDE_COLORS = {
    'A': (0.2, 0.8, 0.2, 1.0),   # green
    'T': (0.9, 0.2, 0.2, 1.0),   # red
    'C': (0.2, 0.4, 0.9, 1.0),   # blue
    'G': (0.95, 0.85, 0.1, 1.0), # yellow
}
NUCLEOTIDE_FALLBACK = (0.4, 0.4, 0.4, 1.0)  # gray for N/unknown


def make_sequence_image(name, sequence):
    """Create a 1D image (width=len(sequence), height=1) with per-base colors."""
    seq = sequence.upper()
    w = max(len(seq), 1)
    img = bpy.data.images.new(name, width=w, height=1, alpha=True)
    pixels = [0.0] * (w * 4)
    for i, base in enumerate(seq):
        c = NUCLEOTIDE_COLORS.get(base, NUCLEOTIDE_FALLBACK)
        pixels[i * 4]     = c[0]
        pixels[i * 4 + 1] = c[1]
        pixels[i * 4 + 2] = c[2]
        pixels[i * 4 + 3] = c[3]
    img.pixels[:] = pixels
    img.pack()
    return img


def make_sequence_texture_material(name, image):
    """Principled BSDF material that samples a 1D nucleotide image along the curve."""
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    tree = mat.node_tree
    tree.nodes.clear()

    output = tree.nodes.new('ShaderNodeOutputMaterial')
    output.location = (600, 0)

    bsdf = tree.nodes.new('ShaderNodeBsdfPrincipled')
    bsdf.location = (300, 0)
    _configure_principled(bsdf)

    tex = tree.nodes.new('ShaderNodeTexImage')
    tex.location = (0, 0)
    tex.image = image
    tex.interpolation = 'Closest'
    tex.extension = 'CLIP'

    # Use Generated coordinates: x goes 0→1 along the curve
    texcoord = tree.nodes.new('ShaderNodeTexCoord')
    texcoord.location = (-400, 0)

    separate = tree.nodes.new('ShaderNodeSeparateXYZ')
    separate.location = (-200, 0)

    combine = tree.nodes.new('ShaderNodeCombineXYZ')
    combine.location = (-50, 0)
    combine.inputs['Y'].default_value = 0.5
    combine.inputs['Z'].default_value = 0.0

    tree.links.new(texcoord.outputs['Generated'], separate.inputs['Vector'])
    tree.links.new(separate.outputs['X'], combine.inputs['X'])
    tree.links.new(combine.outputs['Vector'], tex.inputs['Vector'])
    tree.links.new(tex.outputs['Color'], bsdf.inputs['Base Color'])
    ssc = _bsdf_input(bsdf, 'Subsurface Color')
    if ssc:
        tree.links.new(tex.outputs['Color'], ssc)
    tree.links.new(bsdf.outputs[0], output.inputs[0])

    return mat


# ── DOUBLE-HELIX GEOMETRY ────────────────────────────────────────────────────

def parallel_transport_frames(pts):
    """Compute smoothly-varying coordinate frames along a polyline.
    Returns (tangents, normals, binormals) arrays, each (N, 3)."""
    n = len(pts)
    tangents = np.zeros((n, 3))
    tangents[0] = pts[1] - pts[0]
    tangents[-1] = pts[-1] - pts[-2]
    for i in range(1, n - 1):
        tangents[i] = pts[i + 1] - pts[i - 1]
    norms = np.linalg.norm(tangents, axis=1, keepdims=True)
    tangents /= np.maximum(norms, 1e-8)

    # Initial normal: perpendicular to first tangent
    t0 = tangents[0]
    up = np.array([1.0, 0.0, 0.0]) if abs(t0[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    normal = np.cross(t0, up)
    normal /= np.linalg.norm(normal)

    normals = np.zeros((n, 3))
    binormals = np.zeros((n, 3))
    normals[0] = normal
    binormals[0] = np.cross(tangents[0], normal)

    for i in range(1, n):
        b = np.cross(tangents[i - 1], tangents[i])
        b_len = np.linalg.norm(b)
        if b_len < 1e-8:
            normals[i] = normals[i - 1]
        else:
            b /= b_len
            angle = np.arccos(np.clip(np.dot(tangents[i - 1], tangents[i]), -1, 1))
            cos_a, sin_a = np.cos(angle), np.sin(angle)
            normals[i] = (normals[i - 1] * cos_a
                          + np.cross(b, normals[i - 1]) * sin_a
                          + b * np.dot(b, normals[i - 1]) * (1 - cos_a))
        binormals[i] = np.cross(tangents[i], normals[i])

    return tangents, normals, binormals


def generate_helix_points(centerline, normals, binormals, arc_lengths,
                          helix_radius, helix_pitch, phase=0.0):
    """Offset centerline points along a helical path."""
    angles = 2.0 * np.pi * arc_lengths / helix_pitch + phase
    offset = (np.cos(angles)[:, None] * normals
              + np.sin(angles)[:, None] * binormals) * helix_radius
    return centerline + offset


def make_backbone_curve(name, helix_pts, radius, collection):
    """Create a smooth NURBS bevel curve from helix points."""
    n = len(helix_pts)
    curve_data = bpy.data.curves.new(name, type='CURVE')
    curve_data.dimensions = '3D'
    curve_data.bevel_depth = radius
    curve_data.bevel_resolution = 3
    curve_data.use_fill_caps = True
    curve_data.resolution_u = 12  # smooth interpolation between control points

    spline = curve_data.splines.new('NURBS')
    spline.points.add(n - 1)
    for i in range(n):
        p = helix_pts[i]
        spline.points[i].co = (p[0], p[1], p[2], 1.0)
    spline.use_endpoint_u = True  # curve passes through first/last points
    spline.order_u = 4            # cubic NURBS

    obj = bpy.data.objects.new(name, curve_data)
    collection.objects.link(obj)
    return obj


def build_rung_batch(rung_pairs, rung_colors, collection, mat):
    """Build all base-pair rungs as a single batched mesh with vertex colors."""
    if not rung_pairs:
        return

    bm = bmesh.new()

    # Template: thin cylinder along Z, 6 sides
    seg_count = 6
    template_ring = []
    for i in range(seg_count):
        angle = 2 * math.pi * i / seg_count
        template_ring.append((math.cos(angle) * RUNG_RADIUS,
                              math.sin(angle) * RUNG_RADIUS))

    color_records = []

    for (p0, p1), color in zip(rung_pairs, rung_colors):
        p0 = np.array(p0)
        p1 = np.array(p1)
        axis = p1 - p0
        length = np.linalg.norm(axis)
        if length < 1e-6:
            continue
        axis_n = axis / length

        # Build local frame
        up = np.array([0.0, 0.0, 1.0]) if abs(axis_n[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        side = np.cross(axis_n, up)
        side /= np.linalg.norm(side)
        up2 = np.cross(side, axis_n)

        vert_start = len(bm.verts)

        # Two rings of vertices
        verts_bottom = []
        verts_top = []
        for cx, cy in template_ring:
            v_b = p0 + side * cx + up2 * cy
            v_t = p1 + side * cx + up2 * cy
            verts_bottom.append(bm.verts.new(v_b))
            verts_top.append(bm.verts.new(v_t))

        # Side faces
        for i in range(seg_count):
            j = (i + 1) % seg_count
            try:
                bm.faces.new([verts_bottom[i], verts_bottom[j],
                              verts_top[j], verts_top[i]])
            except ValueError:
                pass

        # Cap faces
        try:
            bm.faces.new(verts_bottom)
        except ValueError:
            pass
        try:
            bm.faces.new(list(reversed(verts_top)))
        except ValueError:
            pass

        vert_count = seg_count * 2
        color_records.append((vert_start, vert_count, color))

    mesh = bpy.data.meshes.new("RungMesh")
    bm.to_mesh(mesh)
    bm.free()

    if hasattr(mesh, 'color_attributes') and color_records:
        ca = mesh.color_attributes.new("Col", 'BYTE_COLOR', 'POINT')
        for vert_start, vert_count, color in color_records:
            for vi in range(vert_start, vert_start + vert_count):
                if vi < len(ca.data):
                    ca.data[vi].color = color

    obj = bpy.data.objects.new("BaseRungs", mesh)
    obj.data.materials.append(mat)
    collection.objects.link(obj)
    return obj


# ── STRAND BUILDING ──────────────────────────────────────────────────────────

def build_strands(data, collection, gc_p5, gc_range, kmer_palette=None):
    """Build double-helix DNA strands: two twisted backbone tubes + base-pair rungs."""
    paths = data.get('paths', [])
    if not paths:
        print("No path data found in layout JSON")
        return

    print(f"Building {len(paths)} double-helix strands...")

    strand_col = bpy.data.collections.new("Strands")
    collection.children.link(strand_col)

    # Shared material for rungs (vertex-colored)
    rung_mat = make_vcol_material("RungMat")

    # Backbone material: muted phosphate-backbone color
    backbone_color_a = (0.85, 0.55, 0.20, 1.0)  # warm orange
    backbone_color_b = (0.50, 0.60, 0.85, 1.0)  # cool blue
    backbone_mat_a = make_flat_material("BackboneMat_A", backbone_color_a)
    backbone_mat_b = make_flat_material("BackboneMat_B", backbone_color_b)

    all_rung_pairs = []
    all_rung_colors = []

    for pi, path_data in enumerate(paths):
        points = path_data['points']
        n_pts = len(points)
        if n_pts < 2:
            continue

        # Subsample centerline — keep density high for smooth helices
        step = max(1, n_pts // 500)
        indices = list(range(0, n_pts, step))
        if indices[-1] != n_pts - 1:
            indices.append(n_pts - 1)
        sampled_pts = np.array([points[i] for i in indices])
        n_sampled = len(sampled_pts)

        if n_sampled < MIN_HELIX_PTS:
            # Too short for helix — single smooth tube fallback
            curve_data = bpy.data.curves.new(f"Strand_{pi}", type='CURVE')
            curve_data.dimensions = '3D'
            curve_data.bevel_depth = BACKBONE_RADIUS * 2
            curve_data.bevel_resolution = 3
            curve_data.use_fill_caps = True
            curve_data.resolution_u = 12
            spline = curve_data.splines.new('NURBS')
            spline.points.add(n_sampled - 1)
            for i in range(n_sampled):
                p = sampled_pts[i]
                spline.points[i].co = (p[0], p[1], p[2], 1.0)
            spline.use_endpoint_u = True
            spline.order_u = min(4, n_sampled)
            obj = bpy.data.objects.new(f"Strand_{pi}_{path_data['name']}", curve_data)
            obj.data.materials.append(backbone_mat_a)
            strand_col.objects.link(obj)
            continue

        # Compute frames
        tangents, normals, binormals = parallel_transport_frames(sampled_pts)

        # Arc lengths
        diffs = np.diff(sampled_pts, axis=0)
        seg_lens = np.linalg.norm(diffs, axis=1)
        arc = np.zeros(n_sampled)
        arc[1:] = np.cumsum(seg_lens)

        # Generate two helical backbone curves (180 degrees apart)
        helix_a = generate_helix_points(sampled_pts, normals, binormals,
                                        arc, HELIX_RADIUS, HELIX_PITCH, phase=0.0)
        helix_b = generate_helix_points(sampled_pts, normals, binormals,
                                        arc, HELIX_RADIUS, HELIX_PITCH, phase=np.pi)

        obj_a = make_backbone_curve(f"Strand_{pi}a_{path_data['name']}",
                                    helix_a, BACKBONE_RADIUS, strand_col)
        obj_b = make_backbone_curve(f"Strand_{pi}b_{path_data['name']}",
                                    helix_b, BACKBONE_RADIUS, strand_col)
        obj_a.data.materials.append(backbone_mat_a)
        obj_b.data.materials.append(backbone_mat_b)

        # Base-pair rungs: place at regular arc-length intervals
        seq = path_data.get('sequence', '')
        total_arc = arc[-1]
        if total_arc < RUNG_SPACING:
            continue

        rung_arcs = np.arange(RUNG_SPACING, total_arc - RUNG_SPACING * 0.5, RUNG_SPACING)
        for ra in rung_arcs:
            # Interpolate position on each helix at this arc length
            idx = np.searchsorted(arc, ra, side='right') - 1
            idx = max(0, min(idx, n_sampled - 2))
            frac = (ra - arc[idx]) / max(seg_lens[idx], 1e-8)
            frac = np.clip(frac, 0, 1)

            pa = helix_a[idx] + frac * (helix_a[idx + 1] - helix_a[idx])
            pb = helix_b[idx] + frac * (helix_b[idx + 1] - helix_b[idx])

            # Color rung by nucleotide at this position
            if seq:
                seq_frac = ra / total_arc
                base_idx = int(seq_frac * len(seq))
                base_idx = min(base_idx, len(seq) - 1)
                base = seq[base_idx].upper()
                color = NUCLEOTIDE_COLORS.get(base, NUCLEOTIDE_FALLBACK)
            else:
                color = NUCLEOTIDE_FALLBACK

            all_rung_pairs.append((pa.tolist(), pb.tolist()))
            all_rung_colors.append(color)

        if (pi + 1) % 50 == 0:
            print(f"  {pi+1}/{len(paths)} strands...")

    # Build all rungs as one batched mesh
    if all_rung_pairs:
        print(f"  Building {len(all_rung_pairs)} base-pair rungs...")
        build_rung_batch(all_rung_pairs, all_rung_colors, strand_col, rung_mat)

    print(f"  {len(paths)} strands built")


# ── JUNCTION NODES ───────────────────────────────────────────────────────────

def build_junctions(data, collection, kmer_palette=None):
    """Build spheres at branch/junction nodes (degree >= 3).
    Uses k-mer coloring if available, falling back to GC coloring."""
    nodes_list = data['nodes']
    junction_nodes = [n for n in nodes_list if n['degree'] >= 4]

    if not junction_nodes:
        print("No junction nodes found")
        return

    print(f"Building {len(junction_nodes)} junction spheres...")

    # Compute overall average GC for fallback
    paths = data.get('paths', [])
    all_gc = [sum(p['gc'])/len(p['gc']) for p in paths if p['gc']] if paths else [0.5]
    default_gc = sum(all_gc) / len(all_gc)

    mat = make_vcol_material("JunctionMat")

    # Build all junction spheres as merged mesh
    bm = bmesh.new()
    template_bm = bmesh.new()
    bmesh.ops.create_icosphere(template_bm, subdivisions=2, radius=1.0)
    template_verts = [(v.co.x, v.co.y, v.co.z) for v in template_bm.verts]
    template_faces = [tuple(v.index for v in f.verts) for f in template_bm.faces]
    template_vcount = len(template_bm.verts)
    template_bm.free()

    color_records = []
    max_deg = max(n['degree'] for n in junction_nodes)
    min_deg = min(n['degree'] for n in junction_nodes)
    deg_range = max_deg - min_deg if max_deg > min_deg else 1

    for n in junction_nodes:
        t = (n['degree'] - min_deg) / deg_range
        r = JUNCTION_RADIUS + (JUNCTION_RADIUS_MAX - JUNCTION_RADIUS) * t
        px, py, pz = n['x'], n['y'], n['z']

        vert_start = len(bm.verts)
        new_verts = []
        for vx, vy, vz in template_verts:
            new_verts.append(bm.verts.new((vx * r + px, vy * r + py, vz * r + pz)))
        for face_indices in template_faces:
            try:
                bm.faces.new([new_verts[fi] for fi in face_indices])
            except ValueError:
                pass

        # Color by k-mer if palette available, else fall back to GC
        color = None
        if kmer_palette and 'kmer_color_index' in n:
            ci = n['kmer_color_index']
            if ci is not None and ci < len(kmer_palette):
                rgb = kmer_palette[ci]
                color = (rgb[0], rgb[1], rgb[2], 1.0)

        if color is None:
            color = gc_color(default_gc)

        # Brighten junction nodes
        color = (min(color[0] * 1.5, 1.0), min(color[1] * 1.5, 1.0),
                 min(color[2] * 1.5, 1.0), 1.0)
        color_records.append((vert_start, template_vcount, color))

    mesh = bpy.data.meshes.new("JunctionMesh")
    bm.to_mesh(mesh)
    bm.free()

    if hasattr(mesh, 'color_attributes'):
        ca = mesh.color_attributes.new("Col", 'BYTE_COLOR', 'POINT')
        for vert_start, vert_count, color in color_records:
            for vi in range(vert_start, vert_start + vert_count):
                if vi < len(ca.data):
                    ca.data[vi].color = color

    obj = bpy.data.objects.new("Junctions", mesh)
    obj.data.materials.append(mat)
    collection.objects.link(obj)
    print(f"  {len(junction_nodes)} junctions ({len(mesh.vertices)} verts)")


# ── SCENE SETUP ──────────────────────────────────────────────────────────────

def setup_scene(data, collection):
    nodes_list = data['nodes']
    xs = [n['x'] for n in nodes_list]
    ys = [n['y'] for n in nodes_list]
    zs = [n['z'] for n in nodes_list]

    cx = (min(xs) + max(xs)) / 2
    cy = (min(ys) + max(ys)) / 2
    cz = (min(zs) + max(zs)) / 2
    span = max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))

    # ── Camera with depth of field ──
    cam_data = bpy.data.cameras.new("GraphCamera3D")
    cam_data.type = 'PERSP'
    cam_data.lens = 50  # tighter lens for shallower DOF look
    cam_data.clip_end = span * 10
    cam_data.dof.use_dof = True
    cam_data.dof.aperture_fstop = 2.8
    cam = bpy.data.objects.new("GraphCamera3D", cam_data)
    cam.location = (cx + span * 0.6, cy - span * 0.6, cz + span * 0.4)
    collection.objects.link(cam)
    bpy.context.scene.camera = cam

    # Track to center, also use as DOF focus target
    target = bpy.data.objects.new("CamTarget", None)
    target.location = (cx, cy, cz)
    collection.objects.link(target)
    constraint = cam.constraints.new('TRACK_TO')
    constraint.target = target
    constraint.track_axis = 'TRACK_NEGATIVE_Z'
    constraint.up_axis = 'UP_Y'
    cam_data.dof.focus_object = target

    # ── Lighting: 3-point setup ──
    # Key light — warm, strong
    key = bpy.data.lights.new("KeyLight", 'AREA')
    key.energy = span * 80
    key.size = span * 0.5
    key.color = (1.0, 0.95, 0.9)
    key_obj = bpy.data.objects.new("KeyLight", key)
    key_obj.location = (cx + span * 0.8, cy - span * 0.4, cz + span * 0.9)
    collection.objects.link(key_obj)
    kc = key_obj.constraints.new('TRACK_TO')
    kc.target = target
    kc.track_axis = 'TRACK_NEGATIVE_Z'
    kc.up_axis = 'UP_Y'

    # Fill light — cool, softer
    fill = bpy.data.lights.new("FillLight", 'AREA')
    fill.energy = span * 30
    fill.size = span * 0.8
    fill.color = (0.8, 0.85, 1.0)
    fill_obj = bpy.data.objects.new("FillLight", fill)
    fill_obj.location = (cx - span * 0.6, cy + span * 0.5, cz + span * 0.3)
    collection.objects.link(fill_obj)
    fc = fill_obj.constraints.new('TRACK_TO')
    fc.target = target
    fc.track_axis = 'TRACK_NEGATIVE_Z'
    fc.up_axis = 'UP_Y'

    # Rim light — accent from behind
    rim = bpy.data.lights.new("RimLight", 'AREA')
    rim.energy = span * 50
    rim.size = span * 0.3
    rim.color = (0.9, 0.95, 1.0)
    rim_obj = bpy.data.objects.new("RimLight", rim)
    rim_obj.location = (cx - span * 0.3, cy + span * 0.2, cz - span * 0.6)
    collection.objects.link(rim_obj)
    rc = rim_obj.constraints.new('TRACK_TO')
    rc.target = target
    rc.track_axis = 'TRACK_NEGATIVE_Z'
    rc.up_axis = 'UP_Y'

    # ── World: dark with subtle volumetric fog ──
    world = bpy.data.worlds.get("World") or bpy.data.worlds.new("World")
    bpy.context.scene.world = world
    world.use_nodes = True
    tree = world.node_tree
    tree.nodes.clear()

    bg = tree.nodes.new('ShaderNodeBackground')
    bg.location = (0, 0)
    bg.inputs['Color'].default_value = (0.005, 0.005, 0.015, 1.0)
    bg.inputs['Strength'].default_value = 0.3

    vol_scatter = tree.nodes.new('ShaderNodeVolumeScatter')
    vol_scatter.location = (0, -200)
    vol_scatter.inputs['Color'].default_value = (0.7, 0.75, 0.85, 1.0)
    vol_scatter.inputs['Density'].default_value = 0.008
    vol_scatter.inputs['Anisotropy'].default_value = 0.3

    output = tree.nodes.new('ShaderNodeOutputWorld')
    output.location = (300, 0)
    tree.links.new(bg.outputs[0], output.inputs['Surface'])
    tree.links.new(vol_scatter.outputs[0], output.inputs['Volume'])

    # ── Render settings: Cycles ──
    scene = bpy.context.scene
    scene.render.engine = 'CYCLES'
    scene.cycles.device = 'GPU'
    scene.cycles.samples = 256
    scene.cycles.use_denoising = True
    scene.cycles.use_adaptive_sampling = True
    scene.cycles.adaptive_threshold = 0.01
    scene.cycles.max_bounces = 8
    scene.cycles.diffuse_bounces = 4
    scene.cycles.glossy_bounces = 4
    scene.cycles.transmission_bounces = 4
    scene.cycles.volume_bounces = 2
    scene.render.resolution_x = 1920
    scene.render.resolution_y = 1080
    scene.render.film_transparent = False

    # Color management — filmic for realistic dynamic range
    try:
        scene.view_settings.view_transform = 'Filmic'
    except TypeError:
        scene.view_settings.view_transform = 'AgX'
    try:
        scene.view_settings.look = 'Medium Contrast'
    except TypeError:
        try:
            scene.view_settings.look = 'AgX - Medium Contrast'
        except TypeError:
            pass
    scene.view_settings.exposure = 0.5

    # ── Compositor: subtle grain + glare ──
    scene.use_nodes = True
    comp = getattr(scene, 'node_tree', None)
    if comp is None:
        print("  Compositor node tree not available, skipping post-processing")
    else:
        comp.nodes.clear()

        rl = comp.nodes.new('CompositorNodeRLayers')
        rl.location = (0, 0)

        # Soft glare (bloom-like in Cycles)
        glare = comp.nodes.new('CompositorNodeGlare')
        glare.location = (200, 0)
        glare.glare_type = 'FOG_GLOW'
        glare.threshold = 1.2
        glare.size = 6
        glare.quality = 'HIGH'

        # Lens distortion for subtle chromatic aberration
        lens = comp.nodes.new('CompositorNodeLensdist')
        lens.location = (400, 0)
        lens.inputs['Dispersion'].default_value = 0.01

        comp_out = comp.nodes.new('CompositorNodeComposite')
        comp_out.location = (600, 0)

        comp.links.new(rl.outputs['Image'], glare.inputs['Image'])
        comp.links.new(glare.outputs['Image'], lens.inputs['Image'])
        comp.links.new(lens.outputs['Image'], comp_out.inputs['Image'])


# ── MAIN ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("3D DNA Strands → Blender")
    print("=" * 60)

    with open(LAYOUT_PATH, 'r') as f:
        data = json.load(f)

    stats = data['stats']
    print(f"  {stats['num_nodes']} nodes, {stats['num_edges']} edges")
    print(f"  {stats.get('num_paths', 0)} paths")

    clear_collection(COLLECTION_NAME)
    col = create_collection(COLLECTION_NAME)

    kmer_palette = data.get('kmer_palette')
    if kmer_palette:
        print(f"  K-mer palette: {len(kmer_palette)} colors")
    else:
        print("  No k-mer palette found, using GC coloring")

    gc_p5, gc_range = compute_gc_normalization(data)
    build_strands(data, col, gc_p5, gc_range, kmer_palette=kmer_palette)
    build_junctions(data, col, kmer_palette=kmer_palette)
    setup_scene(data, col)

    blend_path = os.path.join(SCRIPT_DIR, "assembly_graph_3d.blend")
    bpy.ops.wm.save_as_mainfile(filepath=blend_path)
    print(f"Saved to {blend_path}")
    print("=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
