"""
Headless GPU render script for Blender.

Usage:
  blender --background assembly_graph_3d.blend --python render_gpu.py

Env vars (all optional):
  RENDER_SAMPLES   — cycle samples (default: 512)
  RENDER_RES_X     — width  (default: 1920)
  RENDER_RES_Y     — height (default: 1080)
  RENDER_OUTPUT    — output path (default: ./render.png)
  GPU_BACKEND      — CUDA, OPTIX, HIP, METAL, NONE (default: auto-detect)
"""

import bpy
import os


def enable_gpu(backend=None):
    """Activate all available GPU devices for Cycles."""
    prefs = bpy.context.preferences.addons['cycles'].preferences

    # Try backends in preference order
    backends = [backend] if backend else ['OPTIX', 'CUDA', 'HIP', 'METAL']
    activated = False

    for try_backend in backends:
        try:
            prefs.compute_device_type = try_backend
            prefs.get_devices()
            for dev in prefs.devices:
                if dev.type != 'CPU':
                    dev.use = True
                    activated = True
                    print(f"  Enabled {dev.type} device: {dev.name}")
                else:
                    dev.use = False  # CPU fallback off when GPU available
            if activated:
                print(f"  Using backend: {try_backend}")
                break
        except Exception:
            continue

    if not activated:
        print("  WARNING: No GPU found, falling back to CPU")
        prefs.compute_device_type = 'NONE'
        for dev in prefs.devices:
            dev.use = True


def configure_render():
    scene = bpy.context.scene

    # GPU
    scene.render.engine = 'CYCLES'
    scene.cycles.device = 'GPU'
    backend = os.environ.get('GPU_BACKEND')
    enable_gpu(backend)

    # Samples
    samples = int(os.environ.get('RENDER_SAMPLES', '512'))
    scene.cycles.samples = samples
    scene.cycles.use_adaptive_sampling = True
    scene.cycles.adaptive_threshold = 0.005
    scene.cycles.use_denoising = True

    # Tile size — larger tiles are faster on GPU
    try:
        scene.cycles.tile_size = 2048
    except AttributeError:
        pass

    # Resolution
    scene.render.resolution_x = int(os.environ.get('RENDER_RES_X', '1920'))
    scene.render.resolution_y = int(os.environ.get('RENDER_RES_Y', '1080'))
    scene.render.resolution_percentage = 100

    # Output
    output = os.environ.get('RENDER_OUTPUT', '/render/render.png')
    scene.render.filepath = output
    scene.render.image_settings.file_format = 'PNG'
    scene.render.image_settings.color_depth = '16'
    scene.render.image_settings.compression = 15

    print(f"  Render: {scene.render.resolution_x}x{scene.render.resolution_y}")
    print(f"  Samples: {samples}")
    print(f"  Output: {output}")


if __name__ == "__main__":
    print("=" * 60)
    print("GPU Render")
    print("=" * 60)
    configure_render()
    bpy.ops.render.render(write_still=True)
    print("=" * 60)
    print(f"Done — saved to {bpy.context.scene.render.filepath}")
    print("=" * 60)
