import colorcet as cc
import json
import numpy as np
import pyvista as pv

from jmm.multiple_arrivals import Domain, PointSourceField, \
    ReflectedField, DiffractedField
from jmm.plot import *

from PIL import Image

pv.set_plot_theme('document')

############################################################################
# Load domain

grid = pv.read(f'sethian_shadow/Building.vtu')
with open('sethian_shadow/Building.json', 'r') as f:
    info = json.load(f)

extended_cells = grid.cells.reshape(-1, 5)[:, 1:].astype(np.uintp)
extended_verts = grid.points.astype(np.float64)

num_dom_cells = info['num_dom_cells']
num_dom_verts = info['num_dom_points']

verts = extended_verts[:num_dom_verts]
cells = extended_cells[:num_dom_cells]

domain = Domain(verts, cells, refl_coef=0.7)
mesh = domain.mesh
surf_mesh = mesh.get_surface_mesh()
h = mesh.min_edge_length/2

############################################################################
#

l_int = np.array([_ for _ in range(num_dom_verts) if not domain.mesh.bdv(_)])
src_index = l_int[0]
omega = 3000 # Hz

pt_src_field = PointSourceField(domain, src_index, omega, r=1.5)
pt_src_field.solve()

# refl_index = 3
# refl_field = ReflectedField(
#     domain, refl_index, omega,
#     *pt_src_field._get_reflection_BCs(mesh.get_reflector(refl_index)))
# refl_field.solve()

diff_index = 5
diff_field = DiffractedField(
    domain, diff_index, omega,
    *pt_src_field._get_diffraction_BCs(mesh.get_diffractor(diff_index)))
diff_field.solve()

############################################################################
# Set up different views

def set_view_1(P):
    P.camera.position = (0, 20, -35)
    P.camera.focal_point = (0, 0, 0)
    P.camera.up = (0, 1, 0)
    P.camera.zoom(1.4)

def set_view_2(P):
    P.camera.position = (-30, 40, -30)
    P.camera.focal_point = (0, 0, 0)
    P.camera.up = (0, 1, 0)
    P.camera.zoom(1.3)

def set_view_3(P):
    P.camera.position = (35, 20, 0)
    P.camera.focal_point = (0, 0, 0)
    P.camera.up = (0, 1, 0)
    P.camera.zoom(1.4)

def set_view_4(P):
    P.camera.position = (-1.5, 45, 0)
    P.camera.focal_point = (-1.5, 0, 0)
    # P.camera.up = (0, 1, 0)
    P.camera.zoom(1.07)

def set_view_5(P):
    P.camera.position = (41, 1.5, 0)
    P.camera.focal_point = (1, 1.5, 0)
    P.camera.up = (0, 1, 0)
    P.camera.zoom(1.9)

set_view = {1: set_view_1, 2: set_view_2, 3: set_view_3, 4: set_view_4,
            5: set_view_5}
window_size = {1: (1280, 720), 2: (1280, 1000), 3: (1280, 720), 4:
               (1950, 1650), 5: (1950, 440)}
show_scalar_bar = {1: False, 2: False, 3: False, 4: True, 5: False}
scalar_bar_args = {
    1: None,
    2: None,
    3: None,
    4: {'vertical': True,
        'shadow': True,
        'label_font_size': 48,
        'n_labels': 5,
        'position_x': 0.87,
        'position_y': 0.1,
        'width': 0.075,
        'height': 0.8,
        'shadow': True},
    5: None
}

############################################################################

transparent_background = False

def get_plot_image(view, F, field):
    P = pv.Plotter(off_screen=True, window_size=window_size[view])
    set_view[view](P)
    plot_mesh2(P, surf_mesh, color='gray', opacity=0.25)
    if F is not None:
        if isinstance(F, PointSourceField):
            plot_field_BCs(P, F, color='red', r=3*h)
        elif isinstance(F, ReflectedField):
            plot_field_BCs(P, pt_src_field, color='red', r=3*h)
            plot_field_BCs(P, F, color='cyan', show_edges=True, r=3*h,
                           opacity=0.5)
        elif isinstance(F, DiffractedField):
            plot_field_BCs(P, pt_src_field, color='red', r=3*h)
            plot_field_BCs(P, F, color='cyan', r=6*h, opacity=1)
    plot_scalar_field_kwargs = dict(
        show_scalar_bar=show_scalar_bar[view],
        scalar_bar_args=scalar_bar_args[view],
        render_points_as_spheres=False,
        point_size=10
    )
    try:
        time = F.time
    except:
        time = F.eik.T/type(F).speed_of_sound
    amp = abs(F.amplitude)
    amp[np.isnan(amp)] = 1
    amp = np.clip(amp, 1e-3, 1)
    loudness = 20*np.log10(amp)
    mask = (loudness > -60) & np.isfinite(time)
    if field == 'time':
        plot_scalar_field(P, verts[mask], time[mask], title='',
                          cmap=cc.cm.colorwheel,
                          **plot_scalar_field_kwargs)
    elif field == 'loudness':
        plot_scalar_field(P, verts[mask], loudness[mask],
                          title='', cmap=cc.cm.CET_L3, clim=(-60, 0),
                          **plot_scalar_field_kwargs)
    return P.screenshot(transparent_background=transparent_background)

mode = 'RGBA' if transparent_background else 'RGB'

view = 4

# F, field_index = pt_src_field, src_index
# F, field_index = refl_field, refl_index
F, field_index = diff_field, diff_index

field_type = {
    PointSourceField: 'P',
    ReflectedField: 'R',
    DiffractedField: 'D'
}[type(F)]

field = 'time'
img = Image.fromarray(get_plot_image(view, F, field), mode)
img.save('field_%s_%03d_%s.png' % (field_type, field_index, field))

field = 'loudness'
img = Image.fromarray(get_plot_image(view, F, field), mode)
img.save('field_%s_%03d_%s.png' % (field_type, field_index, field))
