import pyvista as pv
import numpy as np

def VolumeViewer(V, name='Volume'):

	data1 = pv.wrap(np.asarray(V))
	opacity = [0, 0, 0, 0.1, 0.3, 0.6, 1]
	pl = pv.Plotter(lighting='none')
	#pl.add_volume(data1, cmap="bone")
	pl.add_mesh(data1.contour(isosurfaces = 3),color="white")
	light = pv.Light()
	light.set_camera_light()
	#light.set_direction_angle(20, -20)
	pl.add_light(light)
	pl.add_axes()
	pl.add_title(name, font='courier', color='k',font_size=20)
	pl.link_views()
	pl.show()
