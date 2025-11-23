import argparse
import imageio
import numpy as np
import pyvista as pv


parser = argparse.ArgumentParser(
    prog="pp-gif",
    description="A simple program to make a gif from exodus files"
)

parser.add_argument("--exodus-file")
parser.add_argument("--fps", default=1)
parser.add_argument("--gif-file")
parser.add_argument("--view-xy", action="store_true")
args = parser.parse_args()

class GifGenerator:
    def __init__(
        self, 
        exo_file
    ):
        self.reader = pv.get_reader(exo_file)
        self.plotter = pv.Plotter(off_screen=True)

        time_values = self.reader.time_values
        if not time_values:
            # mesh = self.reader.read()
            time_values = [0.]
        
        self.time_values = time_values

    def generate_deforming_mesh_gif(
        self,         
        file_name, 
        displacement_var_name="displ_",
        fps=1,
        view_xy=False
    ):
        n_steps = len(self.time_values)
        print(f"Found {n_steps} time steps")

        # quick check: ensure the nodal variable exists for at least one step
        # read the first time step to inspect arrays
        # self.reader.set_active_time_value(self.time_values[0])
        # mesh0 = self.reader.read()
        # print("Available point (nodal) arrays:", list(mesh0[0][0].point_data.keys()))
        # if var_name not in mesh0[0][0].point_data:
        #     raise ValueError(f"Nodal variable '{var_name}' not found. Available: {list(mesh0.point_data.keys())}")

        # # calculating c lims
        # global_min = np.inf
        # global_max = -np.inf

        # for t in self.time_values:
        #     self.reader.set_active_time_value(t)
        #     mesh = self.reader.read()
        #     data = mesh[0][0].point_data[var_name]

        #     if data is None:
        #         continue

        #     global_min = min(global_min, float(np.nanmin(data)))
        #     global_max = max(global_max, float(np.nanmax(data)))

        #     if not np.isfinite(global_min) or not np.isfinite(global_max):
        #         clim = None
        #         print("Could not compute finite global min/max; using autoscaling.")
        #     else:
        #         clim = (global_min, global_max)

        frames = []
        
        # if component is not None:
        #     var_name = var_name + component
        self.reader.set_active_time_value(self.time_values[0])
        mesh = self.reader.read()
        points = mesh[0][0].points

        for i, t in enumerate(self.time_values):
            print(f"Rendering step {i + 1}/{n_steps} (time value = {t})")
            self.reader.set_active_time_value(t)
            mesh = self.reader.read()
            displ = mesh[0][0].point_data[displacement_var_name]
            # mesh[0][0].points = mesh[0][0].points + displ
            mesh[0][0].points = points + displ
            self.plotter.clear()  # remove previous actors
            actor = self.plotter.add_mesh(
                mesh,
                color="orange",
                show_edges=True
            )

            if view_xy:
                self.plotter.view_xy()

            img = self.plotter.screenshot(return_img=True)
            frames.append(img)

        imageio.mimsave(file_name, frames, fps=fps)


gif_gen = GifGenerator(args.exodus_file)
gif_gen.generate_deforming_mesh_gif(
    args.gif_file,
    fps=int(args.fps),
    view_xy=args.view_xy
)