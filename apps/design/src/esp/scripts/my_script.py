# NOTE this is only hooked up for a single AIM currently

import pyCAPS
import sys

input_file = sys.argv[1]
out_level = True
problem = pyCAPS.Problem(
    problemName = "ESP_Mesh",
    capsFile=input_file,
    outLevel=0
)
aflr2 = problem.analysis.create(aim='aflr2AIM', name='aflr2')
aflr2.input.Mesh_Quiet_Flag = True
aflr2.input.Tess_Params = [problem.geometry.outpmtr.MeshLength, 1.0, 20.0]
aflr2.input.Mesh_Sizing = None

plato = problem.analysis.create(aim='platoAIM', name='plato')
plato.input["Mesh"].link(aflr2.output["Area_Mesh"])
plato.input.Mesh_Morph = False

plato.preAnalysis()
plato.postAnalysis()
