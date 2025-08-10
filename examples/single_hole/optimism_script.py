from jax import grad
from jax import jit
from jax import value_and_grad
import optimism
from optimism import EquationSolver
# from plato_optimism import exodus_writer as ExodusWriter
# from plato_optimism import adjoint_problem_function_space as AdjointFunctionSpace
from optimism import FunctionSpace
from optimism import Interpolants
from optimism import Mechanics
from optimism import Mesh
from optimism import Objective
from optimism import QuadratureRule
from optimism import ReadExodusMesh
from optimism import SparseMatrixAssembler
from optimism import VTKWriter
from optimism.FunctionSpace import DofManager
from optimism.FunctionSpace import EssentialBC
from optimism.inverse import AdjointFunctionSpace
from optimism.material import Neohookean
from typing import Optional

import matplotlib.pyplot as plt

import jax.numpy as np
import numpy as onp

class NodalCoordinateOptimization:

    def __init__(self):
        self.writeOutput = True

        # Target forces for strains of [1.5, 3.0, 4.5, 6.0]
        self.targetSteps = [5, 10, 15, 20] 
        # self.targetForces = [0.03140434, 0.05769101, 0.07248617, 0.07583194] # actual forces
        self.targetForces = [0.035, 0.045, 0.075, 0.076] # artificially perturbed

        self.stateNotStored = True
        self.state = []

        self.quad_rule = QuadratureRule.create_quadrature_rule_on_triangle(degree=2)

        self.ebcs = [
            EssentialBC(nodeSet='yminus_sideset', component=0),
            EssentialBC(nodeSet='yminus_sideset', component=1),
            EssentialBC(nodeSet='yplus_sideset', component=0),
            EssentialBC(nodeSet='yplus_sideset', component=1)
        ]

        props = {
            'elastic modulus': 2. * 49.0e-3 * (1. + 0.48),
            'poisson ratio': 0.48,
            'version': 'adagio'
        }
        self.mat_model = Neohookean.create_material_model_functions(props)
        self.props = Neohookean.create_material_properties(props)
        print(self.props)
        self.eq_settings = EquationSolver.get_settings(
            use_incremental_objective=False,
            max_trust_iters=100,
            tr_size=0.25,
            min_tr_size=1e-15,
            tol=5e-8
        )

        # self.input_mesh = './single_hole.exo'
        self.input_mesh = './single_hole_temp_temp.exo'
        if self.writeOutput:
          self.output_file = 'output.exo'

        self.plot_file = 'disp_control_response.npz'
        self.steps = 20
        self.maxDisp = 1.5

    def reload_mesh(self):
        origMesh = ReadExodusMesh.read_exodus_mesh(self.input_mesh)
        nodeSets = Mesh.create_nodesets_from_sidesets(origMesh)
        self.mesh = Mesh.mesh_with_nodesets(origMesh, nodeSets)
        self.stateNotStored = True

    def run_simulation(self):

        coords = self.mesh.coords

        # setup
        shapeOnRef = Interpolants.compute_shapes(self.mesh.parentElement, self.quad_rule.xigauss)
        func_space = AdjointFunctionSpace.construct_function_space_for_adjoint(coords, shapeOnRef, self.mesh, self.quad_rule)
        mech_funcs = Mechanics.create_mechanics_functions(func_space, mode2D='plane strain', materialModel=self.mat_model)
        dof_manager = DofManager(func_space, 2, self.ebcs)

        # methods defined on the fly
        def get_ubcs(p):
            disp = p[0]
            V = np.zeros(coords.shape)
            index = (self.mesh.nodeSets['yplus_sideset'], 1)
            V = V.at[index].set(disp)
            return dof_manager.get_bc_values(V)

        def create_field(Uu, p):
            return dof_manager.create_field(Uu, get_ubcs(p))

        def energy_function(Uu, p):
            U = create_field(Uu, p)
            internal_variables = p[1]
            return mech_funcs.compute_strain_energy(U, internal_variables, self.props)
        
        def energy_function_alt(U, p):
            internal_variables = p[1]
            return mech_funcs.compute_strain_energy(U, internal_variables, self.props)

        nodal_forces = jit(grad(energy_function_alt, argnums=0))

        def assemble_sparse(Uu, p):
            U = create_field(Uu, p)
            internal_variables = p[1]
            element_stiffnesses = mech_funcs.compute_element_stiffnesses(U, internal_variables, self.props)
            return SparseMatrixAssembler.\
                assemble_sparse_stiffness_matrix(element_stiffnesses, func_space.mesh.conns, dof_manager)
    
        def store_force_displacement(Uu, p, dispval, force, disp):
            U = create_field(Uu, p)
            f = nodal_forces(U, p)

            index = (self.mesh.nodeSets['yplus_sideset'], 1)
            force.append( onp.abs(onp.sum(onp.array(f.at[index].get()))) )

            disp.append( onp.abs(dispval) )

            with open(self.plot_file,'wb') as f:
                np.savez(f, force=force, displacement=disp)

        # only call after calculations are finished
        def plot_solution(plotName, mesh, Uu, p):
            U = create_field(Uu, p)
            writer = VTKWriter.VTKWriter(mesh, baseFileName=plotName)
            writer.add_nodal_field(
                name='displacement',
                nodalData=U,
                fieldType=VTKWriter.VTKFieldType.VECTORS
            )
            writer.write()

        # problem set up
        Uu = dof_manager.get_unknown_values(np.zeros(coords.shape))
        ivs = mech_funcs.compute_initial_state()
        p = Objective.Params(0., ivs)
        precond_strategy = Objective.PrecondStrategy(assemble_sparse)
        objective = Objective.Objective(energy_function, Uu, p, precond_strategy)

        # set up output mesh
        if self.writeOutput:
            plot_solution(f"output_0.vtk", self.mesh, Uu, p)

        # loop over load steps
        disp = 0.
        fd_force = []
        fd_disp = []

        store_force_displacement(Uu, p, disp, fd_force, fd_disp)
        self.state.append((Uu, p))

        for step in range(1, self.steps+1):

            print('--------------------------------------')
            print('LOAD STEP ', step)
            disp = disp - self.maxDisp / self.steps
            p = Objective.param_index_update(p, 0, disp)
            Uu, _ = EquationSolver.nonlinear_equation_solve(objective, Uu, p, self.eq_settings)

            store_force_displacement(Uu, p, disp, fd_force, fd_disp)
            self.state.append((Uu, p))

            if self.writeOutput:
                plot_solution(f"output_{step}.vtk", self.mesh, Uu, p)

        self.stateNotStored = False

    def objective_function(self, coords):
        shapeOnRef = Interpolants.compute_shapes(self.mesh.parentElement, self.quad_rule.xigauss)
        f_space = AdjointFunctionSpace.construct_function_space_for_adjoint(coords, shapeOnRef, self.mesh, self.quad_rule)
        m_funcs = Mechanics.create_mechanics_functions(f_space, mode2D='plane strain', materialModel=self.mat_model)

        dof_manager = DofManager(f_space, 2, self.ebcs)
        def get_ubcs(p):
            disp = p[0]
            V = np.zeros(coords.shape)
            index = (self.mesh.nodeSets['yplus_sideset'], 1)
            V = V.at[index].set(disp)
            return dof_manager.get_bc_values(V)

        endState = self.state[-1]
        U = dof_manager.create_field(endState[0], get_ubcs(endState[1]))

        return m_funcs.compute_strain_energy(U, endState[1][1], self.props)
    
    def get_objective(self):
        if self.stateNotStored:
            self.run_simulation()

        value = -self.objective_function(self.mesh.coords) 
        return onp.array(value).item()        

    def get_gradient(self):
        if self.stateNotStored:
            self.run_simulation()

        gradient = -grad(self.objective_function, argnums=0)(self.mesh.coords)
        return onp.array(gradient)#.flatten().tolist()


if __name__ == '__main__':
    nco = NodalCoordinateOptimization()
    nco.reload_mesh()
    val = nco.get_objective()
    grad = nco.get_gradient()


    print("\n Objective is: ")
    print(val)
    print(grad)

    writer = VTKWriter.VTKWriter(nco.mesh, baseFileName="dfdX_output")
    writer.add_nodal_field(
        name='dfdX',
        nodalData=grad,
        fieldType=VTKWriter.VTKFieldType.VECTORS
    )
    writer.write()

