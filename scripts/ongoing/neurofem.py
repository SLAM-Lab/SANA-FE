# Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS)
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution. Neither the name of the
# copyright holder nor the names of its contributors may be used to endorse or
# promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
from meshpy.triangle import MeshInfo, build, subdivide_facets, refine
from matplotlib.tri import Triangulation
from matplotlib.gridspec import GridSpec
import os
import tqdm
from itertools import product, combinations_with_replacement
import glob

import os
import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.abspath((os.path.join(SCRIPT_DIR, os.pardir, os.pardir)))
sys.path.insert(0, PROJECT_DIR)
import sanafe

def build_disk_mesh(max_V, n_bdry, plot_mesh=False, bdry_steiner=True):
    '''
    Constructs a triangulation of a 2D disk.
    max_V:  The maximum volume (area) of a triangle
    n_bdry: The number of points on the outer boundary of the disk

    Returns;
    tri_points:  Array of (x,y) positions of the vertices (mesh nodes)
    tri_tris:    Array of indices for vertices in each triangle
    tri_bdry:    Boolean array indicating whether a vertex is a boundary vertex.  tri_bdry[idx] = 1 if node <idx> is boundary
    '''
    bdry_pts = []
    bdry_fcs = []

    bdry_idxs = np.linspace(0, 2*np.pi, n_bdry, endpoint=False)
    bdry_pts.extend([(np.cos(x), np.sin(x)) for x in bdry_idxs])
    bdry_fcs.extend([sorted([x, ((x+1) % (n_bdry))]) for x in range(n_bdry)])

    n_bdry_tot = n_bdry

    mesh_info = MeshInfo()
    mesh_info.set_points(bdry_pts)
    mesh_info.set_facets(bdry_fcs)

    mesh = build(mesh_info, verbose=False, min_angle=30, volume_constraints=True, max_volume=max_V, allow_boundary_steiner=bdry_steiner)

    ts = Triangulation([x[0] for x in mesh.points], [x[1] for x in mesh.points], mesh.elements)

    if plot_mesh:
        fig, ax = plt.subplots(1, 1, figsize=(6, 6), dpi=100)
        ax.triplot(ts, color='k')
        fig.savefig('fem_mesh_3.pdf', transparent=True)

    tri_points = np.vstack([x for x in mesh.points])
    tri_tris = np.array([x for x in mesh.elements])
    tri_bdry = np.array([x for x in mesh.point_markers])

    N_boundary = sum(tri_bdry==1)

    N_interior = len(tri_points) - N_boundary

    return tri_points, tri_tris, tri_bdry, N_interior, N_boundary

def save_mesh(mesh_folder, tri_points, tri_tris, tri_bdry, mesh_tag, max_V, N_interior, N_boundary):
    '''
    Save the a triangulation to an npz file
    '''
    triangulation_file = os.path.join(mesh_folder, 'mesh_{}_nmesh-{}_maxv-{:.5}_nbdry-{}.npz'.format(mesh_tag, N_interior, max_V, n_bdry))
    np.savez(triangulation_file, tri_points=tri_points, tri_tris=tri_tris, tri_bdry=tri_bdry, nmesh=N_interior, maxv=max_V)
    return triangulation_file

def calculate_triangle_areas(tri_points, tri_tris, tri_bdry):
    '''
    Calculate the area of each triangle in a triangulation
    '''
    tri_areas = np.zeros(len(tri_tris))

    for idx, triangle in enumerate(tri_tris):
        # compute triangle area:
        v1 = tri_points[triangle[1], :] - tri_points[triangle[0], :]
        v2 = tri_points[triangle[2], :] - tri_points[triangle[0], :]
        A = np.linalg.det(np.vstack([v1, v2]))
        tri_areas[idx] = A

    return tri_areas

def generate_poisson_system_matrix_sparse(tri_points, tri_tris, tri_bdry, tri_areas):
    '''
    Build a sparse matrix representing the Laplace operator on given triangulation
    '''

    N_interior = sum(tri_bdry==0)
    N_boundary = sum(tri_bdry==1)

    A_sys = sp.sparse.lil_array((N_interior, N_interior), dtype=np.float64)

    # for each triangle
    for idx, triangle in enumerate(tri_tris):
        A = tri_areas[idx]
        vtx_vecs = []

        # for each vertex in the triangle, compute the associated vector
        for idx2, vtx in enumerate(triangle):

            # if the vertex is boundary, nothing to do (Dirichlet Boundary conditions)
            if tri_bdry[vtx] == 1:
                continue

            opp_idxs = list(triangle)
            opp_idxs.remove(vtx)
            b = tri_points[opp_idxs[1], :]
            a = tri_points[opp_idxs[0], :]
            vtx_vec = ((-1)**idx2)*np.array([-(b[1]-a[1]), b[0]-a[0]])/np.sqrt(A)
            vtx_vecs.append((vtx, vtx_vec))

        # now compute the dot product of each pair of vertex vectors in this triangle
        for (va, vb) in combinations_with_replacement(vtx_vecs, 2):
            idxa = va[0] - N_boundary
            idxb = vb[0] - N_boundary
            dp = np.dot(va[1], vb[1])
            A_sys[idxa, idxb] += dp
            if idxa != idxb:
                A_sys[idxb, idxa] += dp

    return sp.sparse.csr_array(10*A_sys)


def generate_poisson_mass_matrix_sparse(tri_points, tri_tris, tri_bdry, tri_areas):
    '''
    Generate the mass matrix for the poisson problem for the given triangulation.
    '''
    N_interior = sum(tri_bdry==0)
    N_boundary = sum(tri_bdry==1)
    M = sp.sparse.lil_array((N_interior, N_interior))

    for triangle, area in zip(tri_tris, tri_areas):
        # remove border vertices
        t1 = triangle.copy()
        t1 = [x for x in triangle if tri_bdry[x]==0]
        for vtx1 in t1:
            for vtx2 in t1:
                M[vtx1-N_boundary, vtx2-N_boundary] += area/12
                if vtx1 == vtx2:
                    M[vtx1-N_boundary, vtx2-N_boundary] += area/12
    return sp.sparse.csr_array(M)


def solve_system_sparse_direct(A_sys, M, f):
    '''
    Apply spsolve to solve the FEM problem conventionally
    '''
    soln = sp.sparse.linalg.spsolve(A_sys, M @ f)
    return soln

def generate_rhs_constant(N_interior, tri_points, tri_bdry):
    '''
    Generate a constant forcing function on the whole mesh
    '''
    return -200 * np.ones(N_interior)

def generate_rhs_offset(N_interior, tri_points, tri_bdry):
    '''
    Generate a non-constant forcing function on the mesh
    '''

    px = tri_points[tri_bdry==0, 0]
    py = tri_points[tri_bdry==0, 1]

    # the driving function
    f = -30*(px-0.25)**2 - 30*(py+0.13)**2 + 6
    return -20*f

def compute_const_analytic_soln(tri_points, tri_bdry):
    '''
    The poisson problem on a disk with dirichlet boundary conditions and constant forcing function has an analytic solution
    Evaluate the analytic solution on the given triangulation.
    '''

    analytic_soln = -5*(1 - tri_points[tri_bdry==0, 0]**2 - tri_points[tri_bdry==0, 1]**2)
    return analytic_soln

def plot_soln(tri_points, tri_tris, tri_bdry, soln, title=''):

    fullsoln = np.zeros(len(tri_points))
    fullsoln[tri_bdry==0] = soln

    from scipy.interpolate import LinearNDInterpolator

    interp = LinearNDInterpolator(tri_points, fullsoln)

    X = np.linspace(min(tri_points[:, 0]), max(tri_points[:, 0]), 128)
    Y = np.linspace(min(tri_points[:, 1]), max(tri_points[:, 1]), 128)
    X, Y = np.meshgrid(X, Y)  # 2D grid for interpolation
    Z = interp(X, Y)

    fig, ax = plt.subplots(1,1, figsize=(6, 6))
    plt.pcolormesh(X, Y, Z, shading='auto', cmap=plt.cm.Spectral)
    plt.triplot(tri_points[:,0], tri_points[:,1], tri_tris, alpha=0.75, c='k', linewidth=0.75)

    ax.set_title(title)
    ax.set_frame_on(False)
    ax.set_xticks([])
    ax.set_yticks([])


def generate_gamma_sparse(n_sys, neurons_per_mesh_point, gamma_norm):
    '''
    Generates the Gamma matrix mapping neurons to readout variables.
    n_sys: number of readout variables
    neurons_per_mesh_point: number of neurons to place at each mesh node
    gamma_norm:  Magnitude of the readout kernel for each neuron

    '''

    n_neurons = neurons_per_mesh_point * n_sys

    gamma_data = []
    gamma_row = []
    gamma_col = []
    for pt in range(n_sys):

        start_idx = pt * neurons_per_mesh_point
        half_npm = neurons_per_mesh_point//2

        gamma_row.extend(neurons_per_mesh_point*[pt])
        gamma_col.extend(list(range(start_idx, start_idx + neurons_per_mesh_point)))
        gamma_data.extend(half_npm*[-gamma_norm])
        gamma_data.extend(half_npm*[gamma_norm])

    return sp.sparse.csr_array((np.array(gamma_data), (np.array(gamma_row), np.array(gamma_col))), shape=(n_sys, n_neurons))


def create_spiking_fem_network_sparse(A_sys, gamma, lambda_d=10, lambda_v=20, mu=1e-6, nu=1e-5, tau_A=0.1):
    '''
    Generates sparse weight matrices, thresholds, and reset voltages for a NeuroFEM network implementing the
    linear system A_sys with readout kernel gamma
    '''
    n_sys, n_neurons = gamma.shape
    gTg = (gamma.T) @ gamma
    omega_f = gTg + mu*(lambda_d**2) * sp.sparse.eye(n_neurons)

    GTAG = (gamma.T) @ (A_sys / tau_A) @ gamma  # should be sparse

    threshs = 0.5 * (nu*lambda_d**2 + (gamma_norm**2) * np.ones(n_neurons))
    v_reset = threshs - (gamma_norm**2 + mu*(lambda_d**2)) * np.ones(n_neurons)
    omega_f = omega_f - (gamma_norm**2 + mu*(lambda_d**2)) * sp.sparse.eye(n_neurons)

    return omega_f, GTAG, threshs, v_reset

def save_float_snn_params(model_folder, GTAG, omega_f, lambda_d, lambda_v,
                          threshs, v_reset, dt, gamma, c_in, sigma_v,
                          soln, soln_2, A_sys,
                          mesh_tag, gamma_pow2, N_interior, N_boundary, neurons_per_mesh_point, max_V):
    # Save the model
    float_model_file = os.path.join(model_folder,'fem_network_float_params_{}_gamma-{}_NMESH-{}_NPM-{}_maxV-{}.npz'.format(mesh_tag, gamma_pow2, N_interior, neurons_per_mesh_point, max_V, N_boundary))
    np.savez_compressed(float_model_file, GTAG=GTAG, omega_f=omega_f,
                       lambda_d = lambda_d, lambda_v=lambda_v,
                       thresh=threshs, v_reset=v_reset, dt=dt, gamma=gamma, c_in=c_in, sigma_v=sigma_v,
                       soln_f=soln, soln_f2=soln_2, A_sys=-A_sys)


#### The power of two defining the magnitude of gamma
gamma_pow2 = -6 # In the paper, we used gamma = 2^-6 or gamma = 2^-8
s_gamma = 7 - gamma_pow2
gamma_norm = (2**gamma_pow2 - 2**-s_gamma)

#### The power of two to which we rescale the slow weight matrix
omega_s_pow2 = 1
omega_max = 2**omega_s_pow2 - 2**-(7-omega_s_pow2)


#### network parameters
lambda_d = 8              # Hz; slow variable time constant.  Default = 10 Hz
lambda_v = 16             # Hz; membrane potential time constant. Default = 20 Hz

ki = 16                   # integral control gain
kp = 4                    # proportional control gain


#### simulation parameters
dt = 2**(-12)             # power of two
#n_timesteps = 50000
n_timesteps = 5000


#### place to store the model parameters
model_folder='./runs/neurofem'
mesh_folder = model_folder
mesh_tag = 'runs-example'


##### create the mesh resolutions
# for the paper, we chose max_Vs that gave between ~100 and ~1300 total triangles
# this gave between ~100 and ~10000 mesh nodes

# This generates big files so we leave it off in the example
# max_tris = 1300
# min_tris = 100
# circ_area = np.pi

# tri_sizes = np.linspace(min_tris, max_tris, n_h_steps)
# max_vs = circ_area/tri_sizes

# Instead, we set a single max_V

max_vs = [0.002]

##### set neurons per mesh point
# we tested either 8 or 16 neurons per mesh point

npm_list = np.array([8, 16])

# Whether to apply jacobi preconditioning to the system matrix before generating weight matrices
jacobi_precondition = True #True


for max_V in max_vs:
    if max_V <= 0.0025:
        n_bdry = 128
    else:
        n_bdry = 64

    # build mesh, generate linear system and right hand sides
    tri_points, tri_tris, tri_bdry, N_interior, N_boundary = build_disk_mesh(max_V, n_bdry, plot_mesh=False, bdry_steiner=False)
    tri_areas = calculate_triangle_areas(tri_points, tri_tris, tri_bdry)
    A_sys = generate_poisson_system_matrix_sparse(tri_points, tri_tris, tri_bdry, tri_areas)
    M = generate_poisson_mass_matrix_sparse(tri_points, tri_tris, tri_bdry, tri_areas)
    f1 = generate_rhs_constant(N_interior, tri_points, tri_bdry)
    f2 = generate_rhs_offset(N_interior, tri_points, tri_bdry)

    mesh_file = save_mesh(mesh_folder, tri_points, tri_tris, tri_bdry, mesh_tag, max_V, N_interior, N_boundary)

    if jacobi_precondition:
        Pinv = sp.sparse.diags(1/A_sys.diagonal())
        Pinv_half = sp.sparse.diags(np.sqrt(1/A_sys.diagonal()))
        A_sys =  Pinv @ A_sys
        M = Pinv @ M

    soln_f1 = solve_system_sparse_direct(A_sys, M, f1)
    soln_f2 = solve_system_sparse_direct(A_sys, M, f2)
    soln_const_analytic = compute_const_analytic_soln(tri_points, tri_bdry)

    for neurons_per_mesh_point in npm_list:

        n_neurons = neurons_per_mesh_point * N_interior
        sigma_v = 0.00225 # We kept this constant for all meshes
        tau_A = (gamma_norm**2) * np.amax(np.abs(A_sys))/omega_max
        mu = 0  # spike L2 norm penalization; set to 0
        nu = 0  # spike L1 norm penalization; set to 0

        # create system
        gamma = generate_gamma_sparse(N_interior, neurons_per_mesh_point, gamma_norm)
        omega_f, GTAG, threshs, v_reset = create_spiking_fem_network_sparse(-A_sys, gamma, lambda_d, lambda_v, mu, nu, tau_A)


        # neuron biases
        bias_f1 = gamma.T @ (M @ f1 / tau_A)
        bias_f2 = gamma.T @ (M @ f2 / tau_A)

        # Save the model
        float_model_file = os.path.join(model_folder,'fem_network_float_params_{}_gamma-{}_NMESH-{}_NPM-{}_maxV-{:.4}.npz'.format(mesh_tag, gamma_pow2, N_interior, neurons_per_mesh_point, max_V))
        np.savez_compressed(float_model_file, GTAG=GTAG, omega_f=omega_f,
                           lambda_d = lambda_d, lambda_v=lambda_v,
                           thresh=threshs, v_reset=v_reset, dt=dt, gamma=gamma, bias_f1=bias_f1, bias_f2=bias_f2, f1=f1, f2=f2, M=M, Pinv_half=Pinv_half,
                           sigma_v=sigma_v, soln_f1=soln_f1, soln_f2=soln_f2, soln_analytic=soln_const_analytic, A_sys=-A_sys, max_V=max_V, N_interior=N_interior, npm=neurons_per_mesh_point)
        print(float_model_file)


def convert_model(params_file, model_folder):
    with np.load(params_file, allow_pickle=True) as f:
        gtAg = f['GTAG'][()] #.toarray()
        omega_f = f['omega_f'][()] #.toarray()
        lambda_d = f['lambda_d']
        dt = f['dt']
        threshs = f['thresh']
        v_reset = f['v_reset']
        bias_f1 = f['bias_f1']
        bias_f2 = f['bias_f2']
        sigma_v = f['sigma_v']
        gamma = f['gamma'][()]
        nmesh = f['N_interior']
        npm = f['npm']
        max_V = f['max_V']
        soln_analytic = f['soln_analytic']

    n_sys, n_neurons = gamma.shape
    arch = sanafe.load_arch("arch/neurofem.yaml")  # TODO: create custom arch for this experiment
    snn = sanafe.Network()

    neuron_group = snn.create_neuron_group("solver", n_neurons,
        default_dendrite_hw_name="loihi_neurofem",
        soma_hw_name="loihi_neurofem",
        log_spikes=True,
        model_attributes={
            "lambda_d": lambda_d,
            "lambda_v": lambda_v,
            "ki": ki,
            "kp": kp,
            "dt": dt,
            "sigma_v": sigma_v,
            "reset": v_reset[0],
            "threshold": threshs[0],
        }
    )  # All neurons in single group for now

    print(f"Connecting {n_neurons} neurons")

    gtAg_coo = gtAg.tocoo()
    for j, i, weight_u1 in zip(gtAg_coo.row, gtAg_coo.col, gtAg_coo.data):
        weight_u2 = omega_f[j, i]
        if weight_u1 != 0:
            neuron_group.neurons[int(i)].connect_to_neuron(
                neuron_group.neurons[int(j)],
                {"weight": weight_u1, "compartment": 0})
        if weight_u2 != 0:
            neuron_group.neurons[int(i)].connect_to_neuron(
                neuron_group.neurons[int(j)],
                {"weight": weight_u2, "compartment": 1})

    print("neurons connected")
    # TODO: hardware mapping to cores? figure out actual mapping algorithm
    #  For now just map 1024 neurons per core
    for i, n in enumerate(neuron_group.neurons):
        core = i // 1024
        tile = core // 4
        offset = core % 4
        n.map_to_core(arch.tiles[tile].cores[offset])

    chip = sanafe.SpikingChip(arch)
    chip.load(snn)
    snn.save("runs/neurofem/neurofem.yaml")


    del arch
    del snn

    mapped_neurons = chip.mapped_neuron_groups["solver"]
    for i, n in enumerate(mapped_neurons):
        n.set_attributes(model_attributes={"bias": bias_f1[i]})
    results = chip.sim(n_timesteps // 2, spike_trace=True, processing_threads=16, timing_model="simple")
                       #scheduler_threads=8)
    print("Simulation finished")
    #print(results)

    print("Converting spike format")
    spikes = np.zeros((int(n_neurons), int(n_timesteps)))
    spike_trace = results["spike_trace"]
    for timestep, spike_data in enumerate(spike_trace):
        for spike_address in spike_data:
            spikes[spike_address.neuron_offset, timestep] = 1

    print("Setting up next run")
    for i, n in enumerate(mapped_neurons):
        n.set_attributes(model_attributes={"bias": bias_f2[i]})
    print("Running")

    results = chip.sim((n_timesteps+1) // 2, spike_trace=True, processing_threads=16, timing_model="simple")
    spike_trace = results["spike_trace"]
    for timestep, spike_data in enumerate(spike_trace):
        for spike_address in spike_data:
            spikes[spike_address.neuron_offset, (n_timesteps//2) + timestep] = 1
    #print(results)

    # Calculate output value for the mesh (performed by CPU)
    output = np.zeros((n_sys, n_timesteps))
    for step in range(1, n_timesteps):
        output[:, step] = (1.0-dt*lambda_d)*output[:, step-1] + gamma.dot(spikes[:, step])

    print(output)
    output_file = os.path.join(model_folder, 'sanafe_results_nmesh-{}_npm-{}_nstep-{}.npz'.format(nmesh, npm, n_timesteps))
    np.savez_compressed(output_file, output=output, soln_f1=soln_f1, soln_f2=soln_f2, soln_analytic=soln_analytic, nmesh=nmesh, npm=npm, max_V=max_V, A_sys=A_sys, M=M, f1=f1)

    return output_file


def run_model(params_file, model_folder):
    output_file = ''
    print(params_file)

    # load model parameters
    with np.load(params_file, allow_pickle=True) as f:

        gtAg = f['GTAG'][()] #.toarray()
        omega_f = f['omega_f'][()] #.toarray()
        lambda_d = f['lambda_d']
        lambda_v = f['lambda_v']
        threshs = f['thresh']
        v_reset = f['v_reset']
        dt = f['dt']
        gamma = f['gamma'][()]
        bias_f1 = f['bias_f1']
        bias_f2 = f['bias_f2']
        sigma_v = f['sigma_v']
        soln_f1 = f['soln_f1']
        soln_f2 = f['soln_f2']
        soln_analytic = f['soln_analytic']
        nmesh = f['N_interior']
        npm = f['npm']
        max_V = f['max_V']

        M = f['M'][()]
        A_sys = f['A_sys'][()]
        f1 = f['f1']

    output_file = os.path.join(model_folder, 'neurofem_results_nmesh-{}_npm-{}_nstep-{}.npz'.format(nmesh, npm, n_timesteps))
    # check if output file exists
    if os.path.exists(output_file):
        return output_file

    n_sys, n_neurons = gamma.shape

    V = np.zeros(n_neurons)
    spikes_save = sp.sparse.lil_array((n_neurons, n_timesteps), dtype=np.int8)
    spikes = np.zeros((n_neurons,1), dtype=np.int8)
    output = np.zeros((n_sys, n_timesteps))

    u_int = np.zeros(n_neurons)
    u_err = np.zeros(n_neurons)

    u1 = np.zeros(n_neurons)

    u2 = np.zeros(n_neurons)

    c_in_fixed_gamma = np.copy(bias_f1)

    for step in tqdm.tqdm(range(1, n_timesteps)):

        if step == n_timesteps//2:
            c_in_fixed_gamma = np.copy(bias_f2)

        ######################
        # compute update to u1
        du1 = dt*(-lambda_d*u1) + gtAg @ spikes[:, 0]
        u1 += du1
        ######################

        ######################
        # comptue update to u2
        du2 = dt*(-lambda_d*u2) + lambda_d * omega_f @ spikes[:, 0]
        u2 += du2
        ######################

        #############
        # compute err
        u_err = u1 + c_in_fixed_gamma
        #############

        ##############
        # update u_int
        u_int += dt*u_err
        ##############

        ##########
        # update v
        dv = dt*(-lambda_v*V + kp*u_err + ki*u_int + u2) - omega_f.dot(spikes[:, 0]) + sigma_v*np.random.randn(n_neurons)
        V += dv
        ##########

        spikes[:, 0] = np.greater(V, threshs).astype(np.int8)
        # saving all the spikes makes for large files so we'll skip it for now
        #spikes_save[:, step] = spikes

        V[spikes[:, 0]>0] = v_reset[spikes[:, 0]>0]

        output[:, step] = (1-dt*lambda_d)*output[:, step-1] + gamma.dot(spikes[:, 0])

    # save output
    np.savez_compressed(output_file, output=output, soln_f1=soln_f1, soln_f2=soln_f2, soln_analytic=soln_analytic, nmesh=nmesh, npm=npm, max_V=max_V, A_sys=A_sys, M=M, f1=f1)
    return output_file

float_params_list = []
mesh_list = []


float_params_list = glob.glob(os.path.join(model_folder, 'fem_network_float_params*'))
float_params_list.reverse()

print(float_params_list)

output_file = convert_model(float_params_list[0], model_folder)
#output_file = run_model(float_params_list[0], model_folder)

# load the results

with np.load(output_file, allow_pickle=True) as f:
    output = f['output']
    soln_f1 = f['soln_f1']
    soln_f2 = f['soln_f2']
    soln_analytic = f['soln_analytic']

    A_sys = -f['A_sys'][()]
    M = f['M'][()]
    f1 = f['f1']

fig, ax = plt.subplots(1, 1, figsize=(6, 2.5), tight_layout=True, dpi=300)

rel_resid_1 = np.linalg.norm(output - soln_f1[..., None], axis=0)/np.linalg.norm(soln_f1)
rel_resid_2 = np.linalg.norm(output - soln_f2[..., None], axis=0)/np.linalg.norm(soln_f2)


v1_trace, = ax.plot(rel_resid_1, color='tab:orange', linewidth=2, label='$f_1$')
v2_trace, = ax.plot(rel_resid_2, color='tab:blue', linewidth=2, label='$f_2$')
ax.set_xlabel('Timestep', fontsize=14)
ax.set_ylabel('Relative Error', fontsize=14)
ax.axhline(0, 0, 1, color='k', linestyle='--', linewidth=1, alpha=0.3)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# We average the readout over the last 10000 timesteps
# Compute the relative residuals/errors for the constant forcing term.

n_timesteps_2 = n_timesteps//2
#neurofem_res = np.mean(output[:, n_timesteps_2-10000:n_timesteps_2], axis=1)
neurofem_res = np.mean(output[:, n_timesteps_2-1000:n_timesteps_2], axis=1)

# The relative residual of the linear system
rel_resid = np.linalg.norm(M@f1 - A_sys @ neurofem_res)/np.linalg.norm(M@f1)

# relative error between neurofem and analytic
rel_error_neuro_analytic = np.linalg.norm(soln_analytic - neurofem_res)/np.linalg.norm(soln_analytic)

# relative error between neurofem and solver
rel_error_neuro_solver = np.linalg.norm(soln_f1 - neurofem_res)/np.linalg.norm(soln_f1)

print('Relative residual: {:.6f}'.format(rel_resid))
print('Relative error (NeuroFEM/Analytic): {:.6f}'.format(rel_error_neuro_analytic))
print('Relative error (NeuroFEM/Solver): {:.6f}'.format(rel_error_neuro_solver))

# Plot the neurofem solution against the solver solution

fig, ax = plt.subplots(1,1, figsize=(12, 3))
ax.plot(soln_f1, label='Ground Truth', linewidth=0.5)
ax.plot(neurofem_res, label='Spiking Network', linewidth=0.5)
ax.legend()
ax.set_xlabel('Mesh Node Index', fontsize=14)
ax.set_ylabel('Value', fontsize=14)

# Plot the difference between neurofem and solver

print(neurofem_res)
print(soln_f1)
print(neurofem_res - soln_f1)
#exit()

fig, ax = plt.subplots(1,1, figsize=(12, 3))
ax.plot(neurofem_res - soln_f1, label='NeuroFEM - Solver', linewidth=0.5)
ax.axhline(0, 0, 1, color='k', linestyle='--', linewidth=1, alpha=0.3)
ax.legend()
ax.set_xlabel('Mesh Node Index', fontsize=14)
ax.set_ylabel('Difference', fontsize=14)

# load mesh
with np.load(mesh_file, allow_pickle=True) as f:
    tri_points = f['tri_points']
    tri_tris = f['tri_tris']
    tri_bdry = f['tri_bdry']

neurofem_res_full = np.zeros(len(tri_points))
neurofem_res_full[tri_bdry==0] = neurofem_res

solver_full = np.zeros(len(tri_points))
solver_full[tri_bdry==0] = soln_f1


# SNN 3D plot
fig = plt.figure(figsize=(10, 10), tight_layout=True)

gs = GridSpec(1, 1, figure=fig)
ax_snn = fig.add_subplot(gs[0, 0], projection='3d')

cmap='seismic'
srf = ax_snn.plot_trisurf(tri_points[:, 0], tri_points[:, 1], tri_tris, Z=-neurofem_res_full, linewidth=0.15, antialiased=False, shade=True, cmap='turbo', edgecolor='k',
                           vmin=-5, vmax=5)


ax_snn.set_zlim([-4, 4])
for ax in [ax_snn]:
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_zticklabels([])

ax_snn.set_title('NeuroFEM Solution', fontsize=24)

# SNN 3D plot
fig = plt.figure(figsize=(10, 10), tight_layout=True)

gs = GridSpec(1, 1, figure=fig)
ax_snn = fig.add_subplot(gs[0, 0], projection='3d')

cmap='seismic'
srf = ax_snn.plot_trisurf(tri_points[:, 0], tri_points[:, 1], tri_tris, Z=-neurofem_res_full + solver_full, linewidth=0.15, antialiased=False, shade=True, cmap='turbo', edgecolor='k',
                           vmin=-0.0005, vmax=0.0005)


ax_snn.set_zlim([-0.001, 0.001])
ax_snn.set_title('NeuroFEM - Solver', fontsize=24)
ax_snn.set_zlabel('Difference', fontsize=16)

plt.show()
