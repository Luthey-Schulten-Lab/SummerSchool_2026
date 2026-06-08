import os
import subprocess
import numpy as np

base_dir = os.path.dirname(os.path.abspath(__file__))
dna_monomers_dir = os.path.join(base_dir, "../data/coords/")
output_dir = os.path.join(base_dir, "../data/")
template_dir = os.path.join(base_dir, "./")
btree_chromo_executable = "/Software/btree_chromo/build/apps/btree_chromo"

def create_directives(timestep):
    """
    Creates a directives file for btree_chromo by reading a template
    and replacing placeholders based on the timestep.
    """
    sim_prng_seed = 34
    run_name="summerschool"

    dna_monomers_path = f"{dna_monomers_dir}dna_{run_name}_{timestep}.bin"
    ribos_path = f"{dna_monomers_dir}ribo_{run_name}_{timestep}.bin"
    template_file = template_dir + "template.inp"
    fudge_factor = 160
    sphere_radius = 2000 + fudge_factor # 2160
    sphere_radius += timestep*10 # increase by 10 each minute, from 2000 to 2750
    load_boundary = "spherical_bdry:" + str(sphere_radius) + ".0, 0.0, 0.0, 0.0"
    run_dynamics = "simulator_run_soft_harmonic:20000,1000,10000,append,nofirst"

    if timestep>=60 :
        sphere_radius = 2160 + 600*(timestep - 90)**2/30**2 # from 2760 to 2160
        sphere_height = 30 + np.sqrt(timestep - 60)*2230/np.sqrt(30)  # from 30 to 2260
        
        # Read in monomer coordinates
        with open(dna_monomers_path,'rb') as f:
            DNAbin = np.fromfile(f,dtype=np.float64,count=-1)
        DNAcoords = DNAbin.reshape((3,DNAbin.shape[0]//3),order='F').T
        L = DNAcoords[:54338]
        R = DNAcoords[54338:]
        # Calculate center of mass of both chromosomes
        LCom = np.average(L, axis=0)
        RCom = np.average(R, axis=0)
        # Calculate vector between centers of mass
        n = LCom - RCom
        # Normalize the vector
        n_normalized = n / np.linalg.norm(n)
        # Convert to string with 3 decimal places
        n_string = ', '.join(f"{x:.3f}" for x in n_normalized)
        load_boundary = "overlapping_spheres_bdry:" + str(sphere_height) + ".0, " + str(sphere_radius) + ".0, 0.0, 0.0, 0.0, " + n_string
        BD_steps = 20000 + (timestep-60)*180000/30  # from 20000 to 200000
        run_dynamics = "simulator_run_soft_harmonic:"+str(BD_steps)+",1000,"+str(BD_steps//2)+",append,nofirst"

    load_loops = "load_loops:"+f"{output_dir}loops/loops_{run_name}_{timestep}.txt"
    equilibrate_loops = ""
    input_state = f"input_state:{output_dir}rep_states/rep_state_{run_name}_{timestep}.txt"

    if timestep == 0:
        load_loops = ""
        equilibrate_loops = "translocate:7200,F"
        append_string = "noappend,first"
        input_state = f"input_state:{template_dir}rep_state_initial.txt"
    else:
        append_string = "append,nofirst"

    # Read the template file and replace placeholders
    with open(template_file, 'r') as file:
        directives_content = file.read()

    directives_content = directives_content.format(
        sim_prng_seed=sim_prng_seed,
        base_dir=base_dir,
        run_name=run_name,
        input_state=input_state,
        dna_monomers_path=dna_monomers_path,
        ribos_path=ribos_path,
        output_dir=output_dir,
        timestep=timestep,
        next_timestep=timestep+1,
        dna_monomers_dir=dna_monomers_dir,
        load_boundary=load_boundary,
        load_loops=load_loops,
        equilibrate_loops=equilibrate_loops,
        append_string=append_string,
        run_dynamics=run_dynamics
    )

    # Write to a directives file for this timestep
    directives_filename = f"{output_dir}/btree_chromo_directives_{timestep}_{run_name}.inp"
    with open(directives_filename, 'w') as f:
        f.write(directives_content)

    return directives_filename



def run_btree_chromo(directives_file, env):
    """
    Runs the btree_chromo simulation with a specific directives file.
    """
    try:
        subprocess.run(
            [btree_chromo_executable, directives_file],
            check=True,
            env=env  # Pass the updated environment here
        )
        print(f"Successfully ran btree_chromo with directives: {directives_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running btree_chromo with directives {directives_file}: {e}")
        raise


def main():
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = "/usr/local/lib64:/usr/local/lib:/usr/local/nvidia/lib:/usr/local/nvidia/lib64:/.singularity.d/libs:/Software/LAMMPS/OMP_GPU_Kokkos/lib"

    for step in range(0, 91, 1):
       directives_file = create_directives(step)
       run_btree_chromo(directives_file, env)
       os.remove(directives_file)


if __name__ == "__main__":
    main()
