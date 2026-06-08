# <img src="./figs/STC_logo.png" alt="STC-QCB logo" height="65"> STC-QCB Summer School 2026: Modeling the Minimal Bacterial Cell JCVI-Syn3A

**Welcome to the STC-QCB Summer School 2026!**

These tutorials are part of the 2026 annual summer school ([Full program information here](https://emails.illinois.edu/newsletter/44/374882244.html)) organized by the NSF Science and Technology Center for Quantitative Cell Biology (STC-QCB) at UIUC.

### Schedule

In this computational study of the genetically minimized bacterium, JCVI-Syn3A, you will learn to:

1. Simulate the chromosome dynamics of Syn3A using LAMMPS  
2. Simulate coupled gene expression and metabolism of Syn3A with a spatially homogeneous whole-cell model (WCM) implemented by [**Lattice Microbe**](https://github.com/Luthey-Schulten-Lab/Lattice_Microbes)  
3. Analyze spatially heterogeneous trajectories from the 4DWCM (3D in space plus time) of Syn3A 

### Getting Started: QCB Delta Gateway

The CME and RDME Jupyter tutorials run through the **QCB Delta Gateway** — a JupyterHub front-end on Delta. Set it up once, then open notebooks from any module.

**Step 1: Open an SSH tunnel**

Open a terminal on your laptop and run the following command. Replace `USERNAME` with your NCSA username:

```bash
ssh -L 8000:dt-svc-bbkw01.hsn.cm.delta.internal.ncsa.edu:8000 USERNAME@login.delta.ncsa.illinois.edu
```

You will be prompted for your **NCSA password** and **two-factor authentication (2FA)**. Once you're in, **leave the terminal open** — closing it tears down the tunnel.

**Step 2: Open the gateway in your browser**

Once the SSH tunnel is up, open this URL in any browser on your laptop:

```
https://dt-svc-bbkw01.delta.ncsa.illinois.edu:8000/hub/org/
```

Click on the **QCB Gateway** tab. You should see the **JupyterHub login page** for the QCB Delta Gateway.

<img src="./figs/QCB_Gateway_Login.png" alt="QCB Gateway homepage" width="700">

Click **CI Logon** and sign in with your NCSA Delta credentials.

<img src="./figs/QCB_Gateway_homepage.png" alt="QCB Gateway CI Logon page" width="700">

> [!NOTE]
> If your Gateway account is not approved, please ask the admin, Alfia Parvez, at alfiap@illinois.edu to approve it first.

**Step 3: Allocate resources**

After logging in, the Gateway will ask you to choose compute resources before starting your Jupyter session. Use the settings below:

| Setting | Value |
| --- | --- |
| **Allocation** | **A100 GPU - up to 8 (bgvl-delta-gpu)** — choose the **Batch** (non-interactive) option, not Interactive |
| **Number of CPUs** | **8** |
| **GPU Environment** | **4DCell (LAMMPS/LM)** |
| **Number of GPUs** | **1** |
| **Memory** | **64 GB** |
| **Time limit** | **3 hours** |

<img src="./figs/Resource_Allocation.png" alt="QCB Gateway resource allocation form" width="700">

Click **Start**. 

Jobs may take a short time to queue before your session opens. Once the server is ready, clone the repository from the notebook that starts automatically (Step 4).

**Step 4: Clone the repository**

When your Jupyter session starts, an **Untitled.ipynb** notebook will already be open in JupyterLab.

<img src="./figs/jupyter_env.png" alt="JupyterLab environment with Untitled.ipynb" width="700">

In a **code cell**, run:

```python
!git clone https://github.com/Luthey-Schulten-Lab/SummerSchool_2026.git
%cd SummerSchool_2026
```

Open notebooks from your cloned copy in the Jupyter file browser:

- **CME:** `SummerSchool_2026/CME/`
- **RDME:** `SummerSchool_2026/RDME/`

RDME pre-computed trajectory data (CSVs and VMD files, not in git) is read from the shared bgvl folder:

```
/projects/bgvl/SummerSchool_2026/RDME/data/
/projects/bgvl/SummerSchool_2026/RDME/trajectory/
```

- **Martini:** open notebooks under `SummerSchool_2026/Martini/tutorial_1/` (etc.)

### DNA tutorial (QCB Gateway)

Open **[`DNA/submit_simulation.ipynb`](DNA/submit_simulation.ipynb)** on the Gateway to copy the workshop template into your bgvl folder, then submit the GPU simulation as a **Slurm** job from a Delta SSH login (`sbatch DNA/files/launch_simulation.sh <your-workspace>`) — it runs inside the `DNA_summer2025.sif` Apptainer image. The notebook prints the exact command with your paths (your Gateway username differs from your Delta login name). Monitor the log back on the Gateway. See [`DNA/README.md`](DNA/README.md).
### Useful Links
0. [STC-QCB Official Website](https://qcb.illinois.edu/): Check out the QCB website for more events!  
1. [Delta User Help Document](https://docs.ncsa.illinois.edu/systems/delta/en/latest/): We will use **Delta**, a high-performance computing cluster housed at the National Center for Supercomputing Applications (NCSA), UIUC, to conduct the computational tasks.  
2. [Common Terminal Commands](https://gist.github.com/bradtraversy/cc180de0edee05075a6139e42d5f28ce): You need to know basic terminal commands to navigate directories and launch jobs.  
3. [Jupyter Notebook Basics](https://www.dataquest.io/blog/jupyter-notebook-tutorial/): Jupyter Notebook provides a user-friendly interface to run simulations interactively.  
4. [Try Jupyter Notebook Online](https://jupyter.org/try): Try it out to get familiar with the interface.  
5. [Visual Molecular Dynamics (VMD)](https://www.ks.uiuc.edu/Research/vmd/): **VMD** is a versatile tool that will be used to visualize chromosome dynamics simulated by LAMMPS and 4DWCM trajectories simulated by Lattice Microbe.  

>[!NOTE]  
> These tutorials were written by teaching assistants: Andrew Maytin (DNA), Enguang Fu (CME), and Tianyu Wu (RDME).
