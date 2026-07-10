# STC-QCB Summer School 2026<br><sub>Modeling the Minimal Bacterial Cell JCVI-Syn3A</sub>

<table>
<tr>
<td width="75%" valign="top">

**Welcome to the STC-QCB Summer School 2026!**

These tutorials are part of the 2026 annual summer school ([Full program information here](https://emails.illinois.edu/newsletter/44/374882244.html)) organized by the NSF Science and Technology Center for Quantitative Cell Biology (STC-QCB) at UIUC.

</td>
<td width="25%" align="center" valign="top">

<img src="./figs/STC_logo.png" alt="STC-QCB logo" width="200">

</td>
</tr>
</table>

### Getting started

The workshop is split into modules, where each module focuses on a different method in wholce-cell modelling. Work through them in the order below, opening a module to begin its session:

|   Day   | Module                                   | Topic                                                                                               |
| :------ | :--------------------------------------- | :-------------------------------------------------------------------------------------------------- |
| Mon–Tue | [Martini](Martini/README.md)             | The Martini ecosystem of tools and whole-cell Martini simulations, ending in a toy JCVI-Syn3A cell. |
| Wed AM  | [DNA](DNA/README.md)                     | Simulating a bacterial chromosome.                                                                  |
| Wed PM  | [Lattice Microbes](LM/README.md)         | Stochastic whole-cell simulations (CME and RDME).                                                   |
| Thu AM  | [CME-ODE](CMEODE_WCM/README.md)          | CME-ODE model of the minimal cell.                                                                  |
| Thu PM  | [4DWCM](4DWCM/README.md)                 | 4D whole-cell modeling with Lattice Microbes.                                                       |

Before the first session, set up your Delta access and confirm you can log in:

- **[QCB Gateway setup](QCB_Gateway_setup.md)** — set up and log in to the QCB Delta Gateway. Covers the SSH tunnel, login, resource allocation, and cloning this repository. Start here.
- **[vmd guide](vmd_guide.md)** — log in to the Open OnDemand Desktop for VMD visualization.

### Useful links

- [STC-QCB website](https://qcb.illinois.edu/): the center behind the school, and its other events.
- [Delta documentation](https://docs.ncsa.illinois.edu/systems/delta/en/latest/): the NCSA cluster the tutorials run on.
- [Jupyter basics](https://www.dataquest.io/blog/jupyter-notebook-tutorial/): the notebook interface, for the modules that run in Jupyter.
- [VMD](https://www.ks.uiuc.edu/Research/vmd/): the visualization program used across the modules.

>[!NOTE]
> These tutorials were written by teaching assistants: Abner Apsley (LM), Enguang Fu (4DWCM and CME-ODE), Andrew Maytin (DNA), and Marieke Westendorp and Jan Stevens (Martini). Technical infrastructure and support by Alfia Parvez.

### Schedule

<table>
  <tr>
    <td><sub><b>Time</b></sub></td>
    <td><sub><b>Monday July 13th, 2026</b></sub></td>
    <td><sub><b>Tuesday  July 14th, 2026</b></sub></td>
    <td><sub><b>Wednesday  July 15th, 2026</b></sub></td>
    <td><sub><b>Thursday  July 16th, 2026</b></sub></td>
    <td><sub><b>Friday  July 17th, 2026</b></sub></td>
  </tr>
  <tr>
    <td><sub><b>08:00</b></sub></td>
    <td colspan="5"><sub>Breakfast in Hendrick House dorms</sub></td>
  </tr>
  <tr>
    <td><sub><b>08:30</b></sub></td>
    <td rowspan="2"><sub>Registration<br>(3269 Beckman)</sub></td>
    <td rowspan="6"><sub>The Martini ecosystem of tools</sub></td>
    <td rowspan="6"><sub>Simulating<br>a bacterial chromosome</sub></td>
    <td rowspan="6"><sub>CME-ODE model<br>of the Minimal Cell</sub></td>
    <td rowspan="6"><sub>Presentation Preparation<br>(3269 Beckman)</sub></td>
  </tr>
  <tr>
    <td><sub><b>09:00</b></sub></td>
  </tr>
  <tr>
    <td><sub><b>09:30</b></sub></td>
    <td><sub>Welcome: Sharlene Denos</sub></td>
  </tr>
  <tr>
    <td><sub><b>10:00</b></sub></td>
    <td><sub>Zaida Luthey-Schulten</sub></td>
  </tr>
  <tr>
    <td><sub><b>10:30</b></sub></td>
    <td><sub>Rohit Bhargava</sub></td>
  </tr>
  <tr>
    <td><sub><b>11:00</b></sub></td>
    <td><sub>Group Photo (Beckman Quad)</sub></td>
  </tr>
  <tr>
    <td><sub><b>11:30</b></sub></td>
    <td rowspan="3"><sub>Lunch &amp; TA Posters<br>(Beckman Atrium)</sub></td>
    <td><sub>Lunch (in 612 IGB)</sub></td>
    <td><sub>Lunch (in 612 IGB)</sub></td>
    <td><sub>Lunch (outside 3269 Beckman)</sub></td>
    <td><sub>Lunch (outside 3269 Beckman)</sub></td>
  </tr>
  <tr>
    <td><sub><b>12:00</b></sub></td>
    <td><sub>Shulei Wang (612 IGB)</sub></td>
    <td rowspan="2"><sub>Science Communication Workshop.   Sharlene Denos (612 IGB)</sub></td>
    <td><sub>Hyun Youk (3269 Beckman)</sub></td>
    <td rowspan="2"><sub>Setup student presentations<br>(3269 Beckman)</sub></td>
  </tr>
  <tr>
    <td><sub><b>12:30</b></sub></td>
    <td><sub>Chris Maffeo (612 IGB)</sub></td>
    <td><sub>Jonas Zaehringer (3269 Beckman)</sub></td>
  </tr>
  <tr>
    <td><sub><b>13:00</b></sub></td>
    <td rowspan="10"><sub>Introduction<br>to Martini-GROMACS simulations</sub></td>
    <td rowspan="10"><sub>Whole-cell simulations <br>using Martini-GROMACS</sub></td>
    <td rowspan="10"><sub>Stochastic Simulations<br>in Lattice Microbe</sub></td>
    <td rowspan="10"><sub>4D whole-cell modeling<br>using Lattice Microbes</sub></td>
    <td rowspan="2"><sub>Student Presentations for Theme 1: Gene expression at nanometer resolution using MINFLUX <br>(3269 Beckman)</sub></td>
  </tr>
  <tr>
    <td><sub><b>13:30</b></sub></td>
  </tr>
  <tr>
    <td><sub><b>14:00</b></sub></td>
    <td><sub>Coffee Break</sub></td>
  </tr>
  <tr>
    <td><sub><b>14:30</b></sub></td>
    <td rowspan="2"><sub>Student Presentations for Theme 2: Experimental and Computational Studies of the Minimal Cell<br>(3269 Beckman)</sub></td>
  </tr>
  <tr>
    <td><sub><b>15:00</b></sub></td>
  </tr>
  <tr>
    <td><sub><b>15:30</b></sub></td>
    <td><sub>Coffee Break</sub></td>
  </tr>
  <tr>
    <td><sub><b>16:00</b></sub></td>
    <td rowspan="2"><sub>Student Presentations for Theme 3: Chemical imaging, Spectroscopy, and Computation Module<br>(3269 Beckman)</sub></td>
  </tr>
  <tr>
    <td><sub><b>16:30</b></sub></td>
  </tr>
  <tr>
    <td><sub><b>17:00</b></sub></td>
    <td rowspan="2"><sub>Evaluation Survey<br>(3269 Beckman)</sub></td>
  </tr>
  <tr>
    <td><sub><b>17:30</b></sub></td>
  </tr>
  <tr>
    <td><sub><b>18:00</b></sub></td>
    <td rowspan="6"><sub>Dinner &amp; Tours                                (Minecraft in 3269 Beckman and Microscopy Core tours in IGB basement)</sub></td>
    <td rowspan="2"><sub>Dinner on your own</sub></td>
    <td rowspan="2"><sub>Dinner on your own</sub></td>
    <td rowspan="2"><sub>Dinner @ Riggs</sub></td>
    <td rowspan="6"><sub>Students Depart</sub></td>
  </tr>
  <tr>
    <td><sub><b>18:30</b></sub></td>
  </tr>
  <tr>
    <td><sub><b>19:00</b></sub></td>
    <td rowspan="4"><sub>continuation <br>(if needed)</sub></td>
    <td rowspan="4"><sub>continuation <br>(if needed)</sub></td>
    <td rowspan="4"><sub>continuation <br>(if needed)</sub></td>
  </tr>
  <tr>
    <td><sub><b>19:30</b></sub></td>
  </tr>
  <tr>
    <td><sub><b>20:00</b></sub></td>
  </tr>
  <tr>
    <td><sub><b>20:30</b></sub></td>
  </tr>
</table>
