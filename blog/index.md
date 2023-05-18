@def title = "Notebook Blogs"
@def author = "Stefan Bringuier"


The notebooks are rendered to html via [PlutoSliderServer.jl](https://github.com/JuliaPluto/PlutoSliderServer.jl) and use [Julia version 1.8 ](https://julialang.org/downloads/#current_stable_release). Feel free to post an [issue](https://github.com/stefanbringuier/PlutoNotebookBlogging/issues) if there are errors or problems running notebooks.

{{ get_pluto_tags }}

# Blogs

@@blog-list
 @@blog-item  [Inverse determination of diffusion coefficient](assets/notebooks/BayesianInferenceDiffusion.html) - Using Bayesian inference to inversely determine diffusion coefficient. @@
 @@blog-item  [Simple NN in Julia](assets/notebooks/FluxTutorial.html) - Reworking the [quickstart of Flux.jl](https://fluxml.ai/Flux.jl/stable/models/quickstart/). @@
 @@blog-item [Plotting properties on a periodic table](assets/notebooks/PeriodicTableSchematic.html) - It can be useful at times to color code periodic table to display elemental properties. @@
 @@blog-item [Creating slides using julia](assets/notebooks/slidesviajulia.html) - Demonstration of creating slides using [PPTX.jl](https://github.com/ASML-Labs/PPTX.jl) package. @@
 @@blog-item [Dual numbers and derivatives](assets/notebooks/DualNumbers.html) - Discussing dual numbers and the use for taking derivatives of functions. @@
 @@blog-item [LAMMPS in notebook](assets/notebooks/LAMMPS_Julia.html) - Example running LAMMPS in a notebook and accessing compute result. @@
 @@blog-item [Simulating CVD growth](assets/notebooks/KMC_CVD.html) - Implementation of KMC model to simulate chemical vapor deposition. @@
 @@blog-item [Integrating a function with MC](assets/notebooks/TrivialMC.html) - Basic demo showing how Monte Carlo sampling can be used to integrate a function. @@
 @@blog-item [DFT calculation in a notebook](assets/prerendered_notebooks/AtomicCalculationWorkflow.html) - Demo the use of PyMatgen, ASE, and DFTK.jl in a Pluto.jl notebook. @@
 @@blog-item [Featurization & ML Models](assets/prerendered_notebooks/RemakeBestPracticesPost.html) - A guide through featurization and ML models for materials data. @@
 @@blog-item [Bayesian Inference](assets/notebooks/BMCP_Ch2_1.html) - An example from *Bayesian Modeling and Computation in Python* but in Julia. @@
 @@blog-item [ML in Materials Science](assets/notebooks/ML_MatSci.html) - Plotting the publication activity of Machine Learning use in Materials Science. @@
 @@blog-item [Graphical Intuition of Classification](assets/notebooks/Classifier_GraphicalIntuition.html) - Simple post on graphical intuition behind classification. @@
 @@blog-item [Guassian Processes: Part 1](assets/notebooks/gaussianprocess_part1.html) - Initial workflow showing how Gaussian processes work. @@
 @@blog-item [Density Matrix Renormalization Group](assets/notebooks/dmrg.html) - A walk through basic 1D infinite anti-ferromagnetic spin chain. @@
 @@blog-item [Chemical Kinetics ODE](assets/notebooks/chemkinetics_ode.html) - Solving mass action ODE using very nice Julia modeling packages. @@
 @@blog-item [Welcome](blogpages/welcome) - My hello world. @@
@@
