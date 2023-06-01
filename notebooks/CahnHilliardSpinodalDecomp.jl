### A Pluto.jl notebook ###
# v0.19.25

#> [frontmatter]
#> author = "Stefan Bringuier"
#> Email = "stefanbringuier@gmail.com"
#> title = "Cahn-Hillard Equation & Spinodal Decomposition Simulation"
#> date = "2023-05-31"
#> tags = ["Materials Science", "Simulation", "Simulations", "Numerics", "Differential Equations"]
#> description = "Using the phase field approach solve the Cahn-Hillard equation for standard example"
#> license = "CC-BY-4.0"

using Markdown
using InteractiveUtils

# ╔═╡ 8f199e33-f971-4f96-8883-779f7cc6a946
md"""
# Cahn-Hilliard Equation & Spinodal Decomposition Simulation
**Author: Stefan Bringuier**

!!! note 
    This is an old Jupyter notebook I had which I converted to a pluto notebook using the utility available [here](https://observablehq.com/@amgalanb/pluto-jl-jupyter-conversion). It works pretty well. The notebook is based on my notes from serveral phase field textbooks. One particularly good one is ref. [^1].
"""

# ╔═╡ 1079052f-234f-4698-9bb8-af49b24ef80a
md"""
## Cahn-Hilliard - Part 1
"""

# ╔═╡ 518e6304-8f40-4ef7-925f-e16558080b60
md"""
The Cahn-Hillard model was developed from a phase-separation model focused on spinodal decomposition of a binary alloy. Spinodal decomposition is the process whereby chemical species will move against a compositional gradient. The process is not driven by nucleation and growth process, just diffusion.

The characteristic free-energy functional for spinodial decomposition is:

$$F = \int_{V} \left[ f(c) + \frac{1}{2} \kappa \left(\nabla c\right)^2 \right] dv $$

where $\kappa$ is the gradient energy penalty coefficient and $f(c)$ is the chemical/bulk free energy density represented by:

$$f(c) = Ac^2(1-c)^2$$

This is a simple-double well potential as plotted below ($f(c)$ is intentionally bounded between $0<f(c)<1$). Here $A$ is a positive constant that controls the magnitude of the well-barrier. Other free energy potentials are:

$$f(c) = \frac{1}{4}A\left(1-c^2\right)^2$$

$$f(c) = 4A\left(-\frac{1}{2}c^2 + \frac{1}{4}c^4\right)$$
"""

# ╔═╡ cc84fab4-f9fe-498b-8136-9e997b8a6d09
md"""
### Potential-Well Examples
"""

# ╔═╡ 90c5e28c-3c05-43a3-bb4c-2f0a47f5f339
begin
	A=1.00e0
	f(c)=A*(c^2)*((1-c)^2)
	f1(c)=0.25e0*A*(1-c^2)*(1-c^2)
	f2(c)=4.00e0*A*(-0.5e0*c^2 + 0.25e0*c^4)
	conc = range(-0.20e0,stop=1.20e0,length=120)
	
	dw1 = plot(conc,map(f,conc),title="double-well")
	dw2 = plot(conc,map(f1,conc),title=text("zero-centered well",6))
	dw3 = plot(conc,map(f2,conc),title=text("shifted zero-centered well",6))
	plot(dw1,dw2,dw3,linewidth=2,size = (675,225),layout=(1,3),
	grid=false,legend=false,titlefont = font(10))
end

# ╔═╡ 4432a409-ce78-41bb-863f-dde56f1885f1
md"""
The Cahn-Hillard equation is given by:

$$\frac{\partial c}{\partial t} = \nabla^2 M \left( \frac{\delta F_{ij}}{\delta c} \right)$$

Here $M$ is the mobility and $\frac{\delta F_{ij}}{\delta c}$ is the variational derivative of the free energy functional.

Let us first express the left-hand side of the equation in terms of Explicit Euler time integration:

$$\frac{c_{ij}^{n+1}-c_{ij}^{n}}{\Delta t} = \nabla^2 M \left( \frac{\delta F_{ij}}{\delta c} \right)^n$$

$n$ is the current time interval/step.

The variational derivative of $\frac{\delta F_{ij}}{\delta c}$ is given by:

$$\frac{\delta F_{ij}}{\delta c}^{n} = \mu\left(c_{ij}^n\right) - \kappa \nabla^2 c_{ij}^{n}$$

where $\mu(c_{ij}^{n})$ is the variational derivative of the chemical/bulk free energy and is:

$$\mu(c_{ij}^{n}) = A\left(2c_{ij}^{n}\left(1-c_{ij}^{n}\right)^{2} + 2\left(c_{ij}^{n}\right)^{2}\left(c_{ij}^{n}-1\right)\right)$$

The laplacian operator in the Cahn-Hillard equation can be approximated by the finite-difference 5-point stencil:

![Finite-Difference Stencil](https://github.com/stefanbringuier/PlutoNotebookBlogging/tree/gh-pages/assets/images/finite_diff_grid_fig4.1_SBiner.png)
Finite-difference stencil adapted from [^1]

The discrete energy functional is written as:

$$F_{bulk} = \sum_i^{Nx} \sum_j^{Ny} f(c_{ij}) = Ac_{ij}^2(1-c_{ij})^2$$
$$F_{grad} = \frac{\kappa}{2} \sum_i^{Nx-1} \sum_j^{Ny-1} \left( \left(c_{i+1,j} - c_{ij}\right)^2 + \left(c_{i,j+1}-c_{ij}\right)^2\right)$$
$$F_D = F_{bulk} + F_{grad}$$

Here $F_{grad}$ is calculated using the spatial forward  discrete difference.
"""

# ╔═╡ 45c4ea63-050c-44ac-aeb1-0507340f2675
md"""
A useful tool is to establish non-dimensional characteristic quantities.

For length we can use the ratio of the gradient energy penalty(e.g. J/m$^3$) to the free energy barrier (e.g. J) and then account for the 2D simulation, thus:

$$L^{\prime} = \left(\frac{\kappa}{A}\right)^{\frac{1}{2}}  \text{[Length]}$$

The same can be done for other quantities:

$$F^{\prime} = AL^3 \text{[Energy}-\text{Length}^3]$$

$$t^{\prime} = \frac{L^{\prime 2}\left(c_p^e - c_m^e \right)}{MF^{\prime}} \text{[Time-Mass / Energy-Length]}$$

here $c_p^e$ and $c_m^e$ are the equilibrium concentration of the two phases.

For these simulations the nondimensional time increment was 0.01 with a Euler time integration scheme.

Total number of timesteps is 20,000.

"""

# ╔═╡ a96a12c7-027e-43ba-a4c6-8d013f651899
md"""
# Cahn-Hilliard Finite-Difference Code

### Data Types
Lets first create Julia composite data types that will be useful. We will need a data type that specifies our spatial grid, ```Spatial2D```, a datatype for the concentration field, ```Conc2D```, a simulation time data type,```Time```, and a data type ```Material``` specifying the material parameters.


!!! note
    For Julia data types and function names I will use "CamelCaps" even though this is not a suggested julia style guide.
"""

# ╔═╡ f04ec0c5-fa5a-4467-86d2-9745ffba834a
# Number of grid points, spacing size, and grid area
struct Spatial2D
    nx::Int64;
    ny::Int64;
    dx::Float64;
    dy::Float64;
    size::Int64;
    cellarea::Float64;
    area::Float64;
    function Spatial2D(nx,ny,dx,dy)
        grid = nx*ny;
        cell = dx*dy;
        area =(dx*nx)*(dy*ny);
        new(nx,ny,dx,dy,grid,cell,area);
    end
end

# ╔═╡ 53045f48-b321-4c8c-8fdb-9d1369faa386
#Material parameters
struct Material
    avg_c::Float64;    #average concentration
    mobile::Float64;   #Cahn-Hillard mobility coefficient
    grad_pen::Float64; #Gradient penality coefficent
    pote::Float64;  #Potential well energy barrider
    Material(i,m,g,p) = new(i,m,g,p)
end

# ╔═╡ 67b230cd-77e2-4ecd-83cb-47768ff457c6
struct Field2D
    spatial::Spatial2D #Store spatial info
    material::Material #Store material parameters
    field::Array;
    grad::Array;
    laplacian::Array;
    function Field2D(s2D::Spatial2D,mat::Material)
        nx,ny = s2D.nx,s2D.ny;
        c=zeros(s2D.nx,s2D.ny);
        g=zeros(s2D.nx,s2D.ny);
        l=zeros(s2D.nx,s2D.ny);
        new(s2D,mat,c,g,l);
    end
end

# ╔═╡ 6a8475b7-7ba5-4310-a1f6-d3529cc01651
#Time and screen output values
mutable struct Time
    ntstep::Int64;    #number of timesteps
    prtfreq::Int64;   #print to screen every ntstep
    timestep::Float64;#numerical time step
    tinit::Float64;   #initial simulation time
    Time(nt,pf,ts,ti) = new(nt,pf,ts,ti)
end

# ╔═╡ 39760bb3-05f8-4a9c-89f2-2fcd46cab544
md"""
## Simulation Setup & Energy Functions

Now that we have the datatypes we need to define the functions which will process and operate on our data types. The first thing we need is a function to generate a initial microstructure. We will call it ```GenMircoStruct2D``` since it will generate a two-dimensional microstructure, the function will use the average concentration field and apply some random noise then return a new concentration field which is of type ```Conc2D```.
"""

# ╔═╡ eeb623fd-558d-4c1e-9fe7-677fc704ee91
"""function GenMicroStruc2D(spatial info, material parameters)
        This function is used to generate a 2D concentration profile.

        Input - spatial2D::Spatial2D datatype see defined datatypes.
                conc2D::Conc2D datatype see defined datatypes.
                material::Material datatype see defined datatypes.

        Output - conc2D::Conc2D datatype with modified conc2D.conc values. See defined datatypes
        
        Used variables:
             avg_c - average concentration
             noise*rand() - random noise added to concentration
"""
function GenMicroStruc2D(spatial2D::Spatial2D,material::Material,noise::Float64)
    conc = Field2D(spatial2D,material)
    for i=1:conc.spatial.nx
        for j=1:conc.spatial.ny
            conc.field[i,j]=conc.material.avg_c + noise*(0.50e0-rand())
        end
    end
    return conc
end; #function GenMicroStruc2D

# ╔═╡ c7cc130f-54c8-4c49-8072-818cc4b7d9b3
@mdx """
#### Simulation parameters & microstructure
Now that we have a function to generate an initial microstructure, let us define the simulation details with the data types we created and the ```GenMicroStruct2D``` function. We will plot a contour view of the concentration field.
"""

# ╔═╡ b65e7452-9ae7-4762-8c3a-e1f0b3f2744b
begin
	#Initiate simulation parameters
	nx,ny = 64,64;
	dx,dy = 1.00e0,1.00e0;
	spatial = Spatial2D(nx,ny,dx,dy);
	steps,prnt,dt,ti = 1000,50,1.00e-2,0.00e0;
	time = Time(steps,prnt,dt,ti);
	c̄,M,κ,Aₒ = 0.40e0,1.00e0,0.50e0,1.00e0;
	material = Material(c̄,M,κ,A);
	c = GenMicroStruc2D(spatial,material,0.02e0);
	heatmap(c.field,size=(475,400))
end

# ╔═╡ 0f84a539-dce0-4dc2-b6ff-85bf759b3c99
md"""
#### Boundary Conditions

Here we will assume periodic boundary conditions, which are simple and straightfoward to implement. To do this we will define a function which takens in the grid location $i,j$ and determine if the points lie on the boundary of the simulation cell. If they do then the index will be shifted. The function is going to be named ```PeriodBoundCond```
"""

# ╔═╡ a36fa201-58da-4e78-8ed2-1127fa730d7a
""" function PeriodBoundCond(index i,max i)
         This function determines if neighbor indices need to be shifted based on periodic boundary conditions

         Input - i :: Integers  index of a point on 2D grid
                 mi :: Integers max index of the 2D grid

         Output - ni:: Integers unchanged or shifted index
""" function PeriodBoundCond(i::T,mi::T) where T<:Integer
        iplus,iminus = i+1,i-1;
        if iplus > mi 
            iplus = 1;
        elseif iminus == 0
            iminus = mi;
        end
        return (iplus,iminus)
end; # function PeriodBoundCond

# ╔═╡ f8ff2918-993a-4527-8ec7-00f96b4bfe13
md"""
#### Free Energy functions
"""

# ╔═╡ cf37924a-cfe4-4336-a9b3-88c2c6c268f1
md"""
The next step is to define the total free energy of the system based on our free-energy functional. In this case we use the double-well description to define a function called ```TotalFreeEnergy```. We also need the variational derivative of the functional $\frac{\delta F}{\delta c}$, which is given by the derivative of the integrand. Here I'll define the function ```DervBulkFreeEnergy``` as the derivative at a specific site, $ij$, on the 2D spatial grid
"""

# ╔═╡ 1f7e8f70-eef4-40da-82a0-7050dec53f76
"""function TotalFreeEnergy(concentration field, coefficient for gradient penalty)
        This function calculates the total simulation cell free energy. Assumes standard
        Cahn-Hillard free energy functional given by:
        f(c) + \\frac{1}{2}\\kappa\\left(\\nabla c \\right)^2
        uses forward-differences.

        Input - field::Field2D datatype see defined datatypes.

        Output - free_energy::Real
        
        Question: why don't we include the summation over the discrete volume?

"""
function TotalFreeEnergy(field::Field2D)
    free_energy = 0.00e0;
    for i=1:field.spatial.nx
        iplus, = PeriodBoundCond(i,field.spatial.nx);
        for j=1:field.spatial.ny
            jplus, = PeriodBoundCond(j,field.spatial.nx);
            #Indexing variables
            fij = field.field[i,j];
            fiplus = field.field[iplus,j];
            fjplus = field.field[i,jplus];
            #Energy summation over domain
            enrgy_bulk = fij^2 *(1.00e0-fij)^2;
            enrgy_grad = 0.50e0 * field.material.grad_pen * (fiplus-fij)^2 + (fjplus-fij)^2;
            free_energy = free_energy + enrgy_bulk + enrgy_grad;
        end
    end
    return free_energy
end; #function TotalFreeEnergy

# ╔═╡ e63f6699-7a70-4bd6-ac4c-46bc406cdc23
""" δ(derivative of field at (i,j))
        This function calculates the derivative of the bulk free energy for a double-well
        potential given by:

        \\mu(c_{ij}^{n}) = A\\left(2c_{ij}^{n}\\left(1-c_{ij}^{n}\\right)^{2} + 
        2\\left(c_{ij}^{n}\\right)^{2}\\left(c_{ij}^{n}-1\\right)\\right)

        Input - f ::Float64 value for the field at i,j.

        Output - Real value for derivative of bulk free energy at i,j

""" δij(f::Float64) = 2.00e0 * f * (1.00e0-f)^2 - 2.00e0 * f^2 * (1.00e0-f);

# ╔═╡ 69f55915-f802-4f58-839f-60897a864c97
md"""
## Define Laplacian and concentration derivative functions

We now need to do a few things to evolve the spatial and concentration derivatives. The first thing we need to define a function which takes the Laplacian of a field, in this case the field is the concentration. Let us call this function $\nabla^{2}$ The second thing needed is a routine to calculate the variation derivative of the free-energy with respect to concentration, this funciton will be $\delta F!$. Note that the function modifies the conc arguement and returns a free-energy derivative at every point.
"""

# ╔═╡ f3948247-6b50-4329-aef2-e524c0e781ee
""" function δF!(field to take central difference laplacian of)
                 Input - field :: Contains 2D Array of field to take laplacian of

                 Output - dFdc :: Array of variational derivative of Free-energy on
                                  spatial domain grid.
""" function δF!(field::Field2D)
        nx,ny = field.spatial.nx,field.spatial.ny;
        dx,dy = field.spatial.dx,field.spatial.dy;
        #Array for variational derivative dF/dc (derivative of integrand)
        dFdf = zeros(Float64,nx,ny);
    
        #Solve for δF/δf
        for i=1:nx
        iplus,iminus = PeriodBoundCond(i,nx)
        ii = (iminus,i,iplus);
            for j=1:spatial.ny
                jplus,jminus = PeriodBoundCond(j,ny)
                jj = (jminus,j,jplus);
        
                #Calculate the free energy $\frac{\delta F}{\delta c_{ij}}$
                field.laplacian[i,j] = ∇²(ii,jj,dx,dy,field.field)
            
                #Calculate the bulk free energy derivative in free energy expression
                barrier = field.material.pote;
                dfengry_dfij = barrier*δij(field.field[i,j]);
            
                #Store the variational free energy derivative, 
                # e.g. dF/dcij in right-hand side of Cahn-Hillard eq.
                dFdf[i,j] = dfengry_dfij - field.material.grad_pen * field.laplacian[i,j]          
            end
        end
        return dFdf
end;# function δF!

# ╔═╡ f11e9ecb-6bbd-468d-a0b5-cdb315c6a8ed

""" function Laplacian2D(indices for i, indices for j,
                         x discrete distance, y discrete distance,
                         field to take Laplacian of)
                 Take central difference laplacian of a field

                 Input - i,j :: Tuple of integer indice about point on 2D grid
                         dx,dy :: discrete spatial spacing between grid points
                         field :: 2D Array of field to take laplacian of

                 Output - ∇²f::Laplacian at point i,j of field
""" function ∇²(i::Tuple{T,T,T},j::Tuple{T,T,T},
                dx::S,dy::S,
                field::Array{S,2}) where {T<:Integer,S<:Real}
               
                #i[1] = i-1, i[2] = i, i[3]= i+1
                ∇²f = (field[i[1],j[2]] + field[i[3],j[2]] + 
                field[i[2],j[1]] + field[i[2],j[3]] - 
                4.00e0 * field[i[2],j[2]]) / (dx*dy);
                
                return ∇²f
end; #function ∇²

# ╔═╡ 827cecdf-1345-4149-aaa0-d752c0dff201
md"""
## Spatial Solution to R.H.S of Cahn-Hilliard Equation

Now we want to bring things together by solving for the  $\nabla M\left(\frac{\delta F}{\delta c}\right)$. In this notebook, the material mobility interface parameter is constant so the solution is simplified by taking the Laplacian of the variational term. We will define a function similar to $\delta F!$ that does this for us.
"""

# ╔═╡ c10a9d5b-c92a-452f-8302-5032a9165a06
function ∇²δF(spatial::Spatial2D,
                 dFdf::Array{T,2}) where T<:Real
    
       nx,ny = spatial.nx,spatial.ny;
       dx,dy = spatial.dx,spatial.dy;
       ∇²dFdf = zeros(Float64,nx,ny);
    
       #Solve for laplacian R.H.S. of Cahn-Hillard eq.:
       #$\nabla^{2}\left(\frac{\delta F}{\delta c_{ij}}\right)$
        for i=1:nx
        iplus,iminus = PeriodBoundCond(i,nx)
        ii = (iminus,i,iplus);
        for j=1:ny
           jplus,jminus = PeriodBoundCond(j,ny)
           jj = (jminus,j,jplus);
           ∇²dFdf[i,j] = ∇²(ii,jj,dx,dy,dFdf);
        end
    end
    return ∇²dFdf
end; #function ∇²δF            

# ╔═╡ 161197f4-1076-4fb0-81de-b34c927ca4c2
md"""
## Time evolution function and Run Simulation

With the neccessary spatial and field derivative functions defined that make up the right hand side of the Cahn-Hillard Equation, we can perform numerical integration to time evolve and find the solution of the concentration field at time, $t$, given our initial conditions (e.g., random distribution of concentration) and boundary condtions (e.g., periodic). The function to do this will be called ```TimeEvolveCH!``` and will take in the ```Time``` and ```Field2D``` datatypes.
"""

# ╔═╡ a4bffbdd-a2a1-42ce-92c9-67332738706d
function SimInfo(istep::Integer,energy::Float64,
                 time::Time,field::Field2D;skip=true)
    if skip == true
        return 0.00e0
    else
        #We want to confirm total energy is conserved.
        ΔE = (ototalenergy-TotalFreeEnergy(field)) 
        if  abs(ΔE) > tole
            @warn "Change in total energy is $(ΔE)"
        end
        energy = TotalFreeEnergy(field);
            
        if (mod(istep,time.prtfreq) == 0) || (istep == 1)
            @printf("%5d %14.6e\n",istep*Δt,totalenergy)
        end
        
        return energy
    end
end                

# ╔═╡ 88be3fab-f55b-4ea4-a01c-90c02b7b15a1
""" function TimeEvolveCH(simulation time parameters,
                          field to be evolved under C-H equation)

             Input - time::Time simulation time info.
                     field::Field2D field to evolve under C-H eq.

             Output - None                    
""" function TimeEvolveCH!(time::Time,field::Field2D;tole=1.0e-1)
        
        curtime = time.tinit;
        M=field.material.mobile;
        Δt=time.timestep;
        energy = TotalFreeEnergy(field);
    
        for istep=1:time.ntstep
            curtime += time.timestep;
        
            #Physics
            δFδf = δF!(field);
            ∇²δFδf = ∇²δF(field.spatial,δFδf);
            field.field[:,:] += Δt * M * ∇²δFδf[:,:]; #Integrate CH eq.

                    
            #Simulation checks
            energy = SimInfo(istep,energy,time,field);
          
        end
    
        time.tinit = curtime;
end; #function TimeEvolveCH!

# ╔═╡ d57fca19-7a1e-4189-8672-784ed69819ce
begin
	TimeEvolveCH!(time,c)
	heatmap(c.field,size=(475,400))
end

# ╔═╡ b2f5f523-cea7-4b2a-8e91-44e466cfa2bb
md"""
## Run for longer time 
We can see the continued effect of coarsening by running for longer simulations times.
"""

# ╔═╡ ae2293e5-66e5-4711-b351-a84225570a86
begin
	time.ntstep = 1000
	TimeEvolveCH!(time,c)
	heatmap(c.field,size=(475,400))
end

# ╔═╡ 918532ec-8353-4966-96f3-1e6126f1e916
md"""
Lets run for longer once again.
"""

# ╔═╡ 1a910c0f-1bcc-4755-87a1-63c2496f0551
begin
	time.ntstep = 8000
	TimeEvolveCH!(time,c)
	heatmap(c.field,size=(475,400))
end

# ╔═╡ ba7b2db6-8fd5-4050-8983-f5d58147dbd4
@mdx """
# Re-run
If you want to change the simulation setup parameters go back [here](#simulation-parameters--microstructure)
"""

# ╔═╡ 0e36f365-de6b-4019-8b19-ee77c36dd344
md"""
## Potential Issues

There is some difference in the total free energy as a function of time, this only occurs in the initial stages and then decays. The spinodal coarsening seems to be acurately captured.
"""

# ╔═╡ c4a319b9-c164-4e6d-afbe-da91b7a95846
@mdx """
# References
[^1]:$(DOI("10.1007/978-3-319-41196-5"))
"""

# ╔═╡ fcc0b4f7-32bf-454d-8ec8-3a50139c984d
md"""
## Packages
"""

# ╔═╡ 1c1cd3ea-0ea5-4637-83c2-63ec10d7cbdd
# ╠═╡ show_logs = false
md"""
$(let
	using Pkg
	tb = "| Package | Version |\n| :-- | :-- |\n"
	tb *= "| Julia   | $VERSION |\n"
	for (name, version) in Pkg.installed()
	    tb *= "| $name | $version |\n"
	end
	Markdown.parse(tb)
end)
"""

# ╔═╡ 1b8a9dc3-1cbf-4ab0-93c3-d31eb1bf0353
begin
	using Plots
	plotly()
	using Printf
	using PlutoUI
	using ShortCodes
	using MarkdownLiteral: @mdx
end;

# ╔═╡ 9d5135d6-0ba1-435d-884b-a74dbc906d61
TableOfContents()

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
MarkdownLiteral = "736d6165-7244-6769-4267-6b50796e6954"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
ShortCodes = "f62ebe17-55c5-4640-972f-b59c0dd11ccf"

[compat]
MarkdownLiteral = "~0.1.1"
Plots = "~1.38.14"
PlutoUI = "~0.7.51"
ShortCodes = "~0.3.5"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "5994556fca373a40ea7c86cef4075775332d0858"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "be6ab11021cd29f0344d5c4357b163af05a48cba"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.21.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.CommonMark]]
deps = ["Crayons", "JSON", "PrecompileTools", "URIs"]
git-tree-sha1 = "532c4185d3c9037c0237546d817858b23cf9e071"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.12"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "96d823b94ba8d187a6d8f0826e731195a74b90e9"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "738fec4d684a9a6ee9598a8bfee305b26831f28c"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.2"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "8b8a2fd4536ece6e554168c21860b6820a8a83db"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "19fad9cd9ae44847fe842558a744748084a722d1"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.7+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5e77dbf117412d4f164a464d610ee6050cc75272"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.9.6"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "SnoopPrecompile", "StructTypes", "UUIDs"]
git-tree-sha1 = "84b10656a41ef564c39d2d477d7236966d2b5683"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.12.0"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "099e356f267354f46ba65087981a77da23a279b7"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.0"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MarkdownLiteral]]
deps = ["CommonMark", "HypertextLiteral"]
git-tree-sha1 = "0d3fa2dd374934b62ee16a4721fe68c418b92899"
uuid = "736d6165-7244-6769-4267-6b50796e6954"
version = "0.1.1"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Memoize]]
deps = ["MacroTools"]
git-tree-sha1 = "2b1dfcba103de714d31c033b5dacc2e4a12c7caa"
uuid = "c03570c3-d221-55d1-a50c-7939bbd78826"
version = "0.4.4"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1aa4b74f80b01c6bc2b89992b861b5f210e665b5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.21+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "a5aef8d4a6e8d81f171b2bd4be5265b01384c74c"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.10"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ad59edfb711a4751e0b8271454df47f84a47a29e"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.14"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "b478a748be27bd2f2c73a7690da219d0844db305"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.51"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "259e206946c293698122f63e2b513a7c99a244e8"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.ShortCodes]]
deps = ["Base64", "CodecZlib", "HTTP", "JSON3", "Memoize", "UUIDs"]
git-tree-sha1 = "95479a28f0bb4ad37ec7c7ece7fdbfc400c643e0"
uuid = "f62ebe17-55c5-4640-972f-b59c0dd11ccf"
version = "0.3.5"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "75ebe04c5bed70b91614d684259b661c9e6274a4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.0"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "ca4bccb03acf9faaf4137a9abc1881ed1841aa70"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "ba4aa36b2d5c98d6ed1f149da916b3ba46527b2b"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.14.0"

    [deps.Unitful.extensions]
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.7.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─8f199e33-f971-4f96-8883-779f7cc6a946
# ╟─1079052f-234f-4698-9bb8-af49b24ef80a
# ╟─518e6304-8f40-4ef7-925f-e16558080b60
# ╟─cc84fab4-f9fe-498b-8136-9e997b8a6d09
# ╟─90c5e28c-3c05-43a3-bb4c-2f0a47f5f339
# ╟─4432a409-ce78-41bb-863f-dde56f1885f1
# ╟─45c4ea63-050c-44ac-aeb1-0507340f2675
# ╟─a96a12c7-027e-43ba-a4c6-8d013f651899
# ╠═f04ec0c5-fa5a-4467-86d2-9745ffba834a
# ╠═53045f48-b321-4c8c-8fdb-9d1369faa386
# ╠═67b230cd-77e2-4ecd-83cb-47768ff457c6
# ╠═6a8475b7-7ba5-4310-a1f6-d3529cc01651
# ╟─39760bb3-05f8-4a9c-89f2-2fcd46cab544
# ╠═eeb623fd-558d-4c1e-9fe7-677fc704ee91
# ╟─c7cc130f-54c8-4c49-8072-818cc4b7d9b3
# ╠═b65e7452-9ae7-4762-8c3a-e1f0b3f2744b
# ╟─0f84a539-dce0-4dc2-b6ff-85bf759b3c99
# ╠═a36fa201-58da-4e78-8ed2-1127fa730d7a
# ╟─f8ff2918-993a-4527-8ec7-00f96b4bfe13
# ╟─cf37924a-cfe4-4336-a9b3-88c2c6c268f1
# ╠═1f7e8f70-eef4-40da-82a0-7050dec53f76
# ╠═e63f6699-7a70-4bd6-ac4c-46bc406cdc23
# ╟─69f55915-f802-4f58-839f-60897a864c97
# ╠═f3948247-6b50-4329-aef2-e524c0e781ee
# ╠═f11e9ecb-6bbd-468d-a0b5-cdb315c6a8ed
# ╟─827cecdf-1345-4149-aaa0-d752c0dff201
# ╠═c10a9d5b-c92a-452f-8302-5032a9165a06
# ╟─161197f4-1076-4fb0-81de-b34c927ca4c2
# ╠═a4bffbdd-a2a1-42ce-92c9-67332738706d
# ╠═88be3fab-f55b-4ea4-a01c-90c02b7b15a1
# ╠═d57fca19-7a1e-4189-8672-784ed69819ce
# ╟─b2f5f523-cea7-4b2a-8e91-44e466cfa2bb
# ╠═ae2293e5-66e5-4711-b351-a84225570a86
# ╟─918532ec-8353-4966-96f3-1e6126f1e916
# ╠═1a910c0f-1bcc-4755-87a1-63c2496f0551
# ╟─ba7b2db6-8fd5-4050-8983-f5d58147dbd4
# ╟─0e36f365-de6b-4019-8b19-ee77c36dd344
# ╟─c4a319b9-c164-4e6d-afbe-da91b7a95846
# ╟─fcc0b4f7-32bf-454d-8ec8-3a50139c984d
# ╟─1c1cd3ea-0ea5-4637-83c2-63ec10d7cbdd
# ╠═1b8a9dc3-1cbf-4ab0-93c3-d31eb1bf0353
# ╠═9d5135d6-0ba1-435d-884b-a74dbc906d61
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
