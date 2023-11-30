# 1D_1Group_Diffusion
This is a very primitive (and buggy) project that aims to simulate monoenergetic steady-state neutron diffusion conditions in 1D.

Everything is contained in a single python file. All submethods (gauss-seidel solver, eigenvalue solver, matrix initializer) are relatively short in code length.

You will want to adjust some early parameters depending on the type of problem you're looking to solve.

First, is your diffusion coefficient. This describes how far a typical thermal neutron will travel.given a material density and macroscopic cross section. A great resource for this is https://www.nuclear-power.com/nuclear-power/reactor-physics/neutron-diffusion-theory/diffusion-coefficient/

Next are the cross sections themselves. Right now, I just have two options to allow for heterogeneity in the system. Sigma_A and Sigma_B are what you will want to change depending on your problem. However these are currently dimensionless (assumed barns) but right now I'm working with arbitrary and easy numbers to maximize code functionality.

Sigma_F may also be changed. For all tests, I've left it unchanged. This describes the fission cross section variable.

You'll note a special constant as well -- NU. This is the average number of neutrons per fission (around 2.5). I use this for my F matrix as will be described later.

Material boundary describes where on the 1D x-axis the material changes. Functionality is uncertain as of now. You can change this value and you'll get varying results but accuracy has not been fully tested yet.

L and N are important for this code. L represents the total length of cells with aribitrary units...again for computational ease. N describes the discretization of your system or in easier terms--your resolution. Feel free to bump it up to higher values but the system may break.

There is finally an 'h' term. DO NOT change this. It's a simple tool that calculates mesh spacing for the matrix creation function. 


# SubModules And Their Function

create_matrix_A()

create_matrix_F()

gauss-seidel()

eigenvalue_solve()

plot_flux_distributions()

