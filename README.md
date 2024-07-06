# Component Mode Synthesis

Do you want to run a structural dynamics simulation with full order accuracy, even if it includes some negligible details and takes a long time to complete?

Or would you prefer to simplify the model, achieving high enough accuracy while significantly reducing computation time? If you're lucky, you won't lose any critical information either. 😄

In this repo, I'd like to introduce some model reduction methods for finite element analysis (aka reduced order modeling). Specifically, in the CMS explanation, I've made some modifications that you might easily imagine but take some effort to implement. Don't worry, I've done the hard work for you.

Also, I very strongly recommend you to use conventional finite element solvers, such as Ansys, Adina, Abaqus, Fenics, etc. I used MATLAB's PDE solver, but it is used just for DEMONSTRATION.

## Applications
Dynamic and deformable body

Supramolecular Protein Dynamics


## Theories

1. **Guyan Reduction**
   
   Guyan reduction, also known as static condensation, is a method used to reduce the size of a finite element model by eliminating certain degrees of freedom. This is achieved by assuming that the displacements in the eliminated degrees of freedom can be approximated as linear combinations of the retained degrees of freedom. This method is particularly effective for static or low-frequency dynamic analyses, where higher modes have negligible impact on the response of interest. The primary advantage of Guyan reduction is its simplicity and computational efficiency.

2. **Craig-Bampton Method**
   
   The Craig-Bampton method, also known as Component Mode Synthesis (CMS), is a technique used to reduce the computational cost of dynamic analysis of large structures. This method involves partitioning the structure into substructures, each of which is analyzed independently. The dynamic behavior of each substructure is represented using a combination of fixed-interface normal modes and constraint modes. These modes are then assembled to approximate the global behavior of the entire structure. The Craig-Bampton method is particularly useful for problems involving complex assemblies or where components are reused in different assemblies, as it allows for efficient re-analysis of the system.

Enjoy!
