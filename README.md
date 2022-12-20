# crawlerPlusPlus
C++ version of my mini simulation engine originally written in Julia
I'm writing this to learn more C++ and to practice scientific computing, which I love to do.
I'm purposely restricting myself to the stdlib both to simplify my C++ journey and
to force me to code certain technical things myself (like matrix inversion with Gaussian
elimination algorithm, which is needed for M-SHAKE constraints) rather than relying on libraries.

Everything in this engine is in SI units, including time and everything else.
I tried to make this as OOP as possible/helpful. Writing this made me appreciate OpenMM's object
heirarchy model, which I initially didn't understand/like as a user of the tool. Now it really
makes sense the way they've done things, in terms of being able to reason about the code and maintain it.
For example, I initially didn't understand why Simulation was a class, like why not just have a function
runSimulation() or something for such a high level construct? Now I appreciate why they did it this way,
honestly. It makes sense. I followed the OpenMM model mostly here, except I skipped the System class and
go directly from Topology, Parameters, Coordinates, Velocities, and Forces to Simulation. For such a small project this is fine.
