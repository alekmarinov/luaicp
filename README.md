# LuaICP
Iterative Closest Point algorothm implemented in Lua.
This program is an experiment of applying ICP algorithm for aligning stars in astronomical images.

# Related scripts

```sim-no-weights.```
Simulation with random points without weights

```sim-weights.cmd```
Simulation with random points with assigned weights

```ROZ200_001906.cmd```
Simulation with 100 stars extracted from ROZ200_001906 plate

```icpsim.cmd```
Runner script for the ICP simulation program

```icp.cmd```
Runner script of the ICP program

# Applications

```icpsim.cmd -i```
Runs the simulation program in interactive mode with 100 random points

```icp.cmd data\ROZ200_001906_normals.txt data\ROZ200_001906_noised.txt```
Aligns ROZ200_001906_noised and ROZ200_001906_normals sets

```icpsim.cmd --help```
Show program usage options

## More help:
Pressing 'h' from icpsim in interactive mode will list various control keys

### Even more help:
Read program code. It is pure Lua text.
