# LuaICP
Iterative Closest Point algorothm implemented in Lua.
This program implements Iterative Closest Point algorithm for alignment of two point clouds. It can be used to align astronomical star sets - the one get from standard astronomical catalog and one from arbitrary source, for example stars extracted from astronomical image.

# Source doc
https://alekmarinov.github.io/luaicp/

# Related scripts

```sim-no-weights.cmd```

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

# Demo

![ICP in action](docs/icp.gif?raw=true "ICP demo")

# Acknowledgment
The program is supported by the project "Astroinformatics" with Bulgarian National Science Foundation, DO-02-275
Institute of Information Technologies - Bulgarian Academy of Sciences, http://www.iit.bas.bg/PR_en.html
