# node-to-arc-centric-dbg
Convert a node-centric DBG into an arc-centric one.
Currently the input is the BCALM2 fasta format, and the output is an edge list preceded by the number of nodes.
The arc-centric dbg is represented in "doubled" form, i.e. a directed graph in which each binode is represented by a pair of nodes and each biarc is represented by a pair of arcs.
As an exception, self-complemental arcs are collapsed into a single arc.
The edge list has the columns `<node1> <node2> <weight> <mirror_node1> <mirror_node2> <sequence>`.
The mirror nodes are the nodes corresponding to the reverse complement of an arc.
Note that there may be parallel arcs, so this is not enough to identify the reverse complement arc.