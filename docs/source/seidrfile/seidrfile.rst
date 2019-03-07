.. _seidrfile-label:

SeidrFiles
=================================

Introduction
^^^^^^^^^^^^

Seidr employs it's own file format (called SeidrFile) to store network data. This is done to increase performance, as SeidrFiles are:

* Losslessly compressed using bgzip (to save space)
* Ordered in a lower triangular to enable faster algorithms
* Ranked, so that scores can be rank-aggregated 

The SeidrFile header
^^^^^^^^^^^^^^^^^^^^

A SeidrFile has a **header** that keeps information such as the number of edges, nodes, the node names etc. You can view the header of a SeidrFile with the command::

    seidr view -H <SeidrFile>

The output might look something like this::

    # [G] Nodes: 50
    # [G] Edges: 1225
    # [G] Storage: Dense
    # [G] Algorithms #: 14
    # [G] Supplementary data #: 13
    # [A] ARACNE
    # [A] CLR
    # [A] ELNET
    # [A] MI
    # [A] GENIE3
    # [A] LLR
    # [A] NARROMI
    # [A] PCOR
    # [A] PEARSON
    # [A] PLSNET
    # [A] SPEARMAN
    # [A] SVM
    # [A] TIGRESS
    # [A] irp
    # [S] D1
    # [S] D2
    # [S] D3
    # [S] D4
    # [S] D5
    # [S] D6
    # [S] D7
    # [S] D8
    # [S] D9
    # [S] D10
    # [S] D11
    # [S] D12
    # [S] D13
    # [R] Version: 0.10.0
    # [R] Cmd: seidr aggregate -f -k -m irp aracne.sf clr.sf elnet.sf elranks.sf genie3.sf llr.sf narromi.sf pcor.sf pearson.sf plsnet.sf spearman.sf svm.sf tigress.sf
    # [N] G1
    # [N] G2
    # [N] G3
    # [N] G4
    # [N] G5
    # [N] G6

The SeidrFile body
^^^^^^^^^^^^^^^^^^

In the main body of a SeidrFile, we store the edges of a network. Specifically, for each edge, we have at least four columns:

* Source: For directed edges, this is the originating node, for undirected edges, this is simply one of the partners
* Target: For directed edges, this is the destination node, for undirected edges, this is simply the other partner
* Type: Undirected if the node is undirected, Directed otherwise
* X\_score;X\_rank: This column holds the original score for algorithm "X" as well as its computed rank.

Besides these four mandatory columns, a SeidrFile can hold any number of additional score/rank columns if it is an aggregated or otherwise processed file and and additional supplementary column that annotates the edge with extra information. To view the body of a SeidrFile you can use::

    seidr view <SeidrFile>

Here is the output of a simple imported network::

    G1      G2      Directed        0.004;334084
    G3      G1      Directed        0.334;22729.5
    G1      G4      Directed        0.071;89307
    G4      G2      Directed        0.053;104778
    G3      G4      Directed        0.006;282776


And one that is a little more complex, with 14 score/rank columns and a supplementary column at the end. In aggregated SeidrFiles, the representative score/rank is always the rightmost (last) score/rank column::

    G2  G1  Directed  0.288087;1.30856e+06  nan;nan 1.87357;106802  0.004;334084  -0.018736;243746  0.0904447;42007 0.244;37455.5 0.0128741;202752  -0.159435;202751  1.07712e-05;360264  -0.00225177;1.32058e+06 0.152;26168 nan;nan 0.978291;117022 11

You might notice the columns with nan:nan as score/rank. Seidr uses nan as a placeholder to denote a missing edge. That means this particular edge (G2 -> G1) was not found in the second and thirteenth algorithms.

The SeidrFile index
^^^^^^^^^^^^^^^^^^^

SeidrFiles can be indexed with the command::

    seidr index <SeidrFile>

This will create an index file with the extension .sfi. The index allows us to access edges quickly in a SeidrFile without having to decompress unnecessary data. Some seidr commands therefore need the index. As an example, let's see what happens if we try to pull out a specific edge from a SeidrFile without an index::

    seidr view -n G1000:G3 <SeidrFile>
    [ ERROR   ][ 2018-05-02T21:35:45 ][ seidr ]: SeidrFile index <SeidrFile.sfi> must exist when using --nodelist

Otherwise, if the index exists::

    seidr view -n G1000:G3 ../dream_net1/aggregate/aggregated.sf
    G1000 G3  Undirected  0.388607;611152 nan;nan nan;nan 0.001;581639  -0.0200038;209560 0.00623208;1.16541e+06  0.057;174410  0.00177422;752791 -0.0595161;752789 2.76065e-06;1.11154e+06 -0.0432047;834369 0.031;315583  0.0006;123144 0.507107;458113


