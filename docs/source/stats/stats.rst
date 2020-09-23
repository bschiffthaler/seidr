.. _stats-label:

Getting network statistics
==============================

Centrality
----------


How scores are used in ``seidr stats``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We typically use scores as measures of similarity in ``seidr`` workflows. This means that *higher is better*. As an example in this network::

  A  B  1
  A  C  0.5

the edge A<->B is *stronger* than A<->C. In centrality metrics we often use weights as either similarity, or distance. When we e.g., calculate betweenness centrality, we want to know the shortest path from A to B, therefore weights are usually interpreted as distances here, and therefor *lower is better*.

By default ``seidr`` assumes that weights are *similarities* and handles them as such. When sensible, it will use :math:`\frac{1}{w}`. If your data represents a distance, you must use the flag ``--weight-is-distance``, otherwise your outcome will be wrong. If you set this flag, seidr will use :math:`\frac{1}{w}` for calculations where it assumes the weight indicates a *similarity* (i.e. the behaviour is inverse). See metrics below where the similarity [S] and distance [D] metrics are indicated.

Metrics
^^^^^^^

``seidr`` can calculate a limited number of network centrality statistics on ``SeidrFiles``.

On any ``SeidrFile`` you can run::

  seidr stats seidrfile.sf

to calculate the network centrality statistics. By default all metrics that can be calculated will be. Use the ``-m,--metrics`` option to control this (see above as to the meaning of [S] and [D]):

**For nodes**

* ``PR`` - PageRank [S]
* ``CLO`` - Closeness centrality [D]
* ``BTW`` - Betweenness centrality [D]
* ``STR`` - Strength (weighted degree) centrality [S]
* ``EV`` - Eigenvector centrality [S]
* ``KTZ`` - Katz centrality [S]
* ``LPC`` - Laplacian centrality [S]

**For edges**

* ``SEC`` - Spanning edge centrality [D]
* ``EBC`` - Edge betweenness centrality [D]

To select only few of these run e.g.::

  stats stats -m BTW,CLO seidrfile.sf

Approximate vs exact
^^^^^^^^^^^^^^^^^^^^

By default, ``seidr`` uses approximations where it can to compute centrality statistics. It will sample ``-n,--nsamples`` nodes to do so. If not specified, that number is 10% of nodes. If your network is small, you can turn on exact metrics with ``-e,--exact``.

Viewing stats
^^^^^^^^^^^^^

You can view node level statistics with::

  seidr view --centrality seidrfile.sf

Edge level statistics are stored as edge attributes. You can add tags to see which attributes correspond to which stat::

  seidr view -a seidrfile.sf