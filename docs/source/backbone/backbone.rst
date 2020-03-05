.. _backbone-label:

Calculating a Network Backbone
==============================

``seidr`` implements Coscia and Neffke's backboning algorithm (very neatly described `here <http://www.michelecoscia.com/?p=1236>`_ if you don't feel like handling a lot of math, otherwise here: [Coscia2017]_).

On any ``SeidrFile`` you can run::

  seidr backbone seidrfile.sf

to calculate the network backbone statistics. Not that we are not cutting edges yet. To do that, we need to specify a measure of standard deviations to cut. Essentially, we want to define how extreme an edge has to deviate from its expected value, so that we keep it, the higher, the more stringent. A conservative value, would be 1.28, which corresponds approxiamtely to a P-value of 0.1::

  seidr backbone -F 1.28 -o seidrfile.bb.1.28.sf seidrfile.bb.sf