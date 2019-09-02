from SeidrFile.libseidrfile import SeidrFile_Internal

class SeidrFile(object):
  """docstring for SeidrFile"""
  def __init__(self, path):
    super(SeidrFile, self).__init__()
    self.path = path
    self.sf_ptr = SeidrFile_Internal(path)
    self.attr = self.sf_ptr.attr()
    self.nodes = self.sf_ptr.nodes()
    self.supplementary = self.sf_ptr.supplementary()
    self.algorithms = self.sf_ptr.algorithms()

  def algorithms(self):
    return self.algorithms

  def attr(self):
    return self.attr

  def supplementary(self):
    return self.supplementary

  def nodes(self):
    return self.nodes

  def betweenness(self):
    return self.sf_ptr.betweenness()

  def closeness(self):
    return self.sf_ptr.closeness()

  def pagerank(self):
    return self.sf_ptr.pagerank()

  def eigenvector(self):
    return self.sf_ptr.eigenvector()

  def strength(self):
    return self.sf_ptr.strength()

  def katz(self):
    return self.sf_ptr.katz()

  def edge_scores(self):
    return self.sf_ptr.edge_scores()

  def edge_ranks(self):
    return self.sf_ptr.edge_ranks()

  def edge_supplementary(self):
    return self.sf_ptr.edge_supplementary()

  def edge_index(self):
    return self.sf_ptr.edge_index()

  def edge_flag(self):
    return self.sf_ptr.edge_flag()

  def edge(self):
    return self.sf_ptr.edge()

  def next_edge(self):
    return self.sf_ptr.next_edge()

  def __reset__(self):
    self.sf_ptr.reset()

class edgeIterator(object):
  """docstring for edgeIterator"""
  def __init__(self, seidrfile):
    super(edgeIterator, self).__init__()
    self.seidrfile = seidrfile

  def __len__(self):
    return self.seidrfile.attr.edges

  def __iter__(self):
    return self

  def __next__(self):
    if not self.seidrfile.next_edge():
      raise StopIteration
    return self.seidrfile.edge()

class scoreIterator(object):
  """docstring for scoreIterator"""
  def __init__(self, seidrfile):
    super(scoreIterator, self).__init__()
    self.seidrfile = seidrfile

  def __len__(self):
    return self.seidrfile.attr.edges

  def __iter__(self):
    return self

  def __next__(self):
    if not self.seidrfile.next_edge():
      raise StopIteration
    return self.seidrfile.edge_scores()

class rankIterator(object):
  """docstring for rankIterator"""
  def __init__(self, seidrfile):
    super(rankIterator, self).__init__()
    self.seidrfile = seidrfile

  def __len__(self):
    return self.seidrfile.attr.edges

  def __iter__(self):
    return self

  def __next__(self):
    if not self.seidrfile.next_edge():
      raise StopIteration
    return self.seidrfile.edge_ranks()

class supplementaryIterator(object):
  """docstring for supplementaryIterator"""
  def __init__(self, seidrfile):
    super(supplementaryIterator, self).__init__()
    self.seidrfile = seidrfile

  def __len__(self):
    return self.seidrfile.attr.edges

  def __iter__(self):
    return self

  def __next__(self):
    if not self.seidrfile.next_edge():
      raise StopIteration
    return self.seidrfile.edge_supplementary()

class flagIterator(object):
  """docstring for flagIterator"""
  def __init__(self, seidrfile):
    super(flagIterator, self).__init__()
    self.seidrfile = seidrfile

  def __len__(self):
    return self.seidrfile.attr.edges

  def __iter__(self):
    return self

  def __next__(self):
    if not self.seidrfile.next_edge():
      raise StopIteration
    return self.seidrfile.edge_flag()

class indexIterator(object):
  """docstring for indexIterator"""
  def __init__(self, seidrfile):
    super(indexIterator, self).__init__()
    self.seidrfile = seidrfile

  def __len__(self):
    return self.seidrfile.attr.edges

  def __iter__(self):
    return self

  def __next__(self):
    if not self.seidrfile.next_edge():
      raise StopIteration
    return self.seidrfile.edge_index()
    