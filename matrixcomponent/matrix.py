"""Python object models to be manipulated"""

from dataclasses import dataclass
from typing import List
from sortedcontainers import SortedDict

## Path is all for input files


class Path:
    name: str
    bins: 'SortedDict[int, Path.Bin]'
    links: 'numpy.array'

    def __init__(self, name=''):
        self.name = name
        self.bins = SortedDict()

    @dataclass
    class Bin:
        bin_id: int
        coverage: float
        inversion_rate: float
        first_nucleotide: int
        last_nucleotide: int
        sequence: str = ''

    class LinkEntry:
        def __init__(self, upstream, downstream):
            self.upstream = upstream
            self.downstream = downstream
            # TODO: self.insert_size will require a topology search to find this

    def __contains__(self, item):  # used by " x in Path "
        return item in self.bins

## For Output to RDF  ###########
@dataclass
class LinkColumn:
    upstream: int
    downstream: int
    participants: List[bool]  # in order path_names, true if the individual participates in this LinkColumn
    # participants depends on row ordering of path names, optimized precompute for display


@dataclass
class Bin:
    coverage: float
    inversion: float
    first_nucleotide: int
    last_nucleotide: int


class Component:
    """Block of co-linear variation within a Graph Matrix
        # bin_id and seq are global to column and could be reduced to save memory,
        # careful construction can reuse Bin.sequence memory pointer"""
    first_bin: int
    last_bin: int
    occupants: List[bool]
    matrix: List[List[Bin]]
    arrivals: List[LinkColumn]
    departures: List[LinkColumn]

    def __init__(self, first_bin: int, last_bin: int):
        self.first_bin = first_bin
        self.last_bin = last_bin
        self.occupants = []
        self.matrix = []
        self.arrivals = []  # reverse ordered Links
        self.departures = []  # ordered Links

    def column_count(self):
        """Used to estimate JSON size.  LinkColumns are counted twice because they have a
        participants boolean list."""
        return 2*(len(self.arrivals) + len(self.departures)) + self.last_bin - self.first_bin
