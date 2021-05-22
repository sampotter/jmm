from enum import Enum

class State(Enum):
    Far = 0
    Trial = 1
    Valid = 2
    Boundary = 3
    AdjacentToBoundary = 4
    NewValid = 5
    Shadow = 6

class Ftype(Enum):
    PointSource = 0
    Reflection = 1
    EdgeDiffraction = 2

class Stype(Enum):
    Constant = 0
    NumStype = 1

class Error(Enum):
    Success = 0
    BadArgument = 1
