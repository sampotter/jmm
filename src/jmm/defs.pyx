from enum import Enum

class Stype(Enum):
    Constant = 0

class State(Enum):
    Far = 0
    Trial = 1
    Valid = 2
    Boundary = 3
    AdjacentToBoundary = 4
    NewValid = 5
