# basic objects for calculations

from xml.etree.ElementTree import tostring


class Range:

    def __init__(self):
        self.start  = 0
        self.end    = 0

    def __init__(self, start: float, end: float):
        self.start  = start
        self.end    = end

    def to_str(self) -> str:
        return "(" + str(self.start) + ", " + str(self.end) + ")"

    def contains(self, x: float) -> bool:
        return (x >= self.start) and (x <= self.end)
