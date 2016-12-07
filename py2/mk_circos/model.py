# -*- coding: utf-8 -*-


class Chr:
    """karyotype.txt's chr field."""
    def __init__(self, id_, label, start, end, color):
        self.id_   = id_
        self.label = label
        self.start = start
        self.end   = end
        self.color = color
        self._str  = "\t".join(['chr', '-', id_, label, start, end, color])

    def __str__(self):
        return self._str


class Band:
    """karyotype.txt's band field."""
    def __init__(self, domain, id_, label, strat, end, color):
        self.domain = domain
        self.id_    = id_
        self.label  = label
        self.strat  = strat
        self.end    = end
        self.color  = color
        self._str = "\t".join(['band', domain, id_, label, strat, end, color])

    def __str__(self):
        return self._str


class Link:
    """item of links.txt"""
    def __init__(self, domain1, start1, end1, domain2, start2, end2, link_type):
        self.domain1 = domain1
        self.start1 = start1
        self.end1 = end1
        self.domain2 = domain2
        self.start2 = start2
        self.end2 = end2
        self.type = link_type

    def __str__(self):
        return "\t".join([getattr(self, i) for i in
            ["domain1", "start1", "end1", "domain2", "start2", "end2"]])


class Label:
    """item of band_labels.txt"""
    def __init__(self, domain, start, end, name):
        self.domain = domain
        self.start = start
        self.end = end
        self.name = name

    def __str__(self):
        return "\t".join([getattr(self, i) for i in 
            "domain,start,end,name".split(',')])
