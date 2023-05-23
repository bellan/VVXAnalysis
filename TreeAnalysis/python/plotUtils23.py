import ROOT

class TFileContext(object):
    def __init__(self, *args):
        # print('>>>Opening with args:', args)
        self.tfile = ROOT.TFile(*args)
    def __enter__(self):
        return self.tfile
    def __exit__(self, type, value, traceback):
        # print('<<<Closing TFile "%s"' % (self.tfile.GetName()))
        self.tfile.Close()


def addIfExisting(*args):
    result = None
    for a in [ a for a in args if a is not None ]:
        if result is None:
            result = a
        else:
            result.Add(a)
    return result
