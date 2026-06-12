import operator

class Selector():

    OPS = {
        '>' : operator.gt,
        '<' : operator.lt,
        '>=': operator.ge,
        '<=': operator.le,
        '==': operator.eq,
    }

    def __init__(self, config):
        self.selection = config.get('selection', [])

    def _apply_cut(self, p, attr, op, threshold, opt=''):
        return Selector._apply_cut_static(p, attr, op, threshold, opt)

    @staticmethod
    def passes_cuts(p, cuts):
        return all(
            Selector._apply_cut_static(p, attr, *cut_def)
            for attr, cut_def in cuts.items()
        )

    @staticmethod
    def _apply_cut_static(p, attr, op, threshold, opt=''):
        value = getattr(p, attr)
        if opt == '||':
            value = abs(value)
        return Selector.OPS[op](value, threshold)

    def applySelection(self, particles):
        # One counter for each category
        counts = {cat['name']: 0 for cat in self.selection}

        # Single pass over particles
        for p in particles:
            for cat in self.selection:
                if self.passes_cuts(p, cat['cuts']):
                    counts[cat['name']] += 1

        # Check if the collection passes the selection
        event_ok = True
        for cat in self.selection:
            name  = cat['name']
            min_p = cat.get('min_particles', 0)
            max_p = cat.get('max_particles', float('inf'))
            n     = counts[name]

            cat_ok   = min_p <= n <= max_p
            event_ok = event_ok and cat_ok

        return event_ok
