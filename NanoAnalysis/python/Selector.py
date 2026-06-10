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
        value = getattr(p, attr)
        use_abs = opt == '||'
        if use_abs:
            value = abs(value)
        return Selector.OPS[op](value, threshold)

    def _passes_cuts(self, p, cuts):
        return all(
            self._apply_cut(p, attr, *cut_def)
            for attr, cut_def in cuts.items()
        )

    def applySelection(self, particles):
        # one counter for each category
        counts = {cat['name']: 0 for cat in self.selection}

        # To reduce complexity, just one run over particles
        for p in particles:
            for cat in self.selection:
                if self._passes_cuts(p, cat['cuts']):
                    counts[cat['name']] += 1

        # Now check if the collection passes the selection
        event_ok = True
        for cat in self.selection:
            name  = cat['name']
            min_p = cat.get('min_particles', 0)
            max_p = cat.get('max_particles', float('inf'))
            n     = counts[name]

            cat_ok   = min_p <= n <= max_p
            event_ok = event_ok and cat_ok

        return event_ok

