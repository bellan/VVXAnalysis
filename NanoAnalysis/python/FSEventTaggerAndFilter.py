from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

from VVXAnalysis.NanoAnalysis.Selector import Selector
from VVXAnalysis.NanoAnalysis.Regions import Flags


class FSEventTaggerAndFilter(Module):

    COLLECTIONS = ('leptons', 'photons', 'jets')

    def __init__(self, flags, regions):
        self.regions = regions
        self._plan   = self._build_plan(flags)

    # ------------------------------------------------------------------
    # Plan construction (called once at initialisation)
    # ------------------------------------------------------------------

    def _cuts_signature(self, cuts):
        """Return a hashable key that uniquely identifies a cuts dict,
        ignoring name and min/max_particles.  Two categories that share
        the same signature can be counted with a single particle loop."""
        return tuple(sorted(
            (attr, op, thr, rest[0] if rest else '')
            for attr, (op, thr, *rest) in cuts.items()
        ))

    def _build_plan(self, flags):
        """Pre-compute, for each collection, the minimal set of unique
        cut signatures and the per-flag requirements expressed as
        (counter_index, min_particles, max_particles) triples.

        At analysis time only one particle loop per collection is needed,
        regardless of how many flags share the same cut signature.
        """
        plan = {}

        for coll in self.COLLECTIONS:
            signature_to_id = {}   # cuts signature  -> counter index
            counter_defs    = []   # cuts dict for each unique signature
            flag_reqs       = []   # one entry per flag

            for flag in flags:
                cfg = flag.get(coll)

                if cfg is None:
                    # This flag has no requirement on this collection
                    flag_reqs.append({'flag_name': flag['name'], 'checks': None})
                    continue

                checks = []
                for cat in cfg['selection']:
                    sig = self._cuts_signature(cat['cuts'])

                    if sig not in signature_to_id:
                        signature_to_id[sig] = len(counter_defs)
                        counter_defs.append(cat['cuts'])

                    counter_id = signature_to_id[sig]
                    checks.append((
                        counter_id,
                        cat.get('min_particles', 0),
                        cat.get('max_particles', float('inf'))
                    ))

                flag_reqs.append({'flag_name': flag['name'], 'checks': checks})

            plan[coll] = {
                'counter_defs': counter_defs,  # unique cuts to evaluate
                'flag_reqs'   : flag_reqs,
            }

        return plan

    # ------------------------------------------------------------------
    # Per-event helpers
    # ------------------------------------------------------------------

    def _count_unique(self, particles, counter_defs):
        """Single pass over particles: one counter per unique cut signature.
        Returns a list of integers aligned with counter_defs."""
        counts = [0] * len(counter_defs)

        for p in particles:
            for i, cuts in enumerate(counter_defs):
                if Selector.passes_cuts(p, cuts):
                    counts[i] += 1

        return counts

    def _build_region_word(self, counts_by_coll):
        """Evaluate all flag requirements against pre-computed counts.
        No particle loop here — only integer comparisons."""
        region_word = 0

        for coll, coll_plan in self._plan.items():
            coll_counts = counts_by_coll[coll]

            for freq in coll_plan['flag_reqs']:
                if freq['checks'] is None:
                    continue

                flag_name = freq['flag_name']
                if flag_name not in Flags.__members__:
                    continue

                if all(min_p <= coll_counts[cid] <= max_p
                       for cid, min_p, max_p in freq['checks']):
                    region_word |= Flags[flag_name]

        return region_word

    # ------------------------------------------------------------------
    # Module interface
    # ------------------------------------------------------------------

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("regionWord", "I",
                        title="Word that contains the regions that passed the selection")

    def analyze(self, event):
        """Process event; return True (keep) or False (discard)."""

        raw = {
            'leptons': Collection(event, "Lepton"),
            'photons': Collection(event, "Photon"),
            'jets'   : Collection(event, "Jet"),    # FIXME: AK4 only for now
        }

        # One particle loop per collection over the unique cut signatures
        counts_by_coll = {
            coll: self._count_unique(raw[coll], self._plan[coll]['counter_defs'])
            for coll in self.COLLECTIONS
        }

        region_word = self._build_region_word(counts_by_coll)

        self.out.fillBranch("regionWord", region_word)

        for region in self.regions:
            if region_word & region == region:
                return True

        return False
