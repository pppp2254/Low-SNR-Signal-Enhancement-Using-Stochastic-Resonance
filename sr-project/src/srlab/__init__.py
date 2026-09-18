"""srlab: stochastic resonance simulation library.

All reusable code lives here; experiment scripts in experiments/ are thin and
import from this package (PLAN.md section 3). All randomness goes through
numpy.random.default_rng(seed).
"""

__version__ = "0.1.0"
