"""Tests of the numerical core. No data needed.

python tests/test_core.py         every test prints ok
pytest tests                      the same with pytest
"""

import math
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from rf_ophiolite import originals
from rf_ophiolite.common import is_inside, wrap180
from rf_ophiolite.imaging import moveout_stack, ray_tables
from rf_ophiolite.receiver_functions import hk_best, hk_stack, itdecon, ps_depth
from rf_ophiolite.station_terms import _lsq_terms
from rf_ophiolite.synthetic import phase_times, synthetic_pair


def _pulse(n, dt, t0, f0=0.7):
    t = np.arange(n) * dt
    arg = (math.pi * f0 * (t - t0)) ** 2
    return (1.0 - 2.0 * arg) * np.exp(-arg)


def test_deconvolution_recovers_two_spikes():
    """R = 0.3 Z + 0.1 Z delayed by 4 s gives spikes of 0.3 at 0 s and 0.1 at 4 s."""
    dt, n = 0.05, 2000
    z = _pulse(n, dt, 10.0)
    r = 0.3 * z + 0.1 * _pulse(n, dt, 14.0)
    rf, t, fit = itdecon(r, z, dt, a=2.0, shift=10.0)
    assert fit > 99.0
    assert abs(rf[np.argmin(np.abs(t))] - 0.3) < 0.01
    assert abs(rf[np.argmin(np.abs(t - 4.0))] - 0.1) < 0.01
    assert abs(t[np.argmax(rf)]) < dt


def _synthetic_station(h=35.0, kappa=1.75, vp=6.3, n=25, noise=0.02, seed=3):
    rng = np.random.RandomState(seed)
    dt = 0.05
    rays = np.sort(rng.uniform(0.04, 0.08, n))
    rfs = []
    for p in rays:
        z, r = synthetic_pair(h, kappa, vp, p, dt, rng, noise)
        rf, t, fit = itdecon(r, z, dt, a=2.0, shift=10.0)
        rfs.append(rf[(t >= -5.0) & (t <= 60.0)])
    return np.vstack(rfs), -5.0, dt, rays


def test_hk_stack_recovers_the_crust():
    R, t0, dt, rays = _synthetic_station()
    H = np.arange(20.0, 60.0 + 1e-9, 0.25)
    K = np.arange(1.60, 2.00 + 1e-9, 0.005)
    h, k, _ = hk_best(hk_stack(R, t0, dt, rays, 6.3, H, K), H, K)
    assert abs(h - 35.0) <= 1.0
    assert abs(k - 1.75) <= 0.03


def test_ps_depth_and_moveout():
    R, t0, dt, rays = _synthetic_station()
    t_ps, h, err = ps_depth(R, t0, dt, rays, 6.3, 1.75)
    assert abs(h - 35.0) <= 1.0
    assert err > 0
    # the moveout-corrected stack has its Ps peak at the time of the reference ray parameter
    tt, st = moveout_stack(R, t0, dt, rays, vp=6.3, kappa=1.75, pref=0.06)
    sel = (tt > 2.0) & (tt < 8.0)
    assert abs(tt[sel][np.argmax(st[sel])] - phase_times(35.0, 1.75, 6.3, 0.06)[0]) < 0.15


def test_station_terms_are_recovered():
    """Residuals built as station term + event term + noise give the station terms back."""
    rng = np.random.RandomState(0)
    n_sta, n_evt = 12, 60
    s_true = rng.normal(0.0, 0.3, n_sta)
    s_true -= s_true.mean()
    e_true = rng.normal(0.0, 1.0, n_evt)
    si, ei = np.meshgrid(np.arange(n_sta), np.arange(n_evt), indexing='ij')
    keep = rng.uniform(size=si.size) < 0.7  # 30 % of the picks are missing
    si, ei = si.ravel()[keep], ei.ravel()[keep]
    r = s_true[si] + e_true[ei] + rng.normal(0.0, 0.02, len(si))
    s, e = _lsq_terms(si, ei, r, n_sta, n_evt)
    assert abs(s.sum()) < 1e-6
    assert np.max(np.abs(s - s_true)) < 0.03


def test_ray_tables_constant_velocity():
    z = np.arange(0.0, 60.0 + 1e-9, 0.5)
    vp, vs, p = 6.3, 3.6, 0.06
    tps, xs = ray_tables(p, z, np.full(len(z), vp), np.full(len(z), vs))
    qp, qs = math.sqrt(1 / vp**2 - p**2), math.sqrt(1 / vs**2 - p**2)
    assert np.allclose(tps, z * (qs - qp))
    assert np.allclose(xs, z * p / qs)


def test_small_helpers():
    assert np.allclose(wrap180([190.0, -190.0, 360.0, 45.0]), [-170.0, 170.0, 0.0, 45.0])
    with tempfile.TemporaryDirectory() as d:
        sub = os.path.join(d, 'a', 'b')
        os.makedirs(sub)
        assert is_inside(sub, d)
        assert is_inside(d, d)
        assert not is_inside(d, sub)


def test_originals_manifest_detects_a_changed_file():
    with tempfile.TemporaryDirectory() as d:
        data = os.path.join(d, 'raw', 'Event_1')
        os.makedirs(data)
        for name in ('a.sac', 'b.sac'):
            with open(os.path.join(data, name), 'wb') as fh:
                fh.write(os.urandom(4096))
        manifest = os.path.join(d, 'originals_sha256.txt')
        originals.make([os.path.join(d, 'raw')], manifest, 1, False)
        assert originals.verify(manifest, 1) == 0
        with open(os.path.join(data, 'a.sac'), 'r+b') as fh:  # change one byte
            fh.seek(100)
            fh.write(b'\x00' if fh.read(1) != b'\x00' else b'\x01')
        assert originals.verify(manifest, 1) == 1


if __name__ == '__main__':
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            fn()
            print('ok  ', name)
    print('all tests passed')
