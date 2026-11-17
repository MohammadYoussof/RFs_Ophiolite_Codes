"""Synthetic test of the deconvolution and the H-k stack. No data needed.

A one-layer crust over a half-space. For a set of ray parameters the vertical
component is a band-limited P pulse with noise, and the radial component is built
from delayed copies of it: the direct P, the Ps conversion at the Moho and the two
crustal multiples (PpPs, and PpSs+PsPs with negative polarity). The receiver
functions are recovered with the iterative deconvolution and stacked in H-k, and the
result is compared with the crust that was put in.

    python rfo.py demo --out demo_output

Written to --out: synthetic_test.png, and a line PASSED or FAILED on the screen.
"""

import math
import os
import sys

import numpy as np

from .common import pyplot, save_figure
from .receiver_functions import hk_best, hk_stack_each, itdecon, ps_depth, rf_quality


def add_arguments(p):
    p.add_argument('--out', default='demo_output', help='output folder')
    p.add_argument('--h', type=float, default=35.0, help='Moho depth of the synthetic crust (km)')
    p.add_argument('--kappa', type=float, default=1.75, help='Vp/Vs of the synthetic crust')
    p.add_argument('--vp', type=float, default=6.3, help='crustal Vp (km/s)')
    p.add_argument('--n', type=int, default=40, help='number of synthetic events')
    p.add_argument('--noise', type=float, default=0.03, help='noise level relative to the P pulse')
    p.add_argument('--gauss', type=float, default=2.0, help='Gaussian width parameter of the deconvolution')
    p.add_argument('--seed', type=int, default=1)


def phase_times(h, kappa, vp, p):
    """Arrival times after the direct P of Ps, PpPs and PpSs+PsPs for a layer of thickness h."""
    a = math.sqrt(1.0 / vp**2 - p**2)
    b = math.sqrt((kappa / vp) ** 2 - p**2)
    return h * (b - a), h * (b + a), 2.0 * h * b


def _delay(x, shift, dt):
    """x delayed by shift seconds (Fourier shift, any fraction of a sample)."""
    f = np.fft.rfftfreq(len(x), dt)
    return np.fft.irfft(np.fft.rfft(x) * np.exp(-2j * np.pi * f * shift), len(x))


def synthetic_pair(h, kappa, vp, p, dt, rng, noise, t_p=10.0, length=100.0):
    """Vertical and radial component of one event; the P pulse arrives at t_p."""
    n = int(round(length / dt))
    t = np.arange(n) * dt
    # source: a few pulses of random size within 6 s, as a teleseismic P wave train
    z = np.zeros(n)
    for _ in range(4):
        t0 = t_p + rng.uniform(0.0, 6.0)
        f0 = rng.uniform(0.4, 0.9)
        arg = (math.pi * f0 * (t - t0)) ** 2
        z += rng.uniform(0.4, 1.0) * rng.choice([-1.0, 1.0]) * (1.0 - 2.0 * arg) * np.exp(-arg)
    z /= np.max(np.abs(z))
    t_ps, t_ppps, t_ppss = phase_times(h, kappa, vp, p)
    r = 0.30 * z + 0.12 * _delay(z, t_ps, dt) + 0.07 * _delay(z, t_ppps, dt) - 0.06 * _delay(z, t_ppss, dt)
    z = z + noise * rng.standard_normal(n)
    r = r + 0.3 * noise * rng.standard_normal(n)
    return z, r


def run(a):
    plt = pyplot()
    if not os.path.isdir(a.out):
        os.makedirs(a.out)
    rng = np.random.RandomState(a.seed)
    dt = 0.05
    cfg = dict(min_fit=70.0, direct_min=0.08, direct_max=1.2, pre_ratio=0.5, post_ratio=1.5)
    rays = np.sort(rng.uniform(0.04, 0.08, a.n))
    rfs, kept = [], []
    for p in rays:
        z, r = synthetic_pair(a.h, a.kappa, a.vp, p, dt, rng, a.noise)
        rf, t, fit = itdecon(r, z, dt, a=a.gauss, shift=10.0)
        ok = rf_quality(rf, t, fit, cfg)[0]
        if ok:
            sel = (t >= -5.0) & (t <= 60.0)
            rfs.append(rf[sel])
            kept.append(p)
    if len(rfs) < 5:
        sys.exit('FAILED: only %d of %d synthetic receiver functions were accepted' % (len(rfs), a.n))
    R = np.vstack(rfs)
    p = np.array(kept)
    t0 = -5.0
    H = np.arange(20.0, 60.0 + 1e-9, 0.25)
    K = np.arange(1.60, 2.00 + 1e-9, 0.005)
    G = hk_stack_each(R, t0, dt, p, a.vp, H, K)
    S = G.mean(axis=0).astype(float)
    h_found, k_found, _ = hk_best(S, H, K)
    boot = []
    for _ in range(100):
        idx = rng.randint(0, len(R), len(R))
        boot.append(hk_best(G[idx].mean(axis=0), H, K)[:2])
    boot = np.array(boot)
    t_ps, h_ps, h_ps_err = ps_depth(R, t0, dt, p, a.vp, 1.75)
    print('%d of %d synthetic receiver functions accepted' % (len(R), a.n))
    print('crust put in:   H %.1f km, Vp/Vs %.3f' % (a.h, a.kappa))
    print(
        'H-k stack:      H %.1f +/- %.1f km, Vp/Vs %.3f +/- %.3f'
        % (h_found, boot[:, 0].std(), k_found, boot[:, 1].std())
    )
    print('Ps only:        H %.1f +/- %.1f km with Vp/Vs 1.75, Ps at %.2f s' % (h_ps, h_ps_err, t_ps))

    tt = t0 + np.arange(R.shape[1]) * dt
    fig, axs = plt.subplots(1, 2, figsize=(11, 5), gridspec_kw=dict(width_ratios=[1.5, 1]))
    ax = axs[0]
    scale = 0.8 / np.percentile(np.abs(R), 99)
    for k in range(len(R)):
        y = R[k] * scale + k
        ax.plot(tt, y, 'k', lw=0.3)
        ax.fill_between(tt, k, y, where=y > k, color='0.3', lw=0)
    times = np.array([phase_times(a.h, a.kappa, a.vp, pp) for pp in p])
    for j, (lab, col) in enumerate((('PdS', '0.0'), ('PpPs', '0.4'), ('PpSs+PsPs', '0.6'))):
        ax.plot(times[:, j], np.arange(len(R)), color=col, lw=1.0, ls='--', label=lab)
    ax.set_xlim(-2, 30)
    ax.set_ylim(-1, len(R) + 1)
    ax.set_xlabel('time after P (s)')
    ax.set_ylabel('synthetic receiver functions sorted by ray parameter')
    ax.legend(frameon=False, fontsize=8, loc='upper right', title='predicted times')
    ax.set_title('%d receiver functions recovered by the deconvolution' % len(R), fontsize=9)
    ax = axs[1]
    im = ax.pcolormesh(K, H, S, cmap='Greys', shading='auto')
    ax.contour(K, H, S, levels=[0.9 * S.max()], colors='k', linewidths=0.8)
    ax.scatter(boot[:, 1], boot[:, 0], s=4, c='0.5', alpha=0.5)
    ax.plot(k_found, h_found, 'k+', ms=12, mew=1.5, label='maximum of the stack')
    ax.plot(a.kappa, a.h, 'o', mfc='none', mec='#b2182b', ms=11, mew=1.3, label='crust put in')
    ax.set_xlabel('Vp/Vs')
    ax.set_ylabel('Moho depth H (km)')
    ax.legend(frameon=False, fontsize=8, loc='upper right')
    ax.set_title('H-k stack: %.1f km, Vp/Vs %.2f' % (h_found, k_found), fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.6, label='stack amplitude')
    fig.tight_layout()
    save_figure(fig, a.out, 'synthetic_test.png')
    plt.close(fig)
    good = abs(h_found - a.h) <= 1.0 and abs(k_found - a.kappa) <= 0.03
    print('PASSED' if good else 'FAILED: the H-k stack does not recover the crust within 1 km and 0.03')
    if not good:
        sys.exit(1)
