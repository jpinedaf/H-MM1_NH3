import astropy.units as u
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
plt.ion()

# prepare data
log_N_NH3 = fits.getdata('fits/H-MM1_NH3_Ncol_snr5.fits')
elog_N_NH3 = fits.getdata('fits/H-MM1_NH3_eNcol_snr5.fits')
N_H2 = fits.getdata('data/H-MM1_NH2_6arcsec_aligned.fits')
gd_data = np.isfinite(log_N_NH3 / N_H2)

N_NH3_14 = 10**(log_N_NH3[gd_data] - 14)
eN_NH3_14 = N_NH3_14 * np.log(10.) * elog_N_NH3[gd_data]
N_H2_22 = N_H2[gd_data] * 1e-22

#plt.scatter(N_H2_22, N_NH3_14, alpha=0.4)
plt.errorbar(N_H2_22, N_NH3_14, yerr=eN_NH3_14, 
    fmt='bo', alpha=0.4, zorder=3)

def broken_line(x, a, b, c, x0):
    out = a*x + b
    gd_brk = (x > x0)
    out[gd_brk] = c*x[gd_brk] + b + x0 * (a - c)
    return out

def log_prior(theta):
    a, b, c, x0 = theta
    if 0.1 < a < 20. and -10.0 < b < 10.0 and \
    -1.0 < c < 3.0 and 1.0 < x0 < 4.0:
        return 0.0
    return -np.inf

def log_likelihood(theta, x, y, yerr):
    a, b, c, x0 = theta
    model = a*x + b
    gd_brk = (x > x0)
    model[gd_brk] = c*x[gd_brk] + b + x0 * (a - c)
    sigma2 = yerr ** 2
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)

a0 = 1.975
b0 = 1.120
c0 = -0.110
x00 = 2.575 

xsample = np.arange(0, 6.3, 0.1)
plt.plot(xsample, broken_line(xsample, a0, b0, c0, x00), 
    color='r', zorder=10)

#
# Maximum Likelyhood fit
#
from scipy.optimize import minimize

nll = lambda *args: -log_likelihood(*args)
initial = np.array([a0, b0, c0, x00])
soln = minimize(nll, initial, args=(N_H2_22, N_NH3_14, eN_NH3_14))
a_ml, b_ml, c_ml, x0_ml = soln.x

print("Maximum likelihood estimates:")
print("a = {0:.3f}".format(a_ml))
print("b = {0:.3f}".format(b_ml))
print("c = {0:.3f}".format(c_ml))
print("x0 = {0:.3f}".format(x0_ml))

plt.plot(xsample, broken_line(xsample, a_ml, b_ml, c_ml, x0_ml), 
    color='g', zorder=12)

#
# EMCEE
#
import emcee

pos = soln.x + 1e-4 * np.random.randn(32, 4)
nwalkers, ndim = pos.shape

sampler = emcee.EnsembleSampler(
    nwalkers, ndim, log_probability, 
    args=(N_H2_22, N_NH3_14, eN_NH3_14)
)
sampler.run_mcmc(pos, 5000, progress=True)


fig, axes = plt.subplots(4, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
labels = ["a", "b", "c", "x0"]
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number");
# Autocorrelation
tau = sampler.get_autocorr_time()
print(tau) # ~50

# remove burn-in and thin ~1/2 of autocorrelation 
flat_samples = sampler.get_chain(discard=100, thin=25, flat=True)
print(flat_samples.shape)

import corner

fig = corner.corner(
    flat_samples, labels=labels)


for i in range(ndim):
    mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(mcmc)
    txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
    txt = txt.format(mcmc[1], q[0], q[1], labels[i])
    print(txt)
#
#a = 1.975_{-0.005}^{0.005}
#b = 1.120_{-0.009}^{0.008}
#c = -0.110_{-0.005}^{0.005}
#x0 = 2.575_{-0.004}^{0.003}
#
