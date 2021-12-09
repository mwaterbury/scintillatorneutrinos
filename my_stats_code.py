# Import stats functions
from scipy.stats import poisson
from scipy.stats import norm
from scipy.stats import chi2

# Import special functions
from scipy import special
from scipy.special import entr, logsumexp, betaln, gammaln as gamln, zeta

# Import integration function
from scipy.integrate import quad, dblquad, simps

# Import numpy
import numpy as np

# Import matplotlib
import matplotlib
from matplotlib import pyplot as plt
import matplotlib.image as img

def plot_poisson_pvals(Ns, Nb, alpha = 10**-8):
    lo = int(np.floor(poisson.ppf(alpha, Nb)))
    hi = int(np.ceil(poisson.ppf(1-alpha, Ns + Nb)))
    No = np.arange(lo-1,hi+1)
    like0 = poisson.pmf(No, Nb)
    like1 = poisson.pmf(No, Nb + Ns)
    pvals = poisson.sf(No, Nb)
    Zs = np.arange(1,6)
    pval_sigmas = [1-norm.cdf(Z) for Z in Zs]

    zeros = np.zeros(np.shape(No))

    matplotlib.rcParams.update({'font.size': 14})
    fig = plt.figure(figsize=(8*2,6))

    ax = plt.subplot(1,2,1)
    ax.scatter(No, like0, color='blue', s=8, label='$f(N_o, N_b)$')
    ax.scatter(No, like1, color='orange', s=8, label='$f(N_o, N_b + N_s)$')
    for (Z,pval_sigma) in zip(Zs,pval_sigmas):
        if pval_sigma >= alpha:
            mask = pvals <= pval_sigma
            ax.fill_between(No[mask], zeros[mask], like0[mask], color = 'blue', alpha=0.05)
            ax.text(No[mask][0]+0.5,alpha*10**0.125,f'$Z=${Z}', fontsize='x-small')

    
    ax.set_xlabel('$N_o$')
    ax.set_xlim([lo,hi])
    ax.set_ylim([10**-8,1])
    ax.set_yscale('log')
    ax.legend(loc="upper right")

    ax = plt.subplot(1,2,2)
    ax.scatter(No, pvals, s=8)
    for (Z,pval_sigma) in zip(Zs,pval_sigmas):
        ax.plot([lo-1, hi+1], [pval_sigma, pval_sigma], 'k--')
        ax.text(hi-0.1*(hi-lo), pval_sigma*10**0.125, f'$Z=${Z}')
    ax.set_xlabel('$N_o$')
    ax.set_xlim([lo,hi])
    ax.set_ylim([alpha*10**-0.5,1.5])
    ax.set_yscale('log')
    ax.set_yticks([10**-n for n in [0, 2, 4, 6, 8]]);

def q_stat(N, mu, Ns, Nb):
    loglike0 = np.array(poisson.logpmf(N, mu*Ns+Nb))
    loglike1 = np.array(poisson.logpmf(N, max(N,Nb)))
    maxloglike0 = np.max(loglike0)
    maxloglike1 = np.max(loglike1)
    return -2*(maxloglike0 - maxloglike1)

def plot_poisson_qstat(Ns, Nb, mu=0, alpha=10**-8):
    lo = int(np.floor(poisson.ppf(alpha, Nb)))
    hi = int(np.ceil(poisson.ppf(1-alpha, Ns + Nb)))
    No = np.arange(lo-1,2*hi+1)

    qvals = np.array([q_stat(n, mu, Ns, Nb) for n in No])
    pvals = np.array(poisson.sf(No, Nb))
    chi2s = 0.5*chi2.sf(qvals, df=1)

    dqdx = np.diff(qvals)
    mask = dqdx > 0
    fq = poisson.pmf(No[1:][mask], Nb) / dqdx[mask]

    matplotlib.rcParams.update({'font.size': 14})
    fig = plt.figure(figsize=(8*2,6))

    ax = plt.subplot(1,2,1)
    ax.plot(No, 0.5*chi2.pdf(qvals, df=1), '--', label='$\chi^2(q(N_o))$')
    ax.scatter(No[1:][mask], fq, s=16, label='$\\tilde f(q(N_o))$')
    ax.plot(No, poisson.pmf(No, Nb), '--', label='$f(N)$')
    ax.scatter(No[1:][mask], 0.5*chi2.pdf(qvals[1:][mask], df=1)/fq, s=16, label='$\\tilde f(q(N_o)) / \chi^2(q(N_o))$')

    ax.legend(loc='lower left')
    ax.set_yscale('log')
    ax.set_xlim([Nb,hi])
    ax.set_ylim([alpha * 10**-0.5,10**0.5])
    ax.set_xlabel('$N_o$')

    ax = plt.subplot(1,2,2)
    ax.plot(No, chi2s, '--', label = '$\\tilde{p}(N_o)$')
    ax.scatter(No, pvals, s=16, label = '$\hat{p}(N_o)$')
    ax.scatter(No, pvals/chi2s, s=16, label='$\hat{p}(N_o) / \\tilde p(N_o)$')

    ax.legend(loc='lower left')
    ax.set_yscale('log')
    ax.set_xlim([Nb,hi])
    ax.set_ylim([alpha * 10**-0.5,10**0.5])
    ax.set_xlabel('$N_o$')
    ax.set_yticks([10**-n for n in [0, 2, 4, 6, 8]]);

def marg_like(No, Ns, Nb, db):
    bs = np.linspace(max(Nb-10*db,0), Nb + 10*db)
    integrand = [poisson.pmf(No, Ns + b) * norm.pdf(b, Nb, db)
                 for b in bs]
    return simps(integrand, bs, axis=0)

def marg_loglike(No, Ns, Nb, db):
    return np.log(marg_like(No, Ns, Nb, db))

def marg_q_stat(N, mu, Ns, Nb, db):
    mus = np.linspace(0,Nb+10*db,100) / Ns
    loglike0 = np.array(marg_loglike(N, mu*Ns, Nb, db))
    loglike1 = np.array([marg_loglike(N, mu*Ns, Nb, db)
                for mu in mus])
    maxloglike0 = np.max(loglike0)
    maxloglike1 = np.max(loglike1)
    return -2*(maxloglike0 - maxloglike1)

def plot_marginal_pvals(Ns, Nb, db, alpha=10**-8):
    lo = int(np.floor(poisson.ppf(alpha, Nb)))
    hi = int(np.ceil(poisson.ppf(1-alpha, Ns + Nb)))
    No = np.arange(lo-1,hi+1)
    like0 = marg_like(No, 0, Nb, db)
    like1 = marg_like(No, Ns, Nb, db)
    pvals = np.array([1-np.sum(like0[:n]) for n in range(len(No))])
    Zs = np.arange(1,6)
    pval_sigmas = [1-norm.cdf(Z) for Z in Zs]

    zeros = np.zeros(np.shape(No))

    matplotlib.rcParams.update({'font.size': 14})
    fig = plt.figure(figsize=(8*2,6))

    ax = plt.subplot(1,2,1)
    ax.scatter(No, like0, color='blue', s=8, label='$f(N_o, N_b)$')
    ax.scatter(No, like1, color='orange', s=8, label='$f(N_o, N_b + N_s)$')
    for (Z,pval_sigma) in zip(Zs,pval_sigmas):
        if pval_sigma >= alpha:
            mask = pvals <= pval_sigma
            if np.sum(mask) > 0:
                ax.fill_between(No[mask], zeros[mask], like0[mask], color = 'blue', alpha=0.05)
                ax.text(No[mask][0]+0.5,alpha*10**0.125,f'$Z=${Z}', fontsize='x-small')

    
    ax.set_xlabel('$N_o$')
    ax.set_xlim([lo,hi])
    ax.set_ylim([alpha*10**-0.5,1.5])
    ax.set_yscale('log')
    ax.legend(loc="upper right")

    ax = plt.subplot(1,2,2)
    ax.scatter(No, pvals, s=8)
    for (Z,pval_sigma) in zip(Zs,pval_sigmas):
        ax.plot([lo-1, hi+1], [pval_sigma, pval_sigma], 'k--')
        ax.text(hi-0.1*(hi-lo), pval_sigma*10**0.125, f'$Z=${Z}')
    ax.set_xlabel('$N_o$')
    ax.set_xlim([lo,hi])
    ax.set_ylim([alpha*10**-0.5,1.5])
    ax.set_yscale('log')
    ax.set_yticks([10**-n for n in [0, 2, 4, 6, 8]]);

def plot_marg_qstat(Ns, Nb, db, mu=0, alpha=10**-8):
    lo = int(np.floor(poisson.ppf(alpha, Nb)))
    hi = int(np.ceil(poisson.ppf(1-alpha, Ns + Nb)))
    No = np.arange(lo-1,hi+1)

    qvals = np.array([marg_q_stat(n, mu, Ns, Nb, db) for n in No])
    like0 = marg_like(No, 0, Nb, db)
    like1 = marg_like(No, Ns, Nb, db)
    pvals = np.array([1-np.sum(like0[:n]) for n in range(len(No))])
    chi2s = 0.5*chi2.sf(qvals, df=1)

    dqdx = np.diff(qvals)
    mask = dqdx > 0
    fq = poisson.pmf(No[1:][mask], Nb) / dqdx[mask]

    matplotlib.rcParams.update({'font.size': 14})
    fig = plt.figure(figsize=(8*2,6))

    ax = plt.subplot(1,2,1)
    ax.plot(No, 0.5*chi2.pdf(qvals, df=1), '--', label='$\chi^2(q(N_o))$')
    ax.scatter(No[1:][mask], fq, s=16, label='$\\tilde f(q(N_o))$')
    ax.plot(No, like0, '--', label='$f(N)$')
    ax.scatter(No[1:][mask], 0.5*chi2.pdf(qvals[1:][mask], df=1)/fq, s=16, label='$\\tilde f(q(N_o)) / \chi^2(q(N_o))$')

    ax.legend(loc='lower left')
    ax.set_yscale('log')
    ax.set_xlim([Nb,hi])
    ax.set_ylim([alpha * 10**-0.5,10**0.5])
    ax.set_xlabel('$N_o$')

    ax = plt.subplot(1,2,2)
    ax.plot(No, chi2s, '--', label = '$\\tilde{p}(N_o)$')
    ax.scatter(No, pvals, s=16, label = '$\hat{p}(N_o)$')
    ax.scatter(No, pvals/chi2s, s=16, label='$\hat{p}(N_o) / \\tilde p(N_o)$')

    ax.legend(loc='lower left')
    ax.set_yscale('log')
    ax.set_xlim([Nb,hi])
    ax.set_ylim([alpha * 10**-0.5,10**0.5])
    ax.set_xlabel('$N_o$')
    ax.set_yticks([10**-n for n in [0, 2, 4, 6, 8]]);

def plot_compare_marg(Ns, Nb, db, mu=0, alpha=10**-8):
    yticks = 10**np.arange(np.log10(alpha),0.1,2)
    lo = int(np.floor(poisson.ppf(alpha*10**-6, Nb)))
    hi = int(np.ceil(poisson.ppf(1-alpha*10**-6, Ns + Nb)))
    No = np.arange(lo-1,hi+1)
    like0 = poisson.pmf(No, Nb)
    like1 = poisson.pmf(No, Nb + Ns)
    pvals = poisson.sf(No, Nb)

    marg_like0 = marg_like(No, 0, Nb, db)
    marg_like1 = marg_like(No, Ns, Nb, db)
    marg_pvals = np.array([1-np.sum(marg_like0[:n]) for n in range(len(No))])

    Zs = np.arange(1,8.1)
    pval_sigmas = np.array([1-norm.cdf(Z) for Z in Zs])
    Zs = Zs[pval_sigmas >= alpha]
    pval_sigmas = pval_sigmas[pval_sigmas >= alpha]

    zeros = np.zeros(np.shape(No))

    matplotlib.rcParams.update({'font.size': 14})
    fig = plt.figure(figsize=(8*2,6))

    ax = plt.subplot(1,2,1)
    ax.plot([Ns + Nb, Ns + Nb], [alpha*10**-0.5,10**0.5], 'r--')
    ax.text(Ns+Nb+0.01*(hi-lo),10**0.125,'$N_s + N_b$', color='red')
    ax.scatter(No, like0, color='blue', marker='<', s=16, label='$f(N_o, N_b)$')
    ax.scatter(No, marg_like0, color='orange', marker='<', s=16, label='$f_\mathrm{marg}(N_o, N_b, \sigma_b)$')
    
    ax.scatter(No, like1, color='blue', marker='>', s=16, label='$f(N_o, N_b + N_s)$')
    ax.scatter(No, marg_like1, color='orange', marker='>', s=16, label='$f_\mathrm{marg}(N_o, N_s, N_b, \sigma_b)$')
    for (Z,pval_sigma) in zip(Zs,pval_sigmas):
        if pval_sigma >= alpha:
            mask = pvals <= pval_sigma
            ax.fill_between(No[mask], zeros[mask], like0[mask], color = 'blue', alpha=0.05)
            ax.text(No[mask][0],alpha*10**-0.375,f'{int(Z)}$\sigma$', fontsize='x-small')
            
            mask = marg_pvals <= pval_sigma
            ax.fill_between(No[mask], zeros[mask], marg_like0[mask], color = 'orange', alpha=0.05)
    
    

    ax.legend(loc="upper right")
    
    ax.set_xlabel('$N_o$')
    ax.set_xlim([lo,hi])
    ax.set_ylim([alpha*10**-0.5,10**0.5])
    ax.set_yscale('log')
    ax.set_yticks(yticks);

    ax = plt.subplot(1,2,2)
    ax.plot([Ns + Nb, Ns + Nb], [alpha*10**-0.5,10**0.5], 'r--')
    ax.text(Ns+Nb+0.01*(hi-lo),10**0.125,'$N_s + N_b$', color='red')
    ax.scatter(No, pvals, s=8, label='$\hat{p}(N_o)$')
    ax.scatter(No, marg_pvals, s=8, label='$\hat{p}_\mathrm{marg}(N_o)$')
    for (Z,pval_sigma) in zip(Zs,pval_sigmas):
        ax.plot([lo-1, hi+1], [pval_sigma, pval_sigma], 'k--')
        ax.text(hi-0.06*(hi-lo), pval_sigma*10**0.125, f'{int(Z)}$\sigma$')
    ax.legend(loc="lower left")

    ax.set_xlabel('$N_o$')
    ax.set_xlim([lo,hi])
    ax.set_ylim([alpha*10**-0.5,10**0.5])
    ax.set_yscale('log')
    ax.set_yticks(yticks);