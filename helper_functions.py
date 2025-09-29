import numpy as np

def upsample(x, sps):
    """Insert sps-1 zeros between samples."""
    x = np.asarray(x).reshape(-1)
    y = np.zeros(len(x)*sps, dtype=x.dtype)
    y[::sps] = x
    return y

def rrc_pulse(beta, span, sps):
    """
    Root-raised-cosine (RRC) filter impulse response.
    beta: roll-off, span: in symbols, sps: samples per symbol
    """
    N = span * sps
    # Make it odd length (N+1 taps)
    t = np.arange(-N/2, N/2 + 1) / sps
    h = np.zeros_like(t, dtype=float)
    for i, ti in enumerate(t):
        if np.isclose(ti, 0.0):
            h[i] = 1.0 - beta + (4*beta/np.pi)
        elif np.isclose(abs(ti), 1/(4*beta)):
            # special-case singularities
            h[i] = (beta/np.sqrt(2)) * (
                ((1+2/np.pi)*np.sin(np.pi/(4*beta))) + ((1-2/np.pi)*np.cos(np.pi/(4*beta)))
            )
        else:
            h[i] = (np.sin(np.pi*ti*(1-beta)) + 4*beta*ti*np.cos(np.pi*ti*(1+beta))) / \
                   (np.pi*ti*(1-(4*beta*ti)**2))
    # Normalize for unit energy
    h = h / np.sqrt(np.sum(h**2))
    return h

def qammod(M, symbols):
    """
    Map integers in [0, M-1] to square QAM constellation (Gray-ish), unit average power.
    """
    m = int(np.sqrt(M))
    assert m*m == M, "M must be a perfect square for this qammod."
    # Gray mapping for I and Q roughly
    # Generate PAM levels
    levels = np.arange(-(m-1), m, 2)
    I = levels[(symbols % m).astype(int)]
    Q = levels[(symbols // m).astype(int)]
    const = I + 1j*Q
    # Normalize to unit average power
    const = const / np.sqrt(np.mean(np.abs(const)**2))
    return const

def ecdf(x):
    """Empirical CDF: returns values f, t such that f[i] = P(X <= t[i])"""
    x = np.asarray(x, dtype=float)
    t = np.sort(x)
    f = np.arange(1, len(t)+1, dtype=float)/len(t)
    return f, t

# Fast moving average: O(N)
def movmean(x, L):
    # if L == 1 or L >= x.size:
    #     return np.full_like(x, x.mean(), dtype=complex)
    c = np.cumsum(np.insert(x, 0, 0.0))
    core = (c[L:] - c[:-L]) / L
    # pad to 'same' length (edge replication)
    pad_left = np.full(L//2, core[0])
    pad_right = np.full(x.size - core.size - pad_left.size, core[-1])
    return np.concatenate([pad_left, core, pad_right])

def movvar(x, L):
    """
    Fast moving variance with window length L (centered, same-length output).
    Works for real or complex x. O(N) using cumulative sums.
    Variance for complex signals uses E[|x|^2] - |E[x]|^2.
    Edge handling: replicate edges (centered padding).
    """
    x = np.asarray(x)
    n = x.size
    L = int(L)

    if L <= 1:
        # window of 1 -> zero variance everywhere
        return np.zeros(n, dtype=float)
    if L >= n:
        # fall back to global variance
        v = np.var(x.astype(np.complex128)) if np.iscomplexobj(x) else np.var(x.astype(float))
        return np.full(n, float(v), dtype=float)

    # cumulative sums: x and |x|^2
    # use complex accumulator for x (handles real/complex), float for |x|^2
    csum  = np.cumsum(np.insert(x.astype(np.complex128), 0, 0.0 + 0.0j))
    csum2 = np.cumsum(np.insert(np.abs(x)**2,         0, 0.0        ))

    # windowed sums over length L
    win_sum  = csum[L:]  - csum[:-L]      # complex (sum of x)
    win_sum2 = csum2[L:] - csum2[:-L]     # real    (sum of |x|^2)

    mean = win_sum / L                    # complex mean
    var_core = win_sum2 / L - np.abs(mean)**2   # real variance
    var_core = np.maximum(var_core, 0.0)        # clip tiny negatives

    # Centered padding to return same length
    # core length is n-L+1; total padding needed = L-1
    left = (L - 1) // 2
    right = (L - 1) - left
    pad_left = np.full(left,  var_core[0],  dtype=float)
    pad_right= np.full(right, var_core[-1], dtype=float)

    return np.concatenate([pad_left, var_core, pad_right])