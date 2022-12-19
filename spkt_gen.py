'''
this is the code for generating the spike times for 10seconds of A and C fibers
created by K. Sekiguchi (18th July 20)
'''

import numpy as np
from cfg_mechanical import cfg

def rate_SAI():
    t = [ x for x in np.arange(0, 10001, 1) ]

    rate = []

    for i, key in enumerate(t):
        if i < len(t):
            freq_SAI = (-1.45433609113392e-10 * key ** 3 + 1.340603396708e-6 * key ** 2 - 0.00378224210238498 * key + 4.52737468545426) * 8
            rate.append(freq_SAI)

    return rate, t

def rate_SAII():
    t = [ x for x in np.arange(0, 10001, 1) ]

    rate = []

    for i, key in enumerate(t):
        if i < len(t):
            freq_SAII = (-1.45433609113392e-10 * key ** 3 + 1.340603396708e-6 * key ** 2 - 0.00378224210238498 * key + 4.52737468545426) * 8
            rate.append(freq_SAII)

    return rate, t

def rate_Ad():
    t = [ x for x in np.arange(0, 10001, 1) ]

    rate = []

    for i, key in enumerate(t):
        if i < len(t):
            freq_Ad = ( -7.265297e-15 * key ** 4 + 1.573831e-10 * key ** 3 - 8.201529e-7 * key ** 2 - 0.002539930 * key + 27.18412 ) * cfg.stim_ratios
            rate.append(freq_Ad)

    return rate, t

def rate_C():
    t = [ x for x in np.arange(0, 10001, 1) ]

    rate  = []

    for i, key in enumerate(t):
        if i < len(t):
            freq_C = (  9.630390e-15 * key ** 4 - 2.844577e-10 * key ** 3 + 3.012887e-6 * key ** 2 - 0.01400381 * key + 31.86546 ) * cfg.stim_ratios
            rate.append(freq_C)

    return rate, t

def poisson_generator(rate, t_start=0.0, t_stop=1000.0, seed=None):
    """
    Returns a SpikeTrain whose spikes are a realization of a Poisson process
    with the given rate (Hz) and stopping time t_stop (milliseconds).
    Note: t_start is always 0.0, thus all realizations are as if 
    they spiked at t=0.0, though this spike is not included in the SpikeList.
    Inputs:
    -------
        rate    - the rate of the discharge (in Hz)
        t_start - the beginning of the SpikeTrain (in ms)
        t_stop  - the end of the SpikeTrain (in ms)
        array   - if True, a np array of sorted spikes is returned,
                    rather than a SpikeTrain object.
    Examples:
    --------
        >> gen.poisson_generator(50, 0, 1000)
        >> gen.poisson_generator(20, 5000, 10000, array=True)
    See also:
    --------
        inh_poisson_generator, inh_gamma_generator, inh_adaptingmarkov_generator
    """
    rng = np.random.RandomState(seed)
    #number = int((t_stop-t_start)/1000.0*2.0*rate)
    # less wasteful than double length method above
    n = (t_stop-t_start)/1000.0*rate
    number = np.ceil(n+3*np.sqrt(n))
    if number<100:
        number = min(5+np.ceil(2*n),100)
    if number > 0:
        isi = rng.exponential(1.0/rate, int(number))*1000.0
        if number > 1:
            spikes = np.add.accumulate(isi)
        else:
            spikes = isi
    else:
        spikes = np.array([])
    spikes+=t_start
    i = np.searchsorted(spikes, t_stop)
    extra_spikes = []
    if i==len(spikes):
        # ISI buf overrun
        
        t_last = spikes[-1] + rng.exponential(1.0/rate, 1)[0]*1000.0
        while (t_last<t_stop):
            extra_spikes.append(t_last)
            t_last += rng.exponential(1.0/rate, 1)[0]*1000.0
        
        spikes = np.concatenate((spikes,extra_spikes))
    else:
        spikes = np.resize(spikes,(i,))
        
    return spikes

def inh_poisson_generator(rate, t, t_stop, seed=None):
    """
    Returns a SpikeTrain whose spikes are a realization of an inhomogeneous 
    poisson process (dynamic rate). The implementation uses the thinning 
    method, as presented in the references.
    Inputs:
    -------
        rate   - an array of the rates (Hz) where rate[i] is active on interval 
                    [t[i],t[i+1]]
        t      - an array specifying the time bins (in milliseconds) at which to 
                    specify the rate
        t_stop - length of time to simulate process (in ms)
        array  - if True, a np array of sorted spikes is returned,
                    rather than a SpikeList object.
    Note:
    -----
        t_start=t[0]
    References:
    -----------
    Eilif Muller, Lars Buesing, Johannes Schemmel, and Karlheinz Meier 
    Spike-Frequency Adapting Neural Ensembles: Beyond Mean Adaptation and Renewal Theories
    Neural Comput. 2007 19: 2958-3010.
    Devroye, L. (1986). Non-uniform random variate generation. New York: Springer-Verlag.
    Examples:
    --------
        >> time = arange(0,1000)
        >> stgen.inh_poisson_generator(time,sin(time), 1000)
    See also:
    --------
        poisson_generator, inh_gamma_generator, inh_adaptingmarkov_generator
    """
    rng = np.random.RandomState(seed)
    if np.shape(t)!=np.shape(rate):
        raise ValueError('shape mismatch: t,rate must be of the same shape')
    # get max rate and generate poisson process to be thinned
    rmax = np.max(rate)
    ps = poisson_generator(rate=rmax, t_start=t[0], t_stop=t_stop, seed=None)
    # return empty if no spikes
    if len(ps) == 0:
        np.array([])
        
    # gen uniform rand on 0,1 for each spike
    rn = np.array(rng.uniform(0, 1, len(ps)))
    # instantaneous rate for each spike
    idx = np.searchsorted(t, ps) - 1
    #spike_rate = rate[idx]
    spike_rate = np.array([rate[i] for i in idx])
    # thin and return spikes
    spike_train = ps[rn<spike_rate/rmax]
    return list(spike_train)

