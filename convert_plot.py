
def old_plot(plot_values, order_of_vanishing, plot_delta=0.25, dual_values=None, dual_delta=None):
    """
    The current code in the LMFDB for producing L-function plots

    INPUT:

    - ``plot_values`` -- equally spaced plot values for x >= 0, as taken from the ``plot_values`` column, as a list of floats
    - ``order_of_vanishing`` -- a non-negative integer giving the order of vanishing at the central point
    - ``plot_delta`` -- the delta used in the ``plot_values`` list
    - ``dual_values`` -- ``plot_values`` for the dual L-function (if left as None, indicates self dual)
    - ``dual_delta`` -- the delta used in ``dual_values`` (defaults to the same as ``plot_delta``)

    OUTPUT:

    - the interpolation spline (used in creating the new plot)
    - a list of sampled values from the interpolation, every 0.05
    """
    # First we recreate the current LMFDB code for producing Z-plots from this data
    sign = (-1)**order_of_vanishing
    selfdual = dual_values is None
    plot_delta = RR(plot_delta) # TODO: precision
    if dual_delta is None:
        dual_delta = plot_delta
    # We need to increase precision, since some plot values were in RealField(20), which was not enough for the cubic interpolation we were doing.
    pos_plot = [(j * plot_delta, RR(val)) for j, val in enumerate(plot_values)]
    if selfdual:
        neg_plot = [(-x, sign*y) for (x,y) in pos_plot[1:]]
    else:
        neg_plot = [(-j * dual_delta, val) for j, val in enumerate(dual_values[1:],1)]
    neg_plot.reverse()
    plotpoints = neg_plot + pos_plot

    xmax = min(30, pos_plot[-1][0])
    xmin = max(-30, neg_plot[0][0])
    interpolation = spline(plotpoints)
    F_interp = [(i, interpolation(i)) for i in srange(xmin, xmax, 0.05)] + [(xmax, interpolation(xmax))]

    return interpolation, F_interp

def spline_to_newdata(interpolation, F_interp, pos_zeros, order_of_vanishing, mu_imag, dual_zeros):
    """
    Reconstruct the data to be stored in the new schema from the spline produced by the old data

    INPUT:

    - ``interpolation`` -- the spline object produced by old_plot.
    - ``F_interp`` -- the sampled interpolation (at intervals of 0.05, as in the current plotting code)
    - ``pos_zeros`` -- positive zeros of the L-function, as a list of floats
    - ``order_of_vanishing`` -- a non-negative integer giving the order of vanishing at the central point
    - ``mu_imag`` -- imaginary parts of the mu part of the Selberg data
    - ``dual_zeros`` -- positive zeroes for the dual L-function (if left as None, indicates self dual)

    OUTPUT:

    - ``extrema`` -- a list of (x,y) values giving local mins and maxes
    - ``zeros`` -- the list of (all) zeros
    - ``derivs`` -- the derivative of the L-function, evaluated at each zero
    - ``extra_control`` -- extra control points, as a list of triples (x,y,d)
    """
    xmin = F_interp[0][0]
    xmax = F_interp[-1][0]
    if dual_zeros is None:
        neg_zeros = [-x for x in reversed(pos_zeros)]
    else:
        neg_zeros = [-x for x in reversed(dual_zeros)]
    zero_zeros = [] if order_of_vanishing == 0 else [float(0)]
    zeros = neg_zeros + zero_zeros + pos_zeros
    # It's possible that we have zeros past where we have interpolation; we throw those away
    zeros = [x for x in zeros if xmin <= x <= xmax]
    derivs = [interpolation.derivative(x) for x in zeros]
    # Mins and maxes are a bit annoying
    extrema = []
    for (x0, y0), (x1, y1), (x2, y2) in zip(F_interp[:-2], F_interp[1:-1], F_interp[2:]):
        if y0 <= y1 and y2 <= y1 or y0 >= y1 and y2 >= y1:
            # Local extremum present
            # We just do a black-box approach, rather than using the cubic nature of the spline
            # Start with quadratic approximation for an initial guess, then Newton's method
            denom = (x0 - x1) * (x0 - x2) * (x1 - x2)
            A = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / denom
            B = (x2**2 * (y0 - y1) + x1**2 * (y2 - y0) + x0**2 * (y1 - y2)) / denom
            xguess = -B / (2*A)
            assert x0 <= xguess <= x2
            yderiv = interpolation.derivative(xguess)
            xguesses = [xguess]
            while abs(yderiv) > 0.001:
                yderiv = interpolation.derivative(xguess)
                yderiv2 = interpolation.derivative(xguess, order=2)
                xguess = xguess - yderiv / yderiv2
                xguesses.append(xguess)
                assert len(xguesses) < 10 and x0 <= xguess <= x2
                if xguess > x2:
                    xguess = x2
                elif xguess < x0:
                    xguess = x0
            yguess = interpolation(xguess)
            extrema.append((xguess, interpolation(xguess)))

    control = sorted([x for (x,y) in extrema] + zeros)
    # Add extra control points near 0, and near trivial zeros
    first_pos = min(pos_zeros)
    first_neg = max(neg_zeros)
    near_zero = sorted([first_neg, 0, first_pos] + [x for (x,y) in extrema if first_neg < x < first_pos and abs(x) > 0.01])
    near_zero = [(x0 + x1) / 2 for (x0, x1) in zip(near_zero[:-1], near_zero[1:])]

    # Add control points near trivial zeros
    triv = [RR(-y) for y in mu_imag if xmin < -y < first_neg or first_pos < -y < xmax]
    if triv:
        triv = sum(([x - 2, x - 1, x - 0.5, x, x + 0.5, x + 1, x + 2] for x in triv), [])
        # Remove extra control points that are close to zeros or extrema
        def closest_control(x):
            i = bisect(control, x)
            if i == 0:
                return control[0]
            elif i == len(control):
                return control[-1]
            if x - control[i-1] < control[i] - x:
                return control[i-1]
            return control[i]
        def far_from_control(x):
            return xmin < x < xmax and abs(x - closest_control(x)) > 0.05
        triv = [x for x in triv if far_from_control(x)]

    # Add left and right endpoints
    boundary = [xmin, xmax]
    xs = sorted(boundary + near_zero + triv)
    extra_control = [(x, interpolation(x), interpolation.derivative(x)) for x in xs]
    return extrema, zeros, derivs, extra_control

def convert_plot(plot_values, pos_zeros, order_of_vanishing, mu_imag, plot_delta, dual_values=None, dual_zeros=None, dual_delta=None):
    """
    INPUT:

    - ``plot_values`` -- equally spaced plot values for x >= 0, as taken from the ``plot_values`` column, as a list of floats
    - ``pos_zeros`` -- positive zeros of the L-function, as a list of floats
    - ``order_of_vanishing`` -- a non-negative integer giving the order of vanishing at the central point
    - ``mu_imag`` -- imaginary parts of the mu part of the Selberg data
    - ``plot_delta`` -- the delta used in the ``plot_values`` list
    - ``dual_values`` -- ``plot_values`` for the dual L-function (if left as None, indicates self dual)
    - ``dual_zeros`` -- positive zeroes for the dual L-function (if left as None, indicates self dual)
    - ``dual_delta`` -- the delta used in ``dual_values`` (defaults to the same as ``plot_delta``)

    OUTPUT:

    - ``plot_x`` -- the x-coordinates of extrema
    - ``plot_y`` -- the y-coordinates of extrema
    - ``plot_deriv`` -- the derivatives at the positive zeros
    - ``plot_extras`` -- extra control points, as a list of triples
    """
    interpolation, F_interp = old_plot(plot_values, order_of_vanishing, plot_delta=plot_delta)
    extrema, zeros, derivs, extra_control = spline_to_newdata(interpolation, F_interp, pos_zeros, order_of_vanishing, mu_imag, dual_zeros)
    plot_x = [x for (x,y) in extrema if x >= 0]
    plot_y = [y for (x,y) in extrema if x >= 0]
    plot_deriv = [d for (x,d) in zip(zeros, derivs) if x >= 0]
    plot_extras = [(x,y,d) for (x,y,d) in extra_control if x >= 0]
    return plot_x, plot_y, plot_deriv, plot_extras
