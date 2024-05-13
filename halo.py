

def progress_bar(current, total, time, barLength=10):
    percent = float(current) * 100 / total
    if current > 0:
        eta = '%.2f' % (time * (total - current) / float(current),)
    else:
        eta = '+oo'
    done = '▰' * int(percent / 100 * barLength)
    notdone = '▱' * (barLength - len(done))

    return ' Progress: [%s%s] %.2f%% ETA: %s seconds' % (done, notdone, percent, eta)


def chunkify(l, size=100):
    return [islice(l, elt * size, (elt + 1) * size)
            for elt in range(len(l) // size + 1)]
