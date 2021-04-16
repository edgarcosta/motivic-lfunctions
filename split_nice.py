import sys

def split(fn_in, fn_out, col, n, sep=':'):
    num_lines = sum(1 for line in open(fn_in))
    desired_lines = (num_lines + n - 1)//n
    splits = []
    with open(fn_in) as iF:
        searching = False
        oldcol = None
        begin = 0
        for i, line in enumerate(iF):
            if i != 0 and i % desired_lines == 0:
                searching = True
                oldcol = line.split(sep)[col]
            if searching:
                newcol = line.split(sep)[col]
                if newcol != oldcol:
                    splits.append((begin, i-1))
                    begin = i
                    searching = False
        else:
            splits.append((begin, i))



    with open(fn_in) as iF:
        with open(fn_out) as oF:
            for j, (beg, end) in enumerate(splits):
                with open(fn_in + '.%d' % j, 'w') as iW:
                    with open(fn_out + '.%d' % j, 'w') as oW:
                        for j, (i, o) in enumerate(zip(iF, oF), beg):
                            if j <= end:
                                iW.write(i)
                                oW.write(o)
                            if j == end:
                                break




if __name__ == "__main__":
    try:
        fn_in = sys.argv[1]
        fn_out = sys.argv[2]
        col = int(sys.argv[3])
        n = int(sys.argv[4])
    except ValueError:
        print('Usage: %s input_filename output_filename col n')
        exit()
    split(fn_in, fn_out, col, n)


