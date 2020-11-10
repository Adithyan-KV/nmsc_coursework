def main():
    recurrence(1.24)
    recurrence(1.25)


def recurrence(a):
    print(f'a={a}')
    a_1 = a
    a_2 = a
    sequence = [a_1, a_2]
    errors = [0, 0]

    for i in range(2, 20):
        a_i = 10 * sequence[i - 1] - 9 * sequence[i - 2]
        sequence.append(a_i)
        error = abs(a_i - a) / a
        errors.append(error)
    format_and_print(sequence, errors)


def format_and_print(sequence, errors):
    for i in range(len(sequence)):
        # print('a_{} = {:.16f}; Relative error = {:.16f}'.format(
        #     i, sequence[i], errors[i]))
        print('a_{} & {: .16f} & {: .16f} \\\\'.format(
            i + 1, sequence[i], errors[i]))


if __name__ == "__main__":
    main()
