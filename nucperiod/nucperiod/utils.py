import collections
import itertools


ROMAN_DICT = collections.OrderedDict(sorted({1000: "M", 900: "CM", 500: "D", 400: "CD", 100: "C", 90: "XC", 50: "L", 40: "XL", 10: "X", 9: "IX",
              5: "V", 4: "IV", 1: "I"}.items(), reverse=True))


def roman_num(num):
    for r in ROMAN_DICT.keys():
        x, y = divmod(num, r)
        yield ROMAN_DICT[r] * x
        num -= (r * x)
        if num > 0:
            roman_num(num)
        else:
            break


def int2roman(num):
    return "".join([a for a in roman_num(num)])


def slicing_window(seq, n):
    it = iter(seq)
    result = ''.join(itertools.islice(it, n))

    if len(result) == n:
        yield result

    for elem in it:
        result = result[1:] + elem
        yield result
