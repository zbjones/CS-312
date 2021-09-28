import random


def prime_test(N, k):
    # This is main function, that is connected to the Test button. You don't need to touch it.
    return fermat(N, k), miller_rabin(N, k)


# The time complexity of the the mod_exp function is O(m*n^2) where 'm' is the # of bits of 'y' if m is roughly the
# same size as x and N, we get O(n^3). The complexity of the multiplication and mod step is O(n^2), but we run through
# our recursion log2(y) times, or the length in bits of the input y, or 'm', meaning we can say O(n^3).
# It is a fair assumption to make that 'm' will be of similar size to x and N.
# The space complexity is n recursive calls having n bits of input each, yielding O(n^2)
def mod_exp(x, y, N):
    if y == 0:  # O(1)
        return 1  # O(1)
    z = mod_exp(x, y // 2, N)  # we are calling this log2(y) times, or 'n' times if y is similar to x and N.
    if y % 2 == 0:  # (if it is even), O(1)
        product = z * z  # O(n^2)
        return product % N  # O(n^2)
    else:
        product = x * z * z  # O(n^2)
        return product % N  # O(n^2)


# Space complexity is constant as no values are stored.
# Time complexity is O(k*n^2) but bit shifts mean that instead of multiplying (complexity of O(n^2)) k times,
# we can bit shift in constant time since we are dividing by 2,  yielding an overall complexity of relating
# to the number of bits of k, O(log2(k)), which we have equate to O(n)
def fprobability(k):
    return 1 - (1 / pow(2, k))


# Space complexity is constant as no values are stored.
# Time complexity is O(k*n^2) but bit shifts mean that instead of multiplying (complexity of O(n^2)) k times,
# we can bit shift in constant time since we are dividing by a power of 2, yielding an overall complexity of relating
# to the number of bits of k, O(log2(k)), which we have equated to O(n)
def mprobability(k):
    return 1 - (1 / pow(4, k))


# The time complexity is O(k*n^3).  We are looping k- times and we are running mod_exp, which is O(n^3).  If k is small,
# we can assume that it is simply O(N^3) and drop the k.  But the full expression should be O(k*N^3)
# The space complexity is O(n^2). This is because we are running mod_exp k- times which has a space complexity of
# O(n^2).  As n gets large, the only thing that changes is the number of times that we run mod_exp. We can drop the k
# as well, as there is no need to concurrently store the data from each mod_exp run.
def fermat(N, k):
    for i in range(0, k):  # We loop k- times
        a = random.randint(2, N - 1)
        test_num = mod_exp(a, N - 1, N)  # O(n^3)
        if test_num != 1:
            return 'composite'
    return 'prime' # we assume it is prime if we don't get a negative result (something besides 1 in k-trials)


# The time complexity is O(k*n^4).  We loop k- times through the number of trials.  We also have a while- loop that
# runs at most log2(N) times, which we can simplify to O(n) since log2(y) = n as long as x, y, and N are similar size
# (reasonable assumption). We still call mod_exp, with complexity O(n^3). Therefore, we combine O(n^3) and O(n) and k-
# trials, yielding O(k*n^4).  If k is much smaller than n, we can drop k and assume O(n^4) time.
# Space complexity is O(n^2). This is from mod_exp's recursion. n^2 dominates the increasing size of storing
# random_int, exponent, and the result of mod_exp. Running mod_exp multiple times does not change the complexity as
# each time mod_exp function is called, it can use the same allocated space that grows by O(n^2). (e.g. looping k- times
# does not necessarily increase the space complexity from O(n^2) to O(k*n^2).
def miller_rabin(N, k):
    if N % 2 == 0:
        return 'composite' # we can catch even #'s really quickly.
    for i in range(0, k):  # O(k)
        random_int = random.randint(2, N - 1)
        exponent = N - 1
        while exponent % 2 == 0:  # this step repeats at most log2(N) times, or the number of bits in N. O(n) time.
            z = mod_exp(random_int, exponent, N)  # mod_exp is O(n^3)
            if z == N - 1:  # if mod N is equal to N - 1:
                break # we don't know anything if we get mod N  == -1.  Must get a new random int to check again.
            if z != 1:
                return 'composite'  # we know it's composite if we don't get 1 or -1.
            else:
                exponent = exponent / 2  # Division is O(n^2), but this is a power of 2 so we can do a constant time
                # bit shift that only increases linearly with the number of bits in the exponent. So we can say O(n)
                # here.
    return 'prime'  # we default to prime if we don't prove ourselves wrong before.


# Some useful numbers to check:
# 561, 1105, 29341, 10585, 75361
# 82589933
# 232250619601
# 34830684315505
