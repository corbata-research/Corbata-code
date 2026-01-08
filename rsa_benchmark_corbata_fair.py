import random
import math
import time
import json
from statistics import mean, median, pstdev

try:
    import gmpy2
    HAS_GMPY2 = True
except Exception:
    HAS_GMPY2 = False

try:
    import sympy as sp
    HAS_SYMPY = True
except Exception:
    HAS_SYMPY = False


# ================================================================
# Miller Rabin primality test (probabilistic)
# ================================================================
def is_probable_prime_mr(n, rounds=24, rng=None):
    """Return True if n is probably prime, False if composite (Miller Rabin)."""
    if n < 2:
        return False

    small_primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    if n in small_primes:
        return True
    for p in small_primes:
        if n % p == 0:
            return False

    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    if rng is None:
        rng = random.Random()

    for _ in range(rounds):
        a = rng.randrange(2, n - 1)
        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        witness_found = True
        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                witness_found = False
                break

        if witness_found:
            return False

    return True


def primality_status(n, mr_rounds=24, rng=None):
    """
    Return (is_prime_bool, status_string).
    Uses:
      1 Miller Rabin filter
      2 gmpy2 is_prime if available
      3 sympy isprime if available
    """
    if n < 2:
        return False, "Composite (n < 2)"

    if not is_probable_prime_mr(n, rounds=mr_rounds, rng=rng):
        return False, f"Composite (fails Miller Rabin, rounds={mr_rounds})"

    if HAS_GMPY2:
        res = int(gmpy2.is_prime(n))
        if res == 0:
            return False, "Composite (gmpy2 witness)"
        if res == 2:
            return True, "Prime (gmpy2 proved)"
        return True, "Probably prime (gmpy2 strong test)"

    if HAS_SYMPY:
        if sp.isprime(n):
            return True, "Probably prime (sympy isprime)"
        return False, "Composite (sympy witness)"

    return True, f"Probably prime (Miller Rabin only, rounds={mr_rounds})"


# ================================================================
# Corbata row generator (standalone)
# Based on your generate_corbata_row_local idea
# ================================================================
def generate_corbata_row(r):
    """Return the r-th Corbata row as a list of integers."""
    if r == 1:
        return [1]

    L = (r - 2) ** 2 + 2

    if r == 2:
        return [L]

    row_len = r - 1
    row = [0] * row_len
    row[0] = L

    if r % 2 == 0:
        k = r // 2
        negcount = k - 1
        idx = 0

        for p in range(negcount - 1, -1, -1):
            idx += 1
            diff_val = -2 * (k + p - 1)
            row[idx] = row[idx - 1] + diff_val

        for t in range(1, negcount + 1):
            idx += 1
            diff_val = 2 * (k + t)
            row[idx] = row[idx - 1] + diff_val

    else:
        k = (r - 1) // 2
        idx = 0

        for s in range(k - 2, -1, -1):
            idx += 1
            diff_val = -2 * (k + s)
            row[idx] = row[idx - 1] + diff_val

        idx += 1
        row[idx] = row[idx - 1] + 2

        for t in range(2, k + 1):
            idx += 1
            diff_val = 2 * (k + t)
            row[idx] = row[idx - 1] + diff_val

    if len(row) != r - 1:
        raise ValueError("Row length mismatch at r = " + str(r))

    return row


# ================================================================
# Candidate generation baselines
# ================================================================
def random_odd_candidate(k, rng):
    n = rng.getrandbits(k)
    n |= (1 << (k - 1))
    n |= 1
    return n

def wheel_passes(n, wheel_primes):
    for p in wheel_primes:
        if n % p == 0:
            return False
    return True

def gcd_ok_for_e(n, e):
    return math.gcd(n - 1, e) == 1


# ================================================================
# Corbata residues builder (standalone)
# It does NOT call any external code, only functions above
# ================================================================
def build_corbata_residues(num_rows, modulus_M, mr_rounds, rng, require_coprime_e=False, e=65537):
    """
    Build residues as p mod M where p are primes found inside first num_rows corbata rows.
    Returns residues list and a small report dictionary.
    """
    max_r = num_rows + 1

    used_numbers = set([1])
    all_primes = []
    total_numbers = 0
    total_primes = 0

    for r in range(2, max_r + 1):
        row_vec = generate_corbata_row(r)

        parity = row_vec[0] % 2
        for x in row_vec:
            if (x % 2) != parity:
                raise ValueError("CNS error: row mixes parity at r = " + str(r))

        for v in row_vec:
            if v in used_numbers:
                raise ValueError("CNS error: repeated number " + str(v) + " at r = " + str(r))
            used_numbers.add(v)

        for val in row_vec:
            total_numbers += 1
            ok, _status = primality_status(val, mr_rounds=mr_rounds, rng=rng)
            if ok:
                if require_coprime_e and not gcd_ok_for_e(val, e):
                    continue
                total_primes += 1
                all_primes.append(val)

    all_primes = sorted(set(all_primes))
    residues = sorted({p % modulus_M for p in all_primes})

    report = {
        "rows_used": num_rows,
        "max_r": max_r,
        "total_numbers_scanned": total_numbers,
        "total_primes_found": total_primes,
        "unique_primes_used": len(all_primes),
        "modulus_M": modulus_M,
        "unique_residues": len(residues),
    }
    return residues, report


def corbata_candidate_from_residues(k, modulus_M, residues, rng):
    """
    Sample n in [2^(k-1), 2^k - 1] such that n ≡ residue (mod M).
    This matches the style in your script, but fully self contained.
    """
    low = 1 << (k - 1)
    high = (1 << k) - 1

    base = rng.randrange(low, high + 1)
    target_res = residues[rng.randrange(0, len(residues))]

    n = base - ((base - target_res) % modulus_M)

    if n < low:
        n += modulus_M * ((low - n + modulus_M - 1) // modulus_M)
    elif n > high:
        n -= modulus_M * ((n - high + modulus_M - 1) // modulus_M)

    n |= (1 << (k - 1))
    n |= 1
    return n


# ================================================================
# Prime finders for each method
# ================================================================
def find_prime_random(k, rng, mr_rounds, e, use_wheel, wheel_primes):
    attempts = 0
    while True:
        attempts += 1
        n = random_odd_candidate(k, rng)

        if use_wheel and not wheel_passes(n, wheel_primes):
            continue

        ok, _status = primality_status(n, mr_rounds=mr_rounds, rng=rng)
        if not ok:
            continue

        if not gcd_ok_for_e(n, e):
            continue

        return n, attempts


def find_prime_corbata_guided(k, rng, mr_rounds, e, modulus_M, residues, use_wheel, wheel_primes):
    attempts = 0
    while True:
        attempts += 1
        n = corbata_candidate_from_residues(k, modulus_M, residues, rng)

        if use_wheel and not wheel_passes(n, wheel_primes):
            continue

        ok, _status = primality_status(n, mr_rounds=mr_rounds, rng=rng)
        if not ok:
            continue

        if not gcd_ok_for_e(n, e):
            continue

        return n, attempts


# ================================================================
# Statistics helpers
# ================================================================
def percentile(data, p):
    if not data:
        return None
    xs = sorted(data)
    k = (len(xs) - 1) * (p / 100.0)
    f = int(math.floor(k))
    c = int(math.ceil(k))
    if f == c:
        return xs[f]
    return xs[f] * (c - k) + xs[c] * (k - f)

def summarize_attempts(attempts_list, seconds_list):
    a_mean = mean(attempts_list)
    a_std = pstdev(attempts_list) if len(attempts_list) > 1 else 0.0
    a_med = median(attempts_list)
    out = {
        "trials": len(attempts_list),
        "attempts_mean": a_mean,
        "attempts_std": a_std,
        "attempts_median": a_med,
        "attempts_p25": percentile(attempts_list, 25),
        "attempts_p75": percentile(attempts_list, 75),
        "attempts_p90": percentile(attempts_list, 90),
        "attempts_p99": percentile(attempts_list, 99),
        "time_mean_s": mean(seconds_list),
        "time_std_s": pstdev(seconds_list) if len(seconds_list) > 1 else 0.0,
        "time_median_s": median(seconds_list),
    }
    return out


def run_trials(trials, label, fn):
    attempts_list = []
    seconds_list = []
    for _ in range(trials):
        t0 = time.perf_counter()
        _n, a = fn()
        t1 = time.perf_counter()
        attempts_list.append(a)
        seconds_list.append(t1 - t0)

    stats = summarize_attempts(attempts_list, seconds_list)

    print()
    print(label)
    print("Trials =", stats["trials"])
    print("Attempts mean =", round(stats["attempts_mean"], 6), "std =", round(stats["attempts_std"], 6), "median =", stats["attempts_median"])
    print("Attempts p25 =", stats["attempts_p25"], "p75 =", stats["attempts_p75"], "p90 =", stats["attempts_p90"], "p99 =", stats["attempts_p99"])
    print("Time mean s =", round(stats["time_mean_s"], 6), "std =", round(stats["time_std_s"], 6), "median =", round(stats["time_median_s"], 6))
    return stats


# ================================================================
# Main
# ================================================================
def main():
    print("=== Corbata RSA fair benchmark ===")
    if HAS_GMPY2:
        print("Backend: gmpy2 enabled")
    elif HAS_SYMPY:
        print("Backend: sympy enabled")
    else:
        print("Backend: Miller Rabin only")

    print()
    seed = int(input("Random seed integer, example 12345: ").strip())
    rng = random.Random(seed)

    k = int(input("Bit length k for primes p and q, example 1024 or 2048: ").strip())
    trials = int(input("Number of independent trials M, example 1000 or 5000: ").strip())
    mr_rounds = int(input("Miller Rabin rounds, example 24 or 32: ").strip())

    e_in = input("Public exponent e, press Enter for 65537: ").strip()
    e = 65537 if e_in == "" else int(e_in)

    use_wheel_in = input("Use wheel small prime filter, type yes or no: ").strip().lower()
    use_wheel = use_wheel_in == "yes"
    wheel_primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31] if use_wheel else []

    print()
    print("Now choose modulus M for residue sampling")
    print("Option 1: 2310 = 2*3*5*7*11 (classic)")
    print("Option 2: 30030 = 2*3*5*7*11*13")
    print("Option 3: custom integer")
    opt = input("Choose 1, 2, or 3: ").strip()

    if opt == "1":
        modulus_M = 2310
    elif opt == "2":
        modulus_M = 30030
    else:
        modulus_M = int(input("Enter custom modulus M: ").strip())

    print()
    num_rows = int(input("How many Corbata rows to analyze for residues, example 6000: ").strip())

    print()
    print("Building residues from Corbata")
    residues, residue_report = build_corbata_residues(
        num_rows=num_rows,
        modulus_M=modulus_M,
        mr_rounds=mr_rounds,
        rng=rng,
        require_coprime_e=True,
        e=e
    )

    if not residues:
        print("No residues collected. Increase num_rows or adjust M.")
        return

    print("Residue build report")
    for kkey in residue_report:
        print(kkey, "=", residue_report[kkey])

    print("Residues count =", len(residues))

    # Heuristic baseline expected attempts Ntrad = k ln 2
    Ntrad = k * math.log(2.0)
    print()
    print("Heuristic expected attempts Ntrad = k ln 2 =", round(Ntrad, 6))

    # Define trial wrappers so all methods share same rng and parameters
    def trial_random_plain():
        return find_prime_random(k, rng, mr_rounds, e, use_wheel=False, wheel_primes=wheel_primes)

    def trial_random_wheel():
        return find_prime_random(k, rng, mr_rounds, e, use_wheel=use_wheel, wheel_primes=wheel_primes)

    def trial_corbata_plain():
        return find_prime_corbata_guided(k, rng, mr_rounds, e, modulus_M, residues, use_wheel=False, wheel_primes=wheel_primes)

    def trial_corbata_wheel():
        return find_prime_corbata_guided(k, rng, mr_rounds, e, modulus_M, residues, use_wheel=use_wheel, wheel_primes=wheel_primes)

    # Run benchmarks
    res_random_plain = run_trials(trials, "Method A: random odd sampling", trial_random_plain)

    if use_wheel:
        res_random_wheel = run_trials(trials, "Method B: random sampling with wheel", trial_random_wheel)
    else:
        res_random_wheel = None

    res_corbata_plain = run_trials(trials, "Method C: Corbata guided residue sampling", trial_corbata_plain)

    if use_wheel:
        res_corbata_wheel = run_trials(trials, "Method D: Corbata guided plus wheel", trial_corbata_wheel)
    else:
        res_corbata_wheel = None

    # Efficiency and diagnostic
    def eta(stat):
        return (Ntrad / stat["attempts_mean"]) if stat["attempts_mean"] > 0 else float("inf")

    print()
    print("Efficiency eta = Ntrad divided by mean attempts")
    print("eta A =", round(eta(res_random_plain), 6))
    if res_random_wheel is not None:
        print("eta B =", round(eta(res_random_wheel), 6))
    print("eta C =", round(eta(res_corbata_plain), 6))
    if res_corbata_wheel is not None:
        print("eta D =", round(eta(res_corbata_wheel), 6))

    # Diagnostic: show size of residue set and implied candidate fraction
    # If residues are uniform mod M, then candidate fraction approx len(residues)/M among residues classes
    frac_classes = len(residues) / float(modulus_M)
    print()
    print("Residue class fraction len(residues)/M =", round(frac_classes, 10))
    print("This indicates how much the search space is restricted before primality tests")

    # Save report
    out = {
        "seed": seed,
        "bit_length_k": k,
        "trials_M": trials,
        "mr_rounds": mr_rounds,
        "e": e,
        "use_wheel": use_wheel,
        "wheel_primes": wheel_primes,
        "modulus_M": modulus_M,
        "residue_report": residue_report,
        "residues_count": len(residues),
        "residue_class_fraction": frac_classes,
        "Ntrad_k_ln2": Ntrad,
        "method_A_random": res_random_plain,
        "method_B_random_wheel": res_random_wheel,
        "method_C_corbata": res_corbata_plain,
        "method_D_corbata_wheel": res_corbata_wheel,
        "eta_A": eta(res_random_plain),
        "eta_B": eta(res_random_wheel) if res_random_wheel is not None else None,
        "eta_C": eta(res_corbata_plain),
        "eta_D": eta(res_corbata_wheel) if res_corbata_wheel is not None else None,
    }

    filename = "rsa_benchmark_report.json"
    with open(filename, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)

    print()
    print("Saved report to", filename)
    print("Done")


if __name__ == "__main__":
    main()
