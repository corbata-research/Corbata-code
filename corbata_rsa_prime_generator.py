import random
import math

# ================================================================
#  Optional stronger primality backends
# ================================================================
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
#  Miller Rabin primality test (probabilistic)
# ================================================================
def is_probable_prime(n, k=24):
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

    for _ in range(k):
        a = random.randrange(2, n - 1)
        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False

    return True


def primality_status(n, mr_rounds=24):
    """
    Return (is_prime_bool, status_string).
    Tries:
      1 Miller Rabin (fast filter, probabilistic)
      2 gmpy2 is_prime if available (stronger, still probabilistic but very reliable)
      3 sympy isprime if available (strong, often deterministic in many cases, but treat as strong probable)
    """
    if n < 2:
        return False, "Composite (n < 2)"

    # First, fast MR filter
    if not is_probable_prime(n, k=mr_rounds):
        return False, f"Composite (fails Miller Rabin, rounds={mr_rounds})"

    # If we have gmpy2, use it
    if HAS_GMPY2:
        # gmpy2.is_prime returns:
        # 0 composite, 1 probably prime, 2 definitely prime (for some sizes it can prove)
        res = int(gmpy2.is_prime(n))
        if res == 0:
            return False, "Composite (gmpy2 witness)"
        if res == 2:
            return True, "Prime (gmpy2 proved)"
        return True, "Probably prime (gmpy2 strong test)"

    # Else try sympy
    if HAS_SYMPY:
        # sympy.isprime returns True if prime, False if composite
        # It uses strong methods, but for RSA we still label as strong probable unless you separately prove.
        if sp.isprime(n):
            return True, "Probably prime (sympy isprime)"
        return False, "Composite (sympy witness)"

    # Fall back to MR only
    return True, f"Probably prime (Miller Rabin only, rounds={mr_rounds})"


# ================================================================
#   Generate a single Corbata row (Python version)
# ================================================================
def generate_corbata_row_local(r):
    """Return the r-th Corbata row as a list of integers."""
    if r == 1:
        return [1]

    L = (r - 2)**2 + 2

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

    return row


# ================================================================
#   Build Corbata-based residue pattern modulo M
# ================================================================
def build_corbata_residues(num_rows):
    max_r = num_rows + 1

    used_numbers = set([1])
    all_primes = []

    for r in range(2, max_r + 1):
        row_vec = generate_corbata_row_local(r)

        parity = row_vec[0] % 2
        if any((x % 2) != parity for x in row_vec):
            raise ValueError(f"CNS error: row {r} mixes parity.")

        for v in row_vec:
            if v in used_numbers:
                raise ValueError(f"CNS error: repeated number {v} in row {r}")
            used_numbers.add(v)

        for val in row_vec:
            ok, _status = primality_status(val, mr_rounds=24)
            if ok:
                all_primes.append(val)

    all_primes = sorted(set(all_primes))
    print(f"[Residues] Collected {len(all_primes)} primes (based on tests available).")

    if not all_primes:
        raise ValueError("No primes found. Increase num_rows.")

    M = 2 * 3 * 5 * 7 * 11  # 2310
    residues = sorted({p % M for p in all_primes})

    residues = [
        r for r in residues
        if (r % 2 == 1) and (r % 3 != 0) and (r % 5 != 0)
    ]

    print(f"[Residues] Using {len(residues)} residues modulo {M}.")
    return M, residues


# ================================================================
#   Generate a Corbata-guided prime of bit_len
# ================================================================
def corbata_prime_from_residues(bit_len, M, residues, e=65537):
    low = 2**(bit_len - 1)
    high = 2**bit_len - 1

    small_primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]

    print(f"Searching for a Corbata-style prime ({bit_len} bits)...")
    attempts = 0

    while True:
        attempts += 1

        base = random.randint(low, high)

        target_res = random.choice(residues)
        n = base - ((base - target_res) % M)

        if n < low:
            n += M * math.ceil((low - n) / M)
        elif n > high:
            n -= M * math.ceil((n - high) / M)

        if not (low <= n <= high):
            continue

        if n % 2 == 0:
            continue

        if any(n % sp_ == 0 for sp_ in small_primes):
            continue

        ok, status = primality_status(n, mr_rounds=32)
        if not ok:
            continue

        if math.gcd(n - 1, e) != 1:
            continue

        print(f"Found candidate after {attempts} attempts.")
        print(f"Primality status: {status}")
        return n, attempts, status


# ================================================================
#   Extended Euclidean Algorithm (modular inverse)
# ================================================================
def mod_inverse(a, m):
    old_r, r = a, m
    old_s, s = 1, 0

    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s

    if old_r != 1:
        raise ValueError("No modular inverse exists.")

    return old_s % m


# ================================================================
#   MAIN — generate full RSA keypair using Corbata primes
# ================================================================
def corbata_rsa_prime_generator():
    print("=== Corbata-based RSA key generator ===")
    if HAS_GMPY2:
        print("[Backend] Using gmpy2 for stronger primality checks.")
    elif HAS_SYMPY:
        print("[Backend] Using sympy for stronger primality checks.")
    else:
        print("[Backend] Using Miller Rabin only (probabilistic). Install gmpy2 for best results.")

    while True:
        bit_len = int(input("Bit-length for p and q (>=8): "))
        if bit_len >= 8:
            break

    e_input = input("Public exponent e [default 65537]: ")
    e = 65537 if e_input.strip() == "" else int(e_input)

    while True:
        num_rows = int(input("How many Corbata rows to analyze? "))
        if num_rows > 0:
            break

    M, residues = build_corbata_residues(num_rows)

    print("\nGenerating p...")
    p, attempts_p, status_p = corbata_prime_from_residues(bit_len, M, residues, e)

    print("\nGenerating q...")
    q, attempts_q, status_q = corbata_prime_from_residues(bit_len, M, residues, e)
    while q == p:
        print("q equals p; regenerating q...")
        q, attempts_q, status_q = corbata_prime_from_residues(bit_len, M, residues, e)

    n = p * q
    phi = (p - 1) * (q - 1)
    d = mod_inverse(e, phi)

    print("\n=== Corbata-style RSA Key Generated ===")
    print(f"\nAttempts to find p: {attempts_p}")
    print(f"Attempts to find q: {attempts_q}")

    print(f"\np (bits ~{p.bit_length()}):\n{p}")
    print(f"p status: {status_p}")

    print(f"\nq (bits ~{q.bit_length()}):\n{q}")
    print(f"q status: {status_q}")

    print(f"\nn = p * q (bits ~{n.bit_length()}):\n{n}")
    print(f"\ne:\n{e}")
    print(f"\nd:\n{d}")

    msg = 123456789
    c = pow(msg, e, n)
    m_rec = pow(c, d, n)

    print("\n--- RSA test ---")
    print(f"m  = {msg}")
    print(f"c  = {c}")
    print(f"m' = {m_rec}")
    print("Success?", m_rec == msg)

    approx_bits_n = n.bit_length()
    filename = f"corbata_rsa_key_{approx_bits_n}bits.txt"

    try:
        with open(filename, "w", encoding="utf-8") as f:
            f.write("=== Corbata-style RSA key ===\n\n")
            f.write(f"Bit-length p  ≈ {p.bit_length()}\n")
            f.write(f"Bit-length q  ≈ {q.bit_length()}\n")
            f.write(f"Bit-length n  ≈ {n.bit_length()}\n\n")

            f.write(f"Attempts to find p: {attempts_p}\n")
            f.write(f"Attempts to find q: {attempts_q}\n\n")

            f.write(f"p status: {status_p}\n")
            f.write(f"q status: {status_q}\n\n")

            f.write("p =\n")
            f.write(str(p) + "\n\n")

            f.write("q =\n")
            f.write(str(q) + "\n\n")

            f.write("n = p * q =\n")
            f.write(str(n) + "\n\n")

            f.write("e =\n")
            f.write(str(e) + "\n\n")

            f.write("d =\n")
            f.write(str(d) + "\n\n")

            f.write("--- RSA test ---\n")
            f.write(f"m  = {msg}\n")
            f.write(f"c  = {c}\n")
            f.write(f"m' = {m_rec}\n")
            f.write(f"Success? {m_rec == msg}\n")

        print(f"\nKey saved to file: {filename}")
    except Exception as ex:
        print(f"\n[Warning] Could not save key to file: {ex}")

    print("\nDone.")
    return n, e, d, p, q


if __name__ == "__main__":
    corbata_rsa_prime_generator()
